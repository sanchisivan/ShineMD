# ShineMD main app
#
# This script keeps the workflow fully inside R/Shiny and reads AMBER prmtop +
# NetCDF trajectories from disk with bio3d + ncdf4. The expected workflow is:
# 1) point the app to a project folder,
# 2) choose and order the trajectory segments,
# 3) run the basic analysis,
# 4) launch the more specific membrane / interaction / clustering calculations
#    only when you actually need them.
#
# A small but important detail: ordered NetCDF segments can be combined into a
# continuous time axis, while still allowing the user to run on-demand analyses
# on only a subset of the loaded segments.

suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyjs)
  library(shinyWidgets)
  library(plotly)
  library(ggplot2)
  library(htmlwidgets)
  library(DT)
  library(shinyFiles)
  library(fs)
  library(bio3d)
library(ncdf4)
})

if (!requireNamespace("ncdf4", quietly = TRUE)) {
  stop("Missing dependency: ncdf4. Install with install.packages('ncdf4').")
}

`%||%` <- function(a,b) if (!is.null(a)) a else b

read_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(read.csv(path), error = function(e) NULL)
}

# Keep all exported files in one project-local ShineMD folder.
# That makes the output path predictable and avoids the old results_bio3d label.
project_results_dir <- function(project_dir) file.path(project_dir, "results_ShineMD")
cluster_results_dir <- function(project_dir) file.path(project_results_dir(project_dir), "clustering_rmsd")

normalize_pca_var_df <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) < 1) return(df)
  nms <- names(df)
  ren <- c(prop.var = 'prop_var', cum.var = 'cum_var', PropVar = 'prop_var', CumVar = 'cum_var',
           explained_variance = 'prop_var', cumulative_variance = 'cum_var', component = 'PC')
  for (from in names(ren)) {
    if (from %in% nms && !(ren[[from]] %in% nms)) names(df)[names(df) == from] <- ren[[from]]
  }
  if (!('PC' %in% names(df))) {
    if ('pc' %in% names(df)) names(df)[names(df) == 'pc'] <- 'PC'
  }
  if (!('PC' %in% names(df))) df$PC <- paste0('PC', seq_len(nrow(df)))
  if (!('prop_var' %in% names(df)) && 'variance' %in% names(df)) {
    vv <- suppressWarnings(as.numeric(df$variance))
    s <- sum(vv, na.rm = TRUE)
    if (is.finite(s) && s > 0) df$prop_var <- vv / s
  }
  if (!('cum_var' %in% names(df)) && 'prop_var' %in% names(df)) {
    pv <- suppressWarnings(as.numeric(df$prop_var))
    df$cum_var <- cumsum(replace(pv, !is.finite(pv), 0))
  }
  if ('prop_var' %in% names(df)) df$prop_var <- suppressWarnings(as.numeric(df$prop_var))
  if ('cum_var' %in% names(df)) df$cum_var <- suppressWarnings(as.numeric(df$cum_var))
  pc_raw <- as.character(df$PC)
  pc_idx <- suppressWarnings(as.integer(sub('^PC([0-9]+)$', '\\1', pc_raw)))
  bad <- is.na(pc_idx)
  if (any(bad)) pc_idx[bad] <- seq_len(sum(bad))
  df$PC <- ifelse(grepl('^PC[0-9]+$', pc_raw), pc_raw, paste0('PC', pc_idx))
  df$PC_index <- pc_idx
  df
}

# ─── Custom CSS (preserved from your uploaded app.R) ──────────────────────────
custom_css <- "
@import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;600;700&family=IBM+Plex+Sans:wght@300;400;500;600;700&display=swap');

:root {
  --bg-primary: #0a0e17;
  --bg-secondary: #111827;
  --bg-card: #1a2235;
  --bg-card-hover: #1f2a40;
  --accent-cyan: #06d6a0;
  --accent-blue: #118ab2;
  --accent-warm: #ef476f;
  --accent-gold: #ffd166;
  --text-primary: #e8ecf1;
  --text-muted: #8899aa;
  --border-color: #2a3a52;
  --font-body: 'IBM Plex Sans', sans-serif;
  --font-mono: 'JetBrains Mono', monospace;
}

body, .wrapper, .content-wrapper, .main-sidebar, .left-side,
.skin-blue .main-header .navbar, .skin-blue .main-header .logo {
  font-family: var(--font-body) !important;
}

.skin-blue .main-header .logo {
  background: linear-gradient(135deg, #0d1520, #162035) !important;
  color: var(--accent-cyan) !important;
  font-family: var(--font-mono) !important;
  font-weight: 700 !important;
  font-size: 16px !important;
  letter-spacing: 1px;
  border-bottom: 2px solid var(--accent-cyan) !important;
}
.skin-blue .main-header .logo:hover {
  background: linear-gradient(135deg, #111d2e, #1a2840) !important;
}

.skin-blue .main-header .navbar {
  background: var(--bg-primary) !important;
  border-bottom: 1px solid var(--border-color) !important;
}

.skin-blue .main-sidebar, .skin-blue .left-side {
  background: var(--bg-secondary) !important;
  border-right: 1px solid var(--border-color);
}
.skin-blue .sidebar-menu > li > a {
  color: var(--text-muted) !important;
  border-left: 3px solid transparent;
  font-weight: 500;
  transition: all 0.2s ease;
}
.skin-blue .sidebar-menu > li > a:hover {
  background: var(--bg-card) !important;
  color: var(--text-primary) !important;
  border-left-color: var(--accent-cyan);
}
.skin-blue .sidebar-menu > li.active > a {
  background: var(--bg-card) !important;
  color: var(--accent-cyan) !important;
  border-left-color: var(--accent-cyan);
  font-weight: 600;
}
.skin-blue .sidebar-menu > li > a > .fa {
  color: var(--accent-blue) !important;
}
.skin-blue .sidebar-menu > li.active > a > .fa {
  color: var(--accent-cyan) !important;
}

.content-wrapper, .right-side {
  background: var(--bg-primary) !important;
}
.content-wrapper .content {
  padding: 20px;
}

.box {
  background: var(--bg-card) !important;
  border: 1px solid var(--border-color) !important;
  border-radius: 10px !important;
  box-shadow: 0 4px 20px rgba(0,0,0,0.3) !important;
  color: var(--text-primary) !important;
}
.box .box-header {
  border-bottom: 1px solid var(--border-color) !important;
  padding: 12px 16px !important;
}
.box .box-header .box-title {
  color: var(--text-primary) !important;
  font-weight: 600 !important;
  font-size: 14px !important;
  letter-spacing: 0.3px;
}
.box.box-primary { border-top: 3px solid var(--accent-cyan) !important; }
.box.box-warning { border-top: 3px solid var(--accent-gold) !important; }
.box.box-danger  { border-top: 3px solid var(--accent-warm) !important; }
.box.box-info    { border-top: 3px solid var(--accent-blue) !important; }

h1, h2, h3, h4, h5, h6, p, span, label, .form-group label {
  color: var(--text-primary) !important;
}

.form-control, .selectize-input, .selectize-dropdown {
  background: var(--bg-secondary) !important;
  color: var(--text-primary) !important;
  border: 1px solid var(--border-color) !important;
  border-radius: 6px !important;
  font-family: var(--font-body) !important;
}
/* Selected pills / tags inside multi-select */
.selectize-input .item {
  background: rgba(6,214,160,0.16) !important;
  color: var(--text-primary) !important;
  border: 1px solid rgba(6,214,160,0.38) !important;
  border-radius: 4px !important;
  padding: 2px 6px !important;
  font-size: 12px !important;
  margin: 1px 2px !important;
}
.selectize-input .item.active {
  background: rgba(6,214,160,0.30) !important;
  color: #fff !important;
  border-color: rgba(6,214,160,0.60) !important;
}
/* Remove button (×) inside each pill */
.selectize-input .item .remove {
  color: var(--accent-cyan) !important;
  border-left: 1px solid rgba(6,214,160,0.30) !important;
  margin-left: 4px !important;
  padding-left: 4px !important;
}
.selectize-input .item .remove:hover {
  color: var(--accent-warm) !important;
  background: transparent !important;
}
/* Dropdown option rows */
.selectize-dropdown .option {
  background: var(--bg-secondary) !important;
  color: var(--text-primary) !important;
  padding: 6px 10px !important;
}
.selectize-dropdown .option:hover,
.selectize-dropdown .option.active {
  background: rgba(6,214,160,0.14) !important;
  color: #fff !important;
}
/* Placeholder text */
.selectize-input.has-items input::placeholder,
.selectize-input input::placeholder {
  color: var(--text-muted) !important;
}
.form-control:focus {
  border-color: var(--accent-cyan) !important;
  box-shadow: 0 0 0 2px rgba(6, 214, 160, 0.15) !important;
}

.btn, .btn-default, .btn-primary, .btn-success, .btn-info, .shiny-download-link {
  background: linear-gradient(135deg, var(--accent-blue), var(--accent-cyan)) !important;
  border: 1px solid rgba(6,214,160,0.28) !important;
  color: #fff !important;
  font-weight: 700 !important;
  border-radius: 8px !important;
  padding: 10px 22px !important;
  transition: all 0.18s ease !important;
  text-transform: uppercase;
  font-size: 12px !important;
  letter-spacing: 0.9px;
  box-shadow: 0 6px 18px rgba(6, 214, 160, 0.18) !important;
}
button.btn, a.btn, .btn.action-button, .shiny-download-link.btn {
  display: inline-flex !important;
  align-items: center !important;
  justify-content: center !important;
  gap: 6px;
}
.btn:hover, .btn:focus, .btn-default:hover, .btn-primary:hover, .btn-success:hover,
.shiny-download-link:hover, .shiny-download-link:focus {
  color: #fff !important;
  border-color: rgba(255,255,255,0.18) !important;
  transform: translateY(-1px);
  box-shadow: 0 10px 24px rgba(6, 214, 160, 0.24) !important;
}
.btn:focus, .shiny-download-link:focus {
  outline: none !important;
}
.btn-default, .btn.btn-default, .btn.btn-primary, .btn.btn-success,
.shiny-download-link.btn-default, .shiny-download-link.btn-primary {
  color: #fff !important;
}
.btn.disabled, .btn[disabled], fieldset[disabled] .btn,
.shiny-download-link.disabled, .shiny-download-link[disabled] {
  opacity: 0.5 !important;
  box-shadow: none !important;
  cursor: not-allowed !important;
}

.btn-compact {
  padding: 7px 16px !important;
  font-size: 11px !important;
  letter-spacing: 0.6px;
}

.nav-tabs {
  border-bottom: 1px solid var(--border-color) !important;
}
.nav-tabs > li > a {
  color: var(--text-muted) !important;
  border: none !important;
  font-weight: 500;
}
.nav-tabs > li.active > a, .nav-tabs > li.active > a:hover {
  background: var(--bg-card) !important;
  color: var(--accent-cyan) !important;
  border-bottom: 2px solid var(--accent-cyan) !important;
}

.dataTables_wrapper {
  color: var(--text-primary) !important;
}
table.dataTable {
  color: var(--text-primary) !important;
  background: var(--bg-card) !important;
}
table.dataTable thead th {
  background: var(--bg-secondary) !important;
  color: var(--accent-cyan) !important;
  border-bottom: 1px solid var(--border-color) !important;
  font-family: var(--font-mono) !important;
  font-size: 12px !important;
}
table.dataTable tbody td {
  border-bottom: 1px solid var(--border-color) !important;
  font-size: 13px;
}
table.dataTable tbody tr:hover {
  background: var(--bg-card-hover) !important;
}

.shiny-notification {
  background: var(--bg-card) !important;
  color: var(--text-primary) !important;
  border: 1px solid var(--accent-cyan) !important;
  border-radius: 8px !important;
}

.modal-content, .themed-modal .modal-content {
  background: linear-gradient(180deg, rgba(24, 34, 58, 0.98), rgba(12, 18, 33, 0.98)) !important;
  color: var(--text-primary) !important;
  border: 1px solid rgba(6,214,160,0.35) !important;
  border-radius: 14px !important;
  box-shadow: 0 20px 55px rgba(0,0,0,0.55) !important;
}
.modal-header, .modal-footer,
.themed-modal .modal-header,
.modal-footer,
.themed-modal .modal-footer {
  background: rgba(17,138,178,0.08) !important;
  border-color: rgba(6,214,160,0.18) !important;
}
.modal-header,
.themed-modal .modal-header {
  border-top-left-radius: 14px !important;
  border-top-right-radius: 14px !important;
}
.themed-modal .modal-footer {
  border-bottom-left-radius: 14px !important;
  border-bottom-right-radius: 14px !important;
}
.modal-title, .modal-body, .modal-body p, .modal-body li, .modal-body div, .modal-body span,
.themed-modal .modal-title,
.themed-modal .modal-body,
.themed-modal .modal-body p,
.themed-modal .modal-body li,
.themed-modal .modal-body div,
.themed-modal .modal-body span {
  color: var(--text-primary) !important;
}
.modal .close,
.themed-modal .close {
  color: var(--text-primary) !important;
  opacity: 0.8 !important;
}
.completion-hero {
  background: linear-gradient(135deg, rgba(17,138,178,0.20), rgba(6,214,160,0.14)) !important;
  border: 1px solid rgba(6,214,160,0.26) !important;
  border-radius: 12px !important;
  padding: 16px 18px !important;
  margin-bottom: 16px !important;
}
.completion-hero-title {
  font-size: 17px !important;
  font-weight: 700 !important;
  color: #fff !important;
  margin-bottom: 4px !important;
}
.completion-hero-subtitle {
  color: var(--text-secondary) !important;
  line-height: 1.55 !important;
}
.completion-chip-label {
  margin-top: 10px !important;
  margin-bottom: 6px !important;
  font-size: 11px !important;
  color: var(--text-muted) !important;
  text-transform: uppercase !important;
  letter-spacing: 0.9px !important;
}
.completion-path-chip {
  display: block !important;
  background: rgba(2, 8, 23, 0.45) !important;
  border: 1px solid rgba(6,214,160,0.18) !important;
  border-radius: 10px !important;
  padding: 10px 12px !important;
  font-family: var(--font-mono) !important;
  font-size: 12px !important;
  line-height: 1.5 !important;
  word-break: break-all !important;
}
.completion-steps {
  margin-top: 16px !important;
}
.completion-step {
  display: flex !important;
  gap: 12px !important;
  align-items: flex-start !important;
  padding: 10px 0 !important;
  border-top: 1px solid rgba(6,214,160,0.10) !important;
}
.completion-step:first-child {
  border-top: none !important;
}
.completion-step-index {
  width: 24px !important;
  height: 24px !important;
  min-width: 24px !important;
  border-radius: 999px !important;
  background: linear-gradient(135deg, var(--accent-blue), var(--accent-cyan)) !important;
  color: #fff !important;
  font-weight: 700 !important;
  font-size: 12px !important;
  display: inline-flex !important;
  align-items: center !important;
  justify-content: center !important;
  margin-top: 2px !important;
}
.completion-step-title {
  color: #fff !important;
  font-weight: 600 !important;
  margin-bottom: 2px !important;
}
.completion-step-text {
  color: var(--text-secondary) !important;
  line-height: 1.5 !important;
}

.status-badge {
  display: inline-block;
  padding: 4px 12px;
  border-radius: 20px;
  font-size: 12px;
  font-weight: 600;
  font-family: var(--font-mono);
  letter-spacing: 0.5px;
}
.status-loaded { background: rgba(6,214,160,0.15); color: var(--accent-cyan); }
.status-empty  { background: rgba(136,153,170,0.15); color: var(--text-muted); }

.metric-card {
  background: var(--bg-secondary);
  border: 1px solid var(--border-color);
  border-radius: 8px;
  padding: 16px;
  text-align: center;
  margin-bottom: 12px;
}
.metric-card .metric-value {
  font-family: var(--font-mono);
  font-size: 24px;
  font-weight: 700;
  color: var(--accent-cyan);
}
.metric-card .metric-label {
  font-size: 11px;
  color: var(--text-muted);
  text-transform: uppercase;
  letter-spacing: 1px;
  margin-top: 4px;
}

.section-header {
  color: var(--text-primary);
  font-family: var(--font-mono);
  font-size: 12px;
  text-transform: uppercase;
  letter-spacing: 1.5px;
  margin-bottom: 16px;
  padding-bottom: 8px;
  border-bottom: 1px solid var(--border-color);
}

.sidebar-section { margin-bottom: 20px; }

.skin-blue .main-header .sidebar-toggle { color: var(--text-muted) !important; }
.skin-blue .main-header .sidebar-toggle:hover { color: var(--accent-cyan) !important; }

.treeview-menu { background: var(--bg-primary) !important; }
.treeview-menu > li > a { color: var(--text-muted) !important; font-size: 13px; }
.treeview-menu > li > a:hover { color: var(--text-primary) !important; }


  
  /* Collapsible per-plot styling/export panels */
  details.plot-style {
    border: 1px solid var(--border-color);
    border-radius: 10px;
    padding: 10px 12px;
    background: rgba(0,0,0,0.12);
    margin-top: 8px;
  }
  details.plot-style > summary {
    cursor: pointer;
    font-weight: 600;
    color: var(--accent-gold);
    outline: none;
    list-style: none;
  }
  details.plot-style > summary::-webkit-details-marker { display: none; }
  details.plot-style[open] > summary { margin-bottom: 10px; }
  details.plot-style .form-group { margin-bottom: 10px; }

  /* Download buttons stack (avoid overlap + match spacing) */
  .download-stack {
    display: flex;
    flex-direction: column;
    gap: 10px;
    margin-top: 10px;
  }
  .download-stack .btn,
  .download-stack .shiny-download-link {
    width: 100% !important;
  }

  /* Small info buttons next to plot titles */
  .info-btn.btn {
    background: rgba(255, 209, 102, 0.10) !important;
    border: 1px solid rgba(255, 209, 102, 0.30) !important;
    color: var(--accent-gold) !important;
    border-radius: 999px !important;
    padding: 6px 10px !important;
    font-size: 11px !important;
    letter-spacing: 0.2px !important;
    text-transform: none !important;
    box-shadow: none !important;
    width: auto !important;
    min-width: auto !important;
  }
  .info-btn.btn:hover,
  .info-btn.btn:focus {
    background: rgba(255, 209, 102, 0.16) !important;
    border-color: rgba(255, 209, 102, 0.45) !important;
    transform: translateY(-1px);
  }

  .modal-backdrop.in { opacity: 0.72 !important; }
  .modal-content { background-clip: padding-box !important; }
  .modal-body ul { padding-left: 22px !important; }
  .modal-body li { line-height: 1.5 !important; margin-bottom: 4px !important; }
  .modal-body p, .modal-body li { color: #e8ecf1 !important; }
  .waiter-modal .modal-dialog { margin-top: 12vh; max-width: 560px; }
  .waiter-shell { text-align: center; padding: 14px 10px 8px 10px; }
  .waiter-spinner { font-size: 38px; color: var(--accent-cyan); margin-bottom: 12px; }
  .waiter-title {
    color: var(--text-primary);
    font-size: 22px;
    font-weight: 700;
    margin-bottom: 8px;
  }
  .waiter-subtitle {
    color: var(--text-muted);
    font-size: 14px;
    line-height: 1.6;
    max-width: 440px;
    margin: 0 auto;
  }

  /* Verbatim / pre text (log, session info) */
  pre, .shiny-text-output {
    background: var(--bg-primary) !important;
    color: var(--accent-cyan) !important;
    border: 1px solid var(--border-color) !important;
    border-radius: 8px !important;
    font-family: var(--font-mono) !important;
    font-size: 12px !important;
    padding: 12px !important;
    max-height: 360px;
    overflow-y: auto;
  }

  /* Custom scrollbar */
  ::-webkit-scrollbar { width: 8px; height: 8px; }
  ::-webkit-scrollbar-track { background: var(--bg-secondary); border-radius: 4px; }
  ::-webkit-scrollbar-thumb { background: var(--border-color); border-radius: 4px; }
  ::-webkit-scrollbar-thumb:hover { background: var(--accent-blue); }

  /* File input styling */
  .btn-file { background: var(--bg-secondary) !important; border: 1px solid var(--border-color) !important; color: var(--text-primary) !important; }
  .progress { background: var(--bg-secondary) !important; border-radius: 4px !important; }
  .progress-bar { background: linear-gradient(90deg, var(--accent-blue), var(--accent-cyan)) !important; }

  /* Summary stat row under plots */
  .stat-row {
    display: flex;
    gap: 12px;
    margin-top: 12px;
    flex-wrap: wrap;
  }
  .stat-row .metric-card {
    flex: 1;
    min-width: 120px;
    padding: 10px 12px;
    margin-bottom: 0;
  }
  .stat-row .metric-value { font-size: 18px; }
  .stat-row .metric-label { font-size: 10px; }

  /* Empty plot placeholder */
  .empty-plot-msg {
    display: flex;
    align-items: center;
    justify-content: center;
    height: 300px;
    color: var(--text-muted);
    font-size: 15px;
    font-style: italic;
    border: 1px dashed var(--border-color);
    border-radius: 10px;
    background: rgba(0,0,0,0.08);
  }
"

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────
natural_sort_nc <- function(paths) {
  f <- basename(paths)
  n <- suppressWarnings(as.integer(sub(".*?(\\d+)\\.nc$", "\\1", f, perl=TRUE)))
  n[is.na(n)] <- 1e9
  paths[order(n, f)]
}

detect_inputs <- function(project_dir) {
  prmtops <- dir_ls(project_dir, regexp="\\.(prmtop|parm7)$", type="file")
  prmtop <- if (length(prmtops) >= 1) prmtops[[1]] else NA_character_
  ncs <- dir_ls(project_dir, regexp="\\.nc$", type="file")
  ncs <- natural_sort_nc(ncs)
  list(prmtop = prmtop, trajs = ncs)
}

parse_resno <- function(txt) {
  txt <- gsub("\\s+", "", as.character(txt %||% ""))
  if (!nzchar(txt)) return(integer())
  parts <- unlist(strsplit(txt, ",", fixed = TRUE))
  out <- integer()
  for (p in parts) {
    if (!nzchar(p)) next
    if (grepl(":", p, fixed = TRUE)) {
      ab <- strsplit(p, ":", fixed = TRUE)[[1]]
      if (length(ab) == 2) {
        a <- suppressWarnings(as.integer(ab[1])); b <- suppressWarnings(as.integer(ab[2]))
        if (!is.na(a) && !is.na(b)) out <- c(out, seq(a, b))
      }
    } else {
      a <- suppressWarnings(as.integer(p))
      if (!is.na(a)) out <- c(out, a)
    }
  }
  unique(out)
}


parse_csv_tokens <- function(txt) {
  txt <- gsub("\\s+", "", txt)
  if (!nzchar(txt)) return(character())
  unlist(strsplit(txt, ",", fixed=TRUE))
}

# Compatibility aliases (older blocks may call these)
parse_resno_range <- function(txt) parse_resno(txt)
parse_resid_list  <- function(txt) parse_csv_tokens(txt)

# ─────────────────────────────────────────────────────────────────────────────
# Helper functions used throughout the app. These keep selections NA-safe and avoid
# version-specific surprises when computing RMSD with bio3d.
# ─────────────────────────────────────────────────────────────────────────────

get_prm_natom <- function(prm) {
  pdb_like <- try(as.pdb(prm), silent = TRUE)
  if (!inherits(pdb_like, "try-error") && !is.null(pdb_like$atom) && nrow(pdb_like$atom) > 0) {
    return(nrow(pdb_like$atom))
  }
  if (!is.null(prm$atom) && is.data.frame(prm$atom) && nrow(prm$atom) > 0) {
    return(nrow(prm$atom))
  }
  NA_integer_
}

clean_atom_inds <- function(x, natom = NA_integer_) {
  x <- suppressWarnings(as.integer(x))
  x <- x[!is.na(x)]
  if (!is.na(natom)) x <- x[x >= 1L & x <= natom]
  unique(x)
}

auto_resnos_from_calpha <- function(prm, exclude_resid = character()) {
  pdb_like <- try(as.pdb(prm), silent = TRUE)
  if (inherits(pdb_like, "try-error") || is.null(pdb_like$atom$resno) || is.null(pdb_like$atom$elety)) {
    return(integer())
  }
  a <- pdb_like$atom
  an <- trimws(a$elety %||% "")
  keep <- an == "CA"
  if (!is.null(a$resid) && length(exclude_resid) > 0) keep <- keep & !(a$resid %in% exclude_resid)
  resnos <- unique(a$resno[keep])
  resnos <- resnos[!is.na(resnos)]
  sort(resnos)
}
auto_resnos_from_backbone <- function(prm, exclude_resid = character()) {
  pdb_like <- try(as.pdb(prm), silent = TRUE)
  if (inherits(pdb_like, "try-error") || is.null(pdb_like$atom$resno) || is.null(pdb_like$atom$elety)) {
    return(integer())
  }
  a <- pdb_like$atom
  if (!is.null(a$resid) && length(exclude_resid) > 0) {
    a <- a[!(a$resid %in% exclude_resid), , drop = FALSE]
  }
  an <- trimws(a$elety %||% "")
  resnos <- unique(a$resno[!is.na(a$resno)])
  if (length(resnos) < 1) return(integer())
  good <- vapply(resnos, function(rn) {
    aa <- an[a$resno == rn]
    any(aa == "CA") && (any(aa == "N") || any(aa == "C"))
  }, logical(1))
  sort(resnos[good])
}



# Fast RMSD (Å) against one reference vector. This avoids relying on a bio3d
# calling convention that has changed across versions.
rmsd_to_ref <- function(xyz_mat, ref_vec) {
  if (is.null(dim(xyz_mat))) stop("xyz_mat must be a matrix (frames x coords)")
  if (length(ref_vec) != ncol(xyz_mat)) stop("ref_vec length mismatch with xyz_mat columns")
  n_atoms <- ncol(xyz_mat) / 3
  if (n_atoms %% 1 != 0) stop("Number of coordinate columns is not a multiple of 3")
  dif <- xyz_mat - matrix(ref_vec, nrow = nrow(xyz_mat), ncol = ncol(xyz_mat), byrow = TRUE)
  sqrt(rowSums(dif * dif) / n_atoms)
}

# Kabsch superposition — drop-in replacement for bio3d::fit.xyz which silently
# fails in some bio3d versions.
# fixed      : numeric vector (one reference frame, length = 3*N_total)
# mobile     : numeric matrix  (nframes × 3*N_total)
# fixed.inds : integer vector of COORDINATE-column indices used for fitting
# mobile.inds: integer vector of COORDINATE-column indices used for fitting
# Returns the full mobile matrix with ALL atoms rotated/translated.
kabsch_fit_xyz <- function(fixed, mobile, fixed.inds = NULL, mobile.inds = NULL) {
  if (is.null(dim(mobile))) mobile <- matrix(mobile, nrow = 1)
  if (is.matrix(fixed)) fixed <- as.numeric(fixed[1, ])

  nc <- ncol(mobile)
  if (is.null(fixed.inds))  fixed.inds  <- seq_len(nc)
  if (is.null(mobile.inds)) mobile.inds <- seq_len(nc)

  # Reference fitting atoms as N_fit × 3
  ref_fit <- matrix(fixed[fixed.inds], ncol = 3, byrow = TRUE)
  ref_cen <- colMeans(ref_fit)
  ref_c   <- sweep(ref_fit, 2, ref_cen)          # centered reference

  nF     <- nrow(mobile)
  result <- mobile                                 # output matrix (copy)

  for (i in seq_len(nF)) {
    # Mobile fitting atoms for frame i — N_fit × 3
    mob_fit <- matrix(mobile[i, mobile.inds], ncol = 3, byrow = TRUE)
    mob_cen <- colMeans(mob_fit)
    mob_c   <- sweep(mob_fit, 2, mob_cen)          # centered mobile

    # Cross-covariance & SVD
    H   <- crossprod(mob_c, ref_c)                 # 3×3
    sv  <- svd(H)
    d   <- det(sv$v %*% t(sv$u))
    sgn <- diag(c(1, 1, sign(d)))                  # avoid reflections
    R   <- sv$v %*% sgn %*% t(sv$u)               # 3×3 rotation

    # Apply to ALL atoms of this frame
    all_xyz  <- matrix(mobile[i, ], ncol = 3, byrow = TRUE)   # N_total × 3
    all_xyz  <- sweep(all_xyz, 2, mob_cen)                     # center
    all_xyz  <- all_xyz %*% t(R)                               # rotate
    all_xyz  <- sweep(all_xyz, 2, ref_cen, "+")                # translate

    result[i, ] <- as.vector(t(all_xyz))
  }
  result
}

# Running mean (centered) — returns NA for edges where the window doesn't fit
running_mean <- function(x, k) {
  if (k < 2 || length(x) < k) return(x)
  k <- as.integer(k)
  as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))
}

# Free Energy Landscape from 2D histogram: G(x,y) = -kT * ln(P(x,y))
compute_fel <- function(x, y, nbins = 80, temp_K = 300) {
  kB <- 1.380649e-23  # J/K
  NA_val <- 6.02214076e23  # Avogadro for kcal conversion
  kT_kcal <- kB * temp_K * NA_val / 4184  # kcal/mol

  xr <- range(x, na.rm = TRUE); yr <- range(y, na.rm = TRUE)
  xbreaks <- seq(xr[1], xr[2], length.out = nbins + 1)
  ybreaks <- seq(yr[1], yr[2], length.out = nbins + 1)
  xmid <- (head(xbreaks, -1) + tail(xbreaks, -1)) / 2
  ymid <- (head(ybreaks, -1) + tail(ybreaks, -1)) / 2

  h <- matrix(0, nrow = nbins, ncol = nbins)
  xi <- findInterval(x, xbreaks, all.inside = TRUE)
  yi <- findInterval(y, ybreaks, all.inside = TRUE)
  for (k in seq_along(x)) h[xi[k], yi[k]] <- h[xi[k], yi[k]] + 1L
  h[h == 0] <- NA_real_

  prob <- h / sum(h, na.rm = TRUE)
  G <- -kT_kcal * log(prob)
  G <- G - min(G, na.rm = TRUE)  # shift so minimum = 0

  expand.grid(PC1 = xmid, PC2 = ymid, KEEP.OUT.ATTRS = FALSE) |>
    transform(G_kcal = as.vector(G))
}

# Compute a cpptraj-like RMSD series: fit on one atom set and calculate RMSD on
# the same atom set (or a second set if explicitly requested).
aligned_rmsd_series <- function(xyz_mat, ref_vec, fit_cols, calc_cols = fit_cols) {
  if (is.null(dim(xyz_mat))) stop("xyz_mat must be a matrix (frames x coords)")
  nF <- nrow(xyz_mat)
  if (length(fit_cols) < 3 || length(calc_cols) < 3) return(rep(NA_real_, nF))
  xyz_fit <- kabsch_fit_xyz(fixed = ref_vec, mobile = xyz_mat,
                            fixed.inds = fit_cols, mobile.inds = fit_cols)
  rmsd_to_ref(xyz_fit[, calc_cols, drop = FALSE], ref_vec[calc_cols])
}


make_selections <- function(prm,
                            pep_resno,
                            lipo_resno,
                            membrane_resid,
                            exclude_resid = character(),
                            align_elety = c("N","CA","C","O")) {

  natom <- get_prm_natom(prm)

  # Auto-detect Region A residues if blank (protein/peptide-like backbone first, then CA-only fallback)
  if (length(pep_resno) < 1) {
    pep_resno <- auto_resnos_from_backbone(prm, exclude_resid = unique(c(exclude_resid, membrane_resid)))
    if (length(pep_resno) < 1) {
      pep_resno <- auto_resnos_from_calpha(prm, exclude_resid = unique(c(exclude_resid, membrane_resid)))
    }
  }
  if (length(pep_resno) < 1) {
    stop("Selection A residues are empty and auto-detection failed. Please provide residue numbers for Selection A.")
  }

  # Alignment selection: Region A backbone (preferred), fallback to global backbone-like atoms
  sel_align <- atom.select(prm, resno = pep_resno, elety = align_elety)
  if (length(sel_align$atom) < 3) {
    sel_align <- atom.select(prm, elety = align_elety)
  }
  if (length(sel_align$atom) < 3) {
    # Final fallbacks: use any non-H atoms, then any atoms.
    sel_align <- atom.select(prm, "noh")
  }
  if (length(sel_align$atom) < 3) {
    sel_align <- atom.select(prm, "all")
  }

  # Region A selections
  sel_A_ca    <- atom.select(prm, resno = pep_resno, elety = "CA")
  sel_A_bb    <- atom.select(prm, resno = pep_resno, elety = c("N","CA","C","O"))
  sel_A_bb3   <- atom.select(prm, resno = pep_resno, elety = c("N","CA","C"))
  sel_A_heavy <- atom.select(prm, "noh", resno = pep_resno)
  sel_A_all   <- atom.select(prm, resno = pep_resno)

  # Region B selections (optional)
  sel_B_heavy <- if (length(lipo_resno) > 0) atom.select(prm, "noh", resno = lipo_resno) else list(atom = integer())
  sel_B_all   <- if (length(lipo_resno) > 0) atom.select(prm, resno = lipo_resno) else list(atom = integer())
  sel_B_bb    <- if (length(lipo_resno) > 0) atom.select(prm, resno = lipo_resno, elety = c("N","CA","C","O")) else list(atom = integer())
  sel_B_ca    <- if (length(lipo_resno) > 0) atom.select(prm, resno = lipo_resno, elety = "CA") else list(atom = integer())

  # Combined selections (A + B)
  sel_AplusB_heavy <- list(atom = sort(unique(c(sel_A_heavy$atom, sel_B_heavy$atom))))
  sel_AplusB_bb    <- list(atom = sort(unique(c(sel_A_bb$atom,    sel_B_bb$atom))))
  sel_AplusB_ca    <- list(atom = sort(unique(c(sel_A_ca$atom,    sel_B_ca$atom))))
  sel_AplusB_all   <- list(atom = sort(unique(c(sel_A_all$atom,   sel_B_all$atom))))

  # Membrane selection by residue name (bio3d uses 'resid' for residue name)
  sel_mem <- if (length(membrane_resid) > 0) atom.select(prm, resid = membrane_resid) else list(atom = integer())

  # Clean indices (avoid NA -> NetCDF overflow)
  for (nm in c("atom")) {
    sel_align[[nm]] <- clean_atom_inds(sel_align[[nm]], natom)
    sel_A_ca[[nm]] <- clean_atom_inds(sel_A_ca[[nm]], natom)
    sel_A_bb[[nm]] <- clean_atom_inds(sel_A_bb[[nm]], natom)
    sel_A_bb3[[nm]] <- clean_atom_inds(sel_A_bb3[[nm]], natom)
    sel_A_heavy[[nm]] <- clean_atom_inds(sel_A_heavy[[nm]], natom)
    sel_A_all[[nm]] <- clean_atom_inds(sel_A_all[[nm]], natom)
    sel_B_heavy[[nm]] <- clean_atom_inds(sel_B_heavy[[nm]], natom)
    sel_B_all[[nm]] <- clean_atom_inds(sel_B_all[[nm]], natom)
    sel_B_bb[[nm]] <- clean_atom_inds(sel_B_bb[[nm]], natom)
    sel_B_ca[[nm]] <- clean_atom_inds(sel_B_ca[[nm]], natom)
    sel_AplusB_heavy[[nm]] <- clean_atom_inds(sel_AplusB_heavy[[nm]], natom)
    sel_AplusB_bb[[nm]] <- clean_atom_inds(sel_AplusB_bb[[nm]], natom)
    sel_AplusB_ca[[nm]] <- clean_atom_inds(sel_AplusB_ca[[nm]], natom)
    sel_AplusB_all[[nm]] <- clean_atom_inds(sel_AplusB_all[[nm]], natom)
    sel_mem[[nm]] <- clean_atom_inds(sel_mem[[nm]], natom)
  }

  # Backward compatible keys kept (pep_* / lipo_*) plus generic A/B keys
  list(
    align = sel_align,

    # Region A
    pep_heavy = sel_A_heavy,
    pep_ca = sel_A_ca,
    A_backbone = sel_A_bb,
    A_backbone3 = sel_A_bb3,
    A_heavy = sel_A_heavy,
    A_ca = sel_A_ca,
    A_all = sel_A_all,

    # Region B
    lipo_heavy = sel_B_heavy,
    B_heavy = sel_B_heavy,
    B_ca = sel_B_ca,
    B_backbone = sel_B_bb,
    B_all = sel_B_all,

    # Combined
    AplusB_heavy = sel_AplusB_heavy,
    AplusB_backbone = sel_AplusB_bb,
    AplusB_ca = sel_AplusB_ca,
    AplusB_all = sel_AplusB_all,

    # Membrane
    mem = sel_mem,

    pep_resno = pep_resno,
    lipo_resno = lipo_resno
  )
}



# Convert selection atom indices into column indices of a subset xyz matrix

read_reference_coords_file <- function(path, original_name = NULL, natom_expected = NA_integer_) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop("Reference coordinate file not found.")
  }

  ext <- tolower(tools::file_ext(original_name %||% path))
  lines <- readLines(path, warn = FALSE)
  trim_lines <- trimws(lines)

  if (ext %in% c("pdb", "ent")) {
    ref_pdb <- bio3d::read.pdb(path, ATOM.only = TRUE)
    xyz <- as.numeric(ref_pdb$xyz)
    if (!is.na(natom_expected) && length(xyz) != 3L * natom_expected) {
      stop(sprintf("Reference PDB atom count mismatch: expected %d atoms, found %d.",
                   natom_expected, length(xyz) / 3L))
    }
    return(xyz)
  }

  # Amber ASCII restart / inpcrd / rst7.
  if (ext %in% c("rst7", "rst", "inpcrd")) {
    if (length(lines) < 2) stop("Reference coordinate file is too short.")

    hdr <- suppressWarnings(scan(text = lines[2], what = numeric(), quiet = TRUE))
    nat_file <- if (length(hdr) >= 1) suppressWarnings(as.integer(hdr[1])) else NA_integer_
    nat_use <- if (!is.na(natom_expected)) natom_expected else nat_file

    if (is.na(nat_use) || nat_use < 1) {
      stop("Could not determine atom count for reference coordinate file.")
    }

    nums <- suppressWarnings(scan(
      text = paste(lines[-c(1, 2)], collapse = " "),
      what = numeric(),
      quiet = TRUE
    ))

    need <- 3L * nat_use
    if (length(nums) < need) {
      stop(sprintf("Reference coordinate file does not contain enough coordinates: expected %d values, found %d.",
                   need, length(nums)))
    }

    xyz <- as.numeric(nums[seq_len(need)])
    if (!is.na(nat_file) && nat_file != nat_use) {
      warning(sprintf("Reference coordinate file reports %d atoms but the topology expects %d; using topology atom count.",
                      nat_file, nat_use))
    }
    return(xyz)
  }

  # Generic Amber ASCII coordinate-like fallback. Skip Fortran format lines such
  # as '(8F9.3)' or 'FORMAT(...)'.
  keep <- !grepl("^\\s*(\\(|FORMAT\\()", trim_lines, ignore.case = TRUE)
  nums <- suppressWarnings(scan(
    text = paste(lines[keep], collapse = " "),
    what = numeric(),
    quiet = TRUE
  ))

  nat_use <- natom_expected
  if (is.na(nat_use) || nat_use < 1) {
    stop("Could not determine atom count for generic coordinate file.")
  }

  need <- 3L * nat_use
  if (length(nums) < need) {
    stop(sprintf("Reference coordinate file does not contain enough coordinates: expected %d values, found %d.",
                 need, length(nums)))
  }

  as.numeric(nums[seq_len(need)])
}

subset_xyz_vector <- function(xyz_full, atom_inds) {
  atom_inds <- clean_atom_inds(atom_inds, natom = length(xyz_full) / 3L)
  if (length(atom_inds) < 1) stop("subset_xyz_vector: atom selection is empty.")
  cols <- sort(c(3L * atom_inds - 2L, 3L * atom_inds - 1L, 3L * atom_inds))
  as.numeric(xyz_full[cols])
}

cols_for_atoms <- function(sel_atom, subset_atom) {
  idx <- match(sel_atom, subset_atom)
  if (any(is.na(idx))) stop("Selection atoms are not fully contained in subset atoms (check selections).")
  as.integer(unlist(lapply(idx, function(i) 3L*(i-1L) + 1:3)))
}

# Geometric center from xyz matrix columns (frames x 3M)
center_from_xyz_cols <- function(xyz_mat, cols) {
  sub <- xyz_mat[, cols, drop=FALSE]
  ncol3 <- ncol(sub)
  if (ncol3 %% 3 != 0) stop("Center selection columns not multiple of 3")
  x <- sub[, seq(1, ncol3, by=3), drop=FALSE]
  y <- sub[, seq(2, ncol3, by=3), drop=FALSE]
  z <- sub[, seq(3, ncol3, by=3), drop=FALSE]
  cbind(rowMeans(x), rowMeans(y), rowMeans(z))
}

# Get number of frames in an AMBER NetCDF segment (uses ncdf4 metadata)
get_nc_nframes <- function(nc_path) {
  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  # Prefer explicit dimension if present
  if ("frame" %in% names(nc$dim)) return(nc$dim[["frame"]]$len)
  if ("time"  %in% names(nc$dim)) return(nc$dim[["time"]]$len)

  # If coordinates var exists, find the dim named 'frame' or 'time'
  if ("coordinates" %in% names(nc$var)) {
    v <- nc$var[["coordinates"]]
    dim_names <- vapply(v$dim, function(d) d$name, character(1))
    i <- match("frame", dim_names)
    if (!is.na(i)) return(v$varsize[i])
    i <- match("time", dim_names)
    if (!is.na(i)) return(v$varsize[i])
  }

  # Fallback: choose the largest dimension excluding known non-time dims
  lens <- vapply(nc$dim, function(d) d$len, numeric(1))
  nms <- names(lens)
  keep <- !(nms %in% c("atom", "spatial", "xyz", "x", "y", "z"))
  lens2 <- lens[keep]
  lens2 <- lens2[lens2 > 3]
  if (length(lens2) > 0) return(max(lens2))

  stop("Could not determine number of frames from NetCDF: ", nc_path)
}

# ─────────────────────────────────────────────────────────────────────────────
# Core analysis
# ─────────────────────────────────────────────────────────────────────────────
run_basic_analysis <- function(project_dir,
                               trajs_ordered,
                               combine = TRUE,
                               out_dir = project_results_dir(project_dir),
                               dt_ps = 10,
                               first = NULL, last = NULL, stride = 1,
                               pep_resno = integer(),
                               lipo_resno = integer(),
                               membrane_resid = character(),
                               exclude_resid = c("HOH","WAT","Na+","Cl-"),
                               mem_target = c("B","A"),
                               align_elety = c("N","CA","C","O"),
                               rmsd_ref_mode = c("global_first", "segment_first"),
                               rmsd_ref_path = NULL,
                               rmsd_ref_name = NULL,
                               progress = NULL,
                               logf = NULL,
                               verbose = TRUE) {

  dir_create(out_dir)

  inp <- detect_inputs(project_dir)
  if (is.na(inp$prmtop)) stop("No topology file (.prmtop or .parm7) found in selected folder.")
  if (length(trajs_ordered) < 1) stop("No .nc trajectories selected.")

  prm <- read.prmtop(inp$prmtop)
  mem_target <- match.arg(mem_target)
  rmsd_ref_mode <- match.arg(rmsd_ref_mode)
  sels <- make_selections(prm, pep_resno, lipo_resno, membrane_resid,
                         exclude_resid = exclude_resid,
                         align_elety = align_elety)

  if (length(sels$align$atom) < 3) stop("Alignment selection too small (check Selection A residue numbers).")
  if (length(sels$pep_heavy$atom) < 1) stop("Selection A selection empty (check Selection A residue numbers).")

  has_B   <- length(lipo_resno) > 0
  has_mem <- length(membrane_resid) > 0

  if (has_B && length(sels$lipo_heavy$atom) < 1) stop("Selection B selection empty (check Selection B residue numbers).")
  if (has_mem && length(sels$mem$atom) < 1) stop("Membrane selection empty (check membrane residue names).")

  # IMPORTANT: read.ncdf(at.sel=...) expects ATOM indices as select object
  subset_atom <- sort(unique(c(
    sels$align$atom,
    sels$pep_heavy$atom,
    sels$pep_ca$atom,
    sels$A_backbone$atom,
    sels$A_all$atom,
    sels$lipo_heavy$atom,
    sels$B_backbone$atom,
    sels$B_all$atom,
    sels$mem$atom
  )))
  subset_atom <- clean_atom_inds(subset_atom, natom = get_prm_natom(prm))
  if (length(subset_atom) < 1) stop("Atom selection is empty after cleaning. Check residues / exclude list.")
  at_sel <- as.select(subset_atom)

  # Global reference for alignment / RMSD / PCA. By default this is the first
  # frame of the first selected segment, but the user can override it with an
  # external coordinate file (PDB / rst7 / inpcrd).
  if (!is.null(rmsd_ref_path) && nzchar(rmsd_ref_path)) {
    if (verbose) message("Reading external reference coordinates from: ", rmsd_ref_path)
    ref_full <- read_reference_coords_file(rmsd_ref_path, original_name = rmsd_ref_name, natom_expected = get_prm_natom(prm))
    ref_vec_global <- subset_xyz_vector(ref_full, subset_atom)
  } else {
    ref_seg <- trajs_ordered[[1]]
    if (verbose) message("Reading global reference frame from: ", ref_seg)
    ref_xyz <- read.ncdf(ref_seg, first = 1, last = 1, stride = 1, at.sel = at_sel, verbose = verbose)
    ref_vec_global <- ref_xyz[1, ]
  }

  cols_align   <- cols_for_atoms(sels$align$atom, subset_atom)
  cols_pepH    <- cols_for_atoms(sels$pep_heavy$atom, subset_atom)
  cols_pepCA   <- cols_for_atoms(sels$pep_ca$atom, subset_atom)
  cols_pepBB   <- cols_for_atoms(sels$A_backbone$atom, subset_atom)
  cols_pepAll  <- cols_for_atoms(sels$A_all$atom, subset_atom)
  cols_lipoH   <- if (length(sels$lipo_heavy$atom) > 0) cols_for_atoms(sels$lipo_heavy$atom, subset_atom) else NULL
  cols_lipoBB  <- if (length(sels$B_backbone$atom) > 0) cols_for_atoms(sels$B_backbone$atom, subset_atom) else NULL
  cols_lipoAll <- if (length(sels$B_all$atom) > 0) cols_for_atoms(sels$B_all$atom, subset_atom) else NULL
  cols_mem     <- if (length(sels$mem$atom) > 0) cols_for_atoms(sels$mem$atom, subset_atom) else NULL

  # RMSF grouping by residue for CA atoms (Selection A and B)
  pdb_like <- try(as.pdb(prm), silent = TRUE)
  if (!inherits(pdb_like,"try-error") && !is.null(pdb_like$atom$resno)) {
    rmsf_resno   <- pdb_like$atom$resno[sels$pep_ca$atom]
    rmsf_resno_B <- if (length(sels$B_ca$atom) > 0) pdb_like$atom$resno[sels$B_ca$atom] else integer()
  } else {
    rmsf_resno   <- pep_resno[seq_len(length(sels$pep_ca$atom))]
    rmsf_resno_B <- if (length(sels$B_ca$atom) > 0) lipo_resno[seq_len(length(sels$B_ca$atom))] else integer()
  }

  # Pre-compute segment offsets in frames (true counts), for global time if combine=TRUE
  seg_offsets <- NULL
  if (isTRUE(combine)) {
    nframes <- vapply(trajs_ordered, get_nc_nframes, numeric(1))
    seg_offsets <- c(0, cumsum(nframes))[seq_along(trajs_ordered)]  # 0-based frame offset
    names(seg_offsets) <- basename(trajs_ordered)
  }

  # Storage
  rmsd_pep_all <- list(); rmsd_lipo_all <- list()
  rg_pep_all <- list(); rg_lipo_all <- list()
  dist_all <- list(); z_all <- list()
  ca_accum <- NULL
  ca_B_accum <- NULL

  start_frame_num <- if (is.null(first)) 1L else as.integer(first)

  n_seg <- length(trajs_ordered)

  for (i_seg in seq_along(trajs_ordered)) {
    seg <- trajs_ordered[[i_seg]]
    seg_name <- basename(seg)
    if (!is.null(logf)) logf(sprintf("Segment %d/%d: %s", i_seg, n_seg, seg_name))

    # Progress sub-steps per segment: reading 40%, aligning 20%, metrics 40%
    seg_base <- (i_seg - 1) / max(1, n_seg)
    seg_span <- 1 / max(1, n_seg)
    if (!is.null(progress)) progress(value = seg_base, detail = sprintf("Reading %s (%d/%d)…", seg_name, i_seg, n_seg))
    if (verbose) message("Reading segment: ", seg_name)

    xyz <- read.ncdf(seg, first = first, last = last, stride = stride, at.sel = at_sel, verbose = verbose)
    if (!is.matrix(xyz) || nrow(xyz) < 1) next

    if (!is.null(progress)) progress(value = seg_base + seg_span * 0.40, detail = sprintf("Aligning %s (%d frames)…", seg_name, nrow(xyz)))
    if (!is.null(logf)) logf(sprintf("  Read %d frames, aligning…", nrow(xyz)))

    # Reference used for this segment. External reference overrides segment/global modes.
    ref_vec_seg <- if (!is.null(rmsd_ref_path) && nzchar(rmsd_ref_path)) {
      ref_vec_global
    } else if (identical(rmsd_ref_mode, "segment_first")) {
      as.numeric(xyz[1, ])
    } else {
      ref_vec_global
    }

    xyz_fit <- kabsch_fit_xyz(fixed = ref_vec_seg, mobile = xyz,
                              fixed.inds = cols_align, mobile.inds = cols_align)

    if (!is.null(progress)) progress(value = seg_base + seg_span * 0.60, detail = sprintf("Computing metrics %s (%d/%d)…", seg_name, i_seg, n_seg))

    nF <- nrow(xyz_fit)

    # Local frame numbers (within each file), in ORIGINAL frame units (1-based)
    frame_num <- start_frame_num + (0:(nF-1)) * stride
    time_ps <- (frame_num - 1) * dt_ps
    time_ns <- time_ps / 1000

    # Global time (absolute), based on true segment offsets
    if (isTRUE(combine)) {
      off0 <- seg_offsets[[seg_name]]  # 0-based frames before this segment
      frame_global <- off0 + frame_num
      time_ps_global <- (frame_global - 1) * dt_ps
      time_ns_global <- time_ps_global / 1000
    } else {
      frame_global <- rep(NA_integer_, nF)
      time_ps_global <- rep(NA_real_, nF)
      time_ns_global <- rep(NA_real_, nF)
    }

    # RMSD series computed cpptraj-style: each series fits AND calculates on its
    # own atom set using aligned_rmsd_series(). "Heavy" means non-hydrogen atoms.
    rmsd_pep_backbone <- if (length(cols_pepBB) >= 3) {
      aligned_rmsd_series(xyz, ref_vec_seg, fit_cols = cols_pepBB)
    } else rep(NA_real_, nF)

    rmsd_pep_heavy <- if (length(cols_pepH) >= 3) {
      aligned_rmsd_series(xyz, ref_vec_seg, fit_cols = cols_pepH)
    } else rep(NA_real_, nF)

    rmsd_lipo_backbone <- if (!is.null(cols_lipoBB) && length(cols_lipoBB) >= 3) {
      aligned_rmsd_series(xyz, ref_vec_seg, fit_cols = cols_lipoBB)
    } else rep(NA_real_, nF)

    rmsd_lipo_heavy <- if (!is.null(cols_lipoH) && length(cols_lipoH) >= 3) {
      aligned_rmsd_series(xyz, ref_vec_seg, fit_cols = cols_lipoH)
    } else rep(NA_real_, nF)

    rmsd_pep_all[[seg_name]] <- data.frame(
      segment=seg_name, frame=frame_num, time_ps=time_ps, time_ns=time_ns,
      frame_global=frame_global, time_ps_global=time_ps_global, time_ns_global=time_ns_global,
      rmsd_backbone=rmsd_pep_backbone,
      rmsd_heavy=rmsd_pep_heavy,
      rmsd_total=rmsd_pep_heavy,
      rmsd_A=rmsd_pep_heavy
    )
    rmsd_lipo_all[[seg_name]] <- data.frame(
      segment=seg_name, frame=frame_num, time_ps=time_ps, time_ns=time_ns,
      frame_global=frame_global, time_ps_global=time_ps_global, time_ns_global=time_ns_global,
      rmsd_backbone=rmsd_lipo_backbone,
      rmsd_heavy=rmsd_lipo_heavy,
      rmsd_total=rmsd_lipo_heavy,
      rmsd_A=rmsd_lipo_heavy
    )

    rg_pep <- rgyr(xyz_fit[, cols_pepH, drop=FALSE])
    rg_lipo <- if (!is.null(cols_lipoH) && length(cols_lipoH) >= 3) rgyr(xyz_fit[, cols_lipoH, drop = FALSE]) else rep(NA_real_, nF)

    rg_pep_all[[seg_name]] <- data.frame(
      segment=seg_name, frame=frame_num, time_ps=time_ps, time_ns=time_ns,
      frame_global=frame_global, time_ps_global=time_ps_global, time_ns_global=time_ns_global,
      rg_A=rg_pep
    )
    rg_lipo_all[[seg_name]] <- data.frame(
      segment=seg_name, frame=frame_num, time_ps=time_ps, time_ns=time_ns,
      frame_global=frame_global, time_ps_global=time_ps_global, time_ns_global=time_ns_global,
      rg_A=rg_lipo
    )

    # Membrane metrics: geometric center distance + z proxy (optional)
    cols_target <- if (mem_target == "B" && !is.null(cols_lipoH) && length(cols_lipoH) >= 3) cols_lipoH else cols_pepH

    if (!is.null(cols_mem) && length(cols_mem) >= 3) {
      cen_tar <- center_from_xyz_cols(xyz_fit, cols_target)
      cen_mem <- center_from_xyz_cols(xyz_fit, cols_mem)
      d <- sqrt(rowSums((cen_tar - cen_mem)^2))
      zins <- cen_tar[,3] - cen_mem[,3]
    } else {
      d <- rep(NA_real_, nF)
      zins <- rep(NA_real_, nF)
    }

    dist_all[[seg_name]] <- data.frame(
      segment=seg_name, frame=frame_num, time_ps=time_ps, time_ns=time_ns,
      frame_global=frame_global, time_ps_global=time_ps_global, time_ns_global=time_ns_global,
      dist_A=d
    )
    z_all[[seg_name]] <- data.frame(
      segment=seg_name, frame=frame_num, time_ps=time_ps, time_ns=time_ns,
      frame_global=frame_global, time_ps_global=time_ps_global, time_ns_global=time_ns_global,
      z_A=zins
    )

    # RMSF needs all segments aligned to the SAME global reference regardless of
    # rmsd_ref_mode, so we re-align CA atoms to ref_vec_global before accumulating.
    if (length(cols_pepCA) >= 3) {
      xyz_ca_global <- kabsch_fit_xyz(fixed = ref_vec_global, mobile = xyz,
                                      fixed.inds = cols_align, mobile.inds = cols_align)
      ca_accum <- rbind(ca_accum, xyz_ca_global[, cols_pepCA, drop = FALSE])
      # Selection B CA accumulation (same aligned frame)
      cols_Bca <- if (length(sels$B_ca$atom) > 0) cols_for_atoms(sels$B_ca$atom, subset_atom) else NULL
      if (!is.null(cols_Bca) && length(cols_Bca) >= 3)
        ca_B_accum <- rbind(ca_B_accum, xyz_ca_global[, cols_Bca, drop = FALSE])
    }
  }

  # Write combined CSVs (also write generic filenames for researcher-friendly usage)
  write.csv(do.call(rbind, rmsd_pep_all),  file.path(out_dir,"rmsd_regionA_all.csv"), row.names=FALSE)
  write.csv(do.call(rbind, rg_pep_all),    file.path(out_dir,"rg_regionA_all.csv"), row.names=FALSE)

  # Backward-compatible names (older versions)
  write.csv(do.call(rbind, rmsd_pep_all),  file.path(out_dir,"rmsd_peptide_all.csv"), row.names=FALSE)
  write.csv(do.call(rbind, rg_pep_all),    file.path(out_dir,"rg_peptide_all.csv"), row.names=FALSE)

  if (length(sels$lipo_heavy$atom) > 0) {
    write.csv(do.call(rbind, rmsd_lipo_all), file.path(out_dir,"rmsd_regionB_all.csv"), row.names=FALSE)
    write.csv(do.call(rbind, rg_lipo_all),   file.path(out_dir,"rg_regionB_all.csv"), row.names=FALSE)

    # Backward-compatible
    write.csv(do.call(rbind, rmsd_lipo_all), file.path(out_dir,"rmsd_lipopeptide_all.csv"), row.names=FALSE)
    write.csv(do.call(rbind, rg_lipo_all),   file.path(out_dir,"rg_lipopeptide_all.csv"), row.names=FALSE)
  }

  if (length(sels$mem$atom) > 0) {
    write.csv(do.call(rbind, dist_all), file.path(out_dir,"region_mem_com_dist_all.csv"), row.names=FALSE)
    write.csv(do.call(rbind, z_all),    file.path(out_dir,"region_mem_z_all.csv"), row.names=FALSE)

    # Backward-compatible
    write.csv(do.call(rbind, dist_all), file.path(out_dir,"pep_mem_com_dist_all.csv"), row.names=FALSE)
    write.csv(do.call(rbind, z_all),    file.path(out_dir,"pep_mem_z_all.csv"), row.names=FALSE)
  }

  # RMSF per residue (CA only)

  if (!is.null(ca_accum) && nrow(ca_accum) >= 2) {
    r <- rmsf(ca_accum, grpby = rmsf_resno)
    rmsf_df <- data.frame(resno = unique(rmsf_resno), rmsf_A = as.numeric(r))
    write.csv(rmsf_df, file.path(out_dir,"rmsf_regionA_byres.csv"), row.names=FALSE)
    # Backward-compatible
    write.csv(rmsf_df, file.path(out_dir,"rmsf_peptide_byres.csv"), row.names=FALSE)
  }
  # RMSF Selection B
  if (!is.null(ca_B_accum) && nrow(ca_B_accum) >= 2 && length(rmsf_resno_B) > 0) {
    r_B <- rmsf(ca_B_accum, grpby = rmsf_resno_B)
    rmsf_B_df <- data.frame(resno = unique(rmsf_resno_B), rmsf_B = as.numeric(r_B))
    write.csv(rmsf_B_df, file.path(out_dir,"rmsf_regionB_byres.csv"), row.names=FALSE)
  }

  # PCA / dimensionality reduction on aligned Selection A Cα coordinates
  if (!is.null(ca_accum) && nrow(ca_accum) >= 3 && ncol(ca_accum) >= 3) {
    pca_fit <- prcomp(ca_accum, center = TRUE, scale. = FALSE)
    keep_pc <- min(3, ncol(pca_fit$x))
    pca_scores <- as.data.frame(pca_fit$x[, seq_len(keep_pc), drop = FALSE])
    names(pca_scores) <- paste0("PC", seq_len(keep_pc))

    frame_info <- do.call(rbind, rmsd_pep_all)[, c("segment","frame","time_ps","time_ns","frame_global","time_ps_global","time_ns_global"), drop = FALSE]
    frame_info <- frame_info[seq_len(nrow(pca_scores)), , drop = FALSE]
    pca_scores <- cbind(frame_info, pca_scores)
    write.csv(pca_scores, file.path(out_dir, "pca_regionA_scores.csv"), row.names = FALSE)

    var_expl <- (pca_fit$sdev ^ 2)
    pca_var <- data.frame(
      PC = paste0("PC", seq_along(var_expl)),
      variance = var_expl,
      prop_var = var_expl / sum(var_expl),
      cum_var = cumsum(var_expl / sum(var_expl))
    )
    write.csv(pca_var, file.path(out_dir, "pca_regionA_variance.csv"), row.names = FALSE)
  }

  # Compact run summary (helps reproducible reporting)
  safe_mean <- function(x) if (length(x) < 1 || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  safe_sd <- function(x) if (length(x) < 2 || all(is.na(x))) NA_real_ else sd(x, na.rm = TRUE)

  rmsdA_all <- do.call(rbind, rmsd_pep_all)$rmsd_A
  rgA_all <- do.call(rbind, rg_pep_all)$rg_A

  summary_rows <- data.frame(
    metric = c("rmsd_regionA_mean", "rmsd_regionA_sd", "rg_regionA_mean", "rg_regionA_sd"),
    value = c(safe_mean(rmsdA_all), safe_sd(rmsdA_all), safe_mean(rgA_all), safe_sd(rgA_all)),
    stringsAsFactors = FALSE
  )

  if (length(sels$lipo_heavy$atom) > 0) {
    rmsdB_all <- do.call(rbind, rmsd_lipo_all)$rmsd_A
    rgB_all <- do.call(rbind, rg_lipo_all)$rg_A
    summary_rows <- rbind(
      summary_rows,
      data.frame(
        metric = c("rmsd_regionB_mean", "rmsd_regionB_sd", "rg_regionB_mean", "rg_regionB_sd"),
        value = c(safe_mean(rmsdB_all), safe_sd(rmsdB_all), safe_mean(rgB_all), safe_sd(rgB_all)),
        stringsAsFactors = FALSE
      )
    )
  }

  if (length(sels$mem$atom) > 0) {
    distA_all <- do.call(rbind, dist_all)$dist_A
    zA_all <- do.call(rbind, z_all)$z_A
    summary_rows <- rbind(
      summary_rows,
      data.frame(
        metric = c("membrane_com_distance_mean", "membrane_com_distance_sd", "membrane_dz_mean", "membrane_dz_sd"),
        value = c(safe_mean(distA_all), safe_sd(distA_all), safe_mean(zA_all), safe_sd(zA_all)),
        stringsAsFactors = FALSE
      )
    )
  }

  write.csv(summary_rows, file.path(out_dir, "analysis_summary.csv"), row.names = FALSE)

  invisible(out_dir)
}

# ─────────────────────────────────────────────────────────────────────────────
# RMSD structural clustering (reference-aligned, then RMSD-based distances)
# ─────────────────────────────────────────────────────────────────────────────

# Get element (i,j) from an object of class "dist" (1-based indices)
dist_get <- function(d, i, j) {
  if (i == j) return(0)
  if (i > j) { tmp <- i; i <- j; j <- tmp }
  # dist stores lower triangle by columns
  idx <- (j - 1) * (j - 2) / 2 + i
  d[[as.integer(idx)]]
}

# Medoid index within full distance matrix, plus within-cluster avg pairwise distance
medoid_from_dist <- function(d, idxs) {
  m <- length(idxs)
  if (m == 0) return(list(medoid = NA_integer_, avgpair = NA_real_))
  if (m == 1) return(list(medoid = idxs[1], avgpair = 0))

  sums <- numeric(m)
  pair_sum <- 0
  pair_n <- 0

  for (a in 1:(m - 1)) {
    ia <- idxs[a]
    for (b in (a + 1):m) {
      ib <- idxs[b]
      val <- dist_get(d, ia, ib)
      sums[a] <- sums[a] + val
      sums[b] <- sums[b] + val
      pair_sum <- pair_sum + val
      pair_n <- pair_n + 1
    }
  }

  list(
    medoid = idxs[which.min(sums)],
    avgpair = pair_sum / pair_n
  )
}

rmsd_vec <- function(a, b) {
  sqrt(mean((a - b)^2))
}

# Build a frame index table for the subset that will be read (honors first/last/stride)
build_frame_index <- function(trajs_ordered, dt_ps, first = NULL, last = NULL, stride = 1, combine = TRUE) {
  n_total <- vapply(trajs_ordered, get_nc_nframes, numeric(1))
  seg_offsets <- NULL
  if (isTRUE(combine)) {
    seg_offsets <- c(0, cumsum(n_total))[seq_along(trajs_ordered)]  # 0-based
    names(seg_offsets) <- basename(trajs_ordered)
  }

  out <- list()
  for (i in seq_along(trajs_ordered)) {
    seg <- trajs_ordered[[i]]
    seg_name <- basename(seg)
    nF <- n_total[[i]]

    f1 <- if (is.null(first)) 1L else as.integer(first)
    fL <- if (is.null(last))  nF else as.integer(last)
    f1 <- max(1L, min(nF, f1))
    fL <- max(1L, min(nF, fL))
    if (f1 > fL) next

    frames_vec <- seq(f1, fL, by = as.integer(stride))
    time_ns <- ((frames_vec - 1) * dt_ps) / 1000

    if (isTRUE(combine)) {
      off0 <- seg_offsets[[seg_name]]
      frame_global <- off0 + frames_vec
      time_ns_global <- ((frame_global - 1) * dt_ps) / 1000
    } else {
      frame_global <- rep(NA_integer_, length(frames_vec))
      time_ns_global <- rep(NA_real_, length(frames_vec))
    }

    out[[seg_name]] <- data.frame(
      segment = seg_name,
      frame = frames_vec,
      time_ns = time_ns,
      frame_global = frame_global,
      time_ns_global = time_ns_global,
      stringsAsFactors = FALSE
    )
  }

  df <- do.call(rbind, out)
  if (is.null(df) || nrow(df) < 1) stop("No frames found in the selected window (first/last/stride).")
  df$rowid <- seq_len(nrow(df))
  df
}



# =============================================================================
# Advanced membrane analysis helpers (NetCDF slicing + optional CPPTRAJ)
# =============================================================================

nc_guess_varname <- function(nc, candidates) {
  for (nm in candidates) {
    if (nm %in% names(nc$var)) return(nm)
  }
  # fallback: first var containing a candidate substring
  for (pat in candidates) {
    hit <- names(nc$var)[grepl(pat, names(nc$var), ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[[1]])
  }
  NULL
}

nc_dim_positions <- function(var) {
  dnames <- vapply(var$dim, function(d) d$name, character(1))
  lens   <- vapply(var$dim, function(d) d$len, numeric(1))

  pos_frame <- which(dnames %in% c("frame","time","step"))
  pos_atom  <- which(dnames %in% c("atom","atom_number","atoms","natom"))
  pos_spat  <- which(dnames %in% c("spatial","xyz","coord","coords"))

  if (length(pos_spat) < 1) {
    # Heuristic: a dim of length 3 is usually spatial
    pos_spat <- which(lens == 3)[1]
  }
  list(dnames = dnames, lens = lens, pos_frame = pos_frame[1], pos_atom = pos_atom[1], pos_spat = pos_spat[1])
}

nc_read_coords_atoms <- function(nc_path, atom_inds, frame_inds) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) stop("Package 'ncdf4' is required for NetCDF slicing.")
  if (length(atom_inds) < 1) stop("No atoms selected for NetCDF slicing.")
  if (length(frame_inds) < 1) stop("No frames requested for NetCDF slicing.")

  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  vname <- nc_guess_varname(nc, c("coordinates","coord"))
  if (is.null(vname)) stop("Could not find coordinates variable in NetCDF.")
  var <- nc$var[[vname]]
  dp  <- nc_dim_positions(var)

  if (any(is.na(c(dp$pos_frame, dp$pos_atom, dp$pos_spat)))) {
    stop("Unexpected NetCDF coordinate dimensions (need frame/atom/spatial).")
  }

  frame_min <- min(frame_inds); frame_max <- max(frame_inds)
  atom_min  <- min(atom_inds);  atom_max  <- max(atom_inds)

  start <- rep(1, length(dp$dnames))
  count <- rep(-1, length(dp$dnames))

  start[dp$pos_frame] <- frame_min
  count[dp$pos_frame] <- frame_max - frame_min + 1

  start[dp$pos_atom] <- atom_min
  count[dp$pos_atom] <- atom_max - atom_min + 1

  start[dp$pos_spat] <- 1
  count[dp$pos_spat] <- 3

  arr <- ncdf4::ncvar_get(nc, vname, start = start, count = count)

  # Reorder -> (frame, atom, spatial)
  perm <- c(dp$pos_frame, dp$pos_atom, dp$pos_spat)
  arr <- aperm(arr, perm)

  f_rel <- frame_inds - frame_min + 1
  a_rel <- atom_inds - atom_min + 1

  arr[f_rel, a_rel, , drop = FALSE]
}

nc_read_cell_lengths <- function(nc_path, frame_inds) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) stop("Package 'ncdf4' is required for NetCDF slicing.")
  if (length(frame_inds) < 1) stop("No frames requested for NetCDF slicing.")

  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  vname <- nc_guess_varname(nc, c("cell_lengths","cell_length","cell"))
  if (is.null(vname)) stop("Could not find cell_lengths variable in NetCDF.")
  var <- nc$var[[vname]]
  dp  <- nc_dim_positions(var)

  # Here we expect (frame, spatial=3) or (spatial=3, frame)
  pos_frame <- dp$pos_frame
  pos_spat  <- dp$pos_spat
  if (is.na(pos_spat)) {
    # fallback: dim length 3
    pos_spat <- which(dp$lens == 3)[1]
  }
  if (any(is.na(c(pos_frame, pos_spat)))) stop("Unexpected NetCDF cell_lengths dimensions.")

  frame_min <- min(frame_inds); frame_max <- max(frame_inds)
  start <- rep(1, length(dp$dnames))
  count <- rep(-1, length(dp$dnames))

  start[pos_frame] <- frame_min
  count[pos_frame] <- frame_max - frame_min + 1
  start[pos_spat] <- 1
  count[pos_spat] <- 3

  arr <- ncdf4::ncvar_get(nc, vname, start = start, count = count)

  # Reorder -> (frame, spatial)
  perm <- c(pos_frame, pos_spat)
  arr <- aperm(arr, perm)
  f_rel <- frame_inds - frame_min + 1
  arr[f_rel, , drop = FALSE]
}

detect_total_lipids <- function(prm, membrane_resid) {
  if (length(membrane_resid) < 1) return(NA_integer_)
  sel <- atom.select(prm, resid = membrane_resid)
  if (length(sel$atom) < 1) return(NA_integer_)
  if (is.null(prm$atom) || !("resno" %in% names(prm$atom))) return(NA_integer_)
  length(unique(prm$atom$resno[sel$atom]))
}

compute_bilayer_thickness <- function(trajs_ordered, prm, membrane_resid,
                                      head_atoms = "P",
                                      dt_ps, first = NULL, last = NULL, stride = 1, combine = TRUE) {
  if (length(membrane_resid) < 1) stop("Membrane residue names are empty. Set them in Project tab.")
  head_atoms <- parse_csv_tokens(head_atoms)
  if (length(head_atoms) < 1) head_atoms <- "P"

  sel_head <- atom.select(prm, resid = membrane_resid, elety = head_atoms)
  if (length(sel_head$atom) < 2) stop("No headgroup atoms found. Check membrane residue names and headgroup atom name(s).")

  plan <- build_frame_index(trajs_ordered, dt_ps = dt_ps, first = first, last = last, stride = stride, combine = combine)
  out <- list()

  for (seg in unique(plan$segment)) {
    seg_path <- trajs_ordered[[which(basename(trajs_ordered) == seg)[1]]]
    seg_rows <- which(plan$segment == seg)
    seg_frames <- plan$frame[seg_rows]

    arr <- nc_read_coords_atoms(seg_path, atom_inds = sel_head$atom, frame_inds = seg_frames)
    z <- arr[,,3, drop=TRUE]  # (nframes x natoms)
    if (is.null(dim(z))) z <- matrix(z, ncol = 1)

    z_center <- rowMeans(z)
    upper <- z; upper[upper <= z_center] <- NA_real_
    lower <- z; lower[lower >= z_center] <- NA_real_

    z_up <- rowMeans(upper, na.rm = TRUE)
    z_lo <- rowMeans(lower, na.rm = TRUE)

    thick <- z_up - z_lo

    df <- plan[seg_rows, , drop = FALSE]
    df$thickness_A <- thick
    out[[seg]] <- df
  }
  do.call(rbind, out)
}

compute_area_per_lipid <- function(trajs_ordered, prm, membrane_resid,
                                   n_lip_leaflet = NA_integer_,
                                   dt_ps, first = NULL, last = NULL, stride = 1, combine = TRUE) {
  if (length(membrane_resid) < 1) stop("Membrane residue names are empty. Set them in Project tab.")
  n_total <- detect_total_lipids(prm, membrane_resid)
  if (is.na(n_lip_leaflet)) {
    if (!is.na(n_total)) n_lip_leaflet <- as.integer(round(n_total / 2))
  }
  if (is.na(n_lip_leaflet) || n_lip_leaflet < 1) stop("Could not determine lipids per leaflet. Please enter it manually.")

  plan <- build_frame_index(trajs_ordered, dt_ps = dt_ps, first = first, last = last, stride = stride, combine = combine)
  out <- list()

  for (seg in unique(plan$segment)) {
    seg_path <- trajs_ordered[[which(basename(trajs_ordered) == seg)[1]]]
    seg_rows <- which(plan$segment == seg)
    seg_frames <- plan$frame[seg_rows]

    cell <- nc_read_cell_lengths(seg_path, frame_inds = seg_frames)
    Lx <- cell[,1]; Ly <- cell[,2]
    apl <- (Lx * Ly) / as.numeric(n_lip_leaflet)

    df <- plan[seg_rows, , drop = FALSE]
    df$APL_A2 <- apl
    out[[seg]] <- df
  }

  list(df = do.call(rbind, out), n_total = n_total, n_leaflet = n_lip_leaflet)
}

compute_lipid_enrichment <- function(trajs_ordered, prm, sels, membrane_resid,
                                     target = c("A","B"),
                                     cutoff_A = 8,
                                     head_atom = "P",
                                     by_type = TRUE,
                                     dt_ps, first = NULL, last = NULL, stride = 1, combine = TRUE) {
  target <- match.arg(target)
  if (length(membrane_resid) < 1) stop("Membrane residue names are empty. Set them in Project tab.")
  head_atom <- parse_csv_tokens(head_atom)
  if (length(head_atom) < 1) head_atom <- "P"

  natom <- get_prm_natom(prm)

  sel_head <- atom.select(prm, resid = membrane_resid, elety = head_atom)
  sel_head$atom <- clean_atom_inds(sel_head$atom, natom)
  if (length(sel_head$atom) < 1) stop("No headgroup atoms found. Check membrane residue names and headgroup atom name(s).")

  # Choose target atoms (keep it lightweight)
  target_atoms <- integer()
  if (target == "A") {
    target_atoms <- sels$A_ca$atom
    if (length(target_atoms) < 3) target_atoms <- sels$A_backbone3$atom
    if (length(target_atoms) < 3) target_atoms <- sels$A_heavy$atom
  } else {
    target_atoms <- sels$B_heavy$atom
    if (length(target_atoms) < 3) stop("Selection B is empty; provide Selection B residues in Project tab.")
  }
  target_atoms <- clean_atom_inds(target_atoms, natom)
  if (length(target_atoms) < 1) stop("Target selection is empty; check Selection A/B settings.")

  # Lipid type per headgroup atom (residue name). NAs -> 'UNK'
  lip_type <- rep("LIP", length(sel_head$atom))
  if (!is.null(prm$atom) && ("resid" %in% names(prm$atom))) {
    lip_type <- as.character(prm$atom$resid[sel_head$atom])
  }
  lip_type[is.na(lip_type) | !nzchar(lip_type)] <- "UNK"
  if (!isTRUE(by_type)) lip_type <- rep("All lipids", length(lip_type))

  # Stable level order
  lvls <- sort(unique(lip_type))
  if (length(lvls) < 1) lvls <- "All lipids"

  plan <- build_frame_index(trajs_ordered, dt_ps = dt_ps, first = first, last = last, stride = stride, combine = combine)

  all_seg <- list()
  for (seg in unique(plan$segment)) {
    seg_idx <- which(basename(trajs_ordered) == seg)[1]
    if (is.na(seg_idx)) stop("Could not map segment name to file path: ", seg)
    seg_path <- trajs_ordered[[seg_idx]]

    seg_rows <- which(plan$segment == seg)
    seg_frames <- plan$frame[seg_rows]
    nF <- length(seg_frames)
    if (nF < 1) next

    arr_t <- nc_read_coords_atoms(seg_path, atom_inds = target_atoms, frame_inds = seg_frames)
    arr_h <- nc_read_coords_atoms(seg_path, atom_inds = sel_head$atom, frame_inds = seg_frames)

    # Ensure expected shapes (frames x atoms x 3)
    if (length(dim(arr_t)) != 3L || dim(arr_t)[1] != nF || dim(arr_t)[3] != 3L) {
      stop("Target coordinate array has unexpected shape in segment ", seg, ".")
    }
    if (length(dim(arr_h)) != 3L || dim(arr_h)[1] != nF || dim(arr_h)[3] != 3L) {
      stop("Headgroup coordinate array has unexpected shape in segment ", seg, ".")
    }

    tx <- arr_t[,,1, drop = FALSE]; ty <- arr_t[,,2, drop = FALSE]; tz <- arr_t[,,3, drop = FALSE]
    hx <- arr_h[,,1, drop = FALSE]; hy <- arr_h[,,2, drop = FALSE]; hz <- arr_h[,,3, drop = FALSE]

    # Center-of-geometry of target per frame
    comx <- rowMeans(tx[, , 1, drop = TRUE])
    comy <- rowMeans(ty[, , 1, drop = TRUE])
    comz <- rowMeans(tz[, , 1, drop = TRUE])

    # Distances headgroup->target COM
    dx <- sweep(hx[, , 1, drop = TRUE], 1, comx, "-")
    dy <- sweep(hy[, , 1, drop = TRUE], 1, comy, "-")
    dz <- sweep(hz[, , 1, drop = TRUE], 1, comz, "-")
    dist <- sqrt(dx*dx + dy*dy + dz*dz)  # matrix (nF x nHead)
    if (is.null(dim(dist))) dist <- matrix(dist, nrow = nF)

    within <- dist < as.numeric(cutoff_A)
    if (is.null(dim(within))) within <- matrix(within, nrow = nF)

    # Counts per type
    counts_mat <- matrix(0L, nrow = nF, ncol = length(lvls))
    colnames(counts_mat) <- lvls
    for (j in seq_along(lvls)) {
      cols <- which(lip_type == lvls[[j]])
      if (length(cols) > 0) {
        counts_mat[, j] <- rowSums(within[, cols, drop = FALSE])
      }
    }

    base_df <- plan[seg_rows, , drop = FALSE]
    base_rep <- base_df[rep(seq_len(nF), times = length(lvls)), , drop = FALSE]

    seg_long <- data.frame(
      segment = base_rep$segment,
      frame = base_rep$frame,
      time_ns = base_rep$time_ns,
      time_ps = base_rep$time_ns * 1000,
      frame_global = base_rep$frame_global,
      time_ns_global = base_rep$time_ns_global,
      time_ps_global = base_rep$time_ns_global * 1000,
      lipid_type = rep(lvls, each = nF),
      count = as.integer(as.vector(counts_mat)),
      stringsAsFactors = FALSE
    )

    all_seg[[seg]] <- seg_long
  }

  enrich_df <- do.call(rbind, all_seg)
  if (is.null(enrich_df) || nrow(enrich_df) < 1) stop("No enrichment data were generated (check frame window and selections).")


  sum_df <- aggregate(count ~ lipid_type, data = enrich_df, FUN = function(x) mean(x, na.rm = TRUE))
  list(df = enrich_df, summary = sum_df)
}

get_atom_table <- function(prm) {
  pdb_like <- try(as.pdb(prm), silent = TRUE)
  if (!inherits(pdb_like, "try-error") && !is.null(pdb_like$atom) && is.data.frame(pdb_like$atom) &&
      nrow(pdb_like$atom) > 0) {
    a <- pdb_like$atom
  } else if (!is.null(prm$atom) && is.data.frame(prm$atom) && nrow(prm$atom) > 0) {
    a <- prm$atom
  } else {
    stop("Could not derive atom table from topology object.")
  }
  a$atom_index <- seq_len(nrow(a))
  if (!("resid" %in% names(a))) a$resid <- ""
  if (!("resno" %in% names(a))) a$resno <- seq_len(nrow(a))
  if (!("elety" %in% names(a))) a$elety <- ""
  a
}

sort_chain_atom_df <- function(df) {
  if (is.null(df) || nrow(df) < 2) return(df)
  nm <- trimws(as.character(df$elety))
  digs <- suppressWarnings(as.numeric(gsub("\\D+", "", nm)))
  ord <- order(ifelse(is.na(digs), Inf, digs), nm, seq_len(nrow(df)))
  df[ord, , drop = FALSE]
}

pairwise_contact_frame <- function(Axyz, Bxyz, cutoff_A = 4.5, pair_limit = 5e6, return_hit = FALSE) {
  nA <- nrow(Axyz); nB <- nrow(Bxyz)
  if (nA < 1 || nB < 1) {
    return(list(min_dist = NA_real_, contact_count = 0L,
                A_contact = rep(FALSE, nA), B_contact = rep(FALSE, nB),
                hit_pairs = matrix(integer(0), ncol = 2)))
  }

  cutoff2 <- as.numeric(cutoff_A)^2
  min_d2 <- Inf
  contact_count <- 0L
  A_contact <- rep(FALSE, nA)
  B_contact <- rep(FALSE, nB)
  hit_pairs <- matrix(integer(0), ncol = 2)

  if ((nA * nB) <= pair_limit) {
    dx <- outer(Axyz[,1], Bxyz[,1], "-")
    dy <- outer(Axyz[,2], Bxyz[,2], "-")
    dz <- outer(Axyz[,3], Bxyz[,3], "-")
    d2 <- dx*dx + dy*dy + dz*dz
    min_d2 <- min(d2, na.rm = TRUE)
    hit <- d2 < cutoff2
    contact_count <- sum(hit, na.rm = TRUE)
    A_contact <- rowSums(hit, na.rm = TRUE) > 0
    B_contact <- colSums(hit, na.rm = TRUE) > 0
    if (isTRUE(return_hit) && any(hit, na.rm = TRUE)) {
      hit_pairs <- which(hit, arr.ind = TRUE)
    }
  } else {
    chunk <- max(50L, floor(pair_limit / max(1L, nA)))
    pair_chunks <- list()
    for (j1 in seq(1L, nB, by = chunk)) {
      j2 <- min(nB, j1 + chunk - 1L)
      Bj <- Bxyz[j1:j2, , drop = FALSE]
      dx <- outer(Axyz[,1], Bj[,1], "-")
      dy <- outer(Axyz[,2], Bj[,2], "-")
      dz <- outer(Axyz[,3], Bj[,3], "-")
      d2 <- dx*dx + dy*dy + dz*dz
      min_d2 <- min(min_d2, min(d2, na.rm = TRUE))
      hit <- d2 < cutoff2
      contact_count <- contact_count + sum(hit, na.rm = TRUE)
      A_contact <- A_contact | (rowSums(hit, na.rm = TRUE) > 0)
      B_contact[j1:j2] <- B_contact[j1:j2] | (colSums(hit, na.rm = TRUE) > 0)
      if (isTRUE(return_hit) && any(hit, na.rm = TRUE)) {
        hp <- which(hit, arr.ind = TRUE)
        if (!is.null(hp) && nrow(hp) > 0) {
          hp[,2] <- hp[,2] + j1 - 1L
          pair_chunks[[length(pair_chunks) + 1L]] <- hp
        }
      }
    }
    if (isTRUE(return_hit) && length(pair_chunks) > 0) {
      hit_pairs <- do.call(rbind, pair_chunks)
    }
  }

  list(
    min_dist = if (is.finite(min_d2)) sqrt(min_d2) else NA_real_,
    contact_count = as.integer(contact_count),
    A_contact = A_contact,
    B_contact = B_contact,
    hit_pairs = hit_pairs
  )
}

add_named_counts <- function(counts, keys) {
  if (length(keys) < 1) return(counts)
  keys <- unique(as.character(keys))
  old_names <- names(counts)
  new_keys <- setdiff(keys, old_names)
  if (length(new_keys) > 0) {
    counts <- c(counts, stats::setNames(rep.int(0L, length(new_keys)), new_keys))
  }
  counts[keys] <- counts[keys] + 1L
  counts
}


is_hbond_donor_atom <- function(atom_names) {
  nm <- toupper(trimws(as.character(atom_names)))
  grepl("^(N|O|S)", nm)
}

is_hbond_acceptor_atom <- function(atom_names, include_n = FALSE) {
  nm <- toupper(trimws(as.character(atom_names)))
  pat <- if (isTRUE(include_n)) "^(O|S|N)" else "^(O|S)"
  grepl(pat, nm)
}

compute_hbond_proxy_metrics <- function(trajs_ordered, prm, sels,
                                        cutoff_A = 3.5,
                                        include_n_acceptors = FALSE,
                                        dt_ps, first = NULL, last = NULL, stride = 1, combine = TRUE) {
  atom_df <- get_atom_table(prm)
  natom <- nrow(atom_df)

  A_atoms <- clean_atom_inds(sels[["A_all"]]$atom, natom)
  B_atoms <- clean_atom_inds(sels[["B_all"]]$atom, natom)
  if (length(A_atoms) < 1) stop("Selection A atoms are empty for hydrogen-bond analysis.")
  if (length(B_atoms) < 1) stop("Selection B atoms are empty for hydrogen-bond analysis.")

  atom_names <- trimws(as.character(atom_df$elety))
  A_donor <- A_atoms[is_hbond_donor_atom(atom_names[A_atoms])]
  A_accept <- A_atoms[is_hbond_acceptor_atom(atom_names[A_atoms], include_n = include_n_acceptors)]
  B_donor <- B_atoms[is_hbond_donor_atom(atom_names[B_atoms])]
  B_accept <- B_atoms[is_hbond_acceptor_atom(atom_names[B_atoms], include_n = include_n_acceptors)]

  if (length(A_donor) < 1 && length(B_donor) < 1) {
    stop("No donor-like atoms (N/O/S) were found in Selection A or B.")
  }
  if (length(A_accept) < 1 && length(B_accept) < 1) {
    stop("No acceptor-like atoms (O/S; optionally N if enabled) were found in Selection A or B.")
  }
  if (length(A_donor) < 1 || length(B_accept) < 1) {
    warning("A→B hydrogen-bond direction has no donor/acceptor candidates with the current rules.")
  }
  if (length(B_donor) < 1 || length(A_accept) < 1) {
    warning("B→A hydrogen-bond direction has no donor/acceptor candidates with the current rules.")
  }

  label_res <- function(idx) paste0(atom_df$resid[idx], ":", atom_df$resno[idx])
  label_atom <- function(idx) paste0(label_res(idx), "@", trimws(as.character(atom_df$elety[idx])))

  A_donor_labels <- if (length(A_donor) > 0) label_atom(A_donor) else character()
  A_accept_labels <- if (length(A_accept) > 0) label_atom(A_accept) else character()
  B_donor_labels <- if (length(B_donor) > 0) label_atom(B_donor) else character()
  B_accept_labels <- if (length(B_accept) > 0) label_atom(B_accept) else character()

  plan <- build_frame_index(trajs_ordered, dt_ps = dt_ps, first = first, last = last,
                            stride = stride, combine = combine)
  atom_union <- sort(unique(c(A_donor, A_accept, B_donor, B_accept)))
  map_A_donor  <- match(A_donor, atom_union)
  map_A_accept <- match(A_accept, atom_union)
  map_B_donor  <- match(B_donor, atom_union)
  map_B_accept <- match(B_accept, atom_union)

  ts_rows <- list()
  pair_counts <- integer(0)

  for (seg in unique(plan$segment)) {
    seg_idx <- which(basename(trajs_ordered) == seg)[1]
    seg_path <- trajs_ordered[[seg_idx]]
    seg_rows <- which(plan$segment == seg)
    seg_frames <- plan$frame[seg_rows]
    nF <- length(seg_frames)
    if (nF < 1) next

    arr <- nc_read_coords_atoms(seg_path, atom_inds = atom_union, frame_inds = seg_frames)
    for (ff in seq_len(nF)) {
      frame_row <- plan[seg_rows[ff], , drop = FALSE]
      pair_keys_frame <- character(0)
      count_AtoB <- 0L
      count_BtoA <- 0L
      mins <- numeric(0)

      if (length(map_A_donor) > 0 && length(map_B_accept) > 0) {
        Axyz <- cbind(arr[ff, map_A_donor, 1], arr[ff, map_A_donor, 2], arr[ff, map_A_donor, 3])
        Bxyz <- cbind(arr[ff, map_B_accept, 1], arr[ff, map_B_accept, 2], arr[ff, map_B_accept, 3])
        if (is.null(dim(Axyz))) Axyz <- matrix(Axyz, ncol = 3)
        if (is.null(dim(Bxyz))) Bxyz <- matrix(Bxyz, ncol = 3)
        info_ab <- pairwise_contact_frame(Axyz, Bxyz, cutoff_A = cutoff_A, return_hit = TRUE)
        count_AtoB <- if (!is.null(info_ab$hit_pairs) && nrow(info_ab$hit_pairs) > 0) nrow(info_ab$hit_pairs) else 0L
        mins <- c(mins, info_ab$min_dist)
        if (!is.null(info_ab$hit_pairs) && nrow(info_ab$hit_pairs) > 0) {
          pair_keys_frame <- c(pair_keys_frame,
                               paste0("AtoB||", A_donor_labels[info_ab$hit_pairs[,1]], "||", B_accept_labels[info_ab$hit_pairs[,2]]))
        }
      }

      if (length(map_B_donor) > 0 && length(map_A_accept) > 0) {
        Bxyz <- cbind(arr[ff, map_B_donor, 1], arr[ff, map_B_donor, 2], arr[ff, map_B_donor, 3])
        Axyz <- cbind(arr[ff, map_A_accept, 1], arr[ff, map_A_accept, 2], arr[ff, map_A_accept, 3])
        if (is.null(dim(Bxyz))) Bxyz <- matrix(Bxyz, ncol = 3)
        if (is.null(dim(Axyz))) Axyz <- matrix(Axyz, ncol = 3)
        info_ba <- pairwise_contact_frame(Bxyz, Axyz, cutoff_A = cutoff_A, return_hit = TRUE)
        count_BtoA <- if (!is.null(info_ba$hit_pairs) && nrow(info_ba$hit_pairs) > 0) nrow(info_ba$hit_pairs) else 0L
        mins <- c(mins, info_ba$min_dist)
        if (!is.null(info_ba$hit_pairs) && nrow(info_ba$hit_pairs) > 0) {
          pair_keys_frame <- c(pair_keys_frame,
                               paste0("BtoA||", B_donor_labels[info_ba$hit_pairs[,1]], "||", A_accept_labels[info_ba$hit_pairs[,2]]))
        }
      }

      pair_keys_frame <- unique(pair_keys_frame)
      if (length(pair_keys_frame) > 0) pair_counts <- add_named_counts(pair_counts, pair_keys_frame)
      min_use <- suppressWarnings(min(mins[is.finite(mins)], na.rm = TRUE))
      if (!is.finite(min_use)) min_use <- NA_real_

      ts_rows[[length(ts_rows) + 1L]] <- data.frame(
        segment = frame_row$segment,
        frame = frame_row$frame,
        time_ns = frame_row$time_ns,
        frame_global = frame_row$frame_global,
        time_ns_global = frame_row$time_ns_global,
        hbond_count = as.integer(length(pair_keys_frame)),
        hbond_count_AtoB = as.integer(count_AtoB),
        hbond_count_BtoA = as.integer(count_BtoA),
        min_da_dist = min_use,
        stringsAsFactors = FALSE
      )
    }
  }

  ts_df <- do.call(rbind, ts_rows)
  if (is.null(ts_df) || nrow(ts_df) < 1) stop("No hydrogen-bond proxy data were generated.")

  pair_df <- NULL
  if (length(pair_counts) > 0) {
    pair_df <- data.frame(
      pair_key = names(pair_counts),
      occupancy_frames = as.integer(pair_counts),
      occupancy_frac = as.numeric(pair_counts) / nrow(ts_df),
      occupancy_pct = 100 * as.numeric(pair_counts) / nrow(ts_df),
      stringsAsFactors = FALSE
    )
    parts <- strsplit(pair_df$pair_key, "||", fixed = TRUE)
    pair_df$direction_code <- vapply(parts, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1))
    pair_df$donor_label <- vapply(parts, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
    pair_df$acceptor_label <- vapply(parts, function(x) if (length(x) >= 3) x[3] else NA_character_, character(1))
    pair_df$direction <- ifelse(pair_df$direction_code == "AtoB", "Selection A donor → Selection B acceptor", "Selection B donor → Selection A acceptor")
    pair_df$donor_side <- ifelse(pair_df$direction_code == "AtoB", "Selection A", "Selection B")
    pair_df$acceptor_side <- ifelse(pair_df$direction_code == "AtoB", "Selection B", "Selection A")
    parse_label <- function(x) {
      parts2 <- strsplit(as.character(x), "@", fixed = TRUE)
      residue <- vapply(parts2, function(z) if (length(z) >= 1) z[1] else NA_character_, character(1))
      atom <- vapply(parts2, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1))
      data.frame(residue = residue, atom = atom, stringsAsFactors = FALSE)
    }
    donor_parsed <- parse_label(pair_df$donor_label)
    accept_parsed <- parse_label(pair_df$acceptor_label)
    pair_df$donor_residue <- donor_parsed$residue
    pair_df$donor_atom <- donor_parsed$atom
    pair_df$acceptor_residue <- accept_parsed$residue
    pair_df$acceptor_atom <- accept_parsed$atom
    pair_df$pair_label <- paste0(pair_df$donor_residue, "@", pair_df$donor_atom, " → ", pair_df$acceptor_residue, "@", pair_df$acceptor_atom)
    pair_df <- pair_df[, c("direction","donor_side","acceptor_side","donor_residue","donor_atom","acceptor_residue","acceptor_atom","pair_label","occupancy_frames","occupancy_frac","occupancy_pct")]
    pair_df <- pair_df[order(-pair_df$occupancy_pct, pair_df$direction, pair_df$pair_label), , drop = FALSE]
    rownames(pair_df) <- NULL
  }

  list(ts = ts_df, pairs = pair_df)
}

compute_density_profiles <- function(trajs_ordered, prm, sels, membrane_resid,
                                     head_atoms = "P",
                                     tail_regex = "^C2[0-9]+$|^C3[0-9]+$",
                                     water_resid = c("HOH","WAT","TIP3","TIP3P","SOL"),
                                     ion_resid = c("Na+","K+","Cl-","CLA","SOD","POT","MG2","CA"),
                                     target = c("none","A","B"),
                                     binwidth_A = 1,
                                     zlim_A = NA_real_,
                                     center_membrane = TRUE,
                                     dt_ps, first = NULL, last = NULL, stride = 1, combine = TRUE) {
  target <- match.arg(target)
  if (length(membrane_resid) < 1) stop("Membrane residue names are empty. Set them in Project tab.")
  atom_df <- get_atom_table(prm)
  natom <- nrow(atom_df)

  wrap_axis_to_box <- function(z, boxL) {
    if (!is.finite(boxL) || boxL <= 0) return(z)
    ((z + boxL/2) %% boxL) - boxL/2
  }

  head_atoms <- parse_csv_tokens(head_atoms)
  if (length(head_atoms) < 1) head_atoms <- "P"
  water_resid <- unique(parse_csv_tokens(paste(water_resid, collapse = ",")))
  ion_resid   <- unique(parse_csv_tokens(paste(ion_resid, collapse = ",")))

  mem_sel <- atom.select(prm, resid = membrane_resid)
  head_sel <- atom.select(prm, resid = membrane_resid, elety = head_atoms)

  atom_names <- trimws(as.character(atom_df$elety))
  resid_names <- as.character(atom_df$resid)

  tail_idx <- which(resid_names %in% membrane_resid & grepl(tail_regex, atom_names, perl = TRUE))
  water_idx <- if (length(water_resid) > 0) which(resid_names %in% water_resid) else integer()
  ion_idx <- if (length(ion_resid) > 0) which(resid_names %in% ion_resid) else integer()

  groups <- list(
    Headgroups = clean_atom_inds(head_sel$atom, natom),
    Tails = clean_atom_inds(tail_idx, natom),
    Water = clean_atom_inds(water_idx, natom),
    Ions = clean_atom_inds(ion_idx, natom)
  )
  if (target == "A") groups$Target <- clean_atom_inds(sels$A_heavy$atom, natom)
  if (target == "B") groups$Target <- clean_atom_inds(sels$B_heavy$atom, natom)

  groups <- groups[vapply(groups, length, integer(1)) > 0]
  if (length(groups) < 1) stop("No atoms were selected for density profiles.")

  center_atoms <- clean_atom_inds(if (length(head_sel$atom) > 0) head_sel$atom else mem_sel$atom, natom)
  if (length(center_atoms) < 1) stop("Could not determine membrane centering atoms.")

  plan <- build_frame_index(trajs_ordered, dt_ps = dt_ps, first = first, last = last,
                            stride = stride, combine = combine)

  seg_profiles <- list()
  for (seg in unique(plan$segment)) {
    seg_idx <- which(basename(trajs_ordered) == seg)[1]
    seg_path <- trajs_ordered[[seg_idx]]
    seg_rows <- which(plan$segment == seg)
    seg_frames <- plan$frame[seg_rows]
    nF <- length(seg_frames)
    if (nF < 1) next

    atom_union <- sort(unique(c(center_atoms, unlist(groups, use.names = FALSE))))
    arr <- nc_read_coords_atoms(seg_path, atom_inds = atom_union, frame_inds = seg_frames)
    z_union <- arr[,,3, drop = TRUE]
    if (is.null(dim(z_union))) z_union <- matrix(z_union, nrow = nF)

    idx_center <- match(center_atoms, atom_union)
    z_center <- rowMeans(z_union[, idx_center, drop = FALSE], na.rm = TRUE)

    cell <- tryCatch(nc_read_cell_lengths(seg_path, frame_inds = seg_frames), error = function(e) NULL)
    area <- if (!is.null(cell) && ncol(cell) >= 2) pmax(cell[,1] * cell[,2], 1e-6) else rep(1, nF)
    lz <- if (!is.null(cell) && ncol(cell) >= 3) cell[,3] else rep(NA_real_, nF)

    zlim_use <- suppressWarnings(as.numeric(zlim_A)[1])
    if (!is.finite(zlim_use) || zlim_use <= 0) {
      zbox <- if (any(is.finite(lz) & lz > 0)) stats::median(lz[is.finite(lz) & lz > 0], na.rm = TRUE) / 2 else NA_real_
      if (!is.finite(zbox) || zbox <= 0) {
        all_z <- c()
        for (nm in names(groups)) {
          idxg <- match(groups[[nm]], atom_union)
          zg <- z_union[, idxg, drop = FALSE]
          if (isTRUE(center_membrane)) zg <- sweep(zg, 1, z_center, "-")
          for (ff in seq_len(nF)) {
            zz <- as.numeric(zg[ff, ])
            if (isTRUE(center_membrane) && is.finite(lz[min(ff, length(lz))])) {
              zz <- wrap_axis_to_box(zz, lz[min(ff, length(lz))])
            }
            all_z <- c(all_z, zz)
          }
        }
        zbox <- suppressWarnings(stats::quantile(abs(all_z[is.finite(all_z)]), probs = 0.995, na.rm = TRUE))
      }
      zlim_use <- max(10, as.numeric(zbox))
    }

    breaks <- seq(-zlim_use, zlim_use, by = as.numeric(binwidth_A))
    if (length(breaks) < 3) breaks <- seq(-zlim_use, zlim_use, length.out = 81)
    if (tail(breaks, 1) < zlim_use) breaks <- c(breaks, zlim_use)
    mids <- (head(breaks, -1) + tail(breaks, -1)) / 2
    xmin <- min(breaks, na.rm = TRUE)
    xmax <- max(breaks, na.rm = TRUE)
    eps <- max(.Machine$double.eps * max(1, abs(c(xmin, xmax))), 1e-8)

    seg_long <- list()
    for (grp in names(groups)) {
      idxg <- match(groups[[grp]], atom_union)
      zg <- z_union[, idxg, drop = FALSE]
      if (isTRUE(center_membrane)) zg <- sweep(zg, 1, z_center, "-")
      dens <- numeric(length(mids))
      for (ff in seq_len(nF)) {
        zz <- as.numeric(zg[ff, ])
        zz <- zz[is.finite(zz)]
        if (length(zz) < 1) next
        if (isTRUE(center_membrane) && is.finite(lz[min(ff, length(lz))])) {
          zz <- wrap_axis_to_box(zz, lz[min(ff, length(lz))])
        }
        zz <- pmin(pmax(zz, xmin + eps), xmax - eps)
        h <- hist(zz, breaks = breaks, plot = FALSE, include.lowest = TRUE)
        dens <- dens + h$counts / (area[min(ff, length(area))] * as.numeric(binwidth_A))
      }
      dens <- dens / nF
      seg_long[[grp]] <- data.frame(
        segment = seg,
        group = grp,
        z_mid = mids,
        density = dens,
        n_frames = nF,
        stringsAsFactors = FALSE
      )
    }
    seg_profiles[[seg]] <- do.call(rbind, seg_long)
  }

  prof <- do.call(rbind, seg_profiles)
  if (is.null(prof) || nrow(prof) < 1) stop("No density profile data were generated.")
  prof$wdens <- prof$density * prof$n_frames
  agg <- aggregate(cbind(wdens, n_frames) ~ group + z_mid, data = prof, FUN = sum)
  agg$density <- agg$wdens / pmax(agg$n_frames, 1)
  agg <- agg[order(agg$group, agg$z_mid), c("group","z_mid","density","n_frames")]
  rownames(agg) <- NULL
  agg
}

compute_lipid_tail_order <- function(trajs_ordered, prm, membrane_resid,
                                     chain1_regex = "^C2[0-9]+$",
                                     chain2_regex = "^C3[0-9]+$",
                                     by_type = TRUE,
                                     topology_mode = c("standard", "split"),
                                     include_resid = NULL,
                                     chain1_resid = NULL,
                                     chain2_resid = NULL,
                                     dt_ps, first = NULL, last = NULL, stride = 1, combine = TRUE) {
  if (length(membrane_resid) < 1) stop("Membrane residue names are empty. Set them in Project tab.")
  topology_mode <- match.arg(topology_mode)
  atom_df <- get_atom_table(prm)
  mem_df <- atom_df[atom_df$resid %in% membrane_resid, , drop = FALSE]
  if (nrow(mem_df) < 1) stop("No membrane atoms matched the provided membrane residue names.")

  include_resid <- parse_csv_tokens(include_resid)
  chain1_resid <- parse_csv_tokens(chain1_resid)
  chain2_resid <- parse_csv_tokens(chain2_resid)

  if (topology_mode == "standard") {
    if (length(include_resid) > 0) {
      mem_df <- mem_df[mem_df$resid %in% include_resid, , drop = FALSE]
      if (nrow(mem_df) < 1) {
        stop("No lipid residues matched the tail-order residue filter. Check the residue names for tail-order analysis.")
      }
    }
  } else {
    if (length(chain1_resid) < 1 || length(chain2_resid) < 1) {
      stop("Split-residue tail order requires chain 1 residue names and chain 2 residue names.")
    }
    keep_resid <- unique(c(include_resid, chain1_resid, chain2_resid))
    if (length(keep_resid) > 0) {
      mem_df <- mem_df[mem_df$resid %in% keep_resid, , drop = FALSE]
    }
    if (!any(mem_df$resid %in% chain1_resid)) {
      stop("No atoms matched the chain 1 residue names in split mode.")
    }
    if (!any(mem_df$resid %in% chain2_resid)) {
      stop("No atoms matched the chain 2 residue names in split mode.")
    }
  }

  # Build consecutive tail-bond segments for one chain definition.
  # In split mode the residue names that feed sn1 and sn2 are different, so we
  # keep the originating residue name with each segment and enforce that mapping
  # again later before plotting/exporting.
  add_chain_pairs <- function(df_subset, regex, chain_label, pair_rows) {
    if (is.null(df_subset) || nrow(df_subset) < 2) return(pair_rows)
    for (rn in sort(unique(df_subset$resno))) {
      df_r <- df_subset[df_subset$resno == rn, , drop = FALSE]
      if (nrow(df_r) < 2) next
      source_resid <- as.character(df_r$resid[1])
      lip_type <- if (isTRUE(by_type)) source_resid else "All lipids"
      ch <- sort_chain_atom_df(df_r[grepl(regex, trimws(as.character(df_r$elety)), perl = TRUE), , drop = FALSE])
      if (nrow(ch) < 2) next
      for (i in seq_len(nrow(ch) - 1L)) {
        pair_rows[[length(pair_rows) + 1L]] <- data.frame(
          atom1 = ch$atom_index[i],
          atom2 = ch$atom_index[i + 1L],
          chain = chain_label,
          segment_index = i,
          segment_label = paste0(ch$elety[i], "-", ch$elety[i + 1L]),
          lipid_type = lip_type,
          source_resid = source_resid,
          stringsAsFactors = FALSE
        )
      }
    }
    pair_rows
  }

  pair_rows <- list()
  if (topology_mode == "standard") {
    for (rn in sort(unique(mem_df$resno))) {
      df_r <- mem_df[mem_df$resno == rn, , drop = FALSE]
      if (nrow(df_r) < 2) next
      lip_type <- if (isTRUE(by_type)) as.character(df_r$resid[1]) else "All lipids"

      ch1 <- sort_chain_atom_df(df_r[grepl(chain1_regex, trimws(as.character(df_r$elety)), perl = TRUE), , drop = FALSE])
      ch2 <- sort_chain_atom_df(df_r[grepl(chain2_regex, trimws(as.character(df_r$elety)), perl = TRUE), , drop = FALSE])

      if (nrow(ch1) >= 2) {
        for (i in seq_len(nrow(ch1) - 1L)) {
          pair_rows[[length(pair_rows) + 1L]] <- data.frame(
            atom1 = ch1$atom_index[i],
            atom2 = ch1$atom_index[i + 1L],
            chain = "sn1",
            segment_index = i,
            segment_label = paste0(ch1$elety[i], "-", ch1$elety[i + 1L]),
            lipid_type = lip_type,
            stringsAsFactors = FALSE
          )
        }
      }
      if (nrow(ch2) >= 2) {
        for (i in seq_len(nrow(ch2) - 1L)) {
          pair_rows[[length(pair_rows) + 1L]] <- data.frame(
            atom1 = ch2$atom_index[i],
            atom2 = ch2$atom_index[i + 1L],
            chain = "sn2",
            segment_index = i,
            segment_label = paste0(ch2$elety[i], "-", ch2$elety[i + 1L]),
            lipid_type = lip_type,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  } else {
    pair_rows <- add_chain_pairs(mem_df[mem_df$resid %in% chain1_resid, , drop = FALSE], chain1_regex, "sn1", pair_rows)
    pair_rows <- add_chain_pairs(mem_df[mem_df$resid %in% chain2_resid, , drop = FALSE], chain2_regex, "sn2", pair_rows)
  }

  pair_df <- if (length(pair_rows) > 0) do.call(rbind, pair_rows) else NULL
  if (is.null(pair_df) || nrow(pair_df) < 1) {
    if (topology_mode == "split") {
      stop("No tail atom pairs were built in split mode. Check the chain residue names and regex patterns for your lipid atom names.")
    } else {
      stop("No tail atom pairs were built. Check the chain regex patterns for your lipid atom names.")
    }
  }

  plan <- build_frame_index(trajs_ordered, dt_ps = dt_ps, first = first, last = last,
                            stride = stride, combine = combine)

  pair_sum <- numeric(nrow(pair_df))
  pair_n <- numeric(nrow(pair_df))

  atom_union <- sort(unique(c(pair_df$atom1, pair_df$atom2)))
  map1 <- match(pair_df$atom1, atom_union)
  map2 <- match(pair_df$atom2, atom_union)

  for (seg in unique(plan$segment)) {
    seg_idx <- which(basename(trajs_ordered) == seg)[1]
    seg_path <- trajs_ordered[[seg_idx]]
    seg_rows <- which(plan$segment == seg)
    seg_frames <- plan$frame[seg_rows]
    nF <- length(seg_frames)
    if (nF < 1) next

    arr <- nc_read_coords_atoms(seg_path, atom_inds = atom_union, frame_inds = seg_frames)
    for (j in seq_len(nrow(pair_df))) {
      vx <- arr[, map2[j], 1] - arr[, map1[j], 1]
      vy <- arr[, map2[j], 2] - arr[, map1[j], 2]
      vz <- arr[, map2[j], 3] - arr[, map1[j], 3]
      vn <- sqrt(vx*vx + vy*vy + vz*vz)
      ok <- is.finite(vn) & (vn > 0)
      if (!any(ok)) next
      s <- 0.5 * (3 * (vz[ok] / vn[ok])^2 - 1)
      pair_sum[j] <- pair_sum[j] + sum(abs(s), na.rm = TRUE)
      pair_n[j] <- pair_n[j] + sum(is.finite(s))
    }
  }

  pair_df$sum_absS <- pair_sum
  pair_df$n_obs <- pair_n
  pair_df <- pair_df[pair_df$n_obs > 0, , drop = FALSE]
  if (nrow(pair_df) < 1) stop("No order-parameter values could be computed from the selected tail atoms.")

  if (topology_mode == "split") {
    # Guard against mixed labels in split topologies. Only keep the residue names
    # explicitly assigned to each chain so the plot cannot show PA-sn2 / OL-sn1
    # style cross-combinations.
    valid_chain1 <- pair_df$chain == "sn1" & pair_df$source_resid %in% chain1_resid
    valid_chain2 <- pair_df$chain == "sn2" & pair_df$source_resid %in% chain2_resid
    pair_df <- pair_df[valid_chain1 | valid_chain2, , drop = FALSE]
    if (nrow(pair_df) < 1) {
      stop("No valid tail-parameter segments remained after enforcing the split-residue chain mapping.")
    }
  }

  agg <- aggregate(cbind(sum_absS, n_obs) ~ lipid_type + chain + segment_index + segment_label,
                   data = pair_df, FUN = sum)
  
  agg$absS <- agg$sum_absS / pmax(agg$n_obs, 1)
  agg$series <- if (isTRUE(by_type)) paste(agg$lipid_type, agg$chain, sep = " — ") else agg$chain
  
  agg <- agg[is.finite(agg$absS) & agg$n_obs > 0, , drop = FALSE]
  if (nrow(agg) < 1) {
    stop("No valid tail-order segments remained after aggregation. Check the residue names and regex patterns.")
  }
  
  agg$topology_mode <- topology_mode
  agg <- agg[order(agg$series, agg$segment_index),
             c("lipid_type","chain","series","topology_mode","segment_index","segment_label","absS","n_obs")]
  rownames(agg) <- NULL
  agg
}


compute_interaction_metrics <- function(trajs_ordered, prm, sels,
                                        selA_key = c("A_heavy","A_backbone","A_ca"),
                                        selB_key = c("B_heavy","B_backbone","B_ca"),
                                        cutoff_A = 4.5,
                                        dt_ps, first = NULL, last = NULL, stride = 1, combine = TRUE) {
  selA_key <- match.arg(selA_key)
  selB_key <- match.arg(selB_key)
  atom_df <- get_atom_table(prm)
  natom <- nrow(atom_df)

  A_atoms <- clean_atom_inds(sels[[selA_key]]$atom, natom)
  B_atoms <- clean_atom_inds(sels[[selB_key]]$atom, natom)
  if (length(A_atoms) < 1) stop("Selection A atoms are empty for the chosen interaction selection.")
  if (length(B_atoms) < 1) stop("Selection B atoms are empty for the chosen interaction selection.")

  A_meta <- atom_df[A_atoms, c("resno","resid"), drop = FALSE]
  B_meta <- atom_df[B_atoms, c("resno","resid"), drop = FALSE]
  A_keys <- paste0(A_meta$resid, ":", A_meta$resno)
  B_keys <- paste0(B_meta$resid, ":", B_meta$resno)
  A_levels <- unique(A_keys)
  B_levels <- unique(B_keys)
  A_counts <- setNames(integer(length(A_levels)), A_levels)
  B_counts <- setNames(integer(length(B_levels)), B_levels)
  pair_counts <- integer(0)

  plan <- build_frame_index(trajs_ordered, dt_ps = dt_ps, first = first, last = last,
                            stride = stride, combine = combine)
  ts_rows <- list()

  atom_union <- sort(unique(c(A_atoms, B_atoms)))
  mapA <- match(A_atoms, atom_union)
  mapB <- match(B_atoms, atom_union)

  for (seg in unique(plan$segment)) {
    seg_idx <- which(basename(trajs_ordered) == seg)[1]
    seg_path <- trajs_ordered[[seg_idx]]
    seg_rows <- which(plan$segment == seg)
    seg_frames <- plan$frame[seg_rows]
    nF <- length(seg_frames)
    if (nF < 1) next

    arr <- nc_read_coords_atoms(seg_path, atom_inds = atom_union, frame_inds = seg_frames)
    for (ff in seq_len(nF)) {
      Axyz <- cbind(arr[ff, mapA, 1], arr[ff, mapA, 2], arr[ff, mapA, 3])
      Bxyz <- cbind(arr[ff, mapB, 1], arr[ff, mapB, 2], arr[ff, mapB, 3])
      if (is.null(dim(Axyz))) Axyz <- matrix(Axyz, ncol = 3)
      if (is.null(dim(Bxyz))) Bxyz <- matrix(Bxyz, ncol = 3)

      comA <- colMeans(Axyz)
      comB <- colMeans(Bxyz)
      com_dist <- sqrt(sum((comA - comB)^2))

      info <- pairwise_contact_frame(Axyz, Bxyz, cutoff_A = cutoff_A, return_hit = TRUE)
      frame_row <- plan[seg_rows[ff], , drop = FALSE]
      ts_rows[[length(ts_rows) + 1L]] <- data.frame(
        segment = frame_row$segment,
        frame = frame_row$frame,
        time_ns = frame_row$time_ns,
        frame_global = frame_row$frame_global,
        time_ns_global = frame_row$time_ns_global,
        min_dist_A = info$min_dist,
        com_dist_A = com_dist,
        atom_contacts = info$contact_count,
        stringsAsFactors = FALSE
      )

      if (any(info$A_contact)) {
        hitA <- unique(A_keys[info$A_contact])
        A_counts[hitA] <- A_counts[hitA] + 1L
      }
      if (any(info$B_contact)) {
        hitB <- unique(B_keys[info$B_contact])
        B_counts[hitB] <- B_counts[hitB] + 1L
      }
      if (!is.null(info$hit_pairs) && length(info$hit_pairs) > 0 && nrow(info$hit_pairs) > 0) {
        pair_keys_frame <- unique(paste0(A_keys[info$hit_pairs[,1]], " || ", B_keys[info$hit_pairs[,2]]))
        pair_counts <- add_named_counts(pair_counts, pair_keys_frame)
      }
    }
  }

  ts_df <- do.call(rbind, ts_rows)
  if (is.null(ts_df) || nrow(ts_df) < 1) stop("No interaction time-series data were generated.")

  total_frames <- nrow(ts_df)
  occA <- data.frame(
    side = "Selection A",
    residue = names(A_counts),
    occupancy_frames = as.integer(A_counts),
    occupancy_frac = as.numeric(A_counts) / total_frames,
    occupancy_pct = 100 * as.numeric(A_counts) / total_frames,
    stringsAsFactors = FALSE
  )
  occB <- data.frame(
    side = "Selection B",
    residue = names(B_counts),
    occupancy_frames = as.integer(B_counts),
    occupancy_frac = as.numeric(B_counts) / total_frames,
    occupancy_pct = 100 * as.numeric(B_counts) / total_frames,
    stringsAsFactors = FALSE
  )

  occ_df <- rbind(occA, occB)
  occ_df <- occ_df[order(-occ_df$occupancy_pct, occ_df$side, occ_df$residue), ]
  rownames(occ_df) <- NULL

  pair_df <- NULL
  if (length(pair_counts) > 0) {
    pair_df <- data.frame(
      pair_key = names(pair_counts),
      occupancy_frames = as.integer(pair_counts),
      occupancy_frac = as.numeric(pair_counts) / total_frames,
      occupancy_pct = 100 * as.numeric(pair_counts) / total_frames,
      stringsAsFactors = FALSE
    )
    parts <- strsplit(pair_df$pair_key, " || ", fixed = TRUE)
    pair_df$residue_A <- vapply(parts, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1))
    pair_df$residue_B <- vapply(parts, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
    pair_df$pair_label <- paste0(pair_df$residue_A, " ↔ ", pair_df$residue_B)
    pair_df <- pair_df[, c("residue_A","residue_B","pair_label","occupancy_frames","occupancy_frac","occupancy_pct")]
    pair_df <- pair_df[order(-pair_df$occupancy_pct, pair_df$residue_A, pair_df$residue_B), , drop = FALSE]
    rownames(pair_df) <- NULL
  }

  list(ts = ts_df, occupancy = occ_df, pairs = pair_df)
}

# Choose which frame will be used as the reference for alignment
choose_reference_frame <- function(frames_df, ref_mode = c("initial", "time"), ref_time_ns = 10, combine = TRUE) {
  ref_mode <- match.arg(ref_mode)

  if (ref_mode == "initial") {
    # first available frame in the processed window
    i <- 1L
    return(list(
      segment = frames_df$segment[i],
      frame = frames_df$frame[i],
      time_ns = frames_df$time_ns[i],
      time_ns_global = frames_df$time_ns_global[i]
    ))
  }

  # time-based
  if (isTRUE(combine) && ("time_ns_global" %in% names(frames_df)) && any(!is.na(frames_df$time_ns_global))) {
    tvec <- frames_df$time_ns_global
    i <- which.min(abs(tvec - ref_time_ns))
  } else {
    # without combine, time is ambiguous across segments (it restarts per file) — use first segment only
    first_seg <- frames_df$segment[1]
    sub <- frames_df[frames_df$segment == first_seg, , drop = FALSE]
    i_sub <- which.min(abs(sub$time_ns - ref_time_ns))
    i <- sub$rowid[i_sub]
  }

  list(
    segment = frames_df$segment[i],
    frame = frames_df$frame[i],
    time_ns = frames_df$time_ns[i],
    time_ns_global = frames_df$time_ns_global[i]
  )
}

make_pdb_template_from_prmtop <- function(prm, atom_inds) {
  pdb_full <- try(as.pdb(prm), silent = TRUE)
  if (inherits(pdb_full, "try-error")) stop("Could not build PDB template from prmtop (as.pdb failed).")

  pdb_sub <- try(trim.pdb(pdb_full, inds = atom_inds), silent = TRUE)
  if (!inherits(pdb_sub, "try-error")) return(pdb_sub)

  # fallback manual subset (less feature-complete, but works)
  pdb_full$atom <- pdb_full$atom[atom_inds, , drop = FALSE]
  pdb_full
}

write_pdb_with_xyz <- function(pdb_template, xyz_vec, file) {
  pdb_out <- pdb_template
  xyz_vec <- as.numeric(xyz_vec)

  # Sanity check: XYZ length must match template atom count
  if (!is.null(pdb_out$atom) && nrow(pdb_out$atom) > 0) {
    n_atom <- nrow(pdb_out$atom)
    if (length(xyz_vec) != 3L * n_atom) {
      stop(sprintf(
        "XYZ length (%d) does not match 3*#atoms (%d). Export selection may be inconsistent.",
        length(xyz_vec), 3L * n_atom
      ))
    }
  }

  pdb_out$xyz <- xyz_vec
  write.pdb(pdb_out, file = file)
}

get_export_atom_inds <- function(prm,
                                mode = "sys_heavy",
                                exclude_resid = c("HOH","WAT","NA","Na+","CL","Cl-")) {
  pdb_full <- try(as.pdb(prm), silent = TRUE)
  if (inherits(pdb_full, "try-error") || is.null(pdb_full$atom) || nrow(pdb_full$atom) < 1) {
    stop("Could not derive atom table from prmtop for export selection.")
  }

  a <- pdb_full$atom
  keep <- rep(TRUE, nrow(a))

  # Normalize residue names for matching
  resid <- a$resid
  if (is.null(resid)) resid <- rep("", nrow(a))
  resid_up <- toupper(as.character(resid))
  excl_up  <- toupper(as.character(exclude_resid))

  # "System (no solvent/ions)" excludes residue names listed by the user
  if (startsWith(mode, "sys_") && length(excl_up) > 0) {
    keep <- keep & !(resid_up %in% excl_up)
  }

  # Use element types if available; fallback to atom names
  ele <- a$elety
  if (is.null(ele)) ele <- a$atom
  ele <- as.character(ele)

  # Heavy atoms = remove hydrogens
  if (grepl("_heavy$", mode)) {
    keep <- keep & !grepl("^H", ele)
  }

  # Backbone = N, CA, C, O (+ OXT)
  if (grepl("_backbone$", mode)) {
    bb <- c("N","CA","C","O","OXT")
    keep <- keep & (ele %in% bb)
  }

  which(keep)
}

run_rmsd_structural_clustering <- function(project_dir,
                                          trajs_ordered,
                                          combine = TRUE,
                                          out_dir = cluster_results_dir(project_dir),
                                          dt_ps = 10,
                                          first = NULL, last = NULL, stride = 1,
                                          pep_resno = integer(),
                                          lipo_resno = integer(),
                                          membrane_resid = character(),
                                          exclude_resid = c("HOH","WAT","Na+","Cl-"),
                                          align_elety = c("N","CA","C","O"),
                                          clust_atoms_key = c("A_backbone","A_heavy","A_ca","B_heavy","AplusB_backbone","AplusB_heavy"),
                                          alg = c("hclust","pam"),
                                          ref_mode = c("initial","time"),
                                          ref_time_ns = 10,
                                          ref_coords_path = NULL,
                                          ref_coords_name = NULL,
                                          k = 5,
                                          progress = NULL,
                                          logf = NULL,
                                          verbose = TRUE) {

  dir_create(out_dir)

  inp <- detect_inputs(project_dir)
  if (is.na(inp$prmtop)) stop("No topology file (.prmtop or .parm7) found in selected folder.")
  if (length(trajs_ordered) < 1) stop("No .nc trajectories selected.")

  prm <- read.prmtop(inp$prmtop)
  sels <- make_selections(prm, pep_resno, lipo_resno, membrane_resid,
                         exclude_resid = exclude_resid,
                         align_elety = align_elety)

  clust_atoms_key <- match.arg(clust_atoms_key)
  alg <- match.arg(alg)
  ref_mode <- match.arg(ref_mode)

  sel_align <- sels$align
  sel_clust <- sels[[clust_atoms_key]]

  if (length(sel_align$atom) < 3) stop("Alignment selection too small. Check Selection A residues / exclude list.")
  if (length(sel_clust$atom) < 1) stop("Clustering selection empty. Check Selection A/Selection B residues and atom choice.")

  # atom subset to read (only what is needed for alignment + clustering)
  subset_atom <- sort(unique(c(sel_align$atom, sel_clust$atom)))
  subset_atom <- clean_atom_inds(subset_atom, natom = get_prm_natom(prm))
  if (length(subset_atom) < 1) stop("Clustering atom subset is empty after cleaning. Check selections.")
  at_sel <- as.select(subset_atom)

  # Build frame table (the exact frames that will be processed)
  frames_df <- build_frame_index(trajs_ordered, dt_ps = dt_ps, first = first, last = last,
                                 stride = as.integer(stride), combine = isTRUE(combine))

  # Pick reference frame unless an external coordinate file is provided.
  if (!is.null(ref_coords_path) && nzchar(ref_coords_path)) {
    if (verbose) message("Clustering external reference: ", ref_coords_path)
    ref_info <- list(segment = "external", frame = 1L)
    ref_full <- read_reference_coords_file(ref_coords_path, original_name = ref_coords_name, natom_expected = get_prm_natom(prm))
    ref_vec <- subset_xyz_vector(ref_full, subset_atom)
  } else {
    ref_info <- choose_reference_frame(frames_df, ref_mode = ref_mode, ref_time_ns = ref_time_ns, combine = isTRUE(combine))

    # Read reference coordinates
    trajs_map <- setNames(trajs_ordered, basename(trajs_ordered))
    ref_seg_path <- trajs_map[[ref_info$segment]]
    if (is.null(ref_seg_path) || is.na(ref_seg_path)) stop("Could not locate reference segment file.")
    if (verbose) message("Clustering reference: ", ref_info$segment, " frame ", ref_info$frame)

    ref_xyz <- read.ncdf(ref_seg_path, first = ref_info$frame, last = ref_info$frame, stride = 1,
                        at.sel = at_sel, verbose = verbose)
    ref_vec <- ref_xyz[1, ]
  }

  cols_align <- cols_for_atoms(sel_align$atom, subset_atom)
  cols_clust <- cols_for_atoms(sel_clust$atom, subset_atom)

  # Read and align all frames, accumulating the clustering coordinates
  xyz_all <- NULL

  n_seg <- length(trajs_ordered)

  for (i_seg in seq_along(trajs_ordered)) {
    seg_path <- trajs_ordered[[i_seg]]
    seg_name <- basename(seg_path)
    if (!is.null(logf)) logf(sprintf("Clustering: segment %d/%d: %s", i_seg, n_seg, seg_name))
    seg_base <- 0.15 + 0.55 * (i_seg - 1) / max(1, n_seg)
    seg_span <- 0.55 / max(1, n_seg)
    if (!is.null(progress)) progress(value = seg_base, detail = sprintf("Reading %s (%d/%d)…", seg_name, i_seg, n_seg))
    sub <- frames_df[frames_df$segment == seg_name, , drop = FALSE]
    if (nrow(sub) < 1) next

    f1 <- min(sub$frame)
    fL <- max(sub$frame)

    if (verbose) message("Reading for clustering: ", seg_name, " (", nrow(sub), " frames)")
    xyz <- read.ncdf(seg_path, first = f1, last = fL, stride = as.integer(stride),
                    at.sel = at_sel, verbose = verbose)
    if (!is.matrix(xyz) || nrow(xyz) < 1) next

    if (!is.null(progress)) progress(value = seg_base + seg_span * 0.5, detail = sprintf("Aligning %s (%d frames)…", seg_name, nrow(xyz)))
    if (!is.null(logf)) logf(sprintf("  Read %d frames, aligning…", nrow(xyz)))

    xyz_fit <- kabsch_fit_xyz(fixed = ref_vec, mobile = xyz,
                              fixed.inds = cols_align, mobile.inds = cols_align)

    xyz_all <- rbind(xyz_all, xyz_fit[, cols_clust, drop = FALSE])
  }

  if (is.null(xyz_all) || nrow(xyz_all) < 2) stop("Not enough frames read for clustering.")

  if (!is.null(logf)) logf("Computing RMSD distance matrix (this can take a while for many frames)...")
  if (!is.null(progress)) progress(value = 0.75, detail = "Computing RMSD distance matrix")

  # RMSD distance matrix (frames already aligned to the chosen reference)
  # dist() gives Euclidean distance; convert to RMSD by dividing by sqrt(n_atoms)
  ncoords <- ncol(xyz_all)
  n_atoms_clust <- ncoords / 3L
  d_euc <- dist(xyz_all)
  d_rmsd <- d_euc / sqrt(n_atoms_clust)

  if (!is.null(logf)) logf("Running clustering...")
  if (!is.null(progress)) progress(value = 0.85, detail = "Running clustering")

  # Clustering
  hc <- NULL
  if (alg == "hclust") {
    hc <- hclust(d_rmsd, method = "ward.D2")
    cl_raw <- cutree(hc, k = as.integer(k))
  } else {
    if (!requireNamespace("cluster", quietly = TRUE)) {
      stop("Package 'cluster' is required for PAM (k-medoids). Please install it or use Hierarchical clustering.")
    }
    pam_fit <- cluster::pam(d_rmsd, k = as.integer(k), diss = TRUE)
    cl_raw <- pam_fit$clustering
  }

  # Rename clusters by population: 1 = most populated
  tab <- sort(table(cl_raw), decreasing = TRUE)
  map <- setNames(seq_along(tab), names(tab))
  cl <- as.integer(map[as.character(cl_raw)])

  frames_df$cluster <- cl

  if (!is.null(logf)) logf("Embedding clusters for visualization + computing medoids/centroids...")
  if (!is.null(progress)) progress(value = 0.92, detail = "Embedding + representatives")

  # Embed in 3D for plotting (MDS on RMSD distances; fallback to PCA)
  embed <- try(cmdscale(d_rmsd, k = 3), silent = TRUE)
  if (inherits(embed, "try-error") || is.null(embed)) {
    pc <- prcomp(xyz_all, center = TRUE, scale. = FALSE)
    embed <- pc$x[, 1:3, drop = FALSE]
  }
  plot_df <- data.frame(
    Dim1 = embed[, 1],
    Dim2 = embed[, 2],
    Dim3 = embed[, 3],
    cluster = as.factor(cl),
    segment = frames_df$segment,
    frame = frames_df$frame,
    time_ns = frames_df$time_ns,
    time_ns_global = frames_df$time_ns_global
  )

  # Compute medoids + centroids (on clustering coordinates)
  k_eff <- max(cl)
  medoid_xyz <- vector("list", k_eff)
  centroid_xyz <- vector("list", k_eff)
  centroid_frame_xyz <- vector("list", k_eff)

  sil_df <- NULL
  sil_summary <- data.frame(cluster = seq_len(k_eff), mean_silhouette = rep(NA_real_, k_eff), stringsAsFactors = FALSE)
  silhouette_overall <- NA_real_
  if (length(unique(cl)) > 1 && length(cl) >= 3 && requireNamespace("cluster", quietly = TRUE)) {
    sil_try <- try(cluster::silhouette(cl, d_rmsd), silent = TRUE)
    if (!inherits(sil_try, "try-error")) {
      sil_df <- as.data.frame(unclass(sil_try), stringsAsFactors = FALSE)
      if (ncol(sil_df) >= 3) {
        names(sil_df)[1:3] <- c("cluster", "neighbor_cluster", "sil_width")
        sil_df$cluster <- suppressWarnings(as.integer(sil_df$cluster))
        sil_df$rowid <- seq_len(nrow(sil_df))
        sil_df$segment <- frames_df$segment
        sil_df$frame <- frames_df$frame
        sil_df$time_ns <- frames_df$time_ns
        sil_df$time_ns_global <- frames_df$time_ns_global
        silhouette_overall <- mean(sil_df$sil_width, na.rm = TRUE)
        sil_summary <- stats::aggregate(sil_width ~ cluster, data = sil_df, FUN = function(x) mean(x, na.rm = TRUE))
        names(sil_summary)[2] <- "mean_silhouette"
      }
    }
  }

  summ <- data.frame(
    cluster = seq_len(k_eff),
    n_frames = integer(k_eff),
    pct = numeric(k_eff),
    medoid_rowid = integer(k_eff),
    medoid_segment = character(k_eff),
    medoid_frame = integer(k_eff),
    medoid_time_ns = numeric(k_eff),
    medoid_time_ns_global = numeric(k_eff),
    centroid_frame_rowid = integer(k_eff),
    centroid_frame_segment = character(k_eff),
    centroid_frame = integer(k_eff),
    centroid_frame_time_ns = numeric(k_eff),
    centroid_frame_time_ns_global = numeric(k_eff),
    within_avg_pair_rmsd = numeric(k_eff),
    centroid_rmsd_to_medoid = numeric(k_eff),
    stringsAsFactors = FALSE
  )

  for (cc in seq_len(k_eff)) {
    idxs <- which(cl == cc)
    summ$n_frames[cc] <- length(idxs)
    summ$pct[cc] <- 100 * length(idxs) / length(cl)

    mdat <- medoid_from_dist(d_rmsd, idxs)
    med <- mdat$medoid
    summ$medoid_rowid[cc] <- med
    summ$within_avg_pair_rmsd[cc] <- mdat$avgpair

    medoid_xyz[[cc]] <- xyz_all[med, ]
    centroid_xyz[[cc]] <- colMeans(xyz_all[idxs, , drop = FALSE])

    # centroid-nearest real frame (useful for exporting a physically existing structure)
    if (length(idxs) >= 1) {
      dif <- xyz_all[idxs, , drop = FALSE] - matrix(centroid_xyz[[cc]], nrow = length(idxs), ncol = ncoords, byrow = TRUE)
      d2 <- rowSums(dif * dif)
      near_local <- which.min(d2)
      near <- idxs[near_local]
      centroid_frame_xyz[[cc]] <- xyz_all[near, ]
      summ$centroid_frame_rowid[cc] <- near
      summ$centroid_frame_segment[cc] <- frames_df$segment[near]
      summ$centroid_frame[cc] <- frames_df$frame[near]
      summ$centroid_frame_time_ns[cc] <- frames_df$time_ns[near]
      summ$centroid_frame_time_ns_global[cc] <- frames_df$time_ns_global[near]
    }
    summ$centroid_rmsd_to_medoid[cc] <- rmsd_vec(centroid_xyz[[cc]], medoid_xyz[[cc]])

    summ$medoid_segment[cc] <- frames_df$segment[med]
    summ$medoid_frame[cc] <- frames_df$frame[med]
    summ$medoid_time_ns[cc] <- frames_df$time_ns[med]
    summ$medoid_time_ns_global[cc] <- frames_df$time_ns_global[med]
  }

  # Provide export-ready template info (clustering atoms)
  clust_atom_inds <- sort(sel_clust$atom)
  summ <- merge(summ, sil_summary, by = "cluster", all.x = TRUE, sort = FALSE)
  summ <- summ[order(summ$cluster), , drop = FALSE]

  list(
    frames = frames_df,
    summary = summ,
    plot = plot_df,
    prmtop_path = inp$prmtop,
    trajs_map = trajs_map,
    clust_atom_inds = clust_atom_inds,
    clust_atoms_key = clust_atoms_key,
    alg = alg,
    hc = hc,
    ref_info = ref_info,
    medoid_xyz = medoid_xyz,
    centroid_xyz = centroid_xyz,
    centroid_frame_xyz = centroid_frame_xyz,
    silhouette = sil_df,
    silhouette_summary = sil_summary,
    silhouette_overall = silhouette_overall,
    d_rmsd = d_rmsd
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# UI (preserved layout; added "Combine + Order" controls)
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# Per-plot styling + export UI (local controls next to each plot)
# ─────────────────────────────────────────────────────────────────────────────
plot_style_export_ui <- function(prefix, title = "Plot styling + export",
                                show_points_default = FALSE,
                                base_size_default = 12, line_width_default = 1.2,
                                point_size_default = 3, point_alpha_default = 0.8,
                                fmt_default = "png", width_default = 7, height_default = 4, dpi_default = 300) {
  tags$details(
    class = "plot-style",
    tags$summary(title),
    selectInput(paste0(prefix, "_pub_theme"), "Theme",
                choices = c("Classic" = "classic", "Minimal" = "minimal", "BW" = "bw"),
                selected = "classic"),
    numericInput(paste0(prefix, "_pub_base_size"), "Base font size", value = base_size_default, min = 6),
    numericInput(paste0(prefix, "_line_width"), "Line width (interactive + export)",
                 value = line_width_default, min = 0.2, step = 0.1),
    checkboxInput(paste0(prefix, "_show_points"), "Show points", value = show_points_default),
    numericInput(paste0(prefix, "_point_size"), "Point size", value = point_size_default, min = 1),
    sliderInput(paste0(prefix, "_point_alpha"), "Point opacity", min = 0.1, max = 1,
                value = point_alpha_default, step = 0.1),
    selectInput(paste0(prefix, "_palette"), "Palette",
                choices = c("Default" = "default", "Dark3" = "dark3", "Set2" = "set2", "Viridis" = "viridis"),
                selected = "default"),
    hr(),
    tags$div(class="section-header","Axes"),
    textInput(paste0(prefix, "_x_title"), "X-axis title (optional)", value = ""),
    textInput(paste0(prefix, "_y_title"), "Y-axis title (optional)", value = ""),
    selectInput(paste0(prefix, "_x_scale"), "X scale",
                choices = c("Linear" = "linear", "Log10" = "log"),
                selected = "linear"),
    selectInput(paste0(prefix, "_y_scale"), "Y scale",
                choices = c("Linear" = "linear", "Log10" = "log"),
                selected = "linear"),
    fluidRow(
      column(6, numericInput(paste0(prefix, "_x_min"), "X min (optional)", value = NA)),
      column(6, numericInput(paste0(prefix, "_x_max"), "X max (optional)", value = NA))
    ),
    fluidRow(
      column(6, numericInput(paste0(prefix, "_y_min"), "Y min (optional)", value = NA)),
      column(6, numericInput(paste0(prefix, "_y_max"), "Y max (optional)", value = NA))
    ),
    hr(),
    selectInput(paste0(prefix, "_export_fmt"), "Export format",
                choices = c("PDF" = "pdf", "PNG" = "png", "SVG" = "svg"), selected = fmt_default),
    numericInput(paste0(prefix, "_export_width"), "Export width (in)", value = width_default, min = 3),
    numericInput(paste0(prefix, "_export_height"), "Export height (in)", value = height_default, min = 3),
    numericInput(paste0(prefix, "_export_dpi"), "PNG DPI", value = dpi_default, min = 72)
  )
}

# Small, app-themed Plotly theme toggle (dark default)
plotly_theme_toggle_ui <- function(id, theme = "dark") {
  lab <- if (identical(theme, "light")) "Plot theme: Light" else "Plot theme: Dark"
  ic  <- if (identical(theme, "light")) icon("sun") else icon("moon")
  actionButton(
    id, lab, icon = ic,
    class = "btn-primary btn-compact",
    style = "width:100%; margin-bottom:14px;"
  )
}

info_button_ui <- function(id) {
  actionButton(
    id, label = "", icon = icon("info-circle"),
    class = "info-btn",
    style = "margin-left:10px; vertical-align:middle;"
  )
}

box_title_with_info <- function(title_text, info_id) {
  tagList(
    tags$span(title_text),
    info_button_ui(info_id)
  )
}

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "⚛ ShineMD"),
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Project", tabName = "project", icon = icon("folder-open")),
      menuItem("RMSD", tabName = "rmsd", icon = icon("chart-line")),
      menuItem("RMSF", tabName = "rmsf", icon = icon("wave-square")),
      menuItem("Radius of Gyration", tabName = "rg", icon = icon("circle-notch")),
      menuItem("PCA / DimRed", tabName = "dimred", icon = icon("braille")),
      menuItem("Membrane systems", tabName = "mem", icon = icon("arrows-alt-v")),
      menuItem("Interactions", tabName = "interactions", icon = icon("link")),
      menuItem("Clustering", tabName = "cluster", icon = icon("project-diagram"))
      ,menuItem("Session / Reprod.", tabName = "sessiontab", icon = icon("flask"))
      ,menuItem("About / Contact", tabName = "about", icon = icon("address-card"))
    )
  ),
  dashboardBody(
    useShinyjs(),
    tags$head(
      tags$style(HTML(custom_css)),
      tags$style(HTML("
    .about-inst-logos,
    .about-inst-logos img,
    .about-lab-logo,
    .about-lab-logo img {
      background: transparent !important;
      border: none !important;
      box-shadow: none !important;
      outline: none !important;
    }

    .about-muted {
      color: #93a8c7;
    }

    .about-subtle {
      color: #cfd8e3;
    }

    .about-box-text {
      line-height: 1.75;
    }

    .about-inst-panel {
      background: transparent !important;
      border-radius: 0 !important;
      padding: 0 !important;
      display: inline-flex;
      align-items: center;
      justify-content: center;
      gap: 18px;
      border: none !important;
      box-shadow: none !important;
    }

    .about-inst-panel img {
      background: transparent !important;
      border: none !important;
      box-shadow: none !important;
      outline: none !important;
      display: block;
    }

    .about-inst-sep {
      width: 1px;
      height: 54px;
      background: rgba(255,255,255,0.10);
    }
    .about-shinemd-title {
      display: inline-block;
      margin-top: 0;
      margin-bottom: 14px;
      font-size: 68px;
      font-weight: 800;
      line-height: 1.02;
      letter-spacing: -0.6px;
      background: linear-gradient(90deg, #1ea6c6 0%, #20d29a 100%);
      -webkit-background-clip: text;
      background-clip: text;
      -webkit-text-fill-color: transparent;
      color: transparent;
      text-shadow: 0 0 14px rgba(25, 224, 255, 0.08);
    }
  "))
    ),
    tabItems(
      tabItem(
        tabName = "project",
        fluidRow(
          column(6,
                 box(title = "Project folder (contains .prmtop/.parm7 + .nc segments)", width = 12, status = "primary",
                     div(class="section-header","Select folder"),
                     shinyDirButton("proj_dir_btn", "Browse…", "Select a folder"),
                     textInput("proj_dir", "Or paste folder path", value = ""),
                     uiOutput("proj_status"),
                     hr(),
                     div(class="section-header","Detected inputs"),
                     uiOutput("detected_files")
                 ),
                 box(title = "Status / Log", width = 12, status = "warning",
                     verbatimTextOutput("log")
                 )
          ),
          column(6,
                 box(title = "Run settings", width = 12, status = "info",
                     div(class="section-header","Trajectory/time"),
                     numericInput("dt_ps", "dt per saved frame (ps)", value = 10, min = 0.001),
                     numericInput("stride", "Stride (analyze every N frames)", value = 1, min = 1),
                     numericInput("first", "First frame (1-based; blank = start)", value = NA, min = 1),
                     numericInput("last", "Last frame (1-based; blank = end)", value = NA, min = 1),
                     hr(),
                     div(class="section-header","Combine consecutive segments"),
                     checkboxInput("combine_trajs", "Combine segments into one continuous timeline", value = TRUE),
                     textAreaInput("traj_order", "Trajectory order (one filename per line)",
                                   value = "", rows = 8, placeholder = "seg1.nc
seg2.nc
seg3.nc"),
                     tags$small("Tip: reorder the lines to change the concatenation order."),
                     hr(),
                     div(class="section-header","Selections (A/B)"),
                     textInput("pep_resno", "Selection A residue numbers (protein/peptide; optional; blank = auto macromolecule)", value = "", placeholder = "e.g. 1:535"),
                     textInput("lipo_resno", "Selection B residue numbers (ligand/lipids/other; optional)", value = "", placeholder = "e.g. 536:555"),
                     tags$small(style="color: var(--text-muted);", "Selection A is the primary alignment/analysis selection (typically protein/peptide). Selection B is optional (ligand, lipids, another domain). If left blank, the app auto-detects the macromolecule as Selection A."),
                     textInput("exclude_resid", "Exclude residue names (comma-separated; e.g., HOH,WAT,Na+,Cl-)", value = "HOH,WAT,Na+,Cl-", placeholder = "e.g. HOH,WAT,Na+,Cl-"),
                     selectInput("align_mode", "Alignment selection", choices = c("Selection A backbone (N,CA,C,O)" = "bb",
                                                                                "Selection A Cα only" = "ca"), selected = "bb"),
                     selectInput("rmsd_ref_mode", "RMSD reference mode", choices = c("First frame of selected subset (global)" = "global_first",
                                                                                         "First frame of each segment" = "segment_first"), selected = "global_first"),
                     checkboxInput("use_rmsd_ref_file", "Use external coordinate reference for RMSD / PCA / clustering", value = FALSE),
                     fileInput("rmsd_ref_file", "Reference coordinates (PDB, rst7, inpcrd, crd)", accept = c(".pdb", ".ent", ".rst7", ".inpcrd", ".crd", ".rst")),
                     tags$small(style="color: var(--text-muted); display:block; margin-top:4px;", "Example cpptraj-like workflow: upload the initial tleap coordinates (for example system.rst7). When an external reference is provided, it overrides the internal first-frame RMSD reference."),
                     selectInput("mem_target", "Membrane metrics target", choices = c("Selection B (if provided)" = "B",
                                                                                     "Selection A" = "A"), selected = "B"),
                     textInput("mem_resid", "Membrane residue names (optional; comma-separated)", value = "", placeholder = "e.g. PA,PE,OL,ERG"),
                     checkboxInput("mem_center_profiles", "Center membrane at z = 0 for density/order calculations", value = TRUE),
                     textInput("water_resid", "Water residue names (for membrane profiles)", value = "HOH,WAT,TIP3,TIP3P,SOL", placeholder = "e.g. HOH,WAT,TIP3,TIP3P,SOL"),
                     textInput("ion_resid", "Ion residue names (for membrane profiles)", value = "Na+,K+,Cl-,CLA,SOD,POT,MG2,CA", placeholder = "e.g. Na+,K+,Cl-,CLA,SOD,POT"),
                     tags$small(style="color: var(--text-muted); display:block; margin-top:6px;", "Examples: Selection A = 1:535, Selection B = 536:555, membrane residues = PA,PE,OL,ERG."),
                     tags$small(style="color: var(--text-muted); display:block; margin-top:4px;", "For RMSD/PCA/clustering, trajectories should ideally be whole-molecule and autoimaged before loading."),
                     hr(),
                     actionButton("run_all", "Run basic analysis", class="btn-primary", icon=icon("play"), width="100%")
                 )
          )
        )
      ),

      tabItem(
        tabName = "rmsd",
        fluidRow(
          column(12,
                 box(title=box_title_with_info("RMSD — Selection A", "info_rmsdA"), width=12, status="primary",
                     fluidRow(
                       column(8, plotlyOutput("plot_rmsd_pep", height="360px")),
                       column(4,
                              plotly_theme_toggle_ui("toggle_plotly_theme_rmsd"),
                              selectInput("rmsdA_series", "RMSD series", choices = c("Backbone + heavy atoms" = "both", "Backbone only" = "backbone", "Heavy atoms only" = "heavy"), selected = "both"),
                              numericInput("rmsdA_smooth", "Smoothing window (frames; 0 = off)", value = 0, min = 0, step = 10),
                              plot_style_export_ui("rmsdA", show_points_default = FALSE),
                              tags$br(),
                              downloadButton("dl_rmsdA_plot", "Download plot", class="btn-primary"),
                              downloadButton("dl_rmsdA_data", "Download plot data (CSV)", class="btn-default")
                       )
                     ),
                     uiOutput("rmsdA_stats_row")
)
          )
        ),
        fluidRow(
          column(12,
                 box(title=box_title_with_info("RMSD distribution (Selection A)", "info_rmsdDist"), width=12, status="info",
                     solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                     fluidRow(
                       column(8, plotlyOutput("plot_rmsd_dist", height="320px")),
                       column(4,
                              numericInput("rmsdDist_nbins", "Number of bins", value = 60, min = 10, max = 200, step = 5),
                              plot_style_export_ui("rmsdDist", show_points_default = FALSE),
                              tags$br(),
                              downloadButton("dl_rmsdDist_plot", "Download plot", class="btn-primary")
                       )
                     )
                 )
          )
        ),
        uiOutput("ui_rmsd_B")
      ),

      tabItem(
        tabName = "rmsf",
        fluidRow(
          column(12,
                 box(title=box_title_with_info("RMSF (Selection A, per residue)", "info_rmsf"), width=12, status="warning",
                     fluidRow(
                       column(8, plotlyOutput("plot_rmsf", height="450px")),
                       column(4,
                              plotly_theme_toggle_ui("toggle_plotly_theme_rmsf"),
                              plot_style_export_ui("rmsf", show_points_default = TRUE, line_width_default = 1.0),
                              tags$br(),
                              downloadButton("dl_rmsf_plot", "Download plot", class="btn-primary"),
                              downloadButton("dl_rmsf_data", "Download plot data (CSV)", class="btn-default")
                       )
                     )
                 )
          )
        ),
        uiOutput("ui_rmsf_B")
      ),

      tabItem(
        tabName = "rg",
        fluidRow(
          column(12,
                 box(title=box_title_with_info("Radius of gyration — Selection A", "info_rgA"), width=12, status="danger",
                     fluidRow(
                       column(8, plotlyOutput("plot_rg_pep", height="360px")),
                       column(4,
                              plotly_theme_toggle_ui("toggle_plotly_theme_rg"),
                              numericInput("rgA_smooth", "Smoothing window (frames; 0 = off)", value = 0, min = 0, step = 10),
                              plot_style_export_ui("rgA", show_points_default = FALSE),
                              tags$br(),
                              downloadButton("dl_rgA_plot", "Download plot", class="btn-primary"),
                              downloadButton("dl_rgA_data", "Download plot data (CSV)", class="btn-default")
                       )
                     ),
                     uiOutput("rgA_stats_row")
)
          )
        ),
        uiOutput("ui_rg_B")
      ),


      tabItem(
        tabName = "dimred",
        fluidRow(
          column(12,
                 box(title=box_title_with_info("PCA / Dimensionality reduction (Selection A Cα)", "info_dimred"), width=12, status="primary",
                     fluidRow(
                       column(8, plotlyOutput("plot_dimred", height="480px")),
                       column(4,
                              plotly_theme_toggle_ui("toggle_plotly_theme_dimred"),
                              plot_style_export_ui("dimred", show_points_default = TRUE),
                              tags$br(),
                              actionButton("toggle_dimred_3d", "Toggle 3D PCA plot", icon=icon("cube"),
                                           class="btn-primary", width="100%"),
                              tags$br(),
                              tags$br(),
                              downloadButton("dl_dimred_plot", "Download plot", class="btn-primary"),
                              downloadButton("dl_dimred_data", "Download plot data (CSV)", class="btn-default")
                       )
                     )
                 )
          )
        ),
        fluidRow(
          column(12,
                 box(title=box_title_with_info("Explained variance", "info_pcaVar"), width=12, status="warning",
                     fluidRow(
                       column(8, plotlyOutput("plot_pca_var", height="320px")),
                       column(4,
                              plotly_theme_toggle_ui("toggle_plotly_theme_pcaVar"),
                              plot_style_export_ui("pcaVar", show_points_default = TRUE, line_width_default = 0.8),
                              tags$br(),
                              downloadButton("dl_pcaVar_plot", "Download plot", class="btn-primary"),
                              downloadButton("dl_pcaVar_data", "Download plot data (CSV)", class="btn-default")
                       )
                     )
                 )
          )
        ),
        fluidRow(
          column(12,
                 box(title=box_title_with_info("Free energy landscape (PC1 vs PC2)", "info_fel"), width=12, status="danger",
                     solidHeader = FALSE,
                     fluidRow(
                       column(8, plotlyOutput("plot_fel", height="480px")),
                       column(4,
                              numericInput("fel_nbins", "Grid bins", value = 80, min = 20, max = 200, step = 10),
                              numericInput("fel_temp", "Temperature (K)", value = 300, min = 1),
                              numericInput("fel_maxG", "Max ΔG shown (kcal/mol; blank = auto)", value = NA, min = 0.1, step = 0.5),
                              plot_style_export_ui("fel", show_points_default = FALSE, fmt_default = "png"),
                              tags$br(),
                              downloadButton("dl_fel_plot", "Download plot", class="btn-primary"),
                              downloadButton("dl_fel_data", "Download FEL data (CSV)", class="btn-default")
                       )
                     )
                 )
          )
        )
      ),

      tabItem(
        tabName = "mem",
        fluidRow(
          column(12,
                 box(title=box_title_with_info("Membrane COM distance", "info_memDist"), width=12, status="info",
                     fluidRow(
                       column(8, plotlyOutput("plot_dist", height="360px")),
                       column(4,
                              plotly_theme_toggle_ui("toggle_plotly_theme_memDist"),
                              plot_style_export_ui("memDist", show_points_default = FALSE),
                              tags$br(),
                              downloadButton("dl_memDist_plot", "Download plot", class="btn-primary"),
                              downloadButton("dl_memDist_data", "Download plot data (CSV)", class="btn-default")
                       )
                     )
)
          )
        ),
        fluidRow(
          column(12,
                 box(title=box_title_with_info("Membrane Δz (COM)", "info_memZ"), width=12, status="info",
                     fluidRow(
                       column(8, plotlyOutput("plot_z", height="360px")),
                       column(4,
                              plotly_theme_toggle_ui("toggle_plotly_theme_memZ"),
                              plot_style_export_ui("memZ", show_points_default = FALSE),
                              tags$br(),
                              downloadButton("dl_memZ_plot", "Download plot", class="btn-primary"),
                              downloadButton("dl_memZ_data", "Download plot data (CSV)", class="btn-default")
                       )
                     )
)
          )
        ),

        fluidRow(
          column(6,
                 box(title=box_title_with_info("Bilayer thickness (headgroup z)", "info_memThick"), width=12, status="primary",
                     fluidRow(
                       column(8, plotlyOutput("plot_thickness", height="320px")),
                       column(4,
                              plotly_theme_toggle_ui("toggle_plotly_theme_memThick"),
                              textInput("thick_head_atoms", "Headgroup atom name(s) (comma-separated)", value = "P"),
                              tags$small(style="color: var(--text-muted);", "Tip: use P for phospholipids. Thickness is computed as mean(z_upper) - mean(z_lower) per frame."),
                              actionButton("compute_thickness", "Compute thickness", icon = icon("calculator"),
                                           class="btn-primary", width="100%"),
                              tags$br(), tags$br(),
                              plot_style_export_ui("memThick", show_points_default = FALSE),
                              tags$br(),
                              div(class="download-stack",
                              downloadButton("dl_thick_plot", "Download plot", class="btn-primary"),
                              downloadButton("dl_thick_data", "Download data (CSV)", class="btn-default")
                              )
                       )
                     )
                 )
          ),
          column(6,
                 box(title=box_title_with_info("Area per lipid (APL)", "info_memAPL"), width=12, status="primary",
                     fluidRow(
                       column(8, plotlyOutput("plot_apl", height="320px")),
                       column(4,
                              plotly_theme_toggle_ui("toggle_plotly_theme_memAPL"),
                              numericInput("apl_nlip_leaflet", "Lipids per leaflet (blank = auto total/2)", value = NA, min = 1),
                              verbatimTextOutput("apl_info"),
                              actionButton("compute_apl", "Compute APL", icon = icon("calculator"),
                                           class="btn-primary", width="100%"),
                              tags$br(), tags$br(),
                              plot_style_export_ui("memAPL", show_points_default = FALSE),
                              tags$br(),
                              div(class="download-stack",
                              downloadButton("dl_apl_plot", "Download plot", class="btn-primary"),
                              downloadButton("dl_apl_data", "Download data (CSV)", class="btn-default")
                              )
                       )
                     )
                 )
          )
        ),

        fluidRow(
          column(12,
                 box(title=box_title_with_info("Lipid enrichment around target (headgroups within cutoff)", "info_memEnrich"), width=12, status="warning",
                     fluidRow(
                       column(8, plotlyOutput("plot_enrich", height="420px")),
                       column(4,
                              plotly_theme_toggle_ui("toggle_plotly_theme_memEnrich"),
                              selectInput("enrich_target", "Target region", choices = c("Selection A"="A","Selection B"="B"), selected = "A"),
                              numericInput("enrich_cutoff", "Cutoff (Å)", value = 8, min = 1),
                              textInput("enrich_head_atom", "Headgroup atom (representative)", value = "P"),
                              checkboxInput("enrich_by_type", "Separate by lipid type (residue name)", value = TRUE),
                              actionButton("compute_enrich", "Compute enrichment", icon = icon("bullseye"),
                                           class="btn-primary", width="100%"),
                              tags$br(), tags$br(),
                              plot_style_export_ui("memEnrich", show_points_default = FALSE),
                              tags$br(),
                              downloadButton("dl_enrich_plot", "Download plot", class="btn-primary"),
                              downloadButton("dl_enrich_data", "Download data (CSV)", class="btn-default")
                       )
                     )
                 )
          )
        ),
        fluidRow(
          column(12,
                 box(title=box_title_with_info("Membrane density profiles (headgroups / tails / water / ions / target)", "info_memDensity"), width=12, status="info",
                     fluidRow(
                       column(8, plotlyOutput("plot_mem_density", height="340px")),
                       column(4,
                              plotly_theme_toggle_ui("toggle_plotly_theme_memDensity"),
                              textInput("dens_head_atoms", "Headgroup atom name(s)", value = "P"),
                              textInput("dens_tail_regex", "Tail atom regex", value = "^C2[0-9]+$|^C3[0-9]+$"),
                              selectInput("dens_target", "Include target density", choices = c("No target"="none","Selection A"="A","Selection B"="B"), selected = "none"),
                              numericInput("dens_binwidth", "Bin width (Å)", value = 1, min = 0.1, step = 0.1),
                              numericInput("dens_zlim", "Half-profile range |z| (Å; blank = auto)", value = NA, min = 5),
                              actionButton("compute_density", "Compute density profiles", icon = icon("chart-area"),
                                           class="btn-primary", width="100%"),
                              tags$br(), tags$br(),
                              checkboxGroupInput("dens_exclude_groups",
                                                 "Exclude groups from plot/export:",
                                                 choices  = c("Water", "Ions", "Headgroups", "Tails", "Target"),
                                                 selected = "Water",
                                                 inline   = FALSE),
                              tags$br(),
                              plot_style_export_ui("memDensity", show_points_default = FALSE),
                              tags$br(),
                              div(class="download-stack",
                                  downloadButton("dl_density_plot", "Download plot", class="btn-primary"),
                                  downloadButton("dl_density_data", "Download data (CSV)", class="btn-default"))
                       )
                     )
                 )
          )
        ),
        fluidRow(
          column(12,
                 box(title=box_title_with_info("Lipid tail order profile |S| (pure-R segmental approximation)", "info_memOrder"), width=12, status="warning",
                     fluidRow(
                       column(8, plotlyOutput("plot_mem_order", height="340px")),
                       column(4,
                              plotly_theme_toggle_ui("toggle_plotly_theme_memOrder"),
                              selectInput("order_topology_mode", "Tail-order topology mode",
                                          choices = c("Standard: both chains inside the same lipid residue" = "standard",
                                                      "Split residues: chain 1 and chain 2 are different residue names" = "split"),
                                          selected = "standard"),
                              textInput("order_resid_filter", "Residue names for tail order (blank = all membrane residues)", value = "", placeholder = "e.g. PA,PE,OL"),
                              textInput("order_chain1_resid", "sn1 / chain 1 residue names (for split mode)", value = "", placeholder = "e.g. PA"),
                              textInput("order_chain2_resid", "sn2 / chain 2 residue names (for split mode)", value = "", placeholder = "e.g. OL"),
                              textInput("order_chain1_regex", "sn1 chain atom regex", value = "^C2[0-9]+$", placeholder = "e.g. ^C2[0-9]+$"),
                              textInput("order_chain2_regex", "sn2 chain atom regex", value = "^C3[0-9]+$", placeholder = "e.g. ^C3[0-9]+$"),
                              checkboxInput("order_by_type", "Separate by lipid type", value = TRUE),
                              tags$small(style="color: var(--text-muted); display:block; margin-top:-6px; margin-bottom:10px;", "Tip: sterols such as ERG/CHL usually should be excluded from this acyl-tail approximation. For split topologies like PA + PE + OL, use split mode with chain 1 residue names = PA and chain 2 residue names = OL. In split mode the app now removes invalid cross-labels such as PA — sn2 or OL — sn1."),
                              actionButton("compute_order", "Compute tail order", icon = icon("wave-square"),
                                           class="btn-primary", width="100%"),
                              tags$br(), tags$br(),
                              plot_style_export_ui("memOrder", show_points_default = TRUE),
                              tags$br(),
                              div(class="download-stack",
                                  downloadButton("dl_order_plot", "Download plot", class="btn-primary"),
                                  downloadButton("dl_order_data", "Download data (CSV)", class="btn-default"))
                       )
                     )
                 )
          )
        )
      ),

      tabItem(
        tabName = "interactions",
        fluidRow(
          column(12,
                 box(title=box_title_with_info("Selection A / Selection B interaction time series", "info_intTs"), width=12, status="primary",
                     fluidRow(
                       column(8, plotlyOutput("plot_interact_ts", height="340px")),
                       column(4,
                              plotly_theme_toggle_ui("toggle_plotly_theme_intTs"),
                              selectInput("interact_selA", "Selection A atoms", choices = c("Heavy atoms"="A_heavy","Backbone (N,CA,C,O)"="A_backbone","Cα only"="A_ca"), selected = "A_heavy"),
                              selectInput("interact_selB", "Selection B atoms", choices = c("Heavy atoms"="B_heavy","Backbone (N,CA,C,O)"="B_backbone","Cα only"="B_ca"), selected = "B_heavy"),
                              numericInput("interact_cutoff", "Contact cutoff (Å)", value = 4.5, min = 1, step = 0.1),
                              selectInput("interact_metric", "Metric to display", choices = c("Minimum heavy-atom distance"="min_dist_A","COM distance"="com_dist_A","Atom contacts"="atom_contacts"), selected = "min_dist_A"),
                              actionButton("compute_interactions", "Compute A/B interactions", icon = icon("link"),
                                           class="btn-primary", width="100%"),
                              tags$br(), tags$br(),
                              plot_style_export_ui("intTs", show_points_default = FALSE),
                              tags$br(),
                              div(class="download-stack",
                                  downloadButton("dl_interact_ts_plot", "Download plot", class="btn-primary"),
                                  downloadButton("dl_interact_ts_data", "Download data (CSV)", class="btn-default"))
                       )
                     )
                 )
          )
        ),
        fluidRow(
          column(12,
                 box(title=box_title_with_info("Residue contact occupancy summary", "info_intOcc"), width=12, status="warning",
                     fluidRow(
                       column(8, plotlyOutput("plot_interact_occ", height="360px")),
                       column(4,
                              plotly_theme_toggle_ui("toggle_plotly_theme_intOcc"),
                              selectInput("interact_occ_side", "Show occupancy for", choices = c("Selection A"="Selection A","Selection B"="Selection B"), selected = "Selection A"),
                              numericInput("interact_occ_topn", "Top residues", value = 20, min = 5, max = 100, step = 1),
                              plot_style_export_ui("intOcc", show_points_default = TRUE),
                              tags$br(),
                              div(class="download-stack",
                                  downloadButton("dl_interact_occ_plot", "Download plot", class="btn-primary"),
                                  downloadButton("dl_interact_occ_data", "Download data (CSV)", class="btn-default"))
                       )
                     ),
                     DTOutput("tbl_interact_occ")
                 )
          )
        ),
        fluidRow(
          column(12,
                 box(title=box_title_with_info("A / B residue contact map", "info_intPair"), width=12, status="info",
                     fluidRow(
                       column(8, plotlyOutput("plot_interact_pair", height="430px")),
                       column(4,
                              plotly_theme_toggle_ui("toggle_plotly_theme_intPair"),
                              numericInput("interact_pair_topA", "Top Selection A residues", value = 20, min = 5, max = 100, step = 1),
                              numericInput("interact_pair_topB", "Top Selection B residues", value = 12, min = 3, max = 100, step = 1),
                              numericInput("interact_pair_minpct", "Min pair occupancy (%)", value = 5, min = 0, max = 100, step = 0.5),
                              plot_style_export_ui("intPair", show_points_default = FALSE),
                              tags$br(),
                              div(class="download-stack",
                                  downloadButton("dl_interact_pair_plot", "Download plot", class="btn-primary"),
                                  downloadButton("dl_interact_pair_data", "Download data (CSV)", class="btn-default"))
                       )
                     ),
                     DTOutput("tbl_interact_pairs")
                 )
          )
        ),
        fluidRow(
          column(12,
                 box(title=box_title_with_info("A / B hydrogen-bond summary (pure-R donor/acceptor proxy)", "info_hbTs"), width=12, status="primary",
                     fluidRow(
                       column(8, plotlyOutput("plot_hbond_ts", height="340px")),
                       column(4,
                              plotly_theme_toggle_ui("toggle_plotly_theme_hbTs"),
                              numericInput("hbond_cutoff", "Donor/acceptor cutoff (Å)", value = 3.5, min = 2, max = 5, step = 0.1),
                              checkboxInput("hbond_include_nacc", "Allow N atoms as acceptors (proxy mode)", value = FALSE),
                              selectInput("hbond_metric", "Metric to display", choices = c("Total H-bond proxies" = "hbond_count",
                                                                                             "Selection A donor → B acceptor" = "hbond_count_AtoB",
                                                                                             "Selection B donor → A acceptor" = "hbond_count_BtoA",
                                                                                             "Minimum donor–acceptor distance" = "min_da_dist"), selected = "hbond_count"),
                              actionButton("compute_hbonds", "Compute A/B H-bond proxy", icon = icon("project-diagram"),
                                           class="btn-primary", width = "100%"),
                              tags$br(), tags$br(),
                              plot_style_export_ui("hbTs", show_points_default = FALSE),
                              tags$br(),
                              div(class="download-stack",
                                  downloadButton("dl_hbond_ts_plot", "Download plot", class="btn-primary"),
                                  downloadButton("dl_hbond_ts_data", "Download data (CSV)", class="btn-default"))
                       )
                     )
                 )
          )
        ),
        fluidRow(
          column(12,
                 box(title=box_title_with_info("Persistent A / B hydrogen-bond pairs", "info_hbPair"), width=12, status="warning",
                     fluidRow(
                       column(8, plotlyOutput("plot_hbond_pair", height="360px")),
                       column(4,
                              plotly_theme_toggle_ui("toggle_plotly_theme_hbPair"),
                              selectInput("hbond_pair_direction", "Direction", choices = c("All" = "all",
                                                                                               "Selection A donor → B acceptor" = "Selection A donor → Selection B acceptor",
                                                                                               "Selection B donor → A acceptor" = "Selection B donor → Selection A acceptor"), selected = "all"),
                              numericInput("hbond_pair_topn", "Top H-bond pairs", value = 20, min = 5, max = 100, step = 1),
                              numericInput("hbond_pair_minpct", "Min pair occupancy (%)", value = 5, min = 0, max = 100, step = 0.5),
                              plot_style_export_ui("hbPair", show_points_default = TRUE),
                              tags$br(),
                              div(class="download-stack",
                                  downloadButton("dl_hbond_pair_plot", "Download plot", class="btn-primary"),
                                  downloadButton("dl_hbond_pair_data", "Download data (CSV)", class="btn-default"))
                       )
                     ),
                     DTOutput("tbl_hbond_pairs")
                 )
          )
        )
      ),

      tabItem(
        tabName = "cluster",
        fluidRow(
          column(4,
                 box(title="RMSD structural clustering", width=12, status="primary",
                     tags$small(style="color: var(--text-muted);",
                                "Tip: use the Download plot data (CSV/ZIP) button under 'Cluster plots' to export tables."),
                     tags$br(),
                     tags$br(),
                     div(class="section-header","Clustering setup"),
                     selectInput("clust_alg", "Clustering algorithm",
                                choices = c("Hierarchical (Ward.D2)"="hclust",
                                            "k-medoids (PAM)"="pam"),
                                selected = "hclust"),
                     selectInput("clust_atoms", "Atoms used for RMSD",
                                choices = c("Selection A backbone (N,CA,C,O)"="A_backbone",
                                            "Selection A heavy (no H)"="A_heavy",
                                            "Selection A Cα only"="A_ca",
                                            "Selection B heavy (no H)"="B_heavy",
                                            "Selection A+B backbone (N,CA,C,O)"="AplusB_backbone",
                                            "Selection A+B heavy (no H)"="AplusB_heavy"),
                                selected = "A_heavy"),
                     selectInput("clust_ref_mode", "Reference for alignment",
                                choices = c("Initial structure (t = 0 ns)"="initial",
                                            "Frame at a given time (ns)"="time"),
                                selected = "initial"),
                     fluidRow(
                       column(6, numericInput("clust_ref_time", "Reference time (ns)", value=10, min=0)),
                       column(6, selectInput("clust_ref_quick", "Quick time",
                                             choices=c("10"=10,"100"=100,"250"=250,"500"=500),
                                             selected=10))
                     ),
                     checkboxInput("clust_use_project_window", "Use Project tab first/last/stride window", value=TRUE),
                     numericInput("clust_stride", "Clustering stride (only if override)", value=10, min=1),
                     numericInput("clust_first", "First frame (override; NA=start)", value=NA, min=1),
                     numericInput("clust_last", "Last frame (override; NA=end)", value=NA, min=1),
                     hr(),
                     numericInput("clust_k", "Number of clusters (k)", value=5, min=2),
                     actionButton("run_clust", "Run RMSD clustering", icon=icon("sitemap"),
                                  class="btn-primary", width="100%"),
                     hr(),
                     div(class="section-header","Export (structures + tables)"),
                     numericInput("clust_export_n", "How many clusters to export (0 = all)", value = 0, min = 0),
                     selectInput("clust_export_atoms", "PDB export selection",
                                choices = c(
                                  "System (no solvent/ions), heavy (no H)" = "sys_heavy",
                                  "System (no solvent/ions), all atoms (with H)" = "sys_all",
                                  "System (no solvent/ions), backbone only (N,CA,C,O)" = "sys_backbone",
                                  "Full system, heavy (no H)" = "full_heavy",
                                  "Full system, all atoms (with H)" = "full_all"
                                ),
                                selected = "sys_heavy"),
                     tags$small(style="color: var(--text-muted);",
                                "System (no solvent/ions) excludes residue names from the 'Exclude residue names' field (e.g., HOH, WAT, Na+, Cl-)."),
                     checkboxInput("clust_export_centroid_mean", "Also export centroid mean (non-physical; uses clustering atoms)", value = FALSE),
                     actionButton("export_clust", "Export medoids + centroids", icon=icon("download"),
                                  class="btn-primary", width="100%"),
                     tags$small(style="color: var(--text-muted);",
                                "Exports go to results_ShineMD/clustering_rmsd/"),
                     hr(),
                     div(class="section-header","Export a specific frame as PDB"),
                     tags$small(style="color: var(--text-muted); display:block; margin-bottom:10px;",
                                "Requires a completed clustering run. This export uses the current clustering frame map and clustering window."),
                     selectInput("single_frame_mode", "Pick frame by",
                                 choices = c("Global time (ns)"="time",
                                             "Frame index"="segframe",
                                             "From cluster representative"="clusterrep"),
                                 selected = "time"),
                     conditionalPanel(
                       condition = "input.single_frame_mode == 'time'",
                       numericInput("single_frame_time", "Time (ns)", value = 10, min = 0)
                     ),
                     conditionalPanel(
                       condition = "input.single_frame_mode == 'segframe'",
                       numericInput("single_frame_frame", "Frame index (global if combined)", value = 1, min = 1)
                     ),
                     conditionalPanel(
                       condition = "input.single_frame_mode == 'clusterrep'",
                       numericInput("single_frame_cluster", "Cluster #", value = 1, min = 1),
                       selectInput("single_frame_rep", "Representative",
                                   choices = c("Medoid"="medoid",
                                               "Centroid-nearest frame"="centroid_frame"),
                                   selected = "medoid")
                     ),
                     actionButton("export_single_frame", "Export selected frame", icon=icon("file-export"),
                                  class="btn-primary", width="100%"),
                     tags$small(style="color: var(--text-muted); display:block; margin-top:8px;",
                                "If clustering has not been run yet, compute clustering first in this tab and then return to this export panel.")
                 )
          ),
          column(8,
                 box(title=box_title_with_info("Cluster plots", "info_clusterPlots"), width=12, status="info",
                     tabsetPanel(
                       tabPanel("Distribution",
                                fluidRow(
                                  column(8, plotlyOutput("plot_cluster_dist", height="520px")),
                                  column(4,
                                         plot_style_export_ui("clustDist", show_points_default = TRUE),
                                         tags$div(style="display:flex; align-items:center; gap:10px; margin-top:12px; flex-wrap:wrap;",
                                                  downloadButton("dl_clustDist_plot", "Download plot", class="btn-primary"),
                                                  info_button_ui("info_clustDist")
                                         )
                                  )
                                )
                       ),
                       tabPanel("Clusters vs time",
                                fluidRow(
                                  column(8, plotlyOutput("plot_cluster_time", height="260px")),
                                  column(4,
                                         plot_style_export_ui("clustTime", show_points_default = TRUE, point_size_default = 6),
                                         tags$div(style="display:flex; align-items:center; gap:10px; margin-top:12px; flex-wrap:wrap;",
                                                  downloadButton("dl_clustTime_plot", "Download plot", class="btn-primary"),
                                                  info_button_ui("info_clustTime")
                                         )
                                  )
                                )
                       ),
                       tabPanel("Population",
                                fluidRow(
                                  column(8, plotlyOutput("plot_cluster_pop", height="320px")),
                                  column(4,
                                         plot_style_export_ui("clustPop", show_points_default = TRUE, line_width_default = 0.8),
                                         tags$div(style="display:flex; align-items:center; gap:10px; margin-top:12px; flex-wrap:wrap;",
                                                  downloadButton("dl_clustPop_plot", "Download plot", class="btn-primary"),
                                                  info_button_ui("info_clustPop")
                                         )
                                  )
                                )
                       ),
                       tabPanel("Quality",
                                fluidRow(
                                  column(8, plotlyOutput("plot_cluster_quality", height="320px")),
                                  column(4,
                                         uiOutput("cluster_quality_info"),
                                         plot_style_export_ui("clustQual", show_points_default = TRUE, line_width_default = 0.8),
                                         tags$div(style="display:flex; align-items:center; gap:10px; margin-top:12px; flex-wrap:wrap;",
                                                  downloadButton("dl_clustQual_plot", "Download plot", class="btn-primary"),
                                                  info_button_ui("info_clustQual")
                                         )
                                  )
                                )
                       ),
                       tabPanel("Dendrogram",
                                fluidRow(
                                  column(8, plotOutput("plot_cluster_dend", height="360px")),
                                  column(4,
                                         tags$p(style="color: var(--text-muted); margin-top: 4px;",
                                                "Available when using hierarchical clustering."),
                                         tags$div(style="display:flex; align-items:center; gap:10px; margin-top:12px; flex-wrap:wrap;",
                                                  downloadButton("dl_clustDend_plot", "Download dendrogram (PDF)", class="btn-primary"),
                                                  info_button_ui("info_clustDend")
                                         )
                                  )
                                )
                       ),
                       tabPanel("Pairwise RMSD",
                                fluidRow(
                                  column(8, plotlyOutput("plot_rmsd_heatmap", height="520px")),
                                  column(4,
                                         tags$p(style="color: var(--text-muted); margin-top: 4px;",
                                                "Frame-vs-frame RMSD matrix. Useful for identifying conformational transitions and substates."),
                                         numericInput("rmsd_heatmap_maxval", "Max RMSD for color scale (Å; blank = auto)", value = NA, min = 0.1, step = 0.5),
                                         tags$div(style="display:flex; align-items:center; gap:10px; margin-top:12px; flex-wrap:wrap;",
                                                  downloadButton("dl_rmsdHeatmap_plot", "Download heatmap (PNG)", class="btn-primary"),
                                                  info_button_ui("info_rmsdHeatmap")
                                         )
                                  )
                                )
                       )
                     ),
                     tags$br(),
                     actionButton("toggle_plotly_theme", "Plot theme: Dark", icon=icon("moon"),
                                  class="btn-primary btn-compact", style="width:100%; margin-bottom:14px;"),
                     actionButton("toggle_clust_3d", "Toggle 3D distribution plot", icon=icon("cube"),
                                  class="btn-primary", style="width:100%; margin-bottom:14px;"),
                     downloadButton("dl_cluster_data", "Download plot data (CSV/ZIP)",
                                    class="btn-primary", style="width:100%;"),
                     tags$small(style="color: var(--text-muted);",
                                "Exports: cluster distribution, cluster-vs-time, and cluster summary tables."),
                     hr(),
                     uiOutput("cluster_ref_info")
                 ),
                 box(title=box_title_with_info("Cluster summary", "info_clusterSummary"), width=12, status="warning",
                     DTOutput("tbl_cluster_summary"),
                     tags$small(style="color: var(--text-muted);", "Tip: Use the Download plot data (CSV/ZIP) button to export full per-frame assignments and silhouette-quality tables.")
                 )
          )
        )
      ),

      tabItem(
        tabName = "sessiontab",
        fluidRow(
          box(
            width = 12, title = "Session info (reproducibility)", status = "primary", solidHeader = FALSE,
            tags$p(style = "color: var(--text-muted); margin-bottom: 12px;",
                   "Full R session information for reproducibility. Include this when reporting results or submitting a manuscript."),
            verbatimTextOutput("session_info_text"),
            tags$br(),
            downloadButton("dl_session_info", "Download session info (.txt)", class = "btn-primary btn-compact")
          )
        ),
        fluidRow(
          column(6,
            box(
              width = 12, title = "Installed package versions", status = "info", solidHeader = FALSE,
              DTOutput("pkg_versions_tbl")
            )
          ),
          column(6,
            box(
              width = 12, title = "ShineMD environment", status = "warning", solidHeader = FALSE,
              uiOutput("env_info_panel")
            )
          )
        )
      ),

      tabItem(
        tabName = "about",
        
        fluidRow(
          box(
            width = 8, title = "Overview", status = "primary", solidHeader = FALSE,
            tags$div(
              class = "about-box-text",
              style = "padding: 6px 4px 2px 4px;",
              tags$h1("ShineMD", class = "about-shinemd-title"),
              tags$p(
                "A fully automated, pure-R interactive app for molecular dynamics trajectory analysis.",
                style = "font-size: 20px; margin-bottom: 16px;"
              ),
              tags$p(
                "Developed at the Laboratory of Bioactive Peptides, Department of Organic Chemistry, Faculty of Biochemistry and Biological Sciences, National University of the Littoral, Santa Fe, Argentina.",
                style = "margin-bottom: 0; font-size: 15px;"
              )
            )
          ),
          
          box(
            width = 4, title = "Laboratory", status = "primary", solidHeader = FALSE,
            tags$div(
              class = "about-lab-logo",
              style = "text-align: center; padding: 16px 8px 12px 8px;",
              tags$img(
                src = "lpb_logo_unmatted.png",
                style = "max-width: 325px; width: 98%; height: auto; display: inline-block;"
              ),
              tags$div(
                "Laboratory of Bioactive Peptides",
                style = "margin-top: 16px; font-size: 18px; font-weight: 600; line-height: 1.35; color: #eef4ff;"
              )
            )
          )
        ),
        
        fluidRow(
          box(
            width = 6, title = "Contact", status = "primary", solidHeader = FALSE,
            tags$div(
              class = "about-box-text",
              style = "padding: 4px 2px 2px 2px;",
              
              fluidRow(
                column(
                  width = 6,
                  tags$div(
                    style = "padding-right: 14px;",
                    tags$div(
                      "Prof. Álvaro Sebastián Siano",
                      style = "font-size: 18px; font-weight: 600; color: #eef4ff; margin-bottom: 4px;"
                    ),
                    tags$div(
                      "Corresponding author / scientific contact",
                      style = "font-size: 15px; margin-bottom: 6px;",
                      class = "about-subtle"
                    ),
                    tags$div(
                      "Laboratory of Bioactive Peptides",
                      style = "font-size: 14px; margin-bottom: 6px;",
                      class = "about-muted"
                    ),
                    tags$a("asiano@fbcb.unl.edu.ar", href = "mailto:asiano@fbcb.unl.edu.ar")
                  )
                ),
                
                column(
                  width = 6,
                  tags$div(
                    style = "padding-left: 14px;",
                    tags$div(
                      "Dr. Iván Sanchis",
                      style = "font-size: 18px; font-weight: 600; color: #eef4ff; margin-bottom: 4px;"
                    ),
                    tags$div(
                      "App development / repository maintainer",
                      style = "font-size: 15px; margin-bottom: 6px;",
                      class = "about-subtle"
                    ),
                    tags$div(
                      "Laboratory of Bioactive Peptides",
                      style = "font-size: 14px; margin-bottom: 6px;",
                      class = "about-muted"
                    ),
                    tags$a("sanchisivan@fbcb.unl.edu.ar", href = "mailto:sanchisivan@fbcb.unl.edu.ar")
                  )
                )
              ),
              
              tags$hr(style = "border-color: rgba(255,255,255,0.08); margin: 16px 0 10px 0;"),
              
              tags$div(
                style = "font-size: 13px;",
                tags$span(tags$b("Scientific correspondence: "), "Prof. Álvaro Sebastián Siano"),
                tags$br(),
                tags$span(tags$b("App and repository: "), "Dr. Iván Sanchis")
              )
            )
          ),
          
          box(
            width = 6, title = "Citation & Repository", status = "primary", solidHeader = FALSE,
            tags$div(
              class = "about-box-text",
              style = "padding: 4px 2px;",
              tags$p(
                tags$b("Manuscript title:"), tags$br(),
                "ShineMD: an automated pure-R platform with a graphical user interface for integrated analysis of AMBER molecular dynamics trajectories"
              ),
              tags$p(
                tags$b("Preferred citation:"), tags$br(),
                "Sanchis I, Siano AS. ShineMD: an automated pure-R platform with a graphical user interface for integrated analysis of AMBER molecular dynamics trajectories. ChemRxiv preprint, DOI pending."
              ),
              tags$hr(style = "border-color: rgba(255,255,255,0.08);"),
              tags$p(
                tags$b("GitHub repository: "),
                tags$a("github.com/sanchisivan/ShineMD", href = "https://github.com/sanchisivan/ShineMD", target = "_blank")
              )
            )
          )
        ),
        
        fluidRow(
          box(
            width = 12, title = "Institutions", status = "primary", solidHeader = FALSE,
            tags$div(
            style = "
            display: flex;
            align-items: center;
            justify-content: space-between;
            gap: 36px;
            min-height: 150px;
            padding: 12px 18px 8px 12px;
            flex-wrap: nowrap;
            ",
              
              tags$div(
                class = "about-box-text",
                style = "flex: 1 1 620px;",
                tags$div(
                  "Faculty of Biochemistry and Biological Sciences (FBCB) - National University of the Littoral (UNL)",
                  style = "font-size: 17px; font-weight: 600; margin-bottom: 8px;"
                ),
                tags$div(
                  "National Scientific and Technical Research Council (CONICET)",
                  style = "font-size: 17px; font-weight: 600; margin-bottom: 8px;"
                ),
                tags$div(
                  "Santa Fe, Argentina",
                  style = "font-size: 15px;",
                  class = "about-muted"
                )
              ),
              
            tags$div(
              style = "
    flex: 0 0 440px;
    display: flex;
    align-items: center;
    justify-content: center;
    margin-right: 0;
  ",
              tags$div(
                class = "about-inst-panel",
                style = "
      width: 100%;
      display: flex;
      align-items: center;
      justify-content: center;
      gap: 20px;
    ",
                tags$img(
                  src = "Logo-FBCB.png",
                  style = "height: 84px; width: auto; display: block;"
                ),
                tags$div(class = "about-inst-sep"),
                tags$img(
                  src = "Conicet_Logo_sin_letras.png",
                  style = "height: 84px; width: auto; display: block;"
                )
              )
              )
            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  roots <- c(Home = fs::path_home(), shinyFiles::getVolumes()())
  shinyDirChoose(input, "proj_dir_btn", roots = roots, session = session)

  rv <- reactiveValues(
    proj_dir = NULL,
    prmtop = NULL,
    trajs = character(),
    out_dir = NULL,
    log = "",
    data = list(),
    clust = NULL,
    clust_plot_3d = FALSE,
    dimred_plot_3d = FALSE,
    plotly_theme = "dark",
    prm_obj = NULL
  )

  # Plotly theme toggle available in multiple tabs (all buttons control the same theme)
  plotly_theme_btn_ids <- c(
    "toggle_plotly_theme",           # clustering
    "toggle_plotly_theme_rmsd",
    "toggle_plotly_theme_rmsdB",
    "toggle_plotly_theme_rmsf",
    "toggle_plotly_theme_rmsfB",
    "toggle_plotly_theme_rg",
    "toggle_plotly_theme_rgB",
    "toggle_plotly_theme_dimred",
    "toggle_plotly_theme_pcaVar",
    "toggle_plotly_theme_memDist",
    "toggle_plotly_theme_memZ",
    "toggle_plotly_theme_memThick",
    "toggle_plotly_theme_memAPL",
    "toggle_plotly_theme_memEnrich",
    "toggle_plotly_theme_memDensity",
    "toggle_plotly_theme_memOrder",
    "toggle_plotly_theme_intTs",
    "toggle_plotly_theme_intOcc",
    "toggle_plotly_theme_intPair",
    "toggle_plotly_theme_hbTs",
    "toggle_plotly_theme_hbPair")

  sync_plotly_theme_buttons <- function() {
    is_light <- identical(rv$plotly_theme, "light")
    lab <- if (is_light) "Plot theme: Light" else "Plot theme: Dark"
    ic  <- if (is_light) icon("sun") else icon("moon")
    for (btn_id in plotly_theme_btn_ids) {
      updateActionButton(session, btn_id, label = lab, icon = ic)
    }
  }

  toggle_plotly_theme_now <- function() {
    rv$plotly_theme <- if (identical(rv$plotly_theme, "dark")) "light" else "dark"
    sync_plotly_theme_buttons()
  }

  lapply(plotly_theme_btn_ids, function(btn_id) {
    observeEvent(input[[btn_id]], {
      toggle_plotly_theme_now()
    }, ignoreInit = TRUE)
  })

  observe({
    rv$plotly_theme
    sync_plotly_theme_buttons()
  })

  # ───────────────────────────────────────────────────────────────────────────
  # Per-plot info modals (small ⓘ buttons)
  # ───────────────────────────────────────────────────────────────────────────
  plot_info <- list(
    info_rmsdA = list(
      title = "RMSD — Selection A",
      body = tagList(
        tags$p("Root-mean-square deviation (Å) of Selection A vs the reference. Backbone and heavy-atom series are fit and evaluated on matching atom sets, closer to a cpptraj-style RMSD workflow."),
        tags$ul(
          tags$li("Low/stable RMSD → structural stability or convergence."),
          tags$li("Step changes → conformational transitions or binding/insertion events."),
          tags$li("Use together with clustering/PCA to identify states.")
        )
      )
    ),
    info_rmsdB = list(
      title = "RMSD — Selection B",
      body = tagList(
        tags$p("Root-mean-square deviation (Å) of Selection B vs the reference frame. Mirrors the Selection A panel for your second group of atoms (e.g., lipid headgroups, a second peptide, or a receptor domain)."),
        tags$ul(
          tags$li("Useful for comparing stability between two molecular systems or regions."),
          tags$li("Low/stable RMSD → structural stability; step changes → conformational transitions."),
          tags$li("Use together with Selection A RMSD and the interaction panels to relate structural changes to binding events.")
        )
      )
    ),
    info_rmsf = list(
      title = "RMSF (per residue) — Selection A",
      body = tagList(
        tags$p("Root-mean-square fluctuation (Å) per residue (Cα-based) for Selection A."),
        tags$ul(
          tags$li("Highlights flexible loops/termini and rigid cores."),
          tags$li("Good to compare flexibility between conditions/systems.")
        )
      )
    ),
    info_rmsfB = list(
      title = "RMSF (per residue) — Selection B",
      body = tagList(
        tags$p("Root-mean-square fluctuation (Å) per residue (Cα-based) for Selection B."),
        tags$ul(
          tags$li("Measures per-residue flexibility of the second selection (ligand, second domain, etc.)."),
          tags$li("Compare with Selection A RMSF to identify differences in conformational flexibility between the two regions.")
        )
      )
    ),
    info_rgA = list(
      title = "Radius of gyration — Selection A",
      body = tagList(
        tags$p("Radius of gyration (Å): compactness of Selection A over time."),
        tags$ul(
          tags$li("Decrease → compaction; increase → unfolding/extension."),
          tags$li("Often complements RMSD when there are shape changes.")
        )
      )
    ),
    info_rgB = list(
      title = "Radius of gyration — Selection B",
      body = tagList(
        tags$p("Radius of gyration (Å): compactness of Selection B over time. Mirrors the Selection A panel for your second group of atoms."),
        tags$ul(
          tags$li("Decrease → compaction or insertion; increase → unfolding or extension."),
          tags$li("Compare with Selection A Rg and with the interaction distance panels to see whether compaction correlates with binding.")
        )
      )
    ),
    info_dimred = list(
      title = "PCA / Dimensionality reduction",
      body = tagList(
        tags$p("PCA on aligned Cα coordinates (Selection A). Each point is a frame in a low-dimensional space."),
        tags$ul(
          tags$li("Clusters/regions in PCA space often correspond to conformational states."),
          tags$li("Color (time) helps identify transitions and state lifetimes.")
        )
      )
    ),
    info_pcaVar = list(
      title = "Explained variance",
      body = tagList(
        tags$p("Variance captured by each principal component (PC)."),
        tags$ul(
          tags$li("Helps decide how many PCs are needed to describe most motions."),
          tags$li("Cumulative variance close to 1 indicates good coverage.")
        )
      )
    ),
    info_memDist = list(
      title = "Membrane COM distance",
      body = tagList(
        tags$p("Distance (Å) between the target center-of-mass (A or B) and the membrane center-of-mass."),
        tags$ul(
          tags$li("Decreasing distance → approach/adsorption/insertion."),
          tags$li("Interpret together with Δz and enrichment for membrane systems.")
        )
      )
    ),
    info_memZ = list(
      title = "Membrane Δz (COM)",
      body = tagList(
        tags$p("Δz (Å) between target COM and membrane COM along the bilayer normal (typically z)."),
        tags$ul(
          tags$li("Sign indicates which side/leaflet the target is closer to."),
          tags$li("Useful to detect crossing, insertion depth, or leaflet preference.")
        )
      )
    ),
    info_memThick = list(
      title = "Bilayer thickness",
      body = tagList(
        tags$p("Thickness (Å) estimated as mean(z_upper headgroups) − mean(z_lower headgroups) per frame."),
        tags$ul(
          tags$li("Thinning/thickening can indicate peptide-induced deformation."),
          tags$li("Make sure the headgroup atom name matches your topology (e.g., P or P31).")
        )
      )
    ),
    info_memAPL = list(
      title = "Area per lipid (APL)",
      body = tagList(
        tags$p("Area per lipid (Å²) computed from the box XY area divided by lipids per leaflet."),
        tags$ul(
          tags$li("Higher APL usually means looser packing / more fluid membrane."),
          tags$li("Lower APL usually means tighter packing / more ordered membrane.")
        )
      )
    ),
    info_memEnrich = list(
      title = "Lipid enrichment around target",
      body = tagList(
        tags$p("Counts how many headgroup atoms are within a cutoff (Å) of the target (A or B), per frame."),
        tags$ul(
          tags$li("A contact proxy: higher values → more frequent/stronger surface association."),
          tags$li("If 'separate by lipid type' is enabled, compares lipid preferences (e.g., PA vs PE)."),
          tags$li("Cutoff sensitivity: try 8–12 Å for surface contacts.")
        )
      )
    ),
    info_memDensity = list(
      title = "Membrane density profiles",
      body = tagList(
        tags$p("Density profiles along the membrane normal (z) for headgroups, tails, water, ions, and an optional target selection."),
        tags$ul(
          tags$li("Useful to locate where each component sits across the bilayer thickness."),
          tags$li("Headgroups usually peak near the two membrane interfaces, tails in the hydrophobic core, and water/ions outside the bilayer."),
          tags$li("If membrane centering is enabled, each frame is recentered before profile accumulation to reduce drift artifacts."),
          tags$li("Make sure the headgroup atom names, tail regex, water residue names, and ion residue names match your topology.")
        )
      )
    ),
    info_memOrder = list(
      title = "Lipid tail order profile |S|",
      body = tagList(
        tags$p("Pure-R segmental approximation of lipid tail order along the acyl chains."),
        tags$ul(
          tags$li("Higher |S| values indicate more ordered / aligned chains."),
          tags$li("Lower |S| values indicate more disordered / flexible tails."),
          tags$li("This is useful to compare membrane packing between conditions, regions, or lipid types."),
          tags$li("The quality of the result depends on atom naming consistency for the sn1/sn2 tail regex patterns."),
          tags$li("Use the standard mode when both acyl chains are stored inside the same lipid residue."),
          tags$li("Use the split-residue mode when one lipid is represented by separate residues for chain 1 and chain 2. Example: PA + PE + OL, where PA is one chain, PE is the linker/headgroup, and OL is the other chain."),
          tags$li("In mixed membranes, sterols such as ERG/CHL should usually be excluded, because their carbon names can accidentally match tail regex patterns without representing phospholipid sn1/sn2 chains."),
          tags$li("If one lipid type shows only one point, the regex likely matched only two atoms for that chain in that residue family.")
        )
      )
    ),
    info_intTs = list(
      title = "Selection A / Selection B interaction time series",
      body = tagList(
        tags$p("Tracks simple interaction metrics between the two selections over time."),
        tags$ul(
          tags$li("Minimum distance: closest heavy-atom separation between the two selections."),
          tags$li("COM distance: center-of-mass distance between the selections."),
          tags$li("Atom contacts: number of atom pairs within the chosen cutoff."),
          tags$li("Stable low distances / high contact counts usually indicate a persistent interface.")
        )
      )
    ),
    info_intOcc = list(
      title = "Residue contact occupancy summary",
      body = tagList(
        tags$p("Shows which residues on one side remain in contact with the opposite selection during the trajectory."),
        tags$ul(
          tags$li("Occupancy (%) = fraction of analyzed frames where the residue makes at least one contact."),
          tags$li("High-occupancy residues usually define the most persistent binding/interface patch."),
          tags$li("This is residue-level occupancy, not a guarantee that the exact same atom-atom contact is preserved every frame.")
        )
      )
    ),
    info_intPair = list(
      title = "A/B residue contact map",
      body = tagList(
        tags$p("Heatmap of residue-pair contact occupancy between Selection A and Selection B."),
        tags$ul(
          tags$li("Each cell shows the percentage of frames where that specific residue pair makes at least one contact within the cutoff."),
          tags$li("This is more specific than the side-wise occupancy plot because it resolves which protein residues contact which peptide/ligand residues."),
          tags$li("Use the top-residue filters to focus on the dominant interaction patch and avoid overcrowded maps.")
        )
      )
    ),
    info_hbTs = list(
      title = "A / B hydrogen-bond summary (proxy mode)",
      body = tagList(
        tags$p("Distance-based donor/acceptor proxy for hydrogen bonds between Selection A and Selection B."),
        tags$ul(
          tags$li("Donors are atoms whose names start with N, O, or S."),
          tags$li("Acceptors are O/S by default; you can optionally allow N as a permissive proxy."),
          tags$li("This panel is intentionally labeled as a proxy because it does not enforce a hydrogen-bond angle criterion."),
          tags$li("Use it to find persistent polar contacts quickly in pure R, then confirm key pairs in 3D if needed.")
        )
      )
    ),
    info_hbPair = list(
      title = "Persistent A / B hydrogen-bond pairs",
      body = tagList(
        tags$p("Ranks donor/acceptor atom pairs by the fraction of frames where they satisfy the proxy hydrogen-bond criterion."),
        tags$ul(
          tags$li("High occupancy means a persistent polar interaction across the trajectory."),
          tags$li("Direction matters: A donor → B acceptor is not the same as B donor → A acceptor."),
          tags$li("Pairs here are atom-specific, so they are more detailed than residue-level contact occupancy.")
        )
      )
    ),
    info_clusterPlots = list(
      title = "Cluster plots",
      body = tagList(
        tags$p("Visual summaries of RMSD-based clustering."),
        tags$ul(
          tags$li("Distribution: embedding of frames colored by cluster."),
          tags$li("Clusters vs time: when each state appears."),
          tags$li("Population: how many frames in each cluster."),
          tags$li("Dendrogram: hierarchical relationship (only for hclust).")
        )
      )
    ),
    info_clusterSummary = list(
      title = "Cluster summary",
      body = tagList(
        tags$p("Per-cluster statistics and representative frames (medoids and centroid-nearest frames)."),
        tags$ul(
          tags$li("Use medoids/centroids to export representative PDBs."),
          tags$li("Population and within-cluster RMSD help assess cluster quality.")
        )
      )
    ),
    info_clustDist = list(
      title = "Cluster distribution",
      body = tagList(
        tags$p("2D or 3D embedding of all analyzed frames based on their pairwise RMSD distances (MDS; PCA used as fallback). Each point is a trajectory frame, colored by cluster assignment."),
        tags$ul(
          tags$li("Well-separated clouds → distinct conformational states with low inter-cluster RMSD."),
          tags$li("Overlapping regions → gradual conformational transitions or closely related states."),
          tags$li("Use the 'Toggle 3D distribution' button to switch between 2D and interactive 3D views.")
        )
      )
    ),
    info_clustTime = list(
      title = "Clusters vs time",
      body = tagList(
        tags$p("Scatter plot of cluster assignment over the simulation time. Each point is a frame colored by its cluster."),
        tags$ul(
          tags$li("Long stretches of a single cluster → stable conformational state or slow dynamics."),
          tags$li("Rapid alternation between clusters → fast exchange or sampling of a transition region."),
          tags$li("Use this plot together with the RMSD time series to correlate structural transitions with specific time windows.")
        )
      )
    ),
    info_clustPop = list(
      title = "Cluster population",
      body = tagList(
        tags$p("Bar chart showing the number of frames (and percentage) assigned to each cluster."),
        tags$ul(
          tags$li("Dominant clusters (high %) usually represent the most thermodynamically stable states."),
          tags$li("Very small clusters may reflect transient states or outlier frames — inspect them in a 3D viewer."),
          tags$li("Population is a proxy for relative free energy: more populated states have lower free energy.")
        )
      )
    ),
    info_clustQual = list(
      title = "Cluster quality (silhouette)",
      body = tagList(
        tags$p("Average silhouette width by cluster. This helps judge how well-separated each state is in the RMSD distance space."),
        tags$ul(
          tags$li("Values near 1 mean frames are much closer to their own cluster than to others."),
          tags$li("Values near 0 suggest overlap / weak separation."),
          tags$li("Negative values suggest likely misassignment or too many clusters.")
        )
      )
    ),
    info_clustDend = list(
      title = "Dendrogram",
      body = tagList(
        tags$p("Hierarchical clustering tree built from the pairwise RMSD matrix. Only available when using the hclust method."),
        tags$ul(
          tags$li("Branch height represents the RMSD dissimilarity at which two groups merge."),
          tags$li("Cut the tree at a given height (analogous to choosing k) to define a desired number of clusters."),
          tags$li("Tight subtrees with low merge heights → conformationally similar frame groups.")
        )
      )
    ),
    info_rmsdHeatmap = list(
      title = "Pairwise RMSD matrix",
      body = tagList(
        tags$p("Frame-vs-frame RMSD matrix visualized as a heatmap. Each pixel (i,j) shows the RMSD between frame i and frame j."),
        tags$ul(
          tags$li("Diagonal blocks of low RMSD (blue) indicate conformational substates."),
          tags$li("Off-diagonal transitions (yellow/red) indicate conformational changes."),
          tags$li("This matrix is computed during clustering — run clustering first to populate this panel.")
        )
      )
    ),
    info_fel = list(
      title = "Free energy landscape (FEL)",
      body = tagList(
        tags$p("Gibbs free energy surface projected onto the first two principal components:"),
        tags$p("ΔG(PC1, PC2) = −kT · ln P(PC1, PC2)"),
        tags$ul(
          tags$li("Darker basins correspond to the most populated (lowest free energy) conformational states."),
          tags$li("Barriers between basins indicate conformational transitions."),
          tags$li("The energy is relative: the global minimum is set to 0 kcal/mol."),
          tags$li("Temperature is used for the kT scaling (default 300 K).")
        )
      )
    ),
    info_rmsdDist = list(
      title = "RMSD distribution",
      body = tagList(
        tags$p("Histogram of RMSD values across all analyzed frames. A unimodal distribution suggests a single conformational basin; multimodal distributions indicate distinct substates or conformational transitions during the simulation.")
      )
    )
  )

  lapply(names(plot_info), function(btn_id) {
    observeEvent(input[[btn_id]], {
      info <- plot_info[[btn_id]]
      showModal(modalDialog(
        class = "themed-modal",
        title = tagList(icon("info-circle"), info$title),
        easyClose = TRUE,
        footer = modalButton("Close"),
        info$body
      ))
    }, ignoreInit = TRUE)
  })


  show_themed_completion <- function(title_text, hero_title, subtitle, path_value = NULL,
                                     path_label = "Location", steps = character()) {
    step_nodes <- lapply(seq_along(steps), function(i) {
      tags$div(
        class = "completion-step",
        tags$span(class = "completion-step-index", as.character(i)),
        tags$div(tags$div(class = "completion-step-text", steps[[i]]))
      )
    })

    modal_bits <- list(
      class = "themed-modal",
      title = tagList(icon("check-circle"), title_text),
      easyClose = TRUE,
      footer = modalButton("Close"),
      tags$div(
        class = "completion-hero",
        tags$div(class = "completion-hero-title", hero_title),
        tags$div(class = "completion-hero-subtitle", subtitle)
      )
    )

    if (!is.null(path_value) && nzchar(path_value)) {
      modal_bits <- c(modal_bits, list(
        tags$div(class = "completion-chip-label", path_label),
        tags$div(class = "completion-path-chip", path_value)
      ))
    }

    if (length(steps) > 0) {
      modal_bits <- c(modal_bits, list(
        tags$div(class = "completion-steps", step_nodes)
      ))
    }

    showModal(do.call(modalDialog, modal_bits))
  }

  show_waiter_modal <- function(title_text = "Please wait", subtitle = "Processing data...") {
    showModal(modalDialog(
      class = "themed-modal waiter-modal",
      title = NULL,
      easyClose = FALSE,
      footer = NULL,
      tags$div(
        class = "waiter-shell",
        tags$div(class = "waiter-spinner", icon("spinner", class = "fa-spin")),
        tags$div(class = "waiter-title", title_text),
        tags$div(class = "waiter-subtitle", subtitle)
      )
    ))
  }

  with_waiter_modal <- function(title_text, subtitle = "Processing data...", expr) {
    caller_env <- parent.frame()
    show_waiter_modal(title_text, subtitle)
    on.exit(removeModal(), add = TRUE)
    eval(substitute(expr), envir = caller_env)
  }

  progress_with_waiter <- function(title_text, subtitle = "Processing data...",
                                   message = title_text, value = 0, expr) {
    caller_env <- parent.frame()
    show_waiter_modal(title_text, subtitle)
    on.exit(removeModal(), add = TRUE)
    withProgress(message = message, value = value, {
      eval(substitute(expr), envir = caller_env)
    })
  }

  show_traj_subset_modal <- function(modal_key, title_text, trajs, action_label = "Compute") {
    labs <- basename(trajs)
    showModal(modalDialog(
      class = "themed-modal",
      title = tagList(icon("stream"), title_text),
      easyClose = TRUE,
      size = "m",
      footer = tagList(
        modalButton("Cancel"),
        actionButton(paste0("confirm_", modal_key), action_label, class = "btn-primary")
      ),
      tags$p("Choose which loaded trajectory files will be used for this calculation. All files are selected by default; deselect any segments you want to exclude."),
      selectizeInput(
        inputId = paste0(modal_key, "_trajs"),
        label = "Trajectory files",
        choices = labs,
        selected = labs,
        multiple = TRUE,
        options = list(placeholder = 'Select trajectory files', plugins = list('remove_button'))
      ),
      tags$small(style = "color: var(--text-muted);",
                 "The concatenation order follows the Trajectory order defined in the Project tab. Click the \u00d7 on any tag to deselect a file.")
    ))
  }

  get_modal_trajs <- function(modal_key, trajs_all) {
    chosen <- isolate(input[[paste0(modal_key, "_trajs")]]) %||% character()
    if (!length(chosen)) {
      showNotification("Select at least one trajectory file.", type = "warning")
      return(NULL)
    }
    traj_map <- setNames(trajs_all, basename(trajs_all))
    chosen <- intersect(chosen, names(traj_map))
    if (!length(chosen)) {
      showNotification("The selected trajectory files are not available anymore. Reopen the compute window and try again.", type = "error", duration = 8)
      return(NULL)
    }
    unname(traj_map[chosen])
  }

  observeEvent(input$proj_dir_btn, {
    sel <- parseDirPath(roots, input$proj_dir_btn)
    if (length(sel) && nzchar(sel)) {
      updateTextInput(session, "proj_dir", value = sel)
    }
  })

  observeEvent(input$proj_dir, {
    p <- trimws(input$proj_dir)
    if (!nzchar(p) || !dir_exists(p)) {
      rv$proj_dir <- NULL
      rv$prmtop <- NULL
      rv$trajs <- character()
      return()
    }
    with_waiter_modal(
      "Reading project files",
      "Detecting the topology, indexing trajectory segments, and caching the topology for faster downstream analyses.",
      {
        rv$proj_dir <- p
        inp <- detect_inputs(p)
        rv$prmtop <- if (!is.na(inp$prmtop)) inp$prmtop else NULL
        rv$trajs <- inp$trajs

        # Cache parsed topology for downstream (faster advanced analyses)
        rv$prm_obj <- NULL
        if (!is.null(rv$prmtop) && file.exists(rv$prmtop)) {
          rv$prm_obj <- tryCatch(read.prmtop(rv$prmtop), error = function(e) NULL)
        }

        # prefill order box with detected order (safe default)
        if (length(rv$trajs) > 0) {
          updateTextAreaInput(session, "traj_order", value = paste(basename(rv$trajs), collapse="\n"))
        }
      }
    )
  })

  output$proj_status <- renderUI({
    if (is.null(rv$proj_dir)) {
      tags$span(class="status-badge status-empty", "No project selected")
    } else {
      tags$span(class="status-badge status-loaded", paste0("✓ ", rv$proj_dir))
    }
  })

  output$detected_files <- renderUI({
    if (is.null(rv$proj_dir)) {
      return(tags$p(style="color: var(--text-muted);", "Select a folder to auto-detect prmtop + .nc segments."))
    }
    tagList(
      tags$p(strong("prmtop:"), rv$prmtop %||% "NOT FOUND"),
      tags$p(strong("segments:"), length(rv$trajs)),
      tags$pre(style="color:#8899aa;", paste(basename(rv$trajs), collapse="\n"))
    )
  })

  ordered_trajs <- reactive({
    req(length(rv$trajs) >= 1)
    det <- rv$trajs
    det_names <- basename(det)

    lines <- trimws(unlist(strsplit(input$traj_order %||% "", "\n", fixed = TRUE)))
    lines <- lines[nzchar(lines)]

    if (length(lines) == 0) return(det)

    missing <- setdiff(lines, det_names)
    if (length(missing) > 0) {
      stop(paste0("These files are not in the folder: ", paste(missing, collapse=", ")))
    }

    det[match(lines, det_names)]
  })

  observeEvent(input$run_all, {
    req(rv$proj_dir, rv$prmtop, length(rv$trajs) >= 1)
    rv$log <- "Running analysis...\n"

    # validate order
    trajs <- NULL
    tryCatch({
      trajs <- ordered_trajs()
    }, error=function(e) {
      rv$log <- paste0(rv$log, "ERROR: ", e$message, "\n")
    })
    if (is.null(trajs)) return()

    pep_resno <- parse_resno(input$pep_resno)
    lipo_resno <- parse_resno(input$lipo_resno)
    mem_resid <- parse_csv_tokens(input$mem_resid)
    exclude_resid <- parse_csv_tokens(input$exclude_resid)

    align_elety <- if (identical(input$align_mode, "ca")) "CA" else c("N","CA","C","O")
    mem_target <- input$mem_target %||% "B"

    first <- if (is.na(input$first)) NULL else as.integer(input$first)
    last  <- if (is.na(input$last)) NULL else as.integer(input$last)

    # Planned workload info (helps users estimate waiting time)
    plan_df <- build_frame_index(trajs, dt_ps = input$dt_ps, first = first, last = last,
                               stride = as.integer(input$stride), combine = isTRUE(input$combine_trajs))
    rv$log <- paste0(rv$log, sprintf("Planned workload: %d segments, %d frames after stride.\n",
                                    length(unique(plan_df$segment)), nrow(plan_df)))
    seg_tab <- sort(table(plan_df$segment), decreasing = TRUE)
    if (length(seg_tab) > 0) {
      show <- head(seg_tab, 8)
      rv$log <- paste0(rv$log, "Frames per segment (after stride):\n",
                       paste0(" - ", names(show), ": ", as.integer(show), "\n", collapse = ""),
                       if (length(seg_tab) > 8) " - ...\n" else "")
    }

    tryCatch({
      progress_with_waiter("Running basic analysis",
                           "Reading trajectories, aligning frames, and computing baseline structural descriptors.",
                           message = "Computing metrics...", value = 0, {
        t0 <- Sys.time()
        logf <- function(msg) { rv$log <- paste0(rv$log, msg, "\n") }
        progress_cb <- function(value, detail = NULL) { setProgress(value = value, detail = detail) }
        incProgress(0.02, detail = "Initializing")
        ref_coords_path <- NULL
        ref_coords_name <- NULL
        if (isTRUE(input$use_rmsd_ref_file) && !is.null(input$rmsd_ref_file$datapath) && nzchar(input$rmsd_ref_file$datapath)) {
          ref_coords_path <- input$rmsd_ref_file$datapath
          ref_coords_name <- input$rmsd_ref_file$name
          logf(sprintf("Using external reference coordinates: %s", ref_coords_name %||% basename(ref_coords_path)))
        }
        rv$out_dir <- run_basic_analysis(
          project_dir = rv$proj_dir,
          trajs_ordered = trajs,
          combine = isTRUE(input$combine_trajs),
          dt_ps = input$dt_ps,
          first = first,
          last = last,
          stride = as.integer(input$stride),
          pep_resno = pep_resno,
          lipo_resno = lipo_resno,
          membrane_resid = mem_resid,
          exclude_resid = exclude_resid,
          mem_target = mem_target,
          align_elety = align_elety,
          rmsd_ref_mode = input$rmsd_ref_mode %||% "global_first",
          rmsd_ref_path = ref_coords_path,
          rmsd_ref_name = ref_coords_name,
          progress = progress_cb,
          logf = logf,
          verbose = TRUE
        )
        setProgress(1, detail = "Finalizing")
        dt_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        rv$log <- paste0(rv$log, sprintf("Basic analysis runtime: %.1f s\n", dt_sec))
      })
      rv$log <- paste0(rv$log, "Done. Results in: ", rv$out_dir, "\n")

      rv$data$rmsd_pep  <- read_if_exists(file.path(rv$out_dir, "rmsd_regionA_all.csv")) %||%
        read_if_exists(file.path(rv$out_dir, "rmsd_peptide_all.csv"))
      rv$data$rmsd_lipo <- read_if_exists(file.path(rv$out_dir, "rmsd_regionB_all.csv")) %||%
        read_if_exists(file.path(rv$out_dir, "rmsd_lipopeptide_all.csv"))

      rv$data$rg_pep    <- read_if_exists(file.path(rv$out_dir, "rg_regionA_all.csv")) %||%
        read_if_exists(file.path(rv$out_dir, "rg_peptide_all.csv"))
      rv$data$rg_lipo   <- read_if_exists(file.path(rv$out_dir, "rg_regionB_all.csv")) %||%
        read_if_exists(file.path(rv$out_dir, "rg_lipopeptide_all.csv"))

      rv$data$dist      <- read_if_exists(file.path(rv$out_dir, "region_mem_com_dist_all.csv")) %||%
        read_if_exists(file.path(rv$out_dir, "pep_mem_com_dist_all.csv"))
      rv$data$z         <- read_if_exists(file.path(rv$out_dir, "region_mem_z_all.csv")) %||%
        read_if_exists(file.path(rv$out_dir, "pep_mem_z_all.csv"))

      rmsf_path <- file.path(rv$out_dir, "rmsf_regionA_byres.csv")
      if (!file.exists(rmsf_path)) rmsf_path <- file.path(rv$out_dir, "rmsf_peptide_byres.csv")
      rv$data$rmsf   <- if (file.exists(rmsf_path)) read.csv(rmsf_path) else NULL
      rmsf_B_path <- file.path(rv$out_dir, "rmsf_regionB_byres.csv")
      rv$data$rmsf_B <- if (file.exists(rmsf_B_path)) read.csv(rmsf_B_path) else NULL
      rv$data$pca_scores <- read_if_exists(file.path(rv$out_dir, "pca_regionA_scores.csv"))
      rv$data$pca_var <- normalize_pca_var_df(read_if_exists(file.path(rv$out_dir, "pca_regionA_variance.csv")))
      rv$data$mem_density <- read_if_exists(file.path(rv$out_dir, "membrane_density_profile.csv"))
      rv$data$mem_order <- read_if_exists(file.path(rv$out_dir, "lipid_tail_order_profile.csv"))
      rv$data$interact_ts <- read_if_exists(file.path(rv$out_dir, "ab_interaction_timeseries.csv"))
      rv$data$interact_occ <- read_if_exists(file.path(rv$out_dir, "ab_interaction_occupancy.csv"))
      rv$data$interact_pairs <- read_if_exists(file.path(rv$out_dir, "ab_interaction_pair_occupancy.csv"))
      rv$data$hbond_ts <- read_if_exists(file.path(rv$out_dir, "ab_hbond_proxy_timeseries.csv"))
      rv$data$hbond_pairs <- read_if_exists(file.path(rv$out_dir, "ab_hbond_proxy_pairs.csv"))

      if (!is.null(rv$log) && nzchar(rv$log)) {
        writeLines(rv$log, con = file.path(rv$out_dir, "analysis_log.txt"))
      }

      # Emphasize completion + guide next steps
      if (!is.null(rv$out_dir) && dir.exists(rv$out_dir)) {
        rv$log <- paste0(rv$log,
                         "\n✅ Processing finished successfully.\n\n",
                         "Next steps:\n",
                         "  1) Open the analysis tabs (RMSD / RMSF / Rg / Membrane) to inspect results.\n",
                         "  2) Use the per-plot *Plot styling + export* panel next to each plot to customize and download that specific plot.\n",
                         "  3) If you need representative structures or a specific frame as PDB, go to the Clustering tab and use the export buttons.\n")

        show_themed_completion(
          title_text = "Processing finished",
          hero_title = "Basic analysis completed successfully",
          subtitle = "Your metrics were generated and the results are ready to inspect in the analysis tabs.",
          path_value = rv$out_dir,
          path_label = "Results folder",
          steps = c(
            "Open RMSD, RMSF, Radius of gyration, PCA / DimRed, and Membrane metrics to review the processed trajectory outputs.",
            "Use the plot styling + export panel beside each figure to tune appearance and save publication-ready files.",
            "Go to Clustering to review distribution, cluster populations, and export medoids/centroids or a specific frame as PDB."
          )
        )
      }

    }, error=function(e) {
      rv$log <- paste0(rv$log, "ERROR: ", e$message, "\n")
    })
  })


  # ───────────────────────────────────────────────────────────────────────────
  # Clustering (RMSD-based structural clustering)
  # ───────────────────────────────────────────────────────────────────────────
  observeEvent(input$clust_ref_quick, {
    if (!is.null(input$clust_ref_quick) && nzchar(as.character(input$clust_ref_quick))) {
      updateNumericInput(session, "clust_ref_time", value = as.numeric(input$clust_ref_quick))
    }
  })

  observeEvent(input$toggle_clust_3d, {
    rv$clust_plot_3d <- !isTRUE(rv$clust_plot_3d)
  })

  observeEvent(input$toggle_dimred_3d, {
    rv$dimred_plot_3d <- !isTRUE(rv$dimred_plot_3d)
  })

  # Advanced membrane computations (optional, on-demand)
  observeEvent(input$compute_thickness, {
    req(rv$proj_dir, rv$prmtop, length(rv$trajs) >= 1)
    trajs <- NULL
    tryCatch({ trajs <- ordered_trajs() }, error = function(e) {
      showNotification(paste0("Trajectory order error: ", e$message), type = "error")
    })
    if (is.null(trajs)) return()
    show_traj_subset_modal("thickness", "Compute bilayer thickness — select trajectory files", trajs)
  })

  observeEvent(input$confirm_thickness, {
    req(rv$proj_dir, rv$prmtop, length(rv$trajs) >= 1)
    mem_resid <- parse_csv_tokens(input$mem_resid)
    first <- if (is.na(input$first)) NULL else as.integer(input$first)
    last  <- if (is.na(input$last)) NULL else as.integer(input$last)
    stride <- as.integer(input$stride %||% 1)
    trajs_all <- NULL
    tryCatch({ trajs_all <- ordered_trajs() }, error = function(e) {
      showNotification(paste0("Trajectory order error: ", e$message), type = "error")
    })
    if (is.null(trajs_all)) return()
    trajs <- get_modal_trajs("thickness", trajs_all)
    if (is.null(trajs)) return()
    removeModal()

    rv$log <- paste0(rv$log, "Running bilayer thickness...\n")
    tryCatch({
      progress_with_waiter("Computing bilayer thickness",
                           "Reading the selected trajectories and estimating the bilayer thickness time series.",
                           message = "Computing bilayer thickness...", value = 0, {
        df <- compute_bilayer_thickness(
          trajs_ordered = trajs, prm = (rv$prm_obj %||% read.prmtop(rv$prmtop)), membrane_resid = mem_resid,
          head_atoms = input$thick_head_atoms,
          dt_ps = input$dt_ps, first = first, last = last, stride = stride,
          combine = isTRUE(input$combine_trajs)
        )
        rv$data$thickness <- df
        if (!is.null(rv$out_dir) && dir.exists(rv$out_dir)) {
          write.csv(df, file.path(rv$out_dir, "membrane_thickness.csv"), row.names = FALSE)
        }
        incProgress(1)
      })
      show_themed_completion(
        title_text = "Membrane analysis finished",
        hero_title = "Bilayer thickness computed",
        subtitle = "The bilayer thickness time series is ready. Results were also exported as a CSV file.",
        steps = c(
          "Go to the Membrane systems tab and inspect the Bilayer thickness plot.",
          "Use the Plot styling + export panel to customize the figure and save it as PNG/SVG.",
          "Download the CSV (membrane_thickness.csv) from the export button for further analysis in R or Python."
        )
      )
    }, error = function(e) {
      showNotification(paste0("Thickness failed: ", e$message), type = "error", duration = 10)
      rv$log <- paste0(rv$log, "ERROR (thickness): ", e$message, "\n")
    })
  })

  observeEvent(input$compute_apl, {
    req(rv$proj_dir, rv$prmtop, length(rv$trajs) >= 1)
    trajs <- NULL
    tryCatch({ trajs <- ordered_trajs() }, error = function(e) {
      showNotification(paste0("Trajectory order error: ", e$message), type = "error")
    })
    if (is.null(trajs)) return()
    show_traj_subset_modal("apl", "Compute area per lipid — select trajectory files", trajs)
  })

  observeEvent(input$confirm_apl, {
    req(rv$proj_dir, rv$prmtop, length(rv$trajs) >= 1)
    mem_resid <- parse_csv_tokens(input$mem_resid)
    first <- if (is.na(input$first)) NULL else as.integer(input$first)
    last  <- if (is.na(input$last)) NULL else as.integer(input$last)
    stride <- as.integer(input$stride %||% 1)
    trajs_all <- NULL
    tryCatch({ trajs_all <- ordered_trajs() }, error = function(e) {
      showNotification(paste0("Trajectory order error: ", e$message), type = "error")
    })
    if (is.null(trajs_all)) return()
    trajs <- get_modal_trajs("apl", trajs_all)
    if (is.null(trajs)) return()
    removeModal()

    n_leaf <- input$apl_nlip_leaflet
    if (is.na(n_leaf)) n_leaf <- NA_integer_ else n_leaf <- as.integer(n_leaf)

    rv$log <- paste0(rv$log, "Running area per lipid...\n")
    tryCatch({
      progress_with_waiter("Computing area per lipid",
                           "Reading the selected trajectories and computing leaflet-level area-per-lipid values.",
                           message = "Computing APL...", value = 0, {
        res <- compute_area_per_lipid(
          trajs_ordered = trajs, prm = (rv$prm_obj %||% read.prmtop(rv$prmtop)), membrane_resid = mem_resid,
          n_lip_leaflet = n_leaf,
          dt_ps = input$dt_ps, first = first, last = last, stride = stride,
          combine = isTRUE(input$combine_trajs)
        )
        rv$data$apl <- res$df
        rv$data$apl_meta <- list(n_total = res$n_total, n_leaflet = res$n_leaflet)
        if (!is.null(rv$out_dir) && dir.exists(rv$out_dir)) {
          write.csv(res$df, file.path(rv$out_dir, "membrane_apl.csv"), row.names = FALSE)
        }
        incProgress(1)
      })
      show_themed_completion(
        title_text = "Membrane analysis finished",
        hero_title = "Area per lipid (APL) computed",
        subtitle = "The APL time series is ready. Results were also exported as a CSV file.",
        steps = c(
          "Go to the Membrane systems tab and inspect the Area per lipid plot.",
          "Compare APL values with reference lipid bilayer benchmarks (e.g., ~68 Å² for POPC at 300 K).",
          "Download the CSV (membrane_apl.csv) from the export button for further statistical analysis."
        )
      )
    }, error = function(e) {
      showNotification(paste0("APL failed: ", e$message), type = "error", duration = 10)
      rv$log <- paste0(rv$log, "ERROR (APL): ", e$message, "\n")
    })
  })

  observeEvent(input$compute_enrich, {
    req(rv$proj_dir, rv$prmtop, length(rv$trajs) >= 1)
    trajs <- NULL
    tryCatch({ trajs <- ordered_trajs() }, error = function(e) {
      showNotification(paste0("Trajectory order error: ", e$message), type = "error")
    })
    if (is.null(trajs)) return()
    show_traj_subset_modal("enrich", "Compute lipid enrichment — select trajectory files", trajs)
  })

  observeEvent(input$confirm_enrich, {
    req(rv$proj_dir, rv$prmtop, length(rv$trajs) >= 1)
    mem_resid <- parse_csv_tokens(input$mem_resid)
    exclude_resid <- parse_csv_tokens(input$exclude_resid)
    pep_resno <- parse_resno(input$pep_resno)
    lipo_resno <- parse_resno(input$lipo_resno)

    align_elety <- if (identical(input$align_mode, "ca")) "CA" else c("N","CA","C","O")
    first <- if (is.na(input$first)) NULL else as.integer(input$first)
    last  <- if (is.na(input$last)) NULL else as.integer(input$last)
    stride <- as.integer(input$stride %||% 1)

    trajs_all <- NULL
    tryCatch({ trajs_all <- ordered_trajs() }, error = function(e) {
      showNotification(paste0("Trajectory order error: ", e$message), type = "error")
    })
    if (is.null(trajs_all)) return()
    trajs <- get_modal_trajs("enrich", trajs_all)
    if (is.null(trajs)) return()
    removeModal()

    rv$log <- paste0(rv$log, "Running lipid enrichment...\n")
    tryCatch({
      progress_with_waiter("Computing lipid enrichment",
                           "Scanning membrane contacts around the selected target and aggregating enrichment statistics.",
                           message = "Computing lipid enrichment...", value = 0, {
        sels <- make_selections(
          prm = (rv$prm_obj %||% read.prmtop(rv$prmtop)),
          pep_resno = pep_resno,
          lipo_resno = lipo_resno,
          membrane_resid = mem_resid,
          exclude_resid = exclude_resid,
          align_elety = align_elety
        )

        res <- compute_lipid_enrichment(
          trajs_ordered = trajs, prm = (rv$prm_obj %||% read.prmtop(rv$prmtop)), sels = sels, membrane_resid = mem_resid,
          target = input$enrich_target,
          cutoff_A = as.numeric(input$enrich_cutoff),
          head_atom = input$enrich_head_atom,
          by_type = isTRUE(input$enrich_by_type),
          dt_ps = input$dt_ps, first = first, last = last, stride = stride,
          combine = isTRUE(input$combine_trajs)
        )
        rv$data$enrich <- res$df
        rv$data$enrich_summary <- res$summary
        if (!is.null(rv$out_dir) && dir.exists(rv$out_dir)) {
          write.csv(res$df, file.path(rv$out_dir, "lipid_enrichment.csv"), row.names = FALSE)
          write.csv(res$summary, file.path(rv$out_dir, "lipid_enrichment_summary.csv"), row.names = FALSE)
        }
        incProgress(1)
      })
      show_themed_completion(
        title_text = "Membrane analysis finished",
        hero_title = "Lipid enrichment around target computed",
        subtitle = "Headgroup contact counts per frame are ready. Results were also exported as CSV files.",
        steps = c(
          "Go to the Membrane systems tab and inspect the Lipid enrichment plot.",
          "If 'separate by lipid type' was enabled, compare the enrichment curves per lipid species to identify selectivity.",
          "Download the CSV files (lipid_enrichment.csv / lipid_enrichment_summary.csv) for further analysis."
        )
      )
    }, error = function(e) {
      showNotification(paste0("Enrichment failed: ", e$message), type = "error", duration = 10)
      rv$log <- paste0(rv$log, "ERROR (enrichment): ", e$message, "\n")
    })
  })

  observeEvent(input$compute_density, {
    req(rv$proj_dir, rv$prmtop, length(rv$trajs) >= 1)
    trajs <- NULL
    tryCatch({ trajs <- ordered_trajs() }, error = function(e) {
      showNotification(paste0("Trajectory order error: ", e$message), type = "error")
    })
    if (is.null(trajs)) return()
    show_traj_subset_modal("density", "Compute membrane density profiles — select trajectory files", trajs)
  })

  observeEvent(input$confirm_density, {
    req(rv$proj_dir, rv$prmtop, length(rv$trajs) >= 1)
    mem_resid <- parse_csv_tokens(input$mem_resid)
    exclude_resid <- parse_csv_tokens(input$exclude_resid)
    pep_resno <- parse_resno(input$pep_resno)
    lipo_resno <- parse_resno(input$lipo_resno)
    water_resid <- parse_csv_tokens(input$water_resid)
    ion_resid <- parse_csv_tokens(input$ion_resid)
    align_elety <- if (identical(input$align_mode, "ca")) "CA" else c("N","CA","C","O")
    first <- if (is.na(input$first)) NULL else as.integer(input$first)
    last  <- if (is.na(input$last)) NULL else as.integer(input$last)
    stride <- as.integer(input$stride %||% 1)

    trajs_all <- NULL
    tryCatch({ trajs_all <- ordered_trajs() }, error = function(e) {
      showNotification(paste0("Trajectory order error: ", e$message), type = "error")
    })
    if (is.null(trajs_all)) return()
    trajs <- get_modal_trajs("density", trajs_all)
    if (is.null(trajs)) return()
    removeModal()

    rv$log <- paste0(rv$log, "Running membrane density profiles...\n")
    tryCatch({
      progress_with_waiter("Computing membrane density profiles",
                           "Reading the selected trajectories and building membrane-centered density histograms.",
                           message = "Computing membrane density profiles...", value = 0, {
        prm_use <- (rv$prm_obj %||% read.prmtop(rv$prmtop))
        sels <- make_selections(
          prm = prm_use,
          pep_resno = pep_resno,
          lipo_resno = lipo_resno,
          membrane_resid = mem_resid,
          exclude_resid = exclude_resid,
          align_elety = align_elety
        )
        df <- compute_density_profiles(
          trajs_ordered = trajs, prm = prm_use, sels = sels, membrane_resid = mem_resid,
          head_atoms = input$dens_head_atoms,
          tail_regex = input$dens_tail_regex,
          water_resid = water_resid,
          ion_resid = ion_resid,
          target = input$dens_target,
          binwidth_A = as.numeric(input$dens_binwidth),
          zlim_A = input$dens_zlim,
          center_membrane = isTRUE(input$mem_center_profiles),
          dt_ps = input$dt_ps, first = first, last = last, stride = stride,
          combine = isTRUE(input$combine_trajs)
        )
        rv$data$mem_density <- df
        if (!is.null(rv$out_dir) && dir.exists(rv$out_dir)) {
          write.csv(df, file.path(rv$out_dir, "membrane_density_profile.csv"), row.names = FALSE)
        }
        incProgress(1)
      })
      show_themed_completion(
        title_text = "Membrane analysis finished",
        hero_title = "Membrane density profiles computed",
        subtitle = "Component density profiles along the bilayer normal are ready. Results were also exported as a CSV file.",
        steps = c(
          "Go to the Membrane systems tab and inspect the Density profiles plot.",
          "Headgroups should peak near the two leaflet interfaces; tails in the hydrophobic core; water/ions outside the bilayer.",
          "Download the CSV (membrane_density_profile.csv) for publication-ready profile plots in external tools."
        )
      )
    }, error = function(e) {
      showNotification(paste0("Density profile failed: ", e$message), type = "error", duration = 10)
      rv$log <- paste0(rv$log, "ERROR (density profiles): ", e$message, "\n")
    })
  })

  observeEvent(input$compute_order, {
    req(rv$proj_dir, rv$prmtop, length(rv$trajs) >= 1)
    trajs <- NULL
    tryCatch({ trajs <- ordered_trajs() }, error = function(e) {
      showNotification(paste0("Trajectory order error: ", e$message), type = "error")
    })
    if (is.null(trajs)) return()
    show_traj_subset_modal("order", "Compute lipid tail order — select trajectory files", trajs)
  })

  observeEvent(input$confirm_order, {
    req(rv$proj_dir, rv$prmtop, length(rv$trajs) >= 1)
    mem_resid <- parse_csv_tokens(input$mem_resid)
    first <- if (is.na(input$first)) NULL else as.integer(input$first)
    last  <- if (is.na(input$last)) NULL else as.integer(input$last)
    stride <- as.integer(input$stride %||% 1)

    trajs_all <- NULL
    tryCatch({ trajs_all <- ordered_trajs() }, error = function(e) {
      showNotification(paste0("Trajectory order error: ", e$message), type = "error")
    })
    if (is.null(trajs_all)) return()
    trajs <- get_modal_trajs("order", trajs_all)
    if (is.null(trajs)) return()
    removeModal()

    rv$log <- paste0(rv$log, "Running lipid tail order calculation...\n")
    tryCatch({
      progress_with_waiter("Computing lipid tail order",
                           "Reading the selected trajectories and estimating per-segment order parameters.",
                           message = "Computing lipid tail order...", value = 0, {
        prm_use <- (rv$prm_obj %||% read.prmtop(rv$prmtop))
        df <- compute_lipid_tail_order(
          trajs_ordered = trajs, prm = prm_use, membrane_resid = mem_resid,
          chain1_regex = input$order_chain1_regex,
          chain2_regex = input$order_chain2_regex,
          by_type = isTRUE(input$order_by_type),
          topology_mode = input$order_topology_mode %||% "standard",
          include_resid = input$order_resid_filter,
          chain1_resid = input$order_chain1_resid,
          chain2_resid = input$order_chain2_resid,
          dt_ps = input$dt_ps, first = first, last = last, stride = stride,
          combine = isTRUE(input$combine_trajs)
        )
        rv$data$mem_order <- df
        if (!is.null(rv$out_dir) && dir.exists(rv$out_dir)) {
          write.csv(df, file.path(rv$out_dir, "lipid_tail_order_profile.csv"), row.names = FALSE)
        }
        incProgress(1)
      })
      show_themed_completion(
        title_text = "Membrane analysis finished",
        hero_title = "Lipid tail order (|S|) profile computed",
        subtitle = "Segmental order parameters along the acyl chains are ready. Results were also exported as a CSV file.",
        steps = c(
          "Go to the Membrane systems tab and inspect the Lipid tail order profile plot.",
          "Higher |S| values near the glycerol backbone indicate ordered chains; values decrease toward the terminal methyl.",
          "If only one point appears for a lipid type, check that the tail regex correctly matches the expected atom names.",
          "Download the CSV (lipid_tail_order_profile.csv) to compare order parameters across conditions or systems."
        )
      )
    }, error = function(e) {
      showNotification(paste0("Lipid order failed: ", e$message), type = "error", duration = 10)
      rv$log <- paste0(rv$log, "ERROR (lipid order): ", e$message, "\n")
    })
  })

  observeEvent(input$compute_interactions, {
    req(rv$proj_dir, rv$prmtop, length(rv$trajs) >= 1)
    trajs <- NULL
    tryCatch({ trajs <- ordered_trajs() }, error = function(e) {
      showNotification(paste0("Trajectory order error: ", e$message), type = "error")
    })
    if (is.null(trajs)) return()
    show_traj_subset_modal("interactions", "Compute Selection A / Selection B interactions — select trajectory files", trajs)
  })

  observeEvent(input$confirm_interactions, {
    req(rv$proj_dir, rv$prmtop, length(rv$trajs) >= 1)
    exclude_resid <- parse_csv_tokens(input$exclude_resid)
    pep_resno <- parse_resno(input$pep_resno)
    lipo_resno <- parse_resno(input$lipo_resno)
    mem_resid <- parse_csv_tokens(input$mem_resid)
    align_elety <- if (identical(input$align_mode, "ca")) "CA" else c("N","CA","C","O")
    first <- if (is.na(input$first)) NULL else as.integer(input$first)
    last  <- if (is.na(input$last)) NULL else as.integer(input$last)
    stride <- as.integer(input$stride %||% 1)

    trajs_all <- NULL
    tryCatch({ trajs_all <- ordered_trajs() }, error = function(e) {
      showNotification(paste0("Trajectory order error: ", e$message), type = "error")
    })
    if (is.null(trajs_all)) return()
    trajs <- get_modal_trajs("interactions", trajs_all)
    if (is.null(trajs)) return()
    removeModal()

    rv$log <- paste0(rv$log, "Running A/B interaction analysis...\n")
    tryCatch({
      progress_with_waiter("Computing A/B interactions",
                           "Reading the selected trajectories and evaluating distance, contact, and occupancy metrics.",
                           message = "Computing A/B interactions...", value = 0, {
        prm_use <- (rv$prm_obj %||% read.prmtop(rv$prmtop))
        sels <- make_selections(
          prm = prm_use,
          pep_resno = pep_resno,
          lipo_resno = lipo_resno,
          membrane_resid = mem_resid,
          exclude_resid = exclude_resid,
          align_elety = align_elety
        )
        res <- compute_interaction_metrics(
          trajs_ordered = trajs, prm = prm_use, sels = sels,
          selA_key = input$interact_selA,
          selB_key = input$interact_selB,
          cutoff_A = as.numeric(input$interact_cutoff),
          dt_ps = input$dt_ps, first = first, last = last, stride = stride,
          combine = isTRUE(input$combine_trajs)
        )
        rv$data$interact_ts <- res$ts
        rv$data$interact_occ <- res$occupancy
        rv$data$interact_pairs <- res$pairs
        if (!is.null(rv$out_dir) && dir.exists(rv$out_dir)) {
          write.csv(res$ts, file.path(rv$out_dir, "ab_interaction_timeseries.csv"), row.names = FALSE)
          write.csv(res$occupancy, file.path(rv$out_dir, "ab_interaction_occupancy.csv"), row.names = FALSE)
          if (!is.null(res$pairs) && nrow(res$pairs) > 0) write.csv(res$pairs, file.path(rv$out_dir, "ab_interaction_pair_occupancy.csv"), row.names = FALSE)
        }
        incProgress(1)
      })
      show_themed_completion(
        title_text = "Interaction analysis finished",
        hero_title = "Selection A / Selection B interaction metrics computed",
        subtitle = "Distance, contact, and residue occupancy data are ready. Results were also exported as CSV files.",
        steps = c(
          "Go to the Interactions tab and inspect the time series plots (minimum distance, COM distance, atom contacts).",
          "Check the Residue contact occupancy panel to identify the most persistent interface residues.",
          "Use the A/B residue contact map to pinpoint which pairs form the dominant interaction patch.",
          "Download the CSV tables (ab_interaction_timeseries.csv, ab_interaction_occupancy.csv) for further analysis."
        )
      )
    }, error = function(e) {
      showNotification(paste0("Interaction analysis failed: ", e$message), type = "error", duration = 10)
      rv$log <- paste0(rv$log, "ERROR (interactions): ", e$message, "\n")
    })
  })


  observeEvent(input$compute_hbonds, {
    req(rv$proj_dir, rv$prmtop, length(rv$trajs) >= 1)
    trajs <- NULL
    tryCatch({ trajs <- ordered_trajs() }, error = function(e) {
      showNotification(paste0("Trajectory order error: ", e$message), type = "error")
    })
    if (is.null(trajs)) return()
    show_traj_subset_modal("hbonds", "Compute Selection A / Selection B H-bond proxy — select trajectory files", trajs)
  })

  observeEvent(input$confirm_hbonds, {
    req(rv$proj_dir, rv$prmtop, length(rv$trajs) >= 1)
    exclude_resid <- parse_csv_tokens(input$exclude_resid)
    pep_resno <- parse_resno(input$pep_resno)
    lipo_resno <- parse_resno(input$lipo_resno)
    mem_resid <- parse_csv_tokens(input$mem_resid)
    align_elety <- if (identical(input$align_mode, "ca")) "CA" else c("N","CA","C","O")
    first <- if (is.na(input$first)) NULL else as.integer(input$first)
    last  <- if (is.na(input$last)) NULL else as.integer(input$last)
    stride <- as.integer(input$stride %||% 1)

    trajs_all <- NULL
    tryCatch({ trajs_all <- ordered_trajs() }, error = function(e) {
      showNotification(paste0("Trajectory order error: ", e$message), type = "error")
    })
    if (is.null(trajs_all)) return()
    trajs <- get_modal_trajs("hbonds", trajs_all)
    if (is.null(trajs)) return()
    removeModal()

    rv$log <- paste0(rv$log, "Running A/B hydrogen-bond proxy analysis...\n")
    tryCatch({
      progress_with_waiter("Computing A/B H-bond proxy",
                           "Reading the selected trajectories and evaluating hydrogen-bond proxy contacts.",
                           message = "Computing A/B H-bond proxy...", value = 0, {
        prm_use <- (rv$prm_obj %||% read.prmtop(rv$prmtop))
        sels <- make_selections(
          prm = prm_use,
          pep_resno = pep_resno,
          lipo_resno = lipo_resno,
          membrane_resid = mem_resid,
          exclude_resid = exclude_resid,
          align_elety = align_elety
        )
        res <- compute_hbond_proxy_metrics(
          trajs_ordered = trajs, prm = prm_use, sels = sels,
          cutoff_A = as.numeric(input$hbond_cutoff %||% 3.5),
          include_n_acceptors = isTRUE(input$hbond_include_nacc),
          dt_ps = input$dt_ps, first = first, last = last, stride = stride,
          combine = isTRUE(input$combine_trajs)
        )
        rv$data$hbond_ts <- res$ts
        rv$data$hbond_pairs <- res$pairs
        if (!is.null(rv$out_dir) && dir.exists(rv$out_dir)) {
          write.csv(res$ts, file.path(rv$out_dir, "ab_hbond_proxy_timeseries.csv"), row.names = FALSE)
          if (!is.null(res$pairs) && nrow(res$pairs) > 0) write.csv(res$pairs, file.path(rv$out_dir, "ab_hbond_proxy_pairs.csv"), row.names = FALSE)
        }
        incProgress(1)
      })
      show_themed_completion(
        title_text = "Hydrogen-bond proxy analysis finished",
        hero_title = "A / B hydrogen-bond proxy metrics computed",
        subtitle = "H-bond proxy counts and donor/acceptor pair occupancies are ready. Results were also exported as CSV files.",
        steps = c(
          "Go to the Interactions tab and inspect the H-bond proxy time series plot.",
          "Check the Persistent donor/acceptor pairs panel to identify the most stable polar contacts across the trajectory.",
          "Remember: this is a distance-based proxy without angle filtering — confirm key pairs in a 3D viewer if needed.",
          "Download the CSV tables (ab_hbond_proxy_timeseries.csv, ab_hbond_proxy_pairs.csv) for further analysis."
        )
      )
    }, error = function(e) {
      showNotification(paste0("H-bond proxy analysis failed: ", e$message), type = "error", duration = 10)
      rv$log <- paste0(rv$log, "ERROR (H-bond proxy): ", e$message, "\n")
    })
  })

  observeEvent(input$run_clust, {
    req(rv$proj_dir, rv$prmtop, length(rv$trajs) >= 1)
    trajs <- NULL
    tryCatch({ trajs <- ordered_trajs() }, error=function(e) {
      rv$log <- paste0(rv$log, "ERROR (clustering): ", e$message, "\n")
    })
    if (is.null(trajs)) return()
    show_traj_subset_modal("clust", "Run RMSD clustering — select trajectory files", trajs, action_label = "Run clustering")
  })

  observeEvent(input$confirm_clust, {
    req(rv$proj_dir, rv$prmtop, length(rv$trajs) >= 1)
    rv$log <- paste0(rv$log, "Running RMSD clustering...\n")

    trajs_all <- NULL
    tryCatch({
      trajs_all <- ordered_trajs()
    }, error=function(e) {
      rv$log <- paste0(rv$log, "ERROR (clustering): ", e$message, "\n")
    })
    if (is.null(trajs_all)) return()
    trajs <- get_modal_trajs("clust", trajs_all)
    if (is.null(trajs)) return()
    removeModal()

    pep_resno <- parse_resno(input$pep_resno)
    lipo_resno <- parse_resno(input$lipo_resno)
    mem_resid <- parse_csv_tokens(input$mem_resid)
    exclude_resid <- parse_csv_tokens(input$exclude_resid)

    align_elety <- if (identical(input$align_mode, "ca")) "CA" else c("N","CA","C","O")
    mem_target <- input$mem_target %||% "B"

    if (isTRUE(input$clust_use_project_window)) {
      c_first <- if (is.na(input$first)) NULL else as.integer(input$first)
      c_last  <- if (is.na(input$last))  NULL else as.integer(input$last)
      c_stride <- as.integer(input$stride)
    } else {
      c_first <- if (is.na(input$clust_first)) NULL else as.integer(input$clust_first)
      c_last  <- if (is.na(input$clust_last))  NULL else as.integer(input$clust_last)
      c_stride <- as.integer(input$clust_stride)
    }

    plan_df <- build_frame_index(trajs, dt_ps = input$dt_ps, first = c_first, last = c_last,
                               stride = as.integer(c_stride), combine = isTRUE(input$combine_trajs))
    rv$log <- paste0(rv$log, sprintf("Clustering workload: %d segments, %d frames after stride.\n",
                                    length(unique(plan_df$segment)), nrow(plan_df)))

    tryCatch({
      progress_with_waiter("Running structural clustering",
                           "Reading trajectories, building the RMSD distance matrix, and assigning clusters.",
                           message = "Clustering frames (RMSD)...", value = 0, {
        t0 <- Sys.time()
        logf <- function(msg) { rv$log <- paste0(rv$log, msg, "\n") }
        progress_cb <- function(value, detail = NULL) { setProgress(value = value, detail = detail) }
        incProgress(0.05, detail = "Initializing")
        rv$clust <- run_rmsd_structural_clustering(
          project_dir = rv$proj_dir,
          trajs_ordered = trajs,
          combine = isTRUE(input$combine_trajs),
          dt_ps = input$dt_ps,
          first = c_first,
          last = c_last,
          stride = c_stride,
          pep_resno = pep_resno,
          lipo_resno = lipo_resno,
          membrane_resid = mem_resid,
          exclude_resid = exclude_resid,
          align_elety = align_elety,
          clust_atoms_key = input$clust_atoms,
          alg = input$clust_alg,
          ref_mode = input$clust_ref_mode,
          ref_time_ns = as.numeric(input$clust_ref_time %||% 0),
          ref_coords_path = if (isTRUE(input$use_rmsd_ref_file) && !is.null(input$rmsd_ref_file$datapath) && nzchar(input$rmsd_ref_file$datapath)) input$rmsd_ref_file$datapath else NULL,
          ref_coords_name = if (isTRUE(input$use_rmsd_ref_file) && !is.null(input$rmsd_ref_file$name) && nzchar(input$rmsd_ref_file$name)) input$rmsd_ref_file$name else NULL,
          k = as.integer(input$clust_k),
          progress = progress_cb,
          logf = logf,
          verbose = TRUE
        )
        setProgress(1, detail = "Finalizing")
        dt_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        rv$log <- paste0(rv$log, sprintf("Clustering runtime: %.1f s\n", dt_sec))
      })

      rv$log <- paste0(
        rv$log,
        "Clustering done. ",
        "K=", max(rv$clust$frames$cluster), ", frames=", nrow(rv$clust$frames), ".\n"
      )

      # Save the clustering log next to the exported cluster tables and PDBs.
      writeLines(rv$log, con = file.path(cluster_results_dir(rv$proj_dir), "clustering_log.txt"))

      show_themed_completion(
        title_text = "Clustering finished",
        hero_title = "RMSD clustering completed successfully",
        subtitle = sprintf("Generated %d clusters across %d analyzed frames.", max(rv$clust$frames$cluster), nrow(rv$clust$frames)),
        path_value = cluster_results_dir(rv$proj_dir),
        path_label = "Export folder",
        steps = c(
          "Inspect the Cluster plots panel to review the distribution, cluster-vs-time, population, and dendrogram views.",
          "Use the Toggle 3D distribution plot button inside Cluster plots to switch between 2D and 3D views.",
          "Use the export controls on the left to write medoids, centroids, or a specific frame as PDB."
        )
      )
    }, error=function(e) {
      rv$log <- paste0(rv$log, "ERROR (clustering): ", e$message, "\n")
      rv$clust <- NULL
    })
  })

  
  observeEvent(input$export_clust, {
    req(rv$clust, rv$proj_dir)

    out_dir <- cluster_results_dir(rv$proj_dir)
    dir_create(out_dir, recurse = TRUE)

    tryCatch({
      progress_with_waiter("Exporting cluster structures",
                           "Preparing representative structures and writing clustering export files.",
                           message = "Exporting cluster structures...", value = 0, {
        incProgress(0.05, detail = "Writing tables")

      # Always export tables
      write.csv(rv$clust$frames,  file.path(out_dir, "cluster_assignments.csv"), row.names = FALSE)
      write.csv(rv$clust$summary, file.path(out_dir, "clusters_summary.csv"), row.names = FALSE)
      write.csv(rv$clust$plot,    file.path(out_dir, "cluster_plot_data.csv"), row.names = FALSE)
      if (!is.null(rv$clust$silhouette_summary)) write.csv(rv$clust$silhouette_summary, file.path(out_dir, "cluster_quality_silhouette.csv"), row.names = FALSE)
      if (!is.null(rv$clust$silhouette)) write.csv(rv$clust$silhouette, file.path(out_dir, "cluster_frame_silhouette.csv"), row.names = FALSE)

      # How many clusters?
      k_eff <- nrow(rv$clust$summary)
      n_keep <- as.integer(input$clust_export_n %||% 0)
      if (!is.na(n_keep) && n_keep > 0) k_eff <- min(k_eff, n_keep)

      summ <- rv$clust$summary[seq_len(k_eff), , drop = FALSE]

      # Build export atom selection
      prm <- read.prmtop(rv$clust$prmtop_path)
      incProgress(0.10, detail = "Preparing PDB template")

      pep_resno <- parse_resno(input$pep_resno)
      lipo_resno <- parse_resno(input$lipo_resno)
      exclude_resid <- parse_csv_tokens(input$exclude_resid)
      mem_resid <- parse_csv_tokens(input$mem_resid)

      export_key <- input$clust_export_atoms %||% "sys_heavy"
      natom <- get_prm_natom(prm)
      atom_inds_export <- get_export_atom_inds(prm, mode = export_key, exclude_resid = exclude_resid)
      if (length(atom_inds_export) < 1) {
        stop("Export selection empty. Check exclude list or choose a full-system export option.")
      }
      atom_inds_export <- clean_atom_inds(atom_inds_export, natom)
      if (length(atom_inds_export) < 1) stop("Export atom selection empty after cleaning.")
      pdb_tmpl_export <- make_pdb_template_from_prmtop(prm, atom_inds_export)
      at_sel_export <- as.select(atom_inds_export)

      # Convenience: map segment basename -> full path
      trajs_map <- rv$clust$trajs_map

      # Export medoid + centroid-nearest-frame (REAL frames) for each cluster
      step <- 0.80 / max(1, nrow(summ))
      for (i in seq_len(nrow(summ))) {
        incProgress(step, detail = sprintf("Exporting cluster %d/%d", i, nrow(summ)))
        cc <- summ$cluster[i]

        # Medoid frame
        seg_m <- summ$medoid_segment[i]
        fr_m  <- summ$medoid_frame[i]
        p_m <- trajs_map[[seg_m]]
        if (is.null(p_m)) stop("Could not locate trajectory file for segment: ", seg_m)

        xyz_m <- read.ncdf(p_m, first = fr_m, last = fr_m, stride = 1, at.sel = at_sel_export, verbose = FALSE)
        write_pdb_with_xyz(pdb_tmpl_export, xyz_m[1, ], file.path(out_dir, sprintf("cluster_%02d_medoid.pdb", cc)))

        # Centroid-nearest real frame
        seg_c <- summ$centroid_frame_segment[i]
        fr_c  <- summ$centroid_frame[i]
        if (!is.null(seg_c) && nzchar(seg_c) && !is.na(fr_c)) {
          p_c <- trajs_map[[seg_c]]
          if (is.null(p_c)) stop("Could not locate trajectory file for segment: ", seg_c)

          xyz_c <- read.ncdf(p_c, first = fr_c, last = fr_c, stride = 1, at.sel = at_sel_export, verbose = FALSE)
          write_pdb_with_xyz(pdb_tmpl_export, xyz_c[1, ], file.path(out_dir, sprintf("cluster_%02d_centroid_frame.pdb", cc)))
          # Backward-compatible centroid name
          write_pdb_with_xyz(pdb_tmpl_export, xyz_c[1, ], file.path(out_dir, sprintf("cluster_%02d_centroid.pdb", cc)))
        }

        # Optional: centroid mean coordinates (NON-physical; only valid for the same atoms used for clustering)
        if (isTRUE(input$clust_export_centroid_mean)) {
          if (identical(export_key, rv$clust$clust_atoms_key) && !identical(export_key, "full")) {
            # centroid_xyz is in clustering atom order, which matches sels[[export_key]]$atom
            pdb_tmpl_cent <- pdb_tmpl_export
            write_pdb_with_xyz(
              pdb_tmpl_cent,
              rv$clust$centroid_xyz[[cc]],
              file.path(out_dir, sprintf("cluster_%02d_centroid_mean.pdb", cc))
            )
          } else {
            # Still provide centroid mean on clustering atoms so it is reproducible
            pdb_tmpl_clust <- make_pdb_template_from_prmtop(prm, rv$clust$clust_atom_inds)
            write_pdb_with_xyz(
              pdb_tmpl_clust,
              rv$clust$centroid_xyz[[cc]],
              file.path(out_dir, sprintf("cluster_%02d_centroid_mean_CLUSTATOMS.pdb", cc))
            )
          }
        }
      }

      })

      pdb_written <- list.files(out_dir, pattern = "\\.pdb$", full.names = FALSE)
      rv$log <- paste0(rv$log, "PDB files written: ", length(pdb_written), "\n")
      if (length(pdb_written) > 0) {
        show_list <- head(pdb_written, 12)
        rv$log <- paste0(rv$log, "Examples: ", paste(show_list, collapse = ", "), if (length(pdb_written) > 12) ", ..." else "", "\n")
      }

      rv$log <- paste0(rv$log, "Cluster export done: ", out_dir, "\n")
      writeLines(rv$log, con = file.path(out_dir, "clustering_log.txt"))

      show_themed_completion(
        title_text = "Cluster export finished",
        hero_title = "Representative structures exported",
        subtitle = sprintf("Wrote %d PDB file(s) plus clustering tables.", length(pdb_written)),
        path_value = out_dir,
        path_label = "Export folder",
        steps = c(
          "Open the exported medoid and centroid-frame PDB files in PyMOL, VMD, or your preferred viewer.",
          "The same folder also contains cluster_assignments.csv, clusters_summary.csv, and cluster_plot_data.csv.",
          "If needed, you can now export an additional single frame using the controls just below."
        )
      )
    }, error=function(e) {
      rv$log <- paste0(rv$log, "ERROR (cluster export): ", e$message, "\n")
    })
  })

  observeEvent(input$export_single_frame, {
    req(rv$proj_dir)
    if (is.null(rv$clust)) {
      showNotification(
        "Run clustering first. 'Export selected frame' uses the current clustering results, frame map, and clustering window.",
        type = "warning",
        duration = 8
      )
      return()
    }

    out_dir <- cluster_results_dir(rv$proj_dir)
    dir_create(out_dir, recurse = TRUE)

    tryCatch({
      progress_with_waiter("Exporting selected frame",
                           "Preparing the requested frame and writing the exported PDB file.",
                           message = "Exporting selected frame...", value = 0, {
        incProgress(0.15, detail = "Preparing selection")
      prm <- read.prmtop(rv$clust$prmtop_path)
      natom <- get_prm_natom(prm)

      pep_resno <- parse_resno(input$pep_resno)
      lipo_resno <- parse_resno(input$lipo_resno)
      exclude_resid <- parse_csv_tokens(input$exclude_resid)
      mem_resid <- parse_csv_tokens(input$mem_resid)

      export_key <- input$clust_export_atoms %||% "sys_heavy"
      atom_inds_export <- get_export_atom_inds(prm, mode = export_key, exclude_resid = exclude_resid)
      if (length(atom_inds_export) < 1) stop("Export selection empty. Check exclude list or choose a full-system export option.")
      atom_inds_export <- clean_atom_inds(atom_inds_export, natom)
      if (length(atom_inds_export) < 1) stop("Export atom selection empty.")
      pdb_tmpl_export <- make_pdb_template_from_prmtop(prm, atom_inds_export)
      at_sel_export <- as.select(atom_inds_export)

      trajs_map <- rv$clust$trajs_map

      mode <- input$single_frame_mode %||% "time"
      seg <- NULL; fr <- NULL; t_used <- NA_real_

      if (identical(mode, "time")) {
        df <- rv$clust$frames
        if (isTRUE(input$combine_trajs) && ("time_ns_global" %in% names(df)) && any(!is.na(df$time_ns_global))) {
          df$time_use <- df$time_ns_global
        } else {
          df$time_use <- df$time_ns
        }
        target_t <- as.numeric(input$single_frame_time %||% 0)
        i <- which.min(abs(df$time_use - target_t))
        seg <- df$segment[i]
        fr  <- df$frame[i]
        t_used <- df$time_use[i]
      } else if (identical(mode, "segframe")) {
        df <- rv$clust$frames
        frame_req <- as.integer(input$single_frame_frame %||% 1)
        if (!is.finite(frame_req) || is.na(frame_req) || frame_req < 1) stop("Frame index must be a positive integer.")

        use_global <- isTRUE(input$combine_trajs) && ("frame_global" %in% names(df)) && any(!is.na(df$frame_global))
        cand <- if (use_global) {
          df[df$frame_global == frame_req, , drop = FALSE]
        } else {
          df[df$frame == frame_req, , drop = FALSE]
        }

        if (nrow(cand) < 1) stop("Frame index not found in the current clustering window.")
        if (nrow(cand) > 1) stop("That frame index matches multiple segments. Enable 'Combine segments' or use 'Global time (ns)'.")

        seg <- cand$segment[1]
        fr  <- cand$frame[1]
        t_used <- if (use_global && "time_ns_global" %in% names(cand) && !is.na(cand$time_ns_global[1])) cand$time_ns_global[1] else cand$time_ns[1]
      } else if (identical(mode, "clusterrep")) {
        cc <- as.integer(input$single_frame_cluster %||% 1)
        rep_kind <- input$single_frame_rep %||% "medoid"
        summ <- rv$clust$summary
        if (cc < 1 || cc > nrow(summ)) stop("Cluster # out of range.")
        if (identical(rep_kind, "centroid_frame")) {
          seg <- summ$centroid_frame_segment[cc]
          fr  <- summ$centroid_frame[cc]
          t_used <- summ$centroid_frame_time_ns_global[cc] %||% summ$centroid_frame_time_ns[cc]
        } else {
          seg <- summ$medoid_segment[cc]
          fr  <- summ$medoid_frame[cc]
          t_used <- summ$medoid_time_ns_global[cc] %||% summ$medoid_time_ns[cc]
        }
      }

      if (is.null(seg) || !nzchar(seg) || is.na(fr)) stop("Could not determine the frame to export.")
      seg_path <- trajs_map[[seg]]
      if (is.null(seg_path)) stop("Could not locate trajectory file for segment: ", seg)

      xyz <- read.ncdf(seg_path, first = fr, last = fr, stride = 1, at.sel = at_sel_export, verbose = FALSE)
      incProgress(0.75, detail = "Writing PDB")

      fn <- sprintf("export_frame_%s_frame_%d.pdb", tools::file_path_sans_ext(seg), fr)
      if (!is.na(t_used)) {
        fn <- sprintf("export_frame_%s_frame_%d_t%.3fns.pdb", tools::file_path_sans_ext(seg), fr, t_used)
      }
      write_pdb_with_xyz(pdb_tmpl_export, xyz[1, ], file.path(out_dir, fn))

      })
      rv$log <- paste0(rv$log, "Exported frame PDB: ", file.path(out_dir, fn), "\n")
      writeLines(rv$log, con = file.path(out_dir, "clustering_log.txt"))

      show_themed_completion(
        title_text = "Frame export finished",
        hero_title = "Specific frame exported successfully",
        subtitle = "Your selected frame was written as a PDB file and is ready for visualization.",
        path_value = file.path(out_dir, fn),
        path_label = "Exported file",
        steps = c(
          "Open the PDB in PyMOL, VMD, Chimera, or another structure viewer.",
          "If you need a different frame, just change the frame selector and export again.",
          "You can keep using the same export atom selection menu to switch between heavy-only, backbone-only, or full-system exports."
        )
      )
    }, error=function(e) {
      rv$log <- paste0(rv$log, "ERROR (export frame): ", e$message, "\n")
    })
  })


  output$plot_dimred <- renderPlotly({
    req(rv$data$pca_scores)
    st <- style_for("dimred", defaults = list(show_points = TRUE))
    df <- rv$data$pca_scores
    validate(need(all(c("PC1", "PC2") %in% names(df)), "PCA scores were not generated."))

    if (isTRUE(input$combine_trajs) && ("time_ns_global" %in% names(df)) && any(!is.na(df$time_ns_global))) {
      df$time_use <- df$time_ns_global
    } else {
      df$time_use <- df$time_ns
    }

    hover <- paste0(
      "Segment: ", df$segment,
      "<br>Frame: ", df$frame,
      "<br>Time (ns): ", sprintf("%.3f", df$time_use)
    )

    if (isTRUE(rv$dimred_plot_3d) && ("PC3" %in% names(df))) {
      plot_ly(
        df,
        x = ~PC1, y = ~PC2, z = ~PC3,
        type = "scatter3d", mode = "markers",
        text = hover, hoverinfo = "text",
        marker = list(
          size = st$point_size,
          opacity = st$point_alpha,
          color = df$time_use,
          colorscale = "Viridis",
          showscale = TRUE,
          colorbar = list(title = list(text = "Time (ns)", font = list(size = 12, color = "#e8ecf1")), tickfont = list(color = "#c0ccda", size = 11))
        )
      ) |>
        layout(scene = list(xaxis = list(title = "PC1"),
                            yaxis = list(title = "PC2"),
                            zaxis = list(title = "PC3"))) |> apply_plotly_theme()
    } else {
      plot_ly(
        df,
        x = ~PC1, y = ~PC2,
        type = "scatter", mode = "markers",
        text = hover, hoverinfo = "text",
        marker = list(
          size = st$point_size,
          opacity = st$point_alpha,
          color = df$time_use,
          colorscale = "Viridis",
          showscale = TRUE,
          colorbar = list(title = list(text = "Time (ns)", font = list(size = 12, color = "#e8ecf1")), tickfont = list(color = "#c0ccda", size = 11))
        )
      ) |>
        layout(xaxis = plotly_axis_from_style(st, "x", "PC1"),
               yaxis = plotly_axis_from_style(st, "y", "PC2")) |> apply_plotly_theme()
    }
  })

  output$plot_pca_var <- renderPlotly({
    req(rv$data$pca_var)
    st <- style_for("pcaVar", defaults = list(show_points = TRUE, line_width = 0.8))
    df <- normalize_pca_var_df(head(rv$data$pca_var, 10))
    validate(need(nrow(df) > 0, "PCA variance table is empty."))
    validate(need(all(c("prop_var","cum_var","PC_index","PC") %in% names(df)), "PCA variance table is missing required columns."))

    xax <- plotly_axis_from_style(st, "x", "Principal component")
    xax$tickmode <- "array"
    xax$tickvals <- df$PC_index
    xax$ticktext <- df$PC

    yax <- plotly_axis_from_style(st, "y", "Explained variance")
    yax$tickformat <- ".0%"

    plot_ly(df, x = ~PC_index, y = ~prop_var, type = "bar",
            text = ~paste0(round(prop_var * 100, 1), "%"),
            textposition = "outside", name = "Explained variance") |>
      add_lines(data = df, x = ~PC_index, y = ~cum_var, yaxis = "y2",
                name = "Cumulative variance",
                line = list(width = st$line_width)) |>
      layout(
        xaxis = xax,
        yaxis = yax,
        yaxis2 = list(title = "Cumulative variance", overlaying = "y", side = "right", tickformat = ".0%")
      ) |> apply_plotly_theme()
  })

output$cluster_ref_info <- renderUI({
    req(rv$clust)
    ri <- rv$clust$ref_info
    # Prefer global time if available
    tshow <- if (!is.na(ri$time_ns_global)) ri$time_ns_global else ri$time_ns
    tagList(
      tags$div(class="info-panel",
               tags$b("Reference used for alignment: "),
               tags$span(ri$segment),
               tags$span(" — frame "),
               tags$span(ri$frame),
               tags$span(" — time "),
               tags$span(sprintf("%.3f ns", tshow))
      )
    )
  })

  output$cluster_quality_info <- renderUI({
    req(rv$clust)
    sil_overall <- rv$clust$silhouette_overall %||% NA_real_
    if (!is.finite(sil_overall)) {
      return(tags$div(class = "info-panel",
                      tags$b("Silhouette summary: "),
                      tags$span("not available"),
                      tags$br(),
                      tags$small("Install/use package 'cluster' and at least 2 non-trivial clusters to enable this metric.")))
    }
    qual_lab <- if (sil_overall >= 0.5) "good separation"
    else if (sil_overall >= 0.25) "moderate separation"
    else "weak separation"
    tags$div(class = "info-panel",
             tags$b("Overall mean silhouette: "),
             tags$span(sprintf("%.3f", sil_overall)),
             tags$br(),
             tags$small(sprintf("Interpretation: %s.", qual_lab)))
  })

  output$plot_cluster_dist <- renderPlotly({
    req(rv$clust)
    st <- style_for("clustDist")
    df <- rv$clust$plot

    # choose time label for hover
    t_hover <- if (isTRUE(input$combine_trajs) && any(!is.na(df$time_ns_global))) df$time_ns_global else df$time_ns
    hover <- paste0(
      "Cluster: ", df$cluster,
      "<br>Segment: ", df$segment,
      "<br>Frame: ", df$frame,
      "<br>Time (ns): ", sprintf("%.3f", t_hover)
    )

    if (isTRUE(rv$clust_plot_3d)) {
      plot_ly(
        df,
        x = ~Dim1, y = ~Dim2, z = ~Dim3,
        type = "scatter3d", mode = "markers",
        color = ~cluster, text = hover, hoverinfo = "text",
        marker = list(size = st$point_size, opacity = st$point_alpha)
      ) |>
        layout(scene = list(xaxis = list(title = "Dim1"),
                            yaxis = list(title = "Dim2"),
                            zaxis = list(title = "Dim3"))) |> apply_plotly_theme()
    } else {
      plot_ly(
        df,
        x = ~Dim1, y = ~Dim2,
        type = "scatter", mode = "markers",
        color = ~cluster, text = hover, hoverinfo = "text",
        marker = list(size = st$point_size, opacity = st$point_alpha)
      ) |>
        layout(xaxis = plotly_axis_from_style(st, "x", "Dim1"),
               yaxis = plotly_axis_from_style(st, "y", "Dim2")) |> apply_plotly_theme()
    }
  })

  
  output$plot_cluster_time <- renderPlotly({
    req(rv$clust)
    st <- style_for("clustTime", defaults = list(point_size = 6))
    df <- rv$clust$frames
    validate(need("cluster" %in% names(df), "Run clustering first."))

    # choose time axis: global if segments combined
    if (isTRUE(input$combine_trajs) && ("time_ns_global" %in% names(df)) && any(!is.na(df$time_ns_global))) {
      df$time_use <- df$time_ns_global
    } else {
      df$time_use <- df$time_ns
    }

    df$cluster_f <- factor(df$cluster, levels = sort(unique(df$cluster)))

    plotly::plot_ly(
      df,
      x = ~time_use,
      y = ~as.numeric(cluster),
      type = "scatter",
      mode = "markers",
      color = ~cluster_f,
      marker = list(size = st$point_size, opacity = st$point_alpha),
      text = ~paste0("Cluster: ", cluster,
                     "<br>Segment: ", segment,
                     "<br>Frame: ", frame,
                     "<br>Time (ns): ", sprintf("%.3f", time_use)),
      hoverinfo = "text"
    ) |>
      plotly::layout(
        xaxis = plotly_axis_from_style(st, "x", "Time (ns)"),
        yaxis = modifyList(
          plotly_axis_from_style(st, "y", "Cluster"),
          list(
            tickmode = "array",
            tickvals = sort(unique(as.numeric(df$cluster))),
            ticktext = sort(unique(as.numeric(df$cluster)))
          )
        )
      ) |> apply_plotly_theme()
  })

  output$plot_cluster_pop <- renderPlotly({
    req(rv$clust)
    st <- style_for("clustPop", defaults = list(show_points = TRUE, line_width = 0.8))
    df <- rv$clust$summary
    validate(need(nrow(df) > 0, "Run clustering first."))

    plot_ly(df,
            x = ~factor(cluster),
            y = ~n_frames,
            type = "bar",
            text = ~paste0(round(pct, 1), "%"),
            textposition = "outside") |>
      layout(xaxis = plotly_axis_from_style(st, "x", "Cluster"),
             yaxis = plotly_axis_from_style(st, "y", "Number of frames")) |> apply_plotly_theme()
  })

  output$plot_cluster_quality <- renderPlotly({
    req(rv$clust)
    st <- style_for("clustQual", defaults = list(show_points = TRUE, line_width = 0.8))
    df <- rv$clust$silhouette_summary
    sil_overall <- rv$clust$silhouette_overall %||% NA_real_
    validate(need(!is.null(df) && nrow(df) > 0 && any(is.finite(df$mean_silhouette)),
                  "Silhouette quality is not available for this clustering result."))

    p <- make_cluster_quality_gg(df, sil_overall = sil_overall, title = "Cluster quality (mean silhouette)", style = st)
    ggplotly(p, tooltip = c("x", "y")) |> apply_plotly_theme()
  })

  output$plot_cluster_dend <- renderPlot({
    req(rv$clust)
    validate(need(!is.null(rv$clust$hc), "Dendrogram is available only for hierarchical clustering."))
    hc <- rv$clust$hc
    k  <- max(rv$clust$frames$cluster, na.rm = TRUE)

    # Dark-themed dendrogram with no leaf labels (unreadable at scale)
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    par(bg = "#1a2235",
        col.axis = "#8899aa", col.lab = "#8899aa", col.main = "#e8ecf1",
        fg = "#8899aa",
        mar = c(3, 4, 3, 1))

    # Convert to dendrogram and color branches by cluster
    dend <- as.dendrogram(hc)
    # Suppress individual leaf labels (they overlap when n > ~50)
    dend <- dendrapply(dend, function(n) {
      if (is.leaf(n)) attr(n, "label") <- ""
      n
    })
    plot(dend, main = "Hierarchical clustering dendrogram",
         ylab = "Height (RMSD)", xlab = "",
         leaflab = "none", edgePar = list(col = "#5a7a9a", lwd = 0.6))

    # Draw cluster rectangles
    palette_k <- c("#06d6a0","#ef476f","#ffd166","#118ab2","#073b4c",
                   "#8338ec","#ff006e","#fb5607","#3a86ff","#70e000")
    tryCatch(
      rect.hclust(hc, k = k, border = palette_k[seq_len(min(k, length(palette_k)))]),
      error = function(e) NULL
    )
  }, bg = "#1a2235")

  # ───────────────────────────────────────────────────────────────────────────
  # Pairwise RMSD heatmap (uses the distance matrix from clustering)
  # ───────────────────────────────────────────────────────────────────────────
  output$plot_rmsd_heatmap <- renderPlotly({
    req(rv$clust)
    validate(need(!is.null(rv$clust$d_rmsd), "RMSD distance matrix not available. Run clustering first."))
    d <- rv$clust$d_rmsd
    mat <- as.matrix(d)
    nf <- nrow(mat)

    # Use global time or frame index for axes
    fdf <- rv$clust$frames
    if (isTRUE(input$combine_trajs) && ("time_ns_global" %in% names(fdf)) && any(!is.na(fdf$time_ns_global))) {
      ax_vals <- fdf$time_ns_global
      ax_title <- "Time (ns)"
    } else {
      ax_vals <- seq_len(nf)
      ax_title <- "Frame index"
    }

    maxval <- suppressWarnings(as.numeric(input$rmsd_heatmap_maxval))
    if (!is.finite(maxval) || maxval <= 0) maxval <- max(mat, na.rm = TRUE)
    mat[mat > maxval] <- maxval

    plot_ly(z = mat, x = ax_vals, y = ax_vals, type = "heatmap",
            colorscale = list(c(0, "#0a0e17"), c(0.25, "#118ab2"), c(0.5, "#06d6a0"),
                              c(0.75, "#ffd166"), c(1, "#ef476f")),
            zmin = 0, zmax = maxval,
            colorbar = list(title = list(text = "RMSD (Å)", font = list(size = 12, color = "#e8ecf1")),
                          tickfont = list(color = "#c0ccda", size = 11)),
            hovertemplate = paste0("Frame i: %{x}<br>Frame j: %{y}<br>RMSD: %{z:.3f} Å<extra></extra>")) |>
      layout(xaxis = list(title = ax_title), yaxis = list(title = ax_title)) |>
      apply_plotly_theme()
  })

  output$dl_rmsdHeatmap_plot <- downloadHandler(
    filename = function() "pairwise_rmsd_heatmap.png",
    content = function(file) {
      req(rv$clust, rv$clust$d_rmsd)
      mat <- as.matrix(rv$clust$d_rmsd)
      maxval <- suppressWarnings(as.numeric(input$rmsd_heatmap_maxval))
      if (!is.finite(maxval) || maxval <= 0) maxval <- max(mat, na.rm = TRUE)
      mat[mat > maxval] <- maxval
      grDevices::png(file, width = 2400, height = 2200, res = 300)
      op <- par(no.readonly = TRUE); on.exit(par(op))
      par(mar = c(4, 4, 2, 5))
      image(mat, col = grDevices::colorRampPalette(c("#0a0e17","#118ab2","#06d6a0","#ffd166","#ef476f"))(100),
            zlim = c(0, maxval), xlab = "Frame", ylab = "Frame", main = "Pairwise RMSD (Å)", axes = FALSE)
      axis(1); axis(2)
      grDevices::dev.off()
    }
  )

  # ───────────────────────────────────────────────────────────────────────────
  # RMSD distribution histogram
  # ───────────────────────────────────────────────────────────────────────────
  output$plot_rmsd_dist <- renderPlotly({
    req(rv$data$rmsd_pep)
    st <- style_for("rmsdDist")
    df <- rv$data$rmsd_pep
    which_series <- input$rmsdA_series %||% "both"
    long_df <- rmsd_series_long(df, which_series)
    req(nrow(long_df) > 0)
    nbins <- as.integer(input$rmsdDist_nbins %||% 60)

    p <- plot_ly(alpha = 0.7)
    for (s in unique(long_df$series)) {
      vals <- long_df$rmsd_value[long_df$series == s]
      vals <- vals[is.finite(vals)]
      p <- p |> add_histogram(x = vals, name = s, nbinsx = nbins)
    }
    p |> layout(barmode = "overlay",
                xaxis = plotly_axis_from_style(st, "x", "RMSD (Å)"),
                yaxis = plotly_axis_from_style(st, "y", "Count")) |>
      apply_plotly_theme()
  })

  output$dl_rmsdDist_plot <- downloadHandler(
    filename = function() paste0("rmsd_distribution.", style_for("rmsdDist")$fmt %||% "png"),
    content = function(file) {
      req(rv$data$rmsd_pep)
      st <- style_for("rmsdDist")
      df <- rv$data$rmsd_pep
      which_series <- input$rmsdA_series %||% "both"
      long_df <- rmsd_series_long(df, which_series)
      validate(need(nrow(long_df) > 0, "No RMSD data"))
      nbins <- as.integer(input$rmsdDist_nbins %||% 60)
      p <- ggplot(long_df, aes(x = .data[["rmsd_value"]], fill = .data[["series"]])) +
        geom_histogram(bins = nbins, alpha = 0.7, position = "identity") +
        labs(x = axis_title_or(st$x_title, "RMSD (Å)"),
             y = axis_title_or(st$y_title, "Count"),
             fill = "Series", title = "RMSD distribution (Selection A)") +
        gg_theme_pub(st)
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )

  # ───────────────────────────────────────────────────────────────────────────
  # Summary stats cards (RMSD A / Rg A)
  # ───────────────────────────────────────────────────────────────────────────
  output$rmsdA_stats_row <- renderUI({
    req(rv$data$rmsd_pep)
    df <- rv$data$rmsd_pep
    bb <- if ("rmsd_backbone" %in% names(df)) df$rmsd_backbone else NULL
    hv <- if ("rmsd_heavy" %in% names(df)) df$rmsd_heavy else if ("rmsd_A" %in% names(df)) df$rmsd_A else NULL
    cards <- list()
    if (!is.null(hv) && any(is.finite(hv))) {
      m <- mean(hv, na.rm = TRUE); s <- sd(hv, na.rm = TRUE)
      cards <- c(cards, list(
        div(class = "metric-card",
            div(class = "metric-value", sprintf("%.2f ± %.2f", m, s)),
            div(class = "metric-label", "Heavy atoms RMSD (Å)")),
        div(class = "metric-card",
            div(class = "metric-value", sprintf("%.2f", max(hv, na.rm = TRUE))),
            div(class = "metric-label", "Max RMSD (Å)")),
        div(class = "metric-card",
            div(class = "metric-value", format(length(hv[is.finite(hv)]), big.mark = ",")),
            div(class = "metric-label", "Frames"))
      ))
    }
    if (!is.null(bb) && any(is.finite(bb))) {
      m <- mean(bb, na.rm = TRUE); s <- sd(bb, na.rm = TRUE)
      cards <- c(cards, list(
        div(class = "metric-card",
            div(class = "metric-value", sprintf("%.2f ± %.2f", m, s)),
            div(class = "metric-label", "Backbone RMSD (Å)"))
      ))
    }
    if (length(cards) > 0) div(class = "stat-row", cards) else NULL
  })

  output$rgA_stats_row <- renderUI({
    req(rv$data$rg_pep)
    df <- rv$data$rg_pep
    rg <- df$rg_A
    if (is.null(rg) || !any(is.finite(rg))) return(NULL)
    m <- mean(rg, na.rm = TRUE); s <- sd(rg, na.rm = TRUE)
    mn <- min(rg, na.rm = TRUE); mx <- max(rg, na.rm = TRUE)
    div(class = "stat-row",
      div(class = "metric-card",
          div(class = "metric-value", sprintf("%.2f ± %.2f", m, s)),
          div(class = "metric-label", "Mean ± SD (Å)")),
      div(class = "metric-card",
          div(class = "metric-value", sprintf("%.2f", mn)),
          div(class = "metric-label", "Min Rg (Å)")),
      div(class = "metric-card",
          div(class = "metric-value", sprintf("%.2f", mx)),
          div(class = "metric-label", "Max Rg (Å)")),
      div(class = "metric-card",
          div(class = "metric-value", sprintf("%.2f", mx - mn)),
          div(class = "metric-label", "Range (Å)"))
    )
  })

  # ───────────────────────────────────────────────────────────────────────────
  # Free Energy Landscape (FEL)
  # ───────────────────────────────────────────────────────────────────────────
  fel_data <- reactive({
    req(rv$data$pca_scores)
    df <- rv$data$pca_scores
    validate(need(all(c("PC1", "PC2") %in% names(df)), "PCA scores not available."))
    nbins <- as.integer(input$fel_nbins %||% 80)
    temp  <- as.numeric(input$fel_temp %||% 300)
    compute_fel(df$PC1, df$PC2, nbins = nbins, temp_K = temp)
  })

  output$plot_fel <- renderPlotly({
    fel <- fel_data()
    req(!is.null(fel) && nrow(fel) > 0)
    maxG <- suppressWarnings(as.numeric(input$fel_maxG))
    if (!is.finite(maxG) || maxG <= 0) maxG <- max(fel$G_kcal, na.rm = TRUE)
    fel$G_plot <- pmin(fel$G_kcal, maxG, na.rm = TRUE)

    plot_ly(fel, x = ~PC1, y = ~PC2, z = ~G_plot, type = "heatmap",
            colorscale = list(c(0, "#0a0e17"), c(0.15, "#073b4c"), c(0.3, "#118ab2"),
                              c(0.5, "#06d6a0"), c(0.7, "#ffd166"), c(1, "#ef476f")),
            zmin = 0, zmax = maxG,
            colorbar = list(title = list(text = "ΔG (kcal/mol)", font = list(size = 12, color = "#e8ecf1")),
                          tickfont = list(color = "#c0ccda", size = 11)),
            hovertemplate = "PC1: %{x:.2f}<br>PC2: %{y:.2f}<br>ΔG: %{z:.2f} kcal/mol<extra></extra>") |>
      layout(xaxis = list(title = "PC1", scaleanchor = "y"),
             yaxis = list(title = "PC2")) |>
      apply_plotly_theme()
  })

  output$dl_fel_data <- downloadHandler(
    filename = function() "free_energy_landscape.csv",
    content = function(file) {
      fel <- fel_data()
      req(!is.null(fel))
      write.csv(fel, file, row.names = FALSE)
    }
  )

  output$dl_fel_plot <- downloadHandler(
    filename = function() paste0("free_energy_landscape.", style_for("fel")$fmt %||% "png"),
    content = function(file) {
      fel <- fel_data()
      req(!is.null(fel) && nrow(fel) > 0)
      st <- style_for("fel")
      maxG <- suppressWarnings(as.numeric(input$fel_maxG))
      if (!is.finite(maxG) || maxG <= 0) maxG <- max(fel$G_kcal, na.rm = TRUE)
      fel$G_plot <- pmin(fel$G_kcal, maxG, na.rm = TRUE)
      p <- ggplot(fel[is.finite(fel$G_plot), ], aes(x = .data[["PC1"]], y = .data[["PC2"]], fill = .data[["G_plot"]])) +
        geom_tile() +
        scale_fill_gradientn(
          colours = c("#0a0e17","#073b4c","#118ab2","#06d6a0","#ffd166","#ef476f"),
          na.value = "transparent",
          name = expression(Delta*G~"(kcal/mol)"),
          limits = c(0, maxG)
        ) +
        labs(x = "PC1", y = "PC2", title = "Free energy landscape") +
        coord_fixed() +
        gg_theme_pub(st)
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )

  # ───────────────────────────────────────────────────────────────────────────
  # Session info (reproducibility)
  # ───────────────────────────────────────────────────────────────────────────
  output$session_info_text <- renderText({
    si <- sessionInfo()
    paste(capture.output(print(si)), collapse = "\n")
  })

  output$dl_session_info <- downloadHandler(
    filename = function() paste0("ShineMD_session_info_", Sys.Date(), ".txt"),
    content = function(file) {
      si <- sessionInfo()
      writeLines(capture.output(print(si)), con = file)
    }
  )

  output$pkg_versions_tbl <- renderDT({
    pkgs <- c("shiny","shinydashboard","shinyjs","shinyWidgets","plotly","ggplot2",
              "htmlwidgets","DT","shinyFiles","fs","bio3d","ncdf4","cluster","svglite","viridis")
    vers <- vapply(pkgs, function(p) {
      v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed")
      v
    }, character(1))
    datatable(data.frame(Package = pkgs, Version = vers, stringsAsFactors = FALSE),
              options = list(pageLength = 20, dom = "t"), rownames = FALSE)
  })

  output$env_info_panel <- renderUI({
    r_ver <- paste(R.version$major, R.version$minor, sep = ".")
    platform <- R.version$platform
    os <- Sys.info()[["sysname"]]
    tagList(
      tags$div(class = "metric-card",
        tags$div(class = "metric-value", r_ver),
        tags$div(class = "metric-label", "R version")
      ),
      tags$div(class = "metric-card",
        tags$div(class = "metric-value", style = "font-size: 16px;", platform),
        tags$div(class = "metric-label", "Platform")
      ),
      tags$div(class = "metric-card",
        tags$div(class = "metric-value", os),
        tags$div(class = "metric-label", "Operating system")
      ),
      tags$div(class = "metric-card",
        tags$div(class = "metric-value", format(Sys.time(), "%Y-%m-%d %H:%M")),
        tags$div(class = "metric-label", "Report generated")
      )
    )
  })

output$tbl_cluster_summary <- renderDT({
    req(rv$clust)
    datatable(rv$clust$summary, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$log <- renderText(rv$log)

  time_col <- reactive({
    if (isTRUE(input$combine_trajs)) "time_ns_global" else "time_ns"
  })

  # Per-plot styling/export options (each plot has its own controls in the UI)
  style_for <- function(prefix, defaults = list()) {
    list(
      theme       = input[[paste0(prefix, "_pub_theme")]] %||% defaults$theme %||% "classic",
      base_size   = input[[paste0(prefix, "_pub_base_size")]] %||% defaults$base_size %||% 12,
      line_width  = input[[paste0(prefix, "_line_width")]] %||% defaults$line_width %||% 1.2,
      show_points = isTRUE(input[[paste0(prefix, "_show_points")]] %||% defaults$show_points %||% FALSE),
      point_size  = input[[paste0(prefix, "_point_size")]] %||% defaults$point_size %||% 3,
      point_alpha = input[[paste0(prefix, "_point_alpha")]] %||% defaults$point_alpha %||% 0.8,
      palette     = input[[paste0(prefix, "_palette")]] %||% defaults$palette %||% "default",
      fmt         = tolower(input[[paste0(prefix, "_export_fmt")]] %||% defaults$fmt %||% "png"),
      width       = input[[paste0(prefix, "_export_width")]] %||% defaults$width %||% 7,
      height      = input[[paste0(prefix, "_export_height")]] %||% defaults$height %||% 4,
      dpi         = input[[paste0(prefix, "_export_dpi")]] %||% defaults$dpi %||% 300,
      x_title     = input[[paste0(prefix, "_x_title")]] %||% "",
      y_title     = input[[paste0(prefix, "_y_title")]] %||% "",
      x_scale     = input[[paste0(prefix, "_x_scale")]] %||% "linear",
      y_scale     = input[[paste0(prefix, "_y_scale")]] %||% "linear",
      x_min       = input[[paste0(prefix, "_x_min")]],
      x_max       = input[[paste0(prefix, "_x_max")]],
      y_min       = input[[paste0(prefix, "_y_min")]],
      y_max       = input[[paste0(prefix, "_y_max")]]

    )
  }

  

  # ───────────────────────────────────────────────────────────────────────────
  # Axis customization (interactive Plotly + exported ggplots)
  # ───────────────────────────────────────────────────────────────────────────
  axis_title_or <- function(custom, default) {
    custom <- as.character(custom %||% "")
    if (nzchar(trimws(custom))) trimws(custom) else default
  }

  axis_range_or <- function(minv, maxv) {
    if (is.null(minv) || is.null(maxv)) return(NULL)
    if (!is.finite(minv) || !is.finite(maxv)) return(NULL)
    c(as.numeric(minv), as.numeric(maxv))
  }

  # ggplot coord_cartesian limits: accepts single-sided (NA = use data range).
  # Returns NULL when no limit is set, so coord_cartesian is skipped entirely.
  gg_axis_lim <- function(minv, maxv, scale = "linear") {
    if (!identical(scale %||% "linear", "linear")) return(NULL)
    lo <- if (!is.null(minv) && is.finite(minv)) as.numeric(minv) else NA_real_
    hi <- if (!is.null(maxv) && is.finite(maxv)) as.numeric(maxv) else NA_real_
    if (is.na(lo) && is.na(hi)) return(NULL)
    c(lo, hi)
  }

  plotly_axis_from_style <- function(st, axis = c("x","y"), default_title = "") {
    axis <- match.arg(axis)
    title_use <- axis_title_or(if (axis == "x") st$x_title else st$y_title, default_title)
    type_use  <- if ((if (axis == "x") st$x_scale else st$y_scale) == "log") "log" else "linear"

    r <- if (axis == "x") axis_range_or(st$x_min, st$x_max) else axis_range_or(st$y_min, st$y_max)
    ax <- list(title = list(text = title_use, font = list(size = 13), standoff = 10),
               type = type_use)
    # Apply range only for linear axes (log range is in log10 space, too confusing for users)
    if (!is.null(r) && identical(type_use, "linear")) ax$range <- r
    ax
  }

# ───────────────────────────────────────────────────────────────────────────
  # Plotly dark theme (match the app UI; keep exports publication-grade)
  # ───────────────────────────────────────────────────────────────────────────
  plotly_darkify <- function(p) {
    # base layout
    p <- plotly::layout(
      p,
      paper_bgcolor = "rgba(0,0,0,0)",
      plot_bgcolor  = "rgba(0,0,0,0)",
      font = list(family = "IBM Plex Sans", size = 13, color = "#e8ecf1"),
      legend = list(bgcolor = "rgba(26,34,53,0.85)", font = list(color = "#e8ecf1", size = 12),
                    bordercolor = "rgba(42,58,82,0.5)", borderwidth = 1,
                    x = 1.02, xanchor = "left", y = 1, yanchor = "top"),
      margin = list(l = 70, r = 30, t = 40, b = 60),
      hoverlabel = list(
        bgcolor = "#1a2235",
        bordercolor = "#06d6a0",
        font = list(color = "#e8ecf1", size = 12)
      )
    )

    axis_style <- list(
      gridcolor     = "rgba(42,58,82,0.35)",
      zerolinecolor = "rgba(42,58,82,0.35)",
      linecolor     = "rgba(42,58,82,0.60)",
      tickcolor     = "rgba(42,58,82,0.60)",
      tickfont      = list(color = "#c0ccda", size = 11),
      titlefont     = list(color = "#e8ecf1", size = 13)
    )

    # 2D axes
    if (is.null(p$x$layout)) p$x$layout <- list()
    for (ax in c("xaxis", "yaxis", "yaxis2")) {
      if (!is.null(p$x$layout[[ax]]) || ax %in% c("xaxis", "yaxis")) {
        p$x$layout[[ax]] <- modifyList(p$x$layout[[ax]] %||% list(), axis_style)
      }
    }

    # 3D scene (if present)
    if (!is.null(p$x$layout$scene)) {
      p$x$layout$scene <- modifyList(p$x$layout$scene %||% list(), list(bgcolor = "rgba(0,0,0,0)"))
      for (ax in c("xaxis", "yaxis", "zaxis")) {
        p$x$layout$scene[[ax]] <- modifyList(p$x$layout$scene[[ax]] %||% list(), axis_style)
      }
    }

    plotly::config(p, displaylogo = FALSE, responsive = TRUE)
  }


  plotly_lightify <- function(p) {
    p <- plotly::layout(
      p,
      paper_bgcolor = "#ffffff",
      plot_bgcolor  = "#ffffff",
      font = list(family = "IBM Plex Sans", size = 13, color = "#111827"),
      legend = list(bgcolor = "rgba(255,255,255,0.90)", font = list(color = "#111827", size = 12),
                    bordercolor = "rgba(0,0,0,0.12)", borderwidth = 1,
                    x = 1.02, xanchor = "left", y = 1, yanchor = "top"),
      margin = list(l = 70, r = 30, t = 40, b = 60),
      hoverlabel = list(
        bgcolor = "#ffffff",
        bordercolor = "#118ab2",
        font = list(color = "#111827", size = 12)
      )
    )

    axis_style <- list(
      gridcolor     = "rgba(0,0,0,0.12)",
      zerolinecolor = "rgba(0,0,0,0.12)",
      linecolor     = "rgba(0,0,0,0.25)",
      tickcolor     = "rgba(0,0,0,0.25)",
      tickfont      = list(color = "#374151", size = 11),
      titlefont     = list(color = "#111827", size = 13)
    )

    if (is.null(p$x$layout)) p$x$layout <- list()
    for (ax in c("xaxis", "yaxis", "yaxis2")) {
      if (!is.null(p$x$layout[[ax]]) || ax %in% c("xaxis", "yaxis")) {
        p$x$layout[[ax]] <- modifyList(p$x$layout[[ax]] %||% list(), axis_style)
      }
    }

    if (!is.null(p$x$layout$scene)) {
      p$x$layout$scene <- modifyList(p$x$layout$scene %||% list(), list(bgcolor = "#ffffff"))
      for (ax in c("xaxis", "yaxis", "zaxis")) {
        p$x$layout$scene[[ax]] <- modifyList(p$x$layout$scene[[ax]] %||% list(), axis_style)
      }
    }

    plotly::config(p, displaylogo = FALSE, responsive = TRUE)
  }

  apply_plotly_theme <- function(p) {
    if (identical(rv$plotly_theme, "light")) plotly_lightify(p) else plotly_darkify(p)
  }


  selection_b_active <- reactive({
    has_input_b <- nzchar(trimws(input$lipo_resno %||% ""))
    has_rmsd_b <- !is.null(rv$data$rmsd_lipo) && is.data.frame(rv$data$rmsd_lipo) &&
      nrow(rv$data$rmsd_lipo) > 0 && ("rmsd_A" %in% names(rv$data$rmsd_lipo)) &&
      any(is.finite(rv$data$rmsd_lipo$rmsd_A))
    has_rg_b <- !is.null(rv$data$rg_lipo) && is.data.frame(rv$data$rg_lipo) &&
      nrow(rv$data$rg_lipo) > 0 && ("rg_A" %in% names(rv$data$rg_lipo)) &&
      any(is.finite(rv$data$rg_lipo$rg_A))
    isTRUE(has_input_b && (has_rmsd_b || has_rg_b))
  })

  output$ui_rmsd_B <- renderUI({
    if (!isTRUE(selection_b_active())) return(NULL)
    fluidRow(
      column(12,
             box(title=box_title_with_info("RMSD — Selection B", "info_rmsdB"), width=12, status="primary",
                 fluidRow(
                   column(8, plotlyOutput("plot_rmsd_lipo", height="360px")),
                   column(4,
                          plotly_theme_toggle_ui("toggle_plotly_theme_rmsdB", theme = rv$plotly_theme),
                          selectInput("rmsdB_series", "RMSD series", choices = c("Backbone + heavy atoms" = "both", "Backbone only" = "backbone", "Heavy atoms only" = "heavy"), selected = "both"),
                          numericInput("rmsdB_smooth", "Smoothing window (frames; 0 = off)", value = 0, min = 0, step = 10),
                          plot_style_export_ui("rmsdB", show_points_default = FALSE),
                          tags$br(),
                          downloadButton("dl_rmsdB_plot", "Download plot", class="btn-primary"),
                          downloadButton("dl_rmsdB_data", "Download plot data (CSV)", class="btn-default")
                   )
                 )
             )
      )
    )
  })

  output$ui_rg_B <- renderUI({
    if (!isTRUE(selection_b_active())) return(NULL)
    fluidRow(
      column(12,
             box(title=box_title_with_info("Radius of gyration — Selection B", "info_rgB"), width=12, status="danger",
                 fluidRow(
                   column(8, plotlyOutput("plot_rg_lipo", height="360px")),
                   column(4,
                          plotly_theme_toggle_ui("toggle_plotly_theme_rgB", theme = rv$plotly_theme),
                          plot_style_export_ui("rgB", show_points_default = FALSE),
                          tags$br(),
                          downloadButton("dl_rgB_plot", "Download plot", class="btn-primary"),
                          downloadButton("dl_rgB_data", "Download plot data (CSV)", class="btn-default")
                   )
                 )
             )
      )
    )
  })

    rmsd_series_long <- function(df, which = c("both", "backbone", "heavy")) {
    which <- match.arg(which)
    if (!is.data.frame(df) || nrow(df) < 1) return(data.frame())

    has_backbone <- "rmsd_backbone" %in% names(df) && any(is.finite(df$rmsd_backbone))
    has_heavy <- ("rmsd_heavy" %in% names(df) && any(is.finite(df$rmsd_heavy))) ||
      ("rmsd_total" %in% names(df) && any(is.finite(df$rmsd_total))) ||
      ("rmsd_A" %in% names(df) && any(is.finite(df$rmsd_A)))
    heavy_vec <- if ("rmsd_heavy" %in% names(df)) df$rmsd_heavy else if ("rmsd_total" %in% names(df)) df$rmsd_total else df$rmsd_A

    keep <- switch(which,
                   backbone = c("Backbone"),
                   heavy = c("Heavy atoms"),
                   both = c("Backbone", "Heavy atoms"))

    out <- list()
    if (has_backbone && "Backbone" %in% keep) {
      tmp <- df
      tmp$series <- "Backbone"
      tmp$rmsd_value <- tmp$rmsd_backbone
      out[[length(out) + 1L]] <- tmp
    }
    if (has_heavy && "Heavy atoms" %in% keep) {
      tmp <- df
      tmp$series <- "Heavy atoms"
      tmp$rmsd_value <- heavy_vec
      out[[length(out) + 1L]] <- tmp
    }
    if (length(out) < 1) return(data.frame())
    do.call(rbind, out)
  }

  make_rmsd_overlay_gg <- function(df, which = c("both", "backbone", "heavy"), title = "RMSD", style = NULL, combined = FALSE, smooth_k = 0) {
    which <- match.arg(which)
    style <- style %||% default_plot_style()
    long_df <- rmsd_series_long(df, which)
    req(nrow(long_df) > 0)

    xcol <- if (isTRUE(combined) && "time_ns_global" %in% names(long_df)) "time_ns_global" else if (time_col() %in% names(long_df)) time_col() else "time_ns"
    long_df <- long_df[order(long_df[[xcol]]), , drop = FALSE]

    p <- ggplot(long_df, aes(x = .data[[xcol]], y = .data[["rmsd_value"]], color = .data[["series"]], group = .data[["series"]])) +
      geom_line(linewidth = style$line_width %||% 0.6, alpha = style$line_alpha %||% 1)

    if (isTRUE(style$show_points)) {
      p <- p + geom_point(size = style$point_size %||% 1.4, alpha = style$point_alpha %||% 0.8)
    }

    smooth_k <- as.integer(smooth_k %||% 0)
    if (smooth_k >= 2) {
      sm_rows <- do.call(rbind, lapply(unique(long_df$series), function(s) {
        sub <- long_df[long_df$series == s, , drop = FALSE]
        sm  <- running_mean(sub$rmsd_value, smooth_k)
        data.frame(xcol_val = sub[[xcol]], rmsd_value = sm, series = paste0(s, " (avg ", smooth_k, ")"), stringsAsFactors = FALSE)
      }))
      sm_rows <- sm_rows[is.finite(sm_rows$rmsd_value), , drop = FALSE]
      if (nrow(sm_rows) > 0) {
        p <- p + geom_line(data = sm_rows,
                           aes(x = xcol_val, y = rmsd_value, color = series, group = series),
                           linewidth = (style$line_width %||% 0.6) + 0.6,
                           inherit.aes = FALSE)
      }
    }

    xlim_gg <- gg_axis_lim(style$x_min, style$x_max, style$x_scale)
    ylim_gg <- gg_axis_lim(style$y_min, style$y_max, style$y_scale)
    if (!is.null(xlim_gg) || !is.null(ylim_gg))
      p <- p + coord_cartesian(xlim = xlim_gg, ylim = ylim_gg)

    p <- p + labs(
      title = title,
      x = axis_title_or(style$x_title, if (isTRUE(combined)) "Time (ns)" else "Time (ns)"),
      y = axis_title_or(style$y_title, "RMSD (Å)"),
      color = "Series"
    ) + gg_theme_pub(style)
    apply_palette(p, style)
  }

# RMSD plots
  output$plot_rmsd_pep <- renderPlotly({
    req(rv$data$rmsd_pep)
    st <- style_for("rmsdA")
    df <- rv$data$rmsd_pep
    which_series <- input$rmsdA_series %||% "both"
    long_df <- rmsd_series_long(df, which_series)
    req(nrow(long_df) > 0)
    xcol <- if (isTRUE(input$combine_trajs) && ("time_ns_global" %in% names(long_df))) "time_ns_global" else if (time_col() %in% names(long_df)) time_col() else "time_ns"
    long_df <- long_df[order(long_df[[xcol]]), , drop = FALSE]

    mode_use <- if (isTRUE(st$show_points)) "lines+markers" else "lines"
    p <- plot_ly(long_df, x = long_df[[xcol]], y = ~rmsd_value, color = ~series, type = "scatter", mode = mode_use,
                 line = list(width = st$line_width),
                 marker = list(size = st$point_size, opacity = st$point_alpha),
                 hovertemplate = paste0("Series: %{fullData.name}<br>", if (identical(xcol, "time_ns_global")) "Time (ns): %{x:.3f}<br>" else "Time (ns): %{x:.3f}<br>", "RMSD (Å): %{y:.3f}<extra></extra>"))

    # Smoothing overlay — each series gets its own independent trace
    smooth_k <- as.integer(input$rmsdA_smooth %||% 0)
    if (smooth_k >= 2) {
      for (s in unique(long_df$series)) {
        sub <- long_df[long_df$series == s, , drop = FALSE]
        sm <- running_mean(sub$rmsd_value, smooth_k)
        sm_df <- data.frame(xv = sub[[xcol]], yv = sm, stringsAsFactors = FALSE)
        sm_df <- sm_df[is.finite(sm_df$yv), , drop = FALSE]
        if (nrow(sm_df) > 0) {
          p <- p |> add_trace(data = sm_df, x = ~xv, y = ~yv, type = "scatter", mode = "lines",
                              inherit = FALSE,
                              name = paste0(s, " (avg ", smooth_k, ")"),
                              line = list(width = st$line_width + 1.2),
                              hovertemplate = paste0(s, " avg<br>Time: %{x:.3f}<br>RMSD: %{y:.3f}<extra></extra>"))
        }
      }
    }

    p <- p |> layout(xaxis = plotly_axis_from_style(st, "x", if (identical(xcol, "time_ns_global")) "Time (ns)" else "Time (ns)"),
                     yaxis = plotly_axis_from_style(st, "y", "RMSD (Å)"))
    apply_plotly_theme(p)
  })
  output$plot_rmsd_lipo <- renderPlotly({
    req(rv$data$rmsd_lipo)
    st <- style_for("rmsdB")
    df <- rv$data$rmsd_lipo
    which_series <- input$rmsdB_series %||% "both"
    long_df <- rmsd_series_long(df, which_series)
    req(nrow(long_df) > 0)
    xcol <- if (isTRUE(input$combine_trajs) && ("time_ns_global" %in% names(long_df))) "time_ns_global" else if (time_col() %in% names(long_df)) time_col() else "time_ns"
    long_df <- long_df[order(long_df[[xcol]]), , drop = FALSE]

    mode_use <- if (isTRUE(st$show_points)) "lines+markers" else "lines"
    p <- plot_ly(long_df, x = long_df[[xcol]], y = ~rmsd_value, color = ~series, type = "scatter", mode = mode_use,
                 line = list(width = st$line_width),
                 marker = list(size = st$point_size, opacity = st$point_alpha),
                 hovertemplate = paste0("Series: %{fullData.name}<br>", if (identical(xcol, "time_ns_global")) "Time (ns): %{x:.3f}<br>" else "Time (ns): %{x:.3f}<br>", "RMSD (Å): %{y:.3f}<extra></extra>"))

    smooth_k <- as.integer(input$rmsdB_smooth %||% 0)
    if (smooth_k >= 2) {
      for (s in unique(long_df$series)) {
        sub <- long_df[long_df$series == s, , drop = FALSE]
        sm <- running_mean(sub$rmsd_value, smooth_k)
        sm_df <- data.frame(xv = sub[[xcol]], yv = sm, stringsAsFactors = FALSE)
        sm_df <- sm_df[is.finite(sm_df$yv), , drop = FALSE]
        if (nrow(sm_df) > 0) {
          p <- p |> add_trace(data = sm_df, x = ~xv, y = ~yv, type = "scatter", mode = "lines",
                              inherit = FALSE,
                              name = paste0(s, " (avg ", smooth_k, ")"),
                              line = list(width = st$line_width + 1.2),
                              hovertemplate = paste0(s, " avg<br>Time: %{x:.3f}<br>RMSD: %{y:.3f}<extra></extra>"))
        }
      }
    }

    p <- p |> layout(xaxis = plotly_axis_from_style(st, "x", if (identical(xcol, "time_ns_global")) "Time (ns)" else "Time (ns)"),
                     yaxis = plotly_axis_from_style(st, "y", "RMSD (Å)"))
    apply_plotly_theme(p)
  })
# RMSF
  output$plot_rmsf <- renderPlotly({
    req(rv$data$rmsf)
    st <- style_for("rmsf", defaults = list(show_points = TRUE))
    df <- rv$data$rmsf
    mode_use <- if (isTRUE(st$show_points)) "lines+markers" else "lines"
    plot_ly(df, x=~resno, y=~rmsf_A, type="scatter", mode = mode_use,
            line = list(width = st$line_width),
            marker = list(size = st$point_size, opacity = st$point_alpha)) |>
      layout(xaxis=plotly_axis_from_style(st, "x", "Residue"), yaxis=plotly_axis_from_style(st, "y", "RMSF (Å)")) |> apply_plotly_theme()
  })

  output$ui_rmsf_B <- renderUI({
    if (!isTRUE(selection_b_active())) return(NULL)
    fluidRow(
      column(12,
             box(title=box_title_with_info("RMSF (Selection B, per residue)", "info_rmsfB"), width=12, status="warning",
                 fluidRow(
                   column(8, plotlyOutput("plot_rmsf_B", height="450px")),
                   column(4,
                          plotly_theme_toggle_ui("toggle_plotly_theme_rmsfB", theme = rv$plotly_theme),
                          plot_style_export_ui("rmsfB", show_points_default = TRUE, line_width_default = 1.0),
                          tags$br(),
                          downloadButton("dl_rmsfB_plot", "Download plot", class="btn-primary"),
                          downloadButton("dl_rmsfB_data", "Download plot data (CSV)", class="btn-default")
                   )
                 )
             )
      )
    )
  })

  output$plot_rmsf_B <- renderPlotly({
    req(rv$data$rmsf_B)
    st <- style_for("rmsfB", defaults = list(show_points = TRUE))
    df <- rv$data$rmsf_B
    mode_use <- if (isTRUE(st$show_points)) "lines+markers" else "lines"
    plot_ly(df, x=~resno, y=~rmsf_B, type="scatter", mode = mode_use,
            line = list(width = st$line_width),
            marker = list(size = st$point_size, opacity = st$point_alpha)) |>
      layout(xaxis=plotly_axis_from_style(st, "x", "Residue"), yaxis=plotly_axis_from_style(st, "y", "RMSF (Å)")) |> apply_plotly_theme()
  })
# Rg
  output$plot_rg_pep <- renderPlotly({
    req(rv$data$rg_pep)
    st <- style_for("rgA")
    df <- rv$data$rg_pep
    xcol <- time_col()
    if (!(xcol %in% names(df))) xcol <- "time_ns"

    mode_use <- if (isTRUE(st$show_points)) "lines+markers" else "lines"
    smooth_k <- as.integer(input$rgA_smooth %||% 0)

    if (isTRUE(input$combine_trajs) && ("time_ns_global" %in% names(df))) {
      df <- df[order(df$time_ns_global), ]
      p <- plot_ly(df, x = ~time_ns_global, y = ~rg_A, type="scatter", mode = mode_use,
              name = "Rg",
              line = list(width = st$line_width),
              marker = list(size = st$point_size, opacity = st$point_alpha))
      if (smooth_k >= 2) {
        sm <- running_mean(df$rg_A, smooth_k)
        sm_df <- data.frame(xv = df$time_ns_global, yv = sm, stringsAsFactors = FALSE)
        sm_df <- sm_df[is.finite(sm_df$yv), , drop = FALSE]
        if (nrow(sm_df) > 0) {
          p <- p |> add_trace(data = sm_df, x = ~xv, y = ~yv, type = "scatter", mode = "lines",
                              inherit = FALSE,
                              name = paste0("Avg (", smooth_k, ")"),
                              line = list(width = st$line_width + 1.2))
        }
      }
      p |> layout(xaxis=plotly_axis_from_style(st, "x", "Time (ns)"),
               yaxis=plotly_axis_from_style(st, "y", "Rg (Å)")) |>
      apply_plotly_theme()
    } else {
      p <- plot_ly(df, x = df[[xcol]], y = ~rg_A, color = ~segment, type="scatter", mode = mode_use,
              line = list(width = st$line_width),
              marker = list(size = st$point_size, opacity = st$point_alpha))
      if (smooth_k >= 2) {
        sm <- running_mean(df$rg_A, smooth_k)
        sm_df <- data.frame(xv = df[[xcol]], yv = sm, stringsAsFactors = FALSE)
        sm_df <- sm_df[is.finite(sm_df$yv), , drop = FALSE]
        if (nrow(sm_df) > 0) {
          p <- p |> add_trace(data = sm_df, x = ~xv, y = ~yv, type = "scatter", mode = "lines",
                              inherit = FALSE,
                              name = paste0("Avg (", smooth_k, ")"),
                              line = list(width = st$line_width + 1.2, color = "#ef476f"))
        }
      }
      p |> layout(xaxis=plotly_axis_from_style(st, "x", "Time (ns)"),
               yaxis=plotly_axis_from_style(st, "y", "Rg (Å)")) |> apply_plotly_theme()
    }
  })
  output$plot_rg_lipo <- renderPlotly({
    req(rv$data$rg_lipo)
    st <- style_for("rgB")
    df <- rv$data$rg_lipo
    xcol <- time_col()
    if (!(xcol %in% names(df))) xcol <- "time_ns"

    mode_use <- if (isTRUE(st$show_points)) "lines+markers" else "lines"
    if (isTRUE(input$combine_trajs) && ("time_ns_global" %in% names(df))) {
      df <- df[order(df$time_ns_global), ]
      plot_ly(df, x = ~time_ns_global, y = ~rg_A, type="scatter", mode = mode_use,
              line = list(width = st$line_width),
              marker = list(size = st$point_size, opacity = st$point_alpha)) |>
        layout(xaxis=plotly_axis_from_style(st, "x", "Time (ns)"),
               yaxis=plotly_axis_from_style(st, "y", "Rg (Å)")) |>
      apply_plotly_theme()
    } else {
      plot_ly(df, x = df[[xcol]], y = ~rg_A, color = ~segment, type="scatter", mode = mode_use,
              line = list(width = st$line_width),
              marker = list(size = st$point_size, opacity = st$point_alpha)) |>
        layout(xaxis=plotly_axis_from_style(st, "x", "Time (ns)"),
               yaxis=plotly_axis_from_style(st, "y", "Rg (Å)")) |> apply_plotly_theme()
    }
  })
# Membrane metrics
  output$plot_dist <- renderPlotly({
    req(rv$data$dist)
    st <- style_for("memDist")
    df <- rv$data$dist
    xcol <- time_col()
    if (!(xcol %in% names(df))) xcol <- "time_ns"

    mode_use <- if (isTRUE(st$show_points)) "lines+markers" else "lines"
    if (isTRUE(input$combine_trajs) && ("time_ns_global" %in% names(df))) {
      df <- df[order(df$time_ns_global), ]
      plot_ly(df, x = ~time_ns_global, y = ~dist_A, type="scatter", mode = mode_use,
              line = list(width = st$line_width),
              marker = list(size = st$point_size, opacity = st$point_alpha)) |>
        layout(xaxis=plotly_axis_from_style(st, "x", "Time (ns)"),
               yaxis=plotly_axis_from_style(st, "y", "Distance (Å)")) |>
      apply_plotly_theme()
    } else {
      plot_ly(df, x = df[[xcol]], y = ~dist_A, color = ~segment, type="scatter", mode = mode_use,
              line = list(width = st$line_width),
              marker = list(size = st$point_size, opacity = st$point_alpha)) |>
        layout(xaxis=plotly_axis_from_style(st, "x", "Time (ns)"),
               yaxis=plotly_axis_from_style(st, "y", "Distance (Å)")) |> apply_plotly_theme()
    }
  })
  output$plot_z <- renderPlotly({
    req(rv$data$z)
    st <- style_for("memZ")
    df <- rv$data$z
    xcol <- time_col()
    if (!(xcol %in% names(df))) xcol <- "time_ns"

    mode_use <- if (isTRUE(st$show_points)) "lines+markers" else "lines"
    if (isTRUE(input$combine_trajs) && ("time_ns_global" %in% names(df))) {
      df <- df[order(df$time_ns_global), ]
      plot_ly(df, x = ~time_ns_global, y = ~z_A, type="scatter", mode = mode_use,
              line = list(width = st$line_width),
              marker = list(size = st$point_size, opacity = st$point_alpha)) |>
        layout(xaxis=plotly_axis_from_style(st, "x", "Time (ns)"),
               yaxis=plotly_axis_from_style(st, "y", "Δz (Å)")) |>
      apply_plotly_theme()
    } else {
      plot_ly(df, x = df[[xcol]], y = ~z_A, color = ~segment, type="scatter", mode = mode_use,
              line = list(width = st$line_width),
              marker = list(size = st$point_size, opacity = st$point_alpha)) |>
        layout(xaxis=plotly_axis_from_style(st, "x", "Time (ns)"),
               yaxis=plotly_axis_from_style(st, "y", "Δz (Å)")) |> apply_plotly_theme()
    }
  })


  # Advanced membrane plots
  output$apl_info <- renderText({
    req(rv$prmtop)
    mem_resid <- parse_csv_tokens(input$mem_resid)
    if (length(mem_resid) < 1) return("Membrane residue names not set.")

    prm <- rv$prm_obj
    if (is.null(prm)) {
      prm <- tryCatch(read.prmtop(rv$prmtop), error = function(e) NULL)
    }
    if (is.null(prm)) return("Auto lipid count: NA (failed to read prmtop).")

    n_total <- detect_total_lipids(prm, mem_resid)
    if (is.na(n_total)) return("Auto lipid count: NA (cannot detect from topology).")
    paste0("Auto lipid count (total): ", n_total, "  → per leaflet (approx): ", as.integer(round(n_total/2)))
  })

  output$plot_thickness <- renderPlotly({
    req(rv$data$thickness)
    st <- style_for("memThick")
    df <- rv$data$thickness
    xcol <- time_col()
    if (!(xcol %in% names(df))) xcol <- "time_ns"

    mode_use <- if (isTRUE(st$show_points)) "lines+markers" else "lines"
    if (isTRUE(input$combine_trajs) && ("time_ns_global" %in% names(df))) {
      df <- df[order(df$time_ns_global), ]
      plot_ly(df, x = ~time_ns_global, y = ~thickness_A, type="scatter", mode = mode_use,
              line = list(width = st$line_width),
              marker = list(size = st$point_size, opacity = st$point_alpha)) |>
        layout(xaxis=plotly_axis_from_style(st, "x", "Time (ns)"),
               yaxis=plotly_axis_from_style(st, "y", "Thickness (Å)")) |>
        apply_plotly_theme()
    } else {
      plot_ly(df, x = df[[xcol]], y = ~thickness_A, color = ~segment, type="scatter", mode = mode_use,
              line = list(width = st$line_width),
              marker = list(size = st$point_size, opacity = st$point_alpha)) |>
        layout(xaxis=plotly_axis_from_style(st, "x", "Time (ns)"),
               yaxis=plotly_axis_from_style(st, "y", "Thickness (Å)")) |> apply_plotly_theme()
    }
  })

  output$plot_apl <- renderPlotly({
    req(rv$data$apl)
    st <- style_for("memAPL")
    df <- rv$data$apl
    xcol <- time_col()
    if (!(xcol %in% names(df))) xcol <- "time_ns"

    mode_use <- if (isTRUE(st$show_points)) "lines+markers" else "lines"
    if (isTRUE(input$combine_trajs) && ("time_ns_global" %in% names(df))) {
      df <- df[order(df$time_ns_global), ]
      plot_ly(df, x = ~time_ns_global, y = ~APL_A2, type="scatter", mode = mode_use,
              line = list(width = st$line_width),
              marker = list(size = st$point_size, opacity = st$point_alpha)) |>
        layout(xaxis=plotly_axis_from_style(st, "x", "Time (ns)"),
               yaxis=plotly_axis_from_style(st, "y", "APL (Å² / lipid, per leaflet)")) |>
        apply_plotly_theme()
    } else {
      plot_ly(df, x = df[[xcol]], y = ~APL_A2, color = ~segment, type="scatter", mode = mode_use,
              line = list(width = st$line_width),
              marker = list(size = st$point_size, opacity = st$point_alpha)) |>
        layout(xaxis=plotly_axis_from_style(st, "x", "Time (ns)"),
               yaxis=plotly_axis_from_style(st, "y", "APL (Å² / lipid, per leaflet)")) |> apply_plotly_theme()
    }
  })

  output$plot_enrich <- renderPlotly({
    req(rv$data$enrich)
    st <- style_for("memEnrich")
    df <- rv$data$enrich
    xcol <- if (isTRUE(input$combine_trajs) && ("time_ns_global" %in% names(df))) "time_ns_global" else "time_ns"
    df <- df[order(df[[xcol]]), ]

    mode_use <- if (isTRUE(st$show_points)) "lines+markers" else "lines"
    plot_ly(df, x = df[[xcol]], y = ~count, color = ~lipid_type, type="scatter", mode = mode_use,
            line = list(width = st$line_width),
            marker = list(size = st$point_size, opacity = st$point_alpha)) |>
      layout(xaxis=plotly_axis_from_style(st, "x", if (xcol=="time_ns_global") "Time (ns)" else "Time (ns)"),
             yaxis=plotly_axis_from_style(st, "y", paste0("Headgroups within ", input$enrich_cutoff, " Å"))) |>
      apply_plotly_theme()
  })

  output$plot_mem_density <- renderPlotly({
    req(rv$data$mem_density)
    st <- style_for("memDensity")
    df <- rv$data$mem_density
    excl <- input$dens_exclude_groups
    if (length(excl) > 0) df <- df[!(df$group %in% excl), , drop = FALSE]
    mode_use <- if (isTRUE(st$show_points)) "lines+markers" else "lines"
    plot_ly(df, x = ~z_mid, y = ~density, color = ~group, type = "scatter", mode = mode_use,
            line = list(width = st$line_width),
            marker = list(size = st$point_size, opacity = st$point_alpha)) |>
      layout(xaxis = plotly_axis_from_style(st, "x", "z (Å, membrane-centered)"),
             yaxis = plotly_axis_from_style(st, "y", "Average number density (atoms / Å³)")) |>
      apply_plotly_theme()
  })

  output$plot_mem_order <- renderPlotly({
    req(rv$data$mem_order)
    st <- style_for("memOrder")
    df <- rv$data$mem_order
    
    df$series <- if ("series" %in% names(df)) {
      df$series
    } else if ("lipid_type" %in% names(df)) {
      paste(df$lipid_type, df$chain, sep = " — ")
    } else {
      df$chain
    }
    
    df <- df[is.finite(df$absS) & !is.na(df$series), , drop = FALSE]
    df$series <- droplevels(factor(as.character(df$series)))
    
    mode_use <- if (isTRUE(st$show_points)) "lines+markers" else "lines"
    
    plot_ly(
      df,
      x = ~segment_index, y = ~absS, color = ~series, text = ~segment_label,
      type = "scatter", mode = mode_use,
      line = list(width = st$line_width),
      marker = list(size = st$point_size, opacity = st$point_alpha)
    ) |>
      layout(
        xaxis = plotly_axis_from_style(st, "x", "Tail segment index"),
        yaxis = plotly_axis_from_style(st, "y", "|S|")
      ) |>
      apply_plotly_theme()
  })

  output$plot_interact_ts <- renderPlotly({
    req(rv$data$interact_ts)
    st <- style_for("intTs")
    df <- rv$data$interact_ts
    metric <- input$interact_metric %||% "min_dist_A"
    xcol <- if (isTRUE(input$combine_trajs) && ("time_ns_global" %in% names(df))) "time_ns_global" else "time_ns"
    ylab <- switch(metric,
                   min_dist_A = "Minimum distance (Å)",
                   com_dist_A = "COM distance (Å)",
                   atom_contacts = paste0("Atom contacts < ", input$interact_cutoff, " Å"),
                   metric)
    mode_use <- if (isTRUE(st$show_points)) "lines+markers" else "lines"
    plot_ly(df, x = df[[xcol]], y = df[[metric]], type = "scatter", mode = mode_use,
            line = list(width = st$line_width),
            marker = list(size = st$point_size, opacity = st$point_alpha)) |>
      layout(xaxis = plotly_axis_from_style(st, "x", if (xcol == "time_ns_global") "Time (ns)" else "Time (ns)"),
             yaxis = plotly_axis_from_style(st, "y", ylab)) |>
      apply_plotly_theme()
  })

  output$plot_interact_occ <- renderPlotly({
    req(rv$data$interact_occ)
    st <- style_for("intOcc")
    df <- rv$data$interact_occ
    side_use <- input$interact_occ_side %||% "Selection A"
    topn <- as.integer(input$interact_occ_topn %||% 20)
    df <- df[df$side == side_use & is.finite(df$occupancy_pct), , drop = FALSE]
    validate(need(nrow(df) > 0, paste0("No occupancy rows found for ", side_use, ".")))
    df <- head(df[order(-df$occupancy_pct, df$residue), , drop = FALSE], topn)
    df <- df[order(df$occupancy_pct, df$residue), , drop = FALSE]
    df$residue_plot <- factor(df$residue, levels = df$residue)
    p <- ggplot(df, aes(x = occupancy_pct, y = residue_plot)) +
      geom_col(width = 0.75) +
      labs(x = axis_title_or(st$x_title, "Occupancy (%)"),
           y = axis_title_or(st$y_title, "Residue"),
           title = paste0(side_use, " residue contact occupancy")) +
      gg_theme_pub(st)
    ggplotly(p, tooltip = c("y", "x")) |>
      layout(xaxis = plotly_axis_from_style(st, "x", axis_title_or(st$x_title, "Occupancy (%)")),
             yaxis = plotly_axis_from_style(st, "y", axis_title_or(st$y_title, "Residue")),
             showlegend = FALSE) |>
      apply_plotly_theme()
  })

  output$tbl_interact_occ <- renderDT({
    req(rv$data$interact_occ)
    datatable(rv$data$interact_occ, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$plot_interact_pair <- renderPlotly({
    req(rv$data$interact_pairs)
    st <- style_for("intPair")
    df <- rv$data$interact_pairs
    df <- df[is.finite(df$occupancy_pct) & df$occupancy_pct >= as.numeric(input$interact_pair_minpct %||% 5), , drop = FALSE]
    validate(need(nrow(df) > 0, "No residue-pair contacts passed the current occupancy threshold."))

    topA <- as.integer(input$interact_pair_topA %||% 20)
    topB <- as.integer(input$interact_pair_topB %||% 12)
    rankA <- aggregate(occupancy_pct ~ residue_A, data = df, FUN = sum)
    rankB <- aggregate(occupancy_pct ~ residue_B, data = df, FUN = sum)
    rankA <- rankA[order(-rankA$occupancy_pct, rankA$residue_A), , drop = FALSE]
    rankB <- rankB[order(-rankB$occupancy_pct, rankB$residue_B), , drop = FALSE]
    keepA <- head(rankA$residue_A, topA)
    keepB <- head(rankB$residue_B, topB)
    df <- df[df$residue_A %in% keepA & df$residue_B %in% keepB, , drop = FALSE]
    validate(need(nrow(df) > 0, "No residue pairs remain after applying the top-residue filters."))

    keepA <- rankA$residue_A[rankA$residue_A %in% unique(df$residue_A)]
    keepB <- rankB$residue_B[rankB$residue_B %in% unique(df$residue_B)]
    grid <- expand.grid(residue_A = keepA, residue_B = keepB, stringsAsFactors = FALSE)
    matdf <- merge(grid, df[, c("residue_A","residue_B","occupancy_pct")], by = c("residue_A","residue_B"), all.x = TRUE, sort = FALSE)
    matdf$occupancy_pct[!is.finite(matdf$occupancy_pct)] <- 0
    matdf$residue_A <- factor(matdf$residue_A, levels = rev(keepA))
    matdf$residue_B <- factor(matdf$residue_B, levels = keepB)
    matdf$tip <- paste0("Selection A: ", matdf$residue_A,
                        "<br>Selection B: ", matdf$residue_B,
                        "<br>Pair occupancy: ", sprintf("%.1f", matdf$occupancy_pct), "%")

    p <- ggplot(matdf, aes(x = residue_B, y = residue_A, fill = occupancy_pct, text = tip)) +
      geom_tile(color = "grey80", linewidth = 0.2) +
      labs(x = axis_title_or(st$x_title, "Selection B residue"),
           y = axis_title_or(st$y_title, "Selection A residue"),
           title = "A / B residue contact map") +
      gg_theme_pub(st) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid = element_blank()) +
      scale_fill_gradient(low = "#f0f4fa", high = "#11cbb7", name = "Occupancy (%)")

    ggplotly(p, tooltip = "text") |>
      layout(xaxis = plotly_axis_from_style(st, "x", axis_title_or(st$x_title, "Selection B residue")),
             yaxis = plotly_axis_from_style(st, "y", axis_title_or(st$y_title, "Selection A residue"))) |>
      apply_plotly_theme()
  })

  output$tbl_interact_pairs <- renderDT({
    req(rv$data$interact_pairs)
    datatable(rv$data$interact_pairs, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$plot_hbond_ts <- renderPlotly({
    req(rv$data$hbond_ts)
    st <- style_for("hbTs")
    df <- rv$data$hbond_ts
    metric <- input$hbond_metric %||% "hbond_count"
    xcol <- if (isTRUE(input$combine_trajs) && ("time_ns_global" %in% names(df))) "time_ns_global" else "time_ns"
    ylab <- switch(metric,
                   hbond_count = paste0("H-bond proxies < ", input$hbond_cutoff, " Å"),
                   hbond_count_AtoB = "Selection A donor → B acceptor",
                   hbond_count_BtoA = "Selection B donor → A acceptor",
                   min_da_dist = "Minimum donor–acceptor distance (Å)",
                   metric)
    mode_use <- if (isTRUE(st$show_points)) "lines+markers" else "lines"
    plot_ly(df, x = df[[xcol]], y = df[[metric]], type = "scatter", mode = mode_use,
            line = list(width = st$line_width),
            marker = list(size = st$point_size, opacity = st$point_alpha)) |>
      layout(xaxis = plotly_axis_from_style(st, "x", if (xcol == "time_ns_global") "Time (ns)" else "Time (ns)"),
             yaxis = plotly_axis_from_style(st, "y", ylab)) |>
      apply_plotly_theme()
  })

  output$plot_hbond_pair <- renderPlotly({
    req(rv$data$hbond_pairs)
    st <- style_for("hbPair")
    df <- rv$data$hbond_pairs
    dir_use <- input$hbond_pair_direction %||% "all"
    if (!identical(dir_use, "all")) df <- df[df$direction == dir_use, , drop = FALSE]
    df <- df[is.finite(df$occupancy_pct) & df$occupancy_pct >= as.numeric(input$hbond_pair_minpct %||% 5), , drop = FALSE]
    validate(need(nrow(df) > 0, "No H-bond proxy pairs passed the current filters."))
    topn <- as.integer(input$hbond_pair_topn %||% 20)
    df <- head(df[order(-df$occupancy_pct, df$pair_label), , drop = FALSE], topn)
    df <- df[order(df$occupancy_pct, df$pair_label), , drop = FALSE]
    df$pair_plot <- factor(df$pair_label, levels = df$pair_label)
    p <- ggplot(df, aes(x = occupancy_pct, y = pair_plot, text = paste0(direction, "<br>", pair_label, "<br>Occupancy: ", sprintf("%.1f", occupancy_pct), "%"))) +
      geom_col(width = 0.75) +
      labs(x = axis_title_or(st$x_title, "Occupancy (%)"),
           y = axis_title_or(st$y_title, "Donor → acceptor pair"),
           title = "Persistent A / B hydrogen-bond pairs") +
      gg_theme_pub(st)
    ggplotly(p, tooltip = "text") |>
      layout(xaxis = plotly_axis_from_style(st, "x", axis_title_or(st$x_title, "Occupancy (%)")),
             yaxis = plotly_axis_from_style(st, "y", axis_title_or(st$y_title, "Donor → acceptor pair")),
             showlegend = FALSE) |>
      apply_plotly_theme()
  })

  output$tbl_hbond_pairs <- renderDT({
    req(rv$data$hbond_pairs)
    datatable(rv$data$hbond_pairs, options = list(pageLength = 10, scrollX = TRUE))
  })

# ───────────────────────────────────────────────────────────────────────────
  # Per-plot publication exports (each plot has its own styling + export settings)
  # ───────────────────────────────────────────────────────────────────────────

  gg_theme_pub <- function(style) {
    base_size <- style$base_size %||% 12
    switch(
      style$theme %||% "classic",
      classic = theme_classic(base_size = base_size),
      minimal = theme_minimal(base_size = base_size),
      bw      = theme_bw(base_size = base_size),
      theme_classic(base_size = base_size)
    )
  }

  apply_palette <- function(p, style) {
    pal <- style$palette %||% "default"
    if (identical(pal, "default")) return(p)
    if (pal %in% c("dark3", "set2")) {
      return(p + scale_color_brewer(palette = if (pal == "dark3") "Dark2" else "Set2") +
               scale_fill_brewer(palette = if (pal == "dark3") "Dark2" else "Set2"))
    }
    if (pal == "viridis") {
      if (requireNamespace("viridis", quietly = TRUE)) {
        return(p + viridis::scale_color_viridis(discrete = TRUE) + viridis::scale_fill_viridis(discrete = TRUE))
      } else {
        return(p + scale_color_viridis_d() + scale_fill_viridis_d())
      }
    }
    p
  }

  palette_export_cols <- function(style) {
    pal <- style$palette %||% "default"
    switch(
      pal,
      dark3 = c("#1b9e77", "#d95f02"),
      set2 = c("#66c2a5", "#fc8d62"),
      viridis = c("#440154FF", "#22A884FF"),
      c("#118ab2", "#06d6a0")
    )
  }

  apply_continuous_palette <- function(p, style, aesthetic = c("color", "fill")) {
    aesthetic <- match.arg(aesthetic)
    pal <- style$palette %||% "default"
    cols <- palette_export_cols(style)
    if (pal == "viridis") {
      if (aesthetic == "color") {
        if (requireNamespace("viridis", quietly = TRUE)) return(p + viridis::scale_color_viridis())
        return(p + scale_color_viridis_c())
      }
      if (requireNamespace("viridis", quietly = TRUE)) return(p + viridis::scale_fill_viridis())
      return(p + scale_fill_viridis_c())
    }
    if (aesthetic == "color") return(p + scale_color_gradient(low = cols[1], high = cols[2]))
    p + scale_fill_gradient(low = cols[1], high = cols[2])
  }

  save_plot_file <- function(plot_obj, file, fmt, width_in, height_in, dpi) {
    fmt <- tolower(fmt)
    if (fmt == "pdf") {
      grDevices::pdf(file, width = width_in, height = height_in, useDingbats = FALSE)
      print(plot_obj)
      grDevices::dev.off()
      return(invisible(TRUE))
    }
    if (fmt == "png") {
      grDevices::png(file, width = width_in * dpi, height = height_in * dpi, res = dpi)
      print(plot_obj)
      grDevices::dev.off()
      return(invisible(TRUE))
    }
    if (fmt == "svg") {
      if (requireNamespace("svglite", quietly = TRUE)) {
        svglite::svglite(file, width = width_in, height = height_in)
      } else {
        grDevices::svg(file, width = width_in, height = height_in)
      }
      print(plot_obj)
      grDevices::dev.off()
      return(invisible(TRUE))
    }
    stop("Unsupported export format: ", fmt)
  }

  make_timeseries_gg <- function(df, ycol, ylab, title, style, smooth_k = 0) {
    req(df)
    xcol <- if (isTRUE(input$combine_trajs) && ("time_ns_global" %in% names(df)) && any(!is.na(df$time_ns_global))) "time_ns_global" else "time_ns"
    df <- df[order(df[[xcol]]), ]

    xlab_def <- "Time (ns)"
    xlab_use <- axis_title_or(style$x_title, xlab_def)
    ylab_use <- axis_title_or(style$y_title, ylab)

    single_cols <- palette_export_cols(style)
    p <- ggplot(df, aes(x = .data[[xcol]], y = .data[[ycol]]))
    if (!isTRUE(input$combine_trajs) && ("segment" %in% names(df))) {
      p <- p + aes(color = segment)
      p <- p + geom_line(size = style$line_width %||% 1.2)
      if (isTRUE(style$show_points)) {
        p <- p + geom_point(size = style$point_size %||% 3, alpha = style$point_alpha %||% 0.8)
      }
    } else {
      p <- p + geom_line(size = style$line_width %||% 1.2, color = single_cols[1])
      if (isTRUE(style$show_points)) {
        p <- p + geom_point(size = style$point_size %||% 3, alpha = style$point_alpha %||% 0.8, color = single_cols[1])
      }
    }

    smooth_k <- as.integer(smooth_k %||% 0)
    if (smooth_k >= 2) {
      sm <- running_mean(df[[ycol]], smooth_k)
      sm_df <- data.frame(xcol_val = df[[xcol]], ycol_val = sm, stringsAsFactors = FALSE)
      sm_df <- sm_df[is.finite(sm_df$ycol_val), , drop = FALSE]
      if (nrow(sm_df) > 0) {
        p <- p + geom_line(data = sm_df,
                           aes(x = xcol_val, y = ycol_val),
                           color = "#ef476f",
                           size = (style$line_width %||% 1.2) + 0.4,
                           inherit.aes = FALSE)
      }
    }

    # scales
    if ((style$x_scale %||% "linear") == "log") p <- p + scale_x_log10()
    if ((style$y_scale %||% "linear") == "log") p <- p + scale_y_log10()

    # limits (linear only; log limits are ambiguous); single-sided limits (only min or only max) are supported
    xlim_gg <- gg_axis_lim(style$x_min, style$x_max, style$x_scale)
    ylim_gg <- gg_axis_lim(style$y_min, style$y_max, style$y_scale)
    if (!is.null(xlim_gg) || !is.null(ylim_gg)) p <- p + coord_cartesian(xlim = xlim_gg, ylim = ylim_gg)

    p <- p + labs(x = xlab_use, y = ylab_use, title = title) + gg_theme_pub(style)
    apply_palette(p, style)
  }

  make_rmsf_gg <- function(df, title = "RMSF (Selection A)", style, ycol = "rmsf_A") {
    req(df)

    xlab_use <- axis_title_or(style$x_title, "Residue")
    ylab_use <- axis_title_or(style$y_title, "RMSF (Å)")

    single_cols <- palette_export_cols(style)
    p <- ggplot(df, aes(x = resno, y = .data[[ycol]])) +
      geom_line(size = style$line_width %||% 1.2, color = single_cols[1])

    if (isTRUE(style$show_points)) {
      p <- p + geom_point(size = style$point_size %||% 3, alpha = style$point_alpha %||% 0.8, color = single_cols[1])
    }

    if ((style$x_scale %||% "linear") == "log") p <- p + scale_x_log10()
    if ((style$y_scale %||% "linear") == "log") p <- p + scale_y_log10()

    xlim_gg <- gg_axis_lim(style$x_min, style$x_max, style$x_scale)
    ylim_gg <- gg_axis_lim(style$y_min, style$y_max, style$y_scale)
    if (!is.null(xlim_gg) || !is.null(ylim_gg)) p <- p + coord_cartesian(xlim = xlim_gg, ylim = ylim_gg)

    p <- p + labs(x = xlab_use, y = ylab_use, title = title) + gg_theme_pub(style)
    apply_palette(p, style)
  }

  make_cluster_gg <- function(df, title = "Cluster distribution (MDS)", style) {
    req(df)

    xlab_use <- axis_title_or(style$x_title, "Dim1")
    ylab_use <- axis_title_or(style$y_title, "Dim2")

    p <- ggplot(df, aes(x = Dim1, y = Dim2, color = cluster)) +
      geom_point(size = style$point_size %||% 3, alpha = style$point_alpha %||% 0.8)

    if ((style$x_scale %||% "linear") == "log") p <- p + scale_x_log10()
    if ((style$y_scale %||% "linear") == "log") p <- p + scale_y_log10()

    xlim_gg <- gg_axis_lim(style$x_min, style$x_max, style$x_scale)
    ylim_gg <- gg_axis_lim(style$y_min, style$y_max, style$y_scale)
    if (!is.null(xlim_gg) || !is.null(ylim_gg)) p <- p + coord_cartesian(xlim = xlim_gg, ylim = ylim_gg)

    p <- p + labs(x = xlab_use, y = ylab_use, title = title, color = "Cluster") + gg_theme_pub(style)
    apply_palette(p, style)
  }

  make_cluster_time_gg <- function(df_frames, title = "Clusters vs time", style) {
    req(df_frames)
    df <- df_frames
    if (isTRUE(input$combine_trajs) && ("time_ns_global" %in% names(df)) && any(!is.na(df$time_ns_global))) {
      df$time_use <- df$time_ns_global
    } else {
      df$time_use <- df$time_ns
    }

    xlab_use <- axis_title_or(style$x_title, "Time (ns)")
    ylab_use <- axis_title_or(style$y_title, "Cluster")

    p <- ggplot(df, aes(x = time_use, y = as.numeric(cluster), color = factor(cluster))) +
      geom_point(size = style$point_size %||% 3, alpha = style$point_alpha %||% 0.8) +
      scale_y_continuous(breaks = sort(unique(as.numeric(df$cluster))))

    if ((style$x_scale %||% "linear") == "log") p <- p + scale_x_log10()
    if ((style$y_scale %||% "linear") == "log") p <- p + scale_y_log10()

    xlim_gg <- gg_axis_lim(style$x_min, style$x_max, style$x_scale)
    ylim_gg <- gg_axis_lim(style$y_min, style$y_max, style$y_scale)
    if (!is.null(xlim_gg) || !is.null(ylim_gg)) p <- p + coord_cartesian(xlim = xlim_gg, ylim = ylim_gg)

    p <- p + labs(x = xlab_use, y = ylab_use, title = title, color = "Cluster") + gg_theme_pub(style)
    apply_palette(p, style)
  }

  make_dimred_gg <- function(df_scores, title = "PCA distribution (Selection A)", style) {
    req(df_scores)
    df <- df_scores
    validate(need(all(c("PC1", "PC2") %in% names(df)), "PCA scores were not generated."))
    if (isTRUE(input$combine_trajs) && ("time_ns_global" %in% names(df)) && any(!is.na(df$time_ns_global))) {
      df$time_use <- df$time_ns_global
    } else {
      df$time_use <- df$time_ns
    }

    xlab_use <- axis_title_or(style$x_title, "PC1")
    ylab_use <- axis_title_or(style$y_title, "PC2")

    p <- ggplot(df, aes(x = PC1, y = PC2, color = time_use)) +
      geom_point(size = style$point_size %||% 3, alpha = style$point_alpha %||% 0.8)

    p <- apply_continuous_palette(p, style, aesthetic = "color")

    if ((style$x_scale %||% "linear") == "log") p <- p + scale_x_log10()
    if ((style$y_scale %||% "linear") == "log") p <- p + scale_y_log10()

    xlim_gg <- gg_axis_lim(style$x_min, style$x_max, style$x_scale)
    ylim_gg <- gg_axis_lim(style$y_min, style$y_max, style$y_scale)
    if (!is.null(xlim_gg) || !is.null(ylim_gg)) p <- p + coord_cartesian(xlim = xlim_gg, ylim = ylim_gg)

    p <- p + labs(x = xlab_use, y = ylab_use, color = "Time (ns)", title = title) + gg_theme_pub(style)
    p
  }

  make_pca_var_gg <- function(df_var, title = "Explained variance", style) {
    req(df_var)
    df <- normalize_pca_var_df(head(df_var, 10))
    validate(need(nrow(df) > 0, "PCA variance table is empty."))

    xlab_use <- axis_title_or(style$x_title, "Principal component")
    ylab_use <- axis_title_or(style$y_title, "Explained variance")

    single_cols <- palette_export_cols(style)
    p <- ggplot(df, aes(x = factor(PC, levels = df$PC), y = prop_var, group = 1)) +
      geom_col(fill = single_cols[1]) +
      geom_line(aes(y = cum_var, group = 1), size = style$line_width %||% 0.8, color = single_cols[2]) +
      geom_point(aes(y = cum_var), size = style$point_size %||% 3, alpha = style$point_alpha %||% 0.8, color = single_cols[2])

    if ((style$y_scale %||% "linear") == "log") p <- p + scale_y_log10()

    ylim_gg <- gg_axis_lim(style$y_min, style$y_max, style$y_scale)
    if (!is.null(ylim_gg)) p <- p + coord_cartesian(ylim = ylim_gg)

    p <- p + labs(x = xlab_use, y = ylab_use, title = title) + gg_theme_pub(style)
    p
  }

  make_cluster_pop_gg <- function(df_summary, title = "Cluster populations", style) {
    req(df_summary)

    xlab_use <- axis_title_or(style$x_title, "Cluster")
    ylab_use <- axis_title_or(style$y_title, "Number of frames")

    p <- ggplot(df_summary, aes(x = factor(cluster), y = n_frames, fill = factor(cluster))) +
      geom_col() +
      geom_text(aes(label = paste0(round(pct, 1), "%")), vjust = -0.3, size = 3.5)

    if ((style$y_scale %||% "linear") == "log") p <- p + scale_y_log10()

    ylim_gg <- gg_axis_lim(style$y_min, style$y_max, style$y_scale)
    if (!is.null(ylim_gg)) p <- p + coord_cartesian(ylim = ylim_gg)

    p <- p + labs(x = xlab_use, y = ylab_use, fill = "Cluster", title = title) + gg_theme_pub(style)
    apply_palette(p, style)
  }

  make_cluster_quality_gg <- function(df_quality, sil_overall = NA_real_, title = "Cluster quality (mean silhouette)", style) {
    req(df_quality)
    df <- df_quality
    validate(need(nrow(df) > 0, "Silhouette summary is empty."))

    xlab_use <- axis_title_or(style$x_title, "Cluster")
    ylab_use <- axis_title_or(style$y_title, "Mean silhouette")

    p <- ggplot(df, aes(x = factor(cluster), y = mean_silhouette, fill = factor(cluster),
                        text = paste0("Cluster: ", cluster, "<br>Mean silhouette: ", sprintf("%.3f", mean_silhouette)))) +
      geom_col() +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5)

    if (is.finite(sil_overall)) {
      p <- p + geom_hline(yintercept = sil_overall, linetype = "dotted", linewidth = 0.7)
    }

    ylim_gg <- gg_axis_lim(style$y_min, style$y_max, style$y_scale)
    if (is.null(ylim_gg)) ylim_gg <- c(min(-0.2, min(df$mean_silhouette, na.rm = TRUE) - 0.05),
                                        max(1,    max(df$mean_silhouette, na.rm = TRUE) + 0.05))
    p <- p + coord_cartesian(ylim = ylim_gg)

    p <- p + labs(x = xlab_use, y = ylab_use, fill = "Cluster", title = title) + gg_theme_pub(style)
    apply_palette(p, style)
  }

  # ───────────────────────────────────────────────────────────────────────────
  # Downloads — per plot
  # ───────────────────────────────────────────────────────────────────────────

  # RMSD
  output$dl_rmsdA_data <- downloadHandler(
    filename = function() "rmsd_regionA_plot_data.csv",
    content = function(file) {
      req(rv$data$rmsd_pep)
      write.csv(rv$data$rmsd_pep, file, row.names = FALSE)
    }
  )
  output$dl_rmsdB_data <- downloadHandler(
    filename = function() "rmsd_regionB_plot_data.csv",
    content = function(file) {
      req(rv$data$rmsd_lipo)
      write.csv(rv$data$rmsd_lipo, file, row.names = FALSE)
    }
  )
  output$dl_rmsdA_plot <- downloadHandler(
    filename = function() paste0("rmsd_regionA.", style_for("rmsdA")$fmt %||% "png"),
    content = function(file) {
      req(rv$data$rmsd_pep)
      st <- style_for("rmsdA")
      p <- make_rmsd_overlay_gg(rv$data$rmsd_pep, which = input$rmsdA_series %||% "both", title = "RMSD — Selection A", style = st, combined = isTRUE(input$combine_trajs), smooth_k = input$rmsdA_smooth %||% 0)
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )
  output$dl_rmsdB_plot <- downloadHandler(
    filename = function() paste0("rmsd_regionB.", style_for("rmsdB")$fmt %||% "png"),
    content = function(file) {
      req(rv$data$rmsd_lipo)
      st <- style_for("rmsdB")
      p <- make_rmsd_overlay_gg(rv$data$rmsd_lipo, which = input$rmsdB_series %||% "both", title = "RMSD — Selection B", style = st, combined = isTRUE(input$combine_trajs), smooth_k = input$rmsdB_smooth %||% 0)
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )

  # RMSF
  output$dl_rmsf_data <- downloadHandler(
    filename = function() "rmsf_regionA_plot_data.csv",
    content = function(file) {
      req(rv$data$rmsf)
      write.csv(rv$data$rmsf, file, row.names = FALSE)
    }
  )
  output$dl_rmsf_plot <- downloadHandler(
    filename = function() paste0("rmsf_regionA.", style_for("rmsf")$fmt %||% "png"),
    content = function(file) {
      req(rv$data$rmsf)
      st <- style_for("rmsf", defaults = list(show_points = TRUE))
      p <- make_rmsf_gg(rv$data$rmsf, "RMSF — Selection A", st, ycol = "rmsf_A")
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )
  output$dl_rmsfB_data <- downloadHandler(
    filename = function() "rmsf_regionB_plot_data.csv",
    content = function(file) {
      req(rv$data$rmsf_B)
      write.csv(rv$data$rmsf_B, file, row.names = FALSE)
    }
  )
  output$dl_rmsfB_plot <- downloadHandler(
    filename = function() paste0("rmsf_regionB.", style_for("rmsfB")$fmt %||% "png"),
    content = function(file) {
      req(rv$data$rmsf_B)
      st <- style_for("rmsfB", defaults = list(show_points = TRUE))
      p <- make_rmsf_gg(rv$data$rmsf_B, "RMSF — Selection B", st, ycol = "rmsf_B")
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )

  # PCA / dimensionality reduction
  output$dl_dimred_data <- downloadHandler(
    filename = function() "pca_regionA_scores.csv",
    content = function(file) {
      req(rv$data$pca_scores)
      write.csv(rv$data$pca_scores, file, row.names = FALSE)
    }
  )
  output$dl_pcaVar_data <- downloadHandler(
    filename = function() "pca_regionA_variance.csv",
    content = function(file) {
      req(rv$data$pca_var)
      write.csv(rv$data$pca_var, file, row.names = FALSE)
    }
  )
  output$dl_dimred_plot <- downloadHandler(
    filename = function() paste0("pca_regionA_distribution.", style_for("dimred")$fmt %||% "png"),
    content = function(file) {
      req(rv$data$pca_scores)
      st <- style_for("dimred", defaults = list(show_points = TRUE))
      p <- make_dimred_gg(rv$data$pca_scores, "PCA distribution — Selection A", st)
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )
  output$dl_pcaVar_plot <- downloadHandler(
    filename = function() paste0("pca_regionA_variance.", style_for("pcaVar")$fmt %||% "png"),
    content = function(file) {
      req(rv$data$pca_var)
      st <- style_for("pcaVar", defaults = list(show_points = TRUE, line_width = 0.8))
      p <- make_pca_var_gg(rv$data$pca_var, "Explained variance", st)
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )

  # Radius of gyration
  output$dl_rgA_data <- downloadHandler(
    filename = function() "rg_regionA_plot_data.csv",
    content = function(file) {
      req(rv$data$rg_pep)
      write.csv(rv$data$rg_pep, file, row.names = FALSE)
    }
  )
  output$dl_rgB_data <- downloadHandler(
    filename = function() "rg_regionB_plot_data.csv",
    content = function(file) {
      req(rv$data$rg_lipo)
      write.csv(rv$data$rg_lipo, file, row.names = FALSE)
    }
  )
  output$dl_rgA_plot <- downloadHandler(
    filename = function() paste0("rg_regionA.", style_for("rgA")$fmt %||% "png"),
    content = function(file) {
      req(rv$data$rg_pep)
      st <- style_for("rgA")
      p <- make_timeseries_gg(rv$data$rg_pep, "rg_A", "Rg (Å)", "Radius of gyration — Selection A", st, smooth_k = input$rgA_smooth %||% 0)
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )
  output$dl_rgB_plot <- downloadHandler(
    filename = function() paste0("rg_regionB.", style_for("rgB")$fmt %||% "png"),
    content = function(file) {
      req(rv$data$rg_lipo)
      st <- style_for("rgB")
      p <- make_timeseries_gg(rv$data$rg_lipo, "rg_A", "Rg (Å)", "Radius of gyration — Selection B", st)
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )

  # Membrane metrics
  output$dl_memDist_data <- downloadHandler(
    filename = function() "membrane_com_distance_plot_data.csv",
    content = function(file) {
      req(rv$data$dist)
      write.csv(rv$data$dist, file, row.names = FALSE)
    }
  )
  output$dl_memZ_data <- downloadHandler(
    filename = function() "membrane_dz_plot_data.csv",
    content = function(file) {
      req(rv$data$z)
      write.csv(rv$data$z, file, row.names = FALSE)
    }
  )
  output$dl_memDist_plot <- downloadHandler(
    filename = function() paste0("membrane_com_distance.", style_for("memDist")$fmt %||% "png"),
    content = function(file) {
      req(rv$data$dist)
      st <- style_for("memDist")
      p <- make_timeseries_gg(rv$data$dist, "dist_A", "Distance (Å)", "Membrane COM distance", st)
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )
  output$dl_memZ_plot <- downloadHandler(
    filename = function() paste0("membrane_dz.", style_for("memZ")$fmt %||% "png"),
    content = function(file) {
      req(rv$data$z)
      st <- style_for("memZ")
      p <- make_timeseries_gg(rv$data$z, "z_A", "Δz (Å)", "Membrane Δz", st)
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )

  # Advanced membrane downloads
  output$dl_thick_data <- downloadHandler(
    filename = function() "membrane_thickness.csv",
    content = function(file) {
      req(rv$data$thickness)
      write.csv(rv$data$thickness, file, row.names = FALSE)
    }
  )
  output$dl_thick_plot <- downloadHandler(
    filename = function() paste0("membrane_thickness.", style_for("memThick")$fmt %||% "png"),
    content = function(file) {
      req(rv$data$thickness)
      st <- style_for("memThick")
      p <- make_timeseries_gg(rv$data$thickness, "thickness_A", "Thickness (Å)", "Bilayer thickness", st)
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )

  output$dl_apl_data <- downloadHandler(
    filename = function() "membrane_apl.csv",
    content = function(file) {
      req(rv$data$apl)
      write.csv(rv$data$apl, file, row.names = FALSE)
    }
  )
  output$dl_apl_plot <- downloadHandler(
    filename = function() paste0("membrane_apl.", style_for("memAPL")$fmt %||% "png"),
    content = function(file) {
      req(rv$data$apl)
      st <- style_for("memAPL")
      p <- make_timeseries_gg(rv$data$apl, "APL_A2", "APL (Å² / lipid)", "Area per lipid (per leaflet)", st)
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )

  output$dl_enrich_data <- downloadHandler(
    filename = function() "lipid_enrichment.csv",
    content = function(file) {
      req(rv$data$enrich)
      write.csv(rv$data$enrich, file, row.names = FALSE)
    }
  )
  output$dl_enrich_plot <- downloadHandler(
    filename = function() paste0("lipid_enrichment.", style_for("memEnrich")$fmt %||% "png"),
    content = function(file) {
      req(rv$data$enrich)
      st <- style_for("memEnrich")
      df <- rv$data$enrich
      xcol <- if (isTRUE(input$combine_trajs) && ("time_ns_global" %in% names(df))) "time_ns_global" else "time_ns"
      xlab_def <- if (xcol == "time_ns_global") "Time (ns)" else "Time (ns)"
      ylab_def <- paste0("Headgroups within ", input$enrich_cutoff, " Å")
      xlab_use <- axis_title_or(st$x_title, xlab_def)
      ylab_use <- axis_title_or(st$y_title, ylab_def)

      p <- ggplot(df, aes(x = .data[[xcol]], y = .data[["count"]], color = .data[["lipid_type"]])) +
        geom_line(size = st$line_width) +
        labs(x = xlab_use, y = ylab_use, title = "Lipid enrichment around target") +
        gg_theme_pub(st)

      p <- apply_palette(p, st)

      if ((st$x_scale %||% "linear") == "log") p <- p + scale_x_log10()
      if ((st$y_scale %||% "linear") == "log") p <- p + scale_y_log10()

      xlim <- if (is.finite(st$x_min) && is.finite(st$x_max) && (st$x_scale %||% "linear") == "linear") c(st$x_min, st$x_max) else NULL
      ylim <- if (is.finite(st$y_min) && is.finite(st$y_max) && (st$y_scale %||% "linear") == "linear") c(st$y_min, st$y_max) else NULL
      if (!is.null(xlim) || !is.null(ylim)) p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )



  output$dl_density_data <- downloadHandler(
    filename = function() "membrane_density_profile.csv",
    content = function(file) {
      req(rv$data$mem_density)
      write.csv(rv$data$mem_density, file, row.names = FALSE)
    }
  )
  output$dl_density_plot <- downloadHandler(
    filename = function() paste0("membrane_density_profile.", style_for("memDensity")$fmt %||% "png"),
    content = function(file) {
      req(rv$data$mem_density)
      st <- style_for("memDensity")
      df <- rv$data$mem_density
      excl <- isolate(input$dens_exclude_groups)
      if (length(excl) > 0) df <- df[!(df$group %in% excl), , drop = FALSE]
      p <- ggplot(df, aes(x = z_mid, y = density, color = group)) +
        geom_line(linewidth = st$line_width) +
        labs(x = axis_title_or(st$x_title, "z (Å, membrane-centered)"),
             y = axis_title_or(st$y_title, "Average number density (atoms / Å³)"),
             title = "Membrane density profiles") +
        gg_theme_pub(st)
      p <- apply_palette(p, st)
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )

  output$dl_order_data <- downloadHandler(
    filename = function() "lipid_tail_order_profile.csv",
    content = function(file) {
      req(rv$data$mem_order)
      write.csv(rv$data$mem_order, file, row.names = FALSE)
    }
  )
  output$dl_order_plot <- downloadHandler(
    filename = function() paste0("lipid_tail_order_profile.", style_for("memOrder")$fmt %||% "png"),
    content = function(file) {
      req(rv$data$mem_order)
      st <- style_for("memOrder")
      df <- rv$data$mem_order
      df$series <- if ("series" %in% names(df)) df$series else if ("lipid_type" %in% names(df)) paste(df$lipid_type, df$chain, sep = " — ") else df$chain
      p <- ggplot(df, aes(x = segment_index, y = absS, color = series)) +
        geom_line(linewidth = st$line_width) +
        geom_point(size = st$point_size, alpha = st$point_alpha) +
        labs(x = axis_title_or(st$x_title, "Tail segment index"),
             y = axis_title_or(st$y_title, "|S|"),
             title = "Lipid tail order profile") +
        gg_theme_pub(st)
      p <- apply_palette(p, st)
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )

  output$dl_interact_ts_data <- downloadHandler(
    filename = function() "ab_interaction_timeseries.csv",
    content = function(file) {
      req(rv$data$interact_ts)
      write.csv(rv$data$interact_ts, file, row.names = FALSE)
    }
  )
  output$dl_interact_ts_plot <- downloadHandler(
    filename = function() paste0("ab_interaction_timeseries.", style_for("intTs")$fmt %||% "png"),
    content = function(file) {
      req(rv$data$interact_ts)
      st <- style_for("intTs")
      df <- rv$data$interact_ts
      metric <- input$interact_metric %||% "min_dist_A"
      xcol <- if (isTRUE(input$combine_trajs) && ("time_ns_global" %in% names(df))) "time_ns_global" else "time_ns"
      ylab <- switch(metric,
                     min_dist_A = "Minimum distance (Å)",
                     com_dist_A = "COM distance (Å)",
                     atom_contacts = paste0("Atom contacts < ", input$interact_cutoff, " Å"),
                     metric)
      p <- ggplot(df, aes(x = .data[[xcol]], y = .data[[metric]])) +
        geom_line(linewidth = st$line_width) +
        labs(x = axis_title_or(st$x_title, if (xcol == "time_ns_global") "Time (ns)" else "Time (ns)"),
             y = axis_title_or(st$y_title, ylab),
             title = "Selection A / Selection B interaction time series") +
        gg_theme_pub(st)
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )

  output$dl_interact_occ_data <- downloadHandler(
    filename = function() "ab_interaction_occupancy.csv",
    content = function(file) {
      req(rv$data$interact_occ)
      write.csv(rv$data$interact_occ, file, row.names = FALSE)
    }
  )
  output$dl_interact_occ_plot <- downloadHandler(
    filename = function() paste0("ab_interaction_occupancy.", style_for("intOcc")$fmt %||% "png"),
    content = function(file) {
      req(rv$data$interact_occ)
      st <- style_for("intOcc")
      df <- rv$data$interact_occ
      side_use <- input$interact_occ_side %||% "Selection A"
      topn <- as.integer(input$interact_occ_topn %||% 20)
      df <- df[df$side == side_use & is.finite(df$occupancy_pct), , drop = FALSE]
      validate(need(nrow(df) > 0, paste0("No occupancy rows found for ", side_use, ".")))
      df <- head(df[order(-df$occupancy_pct, df$residue), , drop = FALSE], topn)
      df <- df[order(df$occupancy_pct, df$residue), , drop = FALSE]
      df$residue_plot <- factor(df$residue, levels = df$residue)
      p <- ggplot(df, aes(x = occupancy_pct, y = residue_plot)) +
        geom_col() +
        labs(x = axis_title_or(st$x_title, "Occupancy (%)"),
             y = axis_title_or(st$y_title, "Residue"),
             title = paste0(side_use, " residue contact occupancy")) +
        gg_theme_pub(st)
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )
  output$dl_interact_pair_data <- downloadHandler(
    filename = function() "ab_interaction_pair_occupancy.csv",
    content = function(file) {
      req(rv$data$interact_pairs)
      write.csv(rv$data$interact_pairs, file, row.names = FALSE)
    }
  )
  output$dl_interact_pair_plot <- downloadHandler(
    filename = function() paste0("ab_interaction_pair_map.", style_for("intPair")$fmt %||% "png"),
    content = function(file) {
      req(rv$data$interact_pairs)
      st <- style_for("intPair")
      df <- rv$data$interact_pairs
      df <- df[is.finite(df$occupancy_pct) & df$occupancy_pct >= as.numeric(input$interact_pair_minpct %||% 5), , drop = FALSE]
      validate(need(nrow(df) > 0, "No residue-pair contacts passed the current occupancy threshold."))
      topA <- as.integer(input$interact_pair_topA %||% 20)
      topB <- as.integer(input$interact_pair_topB %||% 12)
      rankA <- aggregate(occupancy_pct ~ residue_A, data = df, FUN = sum)
      rankB <- aggregate(occupancy_pct ~ residue_B, data = df, FUN = sum)
      rankA <- rankA[order(-rankA$occupancy_pct, rankA$residue_A), , drop = FALSE]
      rankB <- rankB[order(-rankB$occupancy_pct, rankB$residue_B), , drop = FALSE]
      keepA <- head(rankA$residue_A, topA)
      keepB <- head(rankB$residue_B, topB)
      df <- df[df$residue_A %in% keepA & df$residue_B %in% keepB, , drop = FALSE]
      validate(need(nrow(df) > 0, "No residue pairs remain after applying the top-residue filters."))
      keepA <- rankA$residue_A[rankA$residue_A %in% unique(df$residue_A)]
      keepB <- rankB$residue_B[rankB$residue_B %in% unique(df$residue_B)]
      grid <- expand.grid(residue_A = keepA, residue_B = keepB, stringsAsFactors = FALSE)
      matdf <- merge(grid, df[, c("residue_A","residue_B","occupancy_pct")], by = c("residue_A","residue_B"), all.x = TRUE, sort = FALSE)
      matdf$occupancy_pct[!is.finite(matdf$occupancy_pct)] <- 0
      matdf$residue_A <- factor(matdf$residue_A, levels = rev(keepA))
      matdf$residue_B <- factor(matdf$residue_B, levels = keepB)
      p <- ggplot(matdf, aes(x = residue_B, y = residue_A, fill = occupancy_pct)) +
        geom_tile(color = "grey80", linewidth = 0.2) +
        labs(x = axis_title_or(st$x_title, "Selection B residue"),
             y = axis_title_or(st$y_title, "Selection A residue"),
             title = "A / B residue contact map") +
        gg_theme_pub(st) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid = element_blank()) +
        scale_fill_gradient(low = "#f0f4fa", high = "#11cbb7", name = "Occupancy (%)")
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )

  output$dl_hbond_ts_data <- downloadHandler(
    filename = function() "ab_hbond_proxy_timeseries.csv",
    content = function(file) {
      req(rv$data$hbond_ts)
      write.csv(rv$data$hbond_ts, file, row.names = FALSE)
    }
  )
  output$dl_hbond_ts_plot <- downloadHandler(
    filename = function() paste0("ab_hbond_proxy_timeseries.", style_for("hbTs")$fmt %||% "png"),
    content = function(file) {
      req(rv$data$hbond_ts)
      st <- style_for("hbTs")
      df <- rv$data$hbond_ts
      metric <- input$hbond_metric %||% "hbond_count"
      xcol <- if (isTRUE(input$combine_trajs) && ("time_ns_global" %in% names(df))) "time_ns_global" else "time_ns"
      ylab <- switch(metric,
                     hbond_count = paste0("H-bond proxies < ", input$hbond_cutoff, " Å"),
                     hbond_count_AtoB = "Selection A donor → B acceptor",
                     hbond_count_BtoA = "Selection B donor → A acceptor",
                     min_da_dist = "Minimum donor–acceptor distance (Å)",
                     metric)
      p <- ggplot(df, aes(x = .data[[xcol]], y = .data[[metric]])) +
        geom_line(linewidth = st$line_width) +
        labs(x = axis_title_or(st$x_title, if (xcol == "time_ns_global") "Time (ns)" else "Time (ns)"),
             y = axis_title_or(st$y_title, ylab),
             title = "A / B hydrogen-bond summary (pure-R donor/acceptor proxy)") +
        gg_theme_pub(st)
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )
  output$dl_hbond_pair_data <- downloadHandler(
    filename = function() "ab_hbond_proxy_pairs.csv",
    content = function(file) {
      req(rv$data$hbond_pairs)
      write.csv(rv$data$hbond_pairs, file, row.names = FALSE)
    }
  )
  output$dl_hbond_pair_plot <- downloadHandler(
    filename = function() paste0("ab_hbond_proxy_pairs.", style_for("hbPair")$fmt %||% "png"),
    content = function(file) {
      req(rv$data$hbond_pairs)
      st <- style_for("hbPair")
      df <- rv$data$hbond_pairs
      dir_use <- input$hbond_pair_direction %||% "all"
      if (!identical(dir_use, "all")) df <- df[df$direction == dir_use, , drop = FALSE]
      df <- df[is.finite(df$occupancy_pct) & df$occupancy_pct >= as.numeric(input$hbond_pair_minpct %||% 5), , drop = FALSE]
      validate(need(nrow(df) > 0, "No H-bond proxy pairs passed the current filters."))
      topn <- as.integer(input$hbond_pair_topn %||% 20)
      df <- head(df[order(-df$occupancy_pct, df$pair_label), , drop = FALSE], topn)
      df <- df[order(df$occupancy_pct, df$pair_label), , drop = FALSE]
      df$pair_plot <- factor(df$pair_label, levels = df$pair_label)
      p <- ggplot(df, aes(x = occupancy_pct, y = pair_plot)) +
        geom_col(width = 0.75) +
        labs(x = axis_title_or(st$x_title, "Occupancy (%)"),
             y = axis_title_or(st$y_title, "Donor → acceptor pair"),
             title = "Persistent A / B hydrogen-bond pairs") +
        gg_theme_pub(st)
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )

  # Clustering
  output$dl_cluster_data <- downloadHandler(
    filename = function() "cluster_plot_data.zip",
    content = function(file) {
      req(rv$clust)
      td <- tempdir()
      f1 <- file.path(td, "cluster_distribution_data.csv")
      f2 <- file.path(td, "cluster_time_data.csv")
      f3 <- file.path(td, "cluster_summary.csv")
      write.csv(rv$clust$plot, f1, row.names = FALSE)
      write.csv(rv$clust$frames, f2, row.names = FALSE)
      write.csv(rv$clust$summary, f3, row.names = FALSE)
      utils::zip(file, files = c(f1, f2, f3), flags = "-j")
    }
  )
  output$dl_clustDist_plot <- downloadHandler(
    filename = function() paste0("cluster_distribution.", style_for("clustDist")$fmt %||% "png"),
    content = function(file) {
      req(rv$clust)
      st <- style_for("clustDist")
      p <- make_cluster_gg(rv$clust$plot, "Cluster distribution (MDS)", st)
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )
  output$dl_clustTime_plot <- downloadHandler(
    filename = function() paste0("cluster_vs_time.", style_for("clustTime")$fmt %||% "png"),
    content = function(file) {
      req(rv$clust)
      st <- style_for("clustTime", defaults = list(point_size = 6))
      p <- make_cluster_time_gg(rv$clust$frames, "Clusters vs time", st)
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )
  output$dl_clustPop_plot <- downloadHandler(
    filename = function() paste0("cluster_population.", style_for("clustPop")$fmt %||% "png"),
    content = function(file) {
      req(rv$clust)
      st <- style_for("clustPop", defaults = list(show_points = TRUE, line_width = 0.8))
      p <- make_cluster_pop_gg(rv$clust$summary, "Cluster populations", st)
      save_plot_file(p, file, st$fmt, st$width, st$height, st$dpi)
    }
  )
  output$dl_clustDend_plot <- downloadHandler(
    filename = function() "cluster_dendrogram.pdf",
    content = function(file) {
      req(rv$clust)
      if (is.null(rv$clust$hc)) stop("Dendrogram is available only for hierarchical clustering.")
      hc <- rv$clust$hc
      k  <- max(rv$clust$frames$cluster, na.rm = TRUE)

      grDevices::pdf(file, width = 10, height = 5, useDingbats = FALSE)
      par(mar = c(3, 4, 3, 1))
      dend <- as.dendrogram(hc)
      dend <- dendrapply(dend, function(n) {
        if (is.leaf(n)) attr(n, "label") <- ""
        n
      })
      plot(dend, main = "Hierarchical clustering dendrogram",
           ylab = "Height (RMSD)", xlab = "",
           leaflab = "none", edgePar = list(lwd = 0.5))
      palette_k <- c("#06d6a0","#ef476f","#ffd166","#118ab2","#073b4c",
                     "#8338ec","#ff006e","#fb5607","#3a86ff","#70e000")
      tryCatch(
        rect.hclust(hc, k = k, border = palette_k[seq_len(min(k, length(palette_k)))]),
        error = function(e) NULL
      )
      grDevices::dev.off()
    }
  )



}

shinyApp(ui, server)
