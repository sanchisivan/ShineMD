# ShineMD

**ShineMD** is an interactive R Shiny application for molecular dynamics (MD) trajectory analysis and visualization. It provides a fully graphical, browser-based interface for processing and exploring MD simulations — no scripting required.

Developed at the **Laboratory of Bioactive Peptides (LPB)**, Faculty of Biochemistry and Biological Sciences, National University of the Littoral (UNL), Santa Fe, Argentina.

---

## Features

| Category | Analyses |
|---|---|
| **Structural** | RMSD, RMSF, Radius of Gyration |
| **Dimensionality reduction** | PCA, Free Energy Landscape (FEL) |
| **Membrane biophysics** | Bilayer thickness, Area per lipid, Density profiles, Lipid tail order \|S\|, Lipid enrichment, COM distances |
| **Interactions** | Contact time series, Residue occupancy, Contact map (heatmap), H-bond proxy analysis |
| **Clustering** | Hierarchical (Ward.D2) and k-medoids (PAM), Structure export (medoids, centroids, single frames) |
| **Reproducibility** | Session info export, package version tracking |

All plots are interactive (pan, zoom, hover), support dark/light themes, and can be exported as PDF, PNG, or SVG. Data tables are downloadable as CSV.

---

## Input

ShineMD works with **AMBER** simulation files:

| File | Extension |
|---|---|
| Topology | `.prmtop` |
| Trajectory | `.nc` (NetCDF) |
| Optional reference structure | `.pdb`, `.rst7`, `.inpcrd`, `.crd` |

Multiple trajectory segments (e.g. `prod_1.nc`, `prod_2.nc`, …) are automatically detected, naturally sorted, and can be concatenated into a continuous timeline.

---

## Requirements

ShineMD runs in **R** (≥ 4.0 recommended). Install all required packages from CRAN:

```r
install.packages(c(
  "shiny", "shinydashboard", "shinyjs", "shinyWidgets", "shinyFiles",
  "plotly", "ggplot2", "DT", "htmlwidgets",
  "bio3d", "ncdf4", "fs"
))
```

---

## Installation and usage

Clone the repository:

```bash
git clone https://github.com/sanchisivan/ShineMD.git
```

Open R in the `ShineMD/` folder and run:

```r
shiny::runApp()
```

The app will open in your default web browser. If not, navigate to the URL shown in the R console (e.g. `http://127.0.0.1:XXXX`).

---

## Documentation

A full **User Manual** is available in this repository:

**[USER_MANUAL.md](USER_MANUAL.md)**

The manual covers every tab, parameter, analysis, and export option in detail, along with tips and troubleshooting guidance.

---

## Repository structure

```
ShineMD/
├── app.R              # Main application (UI + server)
├── README.md
├── USER_MANUAL.md     # Complete user documentation
├── CITATION.cff
├── LICENSE
├── .gitignore
└── www/               # Static assets (logos)
```

---

## Citation

If you use ShineMD in academic work, please cite this repository:

> Sanchis I. *ShineMD: an interactive Shiny application for molecular dynamics trajectory analysis.* GitHub repository. https://github.com/sanchisivan/ShineMD

A `CITATION.cff` file is included for reference managers and GitHub's *Cite this repository* feature.

---

## Contact

<p align="center">
  <img src="www/lpb_logo_unmatted.png" alt="Laboratory of Bioactive Peptides" width="220"/>
</p>

<p align="center">
  <strong>Laboratory of Bioactive Peptides (LPB)</strong><br>
  Faculty of Biochemistry and Biological Sciences (FBCB)<br>
  National University of the Littoral (UNL) · Santa Fe, Argentina
</p>

| | |
|---|---|
| **Dr. Iván Sanchis** (app development) | sanchisivan@fbcb.unl.edu.ar |
| **Prof. Álvaro Sebastián Siano** (group leader) | asiano@fbcb.unl.edu.ar |

---

## License

MIT License — see [LICENSE](LICENSE) for details.
