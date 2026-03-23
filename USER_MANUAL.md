# ShineMD — User Manual

**Version 1.0.0 · March 2026**
Laboratory of Bioactive Peptides (LPB) · Faculty of Biochemistry and Biological Sciences · National University of the Littoral (UNL) · Santa Fe, Argentina

---

## Table of Contents

1. [Overview](#1-overview)
2. [Requirements and Installation](#2-requirements-and-installation)
3. [Launching the Application](#3-launching-the-application)
4. [General Interface](#4-general-interface)
5. [Project Tab — Setup and Configuration](#5-project-tab--setup-and-configuration)
6. [RMSD Tab](#6-rmsd-tab)
7. [RMSF Tab](#7-rmsf-tab)
8. [Radius of Gyration Tab](#8-radius-of-gyration-tab)
9. [PCA / DimRed Tab](#9-pca--dimred-tab)
10. [Membrane Systems Tab](#10-membrane-systems-tab)
11. [Interactions Tab](#11-interactions-tab)
12. [Clustering Tab](#12-clustering-tab)
13. [Session / Reproducibility Tab](#13-session--reproducibility-tab)
14. [About Tab](#14-about-tab)
15. [Export Capabilities](#15-export-capabilities)
16. [Output Folder Structure](#16-output-folder-structure)
17. [Tips, Notes, and Troubleshooting](#17-tips-notes-and-troubleshooting)

---

## 1. Overview

**ShineMD** is an interactive R Shiny application for molecular dynamics (MD) trajectory analysis. It provides a graphical, browser-based interface to explore, analyse, and visualise MD simulations without writing scripts — from basic structural metrics to advanced membrane properties, interaction profiling, and conformational clustering.

### Key capabilities

| Category | Analyses |
|---|---|
| **Structural** | RMSD, RMSF, Radius of Gyration |
| **Dimensionality reduction** | PCA, Free Energy Landscape (FEL) |
| **Membrane biophysics** | Bilayer thickness, Area per lipid, Density profiles, Lipid tail order \|S\|, Lipid enrichment, COM distances |
| **Interactions** | Contact time series, Residue occupancy, Contact map (heatmap), H-bond proxy analysis |
| **Clustering** | Hierarchical (Ward.D2) and k-medoids (PAM), Structure export |
| **Reproducibility** | Full session info export, package version tracking |

### Supported simulation format

ShineMD is designed around **AMBER** output files:

| File type | Extension |
|---|---|
| Topology | `.prmtop` |
| Trajectory | `.nc` (NetCDF binary) |
| Optional reference structure | `.pdb`, `.ent`, `.rst7`, `.inpcrd`, `.crd`, `.rst` |

> **Note:** Multiple trajectory segments (e.g., `prod_1.nc`, `prod_2.nc`, …) are automatically detected and can be combined into a single continuous timeline.

---

## 2. Requirements and Installation

### R and packages

ShineMD runs in **R** (≥ 4.0 recommended). The following packages are required:

```
shiny, shinydashboard, shinyjs, shinyWidgets, shinyFiles,
plotly, ggplot2, DT, htmlwidgets,
bio3d, ncdf4, fs
```

Install missing packages from CRAN:

```r
install.packages(c(
  "shiny", "shinydashboard", "shinyjs", "shinyWidgets", "shinyFiles",
  "plotly", "ggplot2", "DT", "htmlwidgets",
  "bio3d", "ncdf4", "fs"
))
```

### Download ShineMD

Clone or download the repository:

```bash
git clone https://github.com/sanchisivan/ShineMD.git
```

Or download the ZIP archive from the GitHub repository page and extract it.

---

## 3. Launching the Application

Open R (or RStudio) in the `ShineMD/` project folder and run:

```r
shiny::runApp()
```

The application will open in your default web browser. If it does not open automatically, navigate to the URL shown in the R console (typically `http://127.0.0.1:XXXX`).

---

## 4. General Interface

The application uses a **sidebar navigation** layout:

- The **left sidebar** lists all analysis tabs. Click any tab to navigate to it.
- The **main panel** shows the content of the selected tab, organised into collapsible boxes and sub-panels.
- Many panels include an **ⓘ info button** (top-right of each box) that explains what the analysis does and how to interpret results.
- Each plot has a **theme toggle** (dark/light), **export button**, and **data download** button where applicable.
- Interactive plots (powered by Plotly) support **pan, zoom, hover tooltips**, and legend toggling.

---

## 5. Project Tab — Setup and Configuration

The **Project tab** is the starting point. All analyses depend on the configuration defined here.

### 5.1 Loading a project folder

Click **Browse** under *Project folder* to navigate to the folder containing your simulation files. ShineMD will automatically scan for:

- One `.prmtop` topology file
- All `.nc` trajectory segment files (sorted naturally by numeric suffix)

> The detected files are listed so you can confirm the correct ones were found.

### 5.2 Trajectory settings

| Parameter | Description |
|---|---|
| **Time step per frame (ps)** | The time interval in picoseconds between saved frames in your `.nc` files. Used to convert frame indices to physical time (ns). |
| **Stride** | Analyse every N-th frame. Use this to reduce memory usage and computation time on long trajectories. A stride of 1 uses every frame. |
| **First frame** | Start analysis from this frame number (1-based). Leave blank to start from the beginning. |
| **Last frame** | End analysis at this frame number. Leave blank to read to the end. |
| **Combine segments** | When enabled, multiple trajectory segments are concatenated into a single continuous timeline. When disabled, each segment is kept independent. |
| **Trajectory order** | If the segments need a specific concatenation order different from alphabetical/numeric sorting, list them in the desired order in the text box provided. |

### 5.3 Atom selections

ShineMD uses two named selections that determine which atoms are used in each analysis.

#### Selection A (primary)

This is the main structure of interest — typically a protein, peptide, or nucleic acid.

- Enter AMBER-style residue numbers or ranges (e.g., `1-50` or `1,3,5-10`).
- **If left blank**, ShineMD auto-detects the macromolecular chain by looking for backbone atoms (N, CA, C, O).
- Selection A is used for RMSD, RMSF, Rg, PCA, and as one partner in interaction and clustering analyses.

#### Selection B (optional)

A secondary region — typically a ligand, membrane-interacting peptide, a second domain, or a lipid subset.

- Enter residue numbers or ranges as for Selection A.
- If not provided, all Selection B tabs will be hidden.
- Selection B is used for Rg-B, RMSD-B, interaction analyses, and some membrane analyses.

#### Exclude residues

Residue names to exclude from all selections (e.g., water, ions). Default: `HOH,WAT,Na+,Cl-`.

Adjust this list to match your simulation's solvent and ion naming conventions.

### 5.4 Alignment settings

| Option | Description |
|---|---|
| **Backbone (N, CA, C, O)** | Aligns each frame using all backbone heavy atoms of Selection A. More robust, recommended for most cases. |
| **Cα only** | Aligns using only alpha-carbon atoms. Faster; sufficient for most structural analyses. |

### 5.5 RMSD reference structure

Defines the reference frame against which RMSD is calculated:

| Mode | Description |
|---|---|
| **First frame (global)** | RMSD is calculated relative to the very first frame of the entire trajectory (or first combined frame if segments are combined). |
| **First frame of each segment** | Each trajectory segment uses its own first frame as reference. Useful when comparing independent replica runs. |
| **External file** | Upload a `.pdb`, `.rst7`, `.inpcrd`, or `.crd` file to use as a fixed reference structure. This overrides the trajectory-based reference. |

### 5.6 Membrane configuration

These settings are required only if your system contains a lipid bilayer.

| Parameter | Description |
|---|---|
| **Membrane target** | Whether the membrane analysis should track Selection B (e.g., a peptide) or Selection A relative to the bilayer. |
| **Lipid residue names** | Comma-separated AMBER residue names for all lipid types in your bilayer (e.g., `POPC,POPE,CHOL`). |
| **Center at z = 0** | Translates the system so the bilayer midplane is at z = 0 before computing density profiles and tail order parameters. Recommended for membrane systems. |
| **Water residue names** | Names used for water molecules. Default: `HOH,WAT,TIP3,TIP3P,SOL`. |
| **Ion residue names** | Names used for ions. Default: `Na+,K+,Cl-,CLA,SOD,POT,MG2,CA`. |

### 5.7 Running the basic analysis

Once all settings are configured, click **Run basic analysis**. This will:

1. Load and align all trajectory frames.
2. Compute RMSD (Selection A and B), RMSF, and Radius of Gyration.
3. Perform PCA on Selection A Cα atoms.
4. Cache results to `results_ShineMD/` in your project folder.

Progress and any warnings are shown in the **Status / Log** box at the bottom of the Project tab.

> All other analyses (membrane, interactions, clustering) are run on-demand from their respective tabs and may take additional time depending on system size and trajectory length.

---

## 6. RMSD Tab

Root Mean Square Deviation (RMSD) measures how much the structure deviates from the reference over time.

### Sub-panels

#### RMSD — Selection A

- **Time series plot:** RMSD (Å) vs simulation time (ns). Hover over points for exact values.
- **Smoothing window:** Applies a rolling mean over N frames to reduce noise. Set to 0 to disable.
- **Series type:** Choose which atom class to use:
  - *Backbone + heavy atoms* — most complete measure.
  - *Backbone only (N, CA, C, O)* — excludes flexible side chains.
  - *Heavy atoms only* — excludes backbone.
- **Statistics:** Mean, minimum, maximum, and standard deviation are shown below the plot.
- **Export:** Download the plot (PDF/PNG/SVG) and the underlying data (CSV).

#### RMSD distribution — Selection A

- Histogram showing the distribution of RMSD values over all frames.
- Adjust the number of bins (10–200) for finer or coarser resolution.
- A narrow distribution centred at a low value indicates a stable, well-equilibrated system. Multiple peaks may indicate conformational transitions.

#### RMSD — Selection B

- Identical layout to Selection A. Only visible if Selection B was defined.

---

## 7. RMSF Tab

Root Mean Square Fluctuation (RMSF) measures the time-averaged flexibility of each residue.

### RMSF — Selection A (per residue)

- **Bar plot:** RMSF (Å) per residue for all Cα atoms in Selection A.
- Each bar represents the average positional fluctuation of a residue's Cα atom across all trajectory frames.
- Toggle the **points overlay** to visualise individual data points on top of the bars.
- Adjust **line width** for cleaner export.
- High RMSF values indicate flexible regions (e.g., loop regions, termini). Low values indicate rigid or structured regions.
- **Export:** Plot (PDF/PNG/SVG) + data (CSV with residue number and RMSF value).

---

## 8. Radius of Gyration Tab

The radius of gyration (Rg) measures the compactness of a structure — the mass-weighted root mean square distance of atoms from the centre of mass.

### Sub-panels

#### Radius of gyration — Selection A

- **Time series plot:** Rg (Å) vs time (ns).
- **Smoothing window:** Same as in the RMSD tab.
- **Statistics:** Mean, min, max, std.
- An increasing Rg over time may indicate unfolding or expansion. A stable Rg suggests a compact, folded conformation is maintained.
- **Export:** Plot + data (CSV).

#### Radius of gyration — Selection B

- Same layout. Only visible if Selection B is defined.

---

## 9. PCA / DimRed Tab

Principal Component Analysis (PCA) reduces the high-dimensional trajectory to its main modes of motion, allowing you to visualise conformational sampling in a low-dimensional space.

PCA is performed on the Cα atom coordinates of Selection A after alignment.

### Sub-panels

#### PCA scatter plot

- **2D mode:** PC1 vs PC2 scatter plot; each point is one trajectory frame.
- **3D mode:** Toggle to add a PC3 axis. The plot becomes fully interactive (rotate, zoom).
- Points are coloured by trajectory time or segment, allowing you to see how sampling evolves.
- Distinct clusters of points indicate different conformational states visited during the simulation.
- **Export:** Plot + scores data (CSV with PC1, PC2, PC3 values per frame).

#### Explained variance

- Bar chart showing the percentage of total variance explained by each principal component.
- An optional cumulative variance line helps identify how many PCs capture most of the motion.
- If PC1 and PC2 together explain > 60–70% of variance, the 2D projection is a good representation of the system's dynamics.
- **Export:** Plot + variance data (CSV).

#### Free Energy Landscape (FEL)

The FEL is computed from the 2D probability distribution of PC1 vs PC2 values, converted to free energy using:

$$\Delta G(PC_1, PC_2) = -k_B T \ln P(PC_1, PC_2)$$

- **Colour map:** Low energy (blue/dark) = highly populated states (energy minima). High energy (yellow/bright) = rarely visited states.
- **Bins:** Grid resolution for the 2D histogram (20–200). More bins give finer detail but require more frames to populate adequately.
- **Temperature (K):** Used in the Boltzmann conversion. Default: 300 K.
- **Max ΔG ceiling:** Clips the colour scale at a maximum ΔG value for better contrast. Leave blank for automatic scaling.
- **Export:** Plot (PNG recommended for colour fidelity) + data (CSV with PC1, PC2, ΔG values).

---

## 10. Membrane Systems Tab

This tab provides analyses specific to membrane simulations. All membrane panels require that lipid residue names are configured in the Project tab.

Some panels open a **trajectory selection dialog** before computing, allowing you to choose which trajectory segments to include in the calculation.

### 10.1 Membrane COM distance

- Distance (Å) from the centre of mass (COM) of the selected target (Selection A or B) to the membrane COM, as a function of time.
- Useful for tracking whether a peptide or protein stays at the membrane surface, inserts, or desorbs.

### 10.2 Membrane Δz (COM)

- Z-axis displacement of the target COM relative to the membrane midplane, as a function of time.
- Positive Δz values indicate the target is above the midplane; negative values indicate it is below.
- Particularly informative for transmembrane peptides or membrane-active molecules.

### 10.3 Bilayer thickness

- Estimates the bilayer thickness per frame from the mean z-coordinates of headgroup atoms in the upper and lower leaflets.
- **Headgroup atom name:** Default `P` (phosphorus, common for phospholipids). Change to match your lipid force field.
- Thickness = mean(z_upper leaflet headgroup atoms) − mean(z_lower leaflet headgroup atoms).
- Healthy POPC bilayers typically show ~38–40 Å at physiological conditions.

### 10.4 Area per lipid (APL)

- Lateral area available per lipid molecule, calculated from the simulation box XY dimensions.
- **Lipids per leaflet:** Auto-estimated as (total lipid count / 2) by default. Override if your bilayer is asymmetric.
- APL = (box X × box Y) / lipids per leaflet, in Å².
- Typical values: ~65–70 Å² for POPC at 300 K.

### 10.5 Lipid enrichment around target

- Quantifies which lipid types preferentially associate with the target selection.
- For each lipid type, counts how often a representative headgroup atom (default: `P`) is within the cutoff distance from any atom of the target.
- **Cutoff distance (Å):** Default 8 Å. Lipids within this distance are considered "in contact".
- **Group by residue name:** Separates results by lipid type (e.g., POPC vs POPE vs CHOL).
- Enrichment > 1 for a lipid type suggests preferential interaction; < 1 suggests depletion.

### 10.6 Membrane density profiles

- 1D mass density distribution along the z-axis (perpendicular to the bilayer).
- Shows the spatial arrangement of different molecular components.
- Outputs separate density curves for:
  - Headgroup atoms (e.g., phosphorus)
  - Lipid tail atoms
  - Water molecules
  - Ions
  - Target molecule (Selection A, B, or none)

| Parameter | Description |
|---|---|
| **Headgroup atom pattern** | Regex or atom name for headgroup identification. Default: `P`. |
| **Tail atom pattern** | Regex pattern matching tail carbon atom names. Default: `^C2[0-9]+$\|^C3[0-9]+$`. |
| **Target selection** | Include density for Selection A, B, or neither. |
| **Bin width (Å)** | Histogram bin size along z-axis. Default: 1.0 Å. Smaller bins give finer resolution. |
| **Z range (half-range)** | Limits the profile to ±N Å from the bilayer centre. Leave blank to use the full z range. |

### 10.7 Lipid tail order parameter |S|

- Approximates the segmental order parameter |S| for lipid acyl chains, which reflects chain rigidity and packing order.
- Values range from 0 (disordered/fluid) to 0.5 (fully ordered/gel phase).
- Gel-phase lipids typically show |S| > 0.4; fluid-phase lipids show |S| ~ 0.1–0.25.

**Topology modes:**

- **Standard:** Both acyl chains (sn-1 and sn-2) belong to the same residue name (most AMBER lipid force fields).
- **Split:** Each chain is a different residue name (e.g., palmitic acid sn-1 = `PA`, oleoyl sn-2 = `OL`). Provide residue names and atom name regex patterns for each chain separately.

| Parameter | Description |
|---|---|
| **Chain 1 / Chain 2 residue name** | AMBER residue name for each acyl chain (split mode only). |
| **Chain 1 / Chain 2 atom regex** | Regex pattern matching the carbon atom names of each chain. |
| **Residue filter** | Restrict order calculation to specific lipid residues. Leave blank to use all. |
| **Group by lipid type** | Compute and plot separate order profiles for each lipid residue name. |

---

## 11. Interactions Tab

This tab analyses the interactions between **Selection A** and **Selection B**. Both selections must be defined. A trajectory selection dialog is shown before computing; choose the segments to include.

### 11.1 Interaction time series

- Plots one or more interaction metrics between Selection A and B over time.

| Metric | Description |
|---|---|
| **Minimum heavy-atom distance** | Shortest distance (Å) between any heavy atom pair from A and B. |
| **COM distance** | Distance between the centres of mass of A and B. |
| **Atom contacts (count)** | Number of A–B atom pairs within the contact cutoff. |

| Parameter | Description |
|---|---|
| **Atom class** | Which atoms to use: heavy atoms, backbone only, or Cα only. |
| **Contact cutoff (Å)** | Distance threshold for counting contacts. Default: 4.5 Å. |

### 11.2 Residue contact occupancy

- For each residue in Selection A (or B), counts the fraction of frames in which it makes at least one contact with any atom from the other selection.
- Displayed as a bar chart showing the **top N most contacted residues**.
- **Top N residues:** Adjust to show more or fewer residues.
- **View for:** Switch between showing occupancy for Selection A residues or Selection B residues.
- Residues with high occupancy are the key binding/contact interface residues.
- **Export:** Plot + data (CSV) + table (CSV).

### 11.3 A / B residue contact map

- A 2D heatmap where each cell represents the contact occupancy (%) between a residue from Selection A (rows) and a residue from Selection B (columns).
- Darker or more intense cells indicate more persistent contacts.

| Parameter | Description |
|---|---|
| **Top residues A / B** | Show only the N most-contacted residues from each selection. Reduces clutter. |
| **Minimum occupancy (%)** | Only show residue pairs with contact occupancy above this threshold. Default: 5%. |

- **Export:** Plot + data (CSV) + table (CSV).

### 11.4 H-bond proxy time series

- Estimates the number of hydrogen bonds between Selection A and B as a function of time.
- Uses a **distance-based proxy**: donor–acceptor pairs within the cutoff distance are counted, without a geometric angle criterion.
- Outputs four series:
  - **Total H-bond proxies:** All donor–acceptor pairs within cutoff.
  - **A donor → B acceptor:** H-bonds where A provides the hydrogen.
  - **B donor → A acceptor:** H-bonds where B provides the hydrogen.
  - **Minimum donor–acceptor distance:** Closest approach over time.

| Parameter | Description |
|---|---|
| **Donor–acceptor cutoff (Å)** | Distance threshold for D–A pairing. Default: 3.5 Å. Range: 2–5 Å. |
| **Allow N as acceptor** | Advanced option: also include nitrogen atoms as hydrogen-bond acceptors in addition to oxygen atoms. |

> Note: This is a geometric proxy, not a true H-bond calculation (no angle criterion). It is fast and suitable for identifying persistent hydrogen-bonding regions.

### 11.5 Persistent H-bond pairs

- Identifies specific donor–acceptor residue pairs that form hydrogen bonds recurrently throughout the trajectory.
- Each row in the output table represents one unique A–B residue pair and shows its contact occupancy (% of frames).

| Parameter | Description |
|---|---|
| **Direction filter** | Show all pairs, only A→B, or only B→A hydrogen bonds. |
| **Top N pairs** | Display the N most persistent pairs. Default: 20. |
| **Minimum occupancy (%)** | Exclude pairs below this threshold. Default: 5%. |

- **Export:** Plot + data (CSV) + table (CSV).

---

## 12. Clustering Tab

The Clustering tab groups trajectory frames by structural similarity (RMSD-based) and exports representative structures.

A trajectory selection dialog is shown before computing; choose the segments to include.

### 12.1 Clustering controls (left panel)

#### Algorithm

| Algorithm | Description |
|---|---|
| **Hierarchical (Ward.D2)** | Builds a dendrogram by iterative merging of frames. Deterministic. Provides a dendrogram output. |
| **k-medoids (PAM)** | Partitions frames into k clusters around actual medoid frames. More robust to outliers. |

#### Atoms used for RMSD

Defines which atoms are used to compute the pairwise RMSD matrix for clustering:

- Selection A backbone
- Selection A heavy atoms
- Selection A Cα only
- Selection B heavy atoms
- Selection A + B backbone
- Selection A + B heavy atoms

Choose the atom set that best captures the conformational differences relevant to your question.

#### Reference frame for alignment

All frames are aligned prior to clustering. The reference can be:

- **Initial structure (t = 0):** The first frame of the trajectory.
- **Frame at specific time:** A user-defined time point (quick presets: 10, 100, 250, 500 ns).

#### Clustering parameters

- **Number of clusters (k):** The number of structural groups to identify. Minimum: 2. Use the Quality panel (silhouette score) to help choose an optimal k.

Click **Run RMSD clustering** to execute. Results appear in the right panel.

### 12.2 Clustering results (right panel)

#### Distribution plot

- A 2D (or 3D) scatter plot where frames are positioned by MDS/PCA of the RMSD matrix and coloured by cluster ID.
- Compact, well-separated clusters indicate clear structural states. Overlapping clusters suggest high conformational similarity.

#### Clusters vs time

- A time series showing the cluster assignment of each frame.
- Colour-coded by cluster ID.
- Allows you to identify when conformational transitions occur and how stable each cluster is.

#### Population

- Bar chart showing the number of frames assigned to each cluster.
- Large clusters represent highly populated (thermodynamically favoured) conformational states.

#### Quality (silhouette coefficient)

- The silhouette coefficient measures how similar a frame is to its own cluster compared to adjacent clusters. Range: −1 to +1.
- Values close to +1 indicate well-separated clusters. Values near 0 suggest overlapping clusters. Negative values indicate possible misassignment.
- Use this panel to compare different values of k and choose the clustering that gives the best quality.

#### Dendrogram (hierarchical clustering only)

- A tree diagram showing the hierarchical merging of frames.
- The height of each node represents the RMSD distance at which the merge occurred.
- Exportable as a PDF for publication use.

#### Pairwise RMSD heatmap

- A symmetric matrix where cell (i, j) shows the RMSD between frame i and frame j.
- Reveals the overall structural diversity of the trajectory and can identify conformational transitions.
- Adjust the **max RMSD scale** to improve contrast.

### 12.3 Structure export

After clustering, representative structures can be exported:

#### Medoids and centroids

- **Medoid:** The actual trajectory frame that is closest to the cluster centre — a physically meaningful structure.
- **Centroid:** The mean coordinates of all frames in a cluster — not a real conformation but useful as a geometric reference.

**Export atom selection:**

| Option | Contents |
|---|---|
| System (no solvent/ions) — heavy atoms | Protein + ligand heavy atoms only |
| System (no solvent/ions) — all atoms | Protein + ligand including hydrogens |
| System (no solvent/ions) — backbone only | N, CA, C, O of protein only |
| Full system — heavy atoms | All heavy atoms including solvent and ions |
| Full system — all atoms | Complete system |

Click **Export medoids + centroids** to write PDB files to `results_ShineMD/clustering_rmsd/`.

#### Single frame export

Export any individual frame as a PDB file:

| Mode | Description |
|---|---|
| **Global time** | Pick a frame by simulation time (ns). |
| **Frame index** | Pick a frame by its index number. |
| **Cluster representative** | Export the medoid or centroid-nearest frame for a chosen cluster. |

---

## 13. Session / Reproducibility Tab

This tab records the computational environment used to generate your results — important for reproducibility in scientific publications.

### Session info

- Displays the full output of `sessionInfo()` in R: R version, platform, operating system, loaded packages and their versions.
- Click **Download session info** to save a `.txt` file.
- **Recommended:** Include this file as supplementary material when publishing results generated with ShineMD.

### Installed package versions

- An interactive table listing all loaded packages, their versions, and source.
- Useful for checking dependency versions or troubleshooting.

### ShineMD environment

- A summary of the current analysis configuration: project folder path, trajectory files loaded, selections used, time settings, and analysis parameters.

---

## 14. About Tab

The About tab shows authorship, institutional affiliation, contact details, and citation information.

**Contact:**

| Person | Role | Email |
|---|---|---|
| Prof. Álvaro Sebastián Siano | Group leader | asiano@fbcb.unl.edu.ar |
| Dr. Iván Sanchis | App development | sanchisivan@fbcb.unl.edu.ar |

**Suggested citation:**

> Sanchis I. *ShineMD: an interactive Shiny application for molecular dynamics trajectory analysis.* GitHub repository. https://github.com/sanchisivan/ShineMD

---

## 15. Export Capabilities

Every analysis panel in ShineMD provides download buttons. The following formats are available depending on the analysis:

| Output type | Formats |
|---|---|
| **Plots** | PDF, PNG, SVG |
| **Data tables** | CSV |
| **Bulk downloads** | ZIP archive (for large result sets) |
| **Structures** | PDB (medoids, centroids, single frames) |
| **Session info** | Plain text (`.txt`) |

Results are saved to `<project_folder>/results_ShineMD/` when using the Export buttons in the Clustering tab (structure files).

---

## 16. Output Folder Structure

After running analyses and exporting results, ShineMD creates the following folder structure inside your project directory:

```
project_folder/
└── results_ShineMD/
    ├── rmsd_regionA.csv / .pdf
    ├── rmsd_regionB.csv / .pdf          (if Selection B defined)
    ├── rmsf_regionA.csv / .pdf
    ├── rg_regionA.csv / .pdf
    ├── rg_regionB.csv / .pdf            (if Selection B defined)
    ├── pca_regionA_scores.csv
    ├── pca_regionA_variance.csv
    ├── pca_regionA_distribution.pdf
    ├── pca_regionA_variance.pdf
    ├── free_energy_landscape.png / .csv
    ├── membrane_thickness.csv / .pdf
    ├── membrane_apl.csv / .pdf
    ├── membrane_enrichment.csv / .pdf
    ├── membrane_density.csv / .pdf
    ├── membrane_order.csv / .pdf
    ├── ab_interaction_timeseries.csv / .pdf
    ├── ab_occupancy.csv / .pdf
    ├── ab_contactmap.csv / .pdf
    ├── hbond_timeseries.csv / .pdf
    ├── hbond_pairs.csv / .pdf
    └── clustering_rmsd/
        ├── medoid_cluster1.pdb
        ├── medoid_cluster2.pdb
        ├── centroid_cluster1.pdb
        ├── centroid_cluster2.pdb
        ├── cluster_summary.csv
        ├── cluster_dendrogram.pdf
        └── pairwise_rmsd_heatmap.png
```

---

## 17. Tips, Notes, and Troubleshooting

### General recommendations

- **Always check the Project tab log** after running the basic analysis to confirm that all files were loaded correctly and no errors occurred.
- **Use stride wisely:** A stride of 1 uses every frame, which is most accurate but slowest. For initial exploration, a stride of 5–10 is often sufficient.
- **Combine segments** when you want a continuous time axis in all plots. Leave uncombined if you want to compare independent replicas.
- **Export session info** before closing the app if you plan to publish or archive results.

### Atom selection syntax

- Residue numbers are 1-based and match the AMBER topology numbering.
- Use commas and hyphens: `1-50` for a range, `1,5,10` for individual residues, `1-50,75-100` for multiple ranges.
- The **Exclude residues** field uses residue names (not numbers): `HOH,WAT,Na+,Cl-`.

### Membrane analyses

- Make sure **lipid residue names** in the Project tab exactly match the names in your topology. A mismatch will result in empty or incorrect profiles.
- Enable **Center at z = 0** when comparing density profiles or order parameters across different trajectories or conditions.
- For non-standard lipids, adjust the headgroup atom pattern and tail atom regex to match your force field atom naming.

### Clustering

- Start with a small k (e.g., 3–5) and inspect the silhouette score. Increase k if clusters overlap significantly.
- The **pairwise RMSD heatmap** is the most computationally demanding output. For very long trajectories (> 10,000 frames without stride), it may take several minutes.
- Medoid PDB files are physically meaningful and suitable for further analysis (docking, visual inspection, etc.). Centroid files are not real conformations.

### H-bond proxy analysis

- The proxy method is intentionally simple for performance. For rigorous hydrogen bond analysis with angle criteria, use dedicated tools (e.g., CPPTRAJ, VMD) on the exported medoid structures.
- Lowering the donor–acceptor cutoff to 3.0 Å gives stricter (stronger) hydrogen bond candidates. Values above 3.5 Å start to include weak or marginal interactions.

### Free Energy Landscape

- A reliable FEL requires good conformational sampling. If your simulation has not adequately sampled the conformational space, the FEL will reflect only partial coverage and should be interpreted cautiously.
- Increase the number of bins for smoother maps on long, well-sampled trajectories. Use fewer bins when trajectory length is limited.
- Export the FEL as PNG rather than SVG for best colour gradient reproduction.

### Performance

- ShineMD runs entirely in memory. Very large trajectories (many frames × large system) may require significant RAM.
- Running on a machine with 16+ GB RAM is recommended for systems with > 100,000 atoms and > 5,000 frames without stride.
- Results are cached in `results_ShineMD/`. If you re-run the basic analysis after changing settings, the cache is overwritten.

---

*ShineMD User Manual · Version 1.0.0 · March 2026*
*Laboratory of Bioactive Peptides · Faculty of Biochemistry and Biological Sciences · UNL · Santa Fe, Argentina*
