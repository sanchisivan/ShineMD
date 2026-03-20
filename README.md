# ShineMD

**ShineMD** is an interactive **R Shiny** application for molecular dynamics trajectory analysis and visualization.

It provides a user-friendly interface for processing and exploring molecular dynamics simulations, with a focus on structural analysis, dimensionality reduction, clustering, interaction analysis, and publication-ready outputs.

## Main features

- Upload and process molecular dynamics projects
- Support for AMBER topology and trajectory files
- RMSD analysis
- RMSF analysis
- Radius of gyration analysis
- PCA and low-dimensional structural exploration
- Free energy landscape visualization
- Interaction analysis between atom selections
- Membrane-related analyses
- Structural clustering
- Interactive visualization with exportable plots and representative structures

## Input

ShineMD is intended to work with project folders containing simulation files such as:

- topology files (for example, AMBER `.prmtop`)
- trajectory files (for example, NetCDF `.nc`)

## Requirements

ShineMD runs in **R** and uses **Shiny** together with several scientific and visualization packages.

## Run locally

Open R in the project folder and run:

```r
shiny::runApp()
```

## Repository structure

```text
ShineMD/
├── app.R
├── README.md
├── .gitignore
├── LICENSE
├── CITATION.cff
└── www/
```

## Citation

If you use ShineMD in academic work, please cite the associated publication and/or this repository.

Suggested citation:

**Sanchis I. ShineMD: an interactive Shiny application for molecular dynamics trajectory analysis. GitHub repository.**

## Contact

**Iván Sanchis**  
Laboratory of Bioactive Peptides  
Faculty of Biochemistry and Biological Sciences  
National University of the Littoral (UNL)  
Santa Fe, Argentina

GitHub: [sanchisivan](https://github.com/sanchisivan)
