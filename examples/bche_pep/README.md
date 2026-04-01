## bche_pep example

This folder contains a reduced BChE-peptide example dataset for testing the structural and interaction-oriented features of ShineMD.

Files:

- `bche_pep.prmtop`
- `bche_pep.nc`

Provenance:

- Derived from the same BChE-peptide system used as an example in the ShineMD manuscript.
- The peptide was truncated before inclusion here because the full system is still based on unpublished data.
- This dataset is intended for app demonstration and feature validation, especially RMSD, RMSF, radius of gyration, and the Interactions tab.

Recommended Project tab settings:

- `dt per saved frame (ps)`: use the value appropriate for this reduced trajectory
- `Stride`: `1`
- `First frame`: leave blank
- `Last frame`: leave blank
- `Combine segments`: enabled
- `Trajectory order`: `bche_pep.nc`
- `Selection A residue numbers`: `1:535`
- `Selection B residue numbers`: `536:539`
- `Exclude residue names`: `HOH,WAT,Na+,Cl-`
- `Alignment selection`: `Selection A backbone (N,CA,C,O)`
- `RMSD reference mode`: `First frame of selected subset (global)`
- `External coordinate reference`: none

Alternative Selection A for RMSD-focused analyses:

- `Selection A residue numbers`: `5:530`

This trimmed enzyme range can be useful if you want to avoid highly flexible termini when interpreting RMSD.

What this example is good for:

- `RMSD` for Selection A and Selection B
- `RMSF` for Selection A
- `Radius of gyration` for Selection A and Selection B
- all plots in the `Interactions` tab using the `Compute` buttons

Notes:

- Because both Selection A and Selection B are defined, this is the recommended example for testing A/B interaction analyses.
- If you want a more stable enzyme-reference RMSD, prefer `5:530`; if you want the full enzyme selection, use `1:535`.
