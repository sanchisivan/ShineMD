## membrane_pep example

This folder contains a reduced membrane-peptide example dataset for testing the membrane-oriented features of ShineMD.

Files:

- `mem_pep.prmtop`
- `mem_pep.nc`

Provenance:

- Derived from the same membrane system used as an example in the ShineMD manuscript.
- The peptide was truncated to its last 4 residues before inclusion here because the full system is still based on unpublished data.
- This dataset is intended for app demonstration and feature validation, especially membrane density and lipid tail-order analyses.

Recommended Project tab settings:

- `dt per saved frame (ps)`: `2000`
- `Stride`: `1`
- `First frame`: leave blank
- `Last frame`: leave blank
- `Combine segments`: enabled
- `Trajectory order`: `mem_pep.nc`
- `Selection A residue numbers`: `266:269`
- `Selection B residue numbers`: leave blank
- `Exclude residue names`: `HOH,WAT,Na+,Cl-`
- `Alignment selection`: `Selection A backbone (N,CA,C,O)`
- `RMSD reference mode`: `First frame of selected subset (global)`
- `External coordinate reference`: none
- `Membrane metrics target`: `Selection B (if provided)`
- `Membrane residue names`: `PA,PE,OL,PGR`
- `Center membrane at z = 0`: enabled
- `Water residue names`: `HOH,WAT,TIP3,TIP3P,SOL`
- `Ion residue names`: `Na+,K+,Cl-,CLA,SOD,POT,MG2,CA`

Recommended Membrane tab settings:

For density profiles:

- `Headgroup atom name(s)`: `P31`
- `Tail atom regex`: `^C2[0-9]+$|^C3[0-9]+$`
- `Include target density`: `Selection A`
- `Bin width (Å)`: `1`
- `Half-profile range |z| (Å)`: leave blank

For tail order:

- `Tail-order topology mode`: `Split residues: chain 1 and chain 2 are different residue names`
- `Residue names for tail order`: `PA,OL`
- `sn1 / chain 1 residue names`: `PA`
- `sn2 / chain 2 residue names`: `OL`
- `sn1 chain atom regex`: `^C1([2-9]|[1-9][0-9])$`
- `sn2 chain atom regex`: `^C1([2-9]|[1-9][0-9])$`
- `Separate by lipid type`: enabled

Notes:

- This example is meant for membrane-focused analyses, not for reproducing the full unpublished simulation system.
- Because Selection B is left blank in this reduced example, target-dependent membrane plots should use `Selection A` when you want the peptide density overlaid.
