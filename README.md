# ACE2-RBD-energy-analysis
Energy decomposition and interface analysis of the SARS-CoV-2 ACE2–RBD complex.

This repository contains the code and results for the energy analysis of the
SARS-CoV-2 Spike RBD – ACE2 complex (PDB: 6M0J).

## Project structure
- `scripts/`: Python scripts implementing each analysis step
- `data/`: Input structures
- `results/`: Generated plots, tables and PyMOL images
- `report/`: Final PDF report

## Workflow
1. Structure preprocessing and interface definition
2. Interaction energy calculation
3. Optional energy-based interface analysis
4. Alanine scanning
5. Variant analysis and FoldX comparison

## Requirements
- Python 3
- biobb_structure_checking
- BioPython
- NACCESS
- FoldX (for comparison)

## Authors
- Luis Carlos Ospina
- Joseph Nosa
- Nicolás Costa
