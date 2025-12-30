#!/usr/bin/env python3
"""
this script takes a cleaned PDB and outputs a PDBQT with hydrogens + charges
"""

import sys
from biobb_structure_utils.utils.str_check_add_hydrogens import str_check_add_hydrogens

def main():
    if len(sys.argv) != 3:
        print("usage: python3 scripts/step1_pdb_to_pdbqt.py data/input_clean.pdb data/output_cmip.pdbqt")
        sys.exit(1)  

    inp = sys.argv[1]
    out = sys.argv[2]

    props = {
        "mode": "auto",
        "charges": True,  
    }

    str_check_add_hydrogens(
        input_structure_path=inp,
        output_structure_path=out,
        properties=props
    )


if __name__ == "__main__":
    main()
