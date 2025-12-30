#!/usr/bin/env python3

"""
this is a script that takes a pdb with two chains (A and E)
and writes two new pdb files:
- one with only chain A
- one with only chain E
"""

import sys
from Bio.PDB import PDBParser, PDBIO

def save_chain(input_pdb, output_pdb, chain_id):
    """
    this function saves only one chain from the input pdb into a new file
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", input_pdb)
    model = structure[0]

    # this is a new structure object that will contain only the selected chain
    from Bio.PDB import Structure, Model, Chain

    new_structure = Structure.Structure("single_chain")
    new_model = Model.Model(0)
    new_structure.add(new_model)

    # this selects the chain we want
    chain = model[chain_id]

    # this creates a new chain object and copies residues
    new_chain = Chain.Chain(chain_id)
    for res in chain:
        if res.id[0] == " ":  # to skip hetero residues, just in case
            new_chain.add(res.copy())

    new_model.add(new_chain)

    # it writes the new structure to a pdb file
    io = PDBIO()
    io.set_structure(new_structure)
    io.save(output_pdb)

def main():
    """
    the main function that: reads input pdb and writes two pdb files: 
    one for chain A and one for chain E
    """
    if len(sys.argv) != 4:
        print("usage: python3 scripts/step1_split_chains.py data/input_clean.pdb data/output_chainA.pdb data/output_chainE.pdb")
        sys.exit(1)

    input_pdb = sys.argv[1]
    out_A = sys.argv[2]
    out_E = sys.argv[3]

    save_chain(input_pdb, out_A, "A")

    save_chain(input_pdb, out_E, "E")

if __name__ == "__main__":
    main()
