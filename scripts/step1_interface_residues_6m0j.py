#!/usr/bin/env python3

from Bio.PDB import PDBParser, NeighborSearch
import sys

# this script finds interface residues between chain A and chain E
# this is based on a distance cutoff defined by us (7 in this case)

if len(sys.argv) != 3:
    print("usage: python3 scripts/step1_interface_residues_6m0j.py structure.pdb cutoff")
    sys.exit(1)

pdb_file = sys.argv[1]
cutoff = float(sys.argv[2])

# this is the parser that reads the pdb file
parser = PDBParser(QUIET=True)
structure = parser.get_structure("prot", pdb_file)
model = structure[0]

# this selects the two interacting chains
chainA = model["A"]
chainE = model["E"]

def is_standard_residue(res):
    # this checks that the residue is not HETATM
    return res.id[0] == " "

def get_atom_list(chain):
    # this collects all atoms from standard residues
    atoms = []
    for res in chain:
        if is_standard_residue(res):
            atoms.extend(res.get_atoms())
    return atoms

def find_interface(chain1, chain2, cutoff):
    # this finds residues in chain1 that are close to atoms in chain2
    atoms_chain2 = get_atom_list(chain2)
    ns = NeighborSearch(atoms_chain2)
    interface_residues = set()

    for res in chain1:
        if not is_standard_residue(res):
            continue
        for atom in res.get_atoms():
            # this checks if any atom in chain2 is closer than the cutoff
            if ns.search(atom.coord, cutoff):
                interface_residues.add(res)
                break

    return interface_residues

# this gets the final interface lists
iface_A = find_interface(chainA, chainE, cutoff)
iface_E = find_interface(chainE, chainA, cutoff)

print("\nInterface residues in chain A")
for r in sorted(iface_A, key=lambda x: x.id[1]):
    print(f"A {r.resname} {r.id[1]}")

print("\nInterface residues in chain E")
for r in sorted(iface_E, key=lambda x: x.id[1]):
    print(f"E {r.resname} {r.id[1]}")
