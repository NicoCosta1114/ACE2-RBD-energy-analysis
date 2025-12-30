#!/usr/bin/env python3

# this script cleans the 6m0j structure using biobb_structure_checking

import sys
from biobb_structure_checking.structure_checking import StructureChecking
import biobb_structure_checking
import biobb_structure_checking.constants as cts

# this is a simple check for input arguments
if len(sys.argv) != 3:
    print("usage: python3 scripts/step1_preprocess_structure_6m0j.py data/(input.cif|pdb) data/output.pdb")
    sys.exit(1)

input_path = sys.argv[1]
output_path = sys.argv[2]

# this is how we get the base folder of the biobb package
base_dir = biobb_structure_checking.__path__[0]

# here we are loading the default configuration
args = cts.set_defaults(base_dir, {'notebook': False})

# the files to read and to write
args['input_structure_path'] = input_path
args['output_structure_path'] = output_path
args['output_format'] = 'pdb'

# the main structure checker object
st = StructureChecking(base_dir, args)

# to make sure the structure has only one model
st.models()

# this is to keep only ACE2 (A) and RBD (E) chains
st.chains("--select A,E")

# selecting atoms with highest occupancy
st.altloc("occupancy")

# doing this we remove metal ions (HETATM)
st.metals("All")

# and doing this we remove all ligands (HETATM)
st.ligands("All")

# this removes crystallographic water molecules
st.water("yes")

# checking amide assignments
st.amide("all")

# checking chirality
st.chiral()

# this fixes missing backbone atoms or small inconsistencies
st.backbone("--fix_atoms All --fix_chain none --add_caps none")

# this is to detect and re-built missing protein side chains
# add_hydrogen already does it, but we do it to make sure we do it correctly
st.fixside("All")

# this identifies possible S-S bridges
st.getss("all")

# this removes hydrogens if they exist
st.rem_hydrogen("yes")

# this adds hydrogens and reconstructs missing side chain atoms
st.add_hydrogen("auto")

# to save the cleaned structre to output file
st._save_structure(output_path)

# and this is to see the stats of our preprocessing
st.print_stats()

