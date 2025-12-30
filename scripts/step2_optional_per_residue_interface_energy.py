#!/usr/bin/env python3

"""
this script does the the optional part 1 (run it before the script
plots_interface_energies.py, and it also can be run before the
step2_optional_compute_per_residue_full_energies.py).
It computes per-residue interaction energies for interface residues between
chain A and chain E.

it uses:
- pdbqt file with cmip charges and atom types
- vdw parameter file (vdwprm_AD)
- naccess asa files for complex, A-only and E-only
- a text file with interface residues (chain resname resseq)
"""

import sys
import math
from forcefield import VdwParamset

# we reuse the common utilities from energy_utils.py
from energy_utils import (
    distance,
    parse_pdbqt_atom_line,
    read_atoms_from_pdbqt,
    mehler_solmajer_dielectric,
    electrostatic_energy,
    get_vdw_params,
    get_solv_param,
    vdw_energy_pair,
    read_asa_file,
)

# interface helpers 

def load_interface_residues(path):
    """
    reads a text file with interface residues:
    chain resname resseq
    and returns a set of (chain, resseq)
    """
    iface = set()
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            chain_id = parts[0]
            try:
                resseq = int(parts[-1])
            except ValueError:
                continue
            iface.add((chain_id, resseq))
    return iface

def filter_atoms_by_interface(atoms, iface_set):
    """keeps only atoms whose residue is in iface_set"""
    return [a for a in atoms if (a["chain"], a["resseq"]) in iface_set]

# per-residue accumulators 

def init_residue_key(atom):
    """returns a tuple identifying a residue"""
    return (atom["chain"], atom["resseq"], atom["resname"])

# the main function

def main():
    if len(sys.argv) != 7:
        print("usage:")
        print("  python3 scripts/step2_optional_per_residue_interface_energy.py \\")
        print("      data/6m0j_clean_cmip.pdbqt data/vdwprm_AD \\")
        print("      data/6m0j_clean.asa data/6m0j_A_only.asa data/6m0j_E_only.asa \\")
        print("      results/Interface_residues/interface_residues_7A.txt")
        sys.exit(1)

    pdbqt_path      = sys.argv[1]
    vdw_param_file  = sys.argv[2]
    asa_complex_file = sys.argv[3]
    asa_A_file      = sys.argv[4]
    asa_E_file      = sys.argv[5]
    iface_file      = sys.argv[6]

    atoms_A, atoms_E = read_atoms_from_pdbqt(pdbqt_path)

    vdw_pars = VdwParamset(vdw_param_file)

    asa_complex = read_asa_file(asa_complex_file)
    asa_A       = read_asa_file(asa_A_file, chain_override="A")
    asa_E       = read_asa_file(asa_E_file, chain_override="E")

    iface = load_interface_residues(iface_file)

    atoms_A_iface = filter_atoms_by_interface(atoms_A, iface)
    atoms_E_iface = filter_atoms_by_interface(atoms_E, iface)

    print(f"interface atoms A: {len(atoms_A_iface)}")
    print(f"interface atoms E: {len(atoms_E_iface)}")

    E_elec_res = {}
    E_vdw_res  = {}
    G_solv_res = {}

    # pairwise elec + vdw (interface only) 


    for a in atoms_A_iface:
        keyA = init_residue_key(a)
        for b in atoms_E_iface:
            keyB = init_residue_key(b)

            # electrostatic contribution
            e_pair = electrostatic_energy(a, b)
            E_elec_res[keyA] = E_elec_res.get(keyA, 0.0) + e_pair
            E_elec_res[keyB] = E_elec_res.get(keyB, 0.0) + e_pair

            # vdw contribution
            v_pair = vdw_energy_pair(a, b, vdw_pars)
            E_vdw_res[keyA] = E_vdw_res.get(keyA, 0.0) + v_pair
            E_vdw_res[keyB] = E_vdw_res.get(keyB, 0.0) + v_pair

    # solvation per residue (interface only) 


    # we do the same decomposition as for the total ΔG_solv:
    # G_solv_iface = G_solv_complex_iface - G_solv_A_iface - G_solv_E_iface
    # but here we distribute it per atom and then per residue

    def add_solv_for_atoms(atoms, asa_dict, sign, acc_dict):
        """sign = +1 or -1 depending on state"""
        for atom in atoms:
            key_atom = (atom["chain"], atom["resseq"], atom["atom_name"])
            asa_val  = asa_dict.get(key_atom)
            if asa_val is None:
                continue
            sigma = get_solv_param(atom, vdw_pars)
            if sigma is None:
                continue
            res_key = init_residue_key(atom)
            acc_dict[res_key] = acc_dict.get(res_key, 0.0) + sign * sigma * asa_val

    # complex
    add_solv_for_atoms(atoms_A_iface, asa_complex, +1.0, G_solv_res)
    add_solv_for_atoms(atoms_E_iface, asa_complex, +1.0, G_solv_res)

    # isolated A 
    add_solv_for_atoms(atoms_A_iface, asa_A, -1.0, G_solv_res)

    # isolated E 
    add_solv_for_atoms(atoms_E_iface, asa_E, -1.0, G_solv_res)

    # to write a csv 
    out_csv = "results/Energies/per_residue_interface_energies.csv"


    # to collect all residue keys that appear in any component
    all_residues = set(E_elec_res.keys()) | set(E_vdw_res.keys()) | set(G_solv_res.keys())

    with open(out_csv, "w") as fh:
        fh.write("chain,resseq,resname,E_elec,E_vdw,G_solv,E_total\n")
        for chain, resseq, resname in sorted(all_residues, key=lambda x: (x[0], x[1])):
            e_elec = E_elec_res.get((chain, resseq, resname), 0.0)
            e_vdw  = E_vdw_res.get((chain, resseq, resname), 0.0)
            g_solv = G_solv_res.get((chain, resseq, resname), 0.0)
            e_tot  = e_elec + e_vdw + g_solv
            fh.write(f"{chain},{resseq},{resname},{e_elec:.4f},{e_vdw:.4f},{g_solv:.4f},{e_tot:.4f}\n")
            

if __name__ == "__main__":
    main()
