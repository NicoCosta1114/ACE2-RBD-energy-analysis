#!/usr/bin/env python3

"""
this script does the the optional part 1 (run it before the comparison
threholds script). It computes per-residue interaction energies for the full
complex between chain A and chain E (not only the distance-based interface).

it:
- reads a cmip-prepared pdbqt file (with charges and atom types)
- reads vdw parameters (vdwprm_AD)
- reads naccess asa files (complex, A only, E only)
- accumulates electrostatic, vdw and solvation energy per residue
- writes a csv with per-residue energies
- also writes an energy-based interface list
(residues with |E_total| >= threshold)
"""

import sys
import math
import os
from forcefield import VdwParamset
from energy_utils import (
    read_atoms_from_pdbqt,
    read_asa_file,
    electrostatic_energy,
    vdw_energy_pair,
    get_solv_param,
)


# per-residue energy accumulation

def init_residue_dict(atoms_A, atoms_E):
    """
    this function builds a dict for all residues in A and E
    keys: (chain, resseq, resname)
    values: dict with energy components
    """
    residue_energies = {}

    def add_res(atom):
        key = (atom["chain"], atom["resseq"], atom["resname"])
        if key not in residue_energies:
            residue_energies[key] = {
                "E_elec": 0.0,
                "E_vdw": 0.0,
                "G_solv_complex": 0.0,
                "G_solv_A": 0.0,
                "G_solv_E": 0.0,
            }

    for a in atoms_A:
        add_res(a)
    for e in atoms_E:
        add_res(e)

    return residue_energies

def accumulate_pair_energies(atoms_A, atoms_E, vdw_paramset, residue_energies):
    """
    this function loops over all A–E atom pairs
    and adds half of the pair energy to each residue
    """
    for a in atoms_A:
        keyA = (a["chain"], a["resseq"], a["resname"])
        for b in atoms_E:
            keyE = (b["chain"], b["resseq"], b["resname"])

            e_elec = electrostatic_energy(a, b)
            e_vdw = vdw_energy_pair(a, b, vdw_paramset) 

            # we split the pair energy 50/50 between the two residues
            residue_energies[keyA]["E_elec"] += 0.5 * e_elec
            residue_energies[keyE]["E_elec"] += 0.5 * e_elec

            residue_energies[keyA]["E_vdw"] += 0.5 * e_vdw
            residue_energies[keyE]["E_vdw"] += 0.5 * e_vdw
            

def accumulate_solvation(atoms_A, atoms_E,
                         asa_complex, asa_A, asa_E,
                         vdw_paramset, residue_energies):
    """
    this function accumulates solvation energy per residue:
    for each atom we compute sigma * ASA in each state
    """
    all_atoms = atoms_A + atoms_E

    for atom in all_atoms:
        key_res = (atom["chain"], atom["resseq"], atom["resname"])
        res_entry = residue_energies[key_res]

        sigma = get_solv_param(atom, vdw_paramset)
        if sigma is None:
            continue

        key_atom = (atom["chain"], atom["resseq"], atom["atom_name"])

        asa_c = asa_complex.get(key_atom, 0.0)
        asa_a = asa_A.get(key_atom, 0.0)
        asa_e = asa_E.get(key_atom, 0.0)

        res_entry["G_solv_complex"] += sigma * asa_c
        res_entry["G_solv_A"]       += sigma * asa_a
        res_entry["G_solv_E"]       += sigma * asa_e

def write_per_residue_csv(residue_energies, out_csv):
    """this function writes the per-residue energies to a csv file"""
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    with open(out_csv, "w") as f:
        f.write("chain,resname,resseq,E_elec,E_vdw,G_solv,E_total\n")
        for (chain, resseq, resname), vals in sorted(residue_energies.items(),
                                                     key=lambda x: (x[0][0], x[0][1])):
            g_solv = (vals["G_solv_complex"]
                      - vals["G_solv_A"]
                      - vals["G_solv_E"])
            e_total = vals["E_elec"] + vals["E_vdw"] + g_solv
            f.write(f"{chain},{resname},{resseq},"
                    f"{vals['E_elec']:.6f},{vals['E_vdw']:.6f},"
                    f"{g_solv:.6f},{e_total:.6f}\n")

def write_energy_based_interface(residue_energies, threshold, out_txt):
    """
    this function writes an energy-based interface list:
    all residues with |E_total| >= threshold
    """
    with open(out_txt, "w") as f:
        f.write(f"# energy-based interface residues |E_total| >= {threshold:.2f} kcal/mol\n")
        for (chain, resseq, resname), vals in sorted(residue_energies.items(),
                                                     key=lambda x: (x[0][0], x[0][1])):
            g_solv = (vals["G_solv_complex"]
                      - vals["G_solv_A"]
                      - vals["G_solv_E"])
            e_total = vals["E_elec"] + vals["E_vdw"] + g_solv
            if abs(e_total) >= threshold:
                f.write(f"{chain} {resname} {resseq}\n")

# the main function

def main():
    if len(sys.argv) != 8:
        print("usage:")
        print("  python3 scripts/step2_optional_compute_per_residue_full_energies.py \\")
        print("      6m0j_clean_cmip.pdbqt vdwprm_AD \\")
        print("      6m0j_clean.asa 6m0j_A_only.asa 6m0j_E_only.asa \\")
        print("      results/Energies/per_residue_full_energies.csv \\")
        print("      results/Interface_energy/interface_residues_energy_0.5kcal.txt")
        sys.exit(1)

    pdbqt_path     = sys.argv[1]
    vdw_param_file = sys.argv[2]
    asa_complex_fn = sys.argv[3]
    asa_A_fn       = sys.argv[4]
    asa_E_fn       = sys.argv[5]
    out_csv        = sys.argv[6]
    out_iface_txt  = sys.argv[7]

    energy_threshold = 0.5  # this is a threshold in kcal/mol that can be changed
    # but we chose this one because is the one which maximizes F1-score and best
    # matches the geometric interface
    
    atoms_A, atoms_E = read_atoms_from_pdbqt(pdbqt_path)
    print(f"number of atoms in chain A: {len(atoms_A)}")
    print(f"number of atoms in chain E: {len(atoms_E)}")

    vdw_pars = VdwParamset(vdw_param_file)
    print("number of atom types in vdw file:", vdw_pars.ntypes)

    asa_complex = read_asa_file(asa_complex_fn)          # complex A+B
    asa_A       = read_asa_file(asa_A_fn, chain_override="A")  # isolated A
    asa_E       = read_asa_file(asa_E_fn, chain_override="E")  # isolated E

    residue_energies = init_residue_dict(atoms_A, atoms_E)

    accumulate_pair_energies(atoms_A, atoms_E, vdw_pars, residue_energies)

    accumulate_solvation(atoms_A, atoms_E,
                         asa_complex, asa_A, asa_E,
                         vdw_pars, residue_energies)

    write_per_residue_csv(residue_energies, out_csv)

    print(" ", out_iface_txt)
    write_energy_based_interface(residue_energies, energy_threshold, out_iface_txt)


if __name__ == "__main__":
    main()
