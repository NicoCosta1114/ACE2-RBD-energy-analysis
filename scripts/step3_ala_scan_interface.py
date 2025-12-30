#!/usr/bin/env python3
"""
This is a script which does the Ala-scanning on interface residues for an A–E
complex.

It outputs a csv file with the results and a barplot of the top destabilizing
mutations by the ala mutation.
"""

import sys
import os
import csv

import matplotlib.pyplot as plt

from forcefield import VdwParamset
from energy_utils import (
    read_atoms_from_pdbqt,
    electrostatic_energy,
    vdw_energy_pair,
    get_solv_param,
    read_asa_file,
    distance,
)


# atoms that remain after mutating to ALA
ala_keep = {"N", "CA", "C", "O", "CB"}

# residue names that are glycine 
gly_names= {"GLY", "GLYC", "GLH"} 



# functions to read the interface file

def load_interface_residues(path):
    """
    It reads the interface residues file 
    """
    residues = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("Interface"):
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            chain, resname, resseq = parts
            if chain not in {"A", "E"}:
                continue
            try:
                resseq_i = int(resseq)
            except ValueError:
                continue
            residues.append((chain, resseq_i, resname))
    return residues


def filter_atoms_to_interface(atoms, iface_set):
    return [a for a in atoms if (a["chain"], a["resseq"]) in iface_set]


# functions to compute the energy components

def compute_electrostatic_energy(atoms_A, atoms_E):
    e = 0.0
    for a in atoms_A:
        for b in atoms_E:
            e += electrostatic_energy(a, b)
    return e


def compute_vdw_energy(atoms_A, atoms_E, vdw_paramset):
    e = 0.0
    for a in atoms_A:
        for b in atoms_E:
            e += vdw_energy_pair(a, b, vdw_paramset)
    return e


def solvation_energy_state(atoms, asa_dict, vdw_paramset):
    """
    G_solv = sum_i sigma_i * ASA_i
    (we rely on ASA dict already matching heavy atoms; hydrogens typically absent)
    """
    g = 0.0
    for atom in atoms:
        key = (atom["chain"], atom["resseq"], atom["atom_name"].strip())
        asa = asa_dict.get(key)
        if asa is None:
            continue
        sigma = get_solv_param(atom, vdw_paramset)
        g += sigma * asa
    return g


def compute_total_dG(atoms_A, atoms_E, vdw_pars, asa_complex, asa_A, asa_E):
    """
    This is a function to compute the total interaction energy.
    """
    e_elec = compute_electrostatic_energy(atoms_A, atoms_E)
    e_vdw  = compute_vdw_energy(atoms_A, atoms_E, vdw_pars)

    gC = solvation_energy_state(atoms_A + atoms_E, asa_complex, vdw_pars)
    gA = solvation_energy_state(atoms_A,           asa_A,       vdw_pars)
    gE = solvation_energy_state(atoms_E,           asa_E,       vdw_pars)
    dG_solv = gC - gA - gE

    total = e_elec + e_vdw + dG_solv
    return total, e_elec, e_vdw, dG_solv


# Ala mutation by atom removal

def mutate_residue_to_ala_by_filter(atoms, chain, resseq, resname):
    """
    It returns a new atom list where the residue has been mutated to Ala
    by removing sidechain atoms not in ala_keep.

    For Gly we do not add CB (not present). We handle Gly outside.
    """
    out = []
    for a in atoms:
        if a["chain"] != chain or a["resseq"] != resseq:
            out.append(a)
            continue

        an = a["atom_name"].strip()

        # this is to keep backbone + CB, and to drop other sidechain atoms
        if an in ala_keep:
            out.append(a)
        else:
            pass
    return out


# and now, we created here a function to plot a barplot representing the top
# desestabilizing mutations

def make_barplot(res_rows, out_png, top_n):

    # to pick top by ddG
    rows = [r for r in res_rows if r["ddG"] is not None]
    rows.sort(key=lambda r: r["ddG"], reverse=True)
    rows = rows[:min(top_n, len(rows))]

    labels = [f'{r["chain"]} {r["resname"]}{r["resseq"]}' for r in rows]
    values = [r["ddG"] for r in rows]

    plt.figure(figsize=(12, 5))
    plt.bar(labels, values)
    plt.axhline(0.0)
    plt.ylabel("ΔΔG = ΔG_mut - ΔG_wt (kcal/mol)")
    plt.title(f"Top {len(rows)} hotspots residues")
    plt.xticks(rotation=60, ha="right")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


# and this is the main function

def main():
    if len(sys.argv) != 7:
        print("Usage:")
        print("  python3 scripts/step3_ala_scan_interface.py \\")
        print("    data/6m0j_clean_cmip.pdbqt data/vdwprm_AD \\")
        print("    data/6m0j_clean.asa data/6m0j_A_only.asa data/6m0j_E_only.asa \\")
        print("    results/Interface_residues/interface_residues_7A.txt")
        sys.exit(1)

    pdbqt_path     = sys.argv[1]
    vdw_param_file = sys.argv[2]
    asa_complex_fn = sys.argv[3]
    asa_A_fn       = sys.argv[4]
    asa_E_fn       = sys.argv[5]
    iface_file     = sys.argv[6]

    # the output
    out_dir = "results/AlaScan"
    plot_dir = "results/plots"
    out_csv = os.path.join(out_dir, "ala_scan_results.csv")
    out_png = os.path.join(plot_dir, "ala_scan_barplot.png")

    # to load the corresponding atoms
    atoms_A_all, atoms_E_all = read_atoms_from_pdbqt(pdbqt_path)

    # to load the corresponding params + ASA
    vdw_pars   = VdwParamset(vdw_param_file)
    asa_complex = read_asa_file(asa_complex_fn)
    asa_A       = read_asa_file(asa_A_fn, chain_override="A")
    asa_E       = read_asa_file(asa_E_fn, chain_override="E")

    # interface residues + interface atoms
    iface_res_list = load_interface_residues(iface_file)
    iface_set = set((c, r) for (c, r, _rn) in iface_res_list)

    atoms_A = filter_atoms_to_interface(atoms_A_all, iface_set)
    atoms_E = filter_atoms_to_interface(atoms_E_all, iface_set)

    print(f"interface residues: {len(iface_set)}")
    print(f"interface atoms A: {len(atoms_A)}")
    print(f"interface atoms E: {len(atoms_E)}")

    # to compute the total WT interaction energy
    wt_total, wt_elec, wt_vdw, wt_solv = compute_total_dG(
        atoms_A, atoms_E, vdw_pars, asa_complex, asa_A, asa_E
    )
    print("\nWT interface ΔG:")
    print("  electrostatic:", wt_elec)
    print("  vdw         :", wt_vdw)
    print("  solvation   :", wt_solv)
    print("  TOTAL       :", wt_total)

    # this is the loop that iterates through each residue and generates
    # an Ala mutant
    results = []
    for (chain, resseq, resname) in iface_res_list:

        # here we handle Gly separately 
        if resname.strip().upper() in gly_names:
            results.append({
                "chain": chain, "resseq": resseq, "resname": resname,
                "dG_wt": wt_total,
                "dG_mut": None,
                "ddG": None,
            })
            continue

        if chain == "A":
            mut_A = mutate_residue_to_ala_by_filter(atoms_A, "A", resseq, resname)
            mut_E = atoms_E
        else:
            mut_A = atoms_A
            mut_E = mutate_residue_to_ala_by_filter(atoms_E, "E", resseq, resname)

        mut_total, mut_elec, mut_vdw, mut_solv = compute_total_dG(
            mut_A, mut_E, vdw_pars, asa_complex, asa_A, asa_E
        )

        ddG = mut_total - wt_total

        results.append({
            "chain": chain, "resseq": resseq, "resname": resname,
            "dG_wt": wt_total,
            "dG_mut": mut_total,
            "ddG": ddG,
            "wt_elec": wt_elec, "wt_vdw": wt_vdw, "wt_solv": wt_solv,
            "mut_elec": mut_elec, "mut_vdw": mut_vdw, "mut_solv": mut_solv,
        })

    # this to write the CSV
    fieldnames = [
        "chain","resseq","resname",
        "dG_wt","dG_mut","ddG",
        "wt_elec","wt_vdw","wt_solv",
        "mut_elec","mut_vdw","mut_solv",
    ]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in results:
            w.writerow(r)

    # and this is to build plot rows only with ddG
    plot_rows = []
    for r in results:
        if r.get("ddG") is None:
            continue
        plot_rows.append({
            "chain": r["chain"],
            "resseq": r["resseq"],
            "resname": r["resname"],
            "ddG": r["ddG"]
        })
    if plot_rows:
        make_barplot(plot_rows, out_png, top_n=25)
    else:
        print("No ddG values to plot")


if __name__ == "__main__":
    main()
