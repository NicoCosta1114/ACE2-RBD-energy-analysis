"""
this script computes interaction energies between chain A and chain E
(full complex and, optionally, only for interface residues)
"""

#!/usr/bin/env python3

import sys
from forcefield import VdwParamset

# the energy utilities
from energy_utils import (
    read_atoms_from_pdbqt,
    electrostatic_energy,
    vdw_energy_pair,
    get_solv_param,
    read_asa_file,
)


# electrostatic and vdw energies

def compute_electrostatic_energy(atoms_A, atoms_E):
    """this function computes total electrostatic energy between A and E"""
    e_elec = 0.0
    for a in atoms_A:
        for b in atoms_E:
            e_elec += electrostatic_energy(a, b)
    return e_elec


def compute_vdw_energy(atoms_A, atoms_E, vdw_paramset):
    """this function computes total vdw energy between A and E"""
    e_vdw = 0.0
    for a in atoms_A:
        for b in atoms_E:
            e_vdw += vdw_energy_pair(a, b, vdw_paramset)
    return e_vdw


# solvation energy 

def solvation_energy_state(atoms_A, atoms_E, asa_dict, vdw_paramset):
    """
    this function computes solvation energy for one state
    (either complex, or isolated chain A, or isolated chain E)
    """
    e_solv = 0.0
    for atom in atoms_A + atoms_E:
        key = (atom["chain"], atom["resseq"], atom["atom_name"].strip())
        asa_val = asa_dict.get(key)
        if asa_val is None:
            continue

        sigma = get_solv_param(atom, vdw_paramset)
        if sigma is None:
            continue

        e_solv += sigma * asa_val

    return e_solv


# for interface residues

def load_interface_residues(path):
    """
    this function reads a text file with interface residues
    and returns a set of (chain_id, resseq) tuples
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
    """
    this function keeps only atoms whose residue belongs to the interface set
    """
    filtered = []
    for a in atoms:
        key = (a["chain"], a["resseq"])
        if key in iface_set:
            filtered.append(a)
    return filtered



def main():
    """this is the main function"""
    if len(sys.argv) not in (6, 7):
        print("usage:")
        print("  python3 scripts/step2_compute_interaction_energy.py \\")
        print("      data/6m0j_clean_cmip.pdbqt data/vdwprm_AD \\")
        print("      data/6m0j_clean.asa 6m0j_A_only.asa data/6m0j_E_only.asa results/Interface_residues/interface_residues_7A.txt]")
        sys.exit(1)
        
    pdbqt_path    = sys.argv[1]
    vdw_param_file = sys.argv[2]
    asa_complex   = sys.argv[3]
    asa_A         = sys.argv[4]
    asa_E         = sys.argv[5]
    iface_file    = sys.argv[6] if len(sys.argv) == 7 else None

    atoms_A, atoms_E = read_atoms_from_pdbqt(pdbqt_path)
    print(f"number of atoms in chain A: {len(atoms_A)}")
    print(f"number of atoms in chain E: {len(atoms_E)}")

    vdw_pars = VdwParamset(vdw_param_file)
    print("number of atom types in vdw file:", vdw_pars.ntypes)
    
    # FULL COMPLEX
    
    e_elec = compute_electrostatic_energy(atoms_A, atoms_E)
    print("electrostatic energy A–E (kcal/mol):", e_elec)

    e_vdw = compute_vdw_energy(atoms_A, atoms_E, vdw_pars)
    print("van der waals energy A–E (kcal/mol):", e_vdw)

    asa_complex = read_asa_file(asa_complex)          # complex A+B
    asa_A = read_asa_file(asa_A, chain_override="A")  # isolated A
    asa_E = read_asa_file(asa_E, chain_override="E")  # isolated E

    g_solv_complex = solvation_energy_state(atoms_A, atoms_E,
                                            asa_complex, vdw_pars)
    g_solv_A = solvation_energy_state(atoms_A, [], asa_A, vdw_pars)
    g_solv_E = solvation_energy_state([], atoms_E, asa_E, vdw_pars)

    dG_solv_inter = g_solv_complex - g_solv_A - g_solv_E
    print("solvation term ΔG_solv^(A-B) - ΔG_solv^A - ΔG_solv^B (kcal/mol):",
          dG_solv_inter)
    
    print("\nFULL SOLVATION COMPONENTS")
    print("G_solv (complex):", g_solv_complex)
    print("G_solv (A):      ", g_solv_A)
    print("G_solv (E):      ", g_solv_E)
    print("ΔG_solv (full):  ", dG_solv_inter)


    total_dG = e_elec + e_vdw + dG_solv_inter
    print("\nTOTAL interaction energy ΔG^(A-B) (kcal/mol):", total_dG)

    # INTERFACE-ONLY (if file provided)
    
    if iface_file is not None:
        iface = load_interface_residues(iface_file)

        atoms_A_iface = filter_atoms_by_interface(atoms_A, iface)
        atoms_E_iface = filter_atoms_by_interface(atoms_E, iface)

        print(f"number of interface atoms in chain A: {len(atoms_A_iface)}")
        print(f"number of interface atoms in chain E: {len(atoms_E_iface)}")

        e_elec_iface = compute_electrostatic_energy(atoms_A_iface, atoms_E_iface)
        print("electrostatic energy interface (kcal/mol):", e_elec_iface)

        e_vdw_iface = compute_vdw_energy(atoms_A_iface, atoms_E_iface, vdw_pars)
        print("van der waals energy interface (kcal/mol):", e_vdw_iface)
        

        # solvation for complex (A+B interface atoms)
        g_solv_complex_iface = solvation_energy_state(
            atoms_A_iface, atoms_E_iface, asa_complex, vdw_pars
        )

        # solvation for A-only interface atoms
        g_solv_A_iface = solvation_energy_state(
            atoms_A_iface, [], asa_A, vdw_pars
        )

        # solvation for E-only interface atoms
        g_solv_E_iface = solvation_energy_state(
            [], atoms_E_iface, asa_E, vdw_pars
        )

        dG_solv_iface = g_solv_complex_iface - g_solv_A_iface - g_solv_E_iface
        print("solvation energy interface (kcal/mol):", dG_solv_iface)

        print("\nINTERFACE SOLVATION COMPONENTS")
        print("G_solv_iface (complex):", g_solv_complex_iface)
        print("G_solv_iface (A):      ", g_solv_A_iface)
        print("G_solv_iface (E):      ", g_solv_E_iface)
        print("ΔG_solv (iface):       ", dG_solv_iface)

        # full interface energy
        total_iface = e_elec_iface + e_vdw_iface + dG_solv_iface
        total_full  = e_elec + e_vdw + dG_solv_inter

        print("\nTOTAL interface ΔG (elec + vdw + solv):", total_iface)
        print("fraction of full ΔG explained by interface:",
              total_iface / total_full if total_full != 0.0 else "n/a")
        
if __name__ == "__main__":
    main()
