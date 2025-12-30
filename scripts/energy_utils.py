#!/usr/bin/env python3
"""
this module groups common utility functions for interaction-energy calculations:
- distance between 3D points
- parsing of cmip-prepared pdbqt files
- mehler-solmajer dielectric
- electrostatic energy
- vdw and solvation parameter lookup (cmip atom types)
- lennard-jones vdw pair energy
- reading naccess asa files
"""

import math


# basic geometry

def distance(coord1, coord2):
    """this function computes euclidean distance between two 3d points"""
    x1, y1, z1 = coord1
    x2, y2, z2 = coord2
    dx = x1 - x2
    dy = y1 - y2
    dz = z1 - z2
    return math.sqrt(dx*dx + dy*dy + dz*dz)


# pdbqt parsing (cmip-prepared)

def parse_pdbqt_atom_line(line):
    """this function parses one ATOM line from a pdbqt file"""
    atom_name = line[12:16].strip()
    resname   = line[17:20].strip()
    chain_id  = line[21].strip()
    resseq    = int(line[22:26].strip())

    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])

    # charge and atom type at the end of the line
    fields    = line.split()
    atom_type = fields[-1]
    charge    = float(fields[-2])

    return {
        "chain":      chain_id,
        "resname":    resname,
        "resseq":     resseq,
        "atom_name":  atom_name,
        "coord":      (x, y, z),
        "charge":     charge,
        "atom_type":  atom_type,
    }


def read_atoms_from_pdbqt(pdbqt_path):
    """
    this function reads all atoms from a pdbqt file
    and splits them into two lists: chain A and chain E
    """
    atoms_A = []
    atoms_E = []

    with open(pdbqt_path, "r") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom = parse_pdbqt_atom_line(line)
            if atom["chain"] == "A":
                atoms_A.append(atom)
            elif atom["chain"] == "E":
                atoms_E.append(atom)

    return atoms_A, atoms_E


# electrostatic energy (Mehler–Solmajer + Coulomb)

def mehler_solmajer_dielectric(r):
    if r < 0.1:
        r = 0.1
    return 86.9525 / (1.0 - 7.7839 * math.exp(-0.3153 * r)) - 8.5525


def electrostatic_energy(atomA, atomB):
    """
    this computes electrostatic energy (kcal/mol)
    between two atoms from chain A and chain E
    using their charges and the MS dielectric
    """
    q1 = atomA["charge"]
    q2 = atomB["charge"]

    r   = distance(atomA["coord"], atomB["coord"])
    eps = mehler_solmajer_dielectric(r)

    k = 332.16  # coulomb constant in kcal/mol
    return k * q1 * q2 / (eps * r)



def normalize_type(atom_type, atom_name):
    return atom_type


def _require_type_in_paramset(at_type, vdw_paramset, context="vdw/solv"):
    """this helper raises if a normalized type is not in the parameter set"""
    if at_type not in vdw_paramset.at_types:
        raise ValueError(
            f"Missing {context} parameters for atom type '{at_type}'. "
            "Fix the mapping or extend the vdw parameter file."
        )


# vdw parameters and energy 

def get_vdw_params(atom, vdw_paramset):
    """
    this function gets eps and sig for an atom
    based on normalized AD4 atom type
    """
    at = normalize_type(atom["atom_type"], atom["atom_name"])
    _require_type_in_paramset(at, vdw_paramset, context="vdw")
    p = vdw_paramset.at_types[at]
    return p.eps, p.sig


def vdw_energy_pair(a1, a2, vdw_paramset):
    """this function computes lj vdw energy for one atom pair"""
    r = distance(a1["coord"], a2["coord"])
    if r < 1e-6:
        return 0.0

    # this calls the previous function to get the corresponding parameters
    eps_i, sig_i = get_vdw_params(a1, vdw_paramset)
    eps_j, sig_j = get_vdw_params(a2, vdw_paramset)

    # this are lorentz–berthelot combination rules in order to obtain the values
    # for a pair of atoms, and not for only one atom
    eps_ij = (eps_i * eps_j) ** 0.5
    sig_ij = 0.5 * (sig_i + sig_j)

    sr  = sig_ij / r
    sr6 = sr ** 6
    sr12 = sr6 * sr6

    return 4.0 * eps_ij * (sr12 - sr6)



# solvation parameters (fsrf) and NACCESS ASA

def get_solv_param(atom, vdw_paramset):
    """this function gets fsrf (the solvation parameter) for an atom"""
    at = normalize_type(atom["atom_type"], atom["atom_name"])
    _require_type_in_paramset(at, vdw_paramset, context="solvation")
    p = vdw_paramset.at_types[at]
    return p.fsrf


def read_asa_file(asa_path, chain_override=None):
    """
    this function reads a naccess .asa file
    and returns a dict
    """
    asa = {}
    with open(asa_path, "r") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue

            atom_name = line[12:16].strip()
            chain     = line[21].strip()
            resseq    = int(line[22:26])

            fields    = line.split()
            asa_value = float(fields[-2])  # second last number is ASA

            if chain_override is not None:
                chain = chain_override

            key = (chain.strip(), resseq, atom_name.strip())
            asa[key] = asa_value

    return asa
