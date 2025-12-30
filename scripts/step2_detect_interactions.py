#!/usr/bin/env python3

"""
This is a script to detect interface interactions between chain A and chain E,
classify, and plotting them.
"""

import sys
import csv
from collections import defaultdict, Counter
import pandas as pd
import matplotlib.pyplot as plt
from energy_utils import distance 


# definitions for interaction types 

donors = {"N", "NE", "NE1", "NE2", "ND1", "ND2", "NZ", "NH1", "NH2",
    "OG", "OG1", "OH"}

acceptors = {"O", "OXT", "OD1", "OD2", "OE1", "OE2",
    "OG", "OG1", "OH",
    "ND1", "NE2"}

positive_res = {"LYS", "ARG", "HIP"}
negative_res = {"ASP", "GLU"}

aromatic_ring = {"PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "HIS": {"CG", "ND1", "CD2", "CE1", "NE2"},
    "TRP": {"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},}

# the distance thresholds (this distances thresholds are taken because we saw
# that are the usual ones for each corresponding interaction type)
hbond = 3.5
saltbridge = 4.0
vdw = 4.0
aromatic = 5.5

# a broad scan cutoff to avoid useless pair checks
#scan_max = 6.0
scan_max = 7.0
#scan_max = 7.5


# to read atoms from pdbqt

def parse_atom_line_pdbqt(line):
    if not line.startswith("ATOM"):
        return None
    
    atom_name = line[12:16].strip()
    resname = line[17:20].strip()
    chain = line[21]
    resseq = int(line[22:26])
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    
    return {
        "atom": atom_name,
        "resname": resname,
        "chain": chain,
        "resseq": resseq,
        "coord": (x, y, z)
    }

def load_atoms_A_E(pdbqt_path):
    atoms_A, atoms_E = [], []
    
    with open(pdbqt_path, "r") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            at = parse_atom_line_pdbqt(line)
            if at is None:
                continue
            if at["chain"] == "A":
                atoms_A.append(at)
            elif at["chain"] == "E":
                atoms_E.append(at)
                
    return atoms_A, atoms_E


# interface residues from step 1 

def load_interface_residues(path):
    """
    It reads the interface residues file
    """
    iface = set()

    with open(path, "r") as f:
        for line in f:
            line = line.strip()

            # to skip empty lines and headers
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
                resseq = int(resseq)
            except ValueError:
                continue

            iface.add((chain, resseq))

    return iface


def filter_atoms_to_interface(atoms, iface_set):
    return [a for a in atoms if (a["chain"], a["resseq"]) in iface_set]


# interaction classification


def is_salt_bridge(a, b, d):
    if d >= saltbridge:
        return False
    
    return ((a["resname"] in positive_res and b["resname"] in negative_res) or
            (a["resname"] in negative_res and b["resname"] in positive_res))


def is_hbond(a, b, d):
    if d >= hbond:
        return False
    
    a_atom = a["atom"]
    b_atom = b["atom"]
    
    return ((a_atom in donors and b_atom in acceptors) or
            (b_atom in donors and a_atom in acceptors))


def is_aromatic(a, b, d):
    if d >= aromatic:
        return False
        
    ra = aromatic_ring.get(a["resname"])
    rb = aromatic_ring.get(b["resname"])
    
    if ra is None or rb is None:
        return False
    return (a["atom"] in ra) and (b["atom"] in rb)


def is_vdw_hydrophobic(d):
    return d < vdw


def classify_interaction(a, b, d):
    """
    This is a function which calls the other functions and and returns none if
    nocondition is met
    """
    if is_salt_bridge(a, b, d):
        return "salt bridge"
    if is_hbond(a, b, d):
        return "hydrogen bond"
    if is_aromatic(a, b, d):
        return "aromatic stacking"
    if is_vdw_hydrophobic(d):
        return "vdw hydrophobic"
    return None


# detection

def detect_interactions(atoms_A, atoms_E):
    """
    this function scans A–E atom pairs
    and returns a list of detected interactions
    """
    interactions = []
    
    for a in atoms_A:
        for b in atoms_E:
            d = distance(a["coord"], b["coord"])
            if d > scan_max:
                continue
            itype = classify_interaction(a, b, d)
            if itype is None:
                continue
            interactions.append((
                a["resname"], a["resseq"], a["atom"],
                b["resname"], b["resseq"], b["atom"],
                d, itype
            ))
    return interactions


# functions to save the results

def save_interactions_csv(interactions, out_csv):
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "A_resname", "A_resseq", "A_atom",
            "E_resname", "E_resseq", "E_atom",
            "distance_A", "interaction_type"
        ])
        
        for row in interactions:
            Ares, Aseq, Aatom, Eres, Eseq, Eatom, d, itype = row
            w.writerow([Ares, Aseq, Aatom, Eres, Eseq, Eatom, f"{d:.3f}", itype])

def summarize_by_residue(interactions):
    """
    This function returns a dict with residues as keys and a counter each detected
    atom-pair interaction.
    """
    by_res = defaultdict(Counter)
    
    for (Ares, Aseq, Aatom, Eres, Eseq, Eatom, d, itype) in interactions:
        by_res[("A", Aseq, Ares)][itype] += 1
        by_res[("E", Eseq, Eres)][itype] += 1
        
    return by_res


def save_residue_summary_csv(by_res, out_csv):
    cols = ["salt bridge", "hydrogen bond", "aromatic stacking", "vdw hydrophobic"]

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["chain", "resseq", "resname"] + cols + ["total_contacts"])
        
        for (chain, resseq, resname) in sorted(by_res.keys(), key=lambda x: (x[0], x[1], x[2])):
            c = by_res[(chain, resseq, resname)]
            total = sum(c.values())
            w.writerow([chain, resseq, resname] + [c.get(k, 0) for k in cols] + [total])


# and this is the main function

def main():
    if len(sys.argv) != 3:
        print("Usage:")
        print("  python3 scripts/step2_detect_interactions.py data/6m0j_clean_cmip.pdbqt results/Interface_residues/interface_residues_7A.txt")
        sys.exit(1)

    pdbqt_path = sys.argv[1]
    iface_path = sys.argv[2]

    atoms_A, atoms_E = load_atoms_A_E(pdbqt_path)
    print(f"atoms chain A: {len(atoms_A)}")
    print(f"atoms chain E: {len(atoms_E)}")

    iface = load_interface_residues(iface_path)
    print(f"interface residues (chain,resseq): {len(iface)}")

    atoms_A_i = filter_atoms_to_interface(atoms_A, iface)
    atoms_E_i = filter_atoms_to_interface(atoms_E, iface)
    print(f"interface atoms A: {len(atoms_A_i)}")
    print(f"interface atoms E: {len(atoms_E_i)}")

    interactions = detect_interactions(atoms_A_i, atoms_E_i)

    out_csv_pairs = "results/Interactions/interface_interactions.csv"
    out_csv_res   = "results/Interactions/interface_interactions_by_residue.csv"

    save_interactions_csv(interactions, out_csv_pairs)

    by_res = summarize_by_residue(interactions)
    save_residue_summary_csv(by_res, out_csv_res)


    # plot of Wild-type detected interactions
    df = pd.read_csv(out_csv_res)
    
    cols = ["salt bridge", "hydrogen bond", "aromatic stacking", "vdw hydrophobic"]
    totals = df[cols].sum()

    ax = totals.plot(kind="bar")
    ax.set_ylabel("Number of detected contacts")
    ax.set_title("WT interface contacts by type")
    plt.tight_layout()

    out_png = "results/plots/WT_interface_contacts_by_type.png"
    plt.savefig(out_png, dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
