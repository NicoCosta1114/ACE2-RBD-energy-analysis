#!/usr/bin/env python3

"""
this script does the the optional part 1 (run it after the script
compute_per_residue_full_energies.py). It scans different |E_total| thresholds
to define an energy-based interface and compares it to the distance-based
interface (7 Å).
"""

import pandas as pd
import os

def load_geom_interface(path):
    """
    this function reads a text file with interface residues
    and returns a set of (chain, resseq)
    """
    iface = set()
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # to skip comments and header lines
            if line.startswith("#"):
                continue
            if line.lower().startswith("interface residues"):
                continue

            parts = line.split()
            # because we expect at least: chain, resname, resseq
            if len(parts) < 3:
                continue

            chain = parts[0]
            # last token should be the residue number
            try:
                resseq = int(parts[-1])
            except ValueError:
                # if it is not an integer, skip this line
                continue

            iface.add((chain, resseq))
    return iface


def main():
    energies_csv      = "results/Energies/per_residue_full_energies.csv"
    geom_iface_file   = "results/Interface_residues/interface_residues_7A.txt"
    out_csv           = "results/Interface_energy/threshold_comp.csv"

    # load per-residue energies
    df = pd.read_csv(energies_csv)
    # absolute total energy
    df["absE"] = df["E_total"].abs()
    total_absE = df["absE"].sum()

    # load distance-based interface
    geom_iface = load_geom_interface(geom_iface_file)
    n_geom = len(geom_iface)
    print(f"distance-based interface residues (7 Å): {n_geom}")

    # thresholds to test
    thresholds = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0]

    rows = []

    for thr in thresholds:
        sel = df[df["absE"] >= thr].copy()
        # energy-based set of residues 
        energy_set = {(row["chain"], int(row["resseq"])) for _, row in sel.iterrows()}
        n_energy = len(energy_set)

        overlap = energy_set & geom_iface
        n_overlap = len(overlap)

        precision = n_overlap / n_energy if n_energy > 0 else 0.0
        recall    = n_overlap / n_geom   if n_geom   > 0 else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

        frac_energy = sel["absE"].sum() / total_absE if total_absE > 0 else 0.0

        rows.append({
            "threshold_kcal": thr,
            "n_energy_residues": n_energy,
            "n_geom_residues": n_geom,
            "n_overlap": n_overlap,
            "precision_vs_geom": precision,
            "recall_vs_geom": recall,
            "F1_vs_geom": f1,
            "fraction_abs_energy_explained": frac_energy,
        })

    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    out_df = pd.DataFrame(rows)
    out_df.to_csv(out_csv, index=False)

    print(out_df)

    # and this is to compute overlap sets for the best found threshold
    best_row = out_df.loc[out_df["F1_vs_geom"].idxmax()]
    best_thr = float(best_row["threshold_kcal"])
    print(f"\nBest threshold according to F1: {best_thr:.2f} kcal/mol")

    sel_best = df[df["absE"] >= best_thr].copy()
    energy_best = {(row["chain"], int(row["resseq"])) for _, row in sel_best.iterrows()}

    overlap_set   = energy_best & geom_iface
    only_geom_set = geom_iface - energy_best
    only_energy_set = energy_best - geom_iface

    # to save residue lists
    def save_set(res_set, path, header):
        with open(path, "w") as f:
            f.write(f"# {header} (threshold = {best_thr:.2f} kcal/mol)\n")
            for chain, resseq in sorted(res_set, key=lambda x: (x[0], x[1])):
                f.write(f"{chain} {resseq}\n")

    overlap_file = f"results/Interface_energy/overlap_thr{best_thr:.2f}.txt"
    only_geom_file = f"results/Interface_energy/only_geom_thr{best_thr:.2f}.txt"
    only_energy_file = f"results/Interface_energy/only_energy_thr{best_thr:.2f}.txt"

    save_set(overlap_set, overlap_file, "Residues in both (geom and energy)")
    save_set(only_geom_set, only_geom_file, "Residues only in geometric interface")
    save_set(only_energy_set, only_energy_file, "Residues only in energy-based interface")


if __name__ == "__main__":
    main()
