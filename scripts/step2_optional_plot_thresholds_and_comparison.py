#!/usr/bin/env python3

"""
this script does the the optional part 1. (run it after the script
comparison_energy_thresholds.py and compute_per_residue_full_energies.py).

It plots the results of threshold_comp.csv:

- precision, recall and F1 vs |ΔG| threshold
- fraction of absolute interaction energy explained vs threshold
- a Venn diagram comparing: distance-based interface (7 Å) and
energy-based interface with the best found threshold, at 1 kcal/mol in
this case)
"""

import pandas as pd
import matplotlib.pyplot as plt
# for the Venn diagram
from matplotlib_venn import venn2  

scan_csv = "results/Interface_energy/threshold_comp.csv"
df = pd.read_csv(scan_csv)

# Plot 1: precision / recall / F1
plt.figure()
plt.plot(df["threshold_kcal"], df["precision_vs_geom"], marker="o", label="precision")
plt.plot(df["threshold_kcal"], df["recall_vs_geom"],    marker="o", label="recall")
plt.plot(df["threshold_kcal"], df["F1_vs_geom"],        marker="o", label="F1")
plt.xlabel("|ΔG| threshold (kcal/mol)")
plt.ylabel("score")
plt.title("Energy-based interface vs distance-based (7 Å)")
plt.legend()
plt.tight_layout()
plt.savefig("results/plots/threshold_precision_recall.png", dpi=300)
plt.close()

# Plot 2: fraction of energy explained
plt.figure()
plt.plot(df["threshold_kcal"], df["fraction_abs_energy_explained"], marker="o")
plt.xlabel("|ΔG| threshold (kcal/mol)")
plt.ylabel("fraction of |ΔG| explained")
plt.title("Fraction of interaction energy explained vs threshold")
plt.tight_layout()
plt.savefig("results/plots/threshold_energy_fraction.png", dpi=300)
plt.close()

# Plot 3: Venn Diagram for the best threshold = 1 kcal/mol
geom_file   = "results/Interface_energy/only_geom_thr0.50.txt"
energy_file = "results/Interface_energy/only_energy_thr0.50.txt"
overlap_file = "results/Interface_energy/overlap_thr0.50.txt"

def load_set(path):
    """Reads chain + residue number"""
    s = set()
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                chain = parts[0]
                resseq = parts[1]
                s.add(f"{chain}{resseq}")
    return s

only_geom   = load_set(geom_file)
only_energy = load_set(energy_file)
overlap     = load_set(overlap_file)

plt.figure(figsize=(6, 6))
venn2(
    subsets = (
        len(only_geom),
        len(only_energy),
        len(overlap)
    ),
    set_labels = (
        "Distance-based\n(7 Å)",
        "Energy-based\n(0.5 kcal/mol)"
    )
)

plt.title("Comparison of interface definitions")
plt.tight_layout()
plt.savefig("results/plots/venn_interface_overlap.png", dpi=300)
plt.close()

