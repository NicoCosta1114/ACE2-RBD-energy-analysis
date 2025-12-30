#!/usr/bin/env python3

"""
this script does the the optional part 1 (run it after the per_residue_interface
_energy.py script). It reads the per-residue energy table
and generates interesting plots for the report:
- total ΔG per residue (sorted barplot)
- decomposed energies (elec, vdw, solv)
- top 10 stabilizing/destabilizing residues
"""

import pandas as pd
import matplotlib.pyplot as plt

# to load the csv
csv_path = "results/Energies/per_residue_interface_energies.csv"
df = pd.read_csv(csv_path)

# to build a residue label 
df["label"] = df.apply(lambda r: f"{r.chain}:{r.resseq} {r.resname}", axis=1)

# first plot: total ΔG per residue

df_sorted = df.sort_values("E_total", ascending=False)

plt.figure(figsize=(10, 16))
plt.barh(df_sorted["label"], df_sorted["E_total"])
plt.xlabel("ΔG total (kcal/mol)")
plt.title("Per-residue total interaction energy (interface)")
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig("results/plots/total_energy_per_residue.png", dpi=300)
plt.close()

# second plot: stacked energy components 

plt.figure(figsize=(10, 16))
plt.barh(df_sorted["label"], df_sorted["E_elec"], label="electrostatic")
plt.barh(df_sorted["label"], df_sorted["E_vdw"], left=df_sorted["E_elec"], label="vdw")
plt.barh(df_sorted["label"], df_sorted["G_solv"],
         left=df_sorted["E_elec"] + df_sorted["E_vdw"],
         label="solvation")

plt.xlabel("Energy contribution (kcal/mol)")
plt.title("Energy decomposition per residue (interface)")
plt.legend()
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig("results/plots/energy_decomposition.png", dpi=300)
plt.close()

# third plot: top stabilizing and destabilizing residues 

N = 10

top_pos = df_sorted.head(N)       # most destabilizing (positive ΔG)
top_neg = df_sorted.tail(N)       # most stabilizing (negative ΔG)

plt.figure(figsize=(10, 6))
plt.bar(top_pos["label"], top_pos["E_total"], color="red")
plt.xticks(rotation=90)
plt.ylabel("ΔG (kcal/mol)")
plt.title("Top 10 destabilizing residues (interface)")
plt.tight_layout()
plt.savefig("results/plots/top10_destabilizing.png", dpi=300)
plt.close()

plt.figure(figsize=(10, 6))
plt.bar(top_neg["label"], top_neg["E_total"], color="green")
plt.xticks(rotation=90)
plt.ylabel("ΔG (kcal/mol)")
plt.title("Top 10 stabilizing residues (interface)")
plt.tight_layout()
plt.savefig("results/plots/top10_stabilizing.png", dpi=300)
plt.close()

# fourth plot: Full ΔG vs Interface ΔG

# this values are the interaction energies we obtained for each one
# (computate_interaction_energy.py script)
FULL_DG = -99.38257778363081
INTERFACE_DG = -86.65131612440344

plt.figure(figsize=(6, 6))
plt.bar(["Full complex", "Interface only"],
        [FULL_DG, INTERFACE_DG],
        color=["steelblue", "orange"])

plt.ylabel("ΔG (kcal/mol)")
plt.title("Full vs Interface Interaction Energy")

# annotate bars
plt.text(0, FULL_DG + 2, f"{FULL_DG:.1f}", ha="center")
plt.text(1, INTERFACE_DG + 2, f"{INTERFACE_DG:.1f}", ha="center")

plt.tight_layout()
plt.savefig("results/plots/full_vs_interface_energy.png", dpi=300)
plt.close()
