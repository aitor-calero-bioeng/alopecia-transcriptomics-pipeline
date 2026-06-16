

import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt

# 1. Load raw files
print("Loading files...")
matrix = pd.read_csv("GSE90594_series_matrix.txt.gz", sep="\t", comment="!", index_col=0, low_memory=False)
annot = pd.read_csv("GPL17077-17467.txt", sep="\t", comment="#", index_col=0, low_memory=False)

# 2. Merge and clean
print("Processing data matrices...")
merged_data = pd.merge(annot[['GENE_SYMBOL']], matrix, left_index=True, right_index=True)
clean_data = merged_data.dropna(subset=['GENE_SYMBOL']).copy()

alopecia_cols = [f"GSM{i}" for i in range(2407372, 2407386)]
healthy_cols = [f"GSM{i}" for i in range(2407386, 2407400)]
all_patients = alopecia_cols + healthy_cols

for col in all_patients:
    clean_data[col] = pd.to_numeric(clean_data[col], errors='coerce')
clean_data = clean_data.dropna(subset=all_patients).copy()

# 3. Math metrics
print("Calculating expressions...")
clean_data['Mean_Alopecia'] = clean_data[alopecia_cols].mean(axis=1)
clean_data['Mean_Healthy'] = clean_data[healthy_cols].mean(axis=1)
clean_data['Fold_Change'] = clean_data['Mean_Alopecia'] - clean_data['Mean_Healthy']

# 4. Pure Vectorized statistical calculations
print("Running t-tests...")
alopecia_mat = clean_data[alopecia_cols].values.astype(float)
healthy_mat = clean_data[healthy_cols].values.astype(float)

_, p_values = ttest_ind(alopecia_mat, healthy_mat, axis=1, equal_var=False)
clean_data['p_value'] = p_values

# 5. Export results file
print("Saving clean data to CSV...")
clean_data.to_csv("alopecia_differential_expression_results.csv")

# 6. Plot generation block (Zero margin anomalies)
print("Plotting Volcano map...")
clean_data['minus_log10_p'] = -np.log10(clean_data['p_value'])

plt.figure(figsize=(10, 7))

# Plot baseline signatures
plt.scatter(clean_data['Fold_Change'], clean_data['minus_log10_p'], color='lightgrey', alpha=0.5, s=2, label='Not Significant')

# Extract up and down expressions explicitly
up = clean_data[(clean_data['p_value'] < 0.05) & (clean_data['Fold_Change'] > 1)]
down = clean_data[(clean_data['p_value'] < 0.05) & (clean_data['Fold_Change'] < -1)]

plt.scatter(up['Fold_Change'], up['minus_log10_p'], color='crimson', alpha=0.8, s=6, label='Significantly Up-regulated')
plt.scatter(down['Fold_Change'], down['minus_log10_p'], color='royalblue', alpha=0.8, s=6, label='Significantly Down-regulated')

# Axis formatting
plt.title("Volcano Plot: Differential Gene Expression in Androgenic Alopecia", fontsize=14, fontweight='bold', pad=15)
plt.xlabel("Log2 Fold Change (Alopecia - Healthy)", fontsize=12)
plt.ylabel("-Log10 P-Value (Statistical Significance)", fontsize=12)
plt.axvline(x=0, color='black', linestyle='--', linewidth=0.8)
plt.axhline(y=-np.log10(0.05), color='darkgray', linestyle=':', linewidth=1.2)
plt.legend(loc='upper right', fontsize=10)
plt.grid(True, linestyle=':', alpha=0.6)

plt.savefig("volcano_plot.png", dpi=300, bbox_inches='tight')
print("Graph saved as 'volcano_plot.png'!")
plt.show()
