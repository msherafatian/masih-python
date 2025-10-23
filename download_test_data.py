"""
Download a small real single-cell dataset for testing MASIH.

This downloads the PBMC3k dataset (3k peripheral blood mononuclear cells)
which is a standard test dataset in single-cell analysis.
"""

import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

print("Downloading PBMC3k dataset...")
print("This may take a few minutes...")

# Download the preprocessed PBMC3k dataset
adata = sc.datasets.pbmc3k()

print(f"\n✓ Dataset downloaded!")
print(f"  - Cells: {adata.n_obs}")
print(f"  - Genes: {adata.n_vars}")

# Basic preprocessing to make it ready for MASIH
print("\nBasic preprocessing...")

# Filter genes
sc.pp.filter_genes(adata, min_cells=3)

print(f"  - After filtering: {adata.n_vars} genes")

# Save
output_file = 'pbmc3k_raw.h5ad'
adata.write(output_file)

print(f"\n✓ Saved to: {output_file}")
print("\nYou can now upload this file to MASIH!")
print("This dataset contains real immune cells with distinct clusters.")
print("Marker genes should be found successfully with this data.")
