# fix_anndata_for_processing.py - Updated version
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad

# Load the fixed file
adata = ad.read_h5ad(
    "/Users/user/Desktop/Backup/R_projects/Shih/series_fixed.h5ad")

print("Original AnnData:")
print(adata)
print(f"\nRaw data shape: {adata.raw.shape if adata.raw else 'No raw data'}")
print(f"Layers: {list(adata.layers.keys()) if adata.layers else 'None'}")

# Check the data range
print("\nData statistics:")
print(f"X min: {adata.X.min():.4f}, max: {adata.X.max():.4f}")
print(f"X contains NaN: {np.isnan(adata.X).any()}")
print(f"X contains Inf: {np.isinf(adata.X).any()}")

# Since the data appears to be SCT normalized (has negative values from scaling)
# and we don't have counts in layers, we need a different approach

# Option 1: If you want to work with just the 3000 variable genes
print("\n=== Option 1: Working with 3000 variable genes ===")

# Create a copy for processing
adata_3000 = adata.copy()

# Since this is SCT data, we need to handle it differently
# SCT data is already normalized and scaled, so we need to tell the processing pipeline this

# Store the scaled data
adata_3000.layers['scaled_data'] = adata_3000.X.copy()

# For Scanpy, we typically want normalized but not scaled data in X
# Since we have scaled data, we'll try to work with it as is
# But we'll add some metadata to indicate the data is already processed

# Mark that variable genes are already selected
adata_3000.var['highly_variable'] = True  # All 3000 genes are variable

# Save this version
output_path_3000 = "/Users/user/Desktop/Backup/R_projects/Shih/series_3000_ready.h5ad"
adata_3000.write_h5ad(output_path_3000)
print(f"Saved 3000-gene version to: {output_path_3000}")

# Option 2: If you want to use all genes from raw
print("\n=== Option 2: Using all genes from raw data ===")

if adata.raw is not None:
    # Create new AnnData from raw
    adata_full = adata.raw.to_adata()

    # Copy over the metadata
    adata_full.obs = adata.obs.copy()

    # Copy over the reductions (but note they were computed on 3000 genes)
    # So they might not be perfectly valid for the full dataset
    print("Note: PCA/UMAP were computed on 3000 genes, not all genes")

    # Check if raw data is counts or normalized
    print(f"\nRaw data statistics:")
    print(
        f"Raw X min: {adata_full.X.min():.4f}, max: {adata_full.X.max():.4f}")
    is_integer = np.all(adata_full.X.data == adata_full.X.data.astype(int)) if hasattr(
        adata_full.X, 'data') else np.all(adata_full.X == adata_full.X.astype(int))
    print(f"Raw X appears to be counts: {is_integer}")

    # If it's a sparse matrix, we might want to keep it sparse
    if hasattr(adata_full.X, 'toarray'):
        print("Raw data is sparse - keeping it sparse for efficiency")

    # Save this version
    output_path_full = "/Users/user/Desktop/Backup/R_projects/Shih/series_full_ready.h5ad"
    adata_full.write_h5ad(output_path_full)
    print(f"Saved full gene version to: {output_path_full}")

# Option 3: Create a version that your app can process without variable gene selection
print("\n=== Option 3: Pre-marked variable genes version ===")

# This is specifically for your app
adata_app = adata.copy()

# Pre-mark the genes so the app doesn't need to calculate them
adata_app.var['highly_variable'] = True
adata_app.var['highly_variable_nbatches'] = 1
adata_app.var['highly_variable_intersection'] = True
adata_app.var['means'] = np.mean(adata_app.X, axis=0)
adata_app.var['dispersions'] = np.std(adata_app.X, axis=0)
adata_app.var['dispersions_norm'] = adata_app.var['dispersions'] / \
    (adata_app.var['means'] + 0.0001)

# Add a flag to indicate this data is already processed
adata_app.uns['data_already_scaled'] = True
adata_app.uns['hvg_already_selected'] = True

output_path_app = "/Users/user/Desktop/Backup/R_projects/Shih/series_app_ready.h5ad"
adata_app.write_h5ad(output_path_app)
print(f"Saved app-ready version to: {output_path_app}")

print("\n=== RECOMMENDATIONS ===")
print("1. Try 'series_app_ready.h5ad' first - it has pre-marked variable genes")
print("2. If that doesn't work, try 'series_3000_ready.h5ad'")
print("3. For full analysis with all genes, use 'series_full_ready.h5ad'")
print("\nYou may also need to modify your app to skip highly_variable_genes calculation")
print("if 'highly_variable' already exists in adata.var")
