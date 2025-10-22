"""
Create a small test AnnData object for testing MASIH.
"""

from scipy.sparse import csr_matrix
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

# Set random seed
np.random.seed(42)

# Create synthetic data
n_obs = 500  # 500 cells
n_vars = 1000  # 1000 genes

# Generate random count matrix (sparse)
X = np.random.negative_binomial(5, 0.3, (n_obs, n_vars))
X = csr_matrix(X.astype(np.float32))

# Create gene names
var_names = [f'Gene_{i}' for i in range(n_vars)]
# Add some mitochondrial genes
var_names[:20] = [f'MT-Gene{i}' for i in range(20)]

# Create cell barcodes
obs_names = [f'Cell_{i}' for i in range(n_obs)]

# Create AnnData object
adata = ad.AnnData(X=X, obs=pd.DataFrame(index=obs_names),
                   var=pd.DataFrame(index=var_names))

# Add some metadata
adata.obs['sample'] = np.random.choice(['SampleA', 'SampleB'], n_obs)

# DON'T add highly_variable column - let MASIH detect it

# Save
adata.write_h5ad('test_data.h5ad')
print(f"✓ Test data created: test_data.h5ad")
print(f"  - Cells: {adata.n_obs}")
print(f"  - Genes: {adata.n_vars}")
