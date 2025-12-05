# fix_h5ad.py - Fixed version handling groups correctly
import anndata as ad
import pandas as pd
import numpy as np
import os
import h5py


def fix_h5ad_index_issue(input_path, output_path):
    """
    Fix the _index column issue in h5ad files converted from Seurat
    """
    # First check if input file exists
    if not os.path.exists(input_path):
        print(f"Error: Input file does not exist: {input_path}")
        return

    print(f"Loading {input_path}...")
    print(f"File size: {os.path.getsize(input_path)} bytes")

    try:
        with h5py.File(input_path, 'r') as f:
            print("\n=== H5 FILE STRUCTURE ===")
            # First, let's understand the structure

            def print_structure(name, obj):
                if isinstance(obj, h5py.Dataset):
                    print(f"  Dataset: {name}, shape: {obj.shape}")
                else:
                    print(f"  Group: {name}")
                if '_index' in name:
                    print(f"    ⚠️  Found _index at: {name}")
            f.visititems(print_structure)

            print("\n=== EXTRACTING DATA ===")

            # Get X matrix
            X = f['X'][:] if 'X' in f else None
            print(f"X matrix shape: {X.shape if X is not None else 'None'}")

            # Get obs (cell) data
            obs_data = {}
            obs_names = None
            if 'obs' in f:
                obs_group = f['obs']
                if '_index' in obs_group:
                    obs_names = obs_group['_index'][:]
                    if len(obs_names) > 0 and isinstance(obs_names[0], bytes):
                        obs_names = [x.decode('utf-8') for x in obs_names]
                    print(f"Found {len(obs_names)} cell names")

                # Only process datasets, not groups
                for key in obs_group.keys():
                    if key != '_index' and key != '__categories':  # Skip _index and __categories group
                        item = obs_group[key]
                        if isinstance(item, h5py.Dataset):  # Only process datasets
                            data = item[:]
                            if len(data) > 0 and isinstance(data[0], bytes):
                                data = [x.decode('utf-8') for x in data]
                            obs_data[key] = data
                print(f"Obs columns: {list(obs_data.keys())}")

            # Get var (gene) data
            var_data = {}
            var_names = None
            if 'var' in f:
                var_group = f['var']
                if '_index' in var_group:
                    var_names = var_group['_index'][:]
                    if len(var_names) > 0 and isinstance(var_names[0], bytes):
                        var_names = [x.decode('utf-8') for x in var_names]
                    print(f"Found {len(var_names)} gene names")

                # Only process datasets, not groups
                for key in var_group.keys():
                    if key != '_index':
                        item = var_group[key]
                        if isinstance(item, h5py.Dataset):  # Only process datasets
                            data = item[:]
                            if len(data) > 0 and isinstance(data[0], bytes):
                                data = [x.decode('utf-8') for x in data]
                            var_data[key] = data
                print(f"Var columns: {list(var_data.keys())}")

            # Check if raw exists
            raw_X = None
            raw_var_names = None
            if 'raw' in f:
                print("\n⚠️  Found 'raw' slot - checking structure...")
                raw_group = f['raw']

                # Get raw X if it exists
                if 'X' in raw_group:
                    raw_X_group = raw_group['X']
                    if all(k in raw_X_group for k in ['data', 'indices', 'indptr']):
                        print("  Raw X is in sparse format")
                        # We'll handle sparse matrix conversion later if needed

                # Check raw/var for _index
                if 'var' in raw_group and '_index' in raw_group['var']:
                    print(
                        "  ✓ Confirmed: _index found in raw/var - this is the problem!")
                    raw_var_names = raw_group['var']['_index'][:]
                    if len(raw_var_names) > 0 and isinstance(raw_var_names[0], bytes):
                        raw_var_names = [x.decode('utf-8')
                                         for x in raw_var_names]
                    print(f"  Raw has {len(raw_var_names)} gene names")

        # Create clean AnnData
        print("\n=== CREATING CLEAN ANNDATA ===")
        obs_df = pd.DataFrame(obs_data)
        if obs_names is not None:
            obs_df.index = obs_names

        var_df = pd.DataFrame(var_data)
        if var_names is not None:
            var_df.index = var_names

        print(f"Creating AnnData with shape: ({len(obs_df)}, {len(var_df)})")
        adata = ad.AnnData(X=X, obs=obs_df, var=var_df)

        # Copy additional data
        with h5py.File(input_path, 'r') as f:
            # Layers
            if 'layers' in f:
                print("\nCopying layers...")
                for layer_name in f['layers'].keys():
                    layer_data = f['layers'][layer_name][:]
                    adata.layers[layer_name] = layer_data
                    print(f"  - {layer_name}: shape {layer_data.shape}")

            # Obsm (dimensionality reductions)
            if 'obsm' in f:
                print("\nCopying reductions...")
                for key in f['obsm'].keys():
                    obsm_data = f['obsm'][key][:]
                    adata.obsm[key] = obsm_data
                    print(f"  - {key}: shape {obsm_data.shape}")

            # Varm
            if 'varm' in f:
                print("\nCopying varm...")
                for key in f['varm'].keys():
                    varm_data = f['varm'][key][:]
                    adata.varm[key] = varm_data
                    print(f"  - {key}: shape {varm_data.shape}")

            # Handle raw data separately (without _index)
            if 'raw' in f and raw_var_names is not None:
                print("\nProcessing raw data...")
                raw_group = f['raw']

                # Get the sparse matrix components
                if 'X' in raw_group:
                    raw_X_group = raw_group['X']
                    if all(k in raw_X_group for k in ['data', 'indices', 'indptr']):
                        from scipy.sparse import csr_matrix

                        data = raw_X_group['data'][:]
                        indices = raw_X_group['indices'][:]
                        indptr = raw_X_group['indptr'][:]

                        # Create sparse matrix
                        shape = (len(obs_names), len(raw_var_names))
                        raw_X_sparse = csr_matrix(
                            (data, indices, indptr), shape=shape)

                        # Create raw AnnData without _index issues
                        raw_var_df = pd.DataFrame(index=raw_var_names)
                        adata.raw = ad.AnnData(X=raw_X_sparse,
                                               obs=obs_df,
                                               var=raw_var_df)
                        print(
                            f"  Added raw data with shape: {adata.raw.shape}")

        # Save the cleaned file
        print(f"\n=== SAVING CLEANED FILE ===")
        print(f"Output path: {output_path}")
        adata.write_h5ad(output_path)
        print("✓ Successfully saved fixed h5ad file!")

        # Verify it can be loaded
        print("\n=== VERIFYING OUTPUT ===")
        test = ad.read_h5ad(output_path)
        print(f"✓ File can be loaded successfully!")
        print(f"  Shape: {test.shape}")
        print(
            f"  Layers: {list(test.layers.keys()) if test.layers else 'None'}")
        print(
            f"  Reductions: {list(test.obsm.keys()) if test.obsm else 'None'}")
        print(f"  Has raw: {'Yes' if test.raw is not None else 'No'}")

    except Exception as e:
        print(f"\n❌ Error occurred: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    # Paths
    input_file = "/Users/user/Desktop/Backup/R_projects/Shih/series.h5ad"
    output_file = "/Users/user/Desktop/Backup/R_projects/Shih/series_fixed.h5ad"

    print(f"Script will process:")
    print(f"  Input:  {input_file}")
    print(f"  Output: {output_file}")
    print()

    fix_h5ad_index_issue(input_file, output_file)
