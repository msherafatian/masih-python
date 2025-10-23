"""
Trajectory analysis utilities for MASIH application.

This module provides functions for trajectory inference using PAGA
(Partition-based graph abstraction) and diffusion pseudotime.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from typing import Optional, List, Dict, Tuple


def run_trajectory_analysis(adata, cluster_key: str = 'leiden',
                            n_neighbors: int = 15,
                            root_cluster: Optional[str] = None,
                            end_cluster: Optional[str] = None) -> Dict:
    """
    Run trajectory analysis using PAGA and diffusion pseudotime.

    Args:
        adata: AnnData object
        cluster_key: Column containing cluster assignments
        n_neighbors: Number of neighbors for connectivity
        root_cluster: Starting cluster (optional)
        end_cluster: Ending cluster (optional)

    Returns:
        Dictionary with trajectory results
    """
    print("Running trajectory analysis...")

    # CRITICAL: Set matplotlib to non-interactive backend to avoid GUI issues
    import matplotlib
    matplotlib.use('Agg')  # Use non-GUI backend

    # Compute neighborhood graph if not already present
    if 'neighbors' not in adata.uns:
        print(f"Computing neighborhood graph with n_neighbors={n_neighbors}")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X_pca')

    # Run PAGA
    print("Running PAGA...")
    sc.tl.paga(adata, groups=cluster_key)

    # Initialize PAGA positions WITHOUT plotting
    print("Computing PAGA layout...")
    from scipy.sparse import csr_matrix
    import numpy as np

    # Get PAGA connectivity
    connectivities = adata.uns['paga']['connectivities']
    n_groups = connectivities.shape[0]

    # Simple spring layout for PAGA nodes
    np.random.seed(42)
    pos = np.random.random((n_groups, 2))

    # Simple force-directed layout
    for iteration in range(50):
        forces = np.zeros_like(pos)

        # Repulsive forces between all nodes
        for i in range(n_groups):
            for j in range(i + 1, n_groups):
                diff = pos[i] - pos[j]
                dist = np.linalg.norm(diff)
                if dist > 0:
                    force = diff / (dist ** 2 + 0.01)
                    forces[i] += force
                    forces[j] -= force

        # Attractive forces for connected nodes
        rows, cols = connectivities.nonzero()
        for i, j in zip(rows, cols):
            if i != j:
                diff = pos[j] - pos[i]
                dist = np.linalg.norm(diff)
                weight = connectivities[i, j]
                force = diff * weight * dist
                forces[i] += force

        # Update positions
        pos += forces * 0.01

    # Store positions
    adata.uns['paga']['pos'] = pos

    # Determine root cell intelligently
    root_cell_idx = None

    if root_cluster is not None:
        # User specified a root cluster - find a cell in it
        cluster_labels = adata.obs[cluster_key].astype(str)
        root_mask = cluster_labels == str(root_cluster)
        root_indices = np.where(root_mask)[0]

        if len(root_indices) > 0:
            # Pick the cell closest to the cluster center in PCA space
            if 'X_pca' in adata.obsm:
                cluster_pca = adata.obsm['X_pca'][root_mask]
                cluster_center = cluster_pca.mean(axis=0)
                distances = np.linalg.norm(
                    cluster_pca - cluster_center, axis=1)
                closest_idx = distances.argmin()
                root_cell_idx = root_indices[closest_idx]
            else:
                root_cell_idx = root_indices[0]

            print(
                f"Using root cluster: {root_cluster}, root cell index: {root_cell_idx}")
        else:
            print(f"Warning: Root cluster {root_cluster} not found")

    # Auto-detect root if not specified or not found
    if root_cell_idx is None:
        print("Auto-detecting root cell...")

        # Method 1: Use cluster with minimum outgoing connectivity as potential root
        clusters = adata.obs[cluster_key].cat.categories if hasattr(
            adata.obs[cluster_key], 'cat') else sorted(adata.obs[cluster_key].unique())

        # Calculate outgoing connectivity per cluster
        cluster_connectivity = {}
        for i, cluster in enumerate(clusters):
            outgoing = connectivities[i, :].sum(
            ) - connectivities[i, i]  # Exclude self
            cluster_connectivity[cluster] = outgoing

        # Get cluster with minimum outgoing connectivity (likely a start point)
        min_connectivity_cluster = min(
            cluster_connectivity, key=cluster_connectivity.get)
        print(f"Auto-detected root cluster: {min_connectivity_cluster}")

        # Find a cell in this cluster
        cluster_labels = adata.obs[cluster_key].astype(str)
        root_mask = cluster_labels == str(min_connectivity_cluster)
        root_indices = np.where(root_mask)[0]

        if len(root_indices) > 0:
            if 'X_pca' in adata.obsm:
                cluster_pca = adata.obsm['X_pca'][root_mask]
                cluster_center = cluster_pca.mean(axis=0)
                distances = np.linalg.norm(
                    cluster_pca - cluster_center, axis=1)
                closest_idx = distances.argmin()
                root_cell_idx = root_indices[closest_idx]
            else:
                root_cell_idx = root_indices[0]

            print(f"Using auto-detected root cell index: {root_cell_idx}")

    # Set root in adata
    if root_cell_idx is not None:
        adata.uns['iroot'] = root_cell_idx

    # Run diffusion map first (required for DPT)
    print("Computing diffusion map...")
    sc.tl.diffmap(adata, n_comps=10)

    # Run diffusion pseudotime
    print("Computing diffusion pseudotime...")
    if root_cell_idx is not None:
        sc.tl.dpt(adata, n_dcs=10)
    else:
        print("Warning: Could not determine root cell, pseudotime may not be meaningful")
        sc.tl.dpt(adata, n_dcs=10)

    # Verify pseudotime was calculated
    if 'dpt_pseudotime' in adata.obs.columns:
        n_valid = (~adata.obs['dpt_pseudotime'].isna()).sum()
        print(f"✓ Pseudotime calculated for {n_valid}/{adata.n_obs} cells")
    else:
        print("Warning: Pseudotime calculation may have failed")

    # Extract trajectory information
    trajectory_data = {
        'method': 'PAGA + Diffusion Pseudotime',
        'n_neighbors': n_neighbors,
        'root_cluster': root_cluster if root_cluster else str(min_connectivity_cluster) if root_cell_idx else None,
        'end_cluster': end_cluster,
        'root_cell': int(root_cell_idx) if root_cell_idx is not None else None,
        'has_dpt': 'dpt_pseudotime' in adata.obs.columns,
        'has_paga': 'paga' in adata.uns,
        'n_cells': adata.n_obs,
        'clusters': [str(c) for c in adata.obs[cluster_key].unique()]
    }

    print("✓ Trajectory analysis complete")

    return trajectory_data


def prepare_trajectory_subset(adata, cell_types: List[str],
                              annotation_col: str = 'leiden',
                              neotype_col: Optional[str] = None,
                              neotypes: Optional[List[str]] = None):
    """
    Create a subset of data for trajectory analysis.

    Args:
        adata: AnnData object
        cell_types: List of cell types to include
        annotation_col: Column with cell type annotations
        neotype_col: Optional secondary filtering column
        neotypes: Optional list of neotypes to include

    Returns:
        Subset AnnData object
    """
    # Create filter mask
    mask = adata.obs[annotation_col].isin(cell_types)

    # Add neotype filter if provided
    if neotype_col and neotypes:
        neotype_mask = adata.obs[neotype_col].isin(neotypes)
        mask = mask & neotype_mask

    n_cells = mask.sum()
    print(f"Subsetting to {n_cells} cells from {len(cell_types)} cell types")

    if n_cells < 50:
        raise ValueError(
            f"Too few cells selected ({n_cells}). Please select more cell types.")

    # Create subset - preserve the original data structure
    adata_subset = adata[mask].copy()

    print("Preparing subset for trajectory...")

    # Check if data is already processed
    is_logged = False
    if hasattr(adata_subset.X, 'toarray'):
        max_val = adata_subset.X.toarray().max()
    else:
        max_val = adata_subset.X.max()

    # If max value is small (< 20), data is likely already log-normalized
    if max_val < 20:
        is_logged = True
        print("Data appears to be already log-normalized")

    # Only normalize if not already done
    if not is_logged:
        print("Normalizing data...")
        sc.pp.normalize_total(adata_subset, target_sum=1e4)
        sc.pp.log1p(adata_subset)

    # Find highly variable genes with error handling
    try:
        print("Finding highly variable genes...")
        sc.pp.highly_variable_genes(
            adata_subset,
            n_top_genes=min(2000, adata_subset.n_vars),
            flavor='seurat'
        )
    except Exception as e:
        print(f"Warning: Could not find highly variable genes: {e}")
        print("Using all genes instead")
        # If HVG fails, mark all genes as highly variable
        adata_subset.var['highly_variable'] = True

    # Only use highly variable genes for scaling if we found them
    if 'highly_variable' in adata_subset.var.columns:
        n_hvg = adata_subset.var['highly_variable'].sum()
        print(f"Found {n_hvg} highly variable genes")

        if n_hvg > 0:
            # Scale only highly variable genes
            sc.pp.scale(adata_subset, max_value=10)
        else:
            # If no HVGs, scale all genes
            print("No highly variable genes found, scaling all genes")
            sc.pp.scale(adata_subset, max_value=10)
    else:
        # Scale all genes
        sc.pp.scale(adata_subset, max_value=10)

    # Run PCA
    print("Running PCA...")
    n_pcs = min(50, adata_subset.n_obs - 1, adata_subset.n_vars)
    sc.tl.pca(adata_subset, n_comps=n_pcs)

    print("✓ Subset prepared")

    return adata_subset


def get_lineage_paths(adata, cluster_key: str = 'leiden') -> List[List[str]]:
    """
    Extract lineage paths from PAGA results.

    Args:
        adata: AnnData object with PAGA results
        cluster_key: Column containing cluster assignments

    Returns:
        List of lineage paths (each path is a list of cluster names)
    """
    if 'paga' not in adata.uns:
        return []

    # Get connectivity matrix
    connectivities = adata.uns['paga']['connectivities']
    clusters = adata.obs[cluster_key].cat.categories.tolist()

    # Simple lineage extraction: find paths from most connected clusters
    # This is a simplified version - PAGA provides the graph structure
    paths = []

    # Get strongest connections
    threshold = np.percentile(connectivities.data, 75)
    strong_edges = connectivities >= threshold

    # Extract as paths (simplified - could be enhanced)
    for i, cluster in enumerate(clusters):
        connected = [clusters[j]
                     for j in range(len(clusters)) if strong_edges[i, j]]
        if connected:
            paths.append([cluster] + connected)

    return paths


def calculate_trajectory_stats(adata, cluster_key: str = 'leiden',
                               root_cluster: Optional[str] = None,
                               end_cluster: Optional[str] = None) -> Dict:
    """
    Calculate statistics about the trajectory.

    Args:
        adata: AnnData object with trajectory results
        cluster_key: Column containing cluster assignments
        root_cluster: Starting cluster
        end_cluster: Ending cluster

    Returns:
        Dictionary with statistics
    """
    stats = {
        'n_cells': adata.n_obs,
        'n_clusters': len(adata.obs[cluster_key].unique()),
        'clusters': sorted(adata.obs[cluster_key].unique()),
        'root_cluster': root_cluster,
        'end_cluster': end_cluster
    }

    # Pseudotime stats
    if 'dpt_pseudotime' in adata.obs.columns:
        pt = adata.obs['dpt_pseudotime'].values
        pt_valid = pt[~np.isnan(pt)]

        stats['pseudotime_range'] = (
            float(pt_valid.min()), float(pt_valid.max()))
        stats['pseudotime_mean'] = float(pt_valid.mean())
        stats['pseudotime_std'] = float(pt_valid.std())

    # PAGA stats
    if 'paga' in adata.uns:
        connectivities = adata.uns['paga']['connectivities']
        stats['n_connections'] = int((connectivities > 0).sum())
        stats['mean_connectivity'] = float(connectivities.data.mean())

    return stats


def get_pseudotime_by_cluster(adata, cluster_key: str = 'leiden') -> pd.DataFrame:
    """
    Get pseudotime statistics per cluster.

    Args:
        adata: AnnData object
        cluster_key: Column containing cluster assignments

    Returns:
        DataFrame with pseudotime stats per cluster
    """
    if 'dpt_pseudotime' not in adata.obs.columns:
        return pd.DataFrame()

    clusters = sorted(adata.obs[cluster_key].unique())
    stats_list = []

    for cluster in clusters:
        mask = adata.obs[cluster_key] == cluster
        pt = adata.obs.loc[mask, 'dpt_pseudotime'].values
        pt_valid = pt[~np.isnan(pt)]

        if len(pt_valid) > 0:
            stats_list.append({
                'Cluster': cluster,
                'N_cells': mask.sum(),
                'Mean_pseudotime': pt_valid.mean(),
                'Std_pseudotime': pt_valid.std(),
                'Min_pseudotime': pt_valid.min(),
                'Max_pseudotime': pt_valid.max()
            })

    return pd.DataFrame(stats_list)


def export_trajectory_data(adata, trajectory_stats: Dict,
                           cell_types: List[str],
                           annotation_col: str) -> Dict:
    """
    Prepare trajectory data for export.

    Args:
        adata: AnnData object
        trajectory_stats: Trajectory statistics
        cell_types: Selected cell types
        annotation_col: Annotation column name

    Returns:
        Dictionary with exportable data
    """
    export_data = {
        'parameters': {
            'selected_celltypes': cell_types,
            'annotation_column': annotation_col,
            'n_cells': adata.n_obs,
            'method': 'PAGA + Diffusion Pseudotime'
        },
        'statistics': trajectory_stats,
        'pseudotime': adata.obs['dpt_pseudotime'].to_dict() if 'dpt_pseudotime' in adata.obs else {}
    }

    # Add PAGA connectivity
    if 'paga' in adata.uns:
        export_data['paga_connectivities'] = adata.uns['paga']['connectivities'].toarray(
        ).tolist()

    return export_data
