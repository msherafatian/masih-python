"""
Comparative pathway analysis utilities for MASIH application.

This module provides functions for comparing multiple CancerSEA pathways
across clusters and calculating correlations.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional
from scipy import stats


def calculate_pathway_correlation(adata, score_columns: Dict[str, str],
                                  method: str = 'pearson') -> Dict:
    """
    Calculate correlation matrix between pathways.

    Args:
        adata: AnnData object
        score_columns: Dictionary mapping pathway names to score column names
        method: Correlation method ('pearson' or 'spearman')

    Returns:
        Dictionary with correlation matrix and p-values
    """
    # Extract score data
    score_data = pd.DataFrame()
    for pathway, score_col in score_columns.items():
        if score_col in adata.obs.columns:
            score_data[pathway] = adata.obs[score_col]
        else:
            print(
                f"Warning: Score column '{score_col}' not found for pathway '{pathway}'")

    print(f"DEBUG: Score data shape: {score_data.shape}")
    print(f"DEBUG: Score data columns: {score_data.columns.tolist()}")

    if score_data.empty or len(score_data.columns) < 2:
        print(
            f"DEBUG: Not enough data for correlation (need 2+ columns, have {len(score_data.columns)})")
        return {'correlation': pd.DataFrame(), 'p_values': pd.DataFrame()}

    # Calculate correlation
    cor_matrix = score_data.corr(method=method)
    print(f"DEBUG: Correlation matrix:\n{cor_matrix}")

    # Calculate p-values
    n = len(score_data)
    p_matrix = pd.DataFrame(
        np.nan, index=cor_matrix.index, columns=cor_matrix.columns)

    for i, col1 in enumerate(score_data.columns):
        for j, col2 in enumerate(score_data.columns):
            if i != j:  # Skip diagonal
                if method == 'pearson':
                    _, p_val = stats.pearsonr(
                        score_data[col1].dropna(),
                        score_data[col2].dropna()
                    )
                else:
                    _, p_val = stats.spearmanr(
                        score_data[col1].dropna(),
                        score_data[col2].dropna()
                    )
                p_matrix.loc[col1, col2] = p_val

    return {
        'correlation': cor_matrix,
        'p_values': p_matrix
    }


def calculate_pathway_by_cluster(adata, score_columns: Dict[str, str],
                                 cluster_key: str = 'leiden') -> pd.DataFrame:
    """
    Calculate mean pathway scores per cluster.

    Args:
        adata: AnnData object
        score_columns: Dictionary mapping pathway names to score column names
        cluster_key: Column containing cluster assignments

    Returns:
        DataFrame with mean scores per cluster
    """
    clusters = sorted(adata.obs[cluster_key].unique())

    results = []
    for cluster in clusters:
        mask = adata.obs[cluster_key] == cluster

        row = {'Cluster': cluster}
        for pathway, score_col in score_columns.items():
            if score_col in adata.obs.columns:
                scores = adata.obs.loc[mask, score_col]
                row[f'{pathway}_mean'] = scores.mean()
                row[f'{pathway}_se'] = scores.sem()

        results.append(row)

    return pd.DataFrame(results)


def calculate_pathway_statistics(adata, score_columns: Dict[str, str]) -> pd.DataFrame:
    """
    Calculate summary statistics for each pathway.

    Args:
        adata: AnnData object
        score_columns: Dictionary mapping pathway names to score column names

    Returns:
        DataFrame with statistics
    """
    stats_list = []

    for pathway, score_col in score_columns.items():
        if score_col in adata.obs.columns:
            scores = adata.obs[score_col].dropna()

            stats_list.append({
                'Pathway': pathway,
                'Mean': scores.mean(),
                'Median': scores.median(),
                'SD': scores.std(),
                'Min': scores.min(),
                'Max': scores.max(),
                'CV': scores.std() / scores.mean() if scores.mean() != 0 else np.nan
            })

    return pd.DataFrame(stats_list)


def create_pathway_profile(adata, score_columns: Dict[str, str],
                           cluster_key: str = 'leiden') -> pd.DataFrame:
    """
    Create pathway activity profile with summary statistics.

    Args:
        adata: AnnData object
        score_columns: Dictionary mapping pathway names to score column names
        cluster_key: Column containing cluster assignments

    Returns:
        DataFrame with pathway profiles per cluster
    """
    profiles = []

    for cluster in sorted(adata.obs[cluster_key].unique()):
        mask = adata.obs[cluster_key] == cluster

        for pathway, score_col in score_columns.items():
            if score_col in adata.obs.columns:
                scores = adata.obs.loc[mask, score_col].dropna()

                profiles.append({
                    'Cluster': cluster,
                    'Pathway': pathway,
                    'Mean': scores.mean(),
                    'SE': scores.sem(),
                    'Median': scores.median(),
                    'Q1': scores.quantile(0.25),
                    'Q3': scores.quantile(0.75)
                })

    return pd.DataFrame(profiles)


def compare_pathway_distributions(adata, score_column: str,
                                  group_var: str,
                                  test_type: str = 'wilcox') -> Dict:
    """
    Compare pathway distributions between groups statistically.

    Args:
        adata: AnnData object
        score_column: Score column name
        group_var: Grouping variable
        test_type: Statistical test ('wilcox', 't', 'kruskal', 'anova')

    Returns:
        Dictionary with test results
    """
    groups = adata.obs[group_var].unique()

    if len(groups) == 2:
        # Two-group comparison
        group1_scores = adata.obs[adata.obs[group_var]
                                  == groups[0]][score_column].dropna()
        group2_scores = adata.obs[adata.obs[group_var]
                                  == groups[1]][score_column].dropna()

        if test_type == 'wilcox':
            stat, p_val = stats.mannwhitneyu(group1_scores, group2_scores)
        elif test_type == 't':
            stat, p_val = stats.ttest_ind(group1_scores, group2_scores)
        else:
            raise ValueError(f"Invalid test type for 2 groups: {test_type}")

        return {
            'test': test_type,
            'statistic': float(stat),
            'p_value': float(p_val),
            'group1': str(groups[0]),
            'group2': str(groups[1]),
            'mean_diff': float(group1_scores.mean() - group2_scores.mean())
        }

    else:
        # Multiple group comparison
        group_scores = [
            adata.obs[adata.obs[group_var] == group][score_column].dropna()
            for group in groups
        ]

        if test_type == 'kruskal':
            stat, p_val = stats.kruskal(*group_scores)
        elif test_type == 'anova':
            stat, p_val = stats.f_oneway(*group_scores)
        else:
            raise ValueError(
                f"Invalid test type for multiple groups: {test_type}")

        return {
            'test': test_type,
            'statistic': float(stat),
            'p_value': float(p_val),
            'n_groups': len(groups)
        }
