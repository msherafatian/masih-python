"""
Tests for data storage utilities.
"""

import pytest
import numpy as np
import anndata as ad
from pathlib import Path
import shutil

from masih.utils.data_store import DataStoreManager, save_adata, get_adata


@pytest.fixture
def temp_store(tmp_path):
    """Create a temporary data store for testing."""
    store = DataStoreManager(base_dir=str(tmp_path))
    yield store
    # Cleanup after test
    if tmp_path.exists():
        shutil.rmtree(tmp_path)


@pytest.fixture
def sample_adata():
    """Create a sample AnnData object for testing."""
    n_obs = 100
    n_vars = 50

    X = np.random.rand(n_obs, n_vars)
    obs = {'cell_type': [f'type_{i % 3}' for i in range(n_obs)]}
    var = {'gene_name': [f'gene_{i}' for i in range(n_vars)]}

    adata = ad.AnnData(X=X, obs=obs, var=var)
    return adata


def test_generate_session_id(temp_store):
    """Test session ID generation."""
    session_id = temp_store.generate_session_id()
    assert isinstance(session_id, str)
    assert len(session_id) > 0

    # Should generate unique IDs
    session_id2 = temp_store.generate_session_id()
    assert session_id != session_id2


def test_store_and_load_adata(temp_store, sample_adata):
    """Test storing and loading AnnData objects."""
    # Store the data
    metadata = temp_store.store_adata(sample_adata)

    assert 'session_id' in metadata
    assert metadata['n_obs'] == 100
    assert metadata['n_vars'] == 50

    # Load the data
    session_id = metadata['session_id']
    loaded_adata = temp_store.load_adata(session_id)

    assert loaded_adata is not None
    assert loaded_adata.n_obs == sample_adata.n_obs
    assert loaded_adata.n_vars == sample_adata.n_vars
    assert list(loaded_adata.obs.columns) == list(sample_adata.obs.columns)


def test_update_adata(temp_store, sample_adata):
    """Test updating an existing AnnData object."""
    # Store initial data
    metadata = temp_store.store_adata(sample_adata)
    session_id = metadata['session_id']

    # Modify the AnnData
    sample_adata.obs['new_column'] = 'test'

    # Update
    updated_metadata = temp_store.update_adata(session_id, sample_adata)

    # Load and verify
    loaded_adata = temp_store.load_adata(session_id)
    assert 'new_column' in loaded_adata.obs.columns


def test_delete_session(temp_store, sample_adata):
    """Test session deletion."""
    metadata = temp_store.store_adata(sample_adata)
    session_id = metadata['session_id']

    # Verify it exists
    assert temp_store.load_adata(session_id) is not None

    # Delete
    result = temp_store.delete_session(session_id)
    assert result is True

    # Verify it's gone
    assert temp_store.load_adata(session_id) is None


def test_store_additional_data(temp_store, sample_adata):
    """Test storing and loading additional data."""
    metadata = temp_store.store_adata(sample_adata)
    session_id = metadata['session_id']

    # Store some additional data
    marker_data = {'cluster_0': ['gene1', 'gene2'],
                   'cluster_1': ['gene3', 'gene4']}
    result = temp_store.store_additional_data(
        session_id, 'markers', marker_data)
    assert result is True

    # Load it back
    loaded_markers = temp_store.load_additional_data(session_id, 'markers')
    assert loaded_markers == marker_data


def test_convenience_functions(sample_adata):
    """Test the convenience functions."""
    # Use default store location
    metadata = save_adata(sample_adata)
    session_id = metadata['session_id']

    # Load it back
    loaded_adata = get_adata(session_id)
    assert loaded_adata is not None
    assert loaded_adata.n_obs == sample_adata.n_obs

    # Cleanup
    from masih.utils.data_store import data_store
    data_store.delete_session(session_id)


def test_get_session_size(temp_store, sample_adata):
    """Test getting session size."""
    metadata = temp_store.store_adata(sample_adata)
    session_id = metadata['session_id']

    size = temp_store.get_session_size(session_id)
    assert size > 0


def test_load_nonexistent_session(temp_store):
    """Test loading a session that doesn't exist."""
    result = temp_store.load_adata('nonexistent_session_id')
    assert result is None
