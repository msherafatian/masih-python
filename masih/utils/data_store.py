"""
Data storage utilities for MASIH application.

Handles storage and retrieval of AnnData objects using a disk-based cache.
Since dcc.Store has size limitations, we store actual data on disk and only
keep session IDs and metadata in the browser storage.
"""

import os
import pickle
import uuid
from pathlib import Path
from typing import Optional, Dict, Any
import shutil
from datetime import datetime, timedelta

import anndata as ad

from masih.config.settings import AppConfig


class DataStoreManager:
    """
    Manages storage and retrieval of AnnData objects.

    Uses a disk-based cache in a temporary directory with session-based
    organization. Each session gets a unique ID, and data is stored in
    a subdirectory with that ID.
    """

    def __init__(self, base_dir: Optional[str] = None):
        """
        Initialize the data store manager.

        Args:
            base_dir: Base directory for storing data. If None, uses
                     AppConfig.TEMP_DIR
        """
        self.base_dir = Path(base_dir or AppConfig.TEMP_DIR)
        self.base_dir.mkdir(parents=True, exist_ok=True)

    def generate_session_id(self) -> str:
        """
        Generate a unique session ID.

        Returns:
            UUID string for the session
        """
        return str(uuid.uuid4())

    def get_session_path(self, session_id: str) -> Path:
        """
        Get the directory path for a session.

        Args:
            session_id: Session identifier

        Returns:
            Path object for the session directory
        """
        session_path = self.base_dir / session_id
        session_path.mkdir(parents=True, exist_ok=True)
        return session_path

    def store_adata(self, adata: ad.AnnData, session_id: Optional[str] = None) -> Dict[str, Any]:
        """
        Store an AnnData object to disk.

        Args:
            adata: AnnData object to store
            session_id: Session ID (generates new one if None)

        Returns:
            Dictionary with session metadata to store in dcc.Store
        """
        if session_id is None:
            session_id = self.generate_session_id()

        session_path = self.get_session_path(session_id)
        adata_path = session_path / "adata.h5ad"

        # Save AnnData object
        adata.write_h5ad(adata_path)

        # Create metadata for dcc.Store
        metadata = {
            'session_id': session_id,
            'n_obs': adata.n_obs,
            'n_vars': adata.n_vars,
            'obs_keys': list(adata.obs.columns),
            'var_keys': list(adata.var.columns),
            'obsm_keys': list(adata.obsm.keys()),
            'uns_keys': list(adata.uns.keys()),
            'timestamp': datetime.now().isoformat(),
            'file_path': str(adata_path)
        }

        return metadata

    def load_adata(self, session_id: str) -> Optional[ad.AnnData]:
        """
        Load an AnnData object from disk.

        Args:
            session_id: Session identifier

        Returns:
            AnnData object or None if not found
        """
        session_path = self.get_session_path(session_id)
        adata_path = session_path / "adata.h5ad"

        if not adata_path.exists():
            return None

        try:
            adata = ad.read_h5ad(adata_path)
            return adata
        except Exception as e:
            print(f"Error loading AnnData: {e}")
            return None

    def update_adata(self, session_id: str, adata: ad.AnnData) -> Dict[str, Any]:
        """
        Update an existing AnnData object.

        Args:
            session_id: Session identifier
            adata: Updated AnnData object

        Returns:
            Updated metadata dictionary
        """
        return self.store_adata(adata, session_id=session_id)

    def delete_session(self, session_id: str) -> bool:
        """
        Delete all data for a session.

        Args:
            session_id: Session identifier

        Returns:
            True if successful, False otherwise
        """
        session_path = self.get_session_path(session_id)

        if session_path.exists():
            try:
                shutil.rmtree(session_path)
                return True
            except Exception as e:
                print(f"Error deleting session: {e}")
                return False

        return False

    def cleanup_old_sessions(self, max_age_hours: int = 24) -> int:
        """
        Clean up sessions older than specified age.

        Args:
            max_age_hours: Maximum age in hours before cleanup

        Returns:
            Number of sessions cleaned up
        """
        if not self.base_dir.exists():
            return 0

        cutoff_time = datetime.now() - timedelta(hours=max_age_hours)
        cleaned = 0

        for session_dir in self.base_dir.iterdir():
            if session_dir.is_dir():
                # Check modification time
                mtime = datetime.fromtimestamp(session_dir.stat().st_mtime)
                if mtime < cutoff_time:
                    try:
                        shutil.rmtree(session_dir)
                        cleaned += 1
                    except Exception as e:
                        print(f"Error cleaning up {session_dir}: {e}")

        return cleaned

    def get_session_size(self, session_id: str) -> int:
        """
        Get the total size of a session in bytes.

        Args:
            session_id: Session identifier

        Returns:
            Size in bytes
        """
        session_path = self.get_session_path(session_id)

        if not session_path.exists():
            return 0

        total_size = 0
        for file_path in session_path.rglob('*'):
            if file_path.is_file():
                total_size += file_path.stat().st_size

        return total_size

    def store_additional_data(self, session_id: str, key: str, data: Any) -> bool:
        """
        Store additional data (not AnnData) for a session.

        Useful for storing marker results, pathway scores, etc.

        Args:
            session_id: Session identifier
            key: Name for the data (e.g., 'markers', 'pathways')
            data: Data to store (will be pickled)

        Returns:
            True if successful, False otherwise
        """
        session_path = self.get_session_path(session_id)
        data_path = session_path / f"{key}.pkl"

        try:
            with open(data_path, 'wb') as f:
                pickle.dump(data, f)
            return True
        except Exception as e:
            print(f"Error storing additional data: {e}")
            return False

    def load_additional_data(self, session_id: str, key: str) -> Optional[Any]:
        """
        Load additional data for a session.

        Args:
            session_id: Session identifier
            key: Name of the data

        Returns:
            Loaded data or None if not found
        """
        session_path = self.get_session_path(session_id)
        data_path = session_path / f"{key}.pkl"

        if not data_path.exists():
            return None

        try:
            with open(data_path, 'rb') as f:
                return pickle.load(f)
        except Exception as e:
            print(f"Error loading additional data: {e}")
            return None


# Global instance for easy access throughout the application
data_store = DataStoreManager()


# Convenience functions for common operations
def save_adata(adata: ad.AnnData, session_id: Optional[str] = None) -> Dict[str, Any]:
    """
    Save an AnnData object.

    Args:
        adata: AnnData object
        session_id: Optional session ID

    Returns:
        Metadata dictionary for dcc.Store
    """
    return data_store.store_adata(adata, session_id)


def get_adata(session_id: str) -> Optional[ad.AnnData]:
    """
    Retrieve an AnnData object.

    Args:
        session_id: Session identifier

    Returns:
        AnnData object or None
    """
    return data_store.load_adata(session_id)


def update_adata(session_id: str, adata: ad.AnnData) -> Dict[str, Any]:
    """
    Update an existing AnnData object.

    Args:
        session_id: Session identifier
        adata: Updated AnnData object

    Returns:
        Updated metadata dictionary
    """
    return data_store.update_adata(session_id, adata)


def cleanup_sessions(max_age_hours: int = 24) -> int:
    """
    Clean up old sessions.

    Args:
        max_age_hours: Maximum age before cleanup

    Returns:
        Number of sessions cleaned up
    """
    return data_store.cleanup_old_sessions(max_age_hours)
