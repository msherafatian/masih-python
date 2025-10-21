"""
Entry point for running MASIH as a module.

This allows running the application with: python -m masih
"""

from masih.app import run_app

if __name__ == "__main__":
    run_app(debug=True)
