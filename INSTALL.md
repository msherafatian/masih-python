# MASIH-Python Installation Guide

Complete installation instructions for MASIH-Python.

## Table of Contents
1. [System Requirements](#system-requirements)
2. [Installation Methods](#installation-methods)
3. [Docker Installation](#docker-installation)
4. [Troubleshooting](#troubleshooting)
5. [Verification](#verification)

---

## System Requirements

### Minimum Requirements
- **OS**: Windows 10+, macOS 10.15+, or Linux (Ubuntu 20.04+, CentOS 8+)
- **Python**: 3.10 or higher
- **RAM**: 4GB (8GB+ recommended for large datasets)
- **Storage**: 2GB free space
- **Browser**: Chrome, Firefox, Safari, or Edge (latest versions)

### Recommended Requirements
- **RAM**: 16GB+ for datasets with >10,000 cells
- **CPU**: Multi-core processor for faster analysis
- **Storage**: SSD for better performance

---

## Installation Methods

### Method 1: Standard Installation (Recommended for most users)

#### Step 1: Install Python

**Windows:**
1. Download Python from https://www.python.org/downloads/
2. Run installer, **check "Add Python to PATH"**
3. Verify: Open Command Prompt and type `python --version`

**macOS:**
```bash
# Using Homebrew (recommended)
brew install python@3.11

# Or download from python.org
```

**Linux:**
```bash
# Ubuntu/Debian
sudo apt update
sudo apt install python3.11 python3.11-venv python3-pip

# CentOS/RHEL
sudo yum install python311 python311-pip

# Verify
python3 --version
```

#### Step 2: Clone Repository

```bash
# Install git if not already installed
# Windows: https://git-scm.com/download/win
# macOS: brew install git
# Linux: sudo apt install git

# Clone the repository
git clone https://github.com/msherafatian/masih-python.git
cd masih-python
```

#### Step 3: Create Virtual Environment

```bash
# Create virtual environment
python -m venv venv

# Activate it
# On Windows:
venv\Scripts\activate

# On macOS/Linux:
source venv/bin/activate

# You should see (venv) in your prompt
```

#### Step 4: Install Dependencies

```bash
# Upgrade pip
pip install --upgrade pip

# Install all dependencies
pip install -r requirements.txt

# This will take 5-10 minutes depending on your internet speed
```

#### Step 5: Run MASIH

```bash
# Make sure virtual environment is activated
python app.py

# You should see:
# Running on: http://127.0.0.1:8050
```

#### Step 6: Open in Browser

Navigate to: http://127.0.0.1:8050

---

### Method 2: Using Conda (Recommended for bioinformaticians)

#### Step 1: Install Conda

Download and install Miniconda or Anaconda:
- Miniconda (lightweight): https://docs.conda.io/en/latest/miniconda.html
- Anaconda (full): https://www.anaconda.com/products/distribution

#### Step 2: Install MASIH

```bash
# Clone repository
git clone https://github.com/msherafatian/masih-python.git
cd masih-python

# Create conda environment from file
conda env create -f environment.yml

# This creates an environment named 'masih-python'
```

#### Step 3: Activate and Run

```bash
# Activate environment
conda activate masih-python

# Run application
python app.py

# Open browser to http://127.0.0.1:8050
```

#### Step 4: Deactivate When Done

```bash
conda deactivate
```

---

### Method 3: Development Installation

For contributors or those who want to modify the code:

```bash
# Clone and navigate
git clone https://github.com/msherafatian/masih-python.git
cd masih-python

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode
pip install -r requirements-dev.txt
pip install -e .

# This installs the package in "editable" mode
# Changes to code take effect immediately

# Run tests to verify
pytest

# Run the application
python app.py
```

---

## Docker Installation

### Prerequisites
- Docker installed (https://docs.docker.com/get-docker/)
- Docker Compose (included with Docker Desktop)

### Option 1: Using Docker Compose (Easiest)

```bash
# Clone repository
git clone https://github.com/msherafatian/masih-python.git
cd masih-python

# Build and run
docker-compose up -d

# Check status
docker-compose ps

# View logs
docker-compose logs -f masih

# Access at http://localhost:8050
```

### Option 2: Using Docker Directly

```bash
# Clone repository
git clone https://github.com/msherafatian/masih-python.git
cd masih-python

# Build image
docker build -t masih-python:latest .

# Run container
docker run -d \
  --name masih \
  -p 8050:8050 \
  -v $(pwd)/data:/app/data \
  masih-python:latest

# Check status
docker ps

# View logs
docker logs -f masih

# Access at http://localhost:8050
```

### Docker Management

```bash
# Stop container
docker-compose down  # or: docker stop masih

# Restart
docker-compose restart  # or: docker restart masih

# Remove everything
docker-compose down -v
docker rmi masih-python:latest
```

---

## Troubleshooting

### Common Installation Issues

#### Issue 1: "python: command not found"

**Solution:**
```bash
# Try python3 instead
python3 --version
python3 -m venv venv

# Or add Python to PATH (Windows)
# System Properties → Environment Variables → Path → Add Python path
```

#### Issue 2: "pip: command not found"

**Solution:**
```bash
# Install pip
python -m ensurepip --default-pip

# Or use python3
python3 -m pip --version
```

#### Issue 3: Permission Denied (Linux/macOS)

**Solution:**
```bash
# Don't use sudo with pip in virtual environment
# Instead, make sure virtual environment is activated

# If you must install system-wide (not recommended):
pip install --user -r requirements.txt
```

#### Issue 4: "No module named 'scanpy'"

**Solution:**
```bash
# Virtual environment not activated
source venv/bin/activate  # Activate first

# Or reinstall
pip install -r requirements.txt
```

#### Issue 5: Compilation Errors (tables, h5py)

**Solution:**

**On Windows:**
```cmd
# Install Visual C++ Build Tools
# https://visualstudio.microsoft.com/visual-cpp-build-tools/

# Or install pre-compiled wheels
pip install --only-binary :all: tables h5py
```

**On Linux:**
```bash
# Install development packages
sudo apt install build-essential python3-dev libhdf5-dev

# Then reinstall
pip install -r requirements.txt
```

**On macOS:**
```bash
# Install Xcode Command Line Tools
xcode-select --install

# Install HDF5
brew install hdf5

# Then reinstall
pip install -r requirements.txt
```

#### Issue 6: Port 8050 Already in Use

**Solution:**
```bash
# Option 1: Kill process using port
# Windows: netstat -ano | findstr :8050
# Linux/macOS: lsof -ti:8050 | xargs kill

# Option 2: Use different port
export MASIH_PORT=8051  # Linux/macOS
set MASIH_PORT=8051     # Windows
python app.py
```

#### Issue 7: Out of Memory

**Solution:**
```bash
# For large datasets, increase swap or:
# 1. Subsample data before upload
# 2. Use a machine with more RAM
# 3. Process data in batches
```

#### Issue 8: Slow Installation

**Solution:**
```bash
# Use a different PyPI mirror
pip install -r requirements.txt -i https://pypi.tuna.tsinghua.edu.cn/simple

# Or install in stages
pip install numpy pandas scipy
pip install scanpy anndata
pip install dash plotly
pip install -r requirements.txt
```

---

## Verification

### Test Installation

```bash
# Activate environment
source venv/bin/activate  # or: conda activate masih-python

# Check Python version
python --version
# Should be 3.10 or higher

# Check key packages
python -c "import scanpy as sc; print(f'Scanpy {sc.__version__}')"
python -c "import dash; print(f'Dash {dash.__version__}')"
python -c "import plotly; print(f'Plotly {plotly.__version__}')"

# Run tests (if in dev mode)
pytest tests/

# Run the application
python app.py
```

### Expected Output

```
    ╔══════════════════════════════════════════════════════════════╗
    ║                          MASIH                               ║
    ║  Modular Analysis Suite for Interactive Heterogeneity        ║
    ║                      Version 0.1.0                           ║
    ╚══════════════════════════════════════════════════════════════╝
    
    Starting application...
    Running on: http://127.0.0.1:8050
    
    Press CTRL+C to quit
```

### Verify in Browser

1. Open http://127.0.0.1:8050
2. You should see the MASIH interface
3. Try uploading example data (if provided)

---

## Updating MASIH

### Update to Latest Version

```bash
# Navigate to masih-python directory
cd masih-python

# Pull latest changes
git pull origin main

# Update dependencies
source venv/bin/activate  # Activate environment
pip install -r requirements.txt --upgrade

# Or with conda
conda activate masih-python
conda env update -f environment.yml
```

---

## Uninstallation

### Remove Virtual Environment Installation

```bash
# Deactivate environment
deactivate  # or: conda deactivate

# Remove directory
cd ..
rm -rf masih-python  # Linux/macOS
# Windows: rmdir /s masih-python
```

### Remove Conda Installation

```bash
# Remove environment
conda env remove -n masih-python

# Remove directory
rm -rf masih-python
```

### Remove Docker Installation

```bash
docker-compose down -v
docker rmi masih-python:latest
rm -rf masih-python
```

---

## Getting Help

If you encounter issues not covered here:

1. **Check existing issues**: https://github.com/msherafatian/masih-python/issues
2. **Ask a question**: https://github.com/msherafatian/masih-python/discussions
3. **Report a bug**: https://github.com/msherafatian/masih-python/issues/new

When reporting issues, please include:
- Operating system and version
- Python version (`python --version`)
- Complete error message
- Installation method used
- Output of `pip list` or `conda list`

---

## Next Steps

After installation, see:
- [User Guide](docs/user_guide.md) - How to use MASIH
- [API Documentation](docs/api.md) - Programmatic usage
- [Examples](examples/) - Example workflows
- [FAQ](docs/faq.md) - Frequently asked questions

---

**Installation successful? Start analyzing!** 🎉
