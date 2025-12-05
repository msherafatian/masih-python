# MASIH-Python: Modular Analysis Suite for Interactive Heterogeneity <img src="https://raw.githubusercontent.com/msherafatian/masih/main/man/figures/logo.PNG" align="right" width="120" />

[![DOI](https://zenodo.org/badge/1080768114.svg)](https://doi.org/10.5281/zenodo.17824054)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![GitHub issues](https://img.shields.io/github/issues/msherafatian/masih-python.svg)](https://github.com/msherafatian/masih-python/issues)

**A Python/Dash implementation for comprehensive single-cell RNA sequencing analysis, maintaining architectural parity with the [R/Shiny version](https://github.com/msherafatian/masih).**sequencing analysis, maintaining architectural parity with the [R/Shiny version](https://github.com/msherafatian/masih).**

---

## 📋 Table of Contents

- [Overview](#-overview)
- [Key Features](#-key-features)
- [Why MASIH-Python?](#-why-masih-python)
- [Quick Start](#-quick-start)
- [Installation](#-installation)
- [Documentation](#-documentation)
- [Analysis Modules](#-analysis-modules)
- [Supported Data Formats](#-supported-data-formats)
- [Use Cases](#-use-cases)
- [System Requirements](#-system-requirements)
- [Contributing](#-contributing)
- [Citation](#-citation)
- [Roadmap](#-roadmap)
- [Acknowledgments](#-acknowledgments)
- [License](#-license)

---

## 📖 Overview

MASIH-Python provides a user-friendly web interface for analyzing single-cell RNA sequencing (scRNA-seq) data, with a focus on the **multidimensional biological portrait of each cell**.

Built with **Dash** and **Scanpy**, it offers a complete workflow from **raw 10X Genomics data or processed AnnData objects** to **comprehensive biological insights**.

---

## 🌟 Key Features

### 📂 Flexible Input & Preprocessing
- Accepts 10X Genomics outputs, H5AD files, and expression matrices
- Preserves existing analysis results and metadata
- Automatically detects and completes missing analysis steps
- Customizable quality control with interactive parameters

### 🔍 Cellular Landscape Mapping
- **High-Resolution Clustering** — Leiden and Louvain algorithms with resolution optimization
- **Dimensionality Reduction** — PCA, t-SNE, UMAP for revealing data structure
- **Marker Gene Profiling** — Multiple differential expression methods with statistical validation
- **Interactive Visualization** — Plotly-powered plots with real-time exploration

### 🧠 Functional & Dynamic State Inference
- **Functional State Characterization** — CancerSEA-based scoring of 14 cancer-related pathways
- **Trajectory Mapping** — PAGA-based pseudotime analysis for developmental progression
- **Cell Cycle Deconvolution** — G1/S/G2M phase scoring integrated into downstream analyses
- **Comparative Analysis** — Cross-cluster and cross-pathway correlation studies

### 🔤 Insight Sharing & Reporting
- Export publication-quality plots (PNG, PDF, SVG)
- Comprehensive data export (Excel, CSV, H5AD)
- Auto-generate dataset-specific methods text for manuscripts
- Batch export of all visualizations

---

## 🎯 Why MASIH-Python?

- **🐍 Python Ecosystem**: Leverage the power of Scanpy and the Python scientific stack
- **🌐 Web-Based**: No local software installation required (with Docker)
- **🔬 Cancer-Focused**: Specialized tools with CancerSEA pathway integration
- **💻 No Coding Required**: Accessible to all researchers through intuitive interface
- **📊 Comprehensive Workflow**: From raw data to publication-ready figures
- **🔄 Reproducible**: Exports parameters and generates methods text
- **🧩 Modular Design**: Easily extend with new analysis modules
- **🐳 Docker Support**: One-command deployment for easy setup

---

## 🚀 Quick Start

### Installation

```bash
# Option 1: Using Docker (Recommended)
git clone https://github.com/msherafatian/masih-python.git
cd masih-python
docker compose up --build

# Option 2: Using pip
git clone https://github.com/msherafatian/masih-python.git
cd masih-python
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt
python app.py
```

### Launch MASIH-Python

```bash
# If installed locally
python app.py

# If using Docker
docker compose up
```

Access the application at **http://localhost:8050**

### Load Example Data

Download example datasets from the **Upload & QC** tab, or use your own:
- 10X Genomics Cell Ranger output
- AnnData H5AD files
- Expression matrices (CSV/TSV)

---
<!-- 
## 📖 Documentation

### User Guides
- **[Installation Guide](docs/user-guide/installation.md)** - Detailed setup instructions *(coming soon)*
- **[Getting Started](docs/user-guide/getting-started.md)** - Your first analysis *(coming soon)*
- **[Data Preparation](docs/user-guide/data-preparation.md)** - How to prepare your data *(coming soon)*
- **[Troubleshooting](docs/user-guide/troubleshooting.md)** - Common issues and solutions *(coming soon)*

### Tutorials
- **[Basic Analysis](docs/tutorials/basic-analysis.md)** - Complete walkthrough (30 min) *(coming soon)*
- **[Advanced Trajectory Analysis](docs/tutorials/advanced-trajectory.md)** - Pseudotime analysis *(coming soon)*
- **[Comparative Pathway Analysis](docs/tutorials/comparative-pathways.md)** - CancerSEA workflows *(coming soon)*

### For Developers
- **[Contributing Guidelines](CONTRIBUTING.md)** - How to contribute *(coming soon)*
- **[Architecture Overview](docs/developer/architecture.md)** - Technical details *(coming soon)*

---
-->
## 🔬 Analysis Modules

### Core Analysis
- **📤 Data Upload**: Multiple format support (10X, H5AD, CSV/TSV matrices)
- **✅ Quality Control**: Interactive filtering with customizable thresholds
- **🎯 Clustering**: Graph-based clustering with Leiden/Louvain algorithms
- **🧬 Marker Genes**: Statistical testing with multiple methods (Wilcoxon, t-test, logreg)

### Cancer-Specific Features
- **🦠 CancerSEA Integration**: 14 functional state pathways for cancer analysis
- **📊 Pathway Comparison**: Correlation and comparative pathway analysis
- **🔀 Trajectory Analysis**: PAGA-based pseudotime inference
- **🔄 Cell Cycle Scoring**: G1/S/G2M phase identification and integration

### Visualization & Export
- **📈 Interactive Plots**: Real-time exploration with Plotly
- **🖼️ High-Quality Export**: Publication-ready figures (PNG, PDF, SVG)
- **💾 Comprehensive Data Export**: Excel, CSV, and H5AD formats
- **📝 Methods Generation**: Automatic methods text for manuscripts

---

## 📊 Supported Data Formats

- **10X Genomics**: Cell Ranger outputs (matrix directories, H5 files)
- **AnnData Objects**: H5AD files from Scanpy or other Python tools
- **Expression Matrices**: CSV/TSV format (genes × cells or cells × genes)
- **Previous MASIH-Python Sessions**: Reload processed data seamlessly

---

## 🧪 Example Workflows

### Basic Cancer Analysis (30 minutes)
1. 📤 Upload 10X data → 2. ✅ Quality control → 3. 🎯 Clustering → 4. 🦠 CancerSEA analysis → 5. 💾 Export results

### Advanced Trajectory Analysis (45 minutes)
1. 🔄 Basic workflow → 2. 🔍 Cell type selection → 3. 🔀 Trajectory inference → 4. ⏱️ Pseudotime analysis → 5. 📊 Publication figures

### Comparative Study (60 minutes)
1. 📤 Load dataset → 2. 🎯 Clustering analysis → 3. 📊 Comparative pathways → 4. 📈 Statistical analysis → 5. 💾 Export results

---

## 🥼 Use Cases

MASIH-Python is designed for researchers studying:

- **🔬 Tumor Heterogeneity**: Identify and characterize cancer cell subpopulations
- **💊 Treatment Response**: Analyze single-cell responses to therapy
- **📈 Cancer Progression**: Trace developmental trajectories and metastasis
- **⚡ Functional States**: Characterize stemness, invasion, drug resistance
- **🧫 Microenvironment**: Analyze tumor-immune interactions
- **🧬 Gene Expression Patterns**: Discover marker genes and regulatory networks

---

## 📋 System Requirements

### Local Installation
- **Python**: Version 3.10 or higher
- **Operating System**: Windows 10+, macOS 10.14+, or Linux
- **Memory**: 4GB RAM minimum (8GB+ recommended for large datasets)
- **Storage**: 2GB free space for installation and dependencies

### Docker Installation
- **Docker Desktop**: Latest version
- **Memory**: 4GB RAM allocated to Docker (8GB+ recommended)
- **Storage**: 5GB free space for Docker images

---

## 🤝 Contributing

We welcome contributions! Please read our [Contributing Guidelines](CONTRIBUTING.md) for details *(coming soon)*.

### Ways to Contribute
- 🐛 Report bugs and request features
- 📖 Improve documentation
- 🧪 Add new analysis modules
- 🎨 Enhance user interface
- 🧬 Add new pathway databases
- 🐳 Improve Docker deployment

---

## 📄 Citation

If you use MASIH-Python in your research, please cite the Zenodo DOI:

[![DOI](https://zenodo.org/badge/1080768114.svg)](https://doi.org/10.5281/zenodo.17824054)

**Concept DOI (latest MASIH-Python release):**  
https://doi.org/10.5281/zenodo.17824054

**Version-specific DOI (e.g., for the archived release used in a paper):**  
vX.Y.Z → https://doi.org/10.5281/zenodo.17824053

### BibTeX

```bibtex
@software{masih_python_2024,
  author       = {Sherafatian, Masih},
  title        = {MASIH-Python: Modular Analysis Suite for Interactive Heterogeneity},
  year         = 2024,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.17824054},
  url          = {https://doi.org/10.5281/zenodo.17824054}
}

*[Full publication citation will be added upon publication]*

---

## 📞 Support

- **📖 Documentation**: Check our comprehensive guides *(coming soon)*
- **🐛 Issues**: Report bugs on [GitHub Issues](https://github.com/msherafatian/masih-python/issues)
- **✉️ Email**: masihshrftn@gmail.com
- **💬 Discussions**: Join our [GitHub Discussions](https://github.com/msherafatian/masih-python/discussions)

---
<!-- 
## 📈 Roadmap

### Version 0.2.0 (Coming Soon)
- [x] Docker containerization
- [ ] Enhanced documentation and tutorials
- [ ] Command-line interface
- [ ] Batch processing mode
- [ ] Additional clustering algorithms

### Version 0.3.0 (Planned)
- [ ] Cell-cell communication analysis
- [ ] Spatial transcriptomics support
- [ ] Multi-sample integration tools
- [ ] Enhanced pathway databases
- [ ] Cloud deployment options

### Future Features
- [ ] Machine learning cell type prediction
- [ ] Real-time collaboration features
- [ ] REST API for programmatic access
- [ ] Integration with single-cell databases

---
-->
## 🙏 Acknowledgments

MASIH-Python is built on the shoulders of giants:

- **[Scanpy](https://scanpy.readthedocs.io/)** — Single-cell analysis in Python
- **[Dash](https://plotly.com/dash/)** — Web application framework by Plotly
- **[CancerSEA](http://biocc.hrbmu.edu.cn/CancerSEA/)** — Cancer functional state database
- **[CellRank](https://cellrank.readthedocs.io/)** — Trajectory inference
- **[Plotly](https://plotly.com/)** — Interactive visualizations
- **[MASIH R/Shiny](https://github.com/msherafatian/masih)** — Original implementation

Special thanks to all contributors and the single-cell community!

---

## 📜 License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

---

## 🔗 Related Projects

- **[MASIH (R/Shiny)](https://github.com/msherafatian/masih)** — Original R implementation
- **[Scanpy](https://scanpy.readthedocs.io/)** — Single-cell analysis toolkit


---

**MASIH-Python**: Making single-cell cancer analysis accessible to all researchers.

[![Made with ❤️ for Cancer Research](https://img.shields.io/badge/Made%20with%20%E2%9D%A4%EF%B8%8F%20for-Cancer%20Research-red.svg)](https://github.com/msherafatian/masih-python)