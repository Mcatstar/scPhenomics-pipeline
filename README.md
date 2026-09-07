# scPhenomics Pipeline
[![English](https://img.shields.io/badge/English-blue?style=flat-square)](README.md)[![中文](https://img.shields.io/badge/中文-red?style=flat-square)](./docs/docs/README_zh-cn.md)
 **A generalized single‑cell phenomics image analysis pipeline for multi‑channel fluorescence microscopy**
## Overview
  This repository is a **fork and generalized extension** of [CCRMicroscopyCore/kangs](https://github.com/CCRMicroscopyCore/kangs) (Kangs et al., *Nature Metabolism*, 2025).
⚠️ **Note**: This fork is developed for general‑purpose imaging phenomics analysis and is **not regularly synchronized** with the upstream repository.  
For the paper‑exact original implementation, please refer to the [upstream repository](https://github.com/CCRMicroscopyCore/kangs).
  This pipeline provides a workflow for liver lobule quantification and single‑cell phenomics analysis, adapted from Kang S. et al., *Nature Metabolism* 2025.

## Environment Setup
  ### 1. Install Python 3.13
  Ensure that [Python 3.13](https://www.python.org/) is installed on your system.
  ### 2. Install Poetry (Recommended)
  [Poetry](https://python-poetry.org/docs/) is a dependency management and packaging tool for Python. It allows you to declare project dependencies and manage (install/update) them with ease. Poetry also provides a lockfile to ensure reproducible installations.
If you prefer other tools such as `uv`, `pip`, or `conda`, please refer to the [`pyproject.toml`](pyproject.toml) file for the required library versions.

  ### 3. Install VSCode (Optional but Recommended)
  We recommend installing [VSCode](https://code.visualstudio.com/) along with the **Python language support extension**. For data science workflows, additional extensions can be found in the [VSCode Data Science documentation](https://code.visualstudio.com/docs/datascience/overview).

Detailed installation instructions are not provided here. Please refer to the official documentation for each tool.

##  Clone the Repository
  Fork this repository to your own GitHub account (recommended), or clone it directly:
  ```bash
  git clone https://github.com/Mcatstar/scPhenomics-pipeline.git
  cd scPhenomics-pipeline
  ```
## Create Virtual Environment & Install Dependencies
  ### 1. Configure Poetry to create `.venv` in the project root
  ```bash
  poetry config virtualenvs.in-project true
  ```
  ### 2. List available Python versions
  ```bash
  poetry python list
  ```
  ### 3. Specify the Python version and create the virtual environment
  ```bash
  poetry env use /path/to/python3.13
  ```
  ### 4. Install dependencies
  ```bash
  poetry install --no-root
  ```
⚠️ This is a **script‑only pipeline**, not an installable Python package.
The `--no-root` flag prevents Poetry from installing the current project as a package.

  ###  Open in VSCode
  ```bash
  code .
  ```
## Usage Instructions
  1. Use the Jupyter notebooks in the `notebooks/` directory to **validate your data analysis workflow**.
  2. Once the workflow is confirmed, **write production scripts** in the `src/` folder.
### Image Segmentation & Quantification
  Image segmentation and quantification are performed using **Python 3.13**.
  - Segmentation uses custom [Cellpose](https://github.com/MouseLand/cellpose) models.
  - GPU‑enabled HPC clusters are recommended for segmentation, but not required.
  - For image visualization and verification, we recommend [napari](https://github.com/napari/napari) (multi‑dimensional image viewer). Other tools such as [FIJI](https://fiji.sc/) can also be used.
  #### `lobule_analysis_v2.ipynb`
  This notebook performs the following steps:
  - **Load images**
    Define actin, mitochondria, and lipid channels.
  - **Segment cells**
    Using a combination of actin and mitochondria channels.
    → Manual editing of cell labels is supported.
  - **Segment mitochondria**
  - **Segment lipids**
  - **Create overlap map**
    Overlay of mitochondria and lipid segmentations.
  - **Generate Euclidean Distance Transform (EDT) maps**
    - From manually defined central vein
    - From manually defined portal vein
  - **Generate organelle‑based EDT maps**
    - Based on mitochondria
    - Based on lipids
  - **Quantification**
    For cells, mitochondria, lipids, and organelle overlap regions:
    - Channel intensity
    - Area
    - Centroid position
    - Geometric parameters
### Data Analysis
  Data analysis scripts are written in **R (v4.2)**.
  #### `PLIN5_analysis.R`
  - Loads quantification data
  - Calculates **relative central‑portal vein distance** for cells and organelles, normalized to the maximum CV‑PV distance
  - Bins data into **12 regions (R1–R12)** based on relative distance:
    - **R1** → closer to portal vein (PV)
    - **R12** → closer to central vein (CV)
## License
  This project is licensed under the **MIT License**.
  Original copyright belongs to the original authors. See the [LICENSE](https://license/) file for full details.
## Citation
  If you use this pipeline, please cite:
  - **Original method & pretrained Cellpose weights**:
    Kang S. et al., *Nature Metabolism*, 2025
    and the [upstream repository](https://github.com/CCRMicroscopyCore/kangs).
## Contributing
  Contributions are welcome! Please feel free to open issues or submit pull requests.
## Contact
  For questions or suggestions, please open an issue in this repository.