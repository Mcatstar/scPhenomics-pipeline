# scPhenomics image analysis, scPhenomics-pipline

### Overview 
General‑purpose single‑cell phenomics pipeline for multi‑channel fluorescence microscopy, forked & generalized from Kangs et al. (Nature Metabolism 2025).
#### original overview:
Image analysis workflow for liver lobule quantification for scPhenomics analysis. Code associated with Kang S, et al. Nature Metabolism 2025.

### Image segmentation and quantification
Image segmentation and quantifcation utilizes Python (v3.10). Custom [Cellpose](https://github.com/MouseLand/cellpose) models were used for segmentation. Segmentation using GPU-enabled HPC clusters is recommended, but is not essential.
For image visualization and verification, [napari](https://github.com/napari/napari) multi-dimensional image viewer was used, but other software such as [FIJI](https://fiji.sc/) can be substituted.

##### lobule_analysis_v2.ipynb
- Load image
  - Define actin, mitochondria, and lipid channels
- Segment cells using a combination of the actin and mitochondria channel
  - Manual-editing of cell labels
- Segment mitochondria
- Segment lipid
- Create an overlap map of mitochondria and lipids
- Create Euclidean distance map from the central and portal veins
  - Manually define the central vein and create EDT map
  - Manually define the portal vein and create EDT map
- Create Euclidean distance map based on organelles
  - Create EDT map based on mitochondria
  - Create EDT map based on lipids
- Quantification of cells, mitochondria, lipids and organelle overlap regions
  - Measure channel intensity, area, centroid position and geometric parameters.

### Data analysis
Data analysis scripts utilize R (v.4.2)

##### PLIN5_analysis.R
- Load quantification data
- Calculate a relative central - portal vein distance for cells and organelles by normalizing to the maximum CV-PV distance
- Based on relative distance, data was binned into 12 regions R1 - R12, with R1 closer to PV and R12 closer to CV
