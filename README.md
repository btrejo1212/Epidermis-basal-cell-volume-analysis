# Epidermis basal cell volume analysis

Code and example datasets for calculating cell and nuclear volumes in embryonic epidermal basal cells from live imaging and fixed samples.

## Overview

This repository contains three complementary analysis pipelines for quantifying cell volume in the developing epidermis using live imaging and fixed samples.

The workflows are built around segmentation outputs (Cellpose) and downstream analysis of cell shape and volume across time and space.

The three pipelines included are:

1. **Fixed-cell volume analysis**
   - Uses membrane and nucleus masks from fixed samples
   - Computes volume, shape, and embryo-normalized metrics

2. **Live tracked volume (multi-plane)**
   - Uses basal, middle, and apical measurements from QuantifyPolarity
   - Reconstructs volumes across time using manually tracked cell identities

3. **Live single-cell mask volume**
   - Uses full 3D mask stacks for a single tracked cell
   - Computes volume directly from Z-slice segmentation at each timepoint

Each pipeline is independent and can be run separately.

---

## Repository structure

```text
Epidermis-basal-cell-volume-analysis/
├── Data/
│   ├── Fixed volumes of membrane and nucleus/
│   ├── Live volumes from cell tracking/
│   └── Live volumes from single masks/
│
├── src/
│   ├── Fixed volumes of membrane and nucleus/
│   ├── Live volumes from cell tracking/
│   └── Live volumes from single masks/
│
├── results/        (generated after running scripts)
├── requirements.txt
├── README.md

```

## Pipelines
### 1. Fixed-cell volumes (membrane and nucleus)

<ins>Purpose: </ins>

Quantify cell geometry from fixed samples.

<ins>Input:</ins>

-Cellpose segmentation masks for:\
   -membrane\
   -nucleus

<ins>Output:</ins>

-volume (µm³)\
-height, width, aspect ratio\
-embryo-normalized values


### 2. Live volumes from tracked cell identities

<ins>Purpose: </ins>

Track cells over time and reconstruct volume dynamics using QuantifyPolarity measurements.

<ins>Input:</ins>

-QuantifyPolarity outputs (area measurements)\
-manually created tracked_IDs.txt 

<ins>Output:</ins>

-Frustum model volume calculation by combining three planes over time.

### 3. Live volumes from single-cell mask stacks

<ins>Purpose: </ins>

Compute volume directly from 3D mask stacks for individual tracked cells.

<ins>Input: </ins>

-Cellpose masks for a single cell across time\
-organized as Z-stacks per timepoint

<ins>Output:</ins>

-volume timecourse (µm³) with Birth_T and G1exit_T annotations

## Notes

-All pipelines rely on masks generated using **Cellpose**

-Cell cycle timing (G1, S, G2) are not computed by these scripts, but instead manually annotated using FUCCI reporter and/or Geminin staining.

-Each pipeline python code is run independently from the repository root.

-Datasets are reduced in size for repository storage reasons, but should contain enough information to reproduce some outputs & demonstrate file organization.

-Results folder is not present until python code is ran.

## Requirements:
Dependencies:
```
-numpy<2
-pandas
-matplotlib
-tifffile
-imageio
-openpyxl
-xlsxwriter
-scipy
```
## Suggested Citation

If using this code or analysis approach, please cite appropriately and
reference the associated manuscript.

## Contact

For questions or clarification, please reach out.
