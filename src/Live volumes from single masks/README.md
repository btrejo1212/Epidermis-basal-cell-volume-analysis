# Live single-cell volume analysis from segmented mask stacks

## Overview

This script measures volume across time for individually tracked live cells using precomputed Cellpose masks.

Each tracked cell is analyzed as a single object across a series of timepoints. For every timepoint, the script reads the segmented Z-stack for the cell of interest and computes a frustum-based volume in `um^3`.

This workflow is intended for datasets where one tracked cell has been isolated and segmented across time, rather than for full-field tissue segmentation.

---

## What this script does

For each tracked cell, the script:

1. Reads timing annotations from:
   - `Birth/`
   - `G1 exit/`
2. Finds the `*_split` folder containing all timepoints
3. Reads the segmented Z slices from each `T#/Zs/` folder
4. Computes per-timepoint frustum volume
5. Writes one Excel workbook summarizing all tracked cells

---

## Important assumptions

This script does **not** perform:

- denoising
- segmentation
- automated tracking
- phase classification

It assumes that all of those steps were completed upstream.

All masks used by this script were generated earlier using **Cellpose**.

Each mask stack should contain only the single tracked cell being analyzed, typically the cell centered in focus.

---

## Dataset structure

Your data should be organized like this:

```text
Data/
  Live volumes from single masks/
    cell 1/
      Birth/
      G1 exit/
      example_split/
        T1/
          Zs/
        T2/
          Zs/
        T3/
          Zs/
        ...
    cell 2/
      Birth/
      G1 exit/
      example_split/
        T1/
          Zs/
        ...