FIXED MEMBRANE AND NUCLEUS VOLUME ANALYSIS

Overview
--------
This script measures membrane and nucleus geometry from fixed-cell 3D masks
that were generated using Cellpose.

The script processes multiple embryos, multiple positions per embryo, and
multiple cells per position, and outputs per-cell measurements.

Each cell must already be segmented. This script does NOT perform segmentation.

---------------------------------------------------------------------

IMPORTANT: Cell-cycle phase assignment
-------------------------------------
Cell-cycle phase (G1, S, G2) is NOT determined automatically.

Phase is assigned manually before running the script and is taken from the
cell folder name, for example:

    G1 cell 1
    S phase cell 2
    G2 cell 3

These labels were determined by visual inspection of fixed samples using:
    - FUCCI reporters
    - Geminin staining

The script simply reads the folder name and assigns the phase accordingly.

---------------------------------------------------------------------

What the script calculates
-------------------------
For each cell (membrane and nucleus separately):

- Frustum volume (um^3)
- Sum volume (sum of slice areas × Z step) (um^3)
- Height (um)
- Width (um)
- Aspect ratio (Height / Width)
- Number of Z slices

Additionally:
- Embryo-level median normalization is applied to key metrics

---------------------------------------------------------------------

Dataset structure (REQUIRED)
---------------------------
Your data must be organized like this:

```text
Data/
  Fixed volumes of membrane and nucleus/
    Embryo 1/
      Position A/
        G1 cell 1/
          membrane/
          nucleus/ or nuclei/
        S phase cell 1/
          membrane/
          nucleus/
        G2 cell 1/
          membrane/
          nucleus/
      Position B/
        ...
    Embryo 2/
      ...

Notes:
- "membrane" folder must contain membrane masks
- "nucleus" or "nuclei" folder must contain nuclear masks
- Folder names must include "cell <number>"

---------------------------------------------------------------------

Mask formats supported
---------------------
The script accepts Cellpose outputs:

- *_cp_masks.png
- *_cp_masks.tif
- *_cp_masks.tiff
- *_seg.npy

If multiple objects are present in a mask:
- the largest connected component is used (if scipy is installed)

---------------------------------------------------------------------

Width definition (IMPORTANT)
---------------------------
Width is defined as the largest in-plane diameter across all Z slices.

- If scipy is installed:
    uses convex hull distances (more accurate)
- Otherwise:
    uses bounding box diagonal (approximation)

This is NOT a full 3D diameter.

---------------------------------------------------------------------

Normalization
-------------
For each embryo:

- All cells across all positions are pooled
- Median values are computed per metric
- Each cell is normalized:

    normalized_value = value / embryo_median

This allows comparison across positions and phases within the same embryo.

---------------------------------------------------------------------

Output
------
For each position, the script writes:

results/
  Fixed volumes of membrane and nucleus/
    <Embryo>/
      <Position>/
        FrustumVolumes_Membrane_Nucleus_byCell.xlsx

Each row = one cell

Includes:
- metadata (embryo, position, phase, cell number)
- membrane measurements
- nucleus measurements
- normalized values

---------------------------------------------------------------------

Installation
------------
Create a clean environment (recommended), then install:

    pip install -r requirements.txt

If you do not use requirements.txt:

    pip install numpy pandas tifffile imageio openpyxl scipy matplotlib

IMPORTANT:
numpy is restricted to <2 to avoid compatibility issues.

---------------------------------------------------------------------

How to run
----------
From the repository root:

    python "src\\Fixed volumes of membrane and nucleus\\fixed_membrane_nucleus_volumes.py"

Do NOT run from inside the src folder.

---------------------------------------------------------------------

Notes / Limitations
------------------
- Phase labels are manual and depend on folder naming
- Script does not analyze fluorescence to determine phase
- Width is 2D-based, not true 3D diameter
- Results depend on segmentation quality from Cellpose
- Missing or inconsistent masks are skipped with warnings

---------------------------------------------------------------------

Recommended naming
-----------------
Use consistent naming like:

- Embryo 1
- Position A
- G1 cell 1
- S phase cell 1
- G2 cell 1
- membrane
- nucleus

---------------------------------------------------------------------

If something fails
-----------------
Check:
- folder structure matches expected layout
- mask files exist and are readable
- dependencies are installed correctly
- numpy version is <2

---------------------------------------------------------------------
