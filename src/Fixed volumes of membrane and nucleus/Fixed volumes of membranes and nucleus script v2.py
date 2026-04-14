"""
Fixed membrane and nucleus volume analysis from Cellpose-derived masks.

Overview
--------
This script processes fixed-cell 3D segmentation masks to quantify membrane
and nucleus geometry across multiple embryos, positions, and cell-cycle phases.

All masks used in this analysis were generated prior to this script using
Cellpose segmentation. This script assumes that segmentation has already been
completed and that mask files are correctly organized into the expected folder
structure.

For each cell, the script computes the following metrics separately for:
    (1) membrane
    (2) nucleus

- Frustum volume (um^3)
- Sum volume = sum(slice area) * Z step (um^3)
- Height (um)
- Width (um)
- Aspect ratio = Height / Width
- Number of Z slices used

In addition, the script performs embryo-level normalization by dividing each
cell’s value by the median of that metric across all cells within the same embryo.

---------------------------------------------------------------------

Cell-cycle phase assignment (IMPORTANT)
--------------------------------------
Cell-cycle phase (G1, S, G2) is NOT determined automatically by this script.

Instead, phase labels are assigned manually prior to analysis and are encoded
in the cell folder names, for example:

    G1 cell 1
    S phase cell 2
    G2 cell 3

These labels were assigned by visual inspection of fixed-sample fluorescence,
using markers such as:
    - FUCCI reporters
    - Geminin staining

This script extracts phase information ONLY from folder names and does not
analyze fluorescence intensity or classify phase directly.

Example folders for each phase are included in the dataset for reference.

---------------------------------------------------------------------

Supported dataset structure
---------------------------
This script is designed to handle:

- Multiple embryos
- Multiple positions per embryo
- Multiple cells per position
- Multiple cell-cycle phases (G1, S, G2)

Expected structure:

repo/
├── Data/
│   └── Fixed volumes of membrane and nucleus/
│       ├── Embryo 1/
│       │   ├── Position A/
│       │   │   ├── G1 cell 1/
│       │   │   │   ├── membrane/
│       │   │   │   └── nucleus/   or nuclei/
│       │   │   ├── S phase cell 1/
│       │   │   ├── G2 cell 1/
│       │   │   └── ...
│       │   ├── Position B/
│       │   └── ...
│       ├── Embryo 2/
│       └── ...

Each "cell" folder must contain:
    - a membrane mask folder
    - a nucleus (or nuclei) mask folder

---------------------------------------------------------------------

Mask requirements
-----------------
Accepted mask file types:

- *_cp_masks.png
- *_cp_masks.tif
- *_cp_masks.tiff
- *_seg.npy

These are standard outputs from Cellpose.

If multiple connected components are present in a mask:
    - the largest connected component is used (if SciPy is available)

---------------------------------------------------------------------

Measurement definitions
-----------------------

Volume:
    - Frustum volume is computed using adjacent slice areas
    - Sum volume is computed as (sum of slice areas) × Z step

Width:
    - Defined as the maximum in-plane diameter across all Z slices
    - If SciPy is available:
        uses convex hull distances
    - Otherwise:
        uses bounding-box diagonal (approximation)

Height:
    - Number of Z slices × Z step

Aspect ratio:
    - Height / Width

---------------------------------------------------------------------

Normalization
-------------
For each embryo:
    - all cells across all positions are pooled
    - median values are computed per metric
    - normalized columns are added:

        value_normalized = value / embryo_median

This allows comparison across positions and phases within the same embryo.

---------------------------------------------------------------------

Outputs
-------
For each position, the script writes:

results/
  Fixed volumes of membrane and nucleus/
    <Embryo>/
      <Position>/
        FrustumVolumes_Membrane_Nucleus_byCell.xlsx

Each row corresponds to one cell and includes:
    - metadata (embryo, position, phase, cell number)
    - membrane metrics
    - nucleus metrics
    - normalized metrics

---------------------------------------------------------------------

Notes and limitations
---------------------
- Cell-cycle phase labeling is manual and depends on folder naming
- The script does not classify phase from image data
- Width is a 2D per-slice maximum, not a full 3D diameter
- Results depend on segmentation quality from Cellpose
- Multi-component masks are reduced to the largest component (if SciPy is available)
- Missing or inconsistent mask stacks are skipped with warnings

---------------------------------------------------------------------

How to run
----------
From the repository root:

    python "src\\Fixed volumes of membrane and nucleus\\fixed_membrane_nucleus_volumes.py"

Optional input/output paths can be modified in the USER SETTINGS section below.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple
import re

import numpy as np
import pandas as pd
import tifffile
import imageio.v2 as imageio

# Optional SciPy helpers
try:
    from scipy.ndimage import label as cc_label
    from scipy.spatial import ConvexHull
    from scipy.spatial.distance import pdist
    HAVE_SCIPY = True
except Exception:
    cc_label = None
    ConvexHull = None
    pdist = None
    HAVE_SCIPY = False
    print(
        "⚠️ SciPy not fully available. Width will fall back to bounding-box diagonal, "
        "and connected-component filtering warnings may be skipped."
    )


# ============================================================
# USER SETTINGS
# ============================================================

# Optional explicit repo path.
# Leave as None to auto-detect from the script location.
CUSTOM_REPO_DIR: Optional[Path] = None

# Default input folder inside the repo
DEFAULT_INPUT_SUBDIR = Path("Data") / "Fixed volumes of membrane and nucleus"

# Default output folder inside the repo
DEFAULT_OUTPUT_SUBDIR = Path("results") / "Fixed volumes of membrane and nucleus"

# Optional override paths
# Example:
# CUSTOM_INPUT_DIR = Path(r"D:\my_data\Fixed volumes of membrane and nucleus")
# CUSTOM_OUTPUT_DIR = Path(r"D:\my_outputs\Fixed volumes of membrane and nucleus")
CUSTOM_INPUT_DIR: Optional[Path] = None
CUSTOM_OUTPUT_DIR: Optional[Path] = None

# Calibration
PIXEL_SIZE_UM = 0.1726335   # micrometers per pixel in XY
Z_STEP_UM = 0.3000000       # micrometers between Z slices


# ============================================================
# PATH HELPERS
# ============================================================

def resolve_directories() -> Tuple[Path, Path, Path]:
    """
    Resolve repository, input, and output directories.

    Priority:
    1) explicit overrides
    2) infer repo from script location
    3) fall back to current working directory
    """
    if CUSTOM_REPO_DIR is not None:
        repo_dir = Path(CUSTOM_REPO_DIR).resolve()
    else:
        try:
            repo_dir = Path(__file__).resolve().parents[2]
        except NameError:
            repo_dir = Path.cwd().resolve()

    input_dir = Path(CUSTOM_INPUT_DIR).resolve() if CUSTOM_INPUT_DIR is not None else repo_dir / DEFAULT_INPUT_SUBDIR
    output_dir = Path(CUSTOM_OUTPUT_DIR).resolve() if CUSTOM_OUTPUT_DIR is not None else repo_dir / DEFAULT_OUTPUT_SUBDIR

    output_dir.mkdir(parents=True, exist_ok=True)
    return repo_dir, input_dir, output_dir


# ============================================================
# GENERAL HELPERS
# ============================================================

def natural_key(s: str):
    """Sort helper for mixed text/numeric names."""
    parts = re.findall(r"\d+|\D+", str(s))
    return [int(p) if p.isdigit() else p.lower() for p in parts]


_RX_ZIDX = re.compile(r"[_\-]Z(\d+)", re.IGNORECASE)


def is_allowed_mask(path: Path) -> bool:
    """Return True if path is an accepted Cellpose mask file."""
    suffix = path.suffix.lower()
    stem = path.stem.lower()

    if stem.endswith("_cp_masks") and suffix in (".png", ".tif", ".tiff"):
        return True
    if stem.endswith("_seg") and suffix == ".npy":
        return True
    return False


def extract_z_index(path: Path) -> Optional[int]:
    """Extract Z index from filename if present."""
    match = _RX_ZIDX.search(path.stem)
    return int(match.group(1)) if match else None


def is_cell_folder(path: Path) -> bool:
    """
    Return True if folder name looks like a cell folder.
    Examples:
    - G1 cell 1
    - S phase cell 2
    - G2 cell 3
    """
    if not path.is_dir():
        return False
    return re.search(r"\bcell\s*\d+", path.name, flags=re.IGNORECASE) is not None


def extract_cell_number(name: str) -> Optional[int]:
    """Extract numeric cell ID from folder name."""
    match = re.search(r"\bcell\s*(\d+)", name, flags=re.IGNORECASE)
    return int(match.group(1)) if match else None


def extract_phase_from_name(name: str) -> str:
    """
    Infer phase label from folder name.
    """
    lower = name.lower()
    if "g1" in lower:
        return "G1"
    if "g2" in lower:
        return "G2"
    if "s phase" in lower or lower.startswith("s ") or lower.startswith("s_"):
        return "S"
    return name


# ============================================================
# MASK LOADING
# ============================================================

def read_mask_bool(path: Path) -> np.ndarray:
    """
    Read a mask file and return a boolean 2D mask (True = object).

    Supported:
    - *_cp_masks.(png|tif|tiff)
    - *_seg.npy

    If SciPy is available and there are multiple connected components,
    the largest component is kept.
    """
    suffix = path.suffix.lower()
    stem = path.stem.lower()

    if stem.endswith("_cp_masks") and suffix in (".png", ".tif", ".tiff"):
        if suffix in (".tif", ".tiff"):
            arr = tifffile.imread(str(path))
        else:
            arr = imageio.imread(path)

        if arr.ndim == 3:
            arr = np.mean(arr, axis=2)

        mask = arr > 0

    elif stem.endswith("_seg") and suffix == ".npy":
        try:
            arr = np.load(path, allow_pickle=False)
        except Exception:
            arr = np.load(path, allow_pickle=True)

        if isinstance(arr, dict):
            if "masks" in arr:
                arr = arr["masks"]
            else:
                for key in ("labels", "label", "mask", "arr_0"):
                    if key in arr:
                        arr = arr[key]
                        break

        if isinstance(arr, np.ndarray) and arr.dtype == object and arr.size == 1:
            value = arr.item()
            if isinstance(value, dict) and "masks" in value:
                arr = value["masks"]
            else:
                arr = value

        if isinstance(arr, np.ndarray) and arr.ndim > 2:
            arr = arr.squeeze()

        if not isinstance(arr, np.ndarray):
            raise ValueError(f"Unsupported *_seg.npy structure in {path.name}")

        mask = arr if arr.dtype == bool else (arr.astype(float) > 0)

    else:
        raise ValueError(f"Unsupported mask file: {path.name}")

    mask = mask.astype(bool)

    if cc_label is not None:
        labeled, n_components = cc_label(mask)
        if n_components >= 2:
            sizes = np.bincount(labeled.ravel())
            sizes[0] = 0
            component_sizes = sizes[1:].tolist()
            z_idx = extract_z_index(path)
            print(
                f"⚠️ Multi-component mask in {path.name} "
                f"(Z={z_idx if z_idx is not None else 'NA'}, "
                f"components={n_components}, sizes={component_sizes})"
            )
            largest_label = sizes.argmax()
            mask = labeled == largest_label

    return mask


def collect_masks_zordered(folder: Path) -> Tuple[List[int], List[np.ndarray]]:
    """
    Collect one mask per Z slice from a folder, in Z order.

    Priority per Z:
    1) *_cp_masks.(png|tif|tiff)
    2) *_seg.npy
    """
    if not folder.exists():
        return [], []

    files = [p for p in folder.iterdir() if p.is_file() and is_allowed_mask(p)]
    if not files:
        return [], []

    def rank(p: Path) -> int:
        return 0 if p.stem.lower().endswith("_cp_masks") else 1

    by_index: Dict[int, List[Path]] = {}
    without_index: List[Path] = []

    for path in files:
        z_idx = extract_z_index(path)
        if z_idx is not None:
            by_index.setdefault(z_idx, []).append(path)
        else:
            without_index.append(path)

    chosen: Dict[int, Path] = {}

    if by_index:
        for z_idx, candidates in by_index.items():
            chosen[z_idx] = sorted(candidates, key=lambda p: (rank(p), natural_key(p.name)))[0]
    else:
        ordered = sorted(without_index, key=lambda p: (rank(p), natural_key(p.name)))
        for i, path in enumerate(ordered, start=1):
            chosen[i] = path

    z_indices: List[int] = []
    masks: List[np.ndarray] = []

    for z in sorted(chosen.keys()):
        path = chosen[z]
        try:
            mask = read_mask_bool(path)
            z_indices.append(z)
            masks.append(mask)
        except Exception as e:
            print(f"⚠️ Could not read mask {path}: {e}")

    return z_indices, masks


# ============================================================
# SHAPE / VOLUME CALCULATIONS
# ============================================================

def volume_frustum(areas_um2: np.ndarray, dz_um: float) -> float:
    """
    Compute frustum-based volume from per-slice areas in um^2.
    """
    a = np.asarray(areas_um2, dtype=float)

    if a.size < 2:
        if a.size == 1:
            return float(a.sum() * dz_um)
        return 0.0

    return float(
        np.sum(
            (dz_um / 3.0)
            * (a[:-1] + np.sqrt(np.clip(a[:-1] * a[1:], 0, None)) + a[1:])
        )
    )


def max_diameter_px_single_slice(mask_2d: np.ndarray) -> float:
    """
    Estimate the largest in-plane diameter for a 2D mask.

    Preferred (SciPy available):
    - use convex hull vertices and pairwise distances

    Fallback:
    - use bounding-box diagonal
    """
    ys, xs = np.nonzero(mask_2d)
    if ys.size < 2:
        return 0.0

    coords = np.column_stack((ys, xs))

    if HAVE_SCIPY and ConvexHull is not None and pdist is not None:
        try:
            hull = ConvexHull(coords)
            hull_pts = coords[hull.vertices]
        except Exception:
            hull_pts = coords

        if hull_pts.shape[0] < 2:
            return 0.0

        dists = pdist(hull_pts.astype(float))
        return float(dists.max()) if dists.size > 0 else 0.0

    y_min, y_max = ys.min(), ys.max()
    x_min, x_max = xs.min(), xs.max()
    height_px = y_max - y_min + 1
    width_px = x_max - x_min + 1
    return float(np.hypot(height_px, width_px))


def max_diameter_px_from_stack(mask_stack: np.ndarray) -> float:
    """
    Compute the largest in-plane diameter across all Z slices.
    """
    if mask_stack.ndim != 3:
        mask_stack = np.squeeze(mask_stack)

    if mask_stack.ndim != 3:
        return max_diameter_px_single_slice(mask_stack.astype(bool))

    max_d = 0.0
    for z in range(mask_stack.shape[0]):
        d = max_diameter_px_single_slice(mask_stack[z].astype(bool))
        if d > max_d:
            max_d = d
    return max_d


def find_channel_dirs(cell_dir: Path) -> Tuple[Optional[Path], Optional[Path]]:
    """
    Find membrane and nucleus/nuclei subfolders inside a cell folder.
    """
    membrane_dir = None
    nucleus_dir = None

    for sub in cell_dir.iterdir():
        if not sub.is_dir():
            continue

        name = sub.name.lower()
        if "membrane" in name:
            membrane_dir = sub
        if "nucleus" in name or "nuclei" in name:
            nucleus_dir = sub

    return membrane_dir, nucleus_dir


def compute_channel_metrics(channel_dir: Optional[Path]) -> Tuple[float, float, float, float, float, int]:
    """
    Compute metrics for one channel folder.

    Returns:
    - frustum volume (um^3)
    - sum volume (um^3)
    - height (um)
    - width (um)
    - aspect ratio H/W
    - number of slices
    """
    if channel_dir is None or not channel_dir.exists():
        return np.nan, np.nan, np.nan, np.nan, np.nan, 0

    z_indices, masks = collect_masks_zordered(channel_dir)
    if not masks:
        print(f"      ℹ️ No usable masks in {channel_dir}")
        return np.nan, np.nan, np.nan, np.nan, np.nan, 0

    shapes = {m.shape for m in masks}
    if len(shapes) > 1:
        print(f"      ⚠️ Inconsistent slice shapes in {channel_dir}: {shapes}; skipping.")
        return np.nan, np.nan, np.nan, np.nan, np.nan, 0

    mask_stack = np.stack(masks, axis=0)
    n_slices = mask_stack.shape[0]

    areas_px = mask_stack.reshape(n_slices, -1).sum(axis=1).astype(float)
    areas_um2 = areas_px * (PIXEL_SIZE_UM ** 2)

    frustum_vol_um3 = volume_frustum(areas_um2, Z_STEP_UM)
    sum_vol_um3 = float(areas_um2.sum() * Z_STEP_UM)

    height_um = n_slices * Z_STEP_UM
    width_um = max_diameter_px_from_stack(mask_stack) * PIXEL_SIZE_UM
    aspect_hw = (height_um / width_um) if width_um > 0 else np.nan

    return frustum_vol_um3, sum_vol_um3, height_um, width_um, aspect_hw, n_slices


# ============================================================
# MAIN DATASET PROCESSING
# ============================================================

def process_dataset(root_dir: Path, results_root: Path) -> int:
    """
    Process all embryos and positions in the dataset.

    Writes one Excel file per position to:
        results_root / <Embryo> / <Position> / FrustumVolumes_Membrane_Nucleus_byCell.xlsx

    Returns total number of processed cells.
    """
    embryo_dirs = sorted([d for d in root_dir.iterdir() if d.is_dir()], key=lambda p: natural_key(p.name))

    if not embryo_dirs:
        print(f"⚠️ No embryo folders found in {root_dir}")
        return 0

    total_cells_global = 0

    for embryo_dir in embryo_dirs:
        print("\n==============================")
        print(f"Embryo folder: {embryo_dir.name}")
        print("==============================")

        position_dirs = sorted([d for d in embryo_dir.iterdir() if d.is_dir()], key=lambda p: natural_key(p.name))
        if not position_dirs:
            print(f"  ⚠️ No position folders inside {embryo_dir}")
            continue

        embryo_rows: List[Dict[str, object]] = []

        for pos_dir in position_dirs:
            print(f"\n--- Position folder: {pos_dir.name} ---")

            cell_dirs = sorted([d for d in pos_dir.iterdir() if is_cell_folder(d)], key=lambda p: natural_key(p.name))
            if not cell_dirs:
                print(f"  ⚠️ No 'cell' folders found in {pos_dir}")
                continue

            print(f"  Found {len(cell_dirs)} cell folders in {pos_dir.name}")

            for cell_dir in cell_dirs:
                cell_label = cell_dir.name
                cell_num = extract_cell_number(cell_label)
                phase = extract_phase_from_name(cell_label)

                print(f"    → Processing cell folder: {cell_label} (Phase={phase})")

                membrane_dir, nucleus_dir = find_channel_dirs(cell_dir)

                (
                    mem_vol_frustum,
                    mem_vol_sum,
                    mem_height_um,
                    mem_width_um,
                    mem_aspect,
                    mem_n_slices,
                ) = compute_channel_metrics(membrane_dir)

                (
                    nuc_vol_frustum,
                    nuc_vol_sum,
                    nuc_height_um,
                    nuc_width_um,
                    nuc_aspect,
                    nuc_n_slices,
                ) = compute_channel_metrics(nucleus_dir)

                embryo_rows.append({
                    "Embryo": embryo_dir.name,
                    "Position": pos_dir.name,
                    "Phase": phase,
                    "Cell_Number": cell_num,
                    "Cell_Folder": cell_label,

                    "Membrane_Volume_um3": mem_vol_frustum,
                    "Membrane_Volume_Sum_um3": mem_vol_sum,
                    "Membrane_Height_um": mem_height_um,
                    "Membrane_Width_um": mem_width_um,
                    "Membrane_Aspect_H_over_W": mem_aspect,
                    "Membrane_Num_Z_slices": mem_n_slices,

                    "Nucleus_Volume_um3": nuc_vol_frustum,
                    "Nucleus_Volume_Sum_um3": nuc_vol_sum,
                    "Nucleus_Height_um": nuc_height_um,
                    "Nucleus_Width_um": nuc_width_um,
                    "Nucleus_Aspect_H_over_W": nuc_aspect,
                    "Nucleus_Num_Z_slices": nuc_n_slices,
                })

        if not embryo_rows:
            print(f"  ⚠️ No cells processed in embryo {embryo_dir.name}; skipping.")
            continue

        df_embryo = pd.DataFrame(embryo_rows)

        norm_cols = [
            "Membrane_Volume_um3",
            "Membrane_Volume_Sum_um3",
            "Membrane_Height_um",
            "Membrane_Aspect_H_over_W",
            "Nucleus_Volume_um3",
            "Nucleus_Volume_Sum_um3",
            "Nucleus_Height_um",
            "Nucleus_Aspect_H_over_W",
        ]

        embryo_medians = df_embryo[norm_cols].median(skipna=True)

        for col in norm_cols:
            norm_col = f"{col}_NormByEmbryoMedian"
            denom = embryo_medians.get(col, np.nan)

            if pd.isna(denom) or denom == 0:
                df_embryo[norm_col] = np.nan
                print(
                    f"  ⚠️ Embryo {embryo_dir.name}: median for {col} is NaN or 0; "
                    f"{norm_col} set to NaN."
                )
            else:
                df_embryo[norm_col] = df_embryo[col] / denom

        position_names = sorted(df_embryo["Position"].unique(), key=natural_key)

        for pos_name in position_names:
            df_pos = df_embryo[df_embryo["Position"] == pos_name].copy()

            if df_pos.empty:
                print(f"  ⚠️ No cells for position {pos_name} in {embryo_dir.name}; skipping Excel.")
                continue

            df_pos = df_pos.sort_values(
                by=["Phase", "Cell_Number", "Cell_Folder"],
                ignore_index=True,
            )

            out_dir = results_root / embryo_dir.name / pos_name
            out_dir.mkdir(parents=True, exist_ok=True)

            out_xlsx = out_dir / "FrustumVolumes_Membrane_Nucleus_byCell.xlsx"
            df_pos.to_excel(out_xlsx, index=False)

            print(f"  ✅ Done for position {pos_name} in {embryo_dir.name}")
            print(f"     Wrote Excel: {out_xlsx}")
            print(f"     Total rows (cells) in this position = {len(df_pos)}")

            total_cells_global += len(df_pos)

    return total_cells_global


# ============================================================
# ENTRY POINT
# ============================================================

def main() -> None:
    repo_dir, input_dir, output_dir = resolve_directories()

    print("==================================================")
    print("FIXED MEMBRANE / NUCLEUS VOLUME ANALYSIS")
    print("==================================================")
    print(f"Repository directory: {repo_dir}")
    print(f"Input directory:      {input_dir}")
    print(f"Output directory:     {output_dir}")
    print(f"Pixel size (um):      {PIXEL_SIZE_UM}")
    print(f"Z step (um):          {Z_STEP_UM}")
    print("==================================================")

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")

    total_cells = process_dataset(input_dir, output_dir)

    print("\n========================================")
    print("Finished all embryos / positions.")
    print(f"Total cells processed across dataset = {total_cells}")
    print("========================================")


if __name__ == "__main__":
    main()