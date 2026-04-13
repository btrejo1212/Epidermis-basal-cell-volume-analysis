"""
Tracked cell volume pipeline

What this script does
---------------------
1. Finds tracking folders inside:
       Data/Live volumes from cell tracking/Tracked ex cell 1
2. Reads:
       tracked_IDs.txt
       Denoised_Image Seq/All_Timepoints_Combined.xlsx
3. Builds Tracked_Cell_Shapes.xlsx for each valid folder
4. Finds matching apical / middle / basal folders
5. Builds BMA measurement workbooks
6. Computes volume in um^3
7. Generates simple 3D PNG visualizations

Expected repository layout
--------------------------
repo/
├── Data/
│   └── Live volumes from cell tracking/
│       └── Tracked ex cell 1/
│           ├── Apical Z/
│           │   ├── tracked_IDs.txt
│           │   └── Denoised_Image Seq/
│           │       └── All_Timepoints_Combined.xlsx
│           ├── Middle Z/
│           │   ├── tracked_IDs.txt
│           │   └── Denoised_Image Seq/
│           │       └── All_Timepoints_Combined.xlsx
│           ├── Basal Z/
│           │   ├── tracked_IDs.txt
│           │   └── Denoised_Image Seq/
│           │       └── All_Timepoints_Combined.xlsx
│           └── ... optional additional matching groups ...
├── results/
└── src/
    └── Live volumes from cell tracking/
        └── Live volumes from cell tracking script.py

How to run
----------
Run from the repository root:

    python "src\\Live volumes from cell tracking\\Live volumes from cell tracking script.py"

Notes
-----
- The script auto-detects the repository root when run as a .py file.
- By default it uses:
      Data/Live volumes from cell tracking/Tracked ex cell 1
- Outputs are written under:
      results/Live volumes from cell tracking/Tracked ex cell 1
  and also next to some intermediate files in the data folder.

Editable settings are grouped below.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# USER SETTINGS
# ============================================================

# Optional explicit repo path.
# Leave as None to auto-detect from the script location.
CUSTOM_REPO_DIR: Optional[Path] = None

# Default input folder inside the repo.
DEFAULT_INPUT_SUBDIR = Path("Data") / "Live volumes from cell tracking" / "Tracked ex cell 1"

# Default output folder inside the repo.
DEFAULT_OUTPUT_SUBDIR = Path("results") / "Live volumes from cell tracking" / "Tracked ex cell 1"

# Optional overrides.
# Example:
# CUSTOM_INPUT_DIR = Path(r"D:\my_data\Tracked ex cell 1")
# CUSTOM_OUTPUT_DIR = Path(r"D:\my_outputs\Tracked ex cell 1")
CUSTOM_INPUT_DIR: Optional[Path] = None
CUSTOM_OUTPUT_DIR: Optional[Path] = None

# Imaging metadata
PX_PER_UM = 9.2308   # pixels per micron (Fiji)
Z_STEP_UM = 2.5      # microns between slices

# Columns to extract from All_Timepoints_Combined.xlsx
SHAPE_COLUMNS = ["Area", "Perimeter", "Eccentricity"]


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


def natural_key(value: str):
    """Sort strings using human-friendly number-aware ordering."""
    parts = re.findall(r"\d+|\D+", str(value))
    return [int(p) if p.isdigit() else p.lower() for p in parts]


def normalize_group_name(name: str) -> str:
    """
    Normalize a folder name into a grouping key by removing apical/middle/basal
    and generic separators like 'z', underscores, hyphens, and extra spaces.

    Examples:
    - 'Apical Z'  -> 'main'
    - 'D1 Apical' -> 'd1'
    - 'Cell1 basal z' -> 'cell1'
    """
    s = name.lower()
    s = re.sub(r"\b(apical|middle|basal)\b", " ", s)
    s = re.sub(r"\bz\b", " ", s)
    s = re.sub(r"[_\-]+", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s if s else "main"


def detect_plane_type(name: str) -> Optional[str]:
    """Return 'Apical', 'Middle', or 'Basal' based on a folder name."""
    s = name.lower()
    if "apical" in s:
        return "Apical"
    if "middle" in s:
        return "Middle"
    if "basal" in s:
        return "Basal"
    return None


# ============================================================
# STEP 1: TRACKING EXTRACTION
# ============================================================

def extract_track_records(
    track_ids: List[str],
    sheet_names: List[str],
    time_offset: int,
    all_data: Dict[str, pd.DataFrame],
    shape_columns: List[str],
) -> List[List[object]]:
    """
    Extract one tracked branch into a list of records.

    Each output row contains:
    [timepoint_index, cell_id, shape_1, shape_2, ...]

    Missing cells are recorded as None values.
    """
    records: List[List[object]] = []

    for i, cell_id in enumerate(track_ids):
        tp_index = i + time_offset
        if tp_index >= len(sheet_names):
            break

        tp_name = sheet_names[tp_index]

        if cell_id in ("-", ""):
            records.append([tp_index + 1, None, *([None] * len(shape_columns))])
            continue

        df = all_data[tp_name]
        row = df[df["Cell_ID"] == int(cell_id)]

        if not row.empty:
            values = [row.iloc[0].get(col, None) for col in shape_columns]
            records.append([tp_index + 1, int(cell_id), *values])
        else:
            records.append([tp_index + 1, int(cell_id), *([None] * len(shape_columns))])

    return records


def parse_tracking_line(line: str) -> Optional[Tuple[int, str, List[List[str]]]]:
    """
    Parse one line from tracked_IDs.txt.

    Supported examples:
    - START=23: 101.102.103|201.202.203
    - MyTrack: 11.12.13|21.22.23

    Returns:
    - start_time (1-based)
    - sheet name label
    - list of branches, where each branch is a list of cell IDs
    """
    if ":" not in line:
        return None

    prefix, track_str = line.split(":", 1)
    prefix = prefix.strip()
    track_str = track_str.strip()

    if prefix.startswith("START="):
        start_time = int(prefix.replace("START=", "").strip())
        sheet_name = f"T{start_time}"
    else:
        start_time = 1
        sheet_name = prefix

    branches = [
        part.strip().strip(".").split(".")
        for part in track_str.split("|")
    ]

    return start_time, sheet_name, branches


def folder_has_tracking_inputs(folder: Path) -> bool:
    """
    Return True if the folder contains the expected tracking inputs.
    """
    excel_path = folder / "Denoised_Image Seq" / "All_Timepoints_Combined.xlsx"
    txt_path = folder / "tracked_IDs.txt"
    return excel_path.exists() and txt_path.exists()


def process_tracking_folder(base_folder: Path, shape_columns: List[str]) -> bool:
    """
    Build Tracked_Cell_Shapes.xlsx for one tracking folder.

    Required inputs:
    - tracked_IDs.txt
    - Denoised_Image Seq/All_Timepoints_Combined.xlsx
    """
    excel_path = base_folder / "Denoised_Image Seq" / "All_Timepoints_Combined.xlsx"
    txt_path = base_folder / "tracked_IDs.txt"
    output_path = base_folder / "Tracked_Cell_Shapes.xlsx"

    if not excel_path.exists() or not txt_path.exists():
        return False

    print(f"Processing tracking folder: {base_folder}")

    xls = pd.ExcelFile(excel_path)
    sheet_names = xls.sheet_names
    all_data = {sheet: xls.parse(sheet) for sheet in sheet_names}

    with pd.ExcelWriter(output_path, engine="xlsxwriter") as writer:
        with txt_path.open("r", encoding="utf-8") as f:
            lines = [line.strip() for line in f if line.strip()]

        for raw_line in lines:
            parsed = parse_tracking_line(raw_line)
            if parsed is None:
                continue

            start_time, sheet_name, branches = parsed

            parent_branch = branches[0]
            parent_offset_zero_based = start_time - 1
            daughter_offset_zero_based = parent_offset_zero_based + len(parent_branch)

            merged_df: Optional[pd.DataFrame] = None

            for branch_index, branch_ids in enumerate(branches):
                offset = parent_offset_zero_based if branch_index == 0 else daughter_offset_zero_based

                col_suffix = branch_index + 1
                columns = (
                    [f"Timepoint_T{col_suffix}", f"Cell_ID_T{col_suffix}"]
                    + [f"{col}_T{col_suffix}" for col in shape_columns]
                )

                records = extract_track_records(
                    track_ids=branch_ids,
                    sheet_names=sheet_names,
                    time_offset=offset,
                    all_data=all_data,
                    shape_columns=shape_columns,
                )

                branch_df = pd.DataFrame(records, columns=columns)

                if merged_df is None:
                    merged_df = branch_df
                else:
                    merged_df = pd.concat([merged_df, branch_df], axis=1)

            if merged_df is None:
                continue

            merged_df.to_excel(writer, sheet_name=sheet_name[:31], index=False)

    print(f"  Wrote: {output_path}")
    return True


def find_tracking_folders(parent_dir: Path) -> List[Path]:
    """
    Find all immediate subfolders inside parent_dir that contain tracking inputs.
    """
    folders = []
    for subfolder in sorted(parent_dir.iterdir(), key=lambda p: natural_key(p.name)):
        if subfolder.is_dir() and folder_has_tracking_inputs(subfolder):
            folders.append(subfolder)
    return folders


# ============================================================
# STEP 2: BMA GENERATION
# ============================================================

def read_area_series_from_tracked_shapes(folder: Path) -> pd.Series:
    """
    Read the first sheet from Tracked_Cell_Shapes.xlsx and extract a combined
    Time -> Area series across all track branches found in that sheet.
    """
    path = folder / "Tracked_Cell_Shapes.xlsx"
    if not path.exists():
        raise FileNotFoundError(f"Missing tracked shapes workbook: {path}")

    df = pd.read_excel(path, sheet_name=0)

    pairs = [
        (f"Timepoint_T{i}", f"Area_T{i}")
        for i in range(1, 50)
        if f"Timepoint_T{i}" in df.columns and f"Area_T{i}" in df.columns
    ]

    if not pairs:
        raise ValueError(f"No Timepoint/Area columns found in {path}")

    records = []
    for time_col, area_col in pairs:
        sub = df[[time_col, area_col]].dropna().copy()
        if sub.empty:
            continue

        sub[time_col] = (
            sub[time_col]
            .astype(str)
            .str.extract(r"(\d+)")[0]
            .astype(int)
        )
        sub.columns = ["Time", "Area"]
        records.append(sub)

    if not records:
        raise ValueError(f"No valid area records found in {path}")

    combined = (
        pd.concat(records)
        .drop_duplicates(subset="Time")
        .sort_values("Time")
    )

    return pd.Series(combined["Area"].values, index=combined["Time"].values)


def find_plane_groups(parent_dir: Path) -> Dict[str, Dict[str, Path]]:
    """
    Find plane folder groups inside parent_dir.

    This supports both:
    - simple folders like:
        Apical Z / Middle Z / Basal Z
    - grouped folders like:
        D1 Apical / D1 Middle / D1 Basal
        D2 Apical / D2 Middle / D2 Basal

    Returns a dict like:
        {
            "main": {"Apical": Path(...), "Middle": Path(...), "Basal": Path(...)},
            "d1":   {...},
            "d2":   {...},
        }
    """
    groups: Dict[str, Dict[str, Path]] = {}

    for folder in parent_dir.iterdir():
        if not folder.is_dir():
            continue

        plane_type = detect_plane_type(folder.name)
        if plane_type is None:
            continue

        group_name = normalize_group_name(folder.name)
        groups.setdefault(group_name, {})
        groups[group_name][plane_type] = folder

    return groups


def build_bma_workbooks(parent_dir: Path, results_dir: Path) -> List[Path]:
    """
    Build BMA workbooks for each complete apical/middle/basal group.

    Output files are written to results_dir.
    """
    created_files: List[Path] = []
    groups = find_plane_groups(parent_dir)

    if not groups:
        print("No apical/middle/basal folders were detected.")
        return created_files

    for group_name, planes in groups.items():
        required = {"Apical", "Middle", "Basal"}
        if set(planes.keys()) != required:
            print(f"Skipping group '{group_name}': missing one or more of {sorted(required)}")
            continue

        try:
            apical = read_area_series_from_tracked_shapes(planes["Apical"])
            middle = read_area_series_from_tracked_shapes(planes["Middle"])
            basal = read_area_series_from_tracked_shapes(planes["Basal"])
        except Exception as e:
            print(f"Skipping group '{group_name}': {e}")
            continue

        common_times = sorted(set(apical.index) & set(middle.index) & set(basal.index))
        if not common_times:
            print(f"Skipping group '{group_name}': no overlapping timepoints across apical/middle/basal")
            continue

        df = pd.DataFrame({
            "Timepoint_T": common_times,
            "Basal_Area_px2": [float(basal.loc[t]) for t in common_times],
            "Middle_Area_px2": [float(middle.loc[t]) for t in common_times],
            "Apical_Area_px2": [float(apical.loc[t]) for t in common_times],
        })

        out_name = "BMA_measurements.xlsx" if group_name == "main" else f"{group_name}_BMA_measurements.xlsx"
        out_path = results_dir / out_name
        df.to_excel(out_path, index=False)
        created_files.append(out_path)

        print(f"Wrote BMA workbook: {out_path}")

    return created_files


# ============================================================
# STEP 3: VOLUME CALCULATION + 3D VISUALIZATION
# ============================================================

def calc_volume_um3(row: pd.Series, px_per_um: float, z_step_um: float) -> float:
    """
    Compute volume in um^3 from basal, middle, and apical areas in pixel^2.

    Areas are first converted to um^2 using:
        um^2 = px^2 / (px_per_um^2)

    Then the frustum-style 3-plane formula is applied.
    """
    scale_sq = px_per_um ** 2

    area_basal_um2 = row["Basal_Area_px2"] / scale_sq
    area_middle_um2 = row["Middle_Area_px2"] / scale_sq
    area_apical_um2 = row["Apical_Area_px2"] / scale_sq

    volume_um3 = (z_step_um / 3.0) * (
        area_basal_um2
        + 2 * area_middle_um2
        + area_apical_um2
        + np.sqrt(area_basal_um2 * area_middle_um2)
        + np.sqrt(area_middle_um2 * area_apical_um2)
    )

    return float(volume_um3)


def make_disc(radius_um: float, z_um: float, n_points: int = 50) -> np.ndarray:
    """
    Create a circular ring of points in 3D for a simple shape visualization.
    """
    theta = np.linspace(0, 2 * np.pi, n_points)
    x = radius_um * np.cos(theta)
    y = radius_um * np.sin(theta)
    z = np.full_like(x, z_um)
    return np.vstack((x, y, z)).T


def process_volume_and_make_3d_plots(excel_path: Path, px_per_um: float, z_step_um: float) -> Path:
    """
    Add volume calculations to a BMA workbook and generate 3D PNG images.

    Outputs:
    - sibling workbook: *_with_volume.xlsx
    - sibling image folder: *_3D/
    """
    df = pd.read_excel(excel_path)

    required_cols = {"Basal_Area_px2", "Middle_Area_px2", "Apical_Area_px2"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"{excel_path.name} is missing required columns: {sorted(missing)}")

    df["Volume_um3"] = df.apply(calc_volume_um3, axis=1, args=(px_per_um, z_step_um))

    out_folder = excel_path.with_name(f"{excel_path.stem}_3D")
    out_folder.mkdir(parents=True, exist_ok=True)

    for i, row in df.iterrows():
        radius_basal_um = np.sqrt((row["Basal_Area_px2"] / np.pi)) / px_per_um
        radius_middle_um = np.sqrt((row["Middle_Area_px2"] / np.pi)) / px_per_um
        radius_apical_um = np.sqrt((row["Apical_Area_px2"] / np.pi)) / px_per_um

        verts = np.concatenate([
            make_disc(radius_basal_um, 0.0),
            make_disc(radius_middle_um, z_step_um),
            make_disc(radius_apical_um, 2.0 * z_step_um),
        ])

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111, projection="3d")
        ax.plot_trisurf(verts[:, 0], verts[:, 1], verts[:, 2], alpha=0.7)

        ax.set_xlabel("X (um)")
        ax.set_ylabel("Y (um)")
        ax.set_zlabel("Z (um)")
        ax.set_title(f"{excel_path.stem} | Timepoint {i + 1}")

        png_path = out_folder / f"tp_{i + 1:03d}.png"
        plt.savefig(png_path, dpi=200, bbox_inches="tight")
        plt.close(fig)

    out_excel = excel_path.with_name(f"{excel_path.stem}_with_volume.xlsx")
    df.to_excel(out_excel, index=False)

    print(f"Wrote volume workbook: {out_excel}")
    print(f"Wrote 3D plots to:     {out_folder}")

    return out_excel


# ============================================================
# MASTER PIPELINE
# ============================================================

def run_all(parent_dir: Path, results_dir: Path, px_per_um: float, z_step_um: float) -> None:
    """
    Run the full pipeline on one parent directory.
    """
    print("=" * 70)
    print("TRACKED VOLUME PIPELINE")
    print("=" * 70)
    print(f"Input parent directory:  {parent_dir}")
    print(f"Results directory:       {results_dir}")
    print(f"Pixels per micron:       {px_per_um}")
    print(f"Z step (um):             {z_step_um}")
    print()

    if not parent_dir.exists():
        raise FileNotFoundError(f"Input parent directory does not exist: {parent_dir}")

    if not parent_dir.is_dir():
        raise NotADirectoryError(f"Input path is not a directory: {parent_dir}")

    results_dir.mkdir(parents=True, exist_ok=True)

    print("--- STEP 1: TRACKING EXTRACTION ---")
    tracking_folders = find_tracking_folders(parent_dir)

    if not tracking_folders:
        print("No folders with tracked_IDs.txt + All_Timepoints_Combined.xlsx were found.")
    else:
        for folder in tracking_folders:
            process_tracking_folder(folder, SHAPE_COLUMNS)
    print()

    print("--- STEP 2: BUILD BMA WORKBOOKS ---")
    bma_files = build_bma_workbooks(parent_dir, results_dir)
    if not bma_files:
        print("No BMA workbooks were created.")
    print()

    print("--- STEP 3: COMPUTE VOLUME + MAKE 3D PLOTS ---")
    if not bma_files:
        print("Skipping volume and 3D generation because no BMA files were available.")
    else:
        for bma_path in bma_files:
            process_volume_and_make_3d_plots(
                bma_path,
                px_per_um=px_per_um,
                z_step_um=z_step_um,
            )

    print()
    print("Pipeline complete.")


# ============================================================
# ENTRY POINT
# ============================================================

def main() -> None:
    repo_dir, input_dir, output_dir = resolve_directories()

    print(f"Repository directory:    {repo_dir}")
    run_all(
        parent_dir=input_dir,
        results_dir=output_dir,
        px_per_um=PX_PER_UM,
        z_step_um=Z_STEP_UM,
    )


if __name__ == "__main__":
    main()