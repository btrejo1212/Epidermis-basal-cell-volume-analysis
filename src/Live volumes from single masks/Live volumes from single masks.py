"""
Compute frustum-based cell volumes over time from Cellpose mask files.

Default behavior:
- Runs on the example dataset included in the repository.

Optional behavior:
- User can override the input and output directories below to run on their own data.

Expected input structure:
- ROOT_DIR contains one or more cell folders
- each cell folder may contain:
    - a Birth folder
    - a G1 exit folder
    - one *_split folder
- inside the *_split folder:
    - timepoint folders such as T1, T2, ...
- inside each timepoint folder:
    - a Zs/ folder containing mask files

Accepted mask file types:
- *_cp_masks.png
- *_cp_masks.tif
- *_cp_masks.tiff
- *_seg.npy

Output:
- one Excel workbook summarizing volume timecourses
"""

from pathlib import Path
import re
import numpy as np
import pandas as pd
import tifffile
import imageio.v2 as imageio


# ============================================================
# USER SETTINGS
# ============================================================

# Default: use example data inside the repo
DEFAULT_INPUT_SUBDIR = Path("data") / "live_single_masks_example"
DEFAULT_OUTPUT_SUBDIR = Path("results")

# Optional override:
# Leave as None to use the example data/results folders above.
# To run on your own data, replace None with a Path(...), e.g.
# CUSTOM_INPUT_DIR = Path(r"D:\my_data\live_single_masks")
# CUSTOM_OUTPUT_DIR = Path(r"D:\my_outputs")
CUSTOM_INPUT_DIR = None
CUSTOM_OUTPUT_DIR = None

# Imaging metadata
PIXEL_SIZE_UM = 0.11   # micrometers per pixel in XY
Z_STEP_UM = 2.5        # micrometers between Z slices


# ============================================================
# HELPERS
# ============================================================

def _natural_key(s: str):
    """Sort helper for strings with mixed text and numbers."""
    parts = re.findall(r"\d+|\D+", str(s))
    return [int(p) if p.isdigit() else p.lower() for p in parts]


# T in file/folder names, e.g. "t14" or "T 9"
_RX_T_ANY = re.compile(r"(?<!\w)[tT]\s*(\d+)\b")

# Z index in file names, e.g. "_Z002_" or "-Z015"
_RX_ZIDX = re.compile(r"[_\-]Z(\d+)", re.IGNORECASE)


def _is_allowed_mask(p: Path) -> bool:
    """Return True if path is an allowed Cellpose mask file."""
    suf = p.suffix.lower()
    stem = p.stem.lower()

    if stem.endswith("_cp_masks") and suf in (".png", ".tif", ".tiff"):
        return True
    if stem.endswith("_seg") and suf == ".npy":
        return True
    return False


def _z_index(p: Path) -> int | None:
    """Extract Z index from filename, or return None if absent."""
    m = _RX_ZIDX.search(p.stem)
    return int(m.group(1)) if m else None


def _read_mask_bool(path: Path) -> np.ndarray:
    """
    Read a mask file and return a boolean 2D mask (True = cell).

    Supports:
    - *_cp_masks.(png|tif|tiff)
    - *_seg.npy

    If multiple connected components are present, the largest is kept.
    """
    suf = path.suffix.lower()
    stem = path.stem.lower()

    if stem.endswith("_cp_masks") and suf in (".png", ".tif", ".tiff"):
        if suf in (".tif", ".tiff"):
            arr = tifffile.imread(str(path))
        else:
            arr = imageio.imread(path)

        if arr.ndim == 3:
            arr = np.mean(arr, axis=2)

        mask = arr > 0

    elif stem.endswith("_seg") and suf == ".npy":
        try:
            arr = np.load(path, allow_pickle=False)
        except Exception:
            arr = np.load(path, allow_pickle=True)

        # Unwrap common Cellpose / saved-dict structures
        if isinstance(arr, dict):
            if "masks" in arr:
                arr = arr["masks"]
            else:
                for k in ("labels", "label", "mask", "arr_0"):
                    if k in arr:
                        arr = arr[k]
                        break

        if isinstance(arr, np.ndarray) and arr.dtype == object and arr.size == 1:
            v = arr.item()
            if isinstance(v, dict) and "masks" in v:
                arr = v["masks"]
            else:
                arr = v

        if isinstance(arr, np.ndarray) and arr.ndim > 2:
            arr = arr.squeeze()

        if not isinstance(arr, np.ndarray):
            raise ValueError(f"Unsupported *_seg.npy structure in {path.name}")

        mask = arr if arr.dtype == bool else (arr.astype(float) > 0)

    else:
        raise ValueError(f"Unsupported mask file: {path.name}")

    # Keep largest connected component if multiple are present
    try:
        from scipy.ndimage import label as cc_label
        labeled, n = cc_label(mask.astype(bool))
        if n >= 2:
            sizes = np.bincount(labeled.ravel())
            sizes[0] = 0
            mask = labeled == sizes.argmax()
        else:
            mask = mask.astype(bool)
    except Exception:
        mask = mask.astype(bool)

    return mask


def _volume_frustum(areas_um2: np.ndarray, dz_um: float) -> float:
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


def _find_phase_dir(cell_dir: Path, phase: str) -> Path | None:
    """
    Find Birth or G1 exit folder inside a cell directory.

    phase must be one of:
    - "birth"
    - "g1 exit"
    """
    target = phase.lower().replace("_", " ").strip()

    for sub in cell_dir.iterdir():
        if not sub.is_dir():
            continue

        s = sub.name.lower().replace("_", " ").strip()

        if target == "birth":
            if s.startswith("birth"):
                return sub
        elif target == "g1 exit":
            if "g1" in s and "exit" in s:
                return sub

    return None


def _extract_phase_T(phase_dir: Path | None, is_birth: bool) -> int | None:
    """
    Extract T indices from files inside a phase directory.

    Rules:
    - Birth uses min(T)
    - G1 exit uses max(T)
    """
    if phase_dir is None or not phase_dir.exists():
        return None

    Ts = set()
    for p in phase_dir.iterdir():
        if p.is_file():
            m = _RX_T_ANY.search(p.stem)
            if m:
                Ts.add(int(m.group(1)))

    if not Ts:
        return None

    return min(Ts) if is_birth else max(Ts)


def _find_split_dir(cell_dir: Path) -> Path | None:
    """
    Find the *_split folder inside a cell directory.
    """
    for sub in cell_dir.iterdir():
        if sub.is_dir() and sub.name.lower().endswith("_split"):
            return sub
    return None


def _collect_masks_zordered(zs_dir: Path):
    """
    Collect one mask per Z slice from a Zs/ folder, in Z order.

    Priority per Z:
    1) *_cp_masks.(png|tif|tiff)
    2) *_seg.npy

    Returns:
    - z_indices: list[int]
    - ok_masks: list[np.ndarray]
    """
    if not zs_dir.exists():
        return [], []

    files = [p for p in zs_dir.iterdir() if p.is_file() and _is_allowed_mask(p)]
    if not files:
        return [], []

    def rank(pp: Path) -> int:
        return 0 if pp.stem.lower().endswith("_cp_masks") else 1

    by_idx = {}
    others = []

    for p in files:
        idx = _z_index(p)
        if idx is not None:
            by_idx.setdefault(idx, []).append(p)
        else:
            others.append(p)

    chosen = {}

    if by_idx:
        for idx, plist in by_idx.items():
            chosen[idx] = sorted(plist, key=lambda q: (rank(q), _natural_key(q.name)))[0]
    else:
        ordered = sorted(others, key=lambda q: (rank(q), _natural_key(q.name)))
        for k, p in enumerate(ordered, start=1):
            chosen[k] = p

    z_indices = []
    ok_masks = []

    for z in sorted(chosen.keys()):
        p = chosen[z]
        try:
            mask = _read_mask_bool(p)
            z_indices.append(z)
            ok_masks.append(mask)
        except Exception as e:
            print(f"⚠️ Could not read mask {p}: {e}")

    return z_indices, ok_masks


def _is_cell_folder(p: Path) -> bool:
    """
    Return True if directory name looks like a cell folder.
    """
    if not p.is_dir():
        return False
    return re.search(r"\bcell\s*\d+", p.name, flags=re.IGNORECASE) is not None


def _resolve_directories():
    """
    Resolve input and output directories.

    Default:
    - use example input inside the repo
    - use results/ inside the repo

    Optional:
    - use CUSTOM_INPUT_DIR / CUSTOM_OUTPUT_DIR if provided
    """
    repo_dir = Path(__file__).resolve().parents[2]

    root_dir = Path(CUSTOM_INPUT_DIR) if CUSTOM_INPUT_DIR is not None else repo_dir / DEFAULT_INPUT_SUBDIR
    results_dir = Path(CUSTOM_OUTPUT_DIR) if CUSTOM_OUTPUT_DIR is not None else repo_dir / DEFAULT_OUTPUT_SUBDIR

    results_dir.mkdir(parents=True, exist_ok=True)

    return repo_dir, root_dir, results_dir


# ============================================================
# MAIN
# ============================================================

def main():
    repo_dir, root_dir, results_dir = _resolve_directories()

    print("=== Live single-mask volume analysis ===")
    print(f"Repository directory: {repo_dir}")
    print(f"Input directory:      {root_dir}")
    print(f"Output directory:     {results_dir}")
    print(f"Pixel size (um):      {PIXEL_SIZE_UM}")
    print(f"Z step (um):          {Z_STEP_UM}")

    if not root_dir.exists():
        raise FileNotFoundError(
            f"Input directory does not exist:\n{root_dir}\n\n"
            f"Either:\n"
            f"1) place the example dataset in the expected repo location, or\n"
            f"2) set CUSTOM_INPUT_DIR near the top of this script."
        )

    cells_meta = []
    max_T_global = 0
    total_slices_analyzed = 0

    cell_dirs = [d for d in root_dir.iterdir() if _is_cell_folder(d)]
    cell_dirs = sorted(cell_dirs, key=lambda p: _natural_key(p.name))

    if not cell_dirs:
        print(f"⚠️ No 'cell ...' folders found in {root_dir}")

    for cdir in cell_dirs:
        descriptor = cdir.name
        print(f"\n===== Processing cell folder: {descriptor} =====")

        # Birth / G1 exit T detection
        birth_dir = _find_phase_dir(cdir, "birth")
        g1_dir = _find_phase_dir(cdir, "g1 exit")

        birth_T = _extract_phase_T(birth_dir, is_birth=True)
        g1_T = _extract_phase_T(g1_dir, is_birth=False)

        if birth_T is None:
            print(f"  ℹ️ No Birth T found in {birth_dir}")
        else:
            print(f"  • Birth T = {birth_T}")

        if g1_T is None:
            print(f"  ℹ️ No G1 exit T found in {g1_dir}")
        else:
            print(f"  • G1 exit T = {g1_T}")

        # Find *_split folder containing timecourse folders
        split_dir = _find_split_dir(cdir)
        if split_dir is None:
            print(f"  ⚠️ No *_split folder found in {descriptor}, skipping timecourse.")
            cells_meta.append({
                "name": descriptor,
                "birth_t": birth_T,
                "g1_t": g1_T,
                "volumes": {},
            })
            continue

        print(f"  • Using split folder: {split_dir.name}")

        volumes_by_T = {}

        t_dirs = [d for d in split_dir.iterdir() if d.is_dir()]
        t_items = []

        for tdir in t_dirs:
            m = re.search(r"[tT]\s*(\d+)", tdir.name)
            if m:
                t_idx = int(m.group(1))
                t_items.append((t_idx, tdir))

        t_items = sorted(t_items, key=lambda x: x[0])

        if not t_items:
            print(f"  ⚠️ No T# folders found in {split_dir}")
        else:
            for t_idx, tdir in t_items:
                zs_dir = tdir / "Zs"
                if not zs_dir.exists():
                    print(f"    ℹ️ No Zs/ folder in {tdir.name}; skipping this T.")
                    continue

                z_indices, masks = _collect_masks_zordered(zs_dir)
                if not masks:
                    print(f"    ℹ️ No usable masks in {zs_dir}; skipping T{t_idx}.")
                    continue

                shapes = {m.shape for m in masks}
                if len(shapes) > 1:
                    print(f"    ⚠️ Inconsistent slice shapes in {zs_dir}: {shapes}; skipping this T.")
                    continue

                mask_stack = np.stack(masks, axis=0)  # [Z, Y, X]
                num_slices = mask_stack.shape[0]
                total_slices_analyzed += num_slices

                areas_px = mask_stack.reshape(mask_stack.shape[0], -1).sum(axis=1).astype(float)
                areas_um2 = areas_px * (PIXEL_SIZE_UM ** 2)

                frustum_vol = _volume_frustum(areas_um2, Z_STEP_UM)
                volumes_by_T[t_idx] = frustum_vol
                max_T_global = max(max_T_global, t_idx)

                print(
                    f"    T{t_idx:02d}: slices={num_slices}  "
                    f"Zs={z_indices or list(range(1, num_slices + 1))}  "
                    f"Frustum volume={frustum_vol:.1f} um^3"
                )

        cells_meta.append({
            "name": descriptor,
            "birth_t": birth_T,
            "g1_t": g1_T,
            "volumes": volumes_by_T,
        })

    # Build Excel table
    if not cells_meta:
        print("⚠️ No cells processed; nothing to write.")
        return

    if max_T_global <= 0:
        print("⚠️ No timepoints with volumes found; writing metadata only.")
        max_T_global = 0

    time_indices = list(range(1, max_T_global + 1)) if max_T_global > 0 else []
    cell_names = [c["name"] for c in cells_meta]

    if time_indices:
        df_vol = pd.DataFrame(index=time_indices, columns=cell_names, dtype=float)
        for meta in cells_meta:
            name = meta["name"]
            for t_idx, vol in meta["volumes"].items():
                if t_idx in df_vol.index:
                    df_vol.at[t_idx, name] = vol
    else:
        df_vol = pd.DataFrame(columns=cell_names, dtype=float)

    out_xlsx = results_dir / "Cell_Frustum_Volumes_Timecourse.xlsx"

    with pd.ExcelWriter(out_xlsx, engine="xlsxwriter") as writer:
        df_vol.to_excel(
            writer,
            sheet_name="FrustumVolumes",
            startrow=3,
            startcol=1,
            header=False,
            index=False,
        )

        worksheet = writer.sheets["FrustumVolumes"]

        # Row 1: cell names; Row 2: Birth_T; Row 3: G1exit_T
        worksheet.write(0, 0, "CellName")
        worksheet.write(1, 0, "Birth_T")
        worksheet.write(2, 0, "G1exit_T")

        for j, meta in enumerate(cells_meta):
            col = 1 + j
            worksheet.write(0, col, meta["name"])
            worksheet.write(1, col, meta["birth_t"] if meta["birth_t"] is not None else "")
            worksheet.write(2, col, meta["g1_t"] if meta["g1_t"] is not None else "")

        # Column A, from row 4 downward: timepoints
        for i, t in enumerate(time_indices):
            row = 3 + i
            worksheet.write(row, 0, t)

    print("\n✅ Done.")
    print(f"Wrote Excel: {out_xlsx}")
    print(f"Total Z-slices analyzed across all cells and timepoints: {total_slices_analyzed}")


if __name__ == "__main__":
    main()
