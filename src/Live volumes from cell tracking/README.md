# Tracked cell volume analysis (live imaging)

## Overview
This pipeline reconstructs cell volume dynamics across time using manually tracked cells and segmentation-derived measurements.

Segmentation is performed upstream using **Cellpose**, and shape features are extracted using **QuantifyPolarity**. This script links those outputs across time using manually curated cell identities.

---

## Workflow summary

1. Segment each timepoint (Cellpose)
2. Extract cell measurements (QuantifyPolarity)
3. Compile all timepoints into:
   - `All_Timepoints_Combined.xlsx`
4. Manually track cells using their IDs
5. Run this script to reconstruct:
   - time-resolved cell shape
   - volume across time
   - lineage relationships

---

## Tracking format

Tracking is defined in:


tracked_IDs.txt


Each line represents a tracked cell or lineage.

### Format


START=<timepoint>: ID.ID.ID | D1_ID.D1_ID | D2_ID.D2_ID


### Rules

- `.` → next timepoint  
- `|` → cell division branches  
- first segment → parent  
- second → daughter 1 (D1)  
- third → daughter 2 (D2)  
- `-` → missing/untracked  

### Example


START=1: 10.15.18.22. | 50.60.70 | 80.90.100


---

## Data structure

```
Data/
   Live volumes from cell tracking/
      Tracked ex cell 1/
         Apical Z/
            tracked_IDs.txt
            Denoised_Image Seq/
               All_Timepoints_Combined.xlsx
         Middle Z/
         Basal Z/
```

Each plane must contain:
- tracking file
- compiled measurements

---

## Volume calculation

Volume is estimated using three planes:

- basal
- middle
- apical

Areas are converted to µm² and combined using a frustum model.

Output:

Volume (µm³)


---

## Example dataset

The repository includes a reduced dataset:

- Only a subset of timepoints is included in plane folders
- Full tracking and measurements are still present in:
  - `tracked_IDs.txt`
  - `All_Timepoints_Combined.xlsx`

Outputs will still include all timepoints.

---

## Outputs

Generated in:


results/Live volumes from cell tracking/


Includes:
- tracked cell tables
- BVMA measurement tables
- volume calculations
- 3D visualization images

---

## Installation


pip install -r requirements.txt


---

## Run

From repo root:


python "src\Live volumes from cell tracking\Live volumes from cell tracking script.py"


---

## Notes

- Tracking is manual and must be correct
- Script does not perform segmentation or tracking
- Handles:
  - multiple datasets
  - multiple lineages
  - cell divisions
- Visualization is simplified and intended for reference only
