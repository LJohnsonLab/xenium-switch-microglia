# Project: microglia_switch

> **Before any non-trivial work, read the last 3–4 entries of `JOURNAL.md` for recent decisions and unresolved items.**

Xenium in situ spatial transcriptomics of mouse brain under a microglia-specific APOE E4→E2 switch on a 5XFAD amyloid background.

## Pipeline

No analysis code yet. Expected flow once notebooks land:

1. Per-slide ingest of Xenium output (`cells.parquet`, `transcripts.parquet`, `cell_feature_matrix.h5`, `morphology.ome.tif`).
2. Apply per-sample ROIs from `spatial_data/*.geojson` to split each slide's `Region_1` into the seven samples (4s2_F1–F3, 4s2M_F1–F4).
3. Build per-sample Seurat / SpatialExperiment objects, QC, cluster, annotate.
4. Compare microglia and neighbouring cell states between `4s2` (E4 control) and `4s2M` (E2-switched microglia) littermates.

## Data

- `xenium_rawData/` (gitignored except `README.md`): two Xenium runs, both `Region_1` exports from instrument `XETG00118`, panel `mBrain_480g` (480 targets).
  - Slide `0022474` (run 2025-04-18, `4s2_FAD_LG`): samples 4s2_F1, 4s2M_F1, 4s2M_F2. Manifest is `experiment.xenium_slide1.xenium` (non-standard name).
  - Slide `0069080` (run 2026-04-17, `AgeXMetabolism_Run7`): samples 4s2_F2, 4s2_F3, 4s2M_F3, 4s2M_F4.
- `spatial_data/` (tracked): per-slide ROI exports from Xenium Explorer v4.1.1, CSV (point lists) + GeoJSON (polygons), named `20260512-slide_<id>.{csv,geojson}`. These are the authoritative sample masks.
- Source key PDF `4s2MxFAD Xenium_brain key_042126.pdf` lives inside `xenium_rawData/` and is gitignored; layouts on page 2 are needed for any manual ROI redo.

## Environment

R + tidyverse expected (per global preferences); no `renv`, `DESCRIPTION`, `_quarto.yml`, or `requirements.txt` yet. Add them when the first notebook is written.

## Conventions

- Sample IDs follow `xenium_rawData/README.md`: `4s2_Fn` = 5XFAD control (microglia E4), `4s2M_Fn` = 5XFAD tamoxifen-switched (microglia E2). Numbering is continuous across slides within a genotype.
- All mice are 5XFAD. The `4s2` vs `4s2M` distinction is about microglial APOE state, not amyloid background.
- Date-prefixed filenames use `YYYYMMDD-`.
- `xenium_rawData/<slide_dir>/` paths must be referenced verbatim; the two slide directories use different naming patterns.

## Key Decisions

- Raw Xenium output is too large for git; the README inside `xenium_rawData/` is the single tracked artefact in that folder, enforced by an allow-list in `.gitignore`.
- Sample masks are stored as ROI exports (CSV + GeoJSON) committed to `spatial_data/`, not regenerated from the PDF on each run.

## Current State

Day-zero. Raw acquisitions and sample-assignment ROIs are in place. No ingest code, no Seurat / SpatialExperiment objects, no QC, no DE yet. Sample `4s2M_F2` is partial — only a fragment remains on the slide; treat with care in any downstream QC.

## Gotchas

- Slide `0022474` stores its experiment manifest as `experiment.xenium_slide1.xenium`, **not** `experiment.xenium`; loaders that hard-code the standard name will fail on this slide.
- Both slides also carry sections from unrelated experiments inside `Region_1`. The ROIs in `spatial_data/` are the only correct way to isolate this project's seven samples.
- `4s2M_F2` is incomplete — flag it before any cell-count comparison or sample-level QC threshold.
- `.gitignore` is currently untracked (not yet committed by user choice).
