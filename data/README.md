# `data/` — processed analysis artefacts

Everything under `data/` is gitignored except this README. Files here are produced by scripts in `R/` (and later `qmd/`) from the gitignored raw acquisitions under `xenium_rawData/`. They are reproducible from code + raw data; do not commit them.

Filename convention: `<YYYYMMDD>_<short-description>.<ext>`. Use the date the artefact was generated, not the date of the underlying acquisition.

## Current contents

| File | Produced by | Description |
|------|-------------|-------------|
| `20260512_microglia_switch_original_seurat.qs2` | `R/20260512_xenium2seurat.R` | Merged Seurat v5 object across both slides. 480 features (full `mBrain_480g` panel) x 429,980 cells, after the ingest-time filters (cell assigned to one of the seven sample ROIs, `nCount_Xenium > 0`, `median_qv >= 20`). Holds `Xenium` plus the per-feature-class assays (`BlankCodeword`, `ControlCodeword`, `ControlProbe`, etc.) and two FOV slots (one per slide, re-keyed `Xenium_` / `fov2_` by Seurat at merge time). Per-cell metadata: `sample_id`, `slide`, `condition` (control vs switched), `partial_section` (TRUE only for 4s2M_F2), `nCount_raw`, `nFeature_raw`, `median_qv`, `cell_area`. No normalisation, no clustering. Restore with `qs2::qs_read()`. |

## Conventions

- `_original_seurat`: the first Seurat object built straight from `LoadXenium()` + ROI assignment, with QC filters applied but no normalisation or dimensionality reduction. Downstream artefacts (normalised, clustered, annotated) get their own descriptive suffixes.
- Prefer `.qs2` over `.rds` for Seurat objects; `qs2::qs_save()` is much faster and produces smaller files.
