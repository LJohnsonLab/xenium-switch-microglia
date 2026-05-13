# Project: microglia_switch

> **Before any non-trivial work, read the last 3–4 entries of `JOURNAL.md` for recent decisions and unresolved items.**

> **Before producing any prose that explains results (qmd verdicts, JOURNAL entries, chat replies that interpret findings, commit-message bodies that summarise outcomes), read `context/writingStyle_OSmithies.json` and apply that style profile. The intended voice is direct, active, first-person, concise, and quantitative; the JSON is the authoritative spec.**

Xenium in situ spatial transcriptomics of mouse brain under a microglia-specific APOE E4→E2 switch on a 5XFAD amyloid background.

## Pipeline

1. **Ingest** (`R/20260512_xenium2seurat.R`): `LoadXenium` per slide, attach per-sample ROIs from `spatial_data/<sample>_cells_stats.csv`, summarise transcripts via duckplyr, filter to `!is.na(sample_id) & nCount_Xenium > 0 & median_qv >= 20`, merge with `add.cell.ids = c("slide1", "slide2")`, `JoinLayers`, save as `data/20260512_microglia_switch_original_seurat.qs2`.
2. **QC** (`20260512_QC.qmd`): sample-vs-sample and slide-vs-slide diagnostics on per-cell QC and slide `metrics_summary.csv`; verdict in the rendered HTML.
3. **Normalization comparison**: three parallel notebooks at the project root.
   - `20260512_SCT_no_regression.qmd`: SCTransform with no covariate.
   - `20260512_SCT_regressed_by_slide.qmd`: SCTransform with `vars.to.regress = "slide"`.
   - `20260512_SpaNorm.qmd`: per-sample SpaNorm (seven fits), logcounts written back into a v5 `SpaNorm` assay.
   All three end with PCA(30) → FindNeighbors → Louvain at resolution 0.5 → UMAP, plotted by slide and by sample.
4. **DE / spatial comparison** between `4s2` (E4) and `4s2M` (E2) littermates: not started.

## Data

- `xenium_rawData/` (gitignored except `README.md` and the brain-key PDF): two Xenium runs, both `Region_1` exports from instrument `XETG00118`, panel `mBrain_480g` (480 targets).
  - Slide `0022474` (run 2025-04-18, `4s2_FAD_LG`): samples 4s2_F1, 4s2M_F1, 4s2M_F2. Manifest was renamed in-place from `experiment.xenium_slide1.xenium` to `experiment.xenium` on 2026-05-12 (see README).
  - Slide `0069080` (run 2026-04-17, `AgeXMetabolism_Run7`): samples 4s2_F2, 4s2_F3, 4s2M_F3, 4s2M_F4.
- `spatial_data/` (tracked): Xenium Explorer v4.1.1 ROI exports per slide (CSV + GeoJSON, `20260512-slide_<id>.{csv,geojson}`) plus seven `<sample>_cells_stats.csv` per-sample cell tables. These are the authoritative sample masks.
- `data/` (gitignored except `README.md`): processed artefacts. Current files: `20260512_microglia_switch_original_seurat.qs2` (480 × ~430k cells), `20260512_microglia_switch_sct.qs2`, `20260512_microglia_switch_sct_regressed_by_slide.qs2`. The SpaNorm checkpoint (`data/20260512_microglia_switch_spanorm.qs2`) is built by the SpaNorm qmd's `eval: false` chunk on first run.
- Source key PDF `4s2MxFAD Xenium_brain key_042126.pdf` is tracked under `xenium_rawData/`; page 2 holds the sample-layout diagrams needed for any manual ROI redo.

## Environment

R 4.5.2 + tidyverse, native pipe `|>` only. Key packages: Seurat v5, qs2, fs, duckplyr, arrow (>= 24 for the 5.5 GB `transcripts.parquet`; arrow 22 throws "Invalid number of indices: 0"), SpatialExperiment, SpaNorm 1.4.0 (Bioc 3.22, dispatches only on `SpatialExperiment`), scCustomize. No `renv`, `DESCRIPTION`, or `_quarto.yml` yet.

## Conventions

- Sample IDs follow `xenium_rawData/README.md`: `4s2_Fn` = 5XFAD control (microglia E4), `4s2M_Fn` = 5XFAD tamoxifen-switched (microglia E2). Numbering is continuous across slides within a genotype.
- All mice are 5XFAD. The `4s2` vs `4s2M` distinction is about microglial APOE state, not amyloid background.
- Date-prefixed filenames use `YYYYMMDD-`.
- `xenium_rawData/<slide_dir>/` paths must be referenced verbatim; the two slide directories use different naming patterns.

## Key Decisions

- Raw Xenium output is too large for git; only `xenium_rawData/README.md` and the brain-key PDF are tracked, via an allow-list in `.gitignore`. The same allow-list pattern covers `data/` (README only) and ignores `context/` and `notebookLM/` (local-only assets).
- Sample masks are stored as ROI exports (CSV + GeoJSON) committed to `spatial_data/`, not regenerated from the PDF on each run.
- SpaNorm fits run **per sample** (seven fits), not per slide. Spatial smoothness collapses across disconnected tissue sections; running per slide would feed the model multiple unrelated sections at once (Salim, Genome Biology, 2025).
- No Quarto / knitr caching anywhere. Heavy steps save intermediates under `data/*.qs2` and the qmd `eval: false` build chunks reload those checkpoints on render.
- Cached NotebookLM read-throughs live at `notebookLM/<date>_<topic>.md` (gitignored). They back the inline literature citations in the qmds.

## Current State

Ingest, QC, and three normalization notebooks are committed. The merged Seurat object (480 × ~430k cells, 7 samples, 2 slides) holds raw counts; SCT-no-regress and SCT-slide-regress checkpoints exist on disk; the SpaNorm checkpoint is built the first time the SpaNorm qmd is run with the build chunk set to `eval: true`. Next step is the three-way comparison (PCA / UMAP / Louvain) and the DE/spatial analysis between `4s2` (E4) and `4s2M` (E2) littermates. Sample `4s2M_F2` is a partial section; the QC verdict cleared it as non-outlier on per-cell QC, but treat it cautiously in cluster-composition comparisons.

## Gotchas

- Both slides carry sections from unrelated experiments inside `Region_1`. The ROIs in `spatial_data/` are the only correct way to isolate this project's seven samples.
- `4s2M_F2` is a partial section. Flag it in any cluster-composition or cell-count comparison.
- `LoadXenium` names each slide's FOV `fov`; after `merge` they become `fov` and `fov.2`, **not** the slide labels. Build any per-sample FOV lookup from `Images(xen)`, never from a hardcoded slide name.
- SpaNorm 1.4.0 dispatches only on `SpatialExperiment`. Each per-sample build wraps a `SpatialExperiment(counts, spatialCoords)` before calling `SpaNorm()`, then writes the returned `logcounts` back into Seurat as a v5 `SpaNorm` assay via `CreateAssay5Object(data = ...)`.
- Loading `SpatialExperiment` pulls `IRanges`, whose `reduce()` masks `purrr::reduce()`. Qualify the call inside any qmd that uses both packages.
- A custom v5 assay carrying only the `data` layer needs an explicit `ScaleData()` before `RunPCA()`, and `FindVariableFeatures(method = "vst")` fails on it because vst needs counts. Use `VariableFeatures(xen) <- rownames(xen)` for the 480-gene targeted panel.
- Slide effect on per-cell QC (Slide2 vs Slide1): +25 % transcripts/cell, -7 % assignment rate, +33 % control-probe rate. Verdict in `20260512_QC.html` attributes it to acquisition / tissue / instrument-firmware differences, not XOA algorithmic drift (v3.2.0.7 → v3.3.0.1 segmentation changes target boundary-stain only).
