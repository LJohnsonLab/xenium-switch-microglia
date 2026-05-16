# Project: microglia_switch

> **Before any non-trivial work, read the last 3–4 entries of `JOURNAL.md` for recent decisions and unresolved items.**

> **Before producing any prose that explains results (qmd verdicts, JOURNAL entries, chat replies that interpret findings, commit-message bodies that summarise outcomes), read `context/writingStyle_OSmithies.json` and apply that style profile. The intended voice is direct, active, first-person, concise, and quantitative; the JSON is the authoritative spec.**

Xenium in situ spatial transcriptomics of mouse brain under a microglia-specific APOE E4→E2 switch on a 5XFAD amyloid background.

## Pipeline

1. **Ingest** (`R/20260512_xenium2seurat.R`): `LoadXenium` per slide, attach per-sample ROIs from `spatial_data/<sample>_cells_stats.csv`, summarise transcripts via duckplyr, filter to `!is.na(sample_id) & nCount_Xenium > 0 & median_qv >= 20`, merge with `add.cell.ids = c("slide1", "slide2")`, `JoinLayers`, save as `data/20260512_microglia_switch_original_seurat.qs2`.
2. **QC** (`20260512_QC.qmd`): sample-vs-sample and slide-vs-slide diagnostics on per-cell QC and slide `metrics_summary.csv`; verdict in the rendered HTML.
3. **Normalization comparison**: parallel notebooks at the project root.
   - `20260512_SCT_no_regression.qmd`: SCTransform, no covariate.
   - `20260512_SCT_regressed_by_slide.qmd`: SCTransform with `vars.to.regress = "slide"`.
   - `20260512_SpaNorm_sample005.qmd`: per-sample SpaNorm at `sample.p = 0.05`.
   - `20260513_SpaNorm_samplep025.qmd` (**project default**, promoted 2026-05-14): per-sample SpaNorm at `sample.p = 0.25`; cleaner slide balance than 0.05 with identical cluster count.
   All end with PCA(30) → FindNeighbors → Louvain(0.5) → UMAP plotted by slide and by sample.
4. **Marker analysis & cell-type annotation** (every clustering qmd): presto-backed `FindAllMarkers(only.pos = TRUE)` (chunk gated `eval: false` after a one-pass render) writes `results/<date>-find_all_markers_*.csv`; a top-10-per-cluster table (kbl + scroll_box) and a `fiddle()` violin panel render in each qmd. Cell-type calls and contamination flags live in `notes/<date>-cluster_celltype_assignments_*.md`.
5. **Clean SpaNorm re-fit** (`20260516_SpaNorm_samplep025_clean.qmd`): drops the artefactual cluster-20 cells from the parent SpaNorm build, re-fits SpaNorm per sample, brackets Louvain at 0.5 (22 clusters) and 0.7 (30 clusters, **project default**). Saves `data/20260516_microglia_switch_spanorm_p025_clean.qs2` with both `clusters_res05` / `clusters_res07` labels.
6. **Microglia sub-cluster** (`20260517_microglia_subcluster.qmd`): subset cluster 2 of the clean build (~36k microglia), `DietSeurat` to Xenium counts only, re-fit SpaNorm per sample on microglia, brackets Louvain at 0.5 (11) and 0.3 (7, **microglia default**) stashed in `clusters_resXX` columns. Saves `data/20260517_microglia_subcluster.qs2`. Cleaning recipe (wholesale drop of contamination clusters + per-cell `AddModuleScore()` against seven contamination-marker sets) in the assignment note.
7. **DE / spatial comparison** between `4s2` (E4) and `4s2M` (E2) littermates: not started.

## Data

- `xenium_rawData/` (gitignored except `README.md` and the brain-key PDF): two Xenium runs, both `Region_1` exports from instrument `XETG00118`, panel `mBrain_480g` (480 targets).
  - Slide `0022474` (run 2025-04-18, `4s2_FAD_LG`): samples 4s2_F1, 4s2M_F1, 4s2M_F2. Manifest was renamed in-place from `experiment.xenium_slide1.xenium` to `experiment.xenium` on 2026-05-12 (see README).
  - Slide `0069080` (run 2026-04-17, `AgeXMetabolism_Run7`): samples 4s2_F2, 4s2_F3, 4s2M_F3, 4s2M_F4.
- `spatial_data/` (tracked): Xenium Explorer v4.1.1 ROI exports per slide (CSV + GeoJSON, `20260512-slide_<id>.{csv,geojson}`) plus seven `<sample>_cells_stats.csv` per-sample cell tables. These are the authoritative sample masks.
- `data/` (gitignored except `README.md`): processed artefacts. Current files: `20260512_microglia_switch_original_seurat.qs2` (480 × ~430k cells), `20260512_microglia_switch_sct.qs2`, `20260512_microglia_switch_sct_regressed_by_slide.qs2`, `20260512_microglia_switch_spanorm_p025.qs2` (parent SpaNorm), `20260516_microglia_switch_spanorm_p025_clean.qs2` (clean SpaNorm re-fit, project default; both `clusters_res05` 22-cluster and `clusters_res07` 30-cluster labels), `20260517_microglia_subcluster.qs2` (microglia-only SpaNorm re-fit, ~36k cells; both `clusters_res05` 11-cluster and `clusters_res03` 7-cluster labels). The 0.05 variant checkpoint is built on demand by `20260512_SpaNorm_sample005.qmd`.
- Source key PDF `4s2MxFAD Xenium_brain key_042126.pdf` is tracked under `xenium_rawData/`; page 2 holds the sample-layout diagrams needed for any manual ROI redo.

## Environment

R 4.5.2 + tidyverse, native pipe `|>` only. Key packages: Seurat v5, qs2, fs, duckplyr, arrow (>= 24 for the 5.5 GB `transcripts.parquet`; arrow 22 throws "Invalid number of indices: 0"), SpatialExperiment, SpaNorm 1.4.0 (Bioc 3.22, dispatches only on `SpatialExperiment`), scCustomize. No `renv`, `DESCRIPTION`, or `_quarto.yml` yet.

## Conventions

- Sample IDs follow `xenium_rawData/README.md`: `4s2_Fn` = 5XFAD control (microglia E4), `4s2M_Fn` = 5XFAD tamoxifen-switched (microglia E2). Numbering is continuous across slides within a genotype.
- All mice are 5XFAD. The `4s2` vs `4s2M` distinction is about microglial APOE state, not amyloid background.
- Date-prefixed filenames use `YYYYMMDD-`.
- `xenium_rawData/<slide_dir>/` paths must be referenced verbatim; the two slide directories use different naming patterns.
- Shared plotting helpers live at `R/20260514-helper_functions.R`. `fiddle()` draws per-cluster violins with clusters on the x-axis (color-matched, rotated 45°) and genes stacked as left-strip facets, with expression-level tick numbers hidden.
- When running `FindClusters` at multiple resolutions on the same object, stash each result in a named metadata column (`obj$clusters_resXX <- obj$seurat_clusters`). `seurat_clusters` is overwritten by every call; named columns let downstream plots / `Idents()` / `FindAllMarkers` switch resolutions without re-running clustering.

## Key Decisions

- Raw Xenium output is too large for git; only `xenium_rawData/README.md` and the brain-key PDF are tracked, via an allow-list in `.gitignore`. The same allow-list pattern covers `data/` (README only) and ignores `context/` and `notebookLM/` (local-only assets).
- Sample masks are stored as ROI exports (CSV + GeoJSON) committed to `spatial_data/`, not regenerated from the PDF on each run.
- SpaNorm fits run **per sample** (seven fits), not per slide. Spatial smoothness collapses across disconnected tissue sections; running per slide would feed the model multiple unrelated sections at once (Salim, Genome Biology, 2025).
- No Quarto / knitr caching anywhere. Heavy steps save intermediates under `data/*.qs2` and the qmd `eval: false` build chunks reload those checkpoints on render.
- Cached NotebookLM read-throughs live at `notebookLM/<date>_<topic>.md` (gitignored). They back the inline literature citations in the qmds.
- Promoted SpaNorm `sample.p = 0.25` (package default) over the 0.05 value initially suggested for ~400k-cell datasets (Salim, Genome Biology, 2025). At 0.25 the per-slide UMAPs balance more cleanly with identical 25-cluster count; decision recorded in `20260513_SpaNorm_samplep025.qmd`.
- Cell-type calls come from the top-10 FindAllMarkers per cluster ranked by `avg_log2FC`, organized by cell-type similarity (not cluster number) in the assignment md, with explicit contamination flags for off-profile genes.
- Promoted Louvain 0.7 on the clean SpaNorm re-fit as the project working scaffold. 0.7 resolves homeostatic vs DAM microglia (clusters 12 / 6), recovers the discrete T-cell cluster (28: Cd3d / Gzmb / Nkg7 / Itk) that the 0.5 cut lost, and re-splits Sst+ / Vip+ / Pvalb+ GABAergic subclasses.
- Promoted Louvain 0.3 for the microglia sub-cluster (`20260517_microglia_subcluster.qmd`). At 0.5 the OL-contamination and multi-lineage doublet sub-clusters carry a visible slide imbalance; at 0.3 they collapse into broader categories that balance across slides. Three microglia states survive at 0.3: DAM (0), homeostatic (1), IRM + proliferating (6: Cxcl10 / Mki67 / Ifit1 / Ifit2).
- Doublet / segmentation-contamination cleaning uses `AddModuleScore()` against seven contamination-marker sets (OL, excitatory, inhibitory, astrocyte, endothelial, pericyte, ependymal/CP), not DoubletFinder; the 480-probe panel breaks DoubletFinder's synthetic-doublet kNN assumption.
- Custom APOE probes are on the mBrain_480g panel: `rs429358-WT-228` (E3/E4 allele) and `rs7412-ALT-226:T` (E2 allele). They enrich per microglia state at 0.5 (rs429358-WT in DAM Stage-2 cluster 4; rs7412-ALT in low-FC homeostatic cluster 6), giving a per-cell molecular hook for the E4-vs-E2 contrast.
- Cell-type assignment notes live under `notes/`; marker CSVs under `results/`; convergence sidecars under `notes/`.

## Current State

Ingest, QC, four normalization notebooks, the clean SpaNorm re-fit, and a microglia sub-cluster build are committed. Clean SpaNorm at Louvain 0.7 (30 clusters, project default) resolves the homeostatic vs DAM microglia split (12 / 6), recovers the T-cell cluster (28: Cd3d / Gzmb / Nkg7 / Itk) that the parent build at 0.5 had revealed but the 0.5 clean cut lost, and re-splits Sst+ / Vip+ / Pvalb+ GABAergic. The microglia sub-cluster at Louvain 0.3 (7 sub-clusters; microglia working scaffold) carries three robust microglia states: DAM (0), homeostatic (1), IRM + proliferating (6: Cxcl10 / Mki67 / Ifit1 / Ifit2). Custom APOE probes (rs429358-WT for E3/E4, rs7412-ALT for E2) enrich per microglia state and give a per-cell molecular hook for the E4-vs-E2 contrast. Next steps: pseudobulk E4-vs-E2 DE on cleaned microglia clusters 0 / 1 / 6 of the 0.3 sub-cluster (after wholesale drop of contamination clusters + per-cell `AddModuleScore()` filter); sub-cluster the consolidated DAM at higher Louvain to recover Stage-1 vs Stage-2 (visible at 0.5 as cluster 4 antigen-presenting Spp1 / Gpnmb / H2-Eb1 / Cd74 vs 0 lipid-chemokine Lyz2 / Cst7 / Ccl3 / Lpl); spatially inspect the IRM cluster near amyloid plaques; sub-cluster the SpaNorm astrocyte pair (clean-0.7 clusters 2, 8) for an IFN-responsive DAA signature; cross-validate APOE probe per-cell calls against the genotype labels. Sample `4s2M_F2` is a partial section; flag it in cluster-composition comparisons.

## Gotchas

- Both slides carry sections from unrelated experiments inside `Region_1`. The ROIs in `spatial_data/` are the only correct way to isolate this project's seven samples.
- `4s2M_F2` is a partial section. Flag it in any cluster-composition or cell-count comparison.
- `LoadXenium` names each slide's FOV `fov`; after `merge` they become `fov` and `fov.2`, **not** the slide labels. Build any per-sample FOV lookup from `Images(xen)`, never from a hardcoded slide name.
- SpaNorm 1.4.0 dispatches only on `SpatialExperiment`. Each per-sample build wraps a `SpatialExperiment(counts, spatialCoords)` before calling `SpaNorm()`, then writes the returned `logcounts` back into Seurat as a v5 `SpaNorm` assay via `CreateAssay5Object(data = ...)`.
- Loading `SpatialExperiment` pulls `IRanges`, whose `reduce()` masks `purrr::reduce()`. Qualify the call inside any qmd that uses both packages.
- A custom v5 assay carrying only the `data` layer needs an explicit `ScaleData()` before `RunPCA()`, and `FindVariableFeatures(method = "vst")` fails on it because vst needs counts. Use `VariableFeatures(xen) <- rownames(xen)` for the 480-gene targeted panel.
- Slide effect on per-cell QC (Slide2 vs Slide1): +25 % transcripts/cell, -7 % assignment rate, +33 % control-probe rate. Verdict in `20260512_QC.html` attributes it to acquisition / tissue / instrument-firmware differences, not XOA algorithmic drift (v3.2.0.7 → v3.3.0.1 segmentation changes target boundary-stain only).
- `DietSeurat` refuses to drop the currently-default assay. Switch the default to the assay you want to keep first, then call `DietSeurat(..., assays = ...)`.
- Cluster-id filters operate on factor levels as strings; `seurat_clusters %in% c("20", "25")` silently no-ops on levels that do not exist. Confirm the levels are present (`levels(obj$seurat_clusters)`) before relying on the filter.
- The literal slide → FOV map after `merge` is `c(Slide1_0022474 = "fov", Slide2_0069080 = "fov.2")`. Use it directly; do not infer via `Images()`-based introspection.
