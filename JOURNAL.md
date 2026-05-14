# Journal

## 2026-05-12 — Bootstrap: Initial Project Snapshot

**Context:** First `/snapshot` invocation. No prior `CLAUDE.md` existed. Repo was initialised this session with two commits: the spatial_data ROI exports (`f886508`) and the xenium_rawData README (`bece3a1`).

**Baseline:** Data-analysis project for Xenium spatial transcriptomics on mouse brain, 5XFAD background, microglia-specific APOE4→APOE2 switch via `Tmem119-CreERT2`. Two slides (`0022474`, `0069080`) hold seven project samples (4s2_F1–F3 controls, 4s2M_F1–F4 switched) intermixed with unrelated sections. Sample masks live as Xenium Explorer v4.1.1 ROI exports (CSV+GeoJSON) in `spatial_data/`. Raw Xenium acquisitions are gitignored except for `xenium_rawData/README.md`, which encodes the slide↔run↔sample key. A `.gitignore` covering OS/editor/R/Python/Quarto artefacts plus an allow-list for `xenium_rawData/README.md` was authored but not yet committed (user has not authorised it). No analysis code, `renv`, or Quarto config exists yet.

**Ambiguities:**
- Target framework for ingest: Seurat v5 (per global preferences) or SpatialExperiment / Bioconductor?
- Should `4s2M_F2` (partial section) be excluded outright or carried through with a flag?
- Will sample-level objects be built per slide and merged, or all seven samples loaded together from raw via the ROI polygons?
- Is the `.gitignore` ready to commit, or are further entries pending?

**Token Budget:** ~210 tokens.

## 2026-05-12 — Ingest, QC, and three-way normalization scaffold

**What we did.** We built the ingest script `R/20260512_xenium2seurat.R`: `LoadXenium` per slide, attach per-sample ROIs from `spatial_data/<sample>_cells_stats.csv`, summarise transcripts via duckplyr, filter to `nCount_Xenium > 0 & median_qv >= 20`, merge with `JoinLayers`, save `data/20260512_microglia_switch_original_seurat.qs2` (480 × ~430k cells). We wrote `20260512_QC.qmd` (sample-vs-sample and slide-vs-slide diagnostics) and committed three parallel clustering notebooks: SCTransform without regressor, SCTransform with `vars.to.regress = "slide"`, and per-sample SpaNorm (`20260512_SpaNorm.qmd`). All three end with ScaleData → PCA(30) → FindNeighbors → Louvain(0.5) → UMAP. We resolved the SpaNorm 1.4.0 SpatialExperiment-only dispatch by wrapping each per-sample subset in a `SpatialExperiment` before the call, and bringing logcounts back as a v5 `CreateAssay5Object(data = ...)`. We cached the SpaNorm background read-through at `notebookLM/2026-05-12_spanorm-background.md` and cite (Salim, Genome Biology, 2025).

**Slide effect.** Slide2 vs Slide1 shows +25 % transcripts/cell, -7 % assignment rate, +33 % control-probe rate. Verdict (in `20260512_QC.html`) attributes this to acquisition/tissue/firmware drift, not XOA algorithmic change.

**Unresolved.**
- Run the SpaNorm build chunk end-to-end and save `data/20260512_microglia_switch_spanorm.qs2`; the `qs_save` line is commented in the qmd.
- Three-way visual comparison of SCT-no-regress vs SCT-slide-regress vs SpaNorm on PCA/UMAP/Louvain.
- Decide whether to drop `4s2M_F2` (partial section) before DE.

**Token Budget:** ~245 tokens.

## 2026-05-14 — SpaNorm sample.p = 0.25 promoted; marker analysis and cell-type calls landed

**Context & changes.** Two-day session. Built `20260513_SpaNorm_samplep025.qmd` and promoted `sample.p = 0.25` (package default) over `0.05` because the per-slide UMAPs balance more cleanly with identical cluster count; the 0.05 version stays as `20260512_SpaNorm_sample005.qmd`. Imported `R/20260514-helper_functions.R` from the co-registration project (renamed from `20260217-helper_functions.R`) and rewrote `fiddle()`: clusters on x-axis (color-matched, rotated 45°), genes stacked as left-strip facets, expression-level tick numbers hidden. Added a marker-analysis section to both `20260512_SCT_no_regression.qmd` and `20260513_SpaNorm_samplep025.qmd`: presto-backed `FindAllMarkers(only.pos = TRUE)` (gated `eval: false` after one-pass render), kbl + scroll_box top-10 table, and a `fiddle()` panel on top-5. Markers cached to `results/20260514-find_all_markers_{sct,spanorm_p025}_*.csv`. Annotated 30 SCT clusters and 25 SpaNorm clusters by top-10 markers and organized them by cell-type similarity in `results/20260514-cluster_celltype_assignments_{sct,spanorm}.md`.

**Key biological findings.** SpaNorm reveals a discrete cytotoxic-T-cell cluster (24: Cd3d/Gzmb/Nkg7/Itk) that SCT did not separate. SpaNorm loses the SCT IFN-reactive-astrocyte cluster (SCT 28); signal likely moved into the T cells (Ccl5/Tnf) or the reactive-mixed clusters 13/15. SpaNorm removes the Hbb-bs+ RBC contamination cluster (SCT 27) and the Calb2-contaminated OL cluster (SCT 26).

**Unresolved.** Spatially inspect the SpaNorm T-cell cluster (24) near plaques; sub-cluster microglia (SpaNorm 4 / SCT 3) before DE; run DoubletFinder on SpaNorm 13/15/23; sub-cluster SpaNorm astrocytes (5, 8) to recover the SCT cluster 28 DAA signature; start the E4-vs-E2 DE between `4s2` and `4s2M`.

**Token Budget:** ~280 tokens.
