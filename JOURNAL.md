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
