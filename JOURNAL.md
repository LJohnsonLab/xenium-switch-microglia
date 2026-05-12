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
