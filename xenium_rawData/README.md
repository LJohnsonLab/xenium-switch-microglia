# Xenium raw data — sample key and file provenance

Source key: `4s2MxFAD Xenium_brain key_042126.pdf` (dated 2026-04-21).

Mouse brain, fresh-frozen, custom panel `mBrain_480g` (panel design `JRZFVH`, 480 custom targets, 0 predesigned), instrument `XETG00118`.

## Biological model

All mice carry the **5XFAD** transgene (amyloid pathology background). On top of 5XFAD, mice are on an inducible APOE "switch": `APOE4s2flox/flox × Tmem119-CreERT2` (humanized APOE locus). Tamoxifen flips the floxed allele in Tmem119+ microglia from E4 to E2 while leaving E4 expression intact in all other CNS and peripheral cells.

- **`4s2M` / `MC/+`** — 5XFAD; tamoxifen-induced; microglia express E2, rest of body expresses E4. RNA-seq sample prefix `M_*` (e.g. `M_9`, `M_10`).
- **`4s2` / `+/+`** — 5XFAD; Cre-negative or vehicle littermate controls (microglia and body both express E4).

## Subdirectory ↔ slide ↔ run

Acquisition date is extracted from `run_start_time` in each `experiment.xenium` manifest. The two slides were processed on different Xenium software releases; the version gap drives the per-cell QC slide effect documented in `20260512_QC.html` (see Notes).

| Local subdirectory | Slide ID | Run name | Project / run label | Acquisition date | Analysis SW | Instrument SW | Cells (manifest) |
|------------|------------|------------|------------|------------|------------|------------|------------|
| `Slide1_output-XETG00118__0022474__Region_1__20250418__203103/` | `0022474` | `4s2_FAD_LG` (cassette `LG_Slide1`) | APOE 4s2 × 5XFAD — Run 1 | 2025-04-18 | `xenium-3.2.0.7` | `3.2.1.2` | 306,641 |
| `output-XETG00118__0069080__Region_1__20260417__210507/` | `0069080` | `20260417_AgeXMetabolism_Run7` (cassette `20260417_slide2_0069080`) | Aging_Metab — Run 7 | 2026-04-17 | `xenium-3.3.0.1` | `3.4.1.0` | 415,342 |

The two slides export `Region_1` only: a single segmentation spanning the whole capture area. Each slide also carries sections from unrelated experiments inside `Region_1`; this project's seven samples are isolated by ROI masks, not by any slide-level split (see Ingest).

## Samples per slide

Samples relevant to this project (boxed in red on the PDF key).

### Slide `0022474` — APOE 4s2 × 5XFAD Run 1

| Sample ID | Genotype                 | Microglia APOE |
|-----------|--------------------------|----------------|
| `4s2_F1`  | 5XFAD; `+/+` (control)   | E4             |
| `4s2M_F1` | 5XFAD; `MC/+` (switched) | E2             |
| `4s2M_F2` | 5XFAD; `MC/+` (switched) | E2             |

### Slide `0069080` — Aging_Metab Run 7

| Sample ID | Genotype                 | Microglia APOE |
|-----------|--------------------------|----------------|
| `4s2_F2`  | 5XFAD; `+/+` (control)   | E4             |
| `4s2_F3`  | 5XFAD; `+/+` (control)   | E4             |
| `4s2M_F3` | 5XFAD; `MC/+` (switched) | E2             |
| `4s2M_F4` | 5XFAD; `MC/+` (switched) | E2             |

## Sample-ID convention

All samples are 5XFAD; the prefix encodes the microglial APOE state, not amyloid background.

- `4s2_Fn` — 5XFAD; `+/+` (Cre-neg or vehicle), microglia E4. Female `n`.
- `4s2M_Fn` — 5XFAD; `MC/+` (tamoxifen-switched), microglia E2. Female `n`.
- Numbering (`F1`, `F2`, …) is continuous across slides within a genotype.

## Custom APOE allele probes

The `mBrain_480g` panel is fully custom and includes four rs-prefixed APOE allele probes alongside the bulk `Apoe` probe. Probe names as they appear in `gene_panel.json` and in the Seurat feature set:

| Probe | Codon | Allele detected | Use in this project |
|-------|-------|-----------------|---------------------|
| `Apoe` | — | total APOE transcript | bulk APOE expression |
| `rs429358_WT_228` | 112 | T (Cys112) | E2-evidence; this model has no E3, so T uniquely marks E2 |
| `rs7412_ALT_226:T` | 158 | T (Cys158) | E2-evidence; T defines E2 |
| `rs429358_ALT_228:C` | 112 | C (Arg112) | E4 codon-112; insensitive in this assay (mean ≈ 0.21 in 4s2), not used |
| `rs7412_WT_226` | 158 | C (Arg158) | E4/E3 codon-158; does not separate arms, not used |

The per-cell switch+ / switch- call rests on the **two T-allele probes** (`rs429358_WT_228 + rs7412_alt`); see `20260518_microglia_APOE_per_cell_call.qmd` and CLAUDE.md "Custom APOE probes." `allelic_discrimination_slide_0069080.png` (top level, gitignored) is a QC scatter of the two allele axes for slide `0069080`.

## Output-bundle contents

Both slides ship the standard Xenium Onboard Analysis (XOA) bundle; layout is identical between slides. Sizes below are for slide `0022474`; slide `0069080` scales with its larger cell count. Every file in these subdirectories is gitignored (the allow-list in `.gitignore` tracks only this README and the brain-key PDF).

| Path | ~Size | Contents | Used downstream? |
|------|-------|----------|------------------|
| `experiment.xenium` | 4 KB | JSON run manifest (run name, panel, software versions, cell count, paths to the `.zarr` Explorer files). `LoadXenium()` reads this to resolve the bundle. | **Yes** — manifest read by `LoadXenium()`. |
| `cell_feature_matrix.h5` | 53 MB | Cell × feature count matrix (HDF5). | **Yes** — the count assay loaded by `LoadXenium()`. |
| `cell_feature_matrix/` | 248 MB | Same matrix in MEX triplet (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`). Redundant with the `.h5`. | No (`.h5` is used). |
| `cells.parquet` / `cells.csv.gz` | 6 / 18 MB | Per-cell table: centroid x/y, area, transcript counts, nucleus stats. | **Yes** — centroids via `LoadXenium()`. |
| `cell_boundaries.parquet` / `.csv.gz` | 27 / 76 MB | Cell-segmentation polygon vertices. | **Yes** — segmentation FOV via `LoadXenium()`. |
| `nucleus_boundaries.parquet` / `.csv.gz` | 24 / 71 MB | Nucleus-segmentation polygon vertices. | Loaded by `LoadXenium()`; not used in current analysis. |
| `transcripts.parquet` | 5.1 GB | Every decoded transcript: feature, x/y/z, QV, assigned `cell_id`. Requires `arrow >= 24` (arrow 22 throws "Invalid number of indices: 0"). | **Yes** — read separately via `duckplyr` for per-cell QC (`nCount_raw`, `nFeature_raw`, `median_qv`). |
| `gene_panel.json` | 240 KB | Panel definition: the 480 target names (incl. the custom APOE allele probes) and probe metadata. | Reference (probe-name source of truth). |
| `metrics_summary.csv` | 4 KB | One-row slide-level QC: decode rate, transcripts/cell, control-probe rate, segmentation fractions. | **Yes** — slide-vs-slide QC in `20260512_QC.qmd`. |
| `analysis/` | 115 MB | XOA's own clustering / PCA / UMAP / diffexp on the whole slide (graph-based, default params). | No — we re-cluster on ROI-masked, re-normalized data. |
| `analysis_summary.html` | 12 MB | XOA's interactive QC + analysis report. | Reference only. |
| `*.zarr.zip` (`transcripts`, `cells`, `cell_feature_matrix`, `analysis`) | 0.1–5.2 GB | Xenium Explorer viewer files (where the `spatial_data/` ROIs were drawn). | Indirect — ROIs in `spatial_data/` were exported from Explorer reading these. |
| `morphology.ome.tif` | 22 GB | Full-resolution multi-channel morphology image (DAPI + boundary stains). | No (not loaded). |
| `morphology_focus/morphology_focus_0000.ome.tif` | 2.1 GB | Focus-stacked morphology image, single OME-TIFF. | Optional `LoadXenium()` image; not loaded in current ingest. |
| `aux_outputs/` | 438 MB | `overview_scan.png`, FOV-location JSONs, `per_cycle_channel_images/` (60 per-cycle fluorescence TIFFs). | Reference / instrument QC only. |

## Ingest

`R/20260512_xenium2seurat.R` builds `data/20260512_microglia_switch_original_seurat.qs2` (480 × ~430k cells) from these bundles. Transformations applied, in order:

1. **Load per slide.** `LoadXenium(slide_dir)` reads the manifest, `.h5` matrix, `cells.parquet` centroids, and cell/nucleus boundaries. The FOV is named `fov` on each slide; after merge they become `fov` / `fov.2`, never the slide labels.
2. **Per-cell transcript QC.** `transcripts.parquet` is summarized with `duckplyr` (`read_parquet_duckdb`) into `nCount_raw`, `nFeature_raw`, and `median_qv` per `cell_id`.
3. **Sample assignment by ROI.** Each cell's `sample_id` comes from the per-sample masks in `spatial_data/<sample>_cells_stats.csv` (Xenium Explorer v4.1.1 ROI exports), joined on `cell_id`. These ROIs are the only correct way to isolate the seven project samples from the unrelated sections in `Region_1`. `cell_area` is carried over from the ROI table.
4. **Filter.** `subset(!is.na(sample_id) & nCount_Xenium > 0 & median_qv >= 20)` — keeps only ROI-assigned, non-empty, Q20+ cells.
5. **Merge.** `merge(slide1, slide2, add.cell.ids = c("slide1","slide2"))`, then attach `condition` (control / switched) from the sample lookup and flag `partial_section = (sample_id == "4s2M_F2")`. `JoinLayers()` collapses the per-slide count layers.

The slide → FOV map after merge is `c(Slide1_0022474 = "fov", Slide2_0069080 = "fov.2")`.

## Original shared-drive paths (for provenance)

- **APOE 4s2 × 5XFAD Run 1** `07. Xenium Spatial Transcriptomics / 1 APOE4s2 x 5XFAD Xenium 2025-04-21 / 20250418_203052_4s2_FAD_LG-XeniumData / Slide1_output-XETG00118__0022474__Region_1__20250418__203103`
- **Aging_Metab Run 7** `07. Xenium Spatial Transcriptomics / 2 Aging_Metab_Xenium / Run7_20260417__210457__20260417_AgeXMetabolism_Run7 / output-XETG00118__0069080__Region_1__20260417__210507`

## Notes

- Sample-level assignment is **not** part of the XOA export. Both slides export a single `Region_1`. Samples are isolated post hoc via the ROI masks in `spatial_data/` (drawn in Xenium Explorer off the `.zarr.zip` files, using the layouts on page 2 of the PDF key), applied during ingest.
- Slide `0022474` originally exported its manifest as `experiment.xenium_slide1.xenium`; renamed on disk to `experiment.xenium` on 2026-05-12 so `LoadXenium()` resolves it. Slide `0069080` uses the standard `experiment.xenium`.
- Slide effect (Slide2 `0069080` vs Slide1 `0022474`): +25 % transcripts/cell, -7 % assignment rate, +33 % control-probe rate. The QC verdict attributes this to acquisition / tissue / firmware differences (analysis `3.2.0.7` → `3.3.0.1`), not algorithmic drift; the v3.2→v3.3 segmentation changes touch boundary-stain handling only.
- Sample `4s2M_F2` **is not complete**: only a tissue fragment remains on the slide. Flagged `partial_section = TRUE` at ingest; flag it in any cluster-composition or cell-count comparison.
