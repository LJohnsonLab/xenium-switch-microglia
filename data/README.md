# `data/` — processed analysis artefacts


Filename convention: `<YYYYMMDD>_<short-description>.<ext>`. 

## File type

**All seven `.qs2` files are Seurat v5 objects.** . Restore any of them with `qs2::qs_read()`.

## Current contents

Listed in pipeline order. "Built by" is the script/notebook whose `qs_save()` writes the file; "Loaded by" lists the notebooks that `qs_read()` it (as a reload checkpoint, denoted *self*, or as input to the next stage).

| File | Built by | Loaded by | What it holds |
|------|----------|-----------|---------------|
| `20260512_microglia_switch_original_seurat.qs2` | `R/20260512_xenium2seurat.R` | `20260512_QC.qmd`; `20260512_SCT_no_regression.qmd`; `20260512_SCT_regressed_by_slide.qmd`; `20260512_SpaNorm_sample005.qmd`; `20260513_SpaNorm_samplep025.qmd` | Merged Seurat v5 across both slides. 480 features (full `mBrain_480g`) × 429,980 cells after ingest filters (cell in one of the seven sample ROIs, `nCount_Xenium > 0`, `median_qv >= 20`). `Xenium` counts assay + per-feature-class assays (`BlankCodeword`, `ControlCodeword`, `ControlProbe`) and two FOV slots (one per slide). Metadata: `sample_id`, `slide`, `condition`, `partial_section` (TRUE only for `4s2M_F2`), `nCount_raw`, `nFeature_raw`, `median_qv`, `cell_area`. **No normalization, no clustering.** The common ancestor of every downstream build. |
| `20260512_microglia_switch_sct.qs2` | `20260512_SCT_no_regression.qmd` | *self* | `original_seurat` + SCTransform (no covariate) → PCA(30) → Louvain(0.5) → UMAP. Normalization-comparison arm. |
| `20260512_microglia_switch_sct_regressed_by_slide.qs2` | `20260512_SCT_regressed_by_slide.qmd` | *self* | As above but SCTransform with `vars.to.regress = "slide"`. Normalization-comparison arm. |
| `20260512_microglia_switch_spanorm.qs2` | `20260512_SpaNorm_sample005.qmd` | *self* | Per-sample SpaNorm at `sample.p = 0.05` (seven fits) + PCA → clustering → UMAP. The "0.05 variant" checkpoint, built on demand. **Note the name has no `_p025`/`_p005` suffix**; it predates the `sample.p` sweep naming. |
| `20260512_microglia_switch_spanorm_p025.qs2` | `20260513_SpaNorm_samplep025.qmd` | *self*; `20260516_SpaNorm_samplep025_clean.qmd` | Per-sample SpaNorm at `sample.p = 0.25` (package default; promoted 2026-05-14). `SpaNorm` v5 assay + PCA(30) → Louvain(0.5) → UMAP, 25 clusters. **Parent SpaNorm build.** Input to the clean re-fit. |
| `20260516_microglia_switch_spanorm_p025_clean.qs2` | `20260516_SpaNorm_samplep025_clean.qmd` | *self*; `20260517_microglia_subcluster.qmd`; `20260518_microglia_APOE_per_cell_call.qmd` (spatial-context QC) | Parent SpaNorm `p025` with the artefactual cluster-20 cells dropped and SpaNorm re-fit per sample. **Project default scaffold.** Carries both `clusters_res05` (22 clusters) and `clusters_res07` (30 clusters, working default). |
| `20260517_microglia_subcluster.qs2` | `20260517_microglia_subcluster.qmd` | *self*; `20260518_microglia_pseudobulk_DE.qmd`; `20260518_microglia_APOE_per_cell_call.qmd` | Cluster-2 microglia subset (~36k cells) of the clean build, `DietSeurat`'d to Xenium counts, SpaNorm re-fit per sample on microglia. Carries `clusters_res05` (11 clusters) and `clusters_res03` (7 clusters, microglia default). Input to both DE notebooks and the APOE switch call. |

## Lineage

```
xenium_rawData/  ──(R/20260512_xenium2seurat.R)──▶  original_seurat
                                                        │
        ┌───────────────────────────────────┬──────────┼──────────────────────────┐
        ▼                                    ▼          ▼                          ▼
      sct                      sct_regressed_by_slide  spanorm (p=0.05)     spanorm_p025 (p=0.25)
   (comparison)                    (comparison)        (comparison)                │
                                                                                   ▼
                                                                        spanorm_p025_clean  (default scaffold)
                                                                                   │
                                                                                   ▼
                                                                          microglia_subcluster
                                                                                   │
                                                                   ┌───────────────┴───────────────┐
                                                                   ▼                               ▼
                                                        pseudobulk_DE.qmd              APOE_per_cell_call.qmd
                                                        (results/ CSVs)                (results/ CSVs)
```



## Conventions

- `_original_seurat`: first Seurat object built straight from `LoadXenium()` + ROI assignment, QC filters applied, no normalization or dimensionality reduction. Downstream artefacts get descriptive suffixes.
- SpaNorm fits run **per sample** (seven fits), never per slide; the `SpaNorm` v5 assay holds only a `data` layer, so an explicit `ScaleData()` precedes `RunPCA()` and `VariableFeatures()` is set to all 480 genes.
- Prefer `.qs2` over `.rds` for Seurat objects; `qs2::qs_save()` is faster and smaller.
