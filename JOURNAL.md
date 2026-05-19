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

## 2026-05-17 — Clean SpaNorm re-fit, microglia sub-cluster, and the contamination cleaning recipe

**Context & changes.** Built `20260516_SpaNorm_samplep025_clean.qmd` (`data/20260516_microglia_switch_spanorm_p025_clean.qs2`): drops the artefactual cluster-20 cells from the parent SpaNorm build, re-fits SpaNorm per sample, brackets Louvain at 0.5 (22 clusters) and 0.7 (30 clusters). Promoted 0.7 because it resolves homeostatic vs DAM microglia (12 / 6), recovers the discrete T-cell cluster (28: Cd3d / Gzmb / Nkg7 / Itk) that the 0.5 cut lost, and re-splits Sst+ / Vip+ / Pvalb+ GABAergic. The user-requested "cluster 25 removal" silently no-opped (parent only had 0-24); only cluster 20 dropped. Then built `20260517_microglia_subcluster.qmd` (`data/20260517_microglia_subcluster.qs2`): subsets clean cluster 2 (~36k microglia), `DietSeurat` to Xenium only, re-fits SpaNorm per sample on microglia (`df.tps = 6` converges at ~7k cells / section), brackets Louvain at 0.5 (11 sub-clusters) and 0.3 (7). Promoted 0.3 because slide imbalance at 0.5 (OL-contamination, multi-lineage doublet) absorbs at 0.3. Three robust microglia states: DAM (0), homeostatic (1), IRM + proliferating (6: Cxcl10 / Mki67 / Ifit1 / Ifit2). Custom APOE probes rs429358-WT (E3/E4) and rs7412-ALT (E2) enrich per cluster, giving a per-cell molecular hook for the E4-vs-E2 contrast. Cell-type assignment notes moved from `results/` to new `notes/` folder; marker CSVs stay in `results/`. `fiddle()` x-axis text size 6 → 12.

**Cleaning recipe.** Two-strategy contamination filter in `notes/20260517-cluster_celltype_assignments_microglia_subcluster.md`: (1) wholesale drop of contamination clusters 2 / 3 / 4 / 5, keeping 0 / 1 / 6 at Louvain 0.3, then (2) per-cell `AddModuleScore()` against seven contamination-marker sets (OL, excitatory, inhibitory, astrocyte, endothelial, pericyte, ependymal/CP). No DoubletFinder; the 480-probe panel breaks its kNN assumption.

**Unresolved.** Pseudobulk E4-vs-E2 DE on cleaned microglia 0 / 1 / 6 at 0.3. Sub-cluster consolidated DAM to recover Stage-1 vs Stage-2 (visible at 0.5 as cluster 4 Spp1 / Gpnmb / H2-Eb1 / Cd74 vs 0 Lyz2 / Cst7 / Ccl3 / Lpl). Spatially inspect IRM cluster near plaques. Sub-cluster astrocyte pair (clean-0.7 clusters 2, 8) for DAA signature. Cross-validate APOE probe per-cell calls against genotype.

**Token Budget:** ~340 tokens.

## 2026-05-19 — Microglia pseudobulk DE, corrected APOE probe labels, and per-cell ROC-tuned E2 call

**Context & changes.** Extended `20260518_microglia_pseudobulk_DE.qmd` with a full-microglia pool DE section (clusters 0, 2, 4, 6, 10 merged) ahead of the per-state loop, plus per-cluster interpretive commentary that splits each volcano's top-5 hits into expected microglial biology (IFN / DAM / homeostatic / E2-probe signals) and likely segmentation contamination (neuronal Cnr1 / Slc17a7 / Plcb4 / Tunar, OPC Pdgfra, T-cell Ms4a4b). Built and committed the new `20260518_microglia_APOE_per_cell_call.qmd`: probe-validation histograms + summary table against the E4/E4 background, then a logistic-regression score on the two T-allele E2-evidence probes (rs429358-WT-228 + rs7412-ALT-226:T) trained on a per-cell proxy E2-vs-E4 label (4s2M true microglia vs 4s2 true microglia), with pROC AUC + DeLong + bootstrap CI, three operating points (Youden / FPR=5% / FPR=1%), per-cell `e2_call` attached to the microglia object, a genotype-on-cell pseudobulk DE, and an `AddModuleScore`-based E2 UMAP split by condition and by sample.

**Key correction.** Confirmed with JM that `rs429358-WT-228` hybridizes to the T allele at codon 112, so it detects E2 (not E4 as previously labelled across CLAUDE.md, the pseudobulk DE notebook, and the JOURNAL). Both custom E2-evidence probes are now T-allele readouts; the "unexpected E4 signal in switched arm" mystery dissolved. Corrected the labels across all three files and saved a project memory.

**Unresolved.** Sub-cluster the consolidated DAM at higher Louvain (Stage-1 vs Stage-2 split visible at 0.5 cluster 4 Spp1 / Gpnmb / H2-Eb1 / Cd74 vs cluster 0 Lyz2 / Cst7 / Ccl3 / Lpl). Spatially inspect the IRM cluster near amyloid plaques. Sub-cluster the SpaNorm astrocyte pair (clean-0.7 clusters 2, 8) for an IFN-responsive DAA signature. Cross-check the genotype-on-cell DE top hits against the parent five-state DE; recurrent hits across the two designs are the most robust candidates.

**Token Budget:** ~250 tokens.
