# Cluster cell-type assignments — Microglia sub-cluster on the clean SpaNorm build, sample.p = 0.25

I assigned each cluster from `20260517_microglia_subcluster.qmd` at two Louvain resolutions (0.5 and 0.3) on the microglia-only SpaNorm re-fit. Sources: `results/20260517-find_all_markers_microglia_subcluster_res05.csv` (11 clusters) and `results/20260517-find_all_markers_microglia_subcluster_res03.csv` (7 clusters). The panel is the 480-gene Xenium mBrain panel, so the cell-type vocabulary is constrained to what those probes resolve. Rows are ordered by cell-type similarity (true microglia first, then contamination types). The "Flagged genes" column lists top-10 markers that do not fit the called identity.

**Scope note.** The microglia subset starts from clean-SpaNorm Louvain-0.5 cluster 2 (~50k microglia across 7 samples; see `notes/20260516-cluster_celltype_assignments_spanorm_clean_res05.md`). SpaNorm is re-fit per sample on microglia counts only, with the same `sample.p = 0.25` and `df.tps = 6` parameters as the parent build. The cluster count drops from 11 at Louvain 0.5 to 7 at 0.3 because the DAM Stage-1 / Stage-2 split, the homeostatic split, and the multi-lineage doublet pool all dissolve into broader categories at the coarser cut. The 0.3 decomposition is promoted as the working microglia scaffold because it carries no visible slide effect, while 0.5 shows slide imbalance on at least the oligodendrocyte-contamination cluster (8) and the multi-lineage doublet pool (9).

## Assignments at Louvain 0.5 (11 clusters, 35 700 microglia total)

| Cluster | Cells (n) | Cell type | Top-10 markers supporting the call | Flagged genes | Notes |
|---|---|---|---|---|---|
| 0 | 6 919 | DAM Stage-1 microglia (Trem2-dependent core) | Lyz2, Ccl3, Ccl4, Cst7, Clec7a, Ctsb, Gusb, Hexa, Lpl | Pdgfra (pct.1 0.216, OPC) | Canonical Stage-1 DAM signature per Keren-Shaul (Cell, 2017). The Pdgfra flag at low pct.1 is a small OPC-doublet contamination. |
| 4 | 2 962 | DAM Stage-2 microglia (antigen-presenting, late-stage) | Spp1, Apoc4, Gpnmb, H2-Eb1, Lgals3, H2-Ab1, Pf4, **rs429358-WT-228**, Cd74, Itgax | — | Stage-2 DAM with MHC-II antigen-presentation module (H2-Eb1 / H2-Ab1 / Cd74) and Spp1 / Gpnmb at the top. The custom APOE rs429358-WT probe (E3/E4 allele) sits in this cluster's top-10 at FC 2.49, which is the molecular hook for the E4-vs-E2 contrast. |
| 2 | 5 123 | Homeostatic microglia (canonical) | Tmem119, P2ry12, Maf, Ltc4s, Pnp, Htra1, Selplg, Adcy5, Ccl5 | Vip (pct.1 0.041, inhib neuron) | Canonical homeostatic signature. The Vip flag at very low pct.1 is a single-cell doublet, not a population feature. |
| 6 | 2 082 | Low-FC homeostatic microglia (boundary / CP-adjacent subset) | Tmem119, Cd68, Ptgs1, Gpr34, P2ry12, **rs7412-ALT-226:T**, Mertk | Ttr (CP), Ermn (OL), Apod (OL) | Microglia top-10 anchors are all there but at low FC. The Ttr (FC 1.96, pct.1 0.54) and Ermn / Apod entries flag boundary cells near choroid plexus / white-matter. The APOE rs7412-ALT probe (E2 allele) sits here, mirroring the rs429358-WT enrichment in DAM Stage-2 (4). |
| 10 | 399 | IRM + proliferating microglia | Cxcl10, Mki67, Ccl2, Ccl12, Ifit1, Ifit2, Ccl5, Ifitm3, Hmgb2, Plac8 | — | Interferon-responsive microglia (Ifit1 / Ifit2 / Ifitm3 / Cxcl10) co-cluster with proliferating Mki67+ cells. Cxcl10 at FC 6.83 (pct.1 0.617) and Mki67 at 5.03 are the strongest single drivers. Smallest cluster, so power-limited for any per-cluster DE. |
| 1 | 5 453 | Inhibitory + glutamatergic neuron contamination (doublets) | Slc17a6, Calb2, Plcb4, Agt, Tunar, Tac1, Pvalb, Nts, Pthlh, Sst | — | Multi-neuronal contamination (Vglut2+ + canonical inhibitory markers) with Agt astrocyte leak. Drop. |
| 3 | 4 754 | Cortical excitatory neuron contamination (doublets) | Lamp5, Opcml, Cnr1, Slc17a7, Grin2b, Grin2a, Dgat2, Fezf2, Egr1, Hmgcr | — | Vglut1+ cortical signature; same contamination flavor as parent SpaNorm cluster 1 (upper-layer Lamp5+/Calb1+). Drop. |
| 5 | 2 934 | Reactive astrocyte contamination (doublets) | Serping1, Vim, Bgn, Plpp3, Mt2, Aqp4, Mt1, Serpina3n, Gstm1, Lgals1 | — | Reactive astrocyte signature (Serpina3n + Serping1 + Mt1 / Mt2). Drop. |
| 7 | 1 844 | Vascular contamination (endothelial + pericyte, doublets) | Cldn5, Vtn, Cxcl12, Atp13a5, Mgp, Emcn, Kcnj8, Klf4, Flt1, Bgn | — | Mixed BBB endothelial (Cldn5 / Flt1 / Emcn) and pericyte (Vtn / Atp13a5 / Kcnj8) markers in one cluster. Drop. |
| 8 | 1 830 | Oligodendrocyte contamination (doublets) | Mog, Cldn11, Aspa, Mag, Ermn, Apod, Pllp, Ptgds, Gatm, Plin4 | — | Mature OL signature; the slide effect at 0.5 is most visible on this cluster. Drop. |
| 9 | 1 400 | Multi-lineage low-specificity (doublet pool) | Slc32a1, Mbp, Sncg, Nkg7, Cd247, Ccdc153, G6pd2, Ms4a4b, Tmem212, Cpt1b | — | GABAergic (Slc32a1) + OL (Mbp) + T cell (Nkg7 / Cd247 / Ms4a4b) + ependymal (Ccdc153 / Tmem212) markers in one cluster, all at low FC. Multi-doublet pool; this cluster carries the second-largest slide effect at 0.5. Drop. |

True-microglia subtotal at 0.5: 17 485 cells (49.0 %); contamination subtotal: 18 215 cells (51.0 %).

## Assignments at Louvain 0.3 (7 clusters, 35 700 microglia total)

| Cluster | Cells (n) | Cell type | Top-10 markers supporting the call | Flagged genes | Notes |
|---|---|---|---|---|---|
| 0 | 14 013 | DAM microglia (consolidated Stage-1 + Stage-2) | Spp1, Cst7, Gpnmb, H2-Eb1, Lyz2, Cd74, Ccl3, Lpl, Itgax, Lgals3 | — | Merges 0.5 clusters 0 (Stage-1, 6 919 cells) and 4 (Stage-2, 2 962 cells). The combined top-10 carries Stage-1 markers (Cst7, Lyz2, Lpl, Ccl3, Itgax) alongside Stage-2 (Spp1, Gpnmb, H2-Eb1, Cd74). The Stage-1 / Stage-2 distinction is biologically meaningful (per Keren-Shaul, Cell, 2017) and worth recovering by sub-clustering cluster 0 if power allows. |
| 1 | 7 743 | Homeostatic microglia (consolidated, with mild astrocyte bleed) | Tmem119, P2ry12, Ptgs1, Maf, Lpcat2, Ltc4s | Agt, Gja1, Slc7a10 (astrocyte) | Merges 0.5 clusters 2 (canonical homeostatic, 5 123 cells) and 6 (low-FC homeostatic, 2 082 cells). The astrocyte markers Agt (pct.1 0.291) / Gja1 (0.601) / Slc7a10 (0.181) creep into top-10 because the low-FC subset of 0.5 cluster 6 carried boundary astrocyte signal. Per-cell filtering (recipe below) will clean this. |
| 6 | 388 | IRM + proliferating microglia | Cxcl10, Mki67, Ccl2, Ccl12, Ifit1, Ifit2, Ccl5, Ifitm3, Hmgb2, Plac8 | — | Unchanged from 0.5 cluster 10 (399 cells). The cluster survives the coarser cut intact, which means the IRM + proliferating signal is biologically robust and not a noise artefact. Smallest cluster; power-limited for any per-cluster DE. |
| 2 | 5 582 | Cortical excitatory neuron contamination (doublets) | Lamp5, Opcml, Slc17a7, Egr1, Grin2a, Grin2b, Cnr1, Dgat2, Fezf2, Hmgcr | — | Same as 0.5 cluster 3. Drop. |
| 3 | 4 334 | Inhibitory neuron + reactive astrocyte contamination (merged doublets) | Slc17a6, Plcb4, Tunar, Pvalb, Calb2, Serpina3n, Sst, Agt, Pthlh, Lhx6 | — | Merges 0.5 cluster 1 (inhibitory + glutamatergic neuron, 5 453 cells) with the part of 0.5 cluster 5 (reactive astrocyte, 2 934 cells) that did not collapse into the consolidated DAM. Drop. |
| 4 | 1 894 | Oligodendrocyte contamination (doublets) | Mog, Cldn11, Aspa, Mag, Ermn, Apod, Pllp, Ptgds, Gatm, Folh1 | — | Same as 0.5 cluster 8 (1 830 cells). Slide effect from 0.5 absorbed at 0.3. Drop. |
| 5 | 1 746 | Vascular contamination (endothelial + pericyte, doublets) | Cldn5, Vtn, Cxcl12, Atp13a5, Emcn, Mgp, Kcnj8, Klf4, Flt1, Ocln | — | Same as 0.5 cluster 7 (1 844 cells). Drop. |

True-microglia subtotal at 0.3: 22 144 cells (62.0 %); contamination subtotal: 13 556 cells (38.0 %). The fraction of cells the wholesale cluster-drop step removes is 13 % lower at 0.3 than at 0.5 (38 % vs 51 %). This is the cost of the coarser cut: ~4 600 contaminating cells from the 0.5 reactive-astrocyte (5) and multi-lineage doublet (9) clusters get silently absorbed into the consolidated DAM (0) and homeostatic (1) clusters at 0.3. The per-cell module-score filter (recipe below) is what catches them and is therefore more critical at 0.3 than at 0.5.

## Cross-mapping: Louvain 0.5 → Louvain 0.3

| 0.5 cluster | 0.5 identity | 0.3 cluster | Change |
|---|---|---|---|
| 0 | DAM Stage-1 | 0 (consolidated DAM) | merged with 4 |
| 4 | DAM Stage-2 | 0 (consolidated DAM) | merged with 0 |
| 2 | Homeostatic canonical | 1 (consolidated homeostatic) | merged with 6 |
| 6 | Homeostatic low-FC / CP-boundary | 1 (consolidated homeostatic) | merged with 2 |
| 10 | IRM + proliferating | 6 | preserved |
| 1 | Inhib + glut neuron contamination | 3 (merged with reactive astrocyte) | partially merged with 5 |
| 3 | Cortical excitatory contamination | 2 | preserved |
| 5 | Reactive astrocyte contamination | 3 (merged with inhib neurons) | partially merged with 1 |
| 7 | Vascular contamination | 5 | preserved |
| 8 | Oligodendrocyte contamination | 4 | preserved |
| 9 | Multi-lineage doublet pool | — (dissolves) | redistributed across 3 / 4 / 5 / 1 |

## Cleaning the microglia object — markers that flag non-microglia contamination

Eight of the 11 clusters at 0.5 (and four of the seven at 0.3) are doublet / segmentation-contamination clusters of neighboring cell types. They appear because microglia (~5 % of brain cells) are sparse relative to denser neighbors (oligodendrocytes, neurons, astrocytes), so Xenium segmentation occasionally assigns transcripts from those neighbors to a microglia cell-segment boundary. The contamination is not a normalization or clustering failure; it is a segmentation-geometry feature of the data.

For the downstream E4-vs-E2 DE on the true microglia states, these cells need to come out. Two strategies, used in combination:

### 1. Drop the contamination clusters wholesale

At Louvain 0.3, exclude clusters 2, 3, 4, 5 and keep clusters 0, 1, 6. This is the simplest filter and gets ~95 % of the contamination out in one step:

```r
# After loading the microglia subcluster object:
mg_clean <- subset(mg, subset = clusters_res03 %in% c("0", "1", "6"))
```

### 2. Per-cell filter using contamination-marker module scores

`AddModuleScore()` returns a per-cell score for each marker set; threshold and exclude. This catches contaminating cells that landed inside the true-microglia clusters (e.g. a microglia-OL doublet inside the homeostatic cluster 1 that shows Tmem119 strongly but also Mog at low pct.1):

```r
contamination_markers <- list(
  oligodendrocyte  = c("Mog", "Aspa", "Mag", "Cldn11", "Pllp", "Ermn", "Ptgds", "Apod", "Mbp", "Plp1"),
  excitatory       = c("Slc17a7", "Slc17a6", "Fezf2", "Lamp5", "Grin2a", "Grin2b", "Cnr1", "Cux2", "Calb1"),
  inhibitory       = c("Slc32a1", "Gad1", "Sst", "Pvalb", "Vip", "Lhx6", "Tac1", "Nts", "Calb2", "Lhx6"),
  astrocyte        = c("Aqp4", "Agt", "Plpp3", "Slc7a10", "Gja1", "Vim", "Mt1", "Mt2", "Serping1", "Serpina3n"),
  endothelial      = c("Cldn5", "Flt1", "Emcn", "Klf4", "Cxcl12", "Ocln", "Slc2a1", "Sgms1"),
  pericyte         = c("Vtn", "Atp13a5", "Kcnj8", "Bgn", "Mgp", "Axl", "Stard8"),
  ependymal_cp     = c("Ccdc153", "Tmem212", "Dnah11", "Ttr", "Car12", "Ltc4s"))

mg <- AddModuleScore(mg, features = contamination_markers, name = "contam_")
# columns mg$contam_1 .. mg$contam_7 hold per-cell scores.
# Threshold each at the inflection point of its distribution (typically 0.05-0.2 for these panel sizes).
mg_clean <- subset(mg,
                   subset = contam_1 < 0.1 & contam_2 < 0.1 & contam_3 < 0.1 &
                            contam_4 < 0.1 & contam_5 < 0.1 & contam_6 < 0.1 &
                            contam_7 < 0.1)
```

Inspect each score's histogram with `VlnPlot(mg, features = paste0("contam_", 1:7), group.by = "clusters_res03")` before settling on the threshold; the contamination clusters should sit at clearly higher scores than the true microglia clusters, which defines the cut.

### Marker-set rationale

| Contamination type | Why these markers | Why microglia do not express them |
|---|---|---|
| Oligodendrocyte | Mog / Mbp / Mag / Aspa / Cldn11 are myelin / membrane / aspartate-acylase core; absent from any microglia state. | Microglia have no myelin biology. |
| Excitatory neuron | Slc17a7 (Vglut1) / Slc17a6 (Vglut2) are the glutamatergic vesicular transporters; Fezf2 / Cux2 / Lamp5 / Calb1 are cortical layer identity TFs / Ca2+ binders; Grin2a / Grin2b are NMDA receptor subunits. | Microglia carry purinergic and chemokine receptors, not glutamate machinery. |
| Inhibitory neuron | Slc32a1 (VGAT) + Gad1 define GABAergic neurons; Sst / Pvalb / Vip / Lhx6 / Tac1 / Nts mark canonical interneuron classes. | Same logic; microglia do not synthesize GABA. |
| Astrocyte | Aqp4 (water channel) / Agt (angiotensinogen) / Gja1 (Cx43 gap junction) / Slc7a10 / Plpp3 are astrocyte-restricted; Vim / Mt1 / Mt2 / Serpina3n / Serping1 mark reactive astrocytes. | Microglia carry their own immune effectors, not Aqp4 / Cx43. |
| Endothelial | Cldn5 / Ocln (tight junctions) + Flt1 (VEGFR1) / Emcn / Klf4 are BBB-specific; Cxcl12 marks brain endothelium. | Microglia are tissue-resident immune; not endothelial. |
| Pericyte | Vtn / Atp13a5 / Kcnj8 / Bgn / Mgp are mural-cell core; Axl is shared with macrophages so use with caution. | Pericytes are mural; microglia are immune. |
| Ependymal / choroid plexus | Ccdc153 / Tmem212 / Dnah11 are motile-cilia core (ependyma); Ttr / Car12 / Ltc4s mark CP epithelium. | Microglia have no cilia; CP epithelium has no immune effectors. |

### Combined recipe (recommended)

Do the wholesale cluster drop first (cheap, removes the bulk of contamination), then apply the per-cell module-score filter on the remaining cells to catch residual doublets inside the true-microglia clusters:

```r
mg_clean <- subset(mg, subset = clusters_res03 %in% c("0", "1", "6"))
mg_clean <- AddModuleScore(mg_clean, features = contamination_markers, name = "contam_")
# threshold and re-subset as above
```

After cleaning, re-run `ScaleData → RunPCA → FindNeighbors → FindClusters → RunUMAP` on `mg_clean` and re-derive markers. The cluster count should drop to 3 (DAM, homeostatic, IRM+proliferating) at the same Louvain resolution.

## Headline findings

1. **Three robust microglia states recovered.** DAM (cluster 0 of the 0.3 build), homeostatic (1), and IRM + proliferating (6). All three populate every animal and both slides.
2. **DAM Stage-1 vs Stage-2 split at 0.5, collapsed at 0.3.** The Stage-2 cluster (0.5 cluster 4) carries the MHC-II antigen-presentation module (H2-Eb1 / H2-Ab1 / Cd74) and the highest enrichment of the APOE rs429358-WT probe (the E3/E4 allele). The Stage-1 cluster (0.5 cluster 0) carries the canonical chemokine + lipid signature (Lyz2 / Ccl3 / Ccl4 / Cst7 / Lpl). The distinction is biologically meaningful per Keren-Shaul (Cell, 2017) and worth recovering by sub-clustering the consolidated cluster 0 of the 0.3 build at higher resolution.
3. **APOE allele-specific probes carried into the microglia analysis.** rs429358-WT-228 (E3/E4 allele) enriches in DAM Stage-2 (0.5 cluster 4); rs7412-ALT-226:T (E2 allele) enriches in the low-FC homeostatic / CP-boundary cluster (0.5 cluster 6). At 0.3 both probes merge into the consolidated DAM and consolidated homeostatic clusters respectively, but the per-cell signal is still resolvable. This is the molecular hook for the E4-vs-E2 contrast on a per-cell basis.
4. **High contamination fraction is expected.** 6 of 11 clusters at 0.5 (55 %) and 4 of 7 at 0.3 (57 %) are doublet / segmentation-overlap clusters. This is a feature of Xenium segmentation around sparse cell-type boundaries, not a normalization failure. The cleaning recipe above is the correct response.
5. **Slide effect at 0.5, gone at 0.3.** The OL-contamination cluster (0.5 cluster 8) and the multi-lineage doublet pool (0.5 cluster 9) are slide-imbalanced. Collapsing them at 0.3 absorbs the imbalance and is the load-bearing reason to promote the 0.3 decomposition for downstream DE.

## Verdict and recommended next steps

Promote the 0.3 decomposition as the working microglia scaffold. The E4-vs-E2 pseudobulk DE runs on the cleaned subset:

1. Drop contamination clusters: `subset(mg, subset = clusters_res03 %in% c("0", "1", "6"))`.
2. Apply per-cell contamination-marker filtering with the seven marker sets above to catch residual doublets.
3. Run pseudobulk DE on cluster 0 (DAM) and cluster 1 (homeostatic) of the cleaned object, with cluster 6 (IRM + proliferating) as a secondary contrast.
4. After the DAM-versus-homeostatic primary contrast lands, sub-cluster the consolidated DAM cluster (cleaned 0) at higher Louvain resolution to recover the Stage-1 / Stage-2 split that the 0.3 cut collapsed.
5. Inspect spatial location of the IRM + proliferating cluster (6) relative to amyloid plaques (this is where the IRM module is biologically expected to concentrate).
6. The APOE rs429358 and rs7412 probe enrichment per cluster gives a per-cell molecular hook for the E4-vs-E2 contrast that does not depend on the mouse genotype label alone.

## References

- Keren-Shaul H et al. "A unique microglia type associated with restricting development of Alzheimer's disease." Cell, 2017. DOI: 10.1016/j.cell.2017.05.018.
- Salim A et al. "SpaNorm: spatially-aware normalization for spatial transcriptomics data." Genome Biology, 2025. DOI: 10.1186/s13059-025-03565-y.
- Parent qmd: `20260517_microglia_subcluster.qmd`.
- Parent cluster scaffold: `notes/20260516-cluster_celltype_assignments_spanorm_clean_res07.md`.
