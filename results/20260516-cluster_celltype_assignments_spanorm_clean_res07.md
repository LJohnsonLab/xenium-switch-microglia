# Cluster cell-type assignments — SpaNorm clean (cluster 20 removed), sample.p = 0.25, resolution 0.7

I assigned each of the 30 Louvain clusters from `20260516_SpaNorm_samplep025_clean.qmd` at the higher Louvain resolution (0.7) using the top 10 FindAllMarkers per cluster (ranked by `avg_log2FC`; source: `results/20260516-find_all_markers_spanorm_p025_clean_res07.csv`). The panel is the 480-gene Xenium mBrain panel, so the cell-type vocabulary is constrained to what those probes resolve. I ordered the rows by cell-type similarity, not by cluster number, so related populations sit next to each other. The "Flagged genes" column lists top-10 markers that do not fit the called identity and that I would scrutinize as possible doublet or segmentation contamination before downstream DE.

**Scope note.** The clean build dropped only the cluster-20 cells from the parent SpaNorm object (the `seurat_clusters %in% c("20", "25")` filter silently no-opped on "25", which did not exist in a 0-to-24-level parent). The Louvain resolution was then raised from 0.5 to 0.7 on the same clean PCA / SNN graph, producing 30 clusters versus the 22 clusters of the 0.5 decomposition documented in `20260516-cluster_celltype_assignments_spanorm_clean.md`. The 0.5 build is left in the qmd as the intermediate step that motivated this higher-resolution re-run.

## Assignments

| Cluster | Cell type | Top-10 markers supporting the call | Flagged genes (possible contamination) | Notes |
|---|---|---|---|---|
| 13 | Oligodendrocyte precursor cells (OPCs) | Pdgfra, Tnr, Plpp4, Gpt2, Chst11, Dpyd, Pllp, Mki67, Arsb, Lrp1 | — | Pdgfra at FC 6.41 anchors the call; Mki67 at pct.1 0.061 flags the small proliferative OPC fraction. Identical signature to the 0.5 clean cluster 13. |
| 0  | Mature oligodendrocytes (high-FC core) | Mag, Aspa, Mog, Gatm, Ermn, Cldn11, Pllp, Mbp, Plin4, Cers2 | — | Textbook myelinating-OL signature; identical to 0.5 clean cluster 0 and parent cluster 0. The DAM / DAM-like T-cell tracer markers (Nkg7, Ms4a4b) flagged at 0.5 are gone from the top-10 here, consistent with the T-cell signal having been pulled out into cluster 28 by the higher resolution. |
| 9  | Mature oligodendrocytes (lower-FC twin) | Cldn11, Apod, Ermn, Mog, Plin4, Aspa, Gatm, Pllp, Mag, Ptgds | — | Same gene set as cluster 0 but FCs collapse from ~3.5 to ~2.0-3.4. Maps onto 0.5 clean cluster 9 and parent cluster 9. |
| 19 | Mature oligodendrocytes, third state (very low FC) | Mag, Aspa, Mog, Gatm, Pllp, Cldn11, Ermn, Plin4, B3galt5, Cers2 (all FC ≤ 2.06) | — | NEW vs both 0.5 clean and parent. Same myelin signature as clusters 0 and 9 but FCs collapse further (max FC 2.06; pct.1 0.94 vs pct.2 0.19 for Mag). I read this as a third OL maturation or regional state that 0.5 lacked the resolution to separate. Confirm spatially before treating as a distinct biological state. |
| 2  | Astrocytes (gray-matter, Aqp4/Gja1 high) | Plpp3, Htra1, Gja1, Prodh, Slc7a10, Mt2, Cbs, Aqp4, Pla2g7, Hsd11b1 | — | Canonical astrocyte. Cxcl10 dropped out of top-10 here (it was the IFN-responsive hint in the 0.5 clean cluster 4); Hsd11b1 enters as the new top-10 entrant. Maps onto 0.5 clean cluster 4 and parent cluster 5. |
| 8  | Astrocytes (Agt+ subtype) | Agt, Slc7a10, Gja1, Aqp4, Pla2g7, Atp1b2, Cbs, Plpp3, Prodh, Htra1 | — | Agt at FC 5.05 distinguishes from cluster 2. Identical top-10 to 0.5 clean cluster 8 and parent cluster 8. |
| 12 | Microglia (homeostatic) | Tmem119, P2ry12, Selplg, Siglech, Cx3cr1, Il1a, Gpr34, Lpcat2, Ptgs1, Nlrp3 | — | **NEW SPLIT.** Tmem119 (FC 5.43, pct.1 0.881), P2ry12 (5.20), Selplg (4.62), Siglech (4.58) and Cx3cr1 (4.50) define the homeostatic-microglia state cleanly. Neither parent nor 0.5 clean isolated this state on its own. The unexpected Il1a (FC 4.49, pct.1 0.365) and Nlrp3 (3.93) hint at a low-level inflammasome priming inside the homeostatic pool. |
| 6  | Microglia (DAM / disease-associated) | Cst7, Ccl3, Itgax, Clec7a, Ch25h, Lyz2, Bcl2a1b, Ccl4, Cd14, Trem2 | — | **NEW SPLIT.** Textbook DAM signature: Cst7 (FC 6.37), Itgax (5.57), Clec7a (5.57), Ch25h (5.35), Trem2 (4.80), with the Ccl3 / Ccl4 chemokine pair signaling active inflammation. This is the biggest scientific gain of the 0.7 re-clustering and the most actionable cluster for the E4-vs-E2 comparison on a 5XFAD background. Sub-cluster before DE to separate Stage-1 (Trem2-dependent) from Stage-2 (Trem2-independent) DAM if power allows. |
| 28 | T cells / cytotoxic lymphocytes | Cd3d, Ccl5, Itk, Nkg7, Ms4a4b, Tnf, Gzmb, Klrb1c, Cd247, Plac8 | — | **RECOVERED.** Cd3d at FC 7.52 (pct.1 0.72) and Gzmb at 4.93 anchor an unambiguous cytotoxic-T-cell signature, comparable to the parent's cluster 24 (Cd3d FC 7.24) and stronger on Ccl5 (7.37 vs 5.40 in the parent). The 0.5 clean build had lost this cluster into the VLMC neighborhood at pct.1 ≤ 0.03; raising Louvain to 0.7 restores it cleanly. Inspect spatial distribution near plaques as planned for the parent build. |
| 4  | Endothelial cells (BBB) | Cldn5, Cxcl12, Flt1, Emcn, Clec2d, Klf4, Slc2a1, Ocln, Sgms1, Fas | — | Tight-junction + BBB transporter signature; identical to 0.5 clean cluster 5 and parent cluster 6 (with Sgms1 / Slc2a1 in top-10 instead of the parent's Ifitm3). |
| 18 | Pericytes / mural cells | Kcnj8, Atp13a5, Vtn, Bgn, Ifitm1, Stard8, Mgp, Axl, Emcn, Ifitm3 | — | Pericyte-specific markers (Atp13a5 FC 7.33, Kcnj8 FC 7.83). Identical to 0.5 clean cluster 17 and parent cluster 16. |
| 7  | Vascular leptomeningeal cells / fibroblasts (VLMCs) + perivascular macrophages | Tagln, Slc47a1, Lum, Mrc1, Col1a2, Mgp, Pf4, Serping1, Serpinf1, Osr1 | — | Collagen + matrix + smooth-muscle (Tagln) + perivascular-macrophage (Mrc1, Pf4) signature. The T-cell contamination flagged at the 0.5 level (Cd3d / Itk / Nkg7 / Tnf at low pct.1) has been pulled out into cluster 28 here, so the VLMC top-10 is now clean. |
| 20 | Ependymal cells (with mild reactive flavor) | Ccdc153, Tmem212, Dnah11, Plin5, C3, Vim | Lgals3, Ttr, Sncg, Ifitm3 | Motile-cilia genes at FC 7-9 anchor the call. C3 and Lgals3 remain in top-10, hinting at reactive ependyma; Ttr tags adjacent CP at the ventricle wall. Still Slide2-enriched (ventricle sampling); identical to 0.5 clean cluster 18 and parent cluster 17. |
| 22 | Choroid plexus epithelial cells (with MHC-II) | Car12, Ltc4s, Adh1, Tgfbi, H2-Eb1, Ttr, H2-Ab1 | Htr2c, Dnah11, Bgn | Car12 (FC 7.33) and Ttr (FC 4.02) are CP-specific; MHC-II pair (H2-Eb1, H2-Ab1) confirmed in the higher-resolution top-10 as well, supporting the Kolmer-cell / CP-immune-niche reading. Identical to 0.5 clean cluster 19 and matches parent cluster 18 (with MHC-II as the addition). |
| 1  | Excitatory neurons, deep-layer cortex (Fezf2+) | Fezf2, Syt6, Oprk1, Slc17a7, Hmgcr, Fdft1, Egr1, Pfkp, Fut9, Ldha | — | Fezf2 (FC 3.76) + Syt6 mark L5/6; Hmgcr / Fdft1 / Pfkp / Ldha are metabolic, not strict contamination. Identical to 0.5 clean cluster 1 and parent cluster 2. |
| 5  | Excitatory neurons, upper-layer cortex (Cux2+/Lamp5+/Lypd1+/Calb1+) | Lypd1, Lamp5, Otof, Calb1, Pparg, Stard8, Cux2, Slc17a7, Junb, Grin2b | Pparg | Cux2 + Slc17a7 anchor L2/3; Lamp5 here is the glutamatergic flavor (Gad1 absent). The merged-with-inhibitory blend that the 0.5 clean cluster 6 carried is gone here because the inhibitory Lamp5+ / Lypd1+ cells separate into their own clusters (17 Sst+/Lhx6+ and 21 Vip+) at 0.7. Maps onto parent cluster 1. |
| 10 | Excitatory neurons, Vglut2+/Rorb+ (thalamic relay or L4-like) | Slc17a6, Plcb4, Calb2, Tunar, Rorb, Abca7, Calb1, Ldhb, Sptlc2, Adamts2 | — | Rorb on a Vglut2+ background fits thalamic relay > L4. Calb2 still here at pct.1 0.295 (parent cluster 3 contribution still bleeding in). Identical to 0.5 clean cluster 11; maps onto parent cluster 10. |
| 3  | Excitatory neurons, Calb2+/Tac1+/Nts+/Gck+ (mixed regional) | Gck, Nts, Calb2, Tac1, Slc17a6, Htr2c, Tox2, Gaa, Oprk1, Trh | — | Gck (FC 4.25, pct.1 0.146) is unusual outside hepatocytes / hypothalamic glucose-sensing neurons. Identical core to 0.5 clean cluster 3; combines parent cluster 3 (Calb2+/Slc17a6+) with Tac1 / Nts that were elsewhere in the parent. Trh at low pct.1 (0.083) hints at a few Pvalb-cluster-26 doublets but is biologically plausible in hypothalamic populations. |
| 24 | Excitatory neurons, Tshz2+ (recovered) | Tshz2, Dpp4, Htr2c, Pparg, Lamp5, Fut9, Aldob, Slc17a6, H2az1, Dcx | Aldob | **RECOVERED.** Tshz2 (FC 4.22, pct.1 0.964) re-emerges as a discrete cluster; the 0.5 clean build had absorbed it into cluster 15 (Rspo1+/Tcap+/IEG-rich). Maps cleanly onto parent cluster 19 (Tshz2+ deep/hippocampal). Aldob flags a small contamination of hepatic-flavored cells (pct.1 0.075). |
| 25 | Excitatory neurons, Wif1+/Syt6+/Calb2+/Sncg+ (NEW) | Wif1, Syt6, Calb2, Sncg, Tac1, Gys2, Slc17a6, St6galnac4, Ttr, Acly | Gys2 | NEW vs both 0.5 clean and parent. Wif1 (FC 5.80, pct.1 0.335) + Syt6 (5.13) + Calb2 (4.42) + Sncg (4.26) fits an upper-layer or hippocampal CA-region excitatory population with Wnt-pathway antagonism (Wif1) and Calb2 calbindin marking. Gys2 (glycogen synthase 2) is unusual in neurons and warrants checking. |
| 29 | Excitatory neurons, Nts+/Fezf2+/Pou3f1+/Tshz2+ deep-layer (NEW) | Nts, Fezf2, Pou3f1, Bhlhe40, Tshz2, Adamts2, Sorl1, Opcml, Apoc4, Got2 | Apoc4 | NEW vs 0.5 clean; partially traces back to the parent's stressed-neuron + Tshz2 / Nts overlap. Nts at FC 5.89 (pct.1 0.901) is the strongest single marker and identifies a deep-layer projection neuron population. Co-expression of Fezf2 and Tshz2 plus the AD-risk gene Sorl1 (FC 2.06) and the immediate-early gene Bhlhe40 makes this a candidate disease-associated excitatory neuronal cluster. Prioritize for the E4-vs-E2 comparison. |
| 14 | Excitatory neurons, stressed / disease-associated (Pou3f1+/Lpl+) | Pou3f1, Grin2a, Grin2b, Opcml, Chst1, Lpl, Slc17a7, Cnr1, Pfkl, Hmgcr | **Lpl**, Pou3f1 | Neuronal core (Grin2a/b, Opcml, Slc17a7, Cnr1) contaminated by reactive/lipid Lpl. Identical to 0.5 clean cluster 14 and parent cluster 13. Strong candidate for DoubletFinder before DE. |
| 15 | Rspo1+/Tcap+/IEG-rich excitatory | Rspo1, Tcap, Egr1, Whrn, Rorb, Hkdc1, Pparg, Fos, Lamp5, Junb | — | Rspo1 (FC 4.64) + Tcap (4.25) plus the IEG module (Egr1, Fos, Junb). Identical to 0.5 clean cluster 15. Loses the Tshz2 overlap because Tshz2+ cells now sit in cluster 24. Check spatial localization to rule out dissociation-induced IEG response. |
| 16 | Disease-associated reactive (DAM-like + neuronal + Tdo2) | Tdo2, Dgat2, Gpnmb, Plpp4, Grin2a, Jun, Bhlhe41, B4galt1, Lgals3, Serpina3n | Grin2a (neuronal in reactive cluster) | Reactive triad Gpnmb + Lgals3 + Serpina3n; same as 0.5 clean cluster 16 and parent cluster 15. Tdo2 is the parent-cluster-20 residue from neighbours of the dropped cells. With cluster 6 now isolating DAM cleanly, this reactive-mixed cluster should be re-evaluated as either glia-neuron doublets or a distinct DAM-adjacent astrocyte state. |
| 23 | Low-specificity / immediate-early-gene excitatory | Egr1, Lypd1, Fosb, Lamp5, Fezf2, Slc17a7, Fos, Rspo1, Numb, Otof (all FC ≤ 2.23) | — | IEG-dominated layer-marker cluster; FCs are low. Maps onto 0.5 clean cluster 20 and parent cluster 23. Merge or re-resolve. |
| 27 | Low-specificity, multi-lineage (Slc32a1+ / Mbp+) | Slc32a1, Mbp, Cldn11, Adh5, Prdx1, Pllp, Aldh9a1, Hmgcl, Acadl, Degs1 (all FC ≤ 1.52) | — | Slc32a1 (GABAergic vesicular transporter) + Mbp (myelin) in the same top-10 is a multi-lineage doublet signature. All FCs ≤ 1.52. The 0.5 build collapsed this kind of profile into clean cluster 10 (mixed GABAergic); at 0.7 it carves out as a separate low-specificity cluster. Almost certainly doublet-dominated; run DoubletFinder. |
| 17 | Sst+/Lhx6+ GABAergic interneurons (MGE-derived) | Lhx6, Sst, Gad1, Pvalb, Pthlh, Reln, Tcap, Tac1, Calb1, Lgals1 | — | **RECOVERED.** Lhx6 (FC 5.81) + Sst (5.46) + Gad1 (3.44) define MGE-derived inhibitory interneurons; Pvalb at FC 3.41 is the expected Lhx6+ co-expression and not the dedicated Pvalb cluster. Maps onto parent cluster 11. The 0.5 clean build had merged this into the mixed GABAergic cluster 10. |
| 21 | Vip+ GABAergic interneurons (CGE-derived) | Vip, Crh, Ndnf, Cnr1, Reln, Pthlh, Gad1, Macc1, Wif1, Lamp5 | — | **RECOVERED.** Vip (FC 8.50) + Crh (6.08) + Ndnf (4.60) anchor the call; canonical Vip+ CGE interneuron, sharper at 0.7 (Vip FC 9.82 in parent vs 8.50 here, both unambiguous). Maps onto parent cluster 21. |
| 26 | Pvalb+/Trh+ GABAergic interneurons | Trh, Pvalb, Gad1, Spp1, Sst, Whrn, Slc32a1, Gpx1, Serpinf1, Oprk1 | Spp1 | Pvalb (FC 5.71) + Trh (6.65) identifies fast-spiking PV+ interneurons with a Trh+ subset. Spp1 again flags as off-profile (DAM / reactive marker in a PV+ neuron). Identical to 0.5 clean cluster 21 and parent cluster 22. |
| 11 | GABAergic, Tac1+/Adcy5+/Pax6+ (striatal-like) | Ido1, Otof, Tac1, Adcy5, Cd36, Nts, Pax6, Gad1, Fosb, Csgalnact1 | Cd36 | Gad1 + Tac1 + Adcy5 + Pax6 fits striatal medium-spiny / GABAergic projection neurons. Identical to 0.5 clean cluster 12. Partially mapped onto parent cluster 14 (Lamp5+/Lypd1+ inhibitory, which split into excitatory and striatal flavors in the clean build). |

## Summary by lineage

I split the 30 clean-0.7 clusters into:

* **Oligodendrocyte lineage (4 clusters):** OPCs (13) and three mature OL states (0, 9, 19). One more OL cluster than either the parent or the 0.5 clean build; cluster 19 is the new very-low-FC OL state.
* **Astrocytes (2 clusters):** gray-matter Aqp4+ (2) and Agt+ subtype (8). Identical to the 0.5 clean build and parent. The Cxcl10 IFN hint that surfaced inside the 0.5 clean astrocyte cluster drops out of top-10 at 0.7; if a real IFN-responsive astrocyte subpopulation exists, it would need sub-clustering of 2 / 8 specifically.
* **Microglia (2 clusters, NEW SPLIT):** homeostatic (12) and DAM (6). Neither the parent nor the 0.5 clean build resolved this split; raising Louvain to 0.7 is sufficient. Most actionable gain.
* **T cells (1 cluster, RECOVERED):** cluster 28 carries a cleaner cytotoxic-T-cell signature than the parent (Ccl5 / Klrb1c stronger). Top biological priority for spatial inspection near plaques.
* **Ependymal + choroid plexus (2 clusters):** 20 and 22. CP keeps the MHC-II call that surfaced at 0.5.
* **Vascular (3 clusters):** endothelial (4), pericyte (18), VLMC (7). VLMC top-10 is finally clean now that the T-cell contamination is pulled out into cluster 28.
* **Excitatory neurons (8 clusters):** Fezf2+ deep-layer (1), Cux2+/Lamp5+ upper-layer (5), Vglut2+/Rorb+ (10), Calb2+/Tac1+/Nts+/Gck+ mixed (3), Tshz2+ recovered (24), Wif1+/Calb2+/Syt6+/Sncg+ new (25), Nts+/Fezf2+/Pou3f1+/Sorl1+ deep-layer disease-associated new (29), Rspo1+/Tcap+/IEG-rich (15). Eight is the most excitatory clusters of any build so far.
* **GABAergic interneurons (4 clusters, RECOVERED THREE-WAY SPLIT):** Sst+/Lhx6+ (17), Vip+/Crh+ (21), Pvalb+/Trh+ (26), and the striatal Tac1+/Adcy5+/Pax6+ (11). Matches the parent count and identity; the 0.5 clean build had collapsed the Sst / Vip into a single mixed cluster.
* **Disease-associated / stressed (3 clusters):** stressed neurons Pou3f1+/Lpl+ (14), reactive Tdo2+/Gpnmb+/Serpina3n+ (16), and the new Nts+/Fezf2+/Pou3f1+/Sorl1+ deep-layer cluster (29, listed under excitatory too).
* **Low-specificity (2 clusters):** 23 (IEG-dominated) and 27 (multi-lineage Slc32a1+ / Mbp+ doublets). The 0.5 clean build only had one such cluster (clean 20).

## Cross-mapping — parent (25) vs clean 0.5 (22) vs clean 0.7 (30)

| Cell type / population | Parent (res 0.5) | Clean 0.5 (res 0.5) | Clean 0.7 (res 0.7) | Net change |
|---|---|---|---|---|
| Mature OL high-FC | 0 | 0 | 0 | preserved |
| Mature OL low-FC twin | 9 | 9 | 9 | preserved |
| Mature OL third state (very low FC) | absent | absent | **19 (NEW)** | gained at 0.7 |
| OPCs | 12 | 13 | 13 | preserved |
| Astrocyte Aqp4+ | 5 | 4 | 2 | preserved (renumbered) |
| Astrocyte Agt+ | 8 | 8 | 8 | preserved |
| **Microglia homeostatic** | 4 (mixed) | 2 (mixed) | **12 (split)** | resolved at 0.7 |
| **Microglia DAM** | 4 (mixed) | 2 (mixed) | **6 (split)** | resolved at 0.7 |
| **T cells (cytotoxic)** | 24 | scattered in 7 | **28 (recovered)** | lost at clean-0.5, recovered at 0.7 |
| Endothelial BBB | 6 | 5 | 4 | preserved |
| Pericyte | 16 | 17 | 18 | preserved |
| VLMC + perivascular macrophage | 7 | 7 (T-cell contamination) | 7 (clean) | preserved; contamination lifted |
| Ependymal | 17 | 18 | 20 | preserved |
| Choroid plexus (with MHC-II at clean) | 18 | 19 | 22 | preserved; MHC-II from clean onward |
| Excitatory upper-layer Cux2+/Lamp5+/Calb1+ | 1 | 6 (merged with inhibitory Lamp5+) | 5 (un-merged) | clean-0.5 merger reversed at 0.7 |
| Excitatory deep-layer Fezf2+ | 2 | 1 | 1 | preserved |
| Excitatory Calb2+/Slc17a6+ | 3 | 3 + 11 (split) | 3 + 10 (still split) | parent cluster carved across 2 clean clusters |
| Excitatory Vglut2+/Rorb+ | 10 | 11 | 10 | preserved |
| Excitatory Tshz2+ deep/hippocampal | 19 | 15 (absorbed) | **24 (recovered)** | lost at clean-0.5, recovered at 0.7 |
| Excitatory Calb1+/Tdo2+ (region uncertain) | 20 | REMOVED (cell drop) | REMOVED (cell drop) | cluster cells removed by user filter |
| Excitatory Wif1+/Syt6+/Calb2+/Sncg+ | absent | absent | **25 (NEW)** | gained at 0.7 |
| Excitatory Nts+/Fezf2+/Pou3f1+/Sorl1+ deep-layer | absent | absent | **29 (NEW)** | gained at 0.7 |
| Excitatory Rspo1+/Tcap+/IEG | absent | 15 | 15 | gained at clean-0.5 |
| Sst+ GABAergic | 11 | 10 (merged) | **17 (un-merged)** | clean-0.5 collapse reversed at 0.7 |
| Vip+ GABAergic | 21 | 10 (merged) | **21 (un-merged)** | clean-0.5 collapse reversed at 0.7 |
| Pvalb+/Trh+ | 22 | 21 | 26 | preserved |
| Lamp5+/Lypd1+ inhibitory | 14 | split into 6 (excit) + 12 (striatal) | split into 5 (excit) + 11 (striatal) | parent cluster carved across 2 clean clusters at both resolutions |
| Stressed neurons Pou3f1+/Lpl+ | 13 | 14 | 14 | preserved |
| Reactive mixed Gpnmb+/Serpina3n+/Car12+ | 15 | 16 | 16 | preserved |
| Low-specificity IEG | 23 | 20 | 23 | preserved |
| Low-specificity Slc32a1+/Mbp+ (doublets) | absent | absent | **27 (NEW)** | gained at 0.7 (was likely inside parent 23) |

## Differences vs the two earlier calls (key biological observations)

**Resolution 0.7 is a strict improvement over both prior builds on every actionable population except low-specificity / doublet clusters, which is expected:**

1. **Microglia split into homeostatic and DAM.** This is the single most important gain. Cluster 6 (DAM: Cst7 / Itgax / Clec7a / Ch25h / Trem2 / Cd14 / Bcl2a1b) and cluster 12 (homeostatic: Tmem119 / P2ry12 / Selplg / Siglech / Cx3cr1) are the natural target pair for the APOE E4-vs-E2 microglia-switch DE on the 5XFAD background. Neither the parent SpaNorm nor the 0.5 clean build separated these states.
2. **T cells recovered.** Cluster 28 (Cd3d FC 7.52) is comparable to the parent's cluster 24 (Cd3d FC 7.24) and cleaner on Ccl5 / Klrb1c. The 0.5 clean build had lost this signal into the VLMC neighborhood; raising Louvain to 0.7 fully restores it.
3. **GABAergic three-way split recovered.** Sst+/Lhx6+ (17), Vip+/Crh+ (21), and Pvalb+/Trh+ (26) all sit as discrete clusters at 0.7, matching the parent's three-way split. The 0.5 clean collapse is reversed.
4. **Tshz2+ excitatory cluster recovered.** Cluster 24 maps cleanly onto parent cluster 19. The 0.5 clean build had absorbed it into the Rspo1+/Tcap+ cluster 15.
5. **Upper-layer / inhibitory un-merging.** The 0.5 clean build had merged the parent's excitatory Cux2+/Lamp5+/Calb1+ cluster (parent 1) with the inhibitory Lamp5+/Lypd1+ cluster (parent 14) into a single mixed cluster (clean-0.5 cluster 6). At 0.7 the excitatory cluster (now 5) carries Cux2 + Slc17a7 without Gad1, and the inhibitory Lamp5+/Lypd1+ cells redistribute back into the recovered Sst (17) and Vip (21) inhibitory clusters and the striatal (11) cluster.
6. **VLMC top-10 cleaned up.** With the T cells pulled into cluster 28, the VLMC cluster (7) loses the Cd3d / Itk / Nkg7 contamination that the 0.5 clean build flagged.
7. **Three NEW clusters with no parent analogue.** Cluster 19 (third OL state, very low FC), cluster 25 (Wif1+/Syt6+/Calb2+/Sncg+ excitatory), and cluster 29 (Nts+/Fezf2+/Pou3f1+/Sorl1+ deep-layer excitatory) appear only at 0.7. Cluster 29 in particular looks like a disease-associated excitatory neuron population given the AD-risk gene Sorl1 and the immediate-early gene Bhlhe40; it deserves spatial inspection near plaques alongside the DAM cluster 6.

**Things that did NOT improve at 0.7:**

- The reactive Tdo2+/Gpnmb+/Serpina3n+ cluster (16) is unchanged. With DAM now in cluster 6, this reactive cluster reads more clearly as glia-neuron doublets at plaque microenvironments rather than as a separate reactive state. Run DoubletFinder.
- The stressed Pou3f1+/Lpl+ cluster (14) is unchanged.
- The astrocyte clusters do not split into reactive / homeostatic. The Cxcl10 hint seen at 0.5 dropped out of top-10 at 0.7. A reactive astrocyte subpopulation, if present, needs targeted sub-clustering of astrocytes (2 and 8) rather than a global resolution increase.
- The Calb2+ excitatory population (parent cluster 3) is still spread across two clusters (3 and 10) in both clean builds; this carve-up may be real (Calb2 marks multiple cortical and thalamic populations) but should be confirmed.

## Verdict and recommended next steps

The clean-0.7 build is the working SpaNorm scaffold. It strictly improves on the parent build on microglia (two cell types instead of one mixed) and on the 0.5 clean build on T cells and GABAergic resolution. The cell-removal step (cluster 20 dropped) and the resolution increase (0.5 → 0.7) are both load-bearing for the final cell-type map.

I would (1) **start the E4-vs-E2 DE on cluster 6 (DAM) and cluster 12 (homeostatic microglia)** as the primary pseudobulk contrasts; (2) cross-check cluster 28 (T cells) for E4-vs-E2 abundance differences and for spatial proximity to plaques; (3) run DoubletFinder on clusters 14, 16, 23, and 27 before any DE involving those clusters; (4) sub-cluster the astrocyte pair (2 and 8) at higher resolution to test whether a reactive astrocyte subpopulation exists that the 0.7 global resolution missed; (5) confirm the new cluster 19 (third OL state) and cluster 29 (Nts+/Fezf2+/Pou3f1+/Sorl1+ deep-layer) are not dissociation artefacts by inspecting their spatial distribution; and (6) promote `data/20260516_microglia_switch_spanorm_p025_clean.qs2` (with the 0.7 cluster labels written into a dedicated metadata column) as the project-wide working object.
