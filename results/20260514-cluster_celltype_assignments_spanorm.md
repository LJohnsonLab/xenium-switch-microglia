# Cluster cell-type assignments — SpaNorm sample.p = 0.25, resolution 0.5

I assigned each of the 25 Louvain clusters from `20260513_SpaNorm_samplep025.qmd` using the top 10 FindAllMarkers per cluster (ranked by `avg_log2FC`; source: `results/20260514-find_all_markers_spanorm_p025.csv`). The panel is the 480-gene Xenium mBrain panel, so the cell-type vocabulary is constrained to what those probes resolve. I ordered the rows by cell-type similarity, not by cluster number, so related populations sit next to each other. The "Flagged genes" column lists top-10 markers that do not fit the called identity and that I would scrutinize as possible doublet or segmentation contamination before downstream DE.

This SpaNorm build produces five fewer clusters than the SCT no-regression build (25 vs 30) and reveals a discrete T-cell cluster (24) that the SCT clustering did not separate.

## Assignments

| Cluster | Cell type | Top-10 markers supporting the call | Flagged genes (possible contamination) | Notes |
|---|---|---|---|---|
| 12 | Oligodendrocyte precursor cells (OPCs) | Pdgfra, Tnr, Pllp, Chst11, Dpyd, Plpp4 | Abca1, Lpcat2, Rapgef3 | Pdgfra at FC 6.47 anchors the call; Pllp and Tnr are OPC-consistent. The lipid genes are likely OPC-intrinsic, not contamination. |
| 0  | Mature oligodendrocytes (high-FC core) | Mog, Mag, Aspa, Ermn, Cldn11, Apod, Pllp, Ptgds, Folh1, Cers2 | — | Textbook myelinating oligodendrocyte signature; all ten markers fit. The SCT cluster 0 reproduces here cleanly. |
| 9  | Mature oligodendrocytes (lower-FC twin) | Mog, Ermn, Cldn11, Apod, Aspa, Pllp, Mag, Cers2, Phgdh, Ptgds | — | Same gene set as cluster 0 but FCs collapse from ~3.5 to ~2.3. I read this as a second OL state (regional or maturation) rather than a different lineage. SpaNorm absorbed the SCT cluster 26 (Calb2-contaminated OLs) into this clean cluster — the Calb2 contamination flag is no longer present. |
| 5  | Astrocytes (gray-matter, Aqp4/Gja1 high) | Aldh1l1, Gja1, Aqp4, Slc7a10, Plpp3, Mertk, Htra1 | Rapgef3, Prodh, Cbs | Canonical astrocyte. Rapgef3 is the new top-10 entrant under SpaNorm vs SCT cluster 4; it is consistent with astrocyte expression. |
| 8  | Astrocytes (Agt+ subtype) | Agt, Slc7a10, Plpp3, Cbs, Prodh, Gja1, Aldh1l1, Aqp4, Htra1 | Folh1 (oligo-leaning) | Agt at FC 5.04 distinguishes from cluster 5. Folh1 is also expressed in astrocytes but more strongly in OLs. |
| 4  | Microglia (homeostatic + DAM mix) | Cx3cr1, P2ry12, Trem2, C1qa, C1qc, Gpr34, Ly86, Bcl2a1b, Cd68, Hexb | — | All ten markers fit microglia. Co-occurrence of P2ry12 (homeostatic) with Trem2/Cd68 (DAM) suggests this cluster aggregates both states; sub-cluster before E4 vs E2 DE. |
| 24 | T cells / cytotoxic lymphocytes | Cd3d, Cd247, Itk, Nkg7, Gzmb, Ccl5, Ms4a4b, Tnf, Plac8, Clec2d | — | This is the most striking new call vs SCT. Cd3d at FC 7.24 and Gzmb/Nkg7 at FC ~5–6 give an unambiguous cytotoxic-T-cell signature. SpaNorm separates these cells cleanly; SCT had collapsed them into the IFN-reactive astrocyte cluster (SCT cluster 28). I would prioritize this cluster for spatial inspection near plaques. |
| 17 | Ependymal cells | Ccdc153, Tmem212, Dnah11, Vim | Ocln, Slc2a1, Ifitm3, Fas, Sncg, St3gal6 | Motile-cilia genes (FC ~7–9) are unambiguous. The vascular-flavored genes likely reflect neighboring vessels at the ventricle wall rather than true ependymal expression. |
| 18 | Choroid plexus epithelial cells | Ttr, Car12, Ltc4s, Maob, Maf, Ocln | Slc7a10, Htr2c, Ptgds, Whrn | Ttr at FC 4.00 is the gold-standard CP marker. Slc7a10 and Htr2c lean astrocyte/neuron and may indicate adjacent-tissue contamination at the CP interface. |
| 6  | Endothelial cells (BBB) | Cldn5, Flt1, Ocln, Emcn, Cxcl12, Klf4, Ifitm3 | Clec2d, Vim, St3gal6 | Tight-junction + BBB transporter signature. Notably, SpaNorm did not produce an Hbb-bs+ contamination cluster — the SCT cluster 27 RBC carry-over is gone here. |
| 16 | Pericytes / mural cells | Atp13a5, Kcnj8, Vtn, Bgn, Mgp, Stard8, Axl | Cldn5, Flt1, Ifitm3 | Atp13a5 (FC 7.41) and Kcnj8 (FC 7.87) are pericyte-specific. Endothelial markers at low pct.1 reflect physical proximity, not pericyte expression. |
| 7  | Vascular leptomeningeal cells / fibroblasts (VLMCs) | Vtn, Bgn, Col1a2, Serping1, Serpinf1, Mgp, Vim, Ifitm2, Ifitm3 | Apod | Collagen + matrix signature. Apod is the only marker leaning oligodendrocyte; otherwise a clean VLMC call. |
| 1  | Excitatory neurons, upper-layer cortex (Cux2+/Lamp5+/Calb1+) | Cux2, Lamp5, Calb1, Grin2a, Grin2b, Cnr1, Egr1, Hs2st1 | Tnr, Dgat2 | Cux2 anchors L2/3. Lamp5 here is the glutamatergic upper-layer flavor (Gad1/Slc32a1 are absent). Calb1 was not in the top-10 of the SCT counterpart (cluster 2). |
| 2  | Excitatory neurons, deep-layer cortex (Fezf2+) | Fezf2, Grin2a, Grin2b, Cnr1, Egr1, Fut9, Hmgcr, Opcml | Dgat2, Tnr | Fezf2 at FC 3.82 marks L5/6. The lipid/glycolysis genes (Dgat2, Hmgcr) are metabolic, not strictly contamination. |
| 10 | Excitatory neurons, Vglut2+/Rorb+ (thalamic relay or L4-like) | Slc17a6, Plcb4, Rorb, Calb1, Grin2a, Grin2b, Tunar, Tshz2 | Abca7, Tox2 | Rorb on a Vglut2+ background fits thalamic relay neurons more than L4 cortex. |
| 19 | Excitatory neurons, Tshz2+ (Vglut1/Vglut2+, deep/hippocampal) | Tshz2, Slc17a6, Slc17a7, Htr2c, Lamp5, Fut9, Pfkp, Got1 | Pparg, Dcx | Co-expression of Slc17a6 and Slc17a7 is unusual in cortex; suggests a hippocampal or transitional pyramidal population. Pparg in neurons remains atypical and warrants verification. |
| 20 | Excitatory neurons, Calb1+/Tdo2+ (region uncertain) | Calb1, Tdo2, Grin2a, Grin2b, Jun | Car12, Dgat2, Tnr, B3galt5, Ocln | Tdo2 at FC 4.85 is unusual outside hepatocytes and habenula/specific neurons. Calb1 + Grin2a/b argues for a glutamatergic neuron call. Car12 flags possible CP contamination from adjacent tissue. |
| 3  | Excitatory neurons, Calb2+/Slc17a6+ (mixed with GABAergic admixture) | Calb2, Slc17a6, Plcb4, Htr2c, Tox2, Cnr1 | Gad1, Dcx, Gaa, Whrn | Calb2 + Slc17a6 fit a thalamic or hippocampal CA2 glutamatergic population, but Gad1 at FC 1.37 suggests partial doublet contamination with GABAergic neighbors. |
| 11 | Sst+ GABAergic interneurons (MGE-derived) | Sst, Lhx6, Gad1, Slc32a1, Reln, Cnr1, Maf | Tox2, Dcx, Grin2b | Sst at FC 5.10 and Lhx6 at 5.50 are diagnostic. Maps directly onto SCT cluster 17. |
| 22 | Pvalb+ GABAergic interneurons (Trh+ minority) | Pvalb, Gad1, Slc32a1, Trh, Sst | Whrn, Spp1, Serpinf1, Got1, Lpgat1 | Pvalb at FC 5.74 identifies fast-spiking PV+ interneurons. Spp1 remains the strongest off-profile flag — it is a DAM/microglia marker and should not be in a PV+ neuron. |
| 21 | VIP+ GABAergic interneurons (CGE-derived) | Vip, Crh, Gad1, Slc32a1, Cnr1, Calb2, Maf, Pthlh | Htr2c, Plpp4 | Vip at FC 9.82 and Crh at 7.06 anchor the call; canonical VIP+ CGE interneuron. |
| 14 | GABAergic interneurons, Lamp5+/Lypd1+ | Gad1, Lamp5, Lypd1, Otof, Adcy5 | Grin2b, Dcx, Dgat2, Grin2a, Tnr | Gad1 + Lamp5 + Lypd1 mark the inhibitory Lamp5 subclass. Adcy5 may indicate a striatal admixture. |
| 13 | Disease-associated / stressed neurons (mixed reactive + neuronal) | Cnr1, Grin2a, Grin2b, Opcml | **Lpl**, **Serpina3n**, Pou3f1, Tnr, Dgat2, Hmgcr | The neuronal core (Grin2a/b, Opcml, Cnr1) is contaminated by reactive/lipid markers Lpl and Serpina3n. Pou3f1 in this context most likely tags a stressed neuronal state rather than an OL-lineage cell. Strong candidate for DoubletFinder before DE. |
| 15 | Mixed reactive / disease-associated (likely doublets, DAM-like + neuronal) | Serpina3n, Gpnmb, Car12 | Grin2a, Grin2b, Dgat2, Tnr, Cat, St6galnac3, Gad1 | The reactive triad (Serpina3n + Gpnmb + Car12) suggests reactive astrocyte or DAM-like signal, but the neuronal markers Grin2a/Grin2b and the inhibitory Gad1 in the same cluster point to doublets across a plaque microenvironment. Maps onto SCT cluster 12. |
| 23 | Low-specificity / transitional state | Slc32a1, Mbp, Gatm, Cd81 (all FC ≤ 1.11) | Aldh9a1, Adh5, Agt, Acadl, Hmgcl, Degs1 | All FCs are ≤ 1.11 — this cluster is barely distinguishable. Slc32a1 + Mbp + Agt in the same top-10 list points to multi-lineage doublets. Merge or re-resolve. |

## Summary by lineage

I split the 25 clusters into:

* **Oligodendrocyte lineage (3 clusters):** OPCs (12) and two mature OL clusters (0, 9). SpaNorm collapsed the SCT three-way OL split (0, 15, 26) into two, and the Calb2-contamination flag in SCT cluster 26 disappears.
* **Astrocytes (2 clusters):** gray-matter Aqp4+ (5) and Agt+ subtype (8). SpaNorm did not isolate a discrete IFN-responsive reactive-astrocyte cluster — the SCT cluster 28 (Cxcl10+/Ifit1+/Aqp4+) signature is absent here. Some of that signal may have moved into the T-cell cluster (24, which carries Ccl5 and Tnf) or been redistributed into the reactive-mixed clusters (13, 15).
* **Microglia (1 cluster):** cluster 4, still mixing homeostatic and DAM markers. Sub-cluster before DE.
* **T cells (1 cluster, NEW):** cluster 24 — the most consequential SpaNorm-vs-SCT difference. Cd3d/Gzmb/Nkg7/Itk is a clean cytotoxic-T-cell call. I would inspect this cluster's spatial location relative to amyloid plaques.
* **Ependymal + choroid plexus (2 clusters):** 17 and 18.
* **Vascular (3 clusters):** endothelial (6), pericyte (16), VLMC (7). SpaNorm cleaned up the Hbb-bs+ RBC contamination cluster that SCT carried (SCT cluster 27).
* **Excitatory neurons (6 clusters):** L2/3 (1), L5/6 (2), Vglut2+/Rorb+ (10), Tshz2+ deep/hippocampal (19), Calb1+/Tdo2+ regional (20), and Calb2+ mixed (3).
* **GABAergic interneurons (4 clusters):** Sst (11), Pvalb (22), Vip (21), Lamp5+/Lypd1+ (14). SpaNorm collapsed the SCT seven-way GABAergic split into four — the Ndnf+/Reln+ neurogliaform cluster (SCT 22), the Tac1+/Lhx6+ cluster (SCT 13), and the Nts+ neurons (SCT 21) no longer separate at this resolution.
* **Mixed / disease-associated (2 clusters):** 13 and 15. Likely doublets at plaque microenvironments.
* **Low-specificity (1 cluster):** 23.

## Differences vs the SCT no-regression call (key biological observations)

1. **T cells emerge as a discrete cluster.** SpaNorm cluster 24 carries an unambiguous cytotoxic-T-cell signature (Cd3d FC 7.24, Gzmb FC 4.97, Nkg7 FC 6.21). SCT did not separate these cells. This is the most actionable new finding.
2. **The IFN-responsive reactive-astrocyte cluster is gone.** The SCT cluster 28 signature (Cxcl10, Ifit1/2, Aqp4, Ccl12) does not appear as a top-10 driver in any SpaNorm cluster. The reactive astrocyte signal may have redistributed across clusters 13 and 15 (which carry Serpina3n) or collapsed into the broader astrocyte clusters 5 and 8.
3. **Two contamination flags vanish.** SpaNorm produces no Hbb-bs+ RBC cluster (SCT 27) and no Calb2-contaminated OL cluster (SCT 26).
4. **Fewer GABAergic subclasses are resolved.** SpaNorm produces four GABAergic clusters vs SCT's seven; the Ndnf+, Tac1+, and Nts+ subpopulations no longer separate at resolution 0.5.

## Recommended cleanup before downstream DE

I would (1) keep cluster 24 (T cells) intact and inspect its spatial distribution near plaques, (2) sub-cluster microglia (4) to separate homeostatic from DAM, (3) run DoubletFinder on clusters 13, 15, and 23, (4) re-evaluate cluster 23 at a lower Louvain resolution or merge it into the nearest neighbor, (5) check whether the reactive-astrocyte signal that SCT isolated in cluster 28 can be recovered by sub-clustering the SpaNorm astrocyte clusters (5, 8) at higher resolution, and (6) keep the SpaNorm clustering as the primary scaffold for the E4-vs-E2 comparison given the cleaner vascular and oligodendrocyte calls, while cross-validating any DAA-relevant findings against the SCT cluster 28 cells.
