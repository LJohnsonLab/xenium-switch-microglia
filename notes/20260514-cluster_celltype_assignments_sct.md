# Cluster cell-type assignments — SCT no-regression, resolution 0.5

I assigned each of the 30 Louvain clusters from `20260512_SCT_no_regression.qmd` using the top 10 FindAllMarkers per cluster (ranked by `avg_log2FC`; source: `results/20260514-find_all_markers_cst_noRegressed.csv`). The panel is the 480-gene Xenium mBrain panel, so the cell-type vocabulary is constrained to what those probes can resolve. I ordered the rows by cell-type similarity, not by cluster number, so related populations sit next to each other. The "Flagged genes" column lists top-10 markers that do not fit the called identity and that I would scrutinize as possible doublet/segmentation contamination before downstream DE.

## Assignments

| Cluster | Cell type | Top-10 markers supporting the call | Flagged genes (possible contamination) | Notes |
|---|---|---|---|---|
| 10 | Oligodendrocyte precursor cells (OPCs) | Pdgfra, Tnr, Pllp, Chst11, Dpyd, Plpp4 | Abca1, Lpcat2 | Pdgfra at FC 5.88 anchors the call; Pllp/Tnr are OPC-consistent. Lipid genes are likely OPC-intrinsic, not contamination. |
| 0  | Mature oligodendrocytes (high-FC core) | Mog, Mag, Aspa, Ermn, Cldn11, Pllp, Apod, Ptgds, Folh1, Cers2 | — | Textbook myelinating oligodendrocyte signature; all ten markers fit. |
| 15 | Mature oligodendrocytes (lower-FC twin of cluster 0) | Mog, Mag, Ermn, Aspa, Cldn11, Pllp, Apod, Ptgds, Cers2, Folh1 | — | Same gene set as cluster 0 but FCs collapse from ~4.5 to ~1.8. I read this as a second OL state (regional or maturation) rather than a different lineage. |
| 26 | Mature oligodendrocytes with neuronal contamination | Mog, Mag, Aspa, Ermn, Cldn11, Pllp, Apod, Ptgds, Folh1 | **Calb2 (FC 3.06)** | Nine of ten markers are myelin. Calb2 at FC 3.06 is neuronal (calretinin). Likely doublets with Calb2+ neurons or partial segmentation across an OL/interneuron boundary. |
| 4  | Astrocytes (gray-matter, Aqp4/Gja1 high) | Aldh1l1, Gja1, Aqp4, Slc7a10, Plpp3, Mertk, Htra1 | Vcam1 (vascular/reactive), Prodh, Cbs | Canonical astrocyte. Vcam1 may flag reactive astrocytes or endothelial bleed-through; check with cluster 5 proximity. |
| 7  | Astrocytes (Agt+ subtype) | Agt, Slc7a10, Plpp3, Cbs, Prodh, Gja1, Aldh1l1, Aqp4, Htra1 | Folh1 (oligo-leaning) | Agt at FC 4.50 distinguishes from cluster 4. Folh1 is also expressed in astrocytes but more strongly in OLs. |
| 28 | Reactive / interferon-responsive astrocytes (DAA-like) | Aqp4, Serpina3n, Serping1, Ifitm3, Cxcl10, Ifit1, Ifit2, Ccl12, Ccl5, Lgals3bp | — | Strong IFN-stim signature on an Aqp4+ background. Highly relevant to the 5XFAD background; expect enrichment near plaques. |
| 3  | Microglia (homeostatic + DAM mix) | Cx3cr1, P2ry12, Trem2, C1qa, C1qc, Gpr34, Ly86, Bcl2a1b, Cd68, Hexb | — | All ten markers fit microglia. Co-occurrence of P2ry12 (homeostatic) with Trem2/Cd68 (DAM) suggests this cluster aggregates both states; sub-cluster before DE. |
| 18 | Ependymal cells | Ccdc153, Tmem212, Dnah11, Vim | Ocln, Slc2a1, Ifitm3, Fas, Sncg, St3gal6 | Motile-cilia genes (FC ~7–8) are unambiguous. The vascular-flavored genes (Ocln, Slc2a1, Ifitm3) likely reflect neighboring vessels at the ventricle wall rather than true ependymal expression. |
| 20 | Choroid plexus epithelial cells | Ttr, Car12, Ltc4s, Maob, Maf, Hspg2, Ocln, Ptgds | Slc7a10, Htr2c | Ttr is the gold-standard CP marker. Slc7a10/Htr2c are astrocyte/neuron-flavored and may indicate adjacent-tissue contamination at the CP interface. |
| 5  | Endothelial cells (BBB) | Cldn5, Flt1, Ocln, Emcn, Cxcl12, Klf4, Slc2a1, Ifitm3 | Clec2d, Vim | Tight-junction + BBB transporter signature. Vim is broadly mesenchymal so not strictly contamination. |
| 27 | Endothelial cells with red-blood-cell contamination | Flt1, Cldn5, Cxcl12, Emcn, Ocln, Slc2a1, Ifitm3, Fas, Clec2d | **Hbb-bs (FC 8.22)** | Endothelial identity is clear, but Hbb-bs at top rank flags RBC carry-over inside the cell mask. Filter or merge with cluster 5 after removing Hbb-bs+ cells. |
| 14 | Pericytes / mural cells | Atp13a5, Kcnj8, Vtn, Bgn, Mgp, Stard8, Axl | Cldn5, Flt1, Ifitm3 | Atp13a5 + Kcnj8 are pericyte-specific. Endothelial markers at low pct suggest physical proximity, not pericyte expression. |
| 6  | Vascular leptomeningeal cells / fibroblasts (VLMCs) | Vtn, Bgn, Col1a2, Serping1, Serpinf1, Mgp, Vim, Ifitm2, Ifitm3 | Apod | Collagen + matrix signature. Apod is the only marker that leans oligodendrocyte; otherwise a clean VLMC call. |
| 2  | Excitatory neurons, upper-layer cortex (Cux2+/Lamp5+) | Cux2, Lamp5, Grin2a, Grin2b, Cnr1, Bhlhe40, Egr1, Hs2st1 | Pparg, Stard8 | Cux2 anchors L2/3. Lamp5 here is the glutamatergic upper-layer flavor (not the interneuron one — Slc32a1/Gad1 are absent). Pparg in neurons is unusual; check expression specificity. |
| 1  | Excitatory neurons, deep-layer cortex (Fezf2+) | Fezf2, Grin2a, Grin2b, Cnr1, Tnr, Opcml, Egr1 | Dgat2, Hmgcr, Pfkl | Fezf2 marks L5/6. Tnr is perineuronal-net consistent. The lipid/glycolysis genes (Dgat2, Hmgcr, Pfkl) are metabolic, not strictly contamination. |
| 16 | Excitatory neurons, deep-layer / subplate (Fezf2+/Tshz2+) | Fezf2, Tshz2, Grin2a, Grin2b, Cnr1, Lypd1, Htr2c, Otof | Tnr, Dgat2 | Similar to cluster 1 but with Tshz2+Lypd1 shift. Otof is unusual outside hair cells; expression may be panel-specific. |
| 19 | Excitatory neurons, deep-layer (Vglut2+/Tshz2+) | Slc17a6, Tshz2, Fezf2, Grin2a, Grin2b, Htr2c, Lypd1, Fut9 | Dcx, Dgat2 | Slc17a6 swaps in for Slc17a7, suggesting these are thalamic-projecting or non-cortical glutamatergic neurons. Dcx is a juvenile-neuron marker — flag for migrating cells or technical. |
| 11 | Excitatory neurons, Vglut2+/Rorb+ (thalamic relay or L4-like) | Slc17a6, Rorb, Plcb4, Tunar, Grin2a, Grin2b, Tshz2 | Abca7, Ogdhl | Rorb on a Vglut2+ background fits thalamic relay neurons more than L4 cortex. Abca7 (lipid transport) and Ogdhl flagged as off-profile. |
| 9  | Excitatory neurons, Vglut2+/Calb2+ (thalamic/hippocampal) | Slc17a6, Calb2, Calb1, Plcb4 | Tox2, Htr2c, Gaa, Dcx, Ctps2, Gmps | Calb2+Slc17a6 patterns onto specific thalamic nuclei or hippocampal CA2. Several flagged genes are metabolic (Gaa, Ctps2, Gmps) and may reflect a high-activity neuronal state. |
| 29 | Excitatory neurons, Vglut2+/Pvalb+ (thalamic / deep cerebellar nuclei type) | Slc17a6, Pvalb, Calb2, Spp1, Plpp4 | Got1, Gls, Ldhb, Sdha, Adamts2 | Pvalb on a Vglut2+ background is unusual and points away from cortical PV interneurons toward thalamic/cerebellar PV+ glutamatergic neurons. The metabolic enzyme cluster (Gls/Got1/Ldhb/Sdha) sits at ~100% expression in both groups; these are likely housekeeping rather than markers. |
| 21 | Neurotensin+ neurons (thalamic/amygdaloid/hypothalamic) | Nts, Tac1, Syt6, Calb2, Cnr1, Htr2c, Gad1 | Tox2, Gaa, Gmps | Nts at FC 7.11 is the defining call. Mixed Gad1+Tac1 suggests this is not a clean GABAergic or glutamatergic population; could be a regional cluster (CeA, thalamic midline). |
| 17 | Sst+ GABAergic interneurons (MGE-derived) | Sst, Lhx6, Gad1, Slc32a1, Reln, Cnr1, Calb1 | Tox2, Dcx, Grin2b | Sst at FC 6.77 and Lhx6 at 4.65 are diagnostic. Grin2b is broadly expressed in interneurons too. |
| 24 | Pvalb+ GABAergic interneurons (with Trh+ minority) | Pvalb, Gad1, Slc32a1, Trh, Sst | Whrn, Spp1, Plcb4, Serpinf1, Got1 | Pvalb identifies these as fast-spiking cortical/striatal PV+ INs. Spp1 is the strongest contamination flag — it is a DAM/microglia marker and should not be in a PV+ neuron. |
| 13 | GABAergic interneurons, Tac1+/Lhx6+ (MGE-derived) | Tac1, Lhx6, Gad1, Calb1, Calb2, Htr2c | Tox2, Dcx, Grin2b, Gaa | Tac1 here likely marks MGE-derived interneurons rather than striatal D1 MSNs because of co-expressed Lhx6. |
| 22 | CGE-derived interneurons (Ndnf+/Reln+ neurogliaform) | Ndnf, Reln, Lamp5, Cnr1, Maf, Gad1, Slc32a1, Wif1 | Tox2, Gaa | Ndnf+Reln+Lamp5 is the neurogliaform signature. Lamp5 in this cluster is the interneuron flavor (paired with Gad1/Slc32a1), opposite to cluster 2. |
| 25 | VIP+ GABAergic interneurons (CGE-derived) | Vip, Crh, Gad1, Slc32a1, Cnr1, Calb2, Maf, Pthlh | Htr2c, Grin2b | Vip at FC 9.53 and Crh at 6.81 anchor this; canonical VIP+ CGE interneuron. |
| 8  | GABAergic interneurons, Lamp5+/Lypd1+ | Gad1, Slc32a1, Lamp5, Lypd1, Otof | Tox2, Dcx, Htr2c, Grin2b, Gaa | Gad1+Slc32a1+Lamp5 = inhibitory Lamp5 subclass. Otof is unusual but reported in cortical interneurons. |
| 23 | Low-specificity GABAergic / transitional state | Gad1 (FC 0.41) | Whrn, Tox2, St3gal6, Tgfbr1, Gusb, Plpp3, Pla2g12a, B3galt5, Acsl1 | All FCs are ≤0.8 — this cluster is barely distinguishable. I would either merge with a neighboring GABAergic cluster or treat as transitional/poorly resolved at this resolution. |
| 12 | Mixed reactive / disease-associated (likely doublets) | Serpina3n, Gpnmb, Car12 | Grin2a, Grin2b, Dgat2, Tnr, Jun, St6galnac3, Gad1 | The reactive triad (Serpina3n + Gpnmb + Car12) suggests reactive astrocyte or DAM-like signal, but the neuronal markers Grin2a/Grin2b and the inhibitory Gad1 in the same cluster point to doublets across a plaque microenvironment. Flag for careful re-clustering and DoubletFinder. |

## Summary by lineage

I split the 30 clusters into:

* **Oligodendrocyte lineage (4 clusters):** OPCs (10) and three mature OL clusters (0, 15, 26). Cluster 26 carries Calb2 contamination.
* **Astrocytes (3 clusters):** gray-matter Aqp4+ (4), Agt+ subtype (7), and IFN-responsive reactive (28). The reactive astrocyte cluster is the strongest disease-state candidate.
* **Microglia (1 cluster):** cluster 3, mixing homeostatic and DAM markers — worth sub-clustering before E4 vs E2 DE.
* **Ependymal + choroid plexus (2 clusters):** 18 and 20.
* **Vascular (3 clusters):** endothelial (5, 27) and pericyte (14) and VLMC (6). Cluster 27 has Hbb-bs contamination.
* **Excitatory neurons (8 clusters):** L2/3 (2), L5/6 (1, 16, 19), thalamic-relay-like (11, 9, 29), and Nts+ (21).
* **GABAergic interneurons (7 clusters):** Sst (17), Pvalb (24), Tac1/Lhx6 (13), Lamp5/Lypd1 (8), Ndnf/Reln neurogliaform (22), VIP/Crh (25), low-specificity (23).
* **Mixed / disease-associated (1 cluster):** 12. Likely doublets.

## Recommended cleanup before downstream DE

I would (1) drop Hbb-bs+ cells from cluster 27, (2) run DoubletFinder on clusters 12 and 26, (3) sub-cluster microglia (3) to separate homeostatic from DAM, (4) re-evaluate cluster 23 at a lower Louvain resolution or merge it into the nearest GABAergic neighbor, and (5) keep cluster 28 intact — it is the most biologically interesting reactive-astrocyte signature on this 5XFAD background and will be the prime target for the E4-vs-E2 comparison.
