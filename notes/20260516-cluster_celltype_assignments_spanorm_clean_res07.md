# Cluster cell-type assignments — SpaNorm clean (clusters 15 + 20 removed), sample.p = 0.25, Louvain 0.7

I assigned each of the 27 Louvain-0.7 clusters from the rebuilt `20260516_SpaNorm_samplep025_clean.qmd` using the top-10 FindAllMarkers per cluster (ranked by `avg_log2FC`; source: `results/20260516-find_all_markers_spanorm_p025_clean_res07.csv`). The panel is the 480-gene Xenium mBrain panel, so the cell-type vocabulary is constrained to what those probes resolve. I ordered the rows by cell-type similarity, not by cluster number, so related populations sit next to each other. The "Flagged genes" column lists top-10 markers that do not fit the called identity and that I would scrutinize as possible doublet or segmentation contamination before downstream DE.

**Scope note.** Same build as the 0.5 note (`notes/20260516-cluster_celltype_assignments_spanorm_clean_res05.md`); parent clusters 15 (mixed reactive doublets) and 20 (Calb1+/Tdo2+ regional excitatory neurons) were dropped at source from `data/20260512_microglia_switch_spanorm_p025.qs2` before re-fitting SpaNorm per sample. Both Louvain resolutions are stashed on the same Seurat object (`xen$clusters_res05` / `xen$clusters_res07`); the 0.7 labels are the project working scaffold for non-microglia DE.

## Assignments

| Cluster | Cell type | Top-10 markers supporting the call | Flagged genes (possible contamination) | Notes |
|---|---|---|---|---|
| 14 | Oligodendrocyte precursor cells (OPCs) | Pdgfra (6.39), Tnr (3.92), Plpp4 (3.46), Gpt2 (2.40), Chst11 (2.27), Dpyd (2.15), Pllp (1.95), Arsb (1.66) | Mki67 (1.81, pct.1 0.06), Gpnmb (1.72, pct.1 0.20) | Identical top-10 to the 0.5 OPC cluster (12 there); Mki67 at low pct.1 again flags a proliferative fraction. Maps onto 0.5 cluster 12. |
| 0 | Mature oligodendrocytes (high-FC core) | Mag (4.38), Aspa (3.78), Mog (3.55), Gatm (3.55), Ermn (3.43), Cldn11 (3.40), Pllp (3.23), Mbp (3.22), Plin4 (2.53), Cers2 (2.38) | — | Same myelinating-OL signature as 0.5 cluster 0; FCs are slightly lower because the third low-FC OL state (cluster 20 here) carved 5,337 balanced cells out of the high-FC pool. Maps onto 0.5 cluster 0. |
| 7 | Mature oligodendrocytes (lower-FC twin, Slide1-restricted) | Cldn11 (3.26), Ermn (3.21), Apod (3.16), Mog (3.15), Aspa (2.66), Plin4 (2.64), Gatm (2.28), Pllp (2.17), Mag (2.02), Ptgds (1.96) | — | Same identity as 0.5 cluster 9; Louvain at 0.7 did not redistribute these cells. **Slide1-restricted (19,758 cells on Slide1 vs 1,194 on Slide2; ratio 0.06:1).** White-matter sectioning effect; concentrated in 4s2_F1, 4s2M_F1, 4s2M_F2. |
| 20 | Third oligodendrocyte state, low-FC (NEW at 0.7) | Aspa (1.75), Mag (1.72), Mog (1.69), Pllp (1.57), Gatm (1.57), Cldn11 (1.47), Ermn (1.45), Plin4 (1.43), B3galt5 (1.28), Aldh3b1 (1.27) | — | Same myelin-gene set as clusters 0 and 7 but at uniformly low FCs (1.3-1.8) and balanced across slides (2,055 / 3,282). **New split at 0.7** that 0.5 had merged into the lower-FC twin (cluster 9). Possibly a third maturation / regional OL state; worth checking spatial distribution before reading too much into it. |
| 2 | Astrocytes (gray-matter, Aqp4+/Gja1+) | Plpp3 (3.33), Htra1 (3.28), Gja1 (3.18), Prodh (3.01), Aqp4 (2.91), Slc7a10 (2.70), Mt2 (2.70), Cbs (2.62), Pla2g7 (2.59), Aldh1l1 (2.53) | — | Identical top-10 to the 0.5 Aqp4+ astrocyte (4 there). Maps onto 0.5 cluster 4. |
| 8 | Astrocytes (Agt+ subtype) | Agt (4.90), Slc7a10 (3.66), Gja1 (2.88), Aqp4 (2.61), Pla2g7 (2.54), Cbs (2.49), Atp1b2 (2.48), Plpp3 (2.27), Prodh (2.21), Htra1 (2.21) | — | Same identity as 0.5 cluster 8. Maps onto 0.5 cluster 8. |
| 11 | **Homeostatic microglia (NEW at 0.7)** | Tmem119 (5.34), P2ry12 (5.12), Selplg (4.55), Siglech (4.50), Il1a (4.47), Cx3cr1 (4.43), Gpr34 (4.25), Lpcat2 (4.04), Ptgs1 (3.94), Nlrp3 (3.87) | — | Tmem119 / P2ry12 / Selplg / Siglech / Cx3cr1 anchor a clean homeostatic call. Il1a / Nlrp3 / Lpcat2 add a sterile-inflammasome flavor that is biologically plausible inside homeostatic microglia in a 5XFAD context. **This cluster did not exist at 0.5** (merged into cluster 2 with DAM). 14,246 cells (5,131 / 9,115). |
| 5  | **DAM microglia (NEW at 0.7)** | Cst7 (6.29), Ccl3 (5.89), Clec7a (5.54), Itgax (5.51), Ch25h (5.34), Lyz2 (5.21), Bcl2a1b (5.04), Ccl4 (4.97), Cd14 (4.92), Trem2 (4.77) | — | Cst7 / Ccl3 / Clec7a / Itgax / Ch25h / Trem2 / Cd14 give an unambiguous Keren-Shaul-style DAM signature; Ccl4 / Bcl2a1b add chemokine flavour. **This cluster did not exist at 0.5** (merged with homeostatic into cluster 2). **Load-bearing 0.7 gain for the E4-vs-E2 contrast.** 21,143 cells (3,657 / 17,486). The targeted microglia re-fit (`20260517_microglia_subcluster.qmd`) resolves DAM Stage-1 vs Stage-2 inside this cluster. |
| 4 | Endothelial cells (BBB) | Cldn5 (5.49), Cxcl12 (5.19), Flt1 (4.99), Emcn (4.72), Clec2d (4.47), Klf4 (4.32), Slc2a1 (3.94), Ocln (3.86), Sgms1 (3.71) | Ly6g (3.46, pct.1 0.11) | Identical signature to 0.5 cluster 5. Ly6g remains as a low-pct.1 entrant. |
| 6 | VLMC + perivascular macrophages | Tagln (7.87), Slc47a1 (6.88), Mrc1 (6.55), Lum (6.48), Col1a2 (5.69), Mgp (5.65), Pf4 (5.60), Serping1 (5.16), Serpinf1 (4.98), Osr1 (4.72) | — | Same identity as 0.5 cluster 7; T-cell markers no longer redistribute here (T cells are cluster 26). |
| 18 | Pericytes / mural cells | Kcnj8 (7.83), Atp13a5 (7.33), Vtn (6.39), Bgn (3.12), Ifitm1 (3.05), Stard8 (2.77), Mgp (2.44), Axl (2.35), Ifitm2 (1.89) | Emcn (2.01, pct.1 0.44) | Same identity as 0.5 cluster 15. |
| 19 | Ependymal cells (with mild reactive / CP contamination) | Ccdc153 (9.04), Tmem212 (8.94), Dnah11 (7.12), Plin5 (5.45), Vim (3.78) | Lgals3 (4.67), C3 (4.44), Ttr (3.53), Sncg (3.51), Ifitm3 (3.17) | Identical to 0.5 cluster 17. Roughly slide-balanced (3,091 / 2,717). |
| 21 | Choroid plexus epithelial cells (with MHC-II) | Car12 (7.75), Ltc4s (5.56), Adh1 (4.83), Tgfbi (4.71), H2-Eb1 (4.03), Ttr (3.90), H2-Ab1 (3.46), Vcam1 (3.27) | Htr2c (3.97, pct.1 0.98), Dnah11 (3.58, pct.1 0.56) | Same identity as 0.5 cluster 18. MHC-II pair (H2-Eb1 / H2-Ab1) preserved. |
| 1 | Excitatory neurons, deep-layer (Fezf2+) | Fezf2 (3.93), Syt6 (2.59), Slc17a7 (2.19), Oprk1 (2.15), Egr1 (2.01), Hmgcr (1.99), Fdft1 (1.92), Fut9 (1.90), Pfkp (1.79), Dgat2 (1.66) | Dgat2 | Same identity as 0.5 cluster 1; cleaner top-10 (Fut9 enters, Ldha drops). |
| 3 | Excitatory neurons, Calb2+/Nts+/Tac1+ mixed | Gck (4.19), Nts (3.45), Calb2 (3.11), Tac1 (2.64), Slc17a6 (2.37), Tox2 (2.24), Gaa (2.24), Htr2c (2.21), Oprk1 (2.21), Tmx4 (1.77) | Gck | Same identity as 0.5 cluster 3; Sst dropped out of top-10. |
| 9 | Excitatory neurons, upper-layer (Cux2+/Lamp5+) | Lamp5 (3.63), Pparg (3.00), Stard8 (2.83), Otof (2.75), Calb1 (2.68), Cux2 (2.46), Junb (2.27), Fosb (2.10), Chst1 (1.92), Grin2b (1.89) | Pparg, Junb / Fosb | Partial split of 0.5 cluster 6: 0.7 separates the Cux2+/Lamp5+ proper population here from the Rspo1+/Tcap+/IEG-rich tail (cluster 17). Calb1 enters top-10 in the cleaner upper-layer call. |
| 17 | Excitatory neurons, Rspo1+/Tcap+/IEG-rich | Rspo1 (4.44), Tcap (4.28), Egr1 (3.58), Rorb (3.27), Whrn (3.24), Hkdc1 (3.08), Pparg (3.07), Fos (2.53), Junb (2.22), Cux2 (2.12) | Tcap, Rspo1, IEG (Egr1 / Fos / Junb) | The other half of 0.5 cluster 6; Rspo1 / Tcap / Whrn anchor a thalamic or hippocampal sub-population in an activated state. Verify spatial distribution before interpreting; the IEG load can also indicate dissociation / fixation stress. |
| 10 | Excitatory neurons, Vglut2+/Rorb+ (thalamic relay) | Slc17a6 (4.06), Plcb4 (3.42), Calb2 (2.92), Tunar (2.81), Rorb (1.77), Abca7 (1.72), Calb1 (1.68), Ldhb (1.66), Sptlc2 (1.45), Got1 (1.42) | — | Same identity as 0.5 cluster 10. |
| 22 | Excitatory neurons, Tshz2+ (deep / hippocampal) | Tshz2 (4.20), Dpp4 (3.13), Rspo1 (2.74), Pparg (2.48), Htr2c (2.44), Lamp5 (2.27), Fut9 (1.87), H2az1 (1.76), Dcx (1.69) | Aldob (1.96, pct.1 0.07) | Same identity as 0.5 cluster 19. Cleaner Tshz2 FC (4.20 vs 4.12). |
| 24 | Excitatory neurons, Wif1+/Syt6+ | Wif1 (6.08), Syt6 (5.13), Calb2 (4.34), Sncg (4.19), Tac1 (2.66), Plpp4 (1.97), Slc17a6 (1.95), St6galnac4 (1.87), Acly (1.73) | Gys2 (2.61, pct.1 0.03) | Same identity as 0.5 cluster 21. |
| 15 | Stressed excitatory neurons (Pou3f1+/Lpl+) | Pou3f1 (4.10), Grin2a (3.56), Chst1 (2.70), Grin2b (2.67), Opcml (2.63), Slc17a7 (2.18), Lpl (2.06), Pfkl (1.93), Sorl1 (1.87), Hmgcr (1.86) | **Lpl**, Pou3f1 | Same identity as 0.5 cluster 13. **Slide2-enriched (2,193 / 10,264) but driven by sample 4s2M_F3 (4,274 cells); sample-level, not slide-wide.** Filter with `AddModuleScore()` before DE. |
| 13 | Sst+/Lhx6+ GABAergic interneurons (MGE-derived) | Lhx6 (5.57), Sst (5.04), Reln (3.87), Gad1 (3.76), Ndnf (3.51), Pthlh (2.94), Pvalb (2.92), Cnr1 (2.63), Tcap (2.11), Maf (2.05) | Pvalb (2.92, pct.1 0.35) | Same identity as 0.5 cluster 11. |
| 16 | Striatal-flavored GABAergic (Tac1+/Adcy5+/Pax6+) | Ido1 (4.53), Otof (4.18), Tac1 (4.04), Adcy5 (3.16), Cd36 (3.10), Pax6 (2.53), Nts (2.53), Gad1 (2.11), Csgalnact1 (1.98), Fosb (1.87) | Cd36, Ido1 | Same identity as 0.5 cluster 14; Nts enters top-10. |
| 23 | Vip+/Crh+ GABAergic interneurons (CGE-derived) | Vip (9.93), Crh (6.99), Pthlh (5.14), Cnr1 (4.16), Gad1 (3.14), Plpp4 (2.77), Calb2 (2.55), Maf (1.97) | Tph2 (2.18, pct.1 0.07), Macc1 (2.13, pct.1 0.02) | Same identity as 0.5 cluster 20. |
| 25 | Pvalb+/Trh+ GABAergic interneurons | Trh (6.66), Pvalb (5.64), Gad1 (4.36), Spp1 (3.07), Sst (2.68), Whrn (2.52), Slc32a1 (2.47), Gpx1 (2.37), Oprk1 (2.22), Serpinf1 (2.11) | **Spp1**, Sst | Same identity as 0.5 cluster 22. Spp1 in a PV+ neuron remains the standing off-profile flag. |
| 12 | Small Lypd1+/Slc17a7+ neuron-like cluster | Lypd1 (4.61), Slc17a7 (2.22), Reln (2.11), Trh (1.99), Plpp4 (1.90), Otof (1.87), Grin2b (1.68), Dgat2 (1.65), Hmgcr (1.65), Jun (1.65) | Dgat2, Hmgcr (metabolic) | Same identity as 0.5 cluster 16; Slc17a7 enters higher in top-10, reinforcing the excitatory call over an inhibitory one. |
| 26 | T cells / cytotoxic lymphocytes | Ccl5 (7.88), Cd3d (7.61), Itk (6.49), Nkg7 (6.29), Ms4a4b (5.99), Gzmb (4.96), Klrb1c (4.85), Tnf (4.70), Cd247 (4.59), Cxcl10 (4.18) | — | Same identity, same absolute counts (94 / 440) and same ratio as 0.5 cluster 24. Prioritize for spatial inspection near plaques. |

The Slide2-private multi-lineage doublet cluster (0.5 cluster 23, Slc32a1+/Mbp+, all FCs ≤ 1.1) is **absorbed at 0.7 and does not appear as a discrete cluster.** The higher resolution removes the smallest slide-private fragment for free.

## Summary by lineage

I split the 27 clean-0.7 clusters into:

* **Oligodendrocyte lineage (4 clusters):** OPCs (14), mature-OL high-FC core (0), mature-OL lower-FC Slide1-restricted twin (7), third low-FC OL state (20, **NEW at 0.7**).
* **Astrocytes (2 clusters):** gray-matter Aqp4+ (2) and Agt+ subtype (8).
* **Microglia (2 clusters, NEW SPLIT at 0.7):** homeostatic (11: Tmem119 / P2ry12 / Selplg / Siglech / Il1a / Cx3cr1 / Lpcat2 / Ptgs1) and DAM (5: Cst7 / Ccl3 / Clec7a / Itgax / Ch25h / Lyz2 / Bcl2a1b / Cd14 / Trem2).
* **T cells (1 cluster):** cluster 26.
* **Ependymal + choroid plexus (2 clusters):** 19 and 21.
* **Vascular (3 clusters):** endothelial (4), VLMC + perivascular macrophages (6), pericytes (18).
* **Excitatory neurons (7 clusters):** deep-layer Fezf2+ (1), Calb2+/Nts+/Tac1+ mixed (3), Cux2+/Lamp5+ upper-layer (9), Rspo1+/Tcap+/IEG-rich (17), Vglut2+/Rorb+ thalamic relay (10), Tshz2+ deep / hippocampal (22), Wif1+/Syt6+ (24).
* **GABAergic interneurons (4 clusters):** Sst+/Lhx6+ MGE (13), striatal Tac1+/Adcy5+/Pax6+ (16), Vip+/Crh+ CGE (23), Pvalb+/Trh+ (25).
* **Small / atypical (2 clusters):** Pou3f1+ stressed neurons (15), Lypd1+/Slc17a7+ small neuron-like (12).

## Slide differences

Same Slide2 : Slide1 baseline of ~2:1 as in the 0.5 note. Clusters more than ~2x off baseline at 0.7:

| Cluster | Identity | Slide1 / Slide2 | Ratio | Driver |
|---|---|---|---|---|
| 7  | Mature OL lower-FC twin | 19,758 / 1,194 | 0.06:1 | Slide1-restricted; same identity as 0.5 cluster 9. Driver (regional vs transcript-yield residual) not established; see 0.5 note. |
| 0  | Mature OL high-FC | 3,409 / 42,113 | 12.35:1 | Slide2 OL pool concentrates here; not a "more white matter" effect (proportional OL fractions are 18.3 % Slide1 vs 16.5 % Slide2 at 0.5; at 0.7 the third OL state (cluster 20) carves 5,337 balanced cells off the high-FC pool, sharpening the ratio). |
| 15 | Pou3f1+ stressed neurons | 2,193 / 10,264 | 4.68:1 | Sample-driven (4s2M_F3: 4,274 cells). |
| 26 | T cells | 94 / 440 | 4.68:1 | Both slides represented, absolute counts small. |
| 5  | DAM microglia | 3,657 / 17,486 | 4.78:1 | Slide2 captures more amyloid-proximal microglia in absolute numbers, not a slide-effect on the DAM signature itself. |
| 22 | Tshz2+ deep / hippocampal | 915 / 2,596 | 2.84:1 | Regional sampling difference. |
| 12 | Lypd1+ small neuron-like | 3,502 / 10,276 | 2.93:1 | Regional sampling. |

The mature-OL Slide1 / Slide2 split (clusters 7 / 0) survives unchanged from 0.5 and is the load-bearing slide effect. The 0.7 cut absorbs the 0.5 Slide2-private doublet cluster (0.5 cluster 23) and introduces no new slide-private fragment. The DAM enrichment on Slide2 (cluster 5 at 4.78:1) is interpretable as more amyloid-proximal microglia per absolute cell in Slide2 sections; cross-validate against per-cell `AddModuleScore()` for amyloid-pathway markers before reading too much into it.

## Resolution differences vs the 0.5 cut

Two biological splits at 0.7 that 0.5 does not resolve:

1. **Microglia split into homeostatic and DAM.** 0.5 cluster 2 (35,723 cells) splits into 0.7 cluster 11 (homeostatic, 14,246 cells) and 0.7 cluster 5 (DAM, 21,143 cells). **The load-bearing 0.7 gain for the E4-vs-E2 contrast.** The custom APOE probes (`rs429358-WT-228` for E3/E4, `rs7412-ALT-226:T` for E2) enrich differentially across DAM versus homeostatic microglia at the finer microglia sub-cluster resolution (`20260517_microglia_subcluster.qmd`).
2. **Third low-FC oligodendrocyte state.** 0.7 cluster 20 (5,337 cells) separates from the mature-OL pool; 0.5 had merged it into cluster 9. The third state's gene set (Aspa / Mag / Mog / Gatm / Pllp / Cldn11) is identical to clusters 0 and 7 but at FCs 1.3-1.8. Worth checking the spatial distribution (regional vs maturation gradient) before further interpretation.

One small split at 0.7 that 0.5 does not resolve:

3. **Upper-layer Cux2+/Lamp5+ vs Rspo1+/Tcap+/IEG-rich.** 0.5 cluster 6 splits into 0.7 cluster 9 (clean upper-layer Cux2+/Lamp5+/Calb1+) and 0.7 cluster 17 (Rspo1+/Tcap+/IEG-rich activated state). The IEG load suggests this is an activated / dissociation-stressed sub-population of the upper-layer pool; spatial inspection required.

One small cluster from 0.5 that does not survive at 0.7:

4. **Slide2-private multi-lineage doublet cluster** (0.5 cluster 23, Slc32a1+/Mbp+, 1,000 cells, all FCs ≤ 1.1). Absorbed into other clusters at 0.7 and not a discrete island.

All other identities map one-to-one between the two resolutions; the mapping is given in the cluster-composition table inside the qmd (`#### Cluster composition (resolution 0.7)`) and in the per-row "Maps onto 0.5 cluster X" notes above.

## Verdict and recommended cleanup before downstream DE

Louvain 0.7 is the working scaffold for the global E4-vs-E2 DE on non-microglia cell types. The microglia DAM / homeostatic split (clusters 5 / 11) is the primary 0.7 gain; targeted microglia DE flows through `20260517_microglia_subcluster.qmd` rather than directly off these two clusters.

I would (1) drop the residual contamination flag on cluster 15 (Pou3f1+/Lpl+ stressed neurons) using `AddModuleScore()` filtering before DE, (2) confirm cluster 17 (Rspo1+/Tcap+/IEG-rich) is not a dissociation artefact by spatial inspection, (3) sub-cluster astrocytes (2, 8) at higher resolution to test for an IFN-responsive astrocyte sub-population (Cxcl10 sat in the 0.5 cluster-4 tail at pct.1 0.057), (4) spatially inspect cluster 26 (T cells) relative to amyloid plaques, (5) cross-validate the custom APOE probes (`rs429358-WT-228`, `rs7412-ALT-226:T`) against the genotype labels per cell, and (6) check cluster 20 (third OL state) spatially to decide whether it is a regional or maturation sub-state of the mature-OL pool.
