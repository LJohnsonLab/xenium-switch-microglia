# Xenium raw data — sample key

Source key: `4s2MxFAD Xenium_brain key_042126.pdf` (dated 2026-04-21).

Mouse brain, fresh-frozen, custom panel `mBrain_480g` (480 targets, panel design `JRZFVH`), instrument `XETG00118`.

## Biological model

All mice carry the **5XFAD** transgene (amyloid pathology background). On top of 5XFAD, mice are on an inducible APOE "switch": `APOE4s2flox/flox × Tmem119-CreERT2` (humanized APOE locus). Tamoxifen flips the floxed allele in Tmem119+ microglia from E4 to E2 while leaving E4 expression intact in all other CNS and peripheral cells.

- **`4s2M` / `MC/+`** — 5XFAD; tamoxifen-induced; microglia express E2, rest of body expresses E4. RNA-seq sample prefix `M_*` (e.g. `M_9`, `M_10`).
- **`4s2` / `+/+`** — 5XFAD; Cre-negative or vehicle littermate controls (microglia and body both express E4).

## Subdirectory ↔ slide ↔ run

Acquisition date is extracted from `run_start_time` in each `experiment.xenium` manifest.

| Local subdirectory | Slide ID | Run name (`experiment.xenium`) | Project / run label | `run_start_time` | Acquisition date |
|------------|------------|------------|------------|------------|------------|
| `Slide1_output-XETG00118__0022474__Region_1__20250418__203103/` | `0022474` | `4s2_FAD_LG` (cassette `LG_Slide1`) | APOE 4s2 × 5XFAD — Run 1 | `2025-04-18T20:30:52Z` | 2025-04-18 |
| `output-XETG00118__0069080__Region_1__20260417__210507/` | `0069080` | `20260417_AgeXMetabolism_Run7` (cassette `20260417_slide2_0069080`) | Aging_Metab — Run 7 | `2026-04-17T21:04:57Z` | 2026-04-17 |

## Samples per slide

Samples relevant to this project (boxed in red on the PDF key). Each slide also carries additional sections from unrelated experiments that are NOT part of this analysis.

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

## Original shared-drive paths (for provenance)

- **APOE 4s2 × 5XFAD Run 1** `07. Xenium Spatial Transcriptomics / 1 APOE4s2 x 5XFAD Xenium 2025-04-21 / 20250418_203052_4s2_FAD_LG-XeniumData / Slide1_output-XETG00118__0022474__Region_1__20250418__203103`
- **Aging_Metab Run 7** `07. Xenium Spatial Transcriptomics / 2 Aging_Metab_Xenium / Run7_20260417__210457__20260417_AgeXMetabolism_Run7 / output-XETG00118__0069080__Region_1__20260417__210507`

## Notes

- Both slides export `Region_1` only — sample-level segmentation has not yet been applied; sample assignment requires manual ROI definition in Xenium Explorer (or programmatic polygon masks) using the layouts on page 2 of the PDF key.
- Slide `0022474` stores its experiment manifest as `experiment.xenium_slide1.xenium`; slide `0069080` uses the standard `experiment.xenium`.
- Sample `4s2M_F2` **is not complete**, only has a fragment remaining on the slide.