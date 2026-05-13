# SpaNorm convergence log: sample.p = 0.25

- **Source qmd:** `20260513_SpaNorm_samplep025.qmd`
- **Date:** 2026-05-13
- **Run:** 7 per-sample SpaNorm fits at the package-default sample.p = 0.25.

## Verdict

Seven fits, all converged. Every per-sample call closed at outer iter 4 with the log-likelihood flat to six significant figures across iters 3 and 4. No `did not converge` warnings, no negative-binomial degeneracy. The pipeline is healthy at `sample.p = 0.25`.

## Sample-size sanity check

At `sample.p = 0.25` the printed cell count is one quarter of each sample's total. Back-calculating:

| Sample (order in tribble) | Sampled cells | Total ≈ |
|---|---|---|
| 4s2_F1 (Slide1) | 13,750 | 55,000 |
| 4s2M_F1 (Slide1) | 14,757 | 59,028 |
| **4s2M_F2 (Slide1, partial)** | **7,298** | **29,192** |
| 4s2_F2 (Slide2) | 18,306 | 73,224 |
| 4s2_F3 (Slide2) | 18,361 | 73,444 |
| 4s2M_F3 (Slide2) | 18,452 | 73,808 |
| 4s2M_F4 (Slide2) | 16,571 | 66,284 |

Total ≈ 430k cells, matching the merged object.

## Two structural facts the log confirms

1. **4s2M_F2 is half the size of the other Slide1 sections.** That is the partial-section flag we already carry in metadata; the fit still converged on 7.3k cells so SpaNorm is not failing on it.
2. **Slide2 sections carry ~30 % more cells than Slide1 sections** (66k-74k vs 55k-59k). Consistent with the +25 % transcripts/cell and the per-cell QC delta documented in `20260512_QC.html`; Slide2 also fielded denser segmentation.

## Convergence dynamics

Iter 1 takes 6-8 inner NB updates, iter 2 settles in 2-3, iter 3 is already at the fixed point. Iter 4 is just the safety pass. The log-likelihood magnitude scales with cells (~-3.7 M for 7k cells, ~-10.5 M for 18k cells), which is what an additive per-cell likelihood should look like.

## Per-sample iter-1 inner trajectory

| Sample | Iter-1 inner steps to converge | Iter-1 final log-likelihood |
|---|---|---|
| 4s2_F1 | 7 | -6,471,322 |
| 4s2M_F1 | 8 | -7,496,818 |
| 4s2M_F2 | 8 | -3,694,349 |
| 4s2_F2 | 6 | -10,370,380 |
| 4s2_F3 | 7 | -10,495,770 |
| 4s2M_F3 | 7 | -10,249,367 |
| 4s2M_F4 | 7 | -9,324,834 |

## Bottom line

Nothing in the log suggests a bad fit. The verdict on `sample.p = 0.25` vs `sample.p = 0.05` will turn on whether the resulting UMAP looks different from the `sample.p = 0.05` panel, not on the convergence behavior.
