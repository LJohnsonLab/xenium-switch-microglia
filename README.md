# microglia_switch

Xenium in situ spatial transcriptomics of mouse brain under a microglia-specific APOE E4→E2 switch on a 5XFAD amyloid background.

* All mice are 5XFAD; `4s2_Fn` are controls (microglia E4) and `4s2M_Fn` are tamoxifen-switched (microglia E2). 
* The mBrain_480g panel carries custom APOE allele probes that read the genotype per cell.

## Rendered analysis

The notebooks render to a [website published with GitHub Pages](https://ljohnsonlab.github.io/xenium-switch-microglia/) from `docs/`.   
Each notebook reports its own methods, parameters, and results.


## Data availability

Shared plotting helpers live in `R/20260514-helper_functions.R`. Differential-expression result tables are in `results/`

Raw Xenium acquisitions and the processed Seurat objects are too large for git and are not tracked here. The Seurat `.qs2` objects live on [Dropbox](https://www.dropbox.com/scl/fo/iq8rtnvko9lb1sf4h5s20/AFVC-MdhrNk5Jrb6zyFTa88?rlkey=8c063z4elbdy4vqn64v2nc6tq&dl=0):



Place the downloaded `.qs2` files under `data/` to reproduce the notebooks. See `data/README.md` for each object's provenance and `xenium_rawData/README.md` for the raw-bundle layout and sample key.

## Reproduce

Heavy compute chunks are gated `eval: false`; the notebooks reload the `data/*.qs2` checkpoints rather than recomputing.  
Key packages: Seurat v5, SpatialExperiment, SpaNorm, qs2, duckplyr, arrow (>= 24), DESeq2, EnhancedVolcano.


