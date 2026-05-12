library(Seurat)
library(tidyverse)
library(duckplyr)
library(qs2)
library(fs)
library(glue)
library(arrow)


# Project root paths
xenium_dir <- "xenium_rawData"
slide1_dir <- path(xenium_dir, "Slide1_output-XETG00118__0022474__Region_1__20250418__203103")
slide2_dir <- path(xenium_dir, "output-XETG00118__0069080__Region_1__20260417__210507")

# Sample to slide and condition lookup. Single source of truth used by both
# the per-slide ingest and the post-merge metadata join.
samples <- tribble(
  ~sample_id, ~slide,            ~condition,
  "4s2_F1",   "Slide1_0022474",  "control",
  "4s2M_F1",  "Slide1_0022474",  "switched",
  "4s2M_F2",  "Slide1_0022474",  "switched",
  "4s2_F2",   "Slide2_0069080",  "control",
  "4s2_F3",   "Slide2_0069080",  "control",
  "4s2M_F3",  "Slide2_0069080",  "switched",
  "4s2M_F4",  "Slide2_0069080",  "switched")


# Per-cell QC from transcripts.parquet: total count, distinct features,
# median QV across all transcripts assigned to that cell.
summarize_transcripts <- function(slide_dir) {
  read_parquet_duckdb(path(slide_dir, "transcripts.parquet")) |>
    summarize(nCount_raw = n(),
              nFeature_raw = n_distinct(feature_name),
              median_qv = median(qv, na.rm = TRUE),
              .by = cell_id) |>
    collect()}

# Read all per-sample cells_stats CSVs for one slide and bind them with a
# sample_id column derived from the file name.
read_sample_rois <- function(sample_ids) {
  map(sample_ids, \(s) read_csv(glue("spatial_data/{s}_cells_stats.csv"), comment = "#", show_col_types = FALSE)) |>
    set_names(sample_ids) |>
    list_rbind(names_to = "sample_id") |>
    select(sample_id, cell_id = "Cell ID", cell_area = "Area (µm^2)")}

# Attach sample, transcript QC, and slide label to a Seurat meta.data frame.
# Returns a data.frame with the original cell barcodes as rownames.
annotate_meta <- function(meta, rois, transcripts, slide_label) {
  meta |>
    rownames_to_column("cell_id") |>
    left_join(rois, by = "cell_id") |>
    left_join(transcripts, by = "cell_id") |>
    column_to_rownames("cell_id") |>
    mutate(slide = slide_label)}

# Load one slide, annotate it, and filter to ROI-assigned cells with
# nCount_Xenium > 0 and median_qv >= 20.
load_slide <- function(slide_dir, slide_label, sample_ids) {
  obj <- LoadXenium(slide_dir)
  transcripts <- summarize_transcripts(slide_dir)
  rois <- read_sample_rois(sample_ids)
  obj@meta.data <- annotate_meta(obj@meta.data, rois, transcripts, slide_label)
  subset(obj, subset = !is.na(sample_id) & nCount_Xenium > 0 & median_qv >= 20)}


##################################
## Slide1_0022474 ###############
################################
slide1 <- load_slide(slide1_dir,
                     "Slide1_0022474",
                     samples |> filter(slide == "Slide1_0022474") |> pull(sample_id))

##################################
## Slide2_0069080 ###############
################################
slide2 <- load_slide(slide2_dir,
                     "Slide2_0069080",
                     samples |> filter(slide == "Slide2_0069080") |> pull(sample_id))


###############################################
# Merge the two filtered Seurat objects #######
###############################################
xen <- merge(slide1, slide2,
             add.cell.ids = c("slide1", "slide2"),
             project = "microglia_switch")

# Attach condition from the lookup. partial_section flags 4s2M_F2 which only
# has a tissue fragment on the slide; downstream QC decides whether to drop it.
xen@meta.data <- xen@meta.data |>
  rownames_to_column("cell_id") |>
  left_join(samples |> select(sample_id, condition), by = "sample_id") |>
  mutate(partial_section = sample_id == "4s2M_F2") |>
  column_to_rownames("cell_id")

xen <- JoinLayers(xen, assay = "Xenium")


# Save the filtered Seurat object. qs2 is much faster than RDS.
# qs_save(xen, "data/20260512_microglia_switch.qs2")
