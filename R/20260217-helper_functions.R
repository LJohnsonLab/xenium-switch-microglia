### ggplot decoration for the density plots
density_decor <- function(x){
    list(
        geom_density(),
        theme_minimal(),
        theme(axis.text.y = element_blank(),
              axis.text.x = element_text(size=12),
              axis.title.x = element_text(size=14))
    )
}    
######################################################################################################
############################### FIDDLE ###############################################################    
######################################################################################################

# Load libraries only if not already loaded
if (!require(ggplot2, quietly = TRUE)) library(ggplot2)
if (!require(dplyr, quietly = TRUE)) library(dplyr)
if (!require(tidyr, quietly = TRUE)) library(tidyr)
if (!require(forcats, quietly = TRUE)) library(forcats)
if (!require(purrr, quietly = TRUE)) library(purrr)
if (!require(grid, quietly = TRUE)) library(grid)
if (!require(scales, quietly = TRUE)) library(scales)

#' Create a multi-faceted violin plot with colored axis labels
#'
#' @param seurat_obj Seurat object containing the data
#' @param genes_of_interest Character vector of gene names to plot
#' @param classification_col Character string specifying the column name for classification
#' @param assay Character string specifying which assay to use (default: "RNA")
#' @param axis_text_y_size Numeric, size of y-axis text (default: 12)
#' @param axis_text_x_size Numeric, size of x-axis text (default: 6)
#' @param strip_text_size Numeric, size of facet strip text (default: 10.5)
#' @param strip_angle Numeric, angle for strip text (default: 85)
#' @param violin_scale Character, scaling method for violins (default: "width")
#' @param violin_trim Logical, whether to trim violin tails (default: TRUE)
#' @param violin_adjust Numeric, bandwidth adjustment for violins (default: 1)
#' @param facet_nrow Numeric, number of rows for facet_wrap (default: 1)
#'
#' @return Draws the plot; returns invisible NULL
#' @export
fiddle <- function(seurat_obj,
                   genes_of_interest,
                   classification_col,
                   assay = "RNA",
                   axis_text_y_size = 12,
                   axis_text_x_size = 6,
                   strip_text_size = 10.5,
                   strip_angle = 85,
                   violin_scale = "width",
                   violin_trim = TRUE,
                   violin_adjust = 1,
                   facet_nrow = 1) {

    # Validate inputs
    if (!classification_col %in% colnames(seurat_obj@meta.data)) {
        stop(paste("Classification column", classification_col, "not found in Seurat object metadata"))
    }

    # Fetch data from Seurat object
    violin_data <- FetchData(seurat_obj, c(genes_of_interest, classification_col), assay = assay)

    # Reshape data to long format
    violin_data_long <- violin_data |>
        pivot_longer(-all_of(classification_col)) |>
        mutate(seurat_clusters = fct_rev(.data[[classification_col]]),
               name = factor(name, levels = genes_of_interest))

    # Generate color palette matching default ggplot discrete colors
    cluster_levels <- levels(violin_data_long$seurat_clusters)
    palette_colors <- scales::hue_pal()(length(cluster_levels))
    names(palette_colors) <- cluster_levels

    # Build the violin plot
    p <- ggplot(violin_data_long, aes(x = seurat_clusters, y = value, fill = seurat_clusters)) +
        geom_violin(aes(color = seurat_clusters),
                    scale = violin_scale,
                    trim = violin_trim,
                    adjust = violin_adjust,
                    show.legend = FALSE) +
        coord_flip() +
        scale_fill_manual(values = palette_colors) +
        scale_color_manual(values = palette_colors) +
        facet_wrap(~name, scales = "free_x", nrow = facet_nrow) +
        theme_minimal() +
        labs(x = "Cluster", y = "Expression Levels") +
        theme(
            legend.position = "left",
            # Color each y-axis label to match its cluster
            axis.text.y = element_text(
                size = axis_text_y_size,
                color = rev(palette_colors[cluster_levels])
            ),
            axis.text.x = element_text(size = axis_text_x_size),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            strip.text = element_text(size = strip_text_size, face = "italic",
                                      angle = strip_angle, hjust = 1, vjust = 0.7)
        )

    # Convert to grob and remove duplicate y-axes from all facets except the first
    g <- ggplotGrob(p)
    axis_grobs <- which(grepl("axis-l", g$layout$name))
    invisible(
        map(axis_grobs, function(i) {
            if (g$layout[i, "name"] != "axis-l-1-1") {
                g$grobs[[i]] <<- nullGrob()
            }
        })
    )

    # Display the plot
    grid.newpage()
    grid.draw(g)
}


######################################################################################################
# Help function
fiddle_help <- function() {
    cat("fiddle() - Create a multi-faceted violin plot with colored axis labels\n\n")
    cat("Usage:\n")
    cat("  fiddle(seurat_obj, genes_of_interest, classification_col, ...)\n\n")
    cat("Required Parameters:\n")
    cat("  seurat_obj         - Seurat object containing the data\n")
    cat("  genes_of_interest  - Character vector of gene names to plot\n")
    cat("  classification_col - Column name for classification (no default)\n\n")
    cat("Optional Parameters:\n")
    cat("  assay             - Assay to use (default: 'RNA')\n")
    cat("  axis_text_y_size  - Y-axis text size (default: 12)\n")
    cat("  axis_text_x_size  - X-axis text size (default: 6)\n")
    cat("  strip_text_size   - Facet strip text size (default: 10.5)\n")
    cat("  strip_angle       - Strip text angle (default: 85)\n")
    cat("  violin_scale      - Violin scaling method (default: 'width')\n")
    cat("  violin_trim       - Trim violin tails (default: TRUE)\n")
    cat("  violin_adjust     - Bandwidth adjustment (default: 1)\n")
    cat("  facet_nrow        - Number of facet rows (default: 1)\n\n")
    cat("Returns: Draws the plot (invisible NULL)\n\n")
    cat("Example:\n")
    cat("  fiddle(astro, c('Gfap', 'Aqp4'), 'seurat_clusters', assay = 'SCT')\n")
}
######################################################################################################


#' Perform Differential Expression Analysis with Visualization
#'
#' This function performs differential expression analysis between two cell populations
#' using Seurat's FindMarkers function and creates a volcano plot visualization. The function
#' automatically calculates appropriate x-axis limits based on the data range and provides
#' both tabular and visual outputs.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data
#' @param ident1 Character string specifying the first identity group for comparison
#' @param ident2 Character string specifying the second identity group for comparison
#' @param min_pct Numeric value specifying the minimum percentage of cells expressing 
#'   a gene in either group (default: 0.3)
#' @param test_method Character string specifying the statistical test to use 
#'   (default: "DESeq2"). Other options include "wilcox", "t", "negbinom", etc.
#' @param p_adj_threshold Numeric value for adjusted p-value threshold for significance 
#'   (default: 0.05)
#' @param fc_cutoff Numeric value for fold change cutoff in volcano plot (default: 0.001)
#' @param plot_title Character string for the main title prefix 
#'   (default: "Differentially expressed genes")
#'
#' @return A list containing:
#' \describe{
#'   \item{all_results}{Data frame with all differential expression results}
#'   \item{significant_results}{Data frame with only significant results (p_adj < threshold)}
#'   \item{volcano_plot}{EnhancedVolcano plot object}
#'   \item{xlim_used}{Numeric vector of x-axis limits used in the plot}
#'   \item{summary}{List with summary statistics including total genes, significant genes, and comparison info}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Loads required packages (dplyr, EnhancedVolcano) if not already loaded
#'   \item Runs differential expression analysis using Seurat's FindMarkers
#'   \item Filters results for statistical significance
#'   \item Calculates appropriate x-axis limits based on log2 fold change range
#'   \item Generates an enhanced volcano plot with automatic labeling
#'   \item Returns comprehensive results for further analysis
#' }
#'
#' The volcano plot automatically sets x-axis limits to integer values that encompass
#' the full range of log2 fold changes in the data, ensuring all points are visible.
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' results <- anayet(seurat_object, 
#'                   ident1 = "Astrocytes_Adu", 
#'                   ident2 = "Astrocytes_IgG")
#'
#' # Custom parameters
#' results <- anayet(seurat_object,
#'                   ident1 = "Neurons_Treated", 
#'                   ident2 = "Neurons_Control",
#'                   min_pct = 0.25,
#'                   test_method = "wilcox",
#'                   p_adj_threshold = 0.01,
#'                   plot_title = "Top DEGs")
#'
#' # Access results
#' significant_genes <- results$significant_results
#' volcano_plot <- results$volcano_plot
#' print(results$summary)
#' }
#'
#' @seealso 
#' \code{\link[Seurat]{FindMarkers}} for the underlying differential expression analysis
#' \code{\link[EnhancedVolcano]{EnhancedVolcano}} for volcano plot creation
#' \code{\link{anayet_help}} for detailed help and examples
#'
#' @import dplyr
#' @import EnhancedVolcano
#' @export
anayet <- function(seurat_obj,
    ident1,
    ident2,
    min_pct = 0.3,
    test_method = "DESeq2",
    p_adj_threshold = 0.05,
    fc_cutoff = 0.001,
    plot_title = "Differentially expressed genes",
    print_plot = TRUE,
    min_cells_group = 3) {

if (!require(dplyr, quietly = TRUE)) library(dplyr)
if (!require(EnhancedVolcano, quietly = TRUE)) library(EnhancedVolcano)

# --- Run DE ---
de_results <- FindMarkers(seurat_obj,
               min.pct = min_pct,
               ident.1 = ident1,
               ident.2 = ident2,
               test.use = test_method,
               min.cells.group = min_cells_group)

# --- Filter significant ---
de_results_significant <- de_results |>
filter(p_val_adj < p_adj_threshold)

# Color + label significant genes by adjusted p-value
sig_genes <- rownames(de_results_significant)

keyvals <- rep("grey60", nrow(de_results))
names(keyvals) <- rownames(de_results)
keyvals[sig_genes] <- "red2"

fc_range <- range(de_results$avg_log2FC, na.rm = TRUE)
xlim_lower <- floor(fc_range[1])
xlim_upper <- ceiling(fc_range[2])

# Tight y-axis: top out just above the most significant point so empty whitespace
# above the data cloud is removed (EnhancedVolcano otherwise pads y by +5)
y_max       <- max(-log10(de_results$p_val), na.rm = TRUE)
ylim_upper  <- max(ceiling(y_max * 1.05), 1)

cell_type  <- gsub("_.*", "", ident1)
condition1 <- gsub(".*_", "", ident1)
condition2 <- gsub(".*_", "", ident2)
full_title <- paste(plot_title, "in", cell_type)

volcano_plot <- EnhancedVolcano(de_results,
                     lab = rownames(de_results),
                     selectLab = sig_genes,
                     x = "avg_log2FC",
                     y = "p_val",
                     title = full_title,
                     labFace = "italic",
                     caption = "",
                     subtitle = paste(condition1, "vs.", condition2),
                     pCutoff = 1,
                     FCcutoff = fc_cutoff,
                     colCustom = keyvals,
                     legendPosition = "none",
                     drawConnectors = TRUE,
                     directionConnectors = "both",
                     max.overlaps = Inf,
                     xlim = c(xlim_lower, xlim_upper),
                     ylim = c(0, ylim_upper))

if (length(sig_genes) > 0) {
raw_boundary <- max(de_results[sig_genes, "p_val"])
volcano_plot <- volcano_plot +
geom_hline(yintercept = -log10(raw_boundary),
        linetype = "dashed", colour = "black")
}

if (print_plot) print(volcano_plot)

list(
all_results          = de_results,
significant_results  = de_results_significant,
volcano_plot         = volcano_plot,
xlim_used            = c(xlim_lower, xlim_upper),
summary = list(
total_genes       = nrow(de_results),
significant_genes = nrow(de_results_significant),
comparison        = paste(ident1, "vs", ident2)
)
)
}

#### anayet help


#' Display Help Information for the anayet Function
#'
#' This function provides comprehensive help and usage information for the anayet
#' differential expression analysis function. It explains the purpose, parameters,
#' outputs, and provides practical examples.
#'
#' @return NULL (prints help information to console)
#' @export
#' @examples
#' anayet_help()
anayet_help <- function() {
    cat("==========================================================================\n")
    cat("                           ANAYET FUNCTION HELP\n")
    cat("==========================================================================\n\n")
    
    cat("PURPOSE:\n")
    cat("--------\n")
    cat("The anayet() function performs differential expression analysis between two\n")
    cat("cell populations using Seurat and creates publication-ready volcano plots.\n\n")
    
    cat("BASIC USAGE:\n")
    cat("------------\n")
    cat("results <- anayet(seurat_obj, ident1, ident2)\n\n")
    
    cat("REQUIRED PARAMETERS:\n")
    cat("-------------------\n")
    cat("seurat_obj  : Your Seurat object containing single-cell data\n")
    cat("ident1      : First identity group for comparison (e.g., 'Astrocytes_Treated')\n")
    cat("ident2      : Second identity group for comparison (e.g., 'Astrocytes_Control')\n\n")
    
    cat("OPTIONAL PARAMETERS (with defaults):\n")
    cat("------------------------------------\n")
    cat("min_pct = 0.3           : Minimum % of cells expressing a gene\n")
    cat("test_method = 'DESeq2'  : Statistical test ('DESeq2', 'wilcox', 't', etc.)\n")
    cat("p_adj_threshold = 0.05  : Adjusted p-value cutoff for significance\n")
    cat("fc_cutoff = 0.001       : Fold change threshold for volcano plot\n")
    cat("plot_title = 'Differentially expressed genes' : Title prefix for plot\n\n")
    
    cat("WHAT THE FUNCTION DOES:\n")
    cat("----------------------\n")
    cat("1. Loads required packages automatically\n")
    cat("2. Runs differential expression analysis\n")
    cat("3. Filters for statistically significant genes\n")
    cat("4. Calculates appropriate x-axis limits for visualization\n")
    cat("5. Generates a volcano plot with automatic x-axis scaling\n")
    cat("6. Returns comprehensive results for further analysis\n\n")
    
    cat("OUTPUT:\n")
    cat("-------\n")
    cat("The function returns a list containing:\n")
    cat("$all_results        : Complete DE analysis results\n")
    cat("$significant_results: Only significant genes (p_adj < threshold)\n")
    cat("$volcano_plot       : EnhancedVolcano plot object\n")
    cat("$xlim_used         : X-axis limits used in the plot\n")
    cat("$summary           : Summary statistics\n\n")
    
    cat("EXAMPLES:\n")
    cat("---------\n")
    cat("# Basic analysis\n")
    cat("results <- anayet(my_seurat, 'Neurons_KO', 'Neurons_WT')\n\n")
    
    cat("# With custom parameters\n")
    cat("results <- anayet(my_seurat, 'Astrocytes_Treated', 'Astrocytes_Control',\n")
    cat("                  min_pct = 0.25, test_method = 'wilcox', \n")
    cat("                  p_adj_threshold = 0.01, plot_title = 'Top DEGs')\n\n")
    
    cat("# Access specific results\n")
    cat("top_genes <- results$significant_results\n")
    cat("plot <- results$volcano_plot\n")
    cat("print(results$summary)\n\n")
    
    cat("REQUIRED PACKAGES:\n")
    cat("-----------------\n")
    cat("The function will automatically load: dplyr, EnhancedVolcano\n")
    cat("Make sure Seurat is also installed and loaded in your environment.\n\n")
    
    cat("TIPS:\n")
    cat("-----\n")
    cat("- Use descriptive identity names (e.g., 'CellType_Condition')\n")
    cat("- Adjust min_pct for sparse datasets (lower values for rare cell types)\n")
    cat("- Use 'wilcox' test for faster analysis with large datasets\n")
    cat("- The volcano plot x-axis automatically scales to your data range\n\n")
    
    cat("==========================================================================\n")
}



######################################################################################################
############################### scProp ###############################################################    
######################################################################################################



#' Summarize clusters from a Seurat object
#'
#' @description
#' From a Seurat object, compute and plot the number of cells per cluster, 
#' plot per-sample cluster proportions (optionally faceted by treatment), 
#' and build a compareGroups summary table if treatment is available.
#'
#' @param obj A Seurat object.
#' @param cluster_col Character. Column name in obj@meta.data with cluster labels. Default "final_clusters_clean".
#' @param sample_col Character. Column name in obj@meta.data with sample IDs. Default "Mouse.ID".
#' @param treatment_col Character (optional). Column name in obj@meta.data with treatment or group labels. 
#'        If NULL, plots are not faceted and no summary table is created. Default "Treatment.Group".
#' @param counts_title Character. Title for the counts plot. Default "Total cells".
#' @param palette Character. RColorBrewer palette name for fills. Default "Paired".
#' @param percent_digits Integer. Number of digits for percentages. Default 1.
#'
#' @return A list with:
#' \itemize{
#' \item counts_plot: ggplot object with counts per cluster.
#' \item proportions_plot: ggplot object with per-sample cluster proportions (faceted if treatment available).
#' \item summary_table: compareGroups::createTable object if treatment_col exists, otherwise NULL.
#' \item data_counts: data frame of counts per cluster.
#' \item data_props_long: long data frame of per-sample cluster percentages.
#' \item data_props_wide: wide data frame used to build the summary table (if treatment exists).
#' }
#'
#' @examples
#' res <- scProp(microglia)
#' res$counts_plot
#' res$proportions_plot
#' res$summary_table
#'
#' @export
scProp <- function(obj,
                   cluster_col = "final_clusters_clean",
                   sample_col = "Mouse.ID",
                   treatment_col = NULL,
                   counts_title = "Total cells",
                   palette = "Paired",
                   percent_digits = 1) {
    
    # Load packages only if needed
    pkgs <- c("ggplot2","dplyr","tidyr","tibble","forcats")
    for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) stop("Package ", p, " needed.")
    if (!is.null(treatment_col)) {
        if (!requireNamespace("compareGroups", quietly = TRUE)) stop("Package compareGroups needed for summary_table.")
    }
    
    # Metadata
    if (!inherits(obj, "Seurat")) stop("obj must be a Seurat object.")
    md <- obj@meta.data
    req_cols <- c(cluster_col, sample_col)
    if (!is.null(treatment_col)) req_cols <- c(req_cols, treatment_col)
    missing_cols <- setdiff(req_cols, colnames(md))
    if (length(missing_cols) > 0) stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
    
    # Ensure cluster is factor
    cl_vec <- md[[cluster_col]]
    if (!is.factor(cl_vec)) cl_vec <- factor(cl_vec)
    cluster_levels <- levels(cl_vec)
    
    # Counts
    data_counts <- dplyr::count(md, !!rlang::sym(cluster_col), name = "n") |>
        dplyr::arrange(dplyr::desc(n)) |>
        dplyr::mutate(!!cluster_col := forcats::fct_reorder(.data[[cluster_col]], n, .desc = TRUE))
    
    counts_plot <- ggplot2::ggplot(data_counts, ggplot2::aes(x = .data[[cluster_col]], y = n, fill = .data[[cluster_col]])) +
        ggplot2::geom_col(alpha = 0.9) +
        ggplot2::labs(x = NULL, y = NULL, title = counts_title) +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_brewer(palette = palette) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none",
                       axis.text.x = ggplot2::element_text(size = 15, angle = 25),
                       axis.text.y = ggplot2::element_text(size = 15),
                       plot.title = ggplot2::element_text(size = 16, hjust = 0.5))
    
    # Per-sample cluster % 
    data_props_long <- md |>
        dplyr::count(!!rlang::sym(sample_col), !!rlang::sym(cluster_col), name = "n") |>
        dplyr::group_by(!!rlang::sym(sample_col)) |>
        dplyr::mutate(pct = round(100 * n / sum(n), percent_digits)) |>
        dplyr::ungroup() |>
        dplyr::rename(.sample = !!rlang::sym(sample_col), .cluster = !!rlang::sym(cluster_col)) |>
        dplyr::mutate(.cluster = factor(.cluster, levels = cluster_levels))
    
    if (!is.null(treatment_col)) {
        sample_index <- dplyr::distinct(md, !!rlang::sym(sample_col), !!rlang::sym(treatment_col)) |>
            dplyr::rename(.sample = !!rlang::sym(sample_col), .treat = !!rlang::sym(treatment_col))
        data_props_long <- dplyr::left_join(data_props_long, sample_index, by = ".sample")
    }
    
    proportions_plot <- ggplot2::ggplot(data_props_long, ggplot2::aes(x = .sample, y = pct, fill = .cluster)) +
        ggplot2::geom_col() +
        ggplot2::scale_fill_brewer(palette = palette) +
        ggplot2::labs(y = "%", x = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_text(size = 16),
                       axis.text.y = ggplot2::element_text(size = 14))
    
    if (!is.null(treatment_col)) {
        proportions_plot <- proportions_plot + ggplot2::facet_wrap(~.treat) +
            ggplot2::theme(strip.text = ggplot2::element_text(size = 20),
                           panel.grid.minor = ggplot2::element_blank(),
                           panel.grid.major.x = ggplot2::element_blank(),
                           panel.border = ggplot2::element_rect(fill = NA, colour = "grey80"))
    }
    
    # Optional compareGroups
    summary_table <- NULL
    data_props_wide <- NULL
    if (!is.null(treatment_col)) {
        data_props_wide <- data_props_long |>
            dplyr::select(.sample, .cluster, pct) |>
            tidyr::pivot_wider(names_from = .cluster, values_from = pct, values_fill = 0) |>
            dplyr::left_join(sample_index, by = ".sample") |>
            dplyr::rename(!!sample_col := .sample, !!treatment_col := .treat)
        
        frm <- stats::as.formula(paste0(treatment_col, " ~ . - ", sample_col))
        cg <- compareGroups::compareGroups(frm, data = data_props_wide, method = 2)
        summary_table <- compareGroups::createTable(cg)
    }
    
    list(counts_plot = counts_plot,
         proportions_plot = proportions_plot,
         summary_table = summary_table,
         data_counts = data_counts,
         data_props_long = data_props_long,
         data_props_wide = data_props_wide)
}




#' Helper for scProp
#'
#' @description
#' Show required metadata columns, the names of returned elements, and how to print the table.
#'
#' @param cluster_col Character. Cluster column. Default "final_clusters_clean".
#' @param sample_col Character. Sample column. Default "Mouse.ID".
#' @param treatment_col Character or NULL. Treatment column. Default "Treatment.Group".
#'
#' @examples
#' scProp_help()
#'
#' @export
scProp_help <- function(cluster_col = "final_clusters_clean",
                        sample_col = "Mouse.ID",
                        treatment_col = "Treatment.Group") {
    cat("Required metadata columns in obj@meta.data:\n", sep = "")
    cat(" - cluster_col: '", cluster_col, "'\n", sep = "")
    cat(" - sample_col : '", sample_col,  "'\n", sep = "")
    if (!is.null(treatment_col)) {
        cat(" - treatment_col: '", treatment_col, "' (optional)\n", sep = "")
    } else {
        cat(" - treatment_col: <none>\n", sep = "")
    }
    
    cat("\nFunction returns a named list with:\n", sep = "")
    cat(" counts_plot        (ggplot)\n")
    cat(" proportions_plot   (ggplot)\n")
    cat(" summary_table      (compareGroups table or NULL)\n")
    cat(" data_counts        (data.frame)\n")
    cat(" data_props_long    (data.frame)\n")
    cat(" data_props_wide    (data.frame or NULL)\n")
    
    cat("\nPrint the summary table (when present):\n", sep = "")
    cat(" compareGroups::export2md(res$summary_table)\n", sep = "")
    
    cat("\nMinimal usage:\n", sep = "")
    cat(" res <- scProp(obj,\n",
        "   cluster_col = '", cluster_col, "',\n",
        "   sample_col = '",  sample_col,  "',\n",
        "   treatment_col = '", ifelse(is.null(treatment_col), "", treatment_col), "' )\n", sep = "")
}
