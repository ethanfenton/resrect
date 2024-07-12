#' Plot Gene Ranks
#'
#' @param so Seurat object containing the single-cell data.
#' @param plot_genes A list of gene names to highlight in the plots.
#' @param ct1 The cell type to plot.
#' @param output_folder Output folder path for saving plots.
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @importFrom rlang sym
#' @export
plotGeneRanks <- function(so, plot_genes, ct1, output_folder) {
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(rlang)
  require(ggforce)

  # Filter Seurat object by cell type
  so <- subset(so, ct1 == !!ct1)

  # Extract metadata
  metadata <- so@meta.data

  # Helper function to plot and save ranks
  plot_single_rank <- function(data, plot_gene, plot_title, output_file) {
    grey_data <- data %>% filter(gene != plot_gene)
    red_data <- data %>% filter(gene == plot_gene)

    # Calculate median rank for each sample
    median_ranks <- data %>%
      group_by(sample) %>%
      summarise(median_rank = median(rank, na.rm = TRUE))

    p <- ggplot() +
      geom_jitter(data = grey_data, aes(x = sample, y = rank), color = "grey", alpha = 0.2, size = 1.5) +
      geom_jitter(data = red_data, aes(x = sample, y = rank), color = "red", alpha = 1, size = 1.5) +
      geom_point(data = median_ranks, aes(x = sample, y = median_rank), color = "blue", size = 3, shape = 18) +
      labs(x = "Sample", y = "Rank", title = plot_title) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 70, hjust = 1),
        legend.position = "none"
      )
    ggsave(p, filename = output_file, height = 6, width = 8, units = "in", dpi = 300)
  }

  # Create subfolders for original and adjusted ranks
  orig_ranks_folder <- file.path(output_folder, "original_ranks")
  adj_ranks_folder <- file.path(output_folder, "adjusted_ranks")
  dir.create(orig_ranks_folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(adj_ranks_folder, recursive = TRUE, showWarnings = FALSE)

  # Plot the original ranks
  {
    # Extract and process the original rank data
    original_ranks <- as.data.frame(t(so@assays$ranks@counts))
    original_ranks$barcode <- rownames(original_ranks)
    metadata$barcode <- rownames(metadata)
    original_ranks <- merge(original_ranks, metadata, by = "barcode")

    # Filter for the specified cell type
    original_data <- original_ranks %>% filter(!!sym("ct1") == ct1)

    # Melt the data for plotting
    original_data_melted <- original_data %>%
      pivot_longer(cols = -colnames(metadata), names_to = "gene", values_to = "rank")

    # Loop over the genes and create plots
    for (plot_gene in plot_genes) {
      plot_single_rank(
        original_data_melted,
        plot_gene = plot_gene,
        plot_title = paste("Original Ranks:", plot_gene, "Highlighted in", ct1),
        output_file = file.path(orig_ranks_folder, paste0(plot_gene, "_", ct1, "_original_ranks.png"))
      )
    }

    # Clean up
    rm(original_ranks, original_data, original_data_melted)
    gc()
  }

  # Plot the adjusted ranks
  {
    # Extract and process the adjusted rank data
    adjusted_ranks <- as.data.frame(t(so@assays$adj_ranks@counts))  # Note the use of "data" layer for adjusted ranks
    adjusted_ranks$barcode <- rownames(adjusted_ranks)
    adjusted_ranks <- merge(adjusted_ranks, metadata, by = "barcode")

    # Filter for the specified cell type
    adjusted_data <- adjusted_ranks %>% filter(!!sym("ct1") == ct1)

    # Melt the data for plotting
    adjusted_data_melted <- adjusted_data %>%
      pivot_longer(cols = -colnames(metadata), names_to = "gene", values_to = "rank")

    # Loop over the genes and create plots
    for (plot_gene in plot_genes) {
      plot_single_rank(
        adjusted_data_melted,
        plot_gene = plot_gene,
        plot_title = paste("Adjusted Ranks for Gene", plot_gene, "in Cell Type", ct1),
        output_file = file.path(adj_ranks_folder, paste0("adjusted_ranks_", plot_gene, "_", ct1, ".png"))
      )
    }

    # Clean up
    rm(adjusted_ranks, adjusted_data, adjusted_data_melted)
    gc()
  }
}

