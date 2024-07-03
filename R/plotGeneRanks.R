#' Plot Gene Ranks
#'
#' @param so Seurat object containing the single-cell data.
#' @param gene The gene name to highlight in the plot.
#' @param ct1 The cell type to plot.
#' @param output_folder Output folder path for saving plots.
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @importFrom rlang sym
#' @export
plotGeneRanks <- function(so, gene, ct1, output_folder) {
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(rlang)
  require(ggforce)

  # Extract metadata
  metadata <- so@meta.data

  # Helper function to plot and save ranks
  plot_single_rank <- function(data, plot_title, output_file) {
    ggplot(data, aes(x = sample, y = rank)) +
      geom_jitter(aes(color = (gene == gene), alpha = (gene == gene)), size = 1.5) +
      scale_color_manual(values = c("grey", "red")) +
      scale_alpha_manual(values = c(0.2, 1)) +
      labs(x = "Sample", y = "Rank", title = plot_title) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 70, hjust = 1),
        legend.position = "none"
      ) +
      ggsave(filename = output_file, height = 6, width = 8, units = "in", dpi = 300)
  }

  # Plot the original ranks
  {
    # Extract and process the original rank data
    original_ranks <- as.data.frame(t(so@assays$ranks@counts))
    original_ranks$barcode <- rownames(original_ranks)
    original_ranks <- merge(original_ranks, metadata, by = "barcode")

    # Filter for the specified cell type and gene
    original_data <- original_ranks %>% filter(!!sym("ct1") == ct1)

    # Melt the data for plotting
    original_data_melted <- original_data %>%
      pivot_longer(cols = -c(barcode, sample, ct1), names_to = "gene", values_to = "rank") %>%
      filter(gene == gene)

    # Plot the original ranks
    plot_single_rank(
      original_data_melted,
      plot_title = paste("Original Ranks for Gene", gene, "in Cell Type", ct1),
      output_file = file.path(output_folder, paste0("original_ranks_", gene, "_", ct1, ".png"))
    )

    # Clean up
    rm(original_ranks, original_data, original_data_melted)
    gc()
  }

  # Plot the adjusted ranks
  {
    # Extract and process the adjusted rank data
    adjusted_ranks <- as.data.frame(t(so@assays$ranks@adj_ranks))
    adjusted_ranks$barcode <- rownames(adjusted_ranks)
    adjusted_ranks <- merge(adjusted_ranks, metadata, by = "barcode")

    # Filter for the specified cell type and gene
    adjusted_data <- adjusted_ranks %>% filter(!!sym("ct1") == ct1)

    # Melt the data for plotting
    adjusted_data_melted <- adjusted_data %>%
      pivot_longer(cols = -c(barcode, sample, ct1), names_to = "gene", values_to = "rank") %>%
      filter(gene == gene)

    # Plot the adjusted ranks
    plot_single_rank(
      adjusted_data_melted,
      plot_title = paste("Adjusted Ranks for Gene", gene, "in Cell Type", ct1),
      output_file = file.path(output_folder, paste0("adjusted_ranks_", gene, "_", ct1, ".png"))
    )

    # Clean up
    rm(adjusted_ranks, adjusted_data, adjusted_data_melted)
    gc()
  }
}
