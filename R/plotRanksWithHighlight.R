#' Plot Ranks with Highlighted Gene
#'
#' This function plots both the original ranks and the adjusted ranks, highlighting a specified gene.
#'
#' @param so A Seurat object
#' @param grouping1 First grouping column name
#' @param grouping2 Second grouping column name
#' @param output_folder Folder to save the plots
#' @param gene The gene to be highlighted
#' @import dplyr tidyverse ggplot2 Seurat
#' @export
plotRanksWithHighlight <- function(so, grouping1="sample", grouping2="ct1", output_folder="outputs_v2/DE/single_cell/corrected_ranks/resRect/", gene) {
  require(dplyr)
  require(tidyverse)
  require(ggplot2)

  # Function to create plots
  create_plot <- function(data, rank_column, plot_title, filename) {
    summary_stats <- data %>%
      group_by(group, !!sym(grouping2)) %>%
      summarise(median_value = median(!!sym(rank_column), na.rm=TRUE),
                mean_value = mean(!!sym(rank_column), na.rm=TRUE))

    p <- ggplot(data, aes(x=group, y=!!sym(rank_column))) +
      geom_boxplot(aes(fill=orig.ident, color=treatment), outlier.shape=NA) +
      geom_jitter(alpha=0.3, color="grey") +
      geom_point(data = data %>% filter(gene == !!gene), aes(x=group, y=!!sym(rank_column)), color="red", size=3) +
      facet_wrap(stats::as.formula(paste("~", grouping2))) +
      labs(x = "Sample", y = "Rank", title = plot_title) +
      scale_color_manual(values=c("black","blue")) +
      theme(axis.text.x = element_text(angle=70, hjust=1), strip.text = element_text(size=30)) +
      ggtitle(plot_title)

    ggsave(p, file=paste0(output_folder, filename), height=40, width=40, units="in", dpi=600)
    print(p)
  }

  # Prepare data for original ranks
  original_ranks <- as.data.frame(so@assays$ranks@counts)
  metadata <- so@meta.data %>% dplyr::select("treatment", "sex", "genotype", "orig.ident", !!grouping2)
  original_ranks$barcode <- rownames(original_ranks)
  melted_original_data <- original_ranks %>%
    pivot_longer(cols = -barcode, names_to = "gene", values_to = "rank")
  melted_original_data <- melted_original_data[!is.na(melted_original_data$rank),]
  mdata <- metadata %>% dplyr::select(c("treatment", "sex", "genotype", "orig.ident", !!grouping2))
  mdata$barcode <- rownames(mdata)
  melted_original_data <- dplyr::left_join(melted_original_data, mdata, by="barcode")
  melted_original_data$group <- paste(melted_original_data$treatment, melted_original_data$genotype, melted_original_data$sex, melted_original_data$orig.ident, sep="-")
  if (!is.null(group_levels)) {
    melted_original_data$group <- factor(melted_original_data$group, levels = group_levels)
  }

  # Create plot for original ranks
  create_plot(melted_original_data, "rank", "Original Rank Distribution with Highlighted Gene", "highlighted_original_ranks_by_ct.pdf")

  # Prepare data for adjusted ranks
  adjusted_ranks <- as.data.frame(so@assays$median_adj_ranks@counts)
  adjusted_ranks$barcode <- rownames(adjusted_ranks)
  melted_adjusted_data <- adjusted_ranks %>%
    pivot_longer(cols = -barcode, names_to = "gene", values_to = "median_adj_rank")
  melted_adjusted_data <- melted_adjusted_data[!is.na(melted_adjusted_data$median_adj_rank),]
  melted_adjusted_data <- dplyr::left_join(melted_adjusted_data, mdata, by="barcode")
  melted_adjusted_data$group <- paste(melted_adjusted_data$treatment, melted_adjusted_data$genotype, melted_adjusted_data$sex, melted_adjusted_data$orig.ident, sep="-")
  if (!is.null(group_levels)) {
    melted_adjusted_data$group <- factor(melted_adjusted_data$group, levels = group_levels)
  }

  # Create plot for adjusted ranks
  create_plot(melted_adjusted_data, "median_adj_rank", "Adjusted Rank Distribution with Highlighted Gene", "highlighted_adjusted_ranks_by_ct.pdf")
}
