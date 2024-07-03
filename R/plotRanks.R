#' Plot Ranks
#'
#' @param melted_data Dataframe containing the melted data for plotting.
#' @param output_folder Output folder path for saving plots.
#' @param group_levels Levels for grouping samples.
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @importFrom rlang sym
#' @export
plotRanks <- function(melted_data, output_folder, group_levels=NA, adjusted=FALSE) {
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(rlang)
  require(ggforce)

  if(adjusted){
    melted_data$rank <- melted_data$adj_rank
  }
  # Calculate stats of each group
  summary_stats <- melted_data %>%
    dplyr::group_by(group, !!sym("ct1")) %>%
    dplyr::summarise(
      median_value = median(rank, na.rm = TRUE),
      mean_value = mean(rank, na.rm = TRUE)
    )

  # Calculate number of cells in each group
  cell_counts <- melted_data %>%
    group_by(group, !!sym("ct1")) %>%
    summarise(cell_count = n_distinct(barcode))

  summary_stats <- dplyr::left_join(summary_stats, cell_counts, by = c("group", "ct1"))
  melted_data <- dplyr::left_join(melted_data, summary_stats, by = c("group", "ct1"))

  if (!is.na(group_levels)) {
    melted_data$group <- factor(melted_data$group, levels = group_levels)
  }

  # Plotting with ggforce to paginate facets
  num_ct1 <- length(unique(melted_data$ct1))
  for (i in 1:num_ct1) {
    p <- ggplot(melted_data, aes(x = group, y = rank)) +
      geom_boxplot(aes(fill = orig.ident, color = treatment)) +
      ggforce::facet_wrap_paginate(~ct1, ncol = 1, nrow = 1, page = i, scales = "free_y") +
      geom_label(data = summary_stats, aes(label = paste("Median:\n", round(median_value, 2)), y = median_value + 10), vjust = -1, color = "black", size = 3, fill = "white", label.padding = unit(0.5, "lines"), label.r = unit(0.2, "lines")) +
      geom_label(data = summary_stats, aes(label = paste("Count:\n", cell_count), y = median_value - 15), vjust = 2, color = "black", size = 3, fill = "white", label.padding = unit(0.5, "lines"), label.r = unit(0.2, "lines")) +
      labs(x = "Sample", y = "Rank") +
      scale_color_manual(values = c("black", "blue")) +
      theme(
        axis.text.x = element_text(angle = 70, hjust = 1),
        strip.text = element_text(size = 10),  # Adjusted text size for better fitting
        legend.position = "bottom"  # Move legend to the bottom for better space usage
      )

    # Extract the facet label for the current page
    facet_label <- levels(factor(melted_data$ct1))[i]

    # Set the plot title
    if(adjusted){
      p <- p + ggtitle(paste("Adj Rank Distribution by Pool -", facet_label))
      # Save each page as a separate PNG
      ggsave(filename = paste0(output_folder, "/adjusted_ranks_by_ct_", facet_label, ".png"), plot = p, height = 10, width = 10, units = "in", dpi = 300)
    }else{
      p <- p + ggtitle(paste("Orig Rank Distribution by Pool -", facet_label))
      # Save each page as a separate PNG
      ggsave(filename = paste0(output_folder, "/unadjusted_ranks_by_ct_", facet_label, ".png"), plot = p, height = 10, width = 10, units = "in", dpi = 300)
    }
  }

}
