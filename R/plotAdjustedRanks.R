#' Plot adjusted ranks
#'
#' @param so Seurat object
#' @param output_folder Output folder path
#' @param group_levels Group levels
#' @export
plotAdjustedRanks <- function(so, output_folder, group_levels) {
  df_adjusted <- as.data.frame(so@assays[["ranks_adj"]]@counts)
  df_adjusted$barcode <- rownames(df_adjusted)
  metadata <- so@meta.data %>% dplyr::select(c("treatment", "sex", "genotype", "orig.ident", "ct1"))
  metadata$barcode <- rownames(metadata)
  df_adjusted <- dplyr::left_join(df_adjusted, metadata, by="barcode")
  rownames(df_adjusted) <- df_adjusted$barcode
  df_adjusted <- df_adjusted[, -grep('^"', colnames(df_adjusted))]

  melted_data <- tidyr::pivot_longer(df_adjusted, cols = -c(barcode, ct1, orig.ident, genotype, sex, treatment), names_to = "gene", values_to = "rank")
  melted_data <- melted_data[!is.na(melted_data$rank), ]
  melted_data$group <- paste(melted_data$treatment, melted_data$genotype, melted_data$sex, melted_data$orig.ident, sep="-")

  summary_stats <- melted_data %>%
    dplyr::group_by(group, !!rlang::sym("ct1")) %>%
    dplyr::summarise(median_adj_rank = stats::median(rank, na.rm=TRUE),
                     mean_adj_rank = mean(rank, na.rm=TRUE))

  cell_counts <- melted_data %>%
    dplyr::filter(gene == "Actb") %>%
    dplyr::group_by(group, !!rlang::sym("ct1")) %>%
    dplyr::summarise(cell_count = dplyr::n())

  summary_stats <- dplyr::left_join(summary_stats, cell_counts, by=c("group", "ct1"))
  melted_data <- dplyr::left_join(melted_data, summary_stats, by = c("group", "ct1"))
  melted_data$group <- factor(melted_data$group, levels = group_levels)

  p <- ggplot2::ggplot(melted_data, ggplot2::aes(x=group, y=rank)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill=orig.ident, color=treatment)) +
    ggplot2::facet_wrap(~ct1) +
    ggplot2::geom_text(data = summary_stats, ggplot2::aes(label = paste("Median:\n", round(median_adj_rank, 2)), y = median_adj_rank + 10), vjust = -1, color = "black") +
    ggplot2::geom_text(data = summary_stats, ggplot2::aes(label = paste("Mean:\n", round(mean_adj_rank, 2)), y = median_adj_rank), vjust = 1, color = "black") +
    ggplot2::geom_text(data = summary_stats, ggplot2::aes(label = paste("Count:\n", cell_count), y = median_adj_rank - 15), vjust = 2, color = "black") +
    ggplot2::labs(x = "Sample", y = "Adjusted Rank") +
    ggplot2::scale_color_manual(values=c("black", "blue")) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=70, hjust=1), strip.text = ggplot2::element_text(size=30)) +
    ggplot2::ggtitle("Adjusted Rank Distribution by Pool")

  ggplot2::ggsave(filename = paste0(output_folder, "adjusted_ranks_by_ct.pdf"), plot = p, height = 40, width = 40, units = "in", dpi = 600)
}
