#' Adjust ranks and add them to the Seurat object
#'
#' @param so Seurat object
#' @param grouping2 Second grouping variable
#' @return Seurat object with adjusted ranks
#' @export
adjustRanks <- function(so, grouping2="ct1") {
  df_ranked <- as.data.frame(so@assays[["ranks"]]@counts)
  df_ranked$barcode <- rownames(df_ranked)
  metadata <- so@meta.data %>% dplyr::select(c("treatment", "sex", "genotype", "orig.ident", grouping2))
  metadata$barcode <- rownames(metadata)
  df_ranked <- dplyr::left_join(df_ranked, metadata, by="barcode")
  rownames(df_ranked) <- df_ranked$barcode
  df_ranked <- df_ranked[, -grep('^"', colnames(df_ranked))]

  melted_data <- tidyr::pivot_longer(df_ranked, cols = -c(barcode, grouping2, "orig.ident", "genotype", "sex", "treatment"), names_to = "gene", values_to = "rank")
  melted_data <- melted_data[!is.na(melted_data$rank), ]
  melted_data$group <- paste(melted_data$treatment, melted_data$genotype, melted_data$sex, melted_data$orig.ident, sep="-")

  summary_stats <- melted_data %>%
    dplyr::group_by(group, !!rlang::sym(grouping2)) %>%
    dplyr::summarise(median_value = stats::median(rank, na.rm=TRUE),
                     mean_value = mean(rank, na.rm=TRUE))

  melted_data <- dplyr::left_join(melted_data, summary_stats, by = c("group", grouping2))
  melted_data$median_adj_rank <- melted_data$rank - melted_data$median_value
  melted_data$mean_adj_rank <- melted_data$rank - melted_data$mean_value

  summary_stats <- melted_data %>%
    dplyr::group_by(group, !!rlang::sym(grouping2)) %>%
    dplyr::summarise(median_value = stats::median(rank, na.rm=TRUE),
                     median_median_adj_value = stats::median(median_adj_rank, na.rm=TRUE),
                     mean_median_adj_value = mean(median_adj_rank, na.rm=TRUE),
                     mean_value = mean(rank, na.rm=TRUE),
                     median_mean_adj_value = stats::median(mean_adj_rank, na.rm=TRUE),
                     mean_mean_adj_value = mean(mean_adj_rank, na.rm=TRUE))

  cell_counts <- melted_data %>%
    dplyr::filter(gene == "Actb") %>%
    dplyr::group_by(group, !!rlang::sym(grouping2)) %>%
    dplyr::summarise(cell_count = dplyr::n())

  summary_stats <- dplyr::left_join(summary_stats, cell_counts, by=c("group", grouping2))

  df_adjusted <- melted_data %>%
    dplyr::select(barcode, gene, median_adj_rank) %>%
    tidyr::pivot_wider(names_from = "gene", values_from = "median_adj_rank")

  rownames(df_adjusted) <- df_adjusted$barcode
  df_adjusted <- df_adjusted %>% dplyr::select(-barcode)
  df_adjusted_t <- t(df_adjusted)
  new_assay <- Seurat::CreateAssayObject(counts = df_adjusted_t)
  so[["ranks_adj"]] <- new_assay

  return(so)
}
