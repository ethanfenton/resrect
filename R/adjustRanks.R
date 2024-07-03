#' Adjust Ranks for Single-Cell Data
#'
#' This function adjusts ranks for each gene across different cell types based on expression levels.
#'
#' @param melted_data Data frame containing melted data with columns: gene, ct1, sample, rank, and other metadata columns.
#' @param grouping1 Variable in metadata to group samples by (default is "sample").
#' @param grouping2 Variable in metadata to group cell types by (default is "ct1").
#'
#' @return Data frame with an additional column for adjusted ranks.
#'
#' @import dplyr
#' @export
adjustRanks <- function(melted_data, grouping1 = "sample", grouping2 = "ct1") {
  require(dplyr)
  # Calculate the median rank for each cell type (ct1) and sample
  median_ranks <- melted_data %>%
    dplyr::group_by(!!sym(grouping2), !!sym(grouping1)) %>%
    dplyr::summarise(median_rank = median(rank, na.rm = TRUE)) %>%
    ungroup()

  # Join the median ranks with the original data
  melted_data <- melted_data %>%
    dplyr::left_join(median_ranks, by = c(grouping2, grouping1)) %>%
    dplyr::mutate(adj_rank = rank - median_rank)

  return(melted_data)
}
