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
adjustRanks <- function(melted_data, grouping1 = "sample", grouping2 = "ct1", avg_fxn="median") {
  require(dplyr)
  # Convert the function name to a function object
  avg_function <- match.fun(avg_fxn)
  # Calculate the average rank for each cell type (ct1) and sample
  avg_ranks <- melted_data %>%
    dplyr::group_by(!!sym(grouping2), !!sym(grouping1)) %>%
    dplyr::summarise(avg_rank = avg_function(rank, na.rm = TRUE)) %>%
    ungroup()

  # Join the median ranks with the original data
  melted_data <- melted_data %>%
    dplyr::left_join(avg_ranks, by = c(grouping2, grouping1)) %>%
    dplyr::mutate(adj_rank = rank - avg_rank)

  # within each cell type, ensure all ranks are >= 0 by substracting the (presumably negative) lowest rank from each
  melted_data <- melted_data %>%
    dplyr::group_by(!!sym(grouping2)) %>%
    dplyr::mutate(min_ct_rank = min(adj_rank))
  melted_data$adj_rank <- melted_data$adj_rank - melted_data$min_ct_rank

  return(melted_data)
}
