#' Calculate Ranks
#'
#' @param so A Seurat object.
#' @param assay Assay to use.
#' @param layer Data layer to use.
#' @param grouping1 First grouping variable.
#' @param grouping2 Second grouping variable.
#' @param pct_exp_cutoff Percent expression cutoff.
#' @param output_folder Output folder.
#' @return A Seurat object with ranks calculated.
#' @import dplyr
#' @import tidyverse
#' @importFrom Seurat DotPlot CreateAssayObject
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats median
#' @export
calculateRanks <- function(so, assay = "SCT", layer = "data", grouping1 = "sample", grouping2 = "ct1", pct_exp_cutoff = 10, output_folder = "outputs_v2/DE/single_cell/corrected_ranks/resRect/") {
  requireNamespace("dplyr")
  requireNamespace("tidyverse")
  requireNamespace("Seurat")
  requireNamespace("utils")

  DefaultAssay(so) <- assay
  data <- so[[layer]]
  groupings <- data.frame(sample = so[[grouping1]], ct1 = so[[grouping2]])
  groupings <- groupings %>%
    dplyr::mutate(lookup_grouping = paste(sample, ct1, sep = "_"))

  dotplot <- Seurat::DotPlot(so, features = rownames(data), group.by = "lookup_grouping")
  df <- dotplot$data %>%
    dplyr::select(features.plot, pct.exp) %>%
    dplyr::group_by(features.plot) %>%
    dplyr::mutate(mean_pct.exp = mean(pct.exp)) %>%
    dplyr::ungroup()

  df_filtered <- df %>%
    dplyr::filter(mean_pct.exp > pct_exp_cutoff)

  filtered_genes <- unique(df_filtered$features.plot)

  pct_ranks <- data.frame(matrix(nrow = nrow(data), ncol = 0))
  rownames(pct_ranks) <- rownames(data)

  pbar <- utils::txtProgressBar(min = 0, max = length(unique(groupings$lookup_grouping)), style = 3)
  for (i in unique(groupings$lookup_grouping)) {
    idx <- which(groupings$lookup_grouping == i)
    if (length(idx) > 1) {
      pct <- apply(data[filtered_genes, idx, drop = FALSE], 1, function(x) mean(x > 0))
      pct_ranks <- cbind(pct_ranks, rank(-pct, ties.method = "min"))
    }
    utils::setTxtProgressBar(pbar, which(i == unique(groupings$lookup_grouping)))
  }
  close(pbar)

  colnames(pct_ranks) <- unique(groupings$lookup_grouping)

  pct_ranks <- dplyr::mutate(pct_ranks, gene = rownames(pct_ranks))
  pct_ranks_long <- pct_ranks %>%
    tidyr::pivot_longer(-gene, names_to = "lookup_grouping", values_to = "rank")

  pct_ranks_wide <- pct_ranks_long %>%
    tidyr::pivot_wider(names_from = "lookup_grouping", values_from = "rank")

  ranks_assay <- Seurat::CreateAssayObject(counts = as.matrix(pct_ranks_wide))
  so[["ranks"]] <- ranks_assay

  return(so)
}

globalVariables(c("lookup_grouping", "pct.exp", "features.plot"))
