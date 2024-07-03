#' Calculate ranks for single-cell data and generate plots
#'
#' This function calculates ranks for each gene across different cell types based on expression levels,
#' adjusts the ranks, and generates plots.
#'
#' @param so A Seurat object containing single-cell data.
#' @param assay Name of the assay in the Seurat object to use (default is "SCT").
#' @param layer Layer within the assay to use (default is "data").
#' @param grouping1 Variable in metadata to group samples by (default is "sample").
#' @param grouping2 Variable in metadata to group cell types by (default is "ct1").
#' @param pct_exp_cutoff Percentage expression cutoff to filter genes (default is 10).
#' @param output_folder Output folder path for saving plots.
#' @param group_levels Levels for grouping samples.
#'
#' @return Modified Seurat object with ranks added as an assay.
#'
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @importFrom ggforce facet_wrap_paginate
#' @importFrom utils txtProgressBar
#' @importFrom Seurat CreateAssayObject DotPlot
#' @export
calculateRanks <- function(so, assay = "SCT", layer = "data", grouping1 = "sample", grouping2 = "ct1", pct_exp_cutoff = 10, output_folder, group_levels = NA) {
  require(dplyr)
  require(tidyverse)
  require(ggplot2)
  require(ggforce)
  require(utils)
  require(Seurat)
  require(Rcpp)
  require(data.table)
  require(maditr)

  # Source the C++ function
  sourceCpp("./../RESRECT/src/calculate_pct_exp.cpp")

  # Helper Function to replace values with NA based on the threshold
  replace_with_na <- function(data, lookup, lookup_column = "lookup_grouping", grouping_column = "ct1", gene_column = "gene") {
    pb <- utils::txtProgressBar(min = 0, max = nrow(lookup), style = 3, width = 50, char = "=")
    for (i in seq_len(nrow(lookup))) {
      ct_value <- lookup[[lookup_column]][i]
      gene_value <- lookup[[gene_column]][i]
      if (gene_value %in% colnames(data)) {
        data[data[[grouping_column]] == ct_value, gene_value] <- NA
      }
      utils::setTxtProgressBar(pb, i)
    }
    close(pb)
    return(data)
  }
  data <- so@assays[[assay]][layer]
  metadata <- so@meta.data

  # Merge assay data and metadata
  metadata$barcode <- rownames(metadata)
  data <- as.data.frame(t(data))
  data$barcode <- rownames(data)
  data <- dplyr::left_join(data, metadata, by = "barcode") %>% as.data.frame()
  rownames(data) <- data$barcode

  # Determine which genes to include in the calculation of median rank
  features <- dplyr::select(data, -colnames(metadata)) %>% colnames()

  # Call the C++ function to calculate pct.exp
  pct_exp_matrix <- calculatePctExp(as.matrix(data[, features]), as.character(data[[grouping2]]), features)

  lookup <- as.data.frame(pct_exp_matrix) %>%
    rownames_to_column(var = "lookup_grouping") %>%
    pivot_longer(-lookup_grouping, names_to = "gene", values_to = "pct.exp")

  lookup <- lookup %>% dplyr::filter(pct.exp < pct_exp_cutoff)

  # Set the cell-x-gene values of the gene-ct-pairs not passing the expression threshold to NA
  print("Filtering for Expression Level")
  data <- replace_with_na(as.data.frame(data), lookup, grouping_column=grouping2)

  num_nas <- sum(is.na(data[, features]))
  num_non_nas <- sum(!is.na(data[, features]))
  print(paste("Percent of entries passing threshold:", (1 - num_nas / (num_nas + num_non_nas)) * 100, "%"))

  # Calculate the rank for each cell and gene within each cell type
  print("calculating ranks")
  df_ranked <- as.data.table(data)
  rm(data)
  gc()
  df_ranked[, (features) := lapply(.SD, function(x) rank(-x, ties.method = "average", na.last = "keep")), by = grouping2, .SDcols = features]
  setDF(df_ranked)
  rownames(df_ranked) <- df_ranked$barcode

  # Create an assay object with ranked data
  #so[["ranks"]] <- CreateAssayObject(counts = t(df_ranked %>% dplyr::select(-colnames(metadata))))
  so[["ranks"]]$counts <- t(df_ranked %>% dplyr::select(-colnames(metadata)))

  # Reshape data for plotting
  print("reshaping data for plotting")
  gc()
  melted_data <- data.table::melt(df_ranked, id.vars = colnames(metadata), variable.name = "gene", value.name = "rank", na.rm = TRUE)
  rm(df_ranked)
  gc()
  melted_data <- melted_data[!is.na(melted_data$rank), ]
  melted_data$group <- paste(melted_data$treatment, melted_data$genotype, melted_data$sex, melted_data$orig.ident, sep = "-")

  melted_data %>%
    group_by(!!sym(grouping2)) %>%
    summarise(cell_count = n_distinct(barcode)) %>% print(n=50)
  # if I adjust the ranks of all genes and not just the genes with 10% of expression in each ct1, will this be okay?
  # The adjusted median ranks wont be zero any more because there are a diffrerent number of cells in each ct bc of the filtering criteria so the rank range will change
  # but maybe at least the adjustment will still be scaled to the capture the trend of the rank difference

  # Call plotRanks function for original ranks
  orig_ranks_folder <- file.path(output_folder, "orig_ranks")
  if (!dir.exists(orig_ranks_folder)) {
    dir.create(orig_ranks_folder, recursive = TRUE)
  }
  gc()
  print("Plotting Original Ranks")
  # ETHAN UNCOMMENT THE FOLLOWING LINE for final use and delete this comment !!!
  #plotRanks(melted_data, orig_ranks_folder, group_levels)

  # Adjust ranks
  print("Adjusting Ranks")
  melted_data <- adjustRanks(melted_data, grouping1 = grouping1, grouping2 = grouping2)

  # Create subfolders for original and adjusted ranks plots
  adj_ranks_folder <- file.path(output_folder, "adj_ranks")
  if (!dir.exists(adj_ranks_folder)) {
    dir.create(adj_ranks_folder, recursive = TRUE)
  }

  # Plotting adjusted ranks
  print("Plotting Adjusted Ranks")
  # ETHAN UNCOMMENT THE FOLLOWING LINE for final use and delete this comment !!!
  #plotRanks(melted_data, adj_ranks_folder, adjusted=TRUE)
  gc()

  # Add adjusted ranks to the Seurat object
  melted_data <- melted_data %>%
    maditr::dcast(barcode ~ gene, value.var = "adj_rank") %>%
    column_to_rownames("barcode")
  gc()
  #so[["adj_ranks"]] <- CreateAssayObject(counts = t(as.matrix(melted_data)))
  so[["adj_ranks"]]$counts <- t(as.matrix(melted_data))
  return(so)
}
