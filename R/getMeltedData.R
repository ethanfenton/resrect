#' Get Melted DataFrame
#'
#' @param so Seurat object containing the single-cell data.
#' @return A melted dataframe with columns: barcode, gene, ct1, rank_orig, rank_adj, median_rank.
#'
#' @import dplyr
#' @import tidyr
#' @importFrom rlang sym
#' @export
getMeltedData <- function(so) {

  metadata <- so@meta.data
  # Extract and process the original rank data
  original_ranks <- as.data.frame(t(so@assays$ranks@counts))
  original_ranks$barcode <- rownames(original_ranks)
  metadata$barcode <- rownames(metadata)
  original_ranks <- dplyr::left_join(original_ranks, metadata, by = "barcode")
  # Melt the data for plotting
  original_ranks <- original_ranks %>%
    tidyr::pivot_longer(cols = -colnames(metadata), names_to = "gene", values_to = "rank")

  # Extract and process the adjusted rank data
  adj_ranks <- as.data.frame(t(so@assays$adj_rank@counts))
  adj_ranks$barcode <- rownames(adj_ranks)
  # Melt the data for plotting
  adj_ranks <- adj_ranks %>%
    tidyr::pivot_longer(cols=-c(barcode), names_to = "gene", values_to = "adj_rank")

  # combine the original and adjusted ranks
  return(dplyr::left_join(original_ranks, adj_ranks, by="barcode"))

}





getMeltedData <- function(so) {
  library(data.table)
  library(dplyr)
  library(tidyr)
  # Convert metadata to data.table
  metadata <- as.data.table(so@meta.data)
  metadata[, barcode := rownames(metadata)]

  # Extract and process the original rank data
  original_ranks <- as.data.table(t(so@assays$ranks@counts))
  setnames(original_ranks, old = colnames(original_ranks), new = rownames(so@assays$ranks@counts))
  original_ranks[, barcode := rownames(original_ranks)]
  original_ranks <- merge(original_ranks, metadata, by = "barcode")

  # Extract and process the adjusted rank data
  adj_ranks <- as.data.table(t(so@assays$adj_ranks@counts))
  setnames(adj_ranks, old = colnames(adj_ranks), new = rownames(so@assays$adj_ranks@counts))
  adj_ranks[, barcode := rownames(adj_ranks)]
  adj_ranks <- merge(adj_ranks, metadata, by = "barcode")

  # Chunk processing to avoid memory issues
  genes <- setdiff(names(original_ranks), c("barcode", colnames(metadata)))

  # Initialize the melted_data as an empty data.table
  melted_data <- data.table()

  # Process in chunks
  chunk_size <- 1000  # Adjust this size according to your memory capacity
  num_chunks <- ceiling(length(genes) / chunk_size)

  for (i in seq(1, length(genes), by = chunk_size)) {
    chunk_genes <- genes[i:min(i + chunk_size - 1, length(genes))]
    chunk_index <- (i - 1) %/% chunk_size + 1
    print(paste("Processing chunk", chunk_index, "of", num_chunks))

    # Process original ranks
    chunk_original <- original_ranks[, c("barcode", chunk_genes, colnames(metadata)), with = FALSE]
    chunk_melted_orig <- melt(chunk_original, id.vars = c("barcode", colnames(metadata)), variable.name = "gene", value.name = "rank_orig")

    # Process adjusted ranks
    chunk_adj <- adj_ranks[, c("barcode", chunk_genes, colnames(metadata)), with = FALSE]
    chunk_melted_adj <- melt(chunk_adj, id.vars = c("barcode", colnames(metadata)), variable.name = "gene", value.name = "rank_adj")

    # Merge the chunks
    setkey(chunk_melted_orig, barcode, gene)
    setkey(chunk_melted_adj, barcode, gene)
    chunk_combined <- merge(chunk_melted_orig, chunk_melted_adj, by = c("barcode", "gene", colnames(metadata)), all = TRUE)

    # Combine with melted_data
    melted_data <- rbind(melted_data, chunk_combined)

    # Clean up
    rm(chunk_original, chunk_melted_orig, chunk_adj, chunk_melted_adj, chunk_combined)
    gc()  # Collect garbage to free memory

    print(paste("Completed chunk", chunk_index, "of", num_chunks))
  }

  return(melted_data)
}
