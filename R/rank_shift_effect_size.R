#' Calculate Rank Shift Effect Size
#'
#' This function calculates the rank shift effect size as a metric akin to fold change but distinct from it.
#'
#' @param de_df Data frame containing differential expression data.
#' @param so A Seurat object containing single-cell data.
#' @param assay Name of the assay in the Seurat object to use (default is "median_adj_ranks").
#' @param layer Layer within the assay to use (default is "counts").
#'
#' @return A data frame with the average rank data by group.
#'
#' @import dplyr
#' @export
rank_shift_effect_size <- function(rank_de_df, orig_de_df, so, assay = "adj_ranks", layer = "counts",
                                   grouping1 = "treatment",  grouping2="ct1",
                                   cell_type="Dopa-1", contrasts=c("MS","VEH")) {

  cap <- function(orig_values, cap, uncap=FALSE){
    if(!uncap){
      orig_values <- ifelse(abs(orig_values) < abs(cap), orig_values,
                            ifelse(orig_values >= cap, cap,
                                   ifelse(orig_values <= -cap, -cap, NA)))
    }else{
      orig_values <- ifelse(abs(orig_values) > abs(cap), orig_values,
                            ifelse(orig_values <= cap & orig_values > 0, cap,
                                   ifelse(orig_values >= -cap & orig_values < 0, -cap, NA)))
    }

    return(orig_values)
  }

  so <- subset(so, !!sym(grouping2)==cell_type & !!sym(grouping1) %in% contrasts)
  merge_de <- left_join(orig_de_df, rank_de_df, by="Gene.name", suffix=c("_orig","_ranked"))
  merge_de$p_val_adj <-  cap(merge_de$p_val_adj_ranked, cap=0.00000001, uncap=T)

  rank_data <- t(so@assays[[assay]][layer]) %>% as.data.frame()
  features <- colnames(rank_data)
  sample_mdata <- so@meta.data %>% dplyr::select("treatment", "genotype", "ct1")
  sample_mdata$barcode <- rownames(sample_mdata)
  rank_data$barcode <- rownames(rank_data)
  rank_data <- dplyr::left_join(rank_data, sample_mdata, by = "barcode")
  rank_data <- rank_data %>% dplyr::select(-c("barcode"))
  avg_rank_data_by_treatment <- rank_data %>%
    group_by(treatment, genotype, ct1) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>% t()


  colnames(avg_rank_data_by_treatment) <- contrasts
  avg_rank_data_by_treatment <- avg_rank_data_by_treatment[features,] %>% as.data.frame()
  avg_rank_data_by_treatment <- mutate_all(avg_rank_data_by_treatment, as.numeric)
  avg_rank_data_by_treatment$change <- -(avg_rank_data_by_treatment$MS - avg_rank_data_by_treatment$VEH)
  max_ranks <- t(rank_data[features]) %>% as.data.frame()
  max_ranks <- mutate_all(max_ranks, as.numeric) %>% as.matrix()
  max_ranks <- rowMaxs(max_ranks)
  avg_rank_data_by_treatment$max_rank <- max_ranks
  avg_rank_data_by_treatment$change_percentile <- avg_rank_data_by_treatment$change/avg_rank_data_by_treatment$max_rank
  change_pct_hist <- ggplot(avg_rank_data_by_treatment) +
    geom_histogram(aes(x=change_percentile), binwidth=0.01) + xlim(c(-0.1,0.1))
  #ggsave(change_pct_hist)
  avg_rank_data_by_treatment$Gene.name <- rownames(avg_rank_data_by_treatment)
  rank_de_df$Gene.name <- rownames(rank_de_df)
  merge_de <- merge(merge_de, avg_rank_data_by_treatment, by="Gene.name")

  # p_val_adj distribution plots
  p_p_orig <- merge_de %>% ggplot() +
    geom_histogram(aes(x=p_val_adj_orig), binwidth=0.1) +
    coord_cartesian(ylim = c(0, 1000), xlim = c(-0.1, 1.1)) +
    ggtitle(paste(cell_type, "Original Adj P Values"))
  print(p_p_orig)

  p_p <- merge_de %>% ggplot() +
    geom_histogram(aes(x=p_val_adj), binwidth=0.1) +
    coord_cartesian(ylim = c(0, 1000), xlim = c(-0.1, 1.1)) +
    ggtitle(paste(cell_type, "Adjusted Ranks Adj P Values"))
  print(p_p)

  # fold change distribution plots
  p_fc <- merge_de %>% ggplot() +
    geom_histogram(aes(x=avg_log2FC_orig), binwidth=0.1) +
    ggtitle(paste(cell_type, "Original Log2FC - all genes")) +
    theme(plot.title=element_text(size=22))
  print(p_fc)

  p_fc <- merge_de %>% ggplot() +
    geom_histogram(aes(x=change_percentile), binwidth=0.1) +
    ggtitle(paste(cell_type, "Rank shift Effect Size - all genes")) +
    theme(plot.title=element_text(size=22))
  print(p_fc)

  return(merge_de)
}
