# Recorrecting Effect-size and Significance with Rank-based Evaluation of Confounds in Tests

## Installation instructions
```r
# Install devtools package if you haven't already
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install resrect from GitHub
devtools::install_github("ethanfenton/resrect")
```

## Example Usage
```{r}
library(resrect)

# 
## 'so' is a Seurat Object
so <- calculateRanks(so, assay = "SCT", layer = "data", grouping1 = "sample", grouping2 = "ct1", pct_exp_cutoff = 10, output_folder = "outputs_v1/DE/single_cell/corrected_ranks/resrect_mini", group_levels = NA)

## Plotting a gene of interest
plotGeneRanks(so, plot_genes =c("Actb","Scn9a","Fos","Slc17a6"), ct1="Dopa-1", output_folder = "outputs_v1/DE/single_cell/corrected_ranks/resrect_mini/single_gene_plots/")

# Use DE of your choice on rank corrected data
dopa_1_adj_ranks_wt_de_df <- FindMarkers(subset(sud.ranked.so, genotype=="WT" & ct1=="Dopa-1"), ident.1="MS", ident.2="VEH",group.by ="treatment", assay = "adj_ranks", slot="counts", logfc.threshold = 0, min.pct = 0, test.use="wilcox")
dopa_1_adj_ranks_wt_de_df$Gene.name <- rownames(dopa_1_adj_ranks_wt_de_df)

# also do DE on original data for comparrision purposes
dopa_1_orig_wt_de_df <- FindMarkers(subset(sud.ranked.so, genotype=="WT" & ct1=="Dopa-1"), ident.1="MS", ident.2="VEH",group.by ="treatment", assay = "SCT", slot="data", logfc.threshold = 0, min.pct = 0, test.use="wilcox")
dopa_1_orig_wt_de_df$Gene.name <- rownames(dopa_1_orig_wt_de_df)

subset.so <- subset(sud.ranked.so, ct1=="Dopa-1" & treatment %in% c("MS", "VEH") & genotype=="WT")

# Add in rank shift as a metric akin to fold change
dopa_1_wt_de_df <- rank_shift_effect_size(rank_de_df=dopa_1_adj_ranks_wt_de_df, 
                                          orig_de_df = dopa_1_orig_wt_de_df,
                                          so=subset.so,
                                          assay = "adj_ranks", layer = "counts",
                                          grouping1 = "treatment",  grouping2="ct1",
                                          cell_type="Dopa-1", contrasts=c("MS","VEH"))
                                          
# Make volcano plots
# original volcano dopa-1 wt 
dopa_1_wt_de_df$avg_log2FC <- dopa_1_wt_de_df$avg_log2FC_orig
dopa_1_wt_de_df$p_val_adj <- dopa_1_wt_de_df$p_val_adj_orig
labels <- dopa_1_wt_de_df %>% dplyr::filter(p_val_adj_orig < 0.05)
p_volc_orig_dopa_1 <- plot_volcano3(dopa_1_wt_de_df, p_cutoff = 0.05, wilcoxon=T, my_x_lab = "avg L2FC (MS vs VEH)", p_adj = T, labels = labels) + ggtitle("Dopa-1 Original wilcox WT DE") + theme(plot.title=element_text(size=42)) #+ xlim(c(-0.2,0.2))
p_volc_orig_dopa_1

# rank adj volcano dopa-1 wt 
dopa_1_wt_de_df$avg_log2FC <- dopa_1_wt_de_df$change_percentile
dopa_1_wt_de_df$p_val_adj <- dopa_1_wt_de_df$p_val_adj_ranked
labels <- dopa_1_wt_de_df %>% dplyr::filter(p_val_adj < 0.05)
p_volc_rank_adj_dopa_1 <- plot_volcano3(dopa_1_wt_de_df, p_cutoff = 0.05, wilcoxon=T, my_x_lab = "avg rank percentile change (MS vs VEH)", p_adj = T, labels = labels) + ggtitle("Dopa-1 Rank Adj wilcox WT DE") + theme(plot.title=element_text(size=42)) #+ xlim(c(-0.2,0.2))
p_volc_rank_adj_dopa_1


#
```
