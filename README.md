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
so <- calculateRanks(so, assay = "SCT", layer = "data", grouping1 = "sample", grouping2 = "ct1", pct_exp_cutoff = 10, output_folder = "outputs_v1/DE/single_cell/corrected_ranks/resrect_2/", group_levels = NA)

## Plotting a gene of interest
gene_of_interest = "Oprm1"
plotGeneRanks(so, gene=gene_of_interest, output_folder="outputs_v1/DE/single_cell/corrected_ranks/resrect_2/", group_levels = NA)

## Use your de test of choice on the adj_ranks assay. ensure columns include Gene.name, p_adj, and avg_log2FoldChange
bc.MS = so@meta.data$barcode[which(so@meta.data$treatment=="MS" & so@meta.data$genotype=="wt" & so@meta.data$ct1=="Glut-2")]
bc.VEH = so@meta.data$barcode[which(so@meta.data$treatment=="VEH" & so@meta.data$genotype=="wt" & so@meta.data$ct1=="Glut-2")]

de <- FindMarkers(so, cells.1=bc.MS, cells.2=bc.VEH)
plotVolcanos(so, de=my_de)

#
```
