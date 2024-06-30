Recorrecting Effect-size and Significance with Rank-based Evaluation of Confounds in Tests

Installation instructions
```r
# Install devtools package if you haven't already
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install ResRect from GitHub
devtools::install_github("ethanfenton/resrect")
```

```{r}
library(resrect)

# Calculate ranks
## 'so' is a Seurat Object
so <- calculateRanks(so)

# Plot unadjusted ranks
plotRanks(so, output_folder = "path/to/outputs/")

# Adjust ranks
so <- adjustRanks(so)

# Plot adjusted ranks
plotAdjustedRanks(so, output_folder = "path/to/outputs/")
```
