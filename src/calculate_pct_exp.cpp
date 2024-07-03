#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix calculatePctExp(NumericMatrix data, CharacterVector cellTypes, CharacterVector genes) {
  int nCells = data.nrow();
  CharacterVector uniqueCellTypes = sort_unique(cellTypes);
  int nCellTypes = uniqueCellTypes.size();
  int nGenes = genes.size();

  NumericMatrix pctExp(nCellTypes, nGenes);
  colnames(pctExp) = genes;
  rownames(pctExp) = uniqueCellTypes;

  for (int i = 0; i < nCellTypes; i++) {
    for (int j = 0; j < nGenes; j++) {
      int count = 0;
      int total = 0;
      for (int k = 0; k < nCells; k++) {
        if (cellTypes[k] == uniqueCellTypes[i]) {
          total++;
          if (!NumericVector::is_na(data(k, j)) && data(k, j) > 0) {
            count++;
          }
        }
      }
      pctExp(i, j) = (double)count / total * 100;
    }
  }
  return pctExp;
}
