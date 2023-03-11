
#' Population Heterogeneity
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#'
#' @return Numeric vector of gene-wise population heterogeneities
#' @export
#'
#' @examples
#' set.seed(1)
#' CountsMatrix <- Matrix::sparseMatrix(i = sample(10, 100, replace = TRUE),
#' j = sample(10, 100, replace = TRUE),
#' x = 1L,
#' dims = c(10,10))
#' Ig <- Population(CountsMatrix)
Population <- function(CountsMatrix) {

  Total <- Matrix::rowSums(CountsMatrix)
  N <- ncol(CountsMatrix)

  transposeCounts <- Matrix::t(CountsMatrix)

  Indices <- length(transposeCounts@p)-1
  Pop <- vector("numeric", length(Indices))

  for (ind in 1:Indices) {
    freqshrink <- GetFreqShrink(transposeCounts, ind, N, Total)
    LogNfreq <- log(N*freqshrink)
    LogNfreq[LogNfreq == -Inf] <- 0
    Pop[ind] <- t(freqshrink) %*% LogNfreq
  }


  Pop[is.infinite(Pop)] <- 0

  names(Pop) <- colnames(transposeCounts)

  Pop
}

#' Feature Selection by Population Heterogeneity
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#' @param nGenes Number of genes selected
#' @param minHeterogeneity Threshold genes below specified level of heterogeneity (as measured by scEC::Population)
#' @param minCounts Minimum number of transcripts per gene
#'
#' @return Vector of top genes by population heterogeneity
#' @export
#'
FeatureSelection <- function(CountsMatrix, nGenes = 500, minHeterogeneity = 0, minCounts = 100){

  Exp <- rownames(CountsMatrix)[Matrix::rowSums(CountsMatrix) >= minCounts]

  Div <- Population(CountsMatrix[which(rownames(CountsMatrix) %in% Exp),])
  Div <- Div[Div > minHeterogeneity]
  GOI <- names(sort(Div, decreasing = T)[1:nGenes])

  GOI
}
