
#' Differential Expression by Inter-Type Heterogeneity
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#' @param Identity Factor of cell identities
#'
#' @return Numeric vector of gene-wise inter-type heterogeneities
#' @export
#'
#' @examples
#' set.seed(1)
#' CountsMatrix <- Matrix::sparseMatrix(i = sample(10, 100, replace = TRUE),
#' j = sample(10, 100, replace = TRUE),
#' x = 1L,
#' dims = c(10,10))
#' Clustering <- factor(sample(2, 10, replace = TRUE))
#' DifferentialExpression(CountsMatrix, Clustering)
DifferentialExpression <- function(CountsMatrix, Identity) {

  if(length(Identity) != ncol(CountsMatrix)){
    stop("Inconsistent number of cells between objects:\n\tlength(Identity) != ncol(CountsMatrix)")
  }

  Total <- Matrix::rowSums(CountsMatrix)
  N <- ncol(CountsMatrix)

  Ng <- as.vector(table(Identity))

  transposeCounts <- Matrix::t(CountsMatrix)

  Indices <- length(transposeCounts@p)-1
  InterType <- vector("numeric", length(Indices))

  for (ind in 1:Indices) {
    freqshrink <- GetFreqShrink(transposeCounts, ind, N, Total)
    groupedfreqshrink <- tapply(freqshrink, Identity, sum)

    NonZero <- which(groupedfreqshrink != 0)

    InterType[ind] <- t(groupedfreqshrink[NonZero]) %*% log(N*groupedfreqshrink[NonZero] /  Ng[NonZero])
  }

  InterType[is.infinite(InterType)] <- 0

  names(InterType) <- rownames(CountsMatrix)

  InterType
}
