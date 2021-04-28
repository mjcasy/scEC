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
