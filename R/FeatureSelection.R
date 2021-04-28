
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

FeatureSelection <- function(CountsMatrix, minCounts = 100, nGenes = 500){

  Exp <- rownames(CountsMatrix)[Matrix::rowSums(CountsMatrix) > minCounts]

  Div <- Population(CountsMatrix[Exp,])
  GOI <- names(sort(Div, decreasing = T)[1:nGenes])

  GOI
}
