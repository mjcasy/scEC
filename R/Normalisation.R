#' James-Stein Frequency Estimator
#'
#' @param transposeCounts Transposed sparse count matrix
#' @param ind Integer indicating chosen gene (row number in count matrix)
#' @param N N Number of cells
#' @param Total Total Integer of total counts per cell
#'
#' @return Numeric vector of shrinkage cell frequencies
#'
GetFreqShrink <- function(transposeCounts, ind, N, Total){

  tk <- rep(1/N, N)
  tkadj <- tk

  count <- transposeCounts@x[(transposeCounts@p[ind]+1) : transposeCounts@p[ind+1]]
  elements <- transposeCounts@i[(transposeCounts@p[ind]+1) : transposeCounts@p[ind+1]]+1

  freq <- count / Total[ind]

  num <- 1 - sum(freq^2)
  tkadj[elements] <- tkadj[elements] - freq
  den <- (sum(count) - 1)*sum(tkadj^2)
  lambda <- num/den
  lambda[lambda > 1] <- 1

  freqshrink <- lambda*tk
  freqshrink[elements]  <- freqshrink[elements] + (1 - lambda)*freq

  freqshrink
}


#' Normalise Frequencies
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#'
#' @return Feature x cell dense matrix of frequencies
#'
GetFreq <- function(CountsMatrix){

  Total <- Matrix::rowSums(CountsMatrix)
  transposeCounts <- Matrix::t(CountsMatrix)
  Indices <- length(transposeCounts@p)-1
  N <- dim(CountsMatrix)[2]

  freq <- matrix(NA, nrow = nrow(CountsMatrix), ncol = N)
  for (ind in 1:Indices) {
    freq[ind,] <- GetFreqShrink(transposeCounts, ind, N, Total)
  }

  freq
}
