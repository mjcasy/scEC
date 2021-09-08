
#' Single-Cell Entropic Clustering
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#' @param numClus Number of clusters (or maximum number of clusters for greedy)
#' @param Multistart Number of restarts at each step
#' @param Greedy Boolean. Should greedy algorithm be used.
#' @param Seed Seed
#'
#' @return When Greedy = F, vector of integer identities. When Greedy = T,
#' matrix where each column is a vector of integer identities. The Nth column
#' encodes N clusters.
#' @export
#'
#' @examples
Cluster <- function(CountsMatrix, numClus, Multistart = 5, Greedy = F, Seed){

  if(!missing(Seed)){
    set.seed(Seed)
    reticulate::py_set_seed(seed = Seed)
  }

  G <- dim(CountsMatrix)[1]
  N <- dim(CountsMatrix)[2]

  FullFreq <- GetFreq(CountsMatrix)

  if(Greedy == F){
    mu <- PyFunc$multiStartClusterCells(freq = FullFreq, numClusters = numClus, multistart = Multistart)
    Ident <- apply(mu, 1, which.max) - 1

  } else if(Greedy == T){

    Ident <- matrix(NA, nrow = N, ncol = numClus)
    Ident[,1] <- 0

    IdentSplit <- 0
    BoolNewIdent <- matrix(FALSE, nrow = N, ncol = 2*numClus)

    mu <- PyFunc$multiStartSplitCell(freq = FullFreq, multistart = Multistart)

    newIdent <- apply(mu, 1, which.max) - 1
    newBool <- newIdent == 1

    BoolNewIdent[,1] <- newBool

    Score <- sum(PyFunc$intercluster(FullFreq, newIdent)) - sum(PyFunc$intercluster(FullFreq, Ident[,1]))

    k <- 2

    for(i in 2:numClus){

      if(all(Score == 0)){
        message("Cellular Equivalence Reached")
        break()
      }

      Choosen <- which.max(Score)
      Ident[,i] <- Ident[,(i-1)]
      Ident[BoolNewIdent[,Choosen], i] <- i-1

      Idents2Split <- c(IdentSplit[Choosen], i-1)

      Score[Choosen] <- 0

      for(j in Idents2Split){
        BoolCells <- (Ident[,i] == j)

        Freq <- FullFreq[, BoolCells, drop = FALSE]
        Freq <- sweep(Freq, 1, rowSums(Freq), '/')

        if(sum(BoolCells) > 1 & length(table(Freq)) > 1){

          mu <- PyFunc$multiStartSplitCell(freq = Freq, multistart = Multistart)

          newIdent <- apply(mu, 1, which.max) - 1
          newBool <- newIdent == 1

          IdentSplit[k] <- j
          BoolNewIdent[BoolCells, k] <- newBool

          newIdent <- Ident[,i]
          newIdent[BoolNewIdent[,k]] <- i

          Score[k] <- sum(PyFunc$intercluster(FullFreq, newIdent)) -
            sum(PyFunc$intercluster(FullFreq, Ident[,i]))


        } else {

          IdentSplit[k] <- j
          BoolNewIdent[BoolCells, k] <- FALSE
          Score[k] <- 0

        }

        k <- k+1
      }
    }
  }

  Ident
}

#' Comparitive Gain in Inter-Type Heterogeneity
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#' @param IdentMat Matrix output of scEC (Greedy = T)
#'
#' @return Observed versus expected gains in inter-type heterogeneity
#' @export
#'
#' @examples
EntropyGain <- function(CountsMatrix, IdentMat){

  numClus <- ncol(IdentMat)
  Freq <- GetFreq(CountsMatrix)

  Pop <- sum(PyFunc$pop(freq = Freq))
  Unif <- (Pop / log(ncol(CountsMatrix)))*log(1:numClus)
  Exp <- Unif[2:numClus] - Unif[1:(numClus-1)]

  Inter <- c()
  for(i in 1:numClus){
    Inter[i] <- sum(PyFunc$intercluster(freq = Freq, ident = IdentMat[,i]))
  }
  Obs <- Inter[2:numClus] - Inter[1:(numClus-1)]

  ClusNum <- 2:numClus
  Gain <- cbind(ClusNum, Obs, Exp)

  Gain
}
