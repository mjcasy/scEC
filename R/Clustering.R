
#' Single-Cell Entropic Clustering
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#' @param MaxClus Maximum number of clusters
#' @param Multistart Number of restarts at each step
#' @param Seed Seed
#'
#' @return Matrix where each column is a vector of integer identities. The Nth
#' column encodes N clusters.
#' @export
#'
#' @examples
Cluster <- function(CountsMatrix, MaxClus, Multistart = 5, Seed){

  if(!missing(Seed)){
    set.seed(Seed)
    reticulate::py_set_seed(seed = Seed)
  }

  G <- dim(CountsMatrix)[1]
  N <- dim(CountsMatrix)[2]

  FullFreq <- GetFreq(CountsMatrix)

  Ident <- matrix(NA, nrow = N, ncol = MaxClus)
  Ident[,1] <- 0

  IdentSplit <- 0
  BoolNewIdent <- matrix(FALSE, nrow = N, ncol = 2*MaxClus)

  mu <- PyFunc$multiStartSplitCell(freq = FullFreq, multistart = Multistart)

  newIdent <- apply(mu, 1, which.max) - 1
  newBool <- newIdent == 1

  BoolNewIdent[,1] <- newBool

  Score <- sum(PyFunc$intertype(FullFreq, newIdent)) - sum(PyFunc$intertype(FullFreq, Ident[,1]))

  k <- 2

  for(i in 2:MaxClus){

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

        Score[k] <- sum(PyFunc$intertype(FullFreq, newIdent)) -
          sum(PyFunc$intertype(FullFreq, Ident[,i]))


      } else {

        IdentSplit[k] <- j
        BoolNewIdent[BoolCells, k] <- FALSE
        Score[k] <- 0

      }

      k <- k+1
    }
  }

  Ident
}

#' Comparitive Gain in Inter-Type Heterogeneity
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#' @param IdentMat Matrix output of scEC
#'
#' @return Observed versus expected gains in inter-type heterogeneity
#' @export
#'
#' @examples
EntropyGain <- function(CountsMatrix, IdentMat){

  reticulate::source_python("~/OneDrive - University of Southampton/Entropy/scIC/Scripts_for_pkg/Heterogeneity.py")

  MaxClus <- ncol(IdentMat)
  Freq <- GetFreq(CountsMatrix)

  Pop <- sum(PyFunc$pop(freq = Freq))
  Unif <- (Pop / log(ncol(CountsMatrix)))*log(1:MaxClus)
  Exp <- Unif[2:MaxClus] - Unif[1:(MaxClus-1)]

  Inter <- c()
  for(i in 1:MaxClus){
    Inter[i] <- sum(PyFunc$intertype(freq = Freq, ident = IdentMat[,i]))
  }
  Obs <- Inter[2:MaxClus] - Inter[1:(MaxClus-1)]

  ClusNum <- 2:MaxClus
  Gain <- cbind(ClusNum, Obs, Exp)

  Gain
}
