

#' Feature Selection by Inter-cluster Heterogeneity
#'
#' @param ReferenceCountsMatrix Reference count matrix; feature x cell sparse counts matrix of class dgCMatrix
#' @param ReferenceID Factor of reference cell identities
#' @param minCounts Minimum number of transcripts per gene
#' @param nGenes Number of genes selected
#'
#' @return Vector of top genes by inter-cluster heterogeneity
#'
#' @examples
MapFeatureSelection <- function(ReferenceCountsMatrix, ReferenceID, minCounts, nGenes){

  Exp <- rownames(ReferenceCountsMatrix)[Matrix::rowSums(ReferenceCountsMatrix) >= minCounts]

  Div <- DifferentialExpression(ReferenceCountsMatrix[Exp,], ReferenceID)
  GOI <- names(sort(Div, decreasing = T)[1:nGenes])

  GOI
}

#' Mapping of Cells Onto Reference Cell Identities
#'
#' @param MapCountsMatrix Count matrix of cells to be mapped
#' @param ReferenceCountsMatrix Count matrix of cells with known identities
#' @param ReferenceID Factor of reference cell identities
#' @param minCounts Minimum number of transcripts per gene
#' @param nGenes Number of genes selected
#' @param Seed Seed set for both R and Python components
#'
#' @return Mapped cellular identities
#' @export
#'
#' @examples
Map <- function(MapCountsMatrix, ReferenceCountsMatrix, ReferenceID, minCounts = 100, nGenes = 1000, Seed){

  if(!missing(Seed)){
    set.seed(Seed)
    reticulate::py_set_seed(seed = Seed)
  }

  RefN <- ncol(ReferenceCountsMatrix)
  MapN <- ncol(MapCountsMatrix)
  N <- RefN + MapN

  GOI <- MapFeatureSelection(ReferenceCountsMatrix, ReferenceID, minCounts = minCounts, nGenes = nGenes)
  ReferenceCountsMatrix <- ReferenceCountsMatrix[GOI,]
  MapCountsMatrix <- MapCountsMatrix[GOI,]

  RefFreq <- GetFreq(ReferenceCountsMatrix)
  MapFreq <- GetFreq(MapCountsMatrix)

  FullFreq <- cbind((RefN/N) * RefFreq, (MapN/N) * MapFreq)
  RefID <- as.numeric(ReferenceID)-1

  mu <- PyFunc$meld(freq = FullFreq, refID = RefID)
  Ident <- apply(mu, 1, which.max) - 1

  levels(ReferenceID)[Ident+1]

}
