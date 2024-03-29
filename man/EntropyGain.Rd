% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Clustering.R
\name{EntropyGain}
\alias{EntropyGain}
\title{Comparitive Gain in Inter-Type Heterogeneity}
\usage{
EntropyGain(CountsMatrix, IdentMat)
}
\arguments{
\item{CountsMatrix}{Feature x cell sparse counts matrix of class dgCMatrix}

\item{IdentMat}{Matrix output of scEC (Greedy = T)}
}
\value{
Observed versus expected gains in inter-type heterogeneity
}
\description{
Comparitive Gain in Inter-Type Heterogeneity
}
