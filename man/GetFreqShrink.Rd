% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Normalisation.R
\name{GetFreqShrink}
\alias{GetFreqShrink}
\title{James-Stein Frequency Estimator}
\usage{
GetFreqShrink(transposeCounts, ind, N, Total)
}
\arguments{
\item{transposeCounts}{Transposed sparse count matrix}

\item{ind}{Integer indicating chosen gene (row number in count matrix)}

\item{N}{N Number of cells}

\item{Total}{Total Integer of total counts per cell}
}
\value{
Numeric vector of shrinkage cell frequencies
}
\description{
James-Stein Frequency Estimator
}
