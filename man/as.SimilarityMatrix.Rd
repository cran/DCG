% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/similarityMatrix.R
\name{as.SimilarityMatrix}
\alias{as.SimilarityMatrix}
\title{Convert a matrix to a similarity matrix.
\code{as.SimilarityMatrix} convert an adjacency matrix to a similarity matrix.}
\usage{
as.SimilarityMatrix(mat)
}
\arguments{
\item{mat}{a symmetric adjacency matrix}
}
\value{
a similarity matrix.
}
\description{
Convert a matrix to a similarity matrix.
\code{as.SimilarityMatrix} convert an adjacency matrix to a similarity matrix.
}
\examples{
symmetricMatrix <- as.symmetricAdjacencyMatrix(monkeyGrooming, weighted = TRUE, rule = "weak")
similarityMatrix <- as.SimilarityMatrix(symmetricMatrix)

}
