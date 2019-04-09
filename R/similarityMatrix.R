#' convert to a symmetric adjacency matrix
#'
#' \code{as.symmetricAdjacencyMatrix} convert an edgelist or a raw matrix to a symmetric adjacency matrix.
#' @param Data either a dataframe or a matrix, representing raw interactions
#' using either an edgelist or a matrix.
#' Frequency of interactions for each dyad can be represented either
#' by multiple occurrences of the dyad for a 2-column edgelist, or
#' by a third column specifying the frequency of the interaction
#' for a 3-column edgelist.
#' @param weighted If the edgelist is a 3-column edgelist in which weight was
#' specified by frequency, use \code{weighted = TRUE}.
#' @param rule a character vector of length 1, being one of "\code{weak}",
#' "\code{strong}", "\code{upper}", or "\code{lower}".
#' Ways of symmetrizing the matrix. See details for more information.
#' @return a named matrix with the \code{[i,j]}th entry equal to the
#' number of times \code{i} grooms \code{j}.
#'
#' @details There are ways of symmetrizing a matrix.
#' The "\code{weak}" rule symmetrize the matrix by building an edge
#' between nodes \code{[i, j]} and \code{[j, i]} if there is an edge
#' either from \code{i} to \code{j} OR from \code{j} to \code{i}.
#' The "\code{strong}" rule symmetrize the matrix by building an edge
#' between nodes \code{[i, j]} and \code{[j, i]} if there is an edge
#' BOTH from \code{i} to \code{j} AND from \code{j} to \code{i}.
#' The "\code{upper}" and the "\code{lower}" rule symmetrize the matrix
#' by using the "\code{upper}" or the "\code{lower}" triangle respectively.
#'
#'
#' Note, when using a 3-column edgelist (e.g. a weighted edgelist) to
#' represent raw interactions, each dyad must be unique.
#' If more than one rows are found with the same Initiator and recipient,
#' sum of the frequencies will be taken to represent the freqency of
#' interactions between this unique dyad.
#' A warning message will prompt your attention to the accuracy of your
#' raw data when duplicated dyads were found in a three-column edgelist.
#' @examples
#' symmetricMatrix <- as.symmetricAdjacencyMatrix(monkeyGrooming, weighted = TRUE, rule = "weak")
#'
#' @export

as.symmetricAdjacencyMatrix = function(Data, weighted = FALSE, rule = "weak"){
  if (ncol(Data) > 3 & ncol(Data) != nrow(Data)) {
    stop("check your raw data: A edgelist should be of either 2 or 3 columns.
         If it is a matrix, the column number should be equal to row number.")
  }
  # for matrix
  if (ncol(Data) == nrow(Data)){
    # if values on diagonal are not all zeros,
    #   convert to zero, return warnings.
    mat_raw <- as.matrix(Data)
    if (any(diag(mat_raw != 0))) {
      index <- which(diag(mat_raw) != 0)
      diag(mat_raw)[index] <- 0
      warning(paste("check your raw matrix at Row",
                    paste(index, collapse = ","), "and column",
                    paste(index, collapse = ","), ";
                    Non-zero values on diagonal was converted to zeros."))
    }
    subjects <- sort(colnames(mat_raw))
    mat <- mat_raw[subjects, subjects]
  } else {
    mat <- edgelisttomatrix(Data, weighted)
  }
  # symmetrize the matrix
  symmetricMat <- symmetrizeWeightedMatrix(mat, rul = rule)
  return(symmetricMat)
}


#' Convert a matrix to a similarity matrix.
#' \code{as.SimilarityMatrix} convert an adjacency matrix to a similarity matrix.
#' @param mat a symmetric adjacency matrix
#'
#' @return a similarity matrix.
#'
#' @examples
#' symmetricMatrix <- as.symmetricAdjacencyMatrix(monkeyGrooming, weighted = TRUE, rule = "weak")
#' similarityMatrix <- as.SimilarityMatrix(symmetricMatrix)
#'
#' @export


as.SimilarityMatrix <- function(mat){
  # if asymmetric, warning
  if (!isSymmetric(mat)) {
    warning("Asymmetric matrix should not be used.
      DCG is designed to find clusters for undirected network.
      You'll be responsible for interpreting the results on assymetric matrix,
      and use it at your own risks.")
  }
  maxraw <- max(mat)
  Sim<- mat/maxraw
  class(Sim) = c("similarityMatrix", "matrix")
  return(Sim)
}



# symmetrize a weighted matrix
symmetrizeWeightedMatrix <- function(matrix, rul = "weak"){
  # initializing upper and lower matrix
  upperMatrix <- lowerMatrix <- strongMatrix <- matrix
  # keep only the upper or lower triangle in each
  upperMatrix[lower.tri(matrix, diag = TRUE)] <- 0
  lowerMatrix[upper.tri(matrix, diag = TRUE)] <- 0
  tLowerMatrix <- t(lowerMatrix)
  if (rul == "weak") {
    # transpose lower and add it to upper
    upperSymmetricMatrix <- tLowerMatrix  + upperMatrix
    # populate lower symmetric matrix with the
    #    transposed upper symmetric matrix
    symmetricMatrix <- upperSymmetricMatrix + t(upperSymmetricMatrix)
  } else if (rul == "strong") {
    strongMatrix[] <- 0 # assign 0 to all elements
    strongMatrix[upperMatrix < tLowerMatrix] <-
      upperMatrix[upperMatrix < tLowerMatrix]
    strongMatrix[tLowerMatrix < upperMatrix] <-
      tLowerMatrix[tLowerMatrix < upperMatrix]
    strongMatrix[upperMatrix == tLowerMatrix] <-
      upperMatrix[upperMatrix == tLowerMatrix]
    symmetricMatrix <- strongMatrix + t(strongMatrix)
  } else if (rul == "upper") {
    symmetricMatrix <- upperMatrix + t(upperMatrix)
  } else if (rul == "lower") {
    symmetricMatrix <- lowerMatrix + t(lowerMatrix)
  } else {
    stop("Please use one of the following symmetrizing rules:
         'weak', 'strong', 'upper', 'lower'.")
  }
  return(symmetricMatrix)
  }


# transform an edgelist into a matrix
#
# @param edgelist a 2-column (or 3-column for weighted edgelist) dataframe/matrix of edges. The Initiator is in the 1st column by default. For weighted edgelist, the third column should be the weight.
# @param weighted If the edgelist is a 3-column weighted edgelist, use \code{weighted = TRUE}.
# @return a named matrix with \code{[i,j]}th entry equal to the number of times \code{i} initiated interactions over \code{j}.
# It is the matrix representation of the edgelist.
#
# @seealso \code{\link{conductance}}

edgelisttomatrix <- function(edgelist, weighted = FALSE) {

  if (ncol(edgelist) > 3) {
    stop("edgelist should be of 2 column, or 3-column for weighted edgelist")
  }


  if (any(edgelist[,1] == edgelist[,2])) {
    rowIndex <- which(edgelist[,1] == edgelist[,2])
    edgelist <- edgelist[-rowIndex, ]
    warning(paste("check your raw data at row number", paste(rowIndex, collapse = ","), ". The initiator and the recipient are the same. These data were removed"))
  }

  subjects = unique(sort(as.matrix(edgelist[,1:2]))) # work better for IDs of character
  # subjects = sort(unique(c(edgelist[,1], edgelist[,2])))
  N = length(subjects)
  if (N > 10000){
    stop("No more than 10000 unique subjects.")
  }

  mat = matrix(0, N, N)


  if (weighted == TRUE){

    if (ncol(edgelist) != 3){
      stop("Input a matrix or dataframe with three columns, with the third column being Frequency of the interaction")
    }

    if (anyDuplicated(edgelist[,1:2]) != 0) {
      warning(
        "dyads in the weighted edgelist are not unique; the sum of frequencies is taken for duplicated rows."
      )
      edgelist <- sumDuplicate(edgelist)
    }


    # transform the weighted edgelist into a matirx

    for(i in 1:nrow(edgelist)){
      subject1 = which(subjects == edgelist[i,1])
      subject2 = which(subjects == edgelist[i,2])
      mat[subject1, subject2] = edgelist[i, 3]
    }

  } else {

    if (ncol(edgelist) != 2){
      stop("edgelist should be a dataframe or matrix of two columns. If it is a weighted edgelist, it should be a matrix or dataframe of 3 columns and use the argument 'weighted = TRUE'")
    }

    for(i in 1:nrow(edgelist)){
      subject1 = which(subjects == edgelist[i,1])
      subject2 = which(subjects == edgelist[i,2])
      mat[subject1, subject2] = mat[subject1, subject2] + 1
    }
  }

  rownames(mat) = subjects
  colnames(mat) = subjects

  return(mat)
}




#### internal functions


sumDuplicate <- function(weightedEdgelist) {
  uniqueEdgelist <- unique(weightedEdgelist[,1:2])
  for (i in 1:nrow(uniqueEdgelist)){
    uniqueEdgelist[i,3] <-
      sum(
        weightedEdgelist[
          match.2coldf(weightedEdgelist[,1:2],  uniqueEdgelist[i,]),
          3])
  }
  names(uniqueEdgelist) <- names(weightedEdgelist)
  return(uniqueEdgelist)
}


match.2coldf <- function(dataframe, values) {
  # dataframe should be of two columns
  # values should be a vector of length 2, or a row of dataframe of two columns
  rowIndex <- intersect(which(dataframe[,1] == values[[1]]),
                        which(dataframe[,2] == values[[2]]))
  return(rowIndex)
}
