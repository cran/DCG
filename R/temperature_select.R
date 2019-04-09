#' generate temperatures
#' \code{temperatureSample} generate tempatures based on either random or fixed intervals
#'
#' @param start a numeric vector of length 1, indicating the lowest temperature
#' @param end a numeric vector of length 1, indicating the highest temperature
#' @param n an integer between 10 to 30, indicating the number of temperatures (more explanations on what temperatures are).
#' @param method a character vector indicating the method used in selecting temperatures.
#' It should take either 'random' or 'fixedInterval', case-sensitive.
#' @return a numeric vector of length n representing temperatures sampled.
#'
#' @details In using random walks to find community structure, each normalized similarity matrix is evaluated at different temperatures.
#' This allows greater variations in the normalized similarity matrices.
#' It is recommended to try out 20 - 30 temperatures to allow for a thorough exploration of the matrices.
#' A range of temperatures which lead to stable community structures should be considered as reliable. The temperature in the middle of the range should be selected.
#'
#' @seealso \code{\link{getEnsList}}
#' @references
#' Fushing, H., & McAssey, M. P. (2010).
#' Time, temperature, and data cloud geometry.
#' Physical Review E, 82(6), 061110.
#'
#' Chen, C., & Fushing, H. (2012).
#' Multiscale community geometry in a network and its application.
#' Physical Review E, 86(4), 041120.
#'
#' Fushing, H., Wang, H., VanderWaal, K., McCowan, B., & Koehl, P. (2013).
#' Multi-scale clustering by building a robust and self correcting ultrametric topology on data points.
#' PloS one, 8(2), e56259.
#' @examples
#' symmetricMatrix <- as.symmetricAdjacencyMatrix(monkeyGrooming, weighted = TRUE, rule = "weak")
#' Sim <- as.SimilarityMatrix(symmetricMatrix)
#' temperatures <- temperatureSample(start = 0.01, end = 20, n = 20, method = 'random')
#' @export



temperatureSample <- function(start = 0.01, end = 20, n = 20, method = "random")
{
  if (n > 30 | n < 1){
    stop("n should be integers between 1 and 30.")
  }
  if (method == "random") {
    tem_seq <- seq(start, end, (end - start)/1000 + start)
    temperatures <- sample(tem_seq, n)
  } else if (method == "fixedInterval") {
    temperatures <- seq(start, end, (end - start)/n + start)
  } else {
    stop("method should be either 'random' or 'fixedInterval'.")
  }
  temperatures_sorted <- sort(temperatures)
  return(temperatures_sorted)
}

