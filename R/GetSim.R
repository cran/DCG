############################################################################
## For simplicity, suppose a distance matrix D is given.
## The similarity matrix W is calculated at each temperature T.
## The diagonal of W are all 0.
############################################################################
#' \code{GetSim} get similarity matrix from a distance matrix
#' @param D A distance matrix
#' @param T Temperature. \code{\link{temperatureSample}}
#' @details the similarity matrix is calculated at each temperature \code{T}.
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
#' @export

GetSim <- function(D, T) {
  W <- exp(-D/T)
  diag(W) <- 0
  return(W)
}
