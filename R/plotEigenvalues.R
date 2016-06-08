#' plot eigenvalues
#' \code{plotMultiEigenvalues} plot eigenvalues to determine number of communities by finding the elbow point
#'
#' @param Ens_list a list in which elements are numeric vectors representing eigenvalues.
#' @param mfrow A vector of the form \code{c(nr, nc)} passed to \code{\link{par}}.
#' @param mar plotting parameters with useful defaults (\code{\link{par}})
#' @param line plotting parameters with useful defaults (\code{\link{par}})
#' @param cex plotting parameters with useful defaults (\code{\link{par}})
#' @param ... further plotting parameters
#' @details
#' \code{plotMultiEigenvalues} plot multiple eigenvalue plots. The dark blue colored dots indicate eigenvalue greater than 0.
#' Each of the ensemble matrices is decomposed into eigenvalues which is used to determine appropriate number of communities.
#' Plotting out eigenvalues allow us to see where the elbow point is.
#' The curve starting from the elbow point flatten out. The number of points above (excluding) the elbow point indicates number of communities.
#'
#'
#' \code{mfrow} determines the arrangement of multiple plots. It takes the form of
#' \code{c(nr, nc)} with the first parameter being the number of rows and
#' the second parameter being the number of columns. When deciding parameters for mfrow,
#' one should take into considerations size of the plotting device and number of plots.
#' For example, there are 20 plots, mfrow can be set to \code{c(4, 5)} or \code{c(2, 10)}
#' depending on the size and shape of the plotting area.
#' @return a \code{pdf} file in the working directory containing all eigenvalue plots
#' @seealso \code{\link{plotCLUSTERS}}, \code{\link{getEnsList}}
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
#' \dontrun{
#' # for illustration only. skip CRAN check because it ran forever.
#' Ens_list <- getEnsList(Sim, temperatures, MaxIt = 1000, m = 5)
#' }
#' \dontshow{
#' # for CRAN check only
#' Ens_list <- getEnsList(Sim, temperatures, MaxIt = 5, m = 5)
#' }
#'
#' plotMultiEigenvalues(Ens_list = Ens_list, mfrow = c(10, 2), mar = c(1, 1, 1, 1))
#'
#' @export

plotMultiEigenvalues <- function(Ens_list,
                                 mfrow,
                                 mar = c(2, 2, 2, 2),
                                 line = -1.5,
                                 cex = 0.5, ...) {
   eigenvalue_list <- getEigenvalueList(Ens_list)
   op <- par(mfrow = mfrow,
            mar = mar) # set arrangement

  for (i in 1:length(eigenvalue_list)){
    Eigenvalues = eigenvalue_list[[i]]
    n = length(Eigenvalues)
    COLOR = rep("black",n)
    COLOR[which(Eigenvalues > 0)] = "dark blue"
    plot(Eigenvalues,
         type = "b",
         pch = 20,
         col = COLOR,
         cex = cex,
         ylab = "Normalized eigenvalues",
         main = paste0("Eigen-plot", " ", "Ens", i))
  }
  par(op)
}











plotEigenValue <- function(Eigenvalues, ...) {
  n = length(Eigenvalues)
  COLOR = rep("black",n)
  COLOR[which(Eigenvalues > 0)] = "red"
  plot(Eigenvalues,
       type = "b",
       pch = 20,
       col = COLOR,
       cex = 0.5,
       ylab = "Normalized eigenvalues",
       main = "Eigen-plot")
}
