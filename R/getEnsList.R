#' generating a list of ensemble matrices based on the similarity matrix and temperatures
#'
#' \code{getEnsList} get ensemble matrices from given similarity matrix at all temperatures
#'
#' @param simMat a similarity matrix
#' @param temperatures temperatures selected
#' @param MaxIt number of iterations for regulated random walks
#' @param m maxiumnum number of time a node can be visited during random walks
#'
#' @details This step is crucial in finding community structure based on the similarity matrix of the social network.
#' For each \code{temperatures}, the similarity matrix was taken to the power of \code{temperature} as saved as a new similarity matrix.
#' This allows the random walk to explore the similarity matrix at various variations.
#' Random walks are then performed in similarity matrices of various temperatures.
#' In order to prevent random walks being stucked in a locale, the parameter \code{m} was set (to \code{5} by default) to remove a node after \code{m} times of visits of the node.
#' An ensemble matrix is generated at each temperature in which values represent likelihood of two nodes being in the same community.
#'
#'
#' @return a list of ensemble matrices
#'
#' @examples
#' symmetricMatrix <- as.symmetricAdjacencyMatrix(monkeyGrooming, weighted = TRUE, rule = "weak")
#' Sim <- as.SimilarityMatrix(symmetricMatrix)
#' temperatures <- temperatureSample(start = 0.01, end = 20, n = 20, method = 'random')
#' \dontrun{
#' # Note: It takes a while to run the getEnsList example.
#' Ens_list <- getEnsList(Sim, temperatures, MaxIt = 1000, m = 5)
#' }
#' \dontshow{
#' # for CRAN check only
#' Ens_list <- getEnsList(Sim, temperatures, MaxIt = 5, m = 5)
#' }
#'
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
#'
#' @export


getEnsList <- function(simMat, temperatures, MaxIt = 1000, m = 5) {
  if (!isSymmetric(simMat)) {
    warning("Asymmetric matrix should not be used.
            DCG is designed to find clusters for undirected network.
            You'll be responsible for interpreting the results on assymetric matrix,
            and use it at your own risks.")
  }
  ens_list <- lapply(temperatures, function(x) getEns(simMat, x, MaxIt = MaxIt, m = m))
  return(ens_list)
}

#' generate ensemble matrix
#' \code{getEns} get ensemble matrix from given similarity matrix and temperature
#'
#' @param simMat a similarity matrix
#' @param temperature a numeric vector of length 1, indicating the temperature used to transform the similarity matrix to ensemble matrix
#' @param MaxIt number of iterations for regulated random walks
#' @param m maxiumnum number of time a node can be visited during random walks
#' @return a matrix.
#' @details This function involves two steps.
#' It first generate similarity matrices of different variances
#' by taking the raw similarity matrix to the power of each
#' temperature. Then it called the function \code{EstClust} to perform random walks in the network to identify clusters.


getEns <- function(simMat, temperature, MaxIt = 1000, m = 5) {
  data_t <- simMat^temperature
  Ens <- EstClust(data_t, MaxIt, m)
}

