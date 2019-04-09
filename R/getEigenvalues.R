#' generate eigenvalues for all ensemble matrices
#' \code{getEigenvalueList} get eigenvalues from ensemble matrices
#'
#' @param EnsList a list of ensemble matrices
#'
#' @return a list of \code{eigenvalues} for each of the ensemble matrix in the ensemble matrices list.
# #' @examples
# #' Sim <- as.simMat(monkeyGrooming, weighted = TRUE)
# #' temperatures <- temperatureSample(start = 0.01, end = 20, n = 20, method = 'random')
# #' \dontrun{
# #' # for illustration only. skip CRAN check because it ran forever.
# #' Ens_list <- getEnsList(Sim, temperatures, MaxIt = 1000, m = 5)
# #' }
# #' \dontshow{
# #' # for CRAN check only
# #' Ens_list <- getEnsList(Sim, temperatures, MaxIt = 5, m = 5)
# #' }
# #' eigenvalue_list <- getEigenvalueList(Ens_list)




getEigenvalueList <- function(EnsList) {
  eigenvalue_list <- lapply(EnsList, getEigenvalues)
  return(eigenvalue_list)
}




getEigenvalues = function(Ens){
  d = rowSums(Ens)
  n = nrow(Ens)
  Tmp = diag(d^(-1/2))
  NormalizeEns = Tmp %*% Ens %*% Tmp
  Eigenvalues = eigen(NormalizeEns)$values
  return(Eigenvalues)
}
