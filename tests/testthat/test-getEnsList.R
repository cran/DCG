context("get Ensemble matrices for each temperature")

SymAdjacencyMatrix <- as.symmetricAdjacencyMatrix(Data = monkeyGrooming, weighted = TRUE, rule = "weak")  # as.simMat checked.
Sim <- as.SimilarityMatrix(SymAdjacencyMatrix)  # as.simMat checked.
temperatures <- temperatureSample(start = 0.01, end = 20, n = 20, method = 'random')
Ens_list <- getEnsList(Sim, temperatures, MaxIt = 5, m = 5)

test_that("A warning raised if input is not symmetric", {
  assymetricMatrix <- matrix(1:9, 3, 3)
  expect_warning(getEnsList(assymetricMatrix, temperatures, MaxIt = 5, m = 5))
})



test_that("output is a list", {
  expect_is(Ens_list, "list")
})

test_that("the number of elements in the output list is the same as the length of temperatures", {
  expect_equal(length(Ens_list), length(temperatures))
})


test_that("the size of each matrix in the output list is the same as the size of the similarity matrix", {
  expect_equal(dim(Ens_list[[1]])[1], dim(Sim)[1])
  expect_equal(dim(Ens_list[[1]])[2], dim(Sim)[2])
})

