context("Importing Data Tests")


test_that("edgelists of more than 3 columns are not allowed", {
  testEdgelist1 <- data.frame(col1 = letters[1:10], col2 = letters[1:10],
                              col3 = sample(1:10, 10, replace = TRUE),
                              col4 = rnorm(1:10))
  expect_error(as.symmetricAdjacencyMatrix(testEdgelist1))
})

test_that("factors are not allowed in edgelist", {
  testEdgelist2 <- data.frame(col1 = letters[1:10], col2 = letters[11:20],
                              stringsAsFactors = TRUE)
  expect_error(as.symmetricAdjacencyMatrix(testEdgelist2))
})



test_that("only a three-column edgelist is allowed if 'weighted = TRUE'", {
  testEdgelist3 <- data.frame(col1 = letters[1:10],
                              col2 = letters[10:1],
                              stringsAsFactors = FALSE)
  expect_error(as.symmetricAdjacencyMatrix(testEdgelist3, weighted = TRUE))
})

test_that("A warning raised if dyads in weighted edgelist are not unique", {
  testEdgelist4 <- data.frame(col1 = letters[c(1:10, 1)], col2 = letters[c(10:1, 10)],
                              col3 = sample(1:10, 11, replace = TRUE),
                              stringsAsFactors = FALSE)
  expect_warning(as.symmetricAdjacencyMatrix(testEdgelist4, weighted = TRUE))
})


test_that("the initiator and the recipient should not be the same", {
  testEdgelist5 <- data.frame(col1 = letters[1:10],
                              col2 = letters[10:1],
                              stringsAsFactors = FALSE)
  expect_error(as.symmetricAdjacencyMatrix(testEdgelist5, weighted = TRUE))
})


test_that("returns similarity matrix", {
  testEdgelist6 <- data.frame(col1 = letters[1:10],
                              col2 = letters[10:1],
                              stringsAsFactors = FALSE)
  symmetricAdjacencyMatrix <- as.symmetricAdjacencyMatrix(testEdgelist6)
  expect_is(as.SimilarityMatrix(symmetricAdjacencyMatrix), "similarityMatrix")
})


test_that("adjacency matrix is symmetric", {
  testEdgelist6 <- data.frame(col1 = letters[1:10],
                              col2 = letters[10:1],
                              stringsAsFactors = FALSE)
  symmetricAdjacencyMatrix <- as.symmetricAdjacencyMatrix(testEdgelist6)
  expect_true(isSymmetric(symmetricAdjacencyMatrix))
})


test_that("diagonal of the raw win-loss matrix should be zeros", {
  set.seed(1)
  testMatrix7 <- matrix(sample(1:100, 100, TRUE), 10, 10)
  diag(testMatrix7) <- sample(c(0, 1), 10, TRUE, prob = c(0.9, 0.1))
  expect_warning(as.symmetricAdjacencyMatrix(testMatrix7))
})








