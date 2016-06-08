context("sampling Temperatures")

test_that("edgelists of more than 3 columns are not allowed", {
  # select temperatures
  temperatures <- temperatureSample(start = 0.01, end = 20, n = 20, method = 'random')

  # expectation
  expect_false(is.unsorted(temperatures))
})


