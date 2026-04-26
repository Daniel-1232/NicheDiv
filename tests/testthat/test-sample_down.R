test_that("sample.down returns requested number of rows", {
  x <- data.frame(a = 1:10, b = 11:20)
  out <- sample.down(dataframe = x, N.rows = 5, prioritize.NA.poisson = FALSE, seed = 1)
  expect_equal(nrow(out), 5)
  expect_equal(ncol(out), 2)
})

test_that("sample.down returns all rows when N.rows exceeds input", {
  x <- data.frame(a = 1:5, b = 6:10)
  out <- sample.down(dataframe = x, N.rows = 10, prioritize.NA.poisson = FALSE, seed = 1)
  expect_equal(nrow(out), 5)
})
