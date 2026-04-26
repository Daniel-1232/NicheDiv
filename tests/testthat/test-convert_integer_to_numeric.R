test_that("convert.integer.to.numeric converts integer columns only", {
  x <- data.frame(a = 1:3, b = c(1.1, 2.2, 3.3), c = c("x", "y", "z"))
  out <- convert.integer.to.numeric(x)
  expect_true(is.numeric(out$a))
  expect_true(is.numeric(out$b))
  expect_true(is.character(out$c))
  expect_identical(dim(out), dim(x))
})

test_that("convert.integer.to.numeric errors on non-data.frame input", {
  expect_error(convert.integer.to.numeric(1:3), "dataframe must be a data.frame")
})
