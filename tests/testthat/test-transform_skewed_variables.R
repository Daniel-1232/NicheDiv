test_that("transform.skewed.variables returns expected list structure", {
  x <- data.frame(a = c(1, 2, 3, 100), b = c(0.1, 0.2, 0.3, 0.9))
  out <- transform.skewed.variables(data.frame = x, verbose = FALSE)
  expect_true(is.list(out))
  expect_true(all(c("transformed", "summary", "background.transformed") %in% names(out)))
  expect_true(is.data.frame(out$transformed))
  expect_true(is.data.frame(out$summary))
})
