test_that("filter.analogous.variables returns filtered occurrence data.frame", {
  set.seed(1)
  Sp1.background.data <- data.frame(
    id = 1:50,
    var1 = rnorm(50, 0, 1),
    var2 = rnorm(50, 5, 1),
    var3 = rep(1, 50)
  )
  Sp2.background.data <- data.frame(
    id = 1:50,
    var1 = rnorm(50, 0.1, 1),
    var2 = rnorm(50, 5.1, 1),
    var3 = rep(1, 50)
  )
  Sp1.Sp2.occurrence.data <- data.frame(
    id = 1:10,
    var1 = rnorm(10, 0, 1),
    var2 = rnorm(10, 5, 1),
    var3 = rep(1, 10)
  )
  out <- filter.analogous.variables(
    Sp1.background.data = Sp1.background.data,
    Sp2.background.data = Sp2.background.data,
    Sp1.Sp2.occurrence.data = Sp1.Sp2.occurrence.data,
    exclude.cols = "id",
    CV.threshold = 0.01,
    overlap.threshold = 0,
    plot.1D.overlap = FALSE,
    use.parallel = FALSE,
    verbose = FALSE
  )
  expect_true(is.data.frame(out))
  expect_true("id" %in% names(out))
  expect_false("var3" %in% names(out))
})
