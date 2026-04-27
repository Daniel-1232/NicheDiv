test_that("trim.to.analogous.environments errors on non-data.frame input", {
  expect_error(
    trim.to.analogous.environments(
      Sp1.occurrence.data = 1:10,
      Sp2.occurrence.data = data.frame(id = 1:10, var1 = rnorm(10)),
      Sp1.background.data = data.frame(id = 1:10, var1 = rnorm(10)),
      Sp2.background.data = data.frame(id = 1:10, var1 = rnorm(10)),
      exclude.cols = "id",
      verbose = FALSE
    )
  )
})
