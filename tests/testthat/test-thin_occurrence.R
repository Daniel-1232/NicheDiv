test_that("thin.occurrence returns a data.frame with required columns", {
  x <- data.frame(
    Longitude = c(-100, -99.9, -99.8, -99.7),
    Latitude = c(40, 40.1, 40.2, 40.3),
    bio1 = c(1, 2, 3, 4)
  )
  out <- thin.occurrence(
    occurrence.data = x,
    longitude.col = "Longitude",
    latitude.col = "Latitude",
    thinning.dist.km = 1,
    N.thinning.replicates = 2,
    calc.Morans.I = FALSE,
    plot.Morans.I = FALSE,
    verbose = FALSE,
    seed = 1
  )
  expect_true(is.data.frame(out))
  expect_true(all(c("Longitude", "Latitude", "bio1") %in% names(out)))
  expect_lte(nrow(out), nrow(x))
})
