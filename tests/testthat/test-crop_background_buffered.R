test_that("crop.background.buffered returns filtered background data.frame", {
  occurrence.data <- data.frame(
    Longitude = c(-100, -99.9, -99.8),
    Latitude = c(40, 40.1, 40.2)
  )
  background.data <- data.frame(
    Longitude = c(-100, -99.95, -99.7, -120),
    Latitude = c(40, 40.05, 40.3, 60),
    var1 = c(1, 2, 3, 4)
  )
  out <- crop.background.buffered(
    occurrence.data = occurrence.data,
    background.data = background.data,
    latitude.col = "Latitude",
    longitude.col = "Longitude",
    CRS = "EPSG:4326",
    buffer.dist.meters = 50000,
    buffer.method = "bbox",
    verbose = FALSE
  )
  expect_true(is.data.frame(out))
  expect_true(nrow(out) >= 1)
  expect_true(all(names(background.data) %in% names(out)))
})
