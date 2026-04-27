test_that("extract.env.and.background validates missing env.dataset", {
  occ <- data.frame(Longitude = c(-100, -99.9), Latitude = c(40, 40.1))
  expect_error(
    extract.env.and.background(
      occurrence.data = occ,
      latitude.col = "Latitude",
      longitude.col = "Longitude",
      env.dataset = "not_a_real_dataset",
      verbose = FALSE
    )
  )
})
