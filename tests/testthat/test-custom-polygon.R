testthat::test_that("Custom accessible-area polygons filter background correctly", {
  occurrence_data <- data.frame(Longitude = c(0, 0.5, -0.5), Latitude = c(0, 0.5, -0.5))
  background_data <- data.frame(Longitude = c(0, 0.5, 5, -0.5), Latitude = c(0, 0.5, 5, -0.5), Record_ID = c("A", "B", "C", "D"))

  polygon_coordinates <- rbind(c(-1, -1), c(1, -1), c(1, 1), c(-1, 1), c(-1, -1))
  custom_polygon <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(polygon_coordinates)), crs = 4326))

  filtered_data <- crop.background.buffered(occurrence_data, background_data, buffer.method = "custom", custom.polygon = custom_polygon, verbose = FALSE)
  testthat::expect_identical(filtered_data$Record_ID, c("A", "B", "D"))

  projected_polygon <- sf::st_transform(custom_polygon, 3857)
  projected_data <- crop.background.buffered(occurrence_data, background_data, buffer.method = "custom", custom.polygon = projected_polygon, verbose = FALSE)
  testthat::expect_identical(projected_data$Record_ID, c("A", "B", "D"))

  terra_polygon <- terra::vect(custom_polygon)
  terra_data <- crop.background.buffered(occurrence_data, background_data, buffer.method = "custom", custom.polygon = terra_polygon, verbose = FALSE)
  testthat::expect_identical(terra_data$Record_ID, c("A", "B", "D"))

  geopackage_file <- tempfile(fileext = ".gpkg")
  sf::st_write(custom_polygon, geopackage_file, quiet = TRUE)
  geopackage_data <- crop.background.buffered(occurrence_data, background_data, buffer.method = "custom", custom.polygon = geopackage_file, verbose = FALSE)
  testthat::expect_identical(geopackage_data$Record_ID, c("A", "B", "D"))
  unlink(geopackage_file)

  exterior_coordinates <- rbind(c(-2, -2), c(2, -2), c(2, 2), c(-2, 2), c(-2, -2))
  interior_coordinates <- rbind(c(-0.2, -0.2), c(-0.2, 0.2), c(0.2, 0.2), c(0.2, -0.2), c(-0.2, -0.2))
  polygon_with_hole <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(exterior_coordinates, interior_coordinates)), crs = 4326))
  hole_data <- crop.background.buffered(occurrence_data, background_data, buffer.method = "custom", custom.polygon = polygon_with_hole, verbose = FALSE)
  testthat::expect_identical(hole_data$Record_ID, c("B", "D"))

  testthat::expect_error(crop.background.buffered(occurrence_data, background_data, buffer.method = "custom", verbose = FALSE), "custom.polygon must be supplied")
  testthat::expect_error(crop.background.buffered(occurrence_data, background_data, buffer.method = "custom", custom.polygon = sf::st_set_crs(custom_polygon, NA), verbose = FALSE), "defined CRS")
})
