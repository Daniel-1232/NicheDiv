################################################################################
#### Test custom accessible-area polygons
################################################################################
## Create test occurrence and background data
test_occurrences <- data.frame(
  Longitude = c(-500, 500, -500, 500),
  Latitude = c(-500, -500, 500, 500)
)

test_background <- data.frame(
  Longitude = c(0, 500, 1000, 2500, -500, 0),
  Latitude = c(0, 500, 0, 2500, -500, 700),
  Record_ID = c("A", "B", "C", "D", "E", "F")
)

## Create function to construct rectangular polygon coordinates
create.square <- function(xmin, xmax, ymin, ymax) {
  rbind(c(xmin, ymin), c(xmax, ymin), c(xmax, ymax),
        c(xmin, ymax), c(xmin, ymin))
}

## Create test polygon
test_polygon <- sf::st_sf(
  geometry = sf::st_sfc(
    sf::st_polygon(list(create.square(-1000, 1000, -1000, 1000))),
    crs = 3857
  )
)

## Create function to filter background using custom polygons
filter.custom.background <- function(custom_polygon, ...) {
  crop.background.buffered(
    occurrence.data = test_occurrences,
    background.data = test_background,
    CRS = "EPSG:3857",
    buffer.method = "custom",
    custom.polygon = custom_polygon,
    verbose = FALSE,
    ...
  )
}



################################################################################
#### Test custom polygon input formats
################################################################################
testthat::test_that("sf and sfc polygon inputs work", {
  sf_result <- filter.custom.background(test_polygon)
  sfc_result <- filter.custom.background(sf::st_geometry(test_polygon))

  testthat::expect_identical(sf_result$Record_ID, c("A", "B", "C", "E", "F"))
  testthat::expect_identical(sfc_result, sf_result)
})

testthat::test_that("terra SpatVector polygon input works", {
  terra_polygon <- terra::vect(test_polygon)
  terra_result <- filter.custom.background(terra_polygon)

  testthat::expect_identical(terra_result$Record_ID, c("A", "B", "C", "E", "F"))
})

testthat::test_that("GeoPackage input works", {
  geopackage_file <- tempfile(fileext = ".gpkg")
  on.exit(unlink(geopackage_file), add = TRUE)
  sf::st_write(test_polygon, geopackage_file, quiet = TRUE)

  geopackage_result <- filter.custom.background(geopackage_file)
  testthat::expect_identical(geopackage_result$Record_ID, c("A", "B", "C", "E", "F"))
})

testthat::test_that("GeoJSON input works", {
  geojson_file <- tempfile(fileext = ".geojson")
  on.exit(unlink(geojson_file), add = TRUE)
  sf::st_write(test_polygon, geojson_file, quiet = TRUE)

  geojson_result <- filter.custom.background(geojson_file)
  testthat::expect_true(all(c("A", "B", "E", "F") %in% geojson_result$Record_ID))
  testthat::expect_false("D" %in% geojson_result$Record_ID)
})

testthat::test_that("ESRI Shapefile input works", {
  shapefile_directory <- tempfile()
  dir.create(shapefile_directory)
  on.exit(unlink(shapefile_directory, recursive = TRUE), add = TRUE)
  shapefile_path <- file.path(shapefile_directory, "accessible_area.shp")
  sf::st_write(test_polygon, shapefile_path, quiet = TRUE)

  shapefile_result <- filter.custom.background(shapefile_path)
  testthat::expect_identical(shapefile_result$Record_ID, c("A", "B", "C", "E", "F"))
})



################################################################################
#### Test coordinate reference systems
################################################################################
testthat::test_that("Custom polygons are automatically reprojected", {
  geographic_polygon <- sf::st_transform(test_polygon, 4326)
  geographic_result <- filter.custom.background(geographic_polygon)

  testthat::expect_true(all(c("A", "B", "E", "F") %in% geographic_result$Record_ID))
  testthat::expect_false("D" %in% geographic_result$Record_ID)
})

testthat::test_that("Geographic occurrence coordinates work", {
  geographic_occurrences <- data.frame(
    Longitude = c(-84.5, -84.4, -84.6),
    Latitude = c(38, 38.1, 38.05)
  )
  geographic_background <- data.frame(
    Longitude = c(-84.5, -84.45, -80),
    Latitude = c(38, 38.05, 40),
    Record_ID = c("A", "B", "C")
  )
  geographic_polygon <- sf::st_sf(
    geometry = sf::st_sfc(
      sf::st_polygon(list(create.square(-84.7, -84.3, 37.9, 38.2))),
      crs = 4326
    )
  )

  geographic_result <- crop.background.buffered(
    occurrence.data = geographic_occurrences,
    background.data = geographic_background,
    buffer.method = "custom",
    custom.polygon = geographic_polygon,
    verbose = FALSE
  )

  testthat::expect_identical(geographic_result$Record_ID, c("A", "B"))
})



################################################################################
#### Test polygon topology and spatial filtering
################################################################################
testthat::test_that("Points on polygon boundaries are retained", {
  filtered_background <- filter.custom.background(test_polygon)

  testthat::expect_true("C" %in% filtered_background$Record_ID)
  testthat::expect_false("D" %in% filtered_background$Record_ID)
})

testthat::test_that("Polygon holes exclude background points", {
  polygon_with_hole <- sf::st_sf(
    geometry = sf::st_sfc(
      sf::st_polygon(list(
        create.square(-1000, 1000, -1000, 1000),
        rbind(c(-100, -100), c(-100, 100), c(100, 100),
              c(100, -100), c(-100, -100))
      )),
      crs = 3857
    )
  )

  hole_result <- filter.custom.background(polygon_with_hole)
  testthat::expect_identical(hole_result$Record_ID, c("B", "C", "E", "F"))
})

testthat::test_that("Disconnected polygon features are combined correctly", {
  disconnected_polygons <- sf::st_sf(
    geometry = sf::st_sfc(
      sf::st_polygon(list(create.square(-700, -300, -700, -300))),
      sf::st_polygon(list(create.square(300, 700, 300, 700))),
      crs = 3857
    )
  )

  disconnected_result <- filter.custom.background(disconnected_polygons)
  testthat::expect_identical(disconnected_result$Record_ID, c("B", "E"))
})

testthat::test_that("Custom polygons are not buffered", {
  unbuffered_result <- filter.custom.background(test_polygon)
  buffered_result <- filter.custom.background(test_polygon, buffer.dist.meters = 2000)

  testthat::expect_identical(buffered_result, unbuffered_result)
})

testthat::test_that("Background columns and row order are preserved", {
  filtered_background <- filter.custom.background(test_polygon)
  expected_background <- test_background[c(1, 2, 3, 5, 6), , drop = FALSE]

  testthat::expect_identical(filtered_background, expected_background)
})

testthat::test_that("Missing background coordinates are removed", {
  background_with_missing <- rbind(
    test_background,
    data.frame(Longitude = NA_real_, Latitude = NA_real_, Record_ID = "G")
  )

  filtered_background <- crop.background.buffered(
    occurrence.data = test_occurrences,
    background.data = background_with_missing,
    CRS = "EPSG:3857",
    buffer.method = "custom",
    custom.polygon = test_polygon,
    verbose = FALSE
  )

  testthat::expect_identical(filtered_background$Record_ID, c("A", "B", "C", "E", "F"))
})



################################################################################
#### Test invalid custom polygon inputs
################################################################################
testthat::test_that("Missing polygons produce an error", {
  testthat::expect_error(
    crop.background.buffered(
      occurrence.data = test_occurrences,
      background.data = test_background,
      CRS = "EPSG:3857",
      buffer.method = "custom",
      verbose = FALSE
    ),
    "custom.polygon must be supplied"
  )
})

testthat::test_that("Undefined polygon CRS produces an error", {
  polygon_without_crs <- sf::st_set_crs(test_polygon, NA)
  testthat::expect_error(filter.custom.background(polygon_without_crs), "defined CRS")
})

testthat::test_that("Non-polygon geometries are rejected", {
  point_geometry <- sf::st_sf(
    geometry = sf::st_sfc(sf::st_point(c(0, 0)), crs = 3857)
  )

  testthat::expect_error(filter.custom.background(point_geometry), "polygon")
})

testthat::test_that("Missing spatial files are rejected", {
  missing_file <- tempfile(fileext = ".gpkg")
  testthat::expect_error(filter.custom.background(missing_file), "does not exist")
})

testthat::test_that("Custom polygons cannot be used with other methods", {
  testthat::expect_error(
    crop.background.buffered(
      occurrence.data = test_occurrences,
      background.data = test_background,
      CRS = "EPSG:3857",
      buffer.dist.meters = 100,
      buffer.method = "hull",
      custom.polygon = test_polygon,
      verbose = FALSE
    ),
    "custom.polygon requires"
  )
})



################################################################################
#### Test existing accessible-area methods
################################################################################
testthat::test_that("Existing hull method remains functional", {
  hull_result <- crop.background.buffered(
    test_occurrences, test_background, CRS = "EPSG:3857",
    buffer.dist.meters = 500, buffer.method = "hull", verbose = FALSE
  )

  testthat::expect_true("E" %in% hull_result$Record_ID)
  testthat::expect_false("D" %in% hull_result$Record_ID)
})

testthat::test_that("Existing points method remains functional", {
  points_result <- crop.background.buffered(
    test_occurrences, test_background, CRS = "EPSG:3857",
    buffer.dist.meters = 50, buffer.method = "points", verbose = FALSE
  )

  testthat::expect_true(all(c("B", "E") %in% points_result$Record_ID))
  testthat::expect_false("D" %in% points_result$Record_ID)
})

testthat::test_that("Existing bbox method remains functional", {
  bbox_result <- crop.background.buffered(
    test_occurrences, test_background, CRS = "EPSG:3857",
    buffer.dist.meters = 50, buffer.method = "bbox", verbose = FALSE
  )

  testthat::expect_true(all(c("A", "B", "E") %in% bbox_result$Record_ID))
  testthat::expect_false("D" %in% bbox_result$Record_ID)
})

testthat::test_that("Existing positional arguments remain compatible", {
  positional_result <- crop.background.buffered(
    test_occurrences, test_background, "Latitude", "Longitude",
    "EPSG:3857", 50, "points", 3, FALSE
  )

  named_result <- crop.background.buffered(
    occurrence.data = test_occurrences,
    background.data = test_background,
    CRS = "EPSG:3857",
    buffer.dist.meters = 50,
    buffer.method = "points",
    alpha = 3,
    verbose = FALSE
  )

  testthat::expect_identical(positional_result, named_result)
})
