################################################################################
#### Settings
################################################################################
output_dir <- file.path(tempdir(), "NicheDiv_background_test")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

################################################################################
#### Load NicheDiv
################################################################################
pkgload::load_all(".", quiet = TRUE)

################################################################################
#### Prepare small synthetic landmask to avoid GADM download
################################################################################
landmask_dir <- file.path(output_dir, "Env_rasters", "landmask")
dir.create(landmask_dir, recursive = TRUE, showWarnings = FALSE)
ring <- matrix(c(-86, 36, -82, 36, -82, 40, -86, 40, -86, 36), ncol = 2, byrow = TRUE)
land_sf <- sf::st_sf(GID_0 = "USA", geometry = sf::st_sfc(sf::st_polygon(list(ring)), crs = 4326))
land_vect <- terra::vect(land_sf)
saveRDS(terra::wrap(land_vect), file.path(landmask_dir, "gadm36_adm0_r5_pk.rds"))

################################################################################
#### Prepare local synthetic raster
################################################################################
raster <- terra::rast(ncols = 100, nrows = 100, xmin = -86, xmax = -82, ymin = 36, ymax = 40, crs = "EPSG:4326")
terra::values(raster) <- 20
raster_file <- file.path(tempdir(), "NicheDiv_synthetic_climate.tif")
terra::writeRaster(raster, raster_file, overwrite = TRUE)

################################################################################
#### Generate synthetic occurrences and real background samples
################################################################################
occ <- data.frame(Longitude = c(-84.6, -84.0, -84.5, -84.1), Latitude = c(37.9, 38.0, 38.4, 38.6))
rownames(occ) <- paste0("occ_", seq_len(nrow(occ)))
extract.env.and.background(
  occurrence.data = occ,
  longitude.col = "Longitude",
  latitude.col = "Latitude",
  generate.background.data = TRUE,
  N.background.points = 12,
  remove.hydrolakes.background = FALSE,
  buffer.km = 20,
  csv.occurrence.out.file = "synthetic_occurrences.csv",
  csv.background.out.file = "synthetic_background.csv",
  output.dir = output_dir,
  env.datasets = NULL,
  custom.env.rasters = raster_file,
  custom.env.rasters.names = "synthetic",
  custom.env.rasters.variable.names = "synthetic_temperature",
  overwrite = FALSE,
  delete.intermediate.files.folders = FALSE,
  redownload.rasters = FALSE,
  seed = 137,
  verbose = TRUE
)

################################################################################
#### Validate background and occurrence outputs
################################################################################
occ_file <- file.path(output_dir, "synthetic_occurrences.csv")
bg_file <- file.path(output_dir, "synthetic_background.csv")
if (!file.exists(occ_file) || !file.exists(bg_file)) stop("Synthetic background/occurrence CSV was not written")
occ_result <- read.csv(occ_file, check.names = FALSE)
bg_result <- read.csv(bg_file, check.names = FALSE)
if (nrow(occ_result) != 4L) stop("Incorrect occurrence record count")
if (nrow(bg_result) != 12L) stop("Incorrect background record count")
if (!all(c("Latitude", "Longitude") %in% names(bg_result))) stop("Background coordinates are missing")
if (any(!is.finite(bg_result$Latitude)) || any(!is.finite(bg_result$Longitude))) stop("Background coordinates contain invalid values")
synthetic_columns <- grep("synthetic_temperature", names(bg_result), value = TRUE)
if (length(synthetic_columns) != 1L || anyNA(bg_result[[synthetic_columns]])) stop("Background synthetic raster extraction failed")
if (max(abs(bg_result[[synthetic_columns]] - 20)) > 0.001) stop("Incorrect synthetic background raster values")
message("PASS: Background sampling and synthetic raster extraction on ", Sys.info()[["sysname"]])
