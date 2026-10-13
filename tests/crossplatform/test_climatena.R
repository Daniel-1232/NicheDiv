################################################################################
#### Settings
################################################################################
climatena_lib <- file.path(tempdir(), "ClimateNAr_crossplatform_library")
output_dir <- file.path(tempdir(), "ClimateNA_crossplatform_extraction")
expected_checksum <- "376896da27d584f9ad4091965022c46e"

################################################################################
#### Report test platform
################################################################################
print(Sys.info()[c("sysname", "release", "machine")])
print(R.version.string)

################################################################################
#### Install ClimateNAr 3.2.0 from Zenodo
################################################################################
source(file.path("R", "install_nichediv_dependencies.R"))

install_nichediv_dependencies(
  install_climatena = TRUE,
  install_whitebox = FALSE,
  install_data_table = TRUE,
  climatena_lib = climatena_lib,
  timeout = 1800,
  force = TRUE,
  verbose = TRUE
)

################################################################################
#### Verify downloaded archive and installed package
################################################################################
zip_file <- file.path(tempdir(), "ClimateNAr.zip")
if (!file.exists(zip_file)) stop("ClimateNAr ZIP was not downloaded")
checksum <- unname(tools::md5sum(zip_file))
print(checksum)
if (!identical(checksum, expected_checksum)) stop("ClimateNAr ZIP checksum differs from the verified Zenodo release")
if (as.character(utils::packageVersion("ClimateNAr")) != "3.2.0") stop("ClimateNAr version is not 3.2.0")
print(find.package("ClimateNAr"))

################################################################################
#### Load NicheDiv test branch
################################################################################
rm(install_nichediv_dependencies)
pkgload::load_all(".")

################################################################################
#### Create test occurrence data
################################################################################
test_occurrences <- data.frame(
  Longitude = c(-84.5, -84.0),
  Latitude = c(38.0, 38.5)
)
rownames(test_occurrences) <- c("occ_1", "occ_2")

################################################################################
#### Create cached elevation data
################################################################################
intermediate_dir <- file.path(output_dir, "Intermediate_files")
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

test_elevation <- data.frame(
  Elevation = c(300, 350),
  row.names = rownames(test_occurrences)
)
write.csv(test_elevation,
          file.path(intermediate_dir, "Elevation_extracted_occurrence.csv"))

################################################################################
#### Run NicheDiv ClimateNA extraction
################################################################################
oldwd <- getwd()
test_result <- tryCatch(
  extract.env.and.background(
    occurrence.data = test_occurrences,
    longitude.col = "Longitude",
    latitude.col = "Latitude",
    generate.background.data = FALSE,
    buffer.km = 20,
    csv.occurrence.out.file = "ClimateNA_occurrences_test.csv",
    output.dir = output_dir,
    env.datasets = c("elevation", "ClimateNA"),
    CRS.occurrences = "EPSG:4326",
    overwrite = FALSE,
    delete.intermediate.files.folders = FALSE,
    redownload.rasters = FALSE,
    seed = 1,
    verbose = TRUE
  ),
  finally = setwd(oldwd)
)

################################################################################
#### Check output and compare to Windows baseline
################################################################################
result_file <- file.path(output_dir, "ClimateNA_occurrences_test.csv")
if (!file.exists(result_file)) stop("NicheDiv did not generate the occurrence CSV")

climate_result <- read.csv(result_file, check.names = FALSE)
climate_variables <- setdiff(names(climate_result), c("ID", "Latitude", "Longitude", "Elevation"))

print(climate_result[, c("Longitude", "Latitude", "Elevation", "MAT", "MAP")])
print(dim(climate_result))
if (nrow(climate_result) != 2) stop("Incorrect number of occurrence records")
if (ncol(climate_result) != 209) stop("Incorrect number of output columns")
if (length(climate_variables) != 205) stop("Incorrect number of climate variables")
if (anyNA(climate_result[, climate_variables, drop = FALSE])) stop("ClimateNA output contains missing values")
if (anyDuplicated(names(climate_result))) stop("ClimateNA output contains duplicate variable names")
if (max(abs(climate_result$MAT - c(13.4, 12.8))) > 0.1) stop("MAT differs from Windows baseline")
if (max(abs(climate_result$MAP - c(1598, 1507))) > 2) stop("MAP differs from Windows baseline")

reference_file <- file.path("tests", "crossplatform", "reference_climatena_3.2.csv")
if (!file.exists(reference_file)) stop("Committed Windows reference CSV missing")
reference_result <- read.csv(reference_file, check.names = FALSE)
if (!identical(dim(reference_result), dim(climate_result))) stop("ClimateNA output dimensions differ from Windows reference")
if (!identical(names(reference_result), names(climate_result))) stop("ClimateNA variable names or order differ from Windows reference")
matching_columns <- vapply(names(reference_result), function(x) {
  isTRUE(all.equal(reference_result[[x]], climate_result[[x]], tolerance = 1e-6, check.attributes = FALSE))
}, logical(1))
message("Matching columns: ", sum(matching_columns), " / ", length(matching_columns))
if (!all(matching_columns)) {
  print(names(reference_result)[!matching_columns])
  stop("ClimateNA values differ from the Windows reference")
}

message("PASS: ClimateNAr 3.2.0 and all 205 NicheDiv climate variables match Windows on ", Sys.info()[["sysname"]])

