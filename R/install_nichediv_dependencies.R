#' Install external dependencies for NicheDiv
#'
#' Installs optional packages used by selected NicheDiv workflows,
#' including packages not available on CRAN such as ClimateNAr.
#'
#' @param install_climatena Logical; whether to install ClimateNAr.
#' @param install_whitebox Logical; whether to install whitebox.
#' @param install_data_table Logical; whether to install data.table.
#' @param climatena_url Download URL for a ClimateNAr ZIP archive.
#' @param climatena_zip Optional path to a local ClimateNAr ZIP (overrides climatena_url).
#' @param climatena_lib Optional R library for an isolated ClimateNAr installation.
#' @param timeout Download timeout in seconds.
#' @param force Logical; reinstall even if package is already installed.
#' @param verbose Logical; print progress messages.
#'
#' @return Invisibly returns `TRUE` when finished.
#' @export
install_nichediv_dependencies <- function(install_climatena = TRUE,
                                          install_whitebox = TRUE,
                                          install_data_table = TRUE,
                                          climatena_url = "https://zenodo.org/records/23327283/files/ClimateNAr.zip?download=1",
                                          climatena_zip = NULL,
                                          climatena_lib = NULL,
                                          timeout = 800,
                                          force = FALSE,
                                          verbose = TRUE) {

  if (!is.logical(install_climatena) || length(install_climatena) != 1L) stop("install_climatena must be TRUE or FALSE")
  if (!is.logical(install_whitebox) || length(install_whitebox) != 1L) stop("install_whitebox must be TRUE or FALSE")
  if (!is.logical(install_data_table) || length(install_data_table) != 1L) stop("install_data_table must be TRUE or FALSE")
  if (!is.character(climatena_url) || length(climatena_url) != 1L || is.na(climatena_url) || !nzchar(climatena_url)) stop("climatena_url must be a single non-empty character string")
  if (!is.null(climatena_zip) && (!is.character(climatena_zip) || length(climatena_zip) != 1L || is.na(climatena_zip) || !nzchar(climatena_zip))) stop("climatena_zip must be NULL or a single file path")
  if (!is.null(climatena_lib) && (!is.character(climatena_lib) || length(climatena_lib) != 1L || is.na(climatena_lib) || !nzchar(climatena_lib))) stop("climatena_lib must be NULL or a single folder path")
  if (!is.numeric(timeout) || length(timeout) != 1L || !is.finite(timeout) || timeout <= 0) stop("timeout must be a single positive number")
  if (!is.logical(force) || length(force) != 1L) stop("force must be TRUE or FALSE")
  if (!is.logical(verbose) || length(verbose) != 1L) stop("verbose must be TRUE or FALSE")
  install_if_needed <- function(pkg) {
    if (force || !requireNamespace(pkg, quietly = TRUE)) {
      if (verbose) message("Installing ", pkg, " from CRAN")
      utils::install.packages(pkg)
    } else {
      if (verbose) message(pkg, " already installed")
    }
  }
  if (install_whitebox) install_if_needed("whitebox")
  if (install_data_table) install_if_needed("data.table")
  if (install_climatena) {
    if (!is.null(climatena_lib)) {
      dir.create(climatena_lib, recursive = TRUE, showWarnings = FALSE)
      if (!dir.exists(climatena_lib)) stop("Cannot create ClimateNAr library: ", climatena_lib)
      .libPaths(c(normalizePath(climatena_lib, winslash = "/"), .libPaths()))
    }
    target_lib <- .libPaths()[1L]
    if (file.access(target_lib, 2L) != 0L) stop("R library is not writable: ", target_lib)
    for (pkg in c("data.table", "terra", "readr")) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        if (verbose) message("Installing ClimateNAr dependency: ", pkg)
        utils::install.packages(pkg)
      }
    }
    installed_here <- file.exists(file.path(target_lib, "ClimateNAr", "DESCRIPTION"))
    existing_version <- if (requireNamespace("ClimateNAr", quietly = TRUE)) {
      utils::packageVersion("ClimateNAr")
    } else {
      NULL
    }
    needs_install <- force ||
      (!is.null(climatena_lib) && !installed_here) ||
      is.null(existing_version) ||
      existing_version < utils::package_version("3.2.0")
    if (needs_install) {
      if ("ClimateNAr" %in% loadedNamespaces()) stop("Restart R before reinstalling ClimateNAr: its namespace is already loaded")
      dest <- file.path(tempdir(), "ClimateNAr.zip")
      if (!is.null(climatena_zip)) {
        if (!file.exists(climatena_zip)) stop("ClimateNAr ZIP not found: ", climatena_zip)
        if (!identical(normalizePath(climatena_zip, winslash = "/"), normalizePath(dest, winslash = "/", mustWork = FALSE))) {
          if (!file.copy(climatena_zip, dest, overwrite = TRUE)) stop("Failed to copy local ClimateNAr ZIP")
        }
        if (verbose) message("Installing ClimateNAr from local ZIP: ", climatena_zip)
      } else {
        if (verbose) message("Downloading ClimateNAr from: ", climatena_url)
        old_timeout <- getOption("timeout")
        on.exit(options(timeout = old_timeout), add = TRUE)
        options(timeout = max(timeout, old_timeout))
        success <- FALSE
        for (i in 1:3) {
          if (file.exists(dest)) unlink(dest)
          result <- try(utils::download.file(climatena_url, dest, mode = "wb", quiet = !verbose), silent = TRUE)
          if (!inherits(result, "try-error") && is.numeric(result) && length(result) == 1L && result == 0 && file.exists(dest) && file.info(dest)$size > 1e7) {
            success <- TRUE
            break
          }
          if (i < 3) Sys.sleep(3)
        }
        if (!success) stop("Download of ClimateNAr failed from: ", climatena_url)
      }
      if (!file.exists(dest) || !is.finite(file.info(dest)$size) || file.info(dest)$size < 1e7) stop("ClimateNAr archive is missing or unexpectedly small")
      zip_entries <- tryCatch(utils::unzip(dest, list = TRUE)$Name,
                              error = function(e) stop("Invalid ClimateNAr ZIP: ", conditionMessage(e)))
      if (!("ClimateNAr/DESCRIPTION" %in% zip_entries)) stop("ClimateNAr ZIP does not contain ClimateNAr/DESCRIPTION")
      if (any(!startsWith(zip_entries, "ClimateNAr/")) || any(grepl("(^|/)\\.\\.(/|$)", zip_entries))) stop("ClimateNAr ZIP contains unexpected file paths")
      if (.Platform$OS.type == "windows") {
        utils::install.packages(dest, repos = NULL, type = "win.binary", lib = target_lib)
      } else {
        if (any(startsWith(zip_entries, "ClimateNAr/libs/"))) stop("The ZIP contains compiled libraries; a native source build is required")
        extract_dir <- tempfile("ClimateNAr_extract_")
        dir.create(extract_dir)
        on.exit(unlink(extract_dir, recursive = TRUE, force = TRUE), add = TRUE)
        utils::unzip(dest, exdir = extract_dir)
        source_dir <- file.path(extract_dir, "ClimateNAr")
        if (!file.exists(file.path(source_dir, "DESCRIPTION"))) stop("ClimateNAr extraction failed")
        archive_description <- read.dcf(file.path(source_dir, "DESCRIPTION"))
        if (!identical(unname(archive_description[1, "Package"]), "ClimateNAr")) stop("Unexpected package in ClimateNAr archive")
        target_dir <- file.path(target_lib, "ClimateNAr")
        if (dir.exists(target_dir)) unlink(target_dir, recursive = TRUE, force = TRUE)
        if (!isTRUE(file.copy(source_dir, target_lib, recursive = TRUE)) || !file.exists(file.path(target_dir, "DESCRIPTION"))) stop("Failed to install ClimateNAr into: ", target_lib)
      }
      if (!requireNamespace("ClimateNAr", quietly = TRUE)) stop("ClimateNAr was installed, but its namespace could not be loaded; the archive may not be portable to this system")
    } else {
      if (verbose) message("ClimateNAr already installed")
    }
    if (verbose) message("ClimateNAr version: ", as.character(utils::packageVersion("ClimateNAr")))
  }
  invisible(TRUE)
}
