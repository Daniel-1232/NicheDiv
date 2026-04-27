#' Install external dependencies for NicheDiv
#'
#' Installs optional packages used by selected NicheDiv workflows,
#' including packages not available on CRAN such as ClimateNAr.
#'
#' @param install_climatena Logical; whether to install ClimateNAr.
#' @param install_whitebox Logical; whether to install whitebox.
#' @param install_data_table Logical; whether to install data.table.
#' @param climatena_url Download URL for the ClimateNAr Windows binary zip.
#' @param timeout Download timeout in seconds.
#' @param force Logical; reinstall even if package is already installed.
#' @param verbose Logical; print progress messages.
#'
#' @return Invisibly returns `TRUE` when finished.
#' @export
install_nichediv_dependencies <- function(install_climatena = TRUE,
                                          install_whitebox = TRUE,
                                          install_data_table = TRUE,
                                          climatena_url = "https://zenodo.org/records/17401570/files/ClimateNAr.zip?download=1",
                                          timeout = 800,
                                          force = FALSE,
                                          verbose = TRUE) {

  if (!is.logical(install_climatena) || length(install_climatena) != 1L) stop("install_climatena must be TRUE or FALSE")
  if (!is.logical(install_whitebox) || length(install_whitebox) != 1L) stop("install_whitebox must be TRUE or FALSE")
  if (!is.logical(install_data_table) || length(install_data_table) != 1L) stop("install_data_table must be TRUE or FALSE")
  if (!is.character(climatena_url) || length(climatena_url) != 1L) stop("climatena_url must be a single character string")
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
    sysname <- Sys.info()[["sysname"]]
    if (!grepl("Windows", sysname, ignore.case = TRUE)) stop("ClimateNAr installation is currently supported only on Windows")
    if (force || !"ClimateNAr" %in% rownames(utils::installed.packages())) {
      if (verbose) message("Installing ClimateNAr from Zenodo")
      dest <- file.path(tempdir(), "ClimateNAr.zip")
      old_timeout <- getOption("timeout")
      on.exit(options(timeout = old_timeout), add = TRUE)
      options(timeout = max(timeout, old_timeout))
      success <- FALSE
      for (i in 1:3) {
        try({
          utils::download.file(climatena_url, dest, mode = "wb", quiet = !verbose)
          if (file.exists(dest) && file.info(dest)$size > 1e7) {
            success <- TRUE
            break
          }
        }, silent = TRUE)
        Sys.sleep(3)
      }
      if (!success) stop("Download of ClimateNAr failed from: ", climatena_url)
      utils::install.packages(dest, repos = NULL, type = "win.binary")
      if (!"ClimateNAr" %in% rownames(utils::installed.packages())) stop("ClimateNAr installation finished, but the package could not be loaded.")
    } else {
      if (verbose) message("ClimateNAr already installed")
    }
  }
  invisible(TRUE)
}
