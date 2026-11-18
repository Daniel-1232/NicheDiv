#### Cross-platform functional smoke tests ####################################

## Exercise each exported function on small, offline, synthetic inputs.
## A passing smoke test is not a proof of all argument combinations or large-data behavior.

library(NicheDiv)

## Set reproducible environment
set.seed(20261118)
options(warn = 1)

## Track individual probes even when one probe fails
results <- data.frame(
  Function = character(),
  Status = character(),
  Seconds = numeric(),
  Detail = character(),
  stringsAsFactors = FALSE
)

run_probe <- function(name, expr) {
  cat("\n=== ", name, " ===\n", sep = "")
  started <- proc.time()[["elapsed"]]
  error_message <- NULL

  tryCatch(
    eval.parent(substitute(expr)),
    error = function(e) error_message <<- conditionMessage(e)
  )

  elapsed <- proc.time()[["elapsed"]] - started
  passed <- is.null(error_message)
  status <- if (passed) "PASS" else "FAIL"
  detail <- if (passed) "Functional probe succeeded" else error_message

  results[nrow(results) + 1L, ] <<- list(name, status, round(elapsed, 2), detail)
  cat(status, ": ", detail, "\n", sep = "")
  invisible(passed)
}

assert <- function(condition, message = "Unexpected test result") {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

## Ensure no newly exported function silently escapes the test matrix
expected_exports <- c(
  "calc.niche.divergence.metrics",
  "convert.integer.to.numeric",
  "crop.background.buffered",
  "extract.env.and.background",
  "filter.analogous.variables",
  "install_nichediv_dependencies",
  "map.env.variable.names",
  "plot.DAPC.niche.divergence",
  "plot.DAPC.permutation",
  "plot.DAPC.var.contributions",
  "plot.occurrences.map",
  "plot.top.DAPC.predictors",
  "remove.low.CV.vars",
  "run.DAPC.crossval.permutation",
  "sample.down",
  "thin.occurrence",
  "transform.skewed.variables",
  "trim.to.analogous.environments"
)
actual_exports <- getNamespaceExports("NicheDiv")
missing <- setdiff(expected_exports, actual_exports)
new_exports <- setdiff(actual_exports, expected_exports)
if (length(missing) || length(new_exports)) {
  stop("Export list changed; update this test before merging. Missing: ",
       paste(missing, collapse = ", "), "; new: ", paste(new_exports, collapse = ", "))
}
assert(all(vapply(actual_exports, function(x) is.function(getExportedValue("NicheDiv", x)), logical(1))),
       "A declared export is not a callable function")

## Common synthetic occurrence and background environment
n <- 50L
species <- factor(rep(c("Sp1", "Sp2"), each = n), levels = c("Sp1", "Sp2"))
coords <- data.frame(
  Longitude = c(runif(n, -84.8, -84.6), runif(n, -84.7, -84.5)),
  Latitude = c(runif(n, 38.0, 38.2), runif(n, 38.1, 38.3))
)
occ <- data.frame(
  Species = species,
  Env1 = c(rnorm(n, 0, 1), rnorm(n, 0.4, 1)),
  Env2 = c(rnorm(n, 1, 1.1), rnorm(n, 1.3, 1.1)),
  Env3 = c(rnorm(n, -0.5, 1), rnorm(n, -0.3, 1)),
  Env4 = c(rnorm(n, 2, 0.8), rnorm(n, 2.2, 0.8))
)
occ_coord <- cbind(coords, occ)
occ1 <- occ_coord[species == "Sp1", , drop = FALSE]
occ2 <- occ_coord[species == "Sp2", , drop = FALSE]
cols_env <- paste0("Env", 1:4)

## Backgrounds sampled near occurrence values, with overlapping support
bg1 <- occ1[rep(seq_len(nrow(occ1)), each = 3), c("Longitude", "Latitude", cols_env)]
bg2 <- occ2[rep(seq_len(nrow(occ2)), each = 3), c("Longitude", "Latitude", cols_env)]
for (v in cols_env) {
  bg1[[v]] <- bg1[[v]] + rnorm(nrow(bg1), sd = 0.12)
  bg2[[v]] <- bg2[[v]] + rnorm(nrow(bg2), sd = 0.12)
}
rownames(bg1) <- NULL
rownames(bg2) <- NULL

## Helper to check ggplot objects
assert_ggplot <- function(x) {
  assert(inherits(x, "ggplot"), "Expected ggplot result")
  ggplot2::ggplot_build(x)
  invisible(TRUE)
}

## 1. Independent preprocessing functions
run_probe("convert.integer.to.numeric", {
  x <- NicheDiv::convert.integer.to.numeric(data.frame(a = 1:4, b = letters[1:4]))
  assert(is.double(x$a) && identical(x$b, letters[1:4]))
})

run_probe("sample.down", {
  x <- NicheDiv::sample.down(occ_coord, N.rows = 25, prioritize.NA.poisson = FALSE, seed = 7)
  assert(is.data.frame(x) && nrow(x) == 25L)
})

run_probe("map.env.variable.names", {
  x <- NicheDiv::map.env.variable.names(c(Env1 = 1, Env2 = 2), name.length = "short")
  assert(length(x) == 2L && is.numeric(x))
})

run_probe("install_nichediv_dependencies", {
  ## Network-free no-install branch; actual downloads are not exercised here.
  x <- NicheDiv::install_nichediv_dependencies(
    install_climatena = FALSE, install_whitebox = FALSE,
    install_data_table = FALSE, verbose = FALSE
  )
  assert(isTRUE(x))
})

run_probe("crop.background.buffered", {
  x <- NicheDiv::crop.background.buffered(
    occurrence.data = occ1,
    background.data = bg1,
    buffer.method = "hull", buffer.dist.meters = 30000,
    verbose = FALSE
  )
  assert(is.data.frame(x) && nrow(x) > 0L)
})

run_probe("thin.occurrence", {
  x <- NicheDiv::thin.occurrence(
    occurrence.data = occ1,
    thinning.dist.km = 0.1,
    N.thinning.replicates = 2,
    calc.Morans.I = FALSE,
    verbose = FALSE
  )
  assert(is.data.frame(x) && nrow(x) >= 2L)
})

run_probe("transform.skewed.variables", {
  x <- NicheDiv::transform.skewed.variables(
    data.frame = occ,
    background.dataframe = data.frame(Species = species, occ[, cols_env]),
    exclude.cols = "Species", verbose = FALSE
  )
  assert(is.data.frame(x$transformed) && is.data.frame(x$background.transformed))
})

run_probe("remove.low.CV.vars", {
  x <- NicheDiv::remove.low.CV.vars(
    Sp1.occurrence.data = occ1, Sp2.occurrence.data = occ2,
    Sp1.background.data = bg1, Sp2.background.data = bg2,
    exclude.cols = c("Species", "Longitude", "Latitude"),
    CV.threshold = 0.01, verbose = FALSE
  )
  assert(is.list(x) && length(x$kept_variables) >= 1L)
})

run_probe("filter.analogous.variables", {
  x <- NicheDiv::filter.analogous.variables(
    Sp1.background.data = bg1, Sp2.background.data = bg2,
    Sp1.Sp2.occurrence.data = occ,
    exclude.cols = "Species", min.rows = 15,
    overlap.threshold = 0.2, plot.1D.overlap = FALSE,
    use.parallel = FALSE, verbose = FALSE
  )
  assert(is.data.frame(x) && nrow(x) == nrow(occ) && "Species" %in% names(x))
})

run_probe("trim.to.analogous.environments", {
  x <- NicheDiv::trim.to.analogous.environments(
    Sp1.occurrence.data = occ1, Sp2.occurrence.data = occ2,
    Sp1.background.data = bg1, Sp2.background.data = bg2,
    exclude.cols = c("Species", "Longitude", "Latitude"),
    keep.occurrence.cols = c("Species", "Longitude", "Latitude"),
    grid.resolution = 30, analogous.window.size = 5,
    verbose = FALSE
  )
  assert(is.data.frame(x) && nrow(x) > 0L)
})

## 2. Custom raster extraction, with no remote data downloads
run_probe("extract.env.and.background", {
  dir_run <- tempfile("nichediv_ci_extract_")
  dir.create(dir_run, recursive = TRUE)
  raster_file <- file.path(dir_run, "synthetic.tif")
  xras <- terra::rast(nrows = 40, ncols = 40, xmin = -85, xmax = -84,
                      ymin = 37.7, ymax = 38.5, crs = "EPSG:4326")
  terra::values(xras) <- seq_len(terra::ncell(xras))
  names(xras) <- "SyntheticEnv"
  terra::writeRaster(xras, raster_file, overwrite = TRUE)
  input <- occ_coord[1:8, c("Longitude", "Latitude", "Species")]
  rownames(input) <- paste0("occ", seq_len(nrow(input)))
  NicheDiv::extract.env.and.background(
    occurrence.data = input,
    generate.background.data = FALSE,
    buffer.km = 5,
    csv.occurrence.out.file = "extracted.csv",
    output.dir = file.path(dir_run, "output"),
    env.datasets = NULL,
    custom.env.rasters = raster_file,
    custom.env.rasters.names = "synthetic",
    custom.env.rasters.variable.names = "SyntheticEnv",
    verbose = FALSE
  )
  output <- utils::read.csv(file.path(dir_run, "output", "extracted.csv"),
                            check.names = FALSE)
  assert(nrow(output) == 8L && "SyntheticEnv" %in% names(output))
  unlink(dir_run, recursive = TRUE, force = TRUE)
})

## 3. DAPC, niche metrics and visualizations
fit <- NULL
run_probe("run.DAPC.crossval.permutation", {
  fit <- NicheDiv::run.DAPC.crossval.permutation(
    data.input = occ,
    N.crossval.replicates = 2,
    N.permutations = 9,
    N.heldout.replicates = 2,
    N.bootstrap.replicates = 0,
    use.parallel = "none", N.cores = 1,
    seed = 4, verbose = FALSE
  )
  assert(is.list(fit) && is.list(fit$dapc_results))
  assert(is.finite(fit$observed_assign_prop))
})

run_probe("calc.niche.divergence.metrics", {
  assert(!is.null(fit), "DAPC fit was unavailable")
  raw <- NicheDiv::calc.niche.divergence.metrics(fit, verbose = FALSE)
  assert(is.finite(raw$Schoener_D) && raw$Schoener_D >= 0 && raw$Schoener_D <= 1)
  weighted <- NicheDiv::calc.niche.divergence.metrics(
    fit, weight.background = TRUE,
    Sp1.background.data = bg1, Sp2.background.data = bg2,
    verbose = FALSE
  )
  assert(is.finite(weighted$Schoener_D) && weighted$Schoener_D >= 0 && weighted$Schoener_D <= 1)
})

run_probe("plot.DAPC.niche.divergence", {
  assert(!is.null(fit), "DAPC fit was unavailable")
  p <- NicheDiv::plot.DAPC.niche.divergence(fit, show.plot = FALSE, verbose = FALSE)
  assert_ggplot(p)
})

run_probe("plot.DAPC.permutation", {
  assert(!is.null(fit), "DAPC fit was unavailable")
  p <- NicheDiv::plot.DAPC.permutation(fit, show.plot = FALSE, verbose = FALSE)
  assert_ggplot(p)
})

run_probe("plot.DAPC.var.contributions", {
  assert(!is.null(fit), "DAPC fit was unavailable")
  x <- NicheDiv::plot.DAPC.var.contributions(fit, show.plot = FALSE, verbose = FALSE)
  assert(is.data.frame(x) && nrow(x) >= 1L)
})

run_probe("plot.top.DAPC.predictors", {
  assert(!is.null(fit), "DAPC fit was unavailable")
  p <- NicheDiv::plot.top.DAPC.predictors(
    dapc.results = fit,
    predictor.data = occ[, cols_env, drop = FALSE],
    species.labels = occ$Species,
    N.top.variables = 2,
    show.plot = FALSE, verbose = FALSE
  )
  assert_ggplot(p)
})

run_probe("plot.occurrences.map", {
  out <- tempfile(fileext = ".pdf")
  grDevices::pdf(out, width = 7, height = 6)
  tryCatch(
    NicheDiv::plot.occurrences.map(
      coordinates = coords, group.labels = species,
      plot.background.points = FALSE, add.USA.states = FALSE,
      add.USA.counties = FALSE, show.plot = TRUE,
      save = FALSE, verbose = FALSE
    ),
    finally = grDevices::dev.off()
  )
  assert(file.exists(out) && file.info(out)$size > 0)
  unlink(out)
})

## Report every function and fail the GitHub check if any probe failed
report <- results[order(results$Function), , drop = FALSE]
utils::write.csv(report, "function-smoke-report.csv", row.names = FALSE)
print(report, row.names = FALSE)

## Add a summary table to the GitHub Actions run page
summary_file <- Sys.getenv("GITHUB_STEP_SUMMARY", "")
if (nzchar(summary_file)) {
  markdown <- c(
    "## NicheDiv exported-function smoke report",
    "",
    paste0("Platform: `", R.version$platform, "`"),
    "",
    "| Export | Status | Seconds | Detail |",
    "|---|---|---:|---|",
    vapply(seq_len(nrow(report)), function(i) {
      paste0("| `", report$Function[i], "` | ", report$Status[i],
             " | ", report$Seconds[i], " | ",
             gsub("[\r\n|]", " ", report$Detail[i]), " |")
    }, character(1))
  )
  cat(paste(markdown, collapse = "\n"), file = summary_file, append = TRUE, sep = "\n")
}

if (anyDuplicated(results$Function) ||
    !setequal(results$Function, expected_exports)) {
  stop("Smoke test manifest did not exercise every exported function exactly once")
}

failures <- report[report$Status != "PASS", , drop = FALSE]
if (nrow(failures) > 0L) {
  stop(nrow(failures), " of ", nrow(report),
       " exported function probes failed. See function-smoke-report.csv")
}

cat("\nAll ", nrow(report), " exported function smoke probes passed on ",
    R.version$platform, ".\n", sep = "")
