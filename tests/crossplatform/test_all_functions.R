################################################################################
#### Settings
################################################################################
set.seed(137)
report_dir <- file.path("tests", "crossplatform", "results")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
report_file <- file.path(report_dir, paste0("all_functions_", Sys.info()[["sysname"]], "_", R.version$arch, ".csv"))

################################################################################
#### Load NicheDiv
################################################################################
pkgload::load_all(".", quiet = TRUE)
expected_functions <- c(
  "convert.integer.to.numeric",
  "crop.background.buffered",
  "sample.down",
  "thin.occurrence",
  "transform.skewed.variables",
  "remove.low.CV.vars",
  "filter.analogous.variables",
  "trim.to.analogous.environments",
  "run.DAPC.crossval.permutation",
  "plot.DAPC.niche.divergence",
  "plot.DAPC.permutation",
  "calc.niche.divergence.metrics",
  "plot.DAPC.var.contributions",
  "plot.top.DAPC.predictors",
  "plot.occurrences.map",
  "map.env.variable.names",
  "extract.env.and.background",
  "install_nichediv_dependencies"
)
exports <- getNamespaceExports("NicheDiv")
missing_exports <- setdiff(expected_functions, exports)
additional_exports <- setdiff(exports, expected_functions)
if (length(missing_exports)) stop("Expected NicheDiv functions not exported: ", paste(missing_exports, collapse = ", "))
if (length(additional_exports)) stop("Unlisted exported NicheDiv functions need tests: ", paste(additional_exports, collapse = ", "))
if (!all(vapply(expected_functions, function(x) is.function(getExportedValue("NicheDiv", x)), logical(1)))) stop("One or more exports are not functions")

################################################################################
#### Test harness
################################################################################
results <- data.frame(function_name = character(), test_type = character(), status = character(), details = character(), stringsAsFactors = FALSE)
run_case <- function(name, test_type, expr) {
  message("TEST: ", name, " [", test_type, "]")
  error_msg <- tryCatch({
    force(expr)
    NA_character_
  }, error = function(e) conditionMessage(e))
  ok <- is.na(error_msg)
  results[nrow(results) + 1L, ] <<- list(name, test_type, if (ok) "PASS" else "FAIL", if (ok) "" else error_msg)
  if (!ok) message("FAILED: ", error_msg)
  invisible(ok)
}

################################################################################
#### Reproducible synthetic data
################################################################################
n <- 50
background1 <- data.frame(
  MAT = rnorm(n, 14, 3),
  MAP = rnorm(n, 1200, 170),
  CMD = rnorm(n, 250, 45),
  RH = rnorm(n, 65, 9)
)
background2 <- data.frame(
  MAT = rnorm(n, 14.2, 3),
  MAP = rnorm(n, 1205, 170),
  CMD = rnorm(n, 255, 45),
  RH = rnorm(n, 65.2, 9)
)
occurrence1 <- background1[1:35, , drop = FALSE]
occurrence2 <- background2[1:35, , drop = FALSE]
combined_occurrences <- rbind(occurrence1, occurrence2)
combined_occurrences$Species <- factor(rep(c("Sp1", "Sp2"), each = 35))
coords <- data.frame(
  Longitude = c(-84.5, -84.4, -84.3, -84.2, -84.1, -84.0),
  Latitude = c(38.0, 38.1, 38.2, 38.0, 38.1, 38.2)
)

################################################################################
#### Data preprocessing functions
################################################################################
run_case("convert.integer.to.numeric", "calculation", {
  z <- convert.integer.to.numeric(data.frame(a = 1:5, b = letters[1:5]))
  stopifnot(is.double(z$a), identical(z$b, letters[1:5]))
})
run_case("crop.background.buffered", "calculation", {
  bg <- rbind(coords, data.frame(Longitude = -110, Latitude = 25))
  z <- crop.background.buffered(coords, bg, buffer.dist.meters = 25000, buffer.method = "points", verbose = FALSE)
  stopifnot(is.data.frame(z), nrow(z) >= 1L, nrow(z) < nrow(bg))
})
run_case("sample.down", "calculation", {
  z <- sample.down(background1, N.rows = 20, prioritize.NA.poisson = FALSE, seed = 1)
  stopifnot(nrow(z) == 20L, identical(names(z), names(background1)))
})
run_case("thin.occurrence", "calculation", {
  locations <- data.frame(Longitude = seq(-84.5, -83.6, length.out = 40), Latitude = rep(c(38.0, 38.2), 20))
  z <- thin.occurrence(locations, thinning.dist.km = 1, N.thinning.replicates = 3, calc.Morans.I = FALSE, plot.Morans.I = FALSE, verbose = FALSE)
  stopifnot(is.data.frame(z), nrow(z) >= 1L, nrow(z) <= nrow(locations))
})
run_case("transform.skewed.variables", "calculation", {
  z <- transform.skewed.variables(data.frame(MAT = exp(rnorm(50)), MAP = rnorm(50, 1200, 80)), verbose = FALSE)
  stopifnot(!is.null(z))
})
run_case("remove.low.CV.vars", "calculation", {
  z <- remove.low.CV.vars(occurrence1, occurrence2, background1, background2, verbose = FALSE)
  stopifnot(is.list(z), length(z$kept_variables) >= 2L)
})
run_case("filter.analogous.variables", "calculation", {
  z <- filter.analogous.variables(background1, background2, combined_occurrences,
                                 exclude.cols = "Species", overlap.threshold = 0.4,
                                 min.rows = 10, plot.1D.overlap = FALSE,
                                 use.parallel = FALSE, max.pairs = 100, verbose = FALSE)
  stopifnot(!is.null(z))
})
run_case("trim.to.analogous.environments", "calculation", {
  shared_grid <- expand.grid(
    MAT = seq(10, 18, length.out = 17),
    MAP = seq(900, 1500, length.out = 17)
  )
  z <- trim.to.analogous.environments(
    shared_grid, shared_grid, shared_grid, shared_grid,
    grid.resolution = 10, analogous.window.size = 8,
    downsample.equal.sizes = FALSE, verbose = FALSE
  )
  stopifnot(is.data.frame(z), nrow(z) > 0L,
            identical(names(z), names(shared_grid)),
            !anyNA(z))
})
run_case("map.env.variable.names", "calculation", {
  z <- map.env.variable.names(data.frame(MAT = 1:3, MAP = 4:6), name.length = "short")
  stopifnot(is.data.frame(z), nrow(z) == 3L)
})

################################################################################
#### DAPC calculations and plots
################################################################################
dapc <- NULL
run_case("run.DAPC.crossval.permutation", "calculation", {
  dapc <<- run.DAPC.crossval.permutation(combined_occurrences,
                                         N.crossval.replicates = 3,
                                         N.permutations = 10,
                                         fixed.n.pcs = 2,
                                         use.parallel = "none",
                                         seed = 137,
                                         verbose = FALSE)
  stopifnot(is.list(dapc), is.list(dapc$dapc_results), is.finite(dapc$observed_assign_prop))
})
if (!is.null(dapc)) {
  run_case("calc.niche.divergence.metrics", "calculation", {
    z <- calc.niche.divergence.metrics(dapc, combined_occurrences$Species, density.grid.resolution = 256, verbose = FALSE)
    stopifnot(!is.null(z))
  })
  run_case("plot.DAPC.niche.divergence", "plot construction", {
    z <- plot.DAPC.niche.divergence(dapc, show.plot = FALSE, verbose = FALSE)
    stopifnot(!is.null(z))
  })
  run_case("plot.DAPC.permutation", "plot construction", {
    z <- plot.DAPC.permutation(dapc, show.plot = FALSE, verbose = FALSE)
    stopifnot(!is.null(z))
  })
  run_case("plot.DAPC.var.contributions", "plot construction", {
    z <- plot.DAPC.var.contributions(dapc, show.plot = FALSE, verbose = FALSE)
    stopifnot(!is.null(z))
  })
  run_case("plot.top.DAPC.predictors", "plot construction", {
    z <- plot.top.DAPC.predictors(dapc,
                                  predictor.data = combined_occurrences[, c("MAT", "MAP", "CMD", "RH")],
                                  species.labels = combined_occurrences$Species,
                                  N.top.variables = 4,
                                  show.plot = FALSE,
                                  verbose = FALSE)
    stopifnot(!is.null(z))
  })
} else {
  dependent_tests <- c("calc.niche.divergence.metrics", "plot.DAPC.niche.divergence", "plot.DAPC.permutation", "plot.DAPC.var.contributions", "plot.top.DAPC.predictors")
  for (name in dependent_tests) results[nrow(results) + 1L, ] <- list(name, "calculation dependent on DAPC", "BLOCKED", "DAPC failed")
}
run_case("plot.occurrences.map", "plot construction", {
  plot.occurrences.map(coords,
                            group.labels = rep(c("Sp1", "Sp2"), each = 3),
                            plot.background.points = FALSE,
                            add.USA.states = FALSE,
                            show.plot = FALSE,
                            verbose = FALSE)
})

################################################################################
#### ClimateNA export and installer (validated separately)
################################################################################
run_case("install_nichediv_dependencies", "CRAN dependency path; ClimateNAr installation tested separately", {
  install_nichediv_dependencies(install_climatena = FALSE, install_whitebox = FALSE,
                                install_data_table = TRUE, force = FALSE, verbose = FALSE)
  stopifnot(requireNamespace("data.table", quietly = TRUE))
})
results[nrow(results) + 1L, ] <- list("extract.env.and.background", "real-data + synthetic background integration", "SEPARATE", "Executed by test_climatena.R and test_background.R in this workflow")

################################################################################
#### Report results
################################################################################
write.csv(results, report_file, row.names = FALSE)
print(results, row.names = FALSE)
message("Report: ", report_file)
failed <- results$function_name[!results$status %in% c("PASS", "SEPARATE")]
if (length(failed)) stop("Function tests failed/blocked: ", paste(failed, collapse = ", "))
message("PASS: ", sum(results$status == "PASS"), " function test entries completed on ", Sys.info()[["sysname"]])
