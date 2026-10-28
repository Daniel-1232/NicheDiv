#### Verify exported functions ###############################################
test_that("all expected NicheDiv functions are exported", {
  expected_functions <- c(
    "calc.niche.divergence.metrics",
    "convert.integer.to.numeric",
    "crop.background.buffered",
    "extract.env.and.background",
    "filter.analogous.variables",
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
  exported_functions <- getNamespaceExports("NicheDiv")
  expect_true(all(expected_functions %in% exported_functions))
  for (function_name in expected_functions) {
    expect_true(is.function(getExportedValue("NicheDiv", function_name)))
  }
})

#### Verify seed arguments ####################################################
test_that("stochastic functions expose seed arguments", {
  expected_seed_arguments <- c(
    "extract.env.and.background" = "seed",
    "filter.analogous.variables" = "seed",
    "run.DAPC.crossval.permutation" = "seed",
    "sample.down" = "seed",
    "thin.occurrence" = "seed",
    "trim.to.analogous.environments" = "downsampling.seed"
  )
  for (function_name in names(expected_seed_arguments)) {
    function_arguments <- names(formals(getExportedValue("NicheDiv", function_name)))
    expect_true(expected_seed_arguments[[function_name]] %in% function_arguments, info = function_name)
  }
})

#### Test reproducible occurrence sampling ###################################
test_that("sample.down reproduces weighted and uniform sampling", {
  occurrence_data <- data.frame(
    Variable1 = seq_len(100),
    Variable2 = seq_len(100) / 10,
    Variable3 = rep(c(NA_real_, 1, 2, 3), length.out = 100)
  )
  for (weighted_sampling in c(TRUE, FALSE)) {
    first_result <- NicheDiv::sample.down(
      occurrence_data, N.rows = 30,
      prioritize.NA.poisson = weighted_sampling, seed = 123
    )
    second_result <- NicheDiv::sample.down(
      occurrence_data, N.rows = 30,
      prioritize.NA.poisson = weighted_sampling, seed = 123
    )
    expect_identical(first_result, second_result)
    expect_equal(nrow(first_result), 30)
    expect_true(all(first_result$Variable1 %in% occurrence_data$Variable1))
    expect_equal(length(unique(first_result$Variable1)), 30)
  }
})

#### Test reproducible spatial thinning ######################################
test_that("thin.occurrence reproduces thinning", {
  occurrence_data <- data.frame(
    Longitude = rep(seq(-84.5, -84.49, length.out = 12), each = 2),
    Latitude = rep(seq(38.0, 38.01, length.out = 12), each = 2),
    Environment = seq_len(24)
  )
  first_result <- NicheDiv::thin.occurrence(
    occurrence_data, thinning.dist.km = 0.1,
    N.thinning.replicates = 5, calc.Morans.I = FALSE,
    plot.Morans.I = FALSE, seed = 123, verbose = FALSE
  )
  second_result <- NicheDiv::thin.occurrence(
    occurrence_data, thinning.dist.km = 0.1,
    N.thinning.replicates = 5, calc.Morans.I = FALSE,
    plot.Morans.I = FALSE, seed = 123, verbose = FALSE
  )
  expect_identical(first_result, second_result)
  expect_gt(nrow(first_result), 0)
  expect_lt(nrow(first_result), nrow(occurrence_data))
})

#### Test deterministic data conversion ######################################
test_that("convert.integer.to.numeric preserves data", {
  input_data <- data.frame(
    Integer = 1:5,
    Numeric = seq(0.1, 0.5, by = 0.1),
    Group = letters[1:5]
  )
  converted_data <- NicheDiv::convert.integer.to.numeric(input_data)
  expect_type(converted_data$Integer, "double")
  expect_identical(converted_data$Numeric, input_data$Numeric)
  expect_identical(converted_data$Group, input_data$Group)
  expect_identical(dim(converted_data), dim(input_data))
})
