################################################################################
#### Settings
################################################################################
set.seed(137)
pkgload::load_all(".", quiet = TRUE)
backend <- if (.Platform$OS.type == "windows") "Windows PSOCK" else "Unix multicore"
message("Testing parallel backend: ", backend)
################################################################################
#### Reproducible DAPC dataset
################################################################################
n <- 40L
species1 <- data.frame(MAT = rnorm(n, 13, 2), MAP = rnorm(n, 1200, 140),
                       CMD = rnorm(n, 240, 35), RH = rnorm(n, 64, 7))
species2 <- data.frame(MAT = rnorm(n, 15, 2), MAP = rnorm(n, 1280, 140),
                       CMD = rnorm(n, 255, 35), RH = rnorm(n, 67, 7))
occurrences <- rbind(species1, species2)
occurrences$Species <- factor(rep(c("Sp1", "Sp2"), each = n))
################################################################################
#### Test parallel DAPC permutations
################################################################################
result <- run.DAPC.crossval.permutation(
  occurrences, N.crossval.replicates = 3, N.permutations = 8,
  fixed.n.pcs = 2, use.parallel = "auto", N.cores = 2,
  seed = 137, verbose = FALSE
)
stopifnot(is.list(result), is.finite(result$observed_assign_prop),
          length(result$permutation_assign_props) == 8L,
          all(is.finite(result$permutation_assign_props)))
################################################################################
#### Test parallel environmental overlap calculation
################################################################################
overlap <- filter.analogous.variables(
  species1, species2, occurrences,
  exclude.cols = "Species", overlap.threshold = 0.1,
  min.rows = 10, plot.1D.overlap = FALSE,
  use.parallel = TRUE, N.cores = 2, max.pairs = 100,
  seed = 137, verbose = FALSE
)
stopifnot(!is.null(overlap))
message("PASS: Native parallel DAPC and overlap calculations on ", backend)
