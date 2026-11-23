# NicheDiv 1.3.4 (2026-11-23)

## Improvements

- Refined `plot.DAPC.var.contributions()` to display bootstrap contribution distributions as horizontal boxplots with narrower boxes (`width = 0.6`), increased opacity (`alpha = 0.9`), and outliers excluded from the plot.
- Added the `diamond.size` argument (default `1.9`) to control the size of white-filled, black-outlined diamonds marking the original full-data variable contributions.
- Retained ordering of bootstrap boxplots by median contribution and preserved percentile-based 95% bootstrap confidence intervals (`CI_lower` and `CI_upper`) in the returned data.
- Changed the default `equal.color` to `white` for predictors with negligible directional differences according to Cliff's delta.
- Excluded the `Equal` category from the plot legend while retaining its configurable color in the plot.
- Updated the variable-contribution function documentation and retained the original barplot behavior when `N.bootstraps.CI = NULL`.
- Required `ggplot2` version 3.5.0 or newer for the updated bootstrap boxplots.

# NicheDiv 1.3.3 (2026-11-20)

## Improvements

- Updated `plot.DAPC.var.contributions()` to color environmental predictors according to the direction of Cliff's delta rather than differences in species medians.
- Predictors with negligible directional differences (absolute Cliff's delta <= 0.147) retain the `Equal` plotting category. Its color is customizable through `equal.color` (default `#1A1A1A`).
- Corrected alignment between original environmental predictor names and displayed variable names when calculating and plotting Cliff's delta, while preserving the existing DAPC contribution and bootstrap calculations.
- Updated the variable-contribution function documentation.

# NicheDiv 1.3.2 (2026-11-19)

## Improve DAPC variable-contribution plotting

- Replaced bootstrap distribution boxplots with horizontal 95% confidence intervals.
- Retained black diamonds to indicate original full-data variable contributions.
- Fixed ggplot warnings by applying separate color scales to standard bar plots and bootstrap confidence interval plots.

## Update documentation and README

- Updated function documentation to describe the new confidence interval visualization.
- Clarified the input requirements of `plot.DAPC.var.contributions()` and `plot.top.DAPC.predictors()`.
- Updated the README tutorial to demonstrate bootstrap confidence interval plotting.

# NicheDiv 1.3.2 (2026-11-18)

## Improve background-weighted niche divergence metrics

- Removed the reciprocal analogous-environment restriction from `calc.niche.divergence.metrics()` when `weight.background = TRUE`.
- Retained species-specific background correction along the DAPC discriminant axis (LD1), independently correcting occurrence densities for environmental availability.
- Background weighting now follows the full-environment correction configuration (`correct.env = TRUE`, `nae = "YES"`) in `humboldt.niche.similarity()`.
- Retained species-specific occurrence/background density ratios, peak normalization, and the 0.001 threshold for negligible corrected densities.
- Retained safeguards against invalid or zero corrected densities.
- Background-weighted metrics no longer require overlapping background environments between species.


# NicheDiv 1.3.1 (2026-11-17)

## Fix DAPC plotting namespace imports

- Corrected missing ggplot2 imports for `geom_segment`, `geom_point`, and `geom_text` used by `plot.DAPC.permutation()`.
- Imported `geom_boxplot` for the bootstrap confidence interval display in `plot.DAPC.var.contributions()`.
- No changes to DAPC model fitting, permutation statistics, or niche divergence calculations.

# NicheDiv 1.3 (2026-11-17)

## Revise DAPC cross-validation and PC selection

- Replaced repeated random training/validation splits with repeated stratified five-fold cross-validation (80% training and 20% validation per fold), fitting PCA only on the training observations for each fold.
- Introduced the one-standard-error rule: select the smallest evaluated PC count whose mean class-balanced classification error is within one standard error of the minimum; calculate the standard error across five-fold repetition means.
- Revised the two-stage PC search to use exhaustive evaluation for up to 13 PCs and coarse-to-fine adaptive searches for larger candidate ranges. Extend the fine-search boundary one PC at a time when the one-SE-selected count remains at that boundary, reusing previously evaluated results and identical CV splits.
- Reassess the PC count during each label or spatial background permutation and during the new held-out and contribution-bootstrap procedures.
- Changed the default `N.crossval.replicates` from 100 to 10, representing 50 fold evaluations per candidate PC count per stage, and strengthened PC-count constraints based on training-group sizes.

## Add held-out validation and confidence intervals for niche divergence

- Added optional `N.heldout.replicates` for repeated stratified 50/50 training/test splits. Select PCs by inner five-fold cross-validation and fit PCA and discriminant analysis exclusively on each training half; return predictions, discriminant scores, classification accuracy, adjusted Rand index, selected PC counts, and train/test indices for held-out records.
- Extended `calc.niche.divergence.metrics()` to calculate Schoener's D, niche dissimilarity (NDS), niche breadth exclusivity (NE), niche divergence magnitude (ND), and niche divergence angle (theta) for the held-out repetitions, with summaries across repetitions.
- Added optional `N.bootstrap.replicates` for bootstrap uncertainty of mean held-out niche metrics. Resample training and test records separately within groups, refit the selection and discrimination procedure, and calculate percentile-based 95% confidence intervals for D, NDS, NE, ND, and theta.
- Require sufficient successful replicates for bootstrap summaries and retain consistency checks for metric implementation and KDE grid settings. Held-out and bootstrap metric calculations are unweighted.

## Extend permutation inference and refine niche metrics

- Added observed and null-distribution calculations and permutation p-values for Schoener's D and ND, in addition to DAPC assignment accuracy. Use lower-tail tests for D and upper-tail tests for ND and assignment accuracy under label permutation.
- Extended the spatial background permutation workflow to assess D and ND as well as assignment accuracy in each direction, using lower-tail D, upper-tail ND, and two-sided accuracy comparisons.
- Revised numerical integration of niche dissimilarity by splitting grid intervals where the two interpolated density curves cross. Treat theta as undefined when ND equals zero.
- Refined background-environment projection checks and restricted background-corrected metric densities to analogous environmental conditions along the discriminant axis.

## Extend DAPC plots and variable-contribution uncertainty

- Expanded `plot.DAPC.permutation()` from an assignment-accuracy histogram to three panels displaying permutation null distributions for assignment accuracy, Schoener's D, and ND, with observed statistics and empirical p-values.
- Corrected the plot annotation for the smallest attainable empirical p-value based on the number of valid permutation replicates.
- Added optional `N.bootstraps.CI` to `plot.DAPC.var.contributions()` to bootstrap variable contributions, reselect PCs and refit DAPC, and return percentile-based 95% confidence limits (`CI_lower` and `CI_upper`).
- Display contribution-bootstrap distributions as boxplots with full-data estimates marked by black diamonds; when bootstrapping, rank predictors by bootstrap median contribution. Updated group-direction labels to reflect the group with the higher predictor median.

## Improve performance, reproducibility, and input validation

- Reused pooled within-group covariance and Cholesky calculations across nested PC counts; accelerated training-fold PCA where covariance eigendecomposition is suitable and improved standard-deviation calculations.
- Reused Windows parallel CV workers, improved held-out parallel dispatch, and dispatched Unix-parallel label permutations in a single batch only when more than 800 permutations are requested.
- Standardized L'Ecuyer-CMRG random-number generation for held-out sampling across execution backends and restored the caller's original RNG kind and state on function exit.
- Prevented nested Unix parallel CV within parallel held-out repetitions, removed unused held-out contribution calculations, and added explicit checks for failed CV, held-out, bootstrap, and permutation workers.
- Improved handling of factors with unused levels, validated background coordinates and predictors, and included occurrence coordinates in the saved-result signature for spatial background analyses.
- Expanded the contrasting-correlation diagnostic from a minimum of 25 to 12 observations per group and from a maximum of 1,000 to 2,000 predictors, while retaining correlation, distribution-similarity, and multiple-testing criteria.
- Added configurable permutation-progress reporting through `message.N.permutations` and declared the `matrixStats` dependency.

## Changes to function arguments

- Removed the `fixed.n.pcs` override from `run.DAPC.crossval.permutation()`; retained-PC selection is cross-validation based.
- Added `N.heldout.replicates`, `N.bootstrap.replicates`, and `message.N.permutations` to `run.DAPC.crossval.permutation()`.
- Added `N.bootstraps.CI` to `plot.DAPC.var.contributions()`.
- Updated `calc.niche.divergence.metrics()` to obtain group labels from the supplied DAPC result rather than a separate `group.assignment` argument. The former `density.grid.resolution` argument is no longer exposed; the metric grid is fixed internally at 1,024 points. Calls using these removed arguments must be updated.

# NicheDiv 1.2.4 (2026-11-05)
## Add diagnostic for contrasting environmental correlations
- Added a diagnostic to `run.DAPC.crossval.permutation()` that compares species-specific Pearson correlation matrices to identify contrasting relationships between environmental variables that may not be fully represented by DAPC and could affect the accuracy of inferred niche divergence metrics.
- Flag variable pairs with strong opposing correlations despite similar individual distributions, using thresholds for correlation strength, correlation differences, Kolmogorov-Smirnov distance, standardized mean differences, and standard deviation ratios.
- Added approximate Fisher-z tests with Benjamini-Hochberg false discovery rate correction to screen for differences in species-specific correlations.
- Added warnings identifying flagged variable pairs and their correlations, with recommendations to examine the affected variables and assess the sensitivity of niche divergence metrics by rerunning the analysis after excluding implicated variables.
- Skip the diagnostic when fewer than 25 samples per species or more than 1000 predictors are available, and exclude variables with negligible within-species variance.
- Updated function documentation to explain the diagnostic's methodological justification, thresholds, limitations, and recommendations for sensitivity analyses.
- The diagnostic does not modify environmental predictors, DAPC fitting, cross-validation, or permutation testing.

# NicheDiv 1.2.3 (2026-11-05)
## Allow DAPC analyses without permutation testing
- Added support for `N.permutations = NULL` or `0` in `run.DAPC.crossval.permutation()`.
- Allowed DAPC fitting and cross-validation without running standard label permutations or Humboldt-style spatial background permutations.
- Removed background-data and spatial-coordinate requirements when permutation testing is disabled.
- Return an empty permutation null distribution (`numeric(0)`) and an undefined permutation p-value (`NA_real_`) when no permutations are performed.
- Updated `plot.DAPC.permutation()` to report an informative error when plotting is requested without permutation results.
- Updated function documentation and the GitHub tutorial to describe optional permutation testing.

# NicheDiv 1.2.2 (2026-10-29)

## Improve background-weighted niche divergence calculations

- Add configurable missing-data handling to `calc.niche.divergence.metrics()` using `max.NA.prop` and `impute.NA.median`.
- Exclude background records exceeding the missingness threshold and impute remaining missing environmental values using species-specific background medians.
- Improve preprocessing of background environmental data before projection into the fitted PCA and DAPC discriminant space.

# NicheDiv 1.2.1 (2026-10-28)
## Improve stochastic reproducibility

- Added a configurable `downsampling.seed` argument to `trim.to.analogous.environments()` for reproducible occurrence downsampling.
- Forwarded the downsampling seed to both species-specific calls of `sample.down()`.
- Added explicit deterministic seeding of the random-chance baseline in DAPC cross-validation.
- Updated function documentation to describe the new seed argument.

# NicheDiv 1.2.0 (2026-10-28)

## Revise Humboldt-style analogous-environment trimming

- Substantially revised `trim.to.analogous.environments()` to implement reciprocal PCA-based environmental trimming following the Humboldt framework.
- Replaced PCA fitted to combined background environments with a shared PCA fitted to occurrence environments, assigning equal total weight to each species.
- Projected species-specific background environments as supplementary observations into the occurrence-weighted PCA space.
- Added estimation of background environmental densities on a shared two-dimensional PCA grid and calculation of background Schoener's D.
- Implemented overlap-dependent trimming: reciprocal PCA-range restriction when background Schoener's D exceeds 0.97, and additional reciprocal moving-window trimming otherwise.
- Replaced the previous five-iteration trimming procedure with reciprocal PCA-range restriction followed, when required, by four directional moving-window passes across PC1 and PC2.
- Increased the default PCA grid resolution from 50 to 100.
- Improved handling of sf occurrence and background datasets by excluding geometry columns from environmental calculations.
- Added species-specific median imputation for missing background predictors and pooled occurrence-median imputation for the trimming PCA.
- Applied missing-data thresholds to occurrence records and removed predictors with insufficient occurrence variation before fitting PCA.
- Restricted imputed occurrence values to the trimming analysis without replacing environmental values in returned occurrence records.
- Retained optional occurrence downsampling after trimming to balance species sample sizes.
- Improved input validation and reporting of retained occurrence records.

## Implement Humboldt-style spatially structured background permutations

- Substantially revised the optional `background.permutation.test` procedure in `run.DAPC.crossval.permutation()` to implement Humboldt-style spatial background randomization (Brown & Carnaval, 2019).
- Replaced random sampling from the opposite species' background with spatial translations of each species' occurrence configuration within its own geographic background.
- Added preservation and validation of longitude and latitude coordinates while excluding coordinates from environmental predictors used for DAPC.
- Implemented spatial translations in a metric Lambert azimuthal equal-area projection, with displacement limits based on one-third of the corresponding background extent.
- Restricted translated points to the convex hull of the corresponding species' background and matched valid points to nearest background records within 15 km.
- Added `RANN` as a dependency for nearest-neighbor environmental matching.
- Allowed up to 1,000 spatial translation attempts per permutation and species, pooling valid background records across shifts to reconstruct the original occurrence sample size.
- Recalculated PCA and LDA assignment accuracy for each spatial permutation, separately shifting species 1 and species 2 while keeping the other species unchanged.
- Replaced one-sided background permutation p-values with two-sided empirical p-values and required at least 80% valid permutations per direction.
- Renamed directional permutation outputs from `forward` and `reverse` to `species1_shifted` and `species2_shifted`.
- Improved the standard label-permutation p-value calculation by including ties and returning NA when no valid permutations are available.
- Adapted Humboldt-style spatial background randomization to DAPC assignment accuracy rather than Humboldt's original environmental niche-overlap statistics.

## Improve species-specific background correction of niche divergence metrics

- Substantially revised `calc.niche.divergence.metrics()` to implement species-specific Humboldt-style background correction along the DAPC discriminant axis.
- Replaced correction using a pooled background density with separate occurrence-to-background density ratios for each species.
- Added projection of species-specific background environments into the fitted PCA and LD1 space, including discrimination-axis centering and projection-consistency checks.
- Required complete, correctly ordered PCA environmental predictors in each background dataset.
- Extended the shared LD1 density grid to include projected background environments and additional kernel margins.
- Added kernel-tail filtering based on the minimum estimated density at observed occurrence and background records.
- Scaled occurrence and background densities to expected grid-cell counts before calculating species-specific density ratios.
- Added peak normalization, reciprocal analogous-environment masking, and removal of corrected density values below the environmental-support threshold.
- Improved KDE numerical stability by centering LD1 observations before bandwidth selection, resolving translation sensitivity.
- Improved validation of LD1 coordinates, group assignments, grid resolution, background predictors, and logical arguments.
- Preserved unweighted calculations as the default and retained Schoener's D, NDS, NE, ND, and divergence angle outputs.

# NicheDiv 1.1.0 (2026-10-21)

## Accessible-area delineation

- Added user-defined accessible-area polygons to crop.background.buffered().
- Added buffer.method = "custom" and the custom.polygon argument.
- Added support for sf, sfc, terra SpatVector, and spatial polygon files.
- Added automatic polygon CRS transformation and geometry validation.
- Supported multipolygons, disconnected regions, and polygon holes.
- Custom polygons are used directly without additional buffering.
- Preserved the existing hull, points, alpha, and bbox methods.
- Added permanent regression tests for custom polygon handling.

# NicheDiv 1.0.0 (2026-10-21)

## First stable release

- First stable release of the NicheDiv R package for multivariate niche divergence analyses using discriminant analysis of principal components (DAPC).
- Updated ClimateNAr integration to version 3.2.0.
- Added cross-platform ClimateNAr installation and extraction support for Windows, Linux, and macOS.
- Added GitHub Actions validation for Windows, Ubuntu Linux, macOS Intel, and macOS Apple Silicon.
- Improved package portability and cross-platform testing.

