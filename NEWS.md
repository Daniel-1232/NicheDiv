# NicheDiv 1.2.1 (2026-10-28)
## Improve stochastic reproducibility

- Added a configurable `downsampling.seed` argument to `trim.to.analogous.environments()` for reproducible occurrence downsampling.
- Forwarded the downsampling seed to both species-specific calls of `sample.down()`.
- Added explicit deterministic seeding of the random-chance baseline in DAPC cross-validation.
- Updated function documentation to describe the new seed argument.

# NicheDiv 1.2.0 (2026-10-26)

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

