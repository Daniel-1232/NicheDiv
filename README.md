# NicheDiv

NicheDiv is an R package for testing ecological niche divergence between two predefined groups, such as species, lineages, subspecies, populations, or genetic clusters, across highly multivariate environmental space. The approach is described in an associated manuscript that is currently in review.

The package adapts discriminant analysis of principal components (DAPC) to environmental niche data. Environmental variables are first transformed into principal components to reduce dimensionality and collinearity. Discriminant analysis is then used to identify the axis that best separates the two groups. NicheDiv evaluates the significance of this separation with a permutation test and summarizes niche divergence with interpretable metrics, including Schoener’s D, niche dissimilarity, niche breadth exclusivity, niche divergence magnitude, and niche divergence angle.

NicheDiv is designed for workflows that use many environmental predictors, including climatic, topographic, phenological, vegetation, soil, land-cover, anthropogenic, and user-supplied raster layers.

## Main features

* Extract environmental values for occurrence records and background points.
* Generate accessible-area background points from occurrence coordinates.
* Crop background points to group-specific buffered accessible areas.
* Spatially thin occurrence records to reduce spatial autocorrelation.
* Balance sample sizes between groups.
* Transform skewed environmental variables.
* Remove low-information variables based on coefficient of variation.
* Filter predictors to analogous environmental space.
* Run cross-validated DAPC with a permutation test.
* Calculate niche divergence metrics from the discriminant axis.
* Plot DAPC niche divergence, permutation null distributions, variable contributions, top predictors, and occurrence/background maps.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("Daniel-1232/NicheDiv")
```

Load the package:

```r
library(NicheDiv)
```

Several workflows also use common data-manipulation and raster packages:

```r
library(dplyr)
library(terra)
```

Large environmental extraction workflows may require substantial disk space and processing time, especially when downloading high-resolution rasters or generating large background datasets.

## Input data

The main input is an occurrence data frame containing at least:

* one row per occurrence record;
* unique row names or an ID column;
* longitude and latitude columns;
* one grouping column with two groups to compare;
* optional metadata columns;
* optional environmental variables if environmental extraction has already been performed.

Example required columns:

```r
head(occurrence_data[, c("Longitude", "Latitude", "Species")])
```

```text
  Longitude Latitude Species
1   -121.50    38.40     Sp1
2   -121.75    38.55     Sp1
3   -117.25    34.15     Sp2
4   -117.10    34.40     Sp2
```

NicheDiv can either extract environmental data internally from raster layers or use a data frame where environmental variables are already present.

## Recommended workflow

The typical NicheDiv workflow has seven major steps:

1. Set paths, group names, coordinate columns, and analysis parameters.
2. Extract environmental data and generate background points.
3. Prepare occurrence and background data.
4. Crop and downsample group-specific background data.
5. Spatially thin and balance occurrence records.
6. Transform, filter, and restrict predictors to analogous environmental space.
7. Run DAPC, calculate niche divergence metrics, and plot results.

The code below gives a compact version of the full workflow.

## 1. Set working environment and input parameters

```r
#### Set working environment and input #########################################

## Set directories
base_dir <- "path/to/project"
results_dir <- file.path(base_dir, "Results")
figure_dir <- file.path(results_dir, "Figure_files")
intermediate_files_dir_name <- "Intermediate_files"

if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


## Set input files
occurrence_data_file <- file.path(base_dir, "Data/occurrences.csv")

csv_occurrence_out_file <- "Occurrences_env.csv"
csv_background_out_file <- "Background_env.csv"


## Set group and coordinate columns
Sp1_name <- "Group_1"
Sp2_name <- "Group_2"
Sp1_label <- "Group 1"
Sp2_label <- "Group 2"

Species_col <- "Species"
Longitude_col <- "Longitude"
Latitude_col <- "Latitude"
CRS_all <- "EPSG:4326"


## Set analysis parameters
buffer_km <- 5
N_background_points <- 500000
N_background_sample <- 10000
N_permutations <- 1000
CV_threshold <- 0.01
base_colors <- c("#A331A3", "#6CB3A5")

exclude_cols <- c("SampleID", "Locality", "CollectionDate")
set.seed(1)
```

Optional custom rasters can be supplied as GeoTIFF files:

```r
custom_raster_path <- file.path(base_dir, "Data/custom_environmental_layers.tif")
custom_raster_variable_names <- list(names(terra::rast(custom_raster_path)))
```

## 2. Extract environmental data and generate background points

Use `extract.env.and.background()` to extract environmental values for occurrence records and generate background points within the accessible area.

```r
#### Extract environmental data ################################################

## Import occurrences
occurrence_data <- read.csv(occurrence_data_file)


## Extract environmental data and background points
NicheDiv::extract.env.and.background(occurrence.data = occurrence_data,
                                     longitude.col = Longitude_col,
                                     latitude.col = Latitude_col,
                                     generate.background.data = TRUE,
                                     N.background.points = N_background_points,
                                     buffer.km = buffer_km,
                                     remove.hydrolakes.background = FALSE,
                                     landmask.largest.N.pieces = 5,
                                     csv.occurrence.out.file = csv_occurrence_out_file,
                                     csv.background.out.file = csv_background_out_file,
                                     output.dir = results_dir,
                                     intermediate.files.dir = intermediate_files_dir_name,
                                     CRS.occurrences = CRS_all,
                                     overwrite = TRUE,
                                     custom.env.rasters = custom_raster_path,
                                     custom.env.rasters.variable.names = custom_raster_variable_names)
```

If environmental variables have already been extracted, this step can be skipped and the occurrence/background tables can be imported directly.

## 3. Import and prepare extracted data

```r
#### Import and prepare extracted data #########################################

## Import extracted occurrence and background data
Env_data_occurrences <- read.csv(file.path(results_dir, csv_occurrence_out_file))
Env_data_background <- read.csv(file.path(results_dir, csv_background_out_file), check.names = FALSE)


## Remove metadata columns not used as environmental predictors
Env_data_occurrences <- dplyr::select(Env_data_occurrences, -any_of(exclude_cols))
Env_data_background <- dplyr::select(Env_data_background, -any_of(exclude_cols))


## Convert integer columns to numeric
Env_data_occurrences <- NicheDiv::convert.integer.to.numeric(Env_data_occurrences)
Env_data_background <- NicheDiv::convert.integer.to.numeric(Env_data_background)


## Keep the two groups of interest
Env_data_occurrences <- Env_data_occurrences[Env_data_occurrences[[Species_col]] %in% c(Sp1_name, Sp2_name), , drop = FALSE]


## Rename groups for plotting
Env_data_occurrences[[Species_col]] <- dplyr::recode(Env_data_occurrences[[Species_col]],
                                                     !!!setNames(c(Sp1_label, Sp2_label), c(Sp1_name, Sp2_name)))
Env_data_occurrences[[Species_col]] <- factor(Env_data_occurrences[[Species_col]], levels = c(Sp1_label, Sp2_label))


## Split occurrence data by group
Sp1_occurrence_data <- Env_data_occurrences[Env_data_occurrences[[Species_col]] == Sp1_label, , drop = FALSE]
Sp2_occurrence_data <- Env_data_occurrences[Env_data_occurrences[[Species_col]] == Sp2_label, , drop = FALSE]
```

## 4. Crop and downsample background data

Background data should reflect the accessible environmental space for each group. The example below crops the shared background pool to a buffered convex hull around each group’s occurrence records and then downsamples each background to the same target size.

```r
#### Prepare background data ###################################################

## Crop background to each group-specific accessible area
Sp1_background_data <- NicheDiv::crop.background.buffered(occurrence.data = Sp1_occurrence_data,
                                                          background.data = Env_data_background,
                                                          CRS = CRS_all,
                                                          buffer.method = "hull",
                                                          buffer.dist.meters = buffer_km * 1000)

Sp2_background_data <- NicheDiv::crop.background.buffered(occurrence.data = Sp2_occurrence_data,
                                                          background.data = Env_data_background,
                                                          CRS = CRS_all,
                                                          buffer.method = "hull",
                                                          buffer.dist.meters = buffer_km * 1000)


## Downsample background data
Sp1_background_data <- NicheDiv::sample.down(Sp1_background_data, N.rows = N_background_sample)
Sp2_background_data <- NicheDiv::sample.down(Sp2_background_data, N.rows = N_background_sample)
```

Available background geometries are `"hull"`, `"points"`, `"alpha"`, and `"bbox"`. The convex hull is usually a robust default, but point buffers or alpha hulls can be useful for fragmented or spatially complex distributions.

## 5. Spatially thin and balance occurrence records

Spatial thinning reduces clustered sampling and residual spatial autocorrelation. After thinning, both groups are downsampled to the same number of occurrences.

```r
#### Spatial thinning and sample-size balancing ################################

## Thin occurrence records
Sp1_occurrence_thinned <- NicheDiv::thin.occurrence(Sp1_occurrence_data,
                                                    thinning.dist.km = 1)

Sp2_occurrence_thinned <- NicheDiv::thin.occurrence(Sp2_occurrence_data,
                                                    thinning.dist.km = 1)


## Downsample to equal sample size
n_min_occurrence_thinned <- min(nrow(Sp1_occurrence_thinned), nrow(Sp2_occurrence_thinned))

Sp1_occurrence_thinned <- NicheDiv::sample.down(Sp1_occurrence_thinned,
                                                N.rows = n_min_occurrence_thinned)

Sp2_occurrence_thinned <- NicheDiv::sample.down(Sp2_occurrence_thinned,
                                                N.rows = n_min_occurrence_thinned)
```

## 6. Transform and filter environmental variables

Environmental variables may be skewed, uninformative, or non-analogous between accessible areas. NicheDiv includes functions to handle these preprocessing steps before running DAPC.

```r
#### Transform skewed environmental variables ##################################

## Combine occurrence and background datasets
Sp1_Sp2_occurrence_thinned <- dplyr::bind_rows(Sp1_occurrence_thinned, Sp2_occurrence_thinned)

Sp1_background_data[[Species_col]] <- Sp1_label
Sp2_background_data[[Species_col]] <- Sp2_label
Sp1_Sp2_background_data <- dplyr::bind_rows(Sp1_background_data, Sp2_background_data)


## Transform skewed variables
transformation_results <- NicheDiv::transform.skewed.variables(data.frame = Sp1_Sp2_occurrence_thinned,
                                                               exclude.cols = c(Latitude_col, Longitude_col, Species_col, "ID"),
                                                               background.dataframe = Sp1_Sp2_background_data)

Sp1_Sp2_occurrence_transformed <- transformation_results$transformed
Sp1_Sp2_background_transformed <- transformation_results$background.transformed


## Split transformed data by group
Sp1_occurrence_transformed <- Sp1_Sp2_occurrence_transformed[Sp1_Sp2_occurrence_transformed[[Species_col]] == Sp1_label, , drop = FALSE]
Sp2_occurrence_transformed <- Sp1_Sp2_occurrence_transformed[Sp1_Sp2_occurrence_transformed[[Species_col]] == Sp2_label, , drop = FALSE]

Sp1_background_transformed <- Sp1_Sp2_background_transformed[Sp1_Sp2_background_transformed[[Species_col]] == Sp1_label, , drop = FALSE]
Sp2_background_transformed <- Sp1_Sp2_background_transformed[Sp1_Sp2_background_transformed[[Species_col]] == Sp2_label, , drop = FALSE]
```

Remove low-variation variables:

```r
#### Remove low-information variables ##########################################

CV_removal_results <- NicheDiv::remove.low.CV.vars(Sp1.occurrence.data = Sp1_occurrence_transformed,
                                                   Sp2.occurrence.data = Sp2_occurrence_transformed,
                                                   Sp1.background.data = Sp1_background_transformed,
                                                   Sp2.background.data = Sp2_background_transformed,
                                                   exclude.cols = c(Latitude_col, Longitude_col, Species_col, "ID"),
                                                   CV.threshold = CV_threshold)

Sp1_occurrence_filtered <- CV_removal_results$occurrence_Sp1
Sp2_occurrence_filtered <- CV_removal_results$occurrence_Sp2
Sp1_background_filtered <- CV_removal_results$background.Sp1
Sp2_background_filtered <- CV_removal_results$background.Sp2

Sp1_Sp2_occurrence_filtered <- dplyr::bind_rows(Sp1_occurrence_filtered, Sp2_occurrence_filtered)
```

Filter to analogous environmental variables:

```r
#### Filter to analogous environmental variables ###############################

Sp1_Sp2_analogous <- NicheDiv::filter.analogous.variables(Sp1.Sp2.occurrence.data = Sp1_Sp2_occurrence_filtered,
                                                          Sp1.background.data = Sp1_background_filtered,
                                                          Sp2.background.data = Sp2_background_filtered,
                                                          exclude.cols = c(Latitude_col, Longitude_col, Species_col),
                                                          CV.threshold = CV_threshold,
                                                          overlap.threshold = 0.7)
```

This step removes predictors whose background distributions are not sufficiently analogous between groups. This helps reduce bias caused by comparing groups across environmental conditions that are available to one group but not the other.

## 7. Run DAPC with cross-validation and permutation test

```r
#### Run DAPC niche divergence test ############################################

## Extract group assignments
Sp1_Sp2_species_assignment <- factor(Sp1_Sp2_analogous[[Species_col]])


## Set named group colors
Sp1_Sp2_species_colors <- setNames(base_colors[seq_along(levels(Sp1_Sp2_species_assignment))],
                                   levels(Sp1_Sp2_species_assignment))

Sp1_Sp2_species_assignment <- factor(Sp1_Sp2_species_assignment,
                                     levels = names(Sp1_Sp2_species_colors))


## Run cross-validated DAPC with permutation test
DAPC_results <- NicheDiv::run.DAPC.crossval.permutation(data.input = Sp1_Sp2_analogous,
                                                        species.col = Species_col,
                                                        exclude.cols = c(Latitude_col, Longitude_col),
                                                        N.permutations = N_permutations,
                                                        N.crossval.replicates = 300)
```

The permutation test compares the observed DAPC assignment accuracy to a null distribution generated by randomly permuting group labels. A significant result indicates that group separation along the discriminant axis is stronger than expected under random group membership.

## 8. Calculate niche divergence metrics

```r
#### Calculate niche divergence metrics ########################################

Niche_divergence_metrics <- NicheDiv::calc.niche.divergence.metrics(DAPC_results,
                                                                    group.assignment = Sp1_Sp2_species_assignment)

Niche_divergence_metrics
```

Optionally, calculate background-weighted metrics:

```r
Niche_divergence_metrics_weighted <- NicheDiv::calc.niche.divergence.metrics(DAPC_results,
                                                                             weight.background = TRUE,
                                                                             Sp1.background.data = Sp1_background_filtered,
                                                                             Sp2.background.data = Sp2_background_filtered,
                                                                             group.assignment = Sp1_Sp2_species_assignment)

Niche_divergence_metrics_weighted
```

The main metrics are:

* `Schoener_D`: niche overlap; 1 indicates complete overlap and 0 indicates no overlap.
* `Niche_dissimilarity`: density-based divergence along the discriminant axis.
* `Niche_breadth_exclusivity`: range exclusivity along the discriminant axis.
* `Niche_divergence_magnitude`: combined divergence magnitude in the niche divergence plane.
* `Niche_divergence_angle`: relative contribution of density-based versus range-based divergence.

## 9. Plot DAPC results

Plot the discriminant-axis density distributions:

```r
#### Plot DAPC niche divergence ################################################

NicheDiv::plot.DAPC.niche.divergence(DAPC_results,
                                     group.colors = Sp1_Sp2_species_colors,
                                     save = TRUE,
                                     overwrite = TRUE,
                                     type = "svg",
                                     output.dir = figure_dir,
                                     filename = "DAPC_niche_divergence",
                                     width = 16,
                                     height = 12)
```

Plot the permutation null distribution:

```r
#### Plot permutation test #####################################################

NicheDiv::plot.DAPC.permutation(DAPC_results,
                                save = TRUE,
                                overwrite = TRUE,
                                type = "svg",
                                output.dir = figure_dir,
                                filename = "DAPC_permutation_test",
                                width = 16,
                                height = 9)
```

Plot environmental variable contributions:

```r
#### Plot variable contributions ##############################################

DAPC_results_short_names <- DAPC_results
DAPC_results_short_names$dapc_results$var.contr <- NicheDiv::map.env.variable.names(DAPC_results_short_names$dapc_results$var.contr, "short")
DAPC_results_short_names$dapc_results$var.load <- NicheDiv::map.env.variable.names(DAPC_results_short_names$dapc_results$var.load, "short")

DAPC_var_contr <- NicheDiv::plot.DAPC.var.contributions(DAPC_results_short_names,
                                                        group.colors = Sp1_Sp2_species_colors,
                                                        save = TRUE,
                                                        overwrite = TRUE,
                                                        type = "svg",
                                                        output.dir = figure_dir,
                                                        filename = "DAPC_variable_contributions",
                                                        width = 16,
                                                        height = 10)

head(DAPC_var_contr)
```

Plot raw distributions of the top contributing predictors:

```r
#### Plot top predictors #######################################################

Sp1_Sp2_analogous_short_names <- NicheDiv::map.env.variable.names(Sp1_Sp2_analogous, "short")

NicheDiv::plot.top.DAPC.predictors(dapc.results = DAPC_results_short_names,
                                   predictor.data = Sp1_Sp2_analogous_short_names,
                                   species.labels = Sp1_Sp2_species_assignment,
                                   group.colors = Sp1_Sp2_species_colors,
                                   save = TRUE,
                                   overwrite = TRUE,
                                   type = "svg",
                                   output.dir = figure_dir,
                                   filename = "DAPC_top_predictors",
                                   width = 16,
                                   height = 10)
```

## 10. Plot occurrences and background points

```r
#### Plot occurrence and background map ########################################

background_labels <- factor(c(rep(levels(Sp1_Sp2_species_assignment)[1], nrow(Sp1_background_data)),
                              rep(levels(Sp1_Sp2_species_assignment)[2], nrow(Sp2_background_data))),
                            levels = levels(Sp1_Sp2_species_assignment))

background_data_combined <- dplyr::bind_rows(Sp1_background_data, Sp2_background_data)

NicheDiv::plot.occurrences.map(coordinates = Sp1_Sp2_analogous,
                               group.labels = Sp1_Sp2_species_assignment,
                               group.colors = unname(Sp1_Sp2_species_colors),
                               plot.background.points = TRUE,
                               background.coords = background_data_combined,
                               background.group.labels = background_labels,
                               legend.group.names = c(Sp1_label, Sp2_label),
                               save = TRUE,
                               overwrite = TRUE,
                               type = "svg",
                               output.dir = figure_dir,
                               filename = "Occurrence_background_map",
                               width = 16,
                               height = 12)
```

## Optional: Brown and Carnaval-style analogous trimming

In addition to variable-level analogy filtering, NicheDiv includes `trim.to.analogous.environments()` to remove occurrence records from non-analogous environmental conditions following the logic of Brown and Carnaval-style environmental analogy correction.

```r
#### Optional Brown and Carnaval-style correction ##############################

Sp1_Sp2_analogous_trimmed <- NicheDiv::trim.to.analogous.environments(Sp1.occurrence.data = Sp1_occurrence_filtered,
                                                                      Sp2.occurrence.data = Sp2_occurrence_filtered,
                                                                      Sp1.background.data = Sp1_background_filtered,
                                                                      Sp2.background.data = Sp2_background_filtered,
                                                                      exclude.cols = c(Latitude_col, Longitude_col, Species_col),
                                                                      keep.occurrence.cols = c(Latitude_col, Longitude_col, Species_col))
```

The trimmed dataset can then be passed to `run.DAPC.crossval.permutation()` using the same DAPC workflow shown above.

## Optional: run DAPC without analogous-variable filtering

For comparison, users may also run the DAPC test on the filtered occurrence data before analogous-variable filtering. This can help evaluate how much non-analogous environmental space affects the final result.

```r
#### Optional DAPC without analogous-variable filtering ########################

Sp1_Sp2_species_assignment_no_analogy <- factor(Sp1_Sp2_occurrence_filtered[[Species_col]])

Sp1_Sp2_species_colors_no_analogy <- setNames(base_colors[seq_along(levels(Sp1_Sp2_species_assignment_no_analogy))],
                                              levels(Sp1_Sp2_species_assignment_no_analogy))

Sp1_Sp2_species_assignment_no_analogy <- factor(Sp1_Sp2_species_assignment_no_analogy,
                                                levels = names(Sp1_Sp2_species_colors_no_analogy))

DAPC_results_no_analogy <- NicheDiv::run.DAPC.crossval.permutation(data.input = Sp1_Sp2_occurrence_filtered,
                                                                   species.col = Species_col,
                                                                   exclude.cols = c(Latitude_col, Longitude_col),
                                                                   N.permutations = N_permutations,
                                                                   N.crossval.replicates = 300)

Niche_divergence_metrics_no_analogy <- NicheDiv::calc.niche.divergence.metrics(DAPC_results_no_analogy,
                                                                               group.assignment = Sp1_Sp2_species_assignment_no_analogy)
```

## Interpreting results

A complete NicheDiv analysis should be interpreted using three complementary outputs.

First, the permutation test evaluates whether group separation is stronger than expected under random group labels. A significant result indicates that the groups are more environmentally separable than expected by chance.

Second, the niche divergence metrics quantify the amount and type of divergence along the discriminant axis. Schoener’s D describes overlap, niche dissimilarity describes density-based divergence, niche breadth exclusivity describes range exclusivity, niche divergence magnitude summarizes overall divergence in the niche divergence plane, and niche divergence angle describes whether divergence is driven more by density differences, range exclusivity, or both.

Third, the variable-contribution plots identify environmental predictors most associated with the discriminant separation. These contributions should be interpreted as environmental associations with the discriminant axis, not as direct evidence of causal ecological mechanisms.

## Main functions

| Function                           | Purpose                                                                          |
| ---------------------------------- | -------------------------------------------------------------------------------- |
| `extract.env.and.background()`     | Extract environmental variables and generate background data.                    |
| `convert.integer.to.numeric()`     | Convert integer columns to numeric.                                              |
| `crop.background.buffered()`       | Crop background points to buffered accessible areas.                             |
| `sample.down()`                    | Downsample occurrence or background records.                                     |
| `thin.occurrence()`                | Spatially thin occurrence records and evaluate residual spatial autocorrelation. |
| `transform.skewed.variables()`     | Transform skewed environmental variables.                                        |
| `remove.low.CV.vars()`             | Remove variables with low coefficient of variation.                              |
| `filter.analogous.variables()`     | Retain predictors with analogous background distributions.                       |
| `trim.to.analogous.environments()` | Remove occurrence records from non-analogous environmental conditions.           |
| `run.DAPC.crossval.permutation()`  | Run cross-validated DAPC and permutation testing.                                |
| `calc.niche.divergence.metrics()`  | Calculate Schoener’s D and niche divergence plane metrics.                       |
| `plot.DAPC.niche.divergence()`     | Plot density distributions along the DAPC discriminant axis.                     |
| `plot.DAPC.permutation()`          | Plot the permutation null distribution.                                          |
| `plot.DAPC.var.contributions()`    | Plot variable contributions to the discriminant axis.                            |
| `plot.top.DAPC.predictors()`       | Plot raw distributions of top contributing predictors.                           |
| `plot.occurrences.map()`           | Plot occurrence and background records on a map.                                 |
| `map.env.variable.names()`         | Convert environmental variable names to shorter or more readable labels.         |

## Notes and recommendations

Use biologically meaningful groups. NicheDiv tests divergence between predefined groups and does not infer species limits by itself.

Use a biologically justified accessible area. Background data should represent the environmental conditions plausibly available to each group.

Balance occurrence sample sizes before DAPC. Unequal sample sizes can affect discrimination and interpretation.

Inspect spatial thinning diagnostics. Thinning reduces spatial clustering but may not remove all spatial autocorrelation.

Use analogous-environment filtering when comparing groups from different accessible areas. This helps reduce bias from environmental conditions available to only one group.

Do not interpret variable contributions as proof of causation. They identify predictors associated with niche divergence and should be evaluated with biological knowledge, natural history, and independent evidence.

## Example output files

A standard analysis can produce:

```text
Results/
├── Occurrences_env.csv
├── Background_env.csv
└── Figure_files/
    ├── DAPC_niche_divergence.svg
    ├── DAPC_permutation_test.svg
    ├── DAPC_variable_contributions.svg
    ├── DAPC_top_predictors.svg
    └── Occurrence_background_map.svg
```

## Citation

The methodological framework implemented in NicheDiv is described in an associated manuscript that is currently in review.

Please cite NicheDiv and the associated manuscript when using the package:

```r
citation("NicheDiv")
```

If the package citation is not yet installed, cite the manuscript as:

```text
Schönberger, D., MacDonald, Z. G., Tuttle, J. P., Schmidt, B. C., & Dupuis, J. R. NicheDiv: A DAPC framework to quantify niche divergence across highly multivariate environmental space. Manuscript in review.
```

Please also include the package version or GitHub commit used for the analysis.

## License

Add license information here before public release.
