# NicheDiv R package

NicheDiv is an R package for testing niche divergence across highly multivariate environmental space between two predefined groups (e.g., species, lineages, populations, or genetic clusters). This is done by adapting discriminant analysis of principal components (DAPC) to environmental niche data. Environmental variables are first transformed into principal components (PCs) to reduce dimensionality and collinearity. Discriminant analysis is then used to identify the axis that best separates the two groups. NicheDiv summarizes niche divergence with easily interpretable metrics.


## Main advantages of the approach

- Automatically extracts environmental values for occurrence records and background points from implemented and user-supplied environmental layers.
- Implemented global environmental layers cover monthly to seasonal climate, topography, phenology, hyrdology, vegetation, soil, land cover, and anthropogenic variables.
- Mitigates common biases in niche divergence testing by incorporating background environments, analogous-environment filtering, variable transformation, low-variation filtering, and occurrence thinning.
- Can handle hundreds of correlated environmental variables.
- Identifies environmental variables contributing most to niche separation.
- Visualizes DAPC-based niche divergence along a single discriminant axis, making results easy to interpret.

## Development status
NicheDiv is under active development. The methodological framework is described in an associated preprint (https://doi.org/10.64898/2026.06.19.733388), and the manuscript is currently in review.  

For bug reports, feedback, or questions, please contact me: daniel.schoenberger@uky.edu.


## Installation

Install the R package from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("Daniel-1232/NicheDiv")
```

Load the package:

```r
library(NicheDiv)
```


## Input data

The approach only requires a data frame with occurrence records:

* one row per occurrence record;
* unique row names or an ID column;
* longitude and latitude columns;
* one grouping column with two or more groups to compare;

Example:

```r
head(occurrence_data[, c("ID", "Longitude", "Latitude", "Species")])
```

```text
     ID Longitude Latitude Species
1 ID_001  -121.50   38.40   Sp1
2 ID_002  -121.75   38.55   Sp1
3 ID_003  -117.25   34.15   Sp2
4 ID_004  -117.10   34.40   Sp2
```

## Recommended workflow

The typical NicheDiv workflow has seven major steps.
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
N_permutations <- 1000
CV_threshold <- 0.01
base_colors <- c("#A331A3", "#6CB3A5")

exclude_cols <- c("ID", "Locality", "CollectionDate") #columns to ignore throughout
set.seed(1)
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
                                     remove.hydrolakes.background = TRUE,
                                     csv.occurrence.out.file = csv_occurrence_out_file,
                                     csv.background.out.file = csv_background_out_file,
                                     output.dir = results_dir,
                                     intermediate.files.dir = intermediate_files_dir_name,
                                     CRS.occurrences = CRS_all,
                                     env.datasets = c("elevation", "ClimateNA", "EVI", "terrain",
                                                      "ENVIREM", "footprint", "landcover", "soil",
                                                      "forest_height", "atmosphere", "nightlight",
                                                      "burned_area", "snow_water_equivalent",
                                                      "daylength", "soil_moisture"))
```

In the example above, all implemented environmental datasets are used. This can take several hours because the raster layers need to be downloaded.
Some datasets are only available for North America ("ClimateNA", "daylength", "snow_water_equivalent")

Optional custom rasters can also be supplied by the user as one or more GeoTIFF files:

```r
custom_raster_path <- file.path(base_dir, "Data/custom_environmental_layers.tif")
custom_raster_variable_names <- names(terra::rast(custom_raster_path))

NicheDiv::extract.env.and.background(occurrence.data = occurrence_data,
                                     longitude.col = Longitude_col,
                                     latitude.col = Latitude_col,
                                     generate.background.data = TRUE,
                                     N.background.points = N_background_points,
                                     buffer.km = buffer_km,
                                     remove.hydrolakes.background = TRUE,
                                     csv.occurrence.out.file = csv_occurrence_out_file,
                                     csv.background.out.file = csv_background_out_file,
                                     output.dir = results_dir,
                                     intermediate.files.dir = intermediate_files_dir_name,
                                     CRS.occurrences = CRS_all,
                                     env.datasets = c("elevation", "ClimateNA"),
                                     custom.env.rasters = custom_raster_path,
                                     custom.env.rasters.variable.names = custom_raster_variable_names)
```

## 3. Import and prepare extracted data

```r
#### Import and prepare extracted data #########################################

## Import extracted occurrence and background data
Env_data_occurrences <- read.csv(file.path(results_dir, csv_occurrence_out_file))
Env_data_background <- read.csv(file.path(results_dir, csv_background_out_file), check.names = FALSE)


## Remove metadata columns not used as environmental predictors
Env_data_occurrences <- Env_data_occurrences[, setdiff(colnames(Env_data_occurrences), exclude_cols), drop = FALSE]
Env_data_background <- Env_data_background[, setdiff(colnames(Env_data_background), exclude_cols), drop = FALSE]


## Convert integer columns to numeric
Env_data_occurrences <- NicheDiv::convert.integer.to.numeric(Env_data_occurrences)
Env_data_background <- NicheDiv::convert.integer.to.numeric(Env_data_background)


## Keep the two groups of interest
Env_data_occurrences <- Env_data_occurrences[Env_data_occurrences[[Species_col]] %in% c(Sp1_name, Sp2_name), , drop = FALSE]


## Rename groups for plotting
group_name_map <- setNames(c(Sp1_label, Sp2_label), c(Sp1_name, Sp2_name))
Env_data_occurrences[[Species_col]] <- group_name_map[as.character(Env_data_occurrences[[Species_col]])]
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
Sp1_background_data <- NicheDiv::sample.down(Sp1_background_data, N.rows = 10000)
Sp2_background_data <- NicheDiv::sample.down(Sp2_background_data, N.rows = 10000)
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
Sp1_Sp2_occurrence_thinned <- rbind(Sp1_occurrence_thinned, Sp2_occurrence_thinned)
Sp1_background_data[[Species_col]] <- Sp1_label
Sp2_background_data[[Species_col]] <- Sp2_label
Sp1_Sp2_background_data <- rbind(Sp1_background_data, Sp2_background_data)


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

Sp1_Sp2_occurrence_filtered <- rbind(Sp1_occurrence_filtered, Sp2_occurrence_filtered)
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

## 7. Run DAPC

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


* `Schoener_D (D)`: niche overlap between the two groups along the discriminant axis. Values range from 0 to 1, where 1 indicates complete overlap and 0 indicates no overlap. Lower values therefore indicate stronger niche differentiation.
* `Niche_dissimilarity (NDS)`: density-based niche divergence along the discriminant axis. This metric captures how strongly the two groups differ in the distribution of their occurrence densities, even when their total occupied ranges overlap.
* `Niche_breadth_exclusivity (NE)`: range-based niche divergence along the discriminant axis. This metric captures how much of each group’s occupied environmental range is exclusive rather than shared with the other group.
* `Niche_divergence_magnitude (ND)`: combined divergence magnitude in the niche divergence plane. This summarizes the joint strength of density-based divergence and range-based exclusivity in a single value.
* `Niche_divergence_angle (θ)`: relative contribution of density-based versus range-based divergence. Angles closer to the density-based axis indicate that divergence is mainly driven by differences in occurrence density, whereas angles closer to the range-based axis indicate that divergence is mainly driven by exclusive environmental ranges.

The most important summary metrics are `D` and `ND`. Stronger niche divergence is indicated by lower `D` values and higher `ND` values. As a general rule of thumb: `D` values around or below 0.4 and `ND` values around or above 0.9 indicate strong divergence in the current framework.

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

background_data_combined <- rbind(Sp1_background_data, Sp2_background_data)

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

* Use biologically meaningful groups. NicheDiv tests divergence between predefined groups and does not infer species limits by itself.

Use a biologically justified accessible area. Background data should represent the environmental conditions plausibly available to each group.

Balance occurrence sample sizes before DAPC. Unequal sample sizes can affect discrimination and interpretation.

Inspect spatial thinning diagnostics. Thinning reduces spatial clustering but may not remove all spatial autocorrelation.

Use analogous-environment filtering when comparing groups from different accessible areas. This helps reduce bias from environmental conditions available to only one group.

* Do not interpret variable contributions as proof of causation. They identify predictors associated with niche divergence and should be evaluated with biological knowledge, natural history, and independent evidence.


## Citation

Please cite the NicheDiv framework as follows:

Schönberger, D., MacDonald, Z. G., Schmidt, B. C., & Dupuis, J. R. NicheDiv: A DAPC framework to quantify niche divergence across highly multivariate environmental space. bioRxiv. https://doi.org/10.64898/2026.06.19.733388 



## License

NicheDiv is released under the MIT License. See the `LICENSE` file for details.
