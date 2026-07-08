# *NicheDiv* *R* package

*NicheDiv* is an *R* package for testing pairwise niche divergence across highly multivariate environmental space.

This is done by adapting discriminant analysis of principal components (DAPC) to environmental niche data. DAPC was originally developed for population genetics (Jombart et al. 2010) but is well suited for numerous correlated environmental variables. Environmental variables are first transformed into principal components (PCs) to reduce dimensionality and collinearity. Discriminant analysis (Fisher 1936, Lachenbruch & Goldstein 1979) is then used to identify the axis that best separates the two groups. Our method summarizes niche divergence with easily interpretable metrics and density plots.

The motivation behind our approach is that ecological niches are highly multidimensional (Hutchinson, 1957) and are rarely captured completely by the commonly used annual climate variables alone (Elith & Leathwick 2009; Kearney & Porter 2009; Soberón, 2007), such as WorldClim’s BIO1–BIO19 variables (Hijmans et al. 2005). Seasonal and monthly variables are often required to capture phenology, resource availability, physiological stress, and other time-dependent ecological processes that may be obscured by annual averages (Prajzlerová et al. 2025; Zimmermann et al. 2009). *NicheDiv* tackles this problem in two ways: first, by automatically extracting environmental values from a broad set of implemented GIS layers covering both abiotic and biotic environmental dimensions; and second, by making it possible to test niche divergence across this high-dimensional and correlated environmental space using our DAPC-based framework.

## Main advantages of the approach

- Requires only occurrence data as input.
- Automatically extracts environmental values for occurrence records and background points from implemented and user-supplied environmental GIS layers. Implemented GIS layers cover monthly to seasonal climate, topography, phenology, hydrology, vegetation, soil, land cover, and anthropogenic variables (most at global extent).
- Implements a preprocessing pipeline that reduces common biases: delimiting accessible background space, spatially thinning occurrences, balancing sample sizes, filtering low-information variables, and screening predictors for between-group environmental analogy (e.g., Dormann et al. 2013, Soberón 2007, Brown & Carnaval 2019).
- Can handle hundreds of correlated environmental variables.
- Identifies environmental variables contributing most to niche separation.
- Visualizes multivariate niche divergence along a single discriminant axis, making results easy to interpret.
- Can distinguish different forms of niche divergence (weighted, nested, soft, and hard niche divergence)
- Compared with alternative divergence tests, *NicheDiv* generally retains more variation, and scales more consistently with increasing divergence.

## Development status
The framework is described in a preprint (https://doi.org/10.64898/2026.06.19.733388) and the manuscript is currently in review.  

Current *R* package version: 0.1.0

For bug reports, feedback, or questions, please contact me: daniel.schoenberger@uky.edu.


# Tutorial

## Installation

Install and load the *NicheDiv* *R* package:

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("Daniel-1232/NicheDiv")

library(NicheDiv)
packageVersion("NicheDiv")
```


## Input data

The approach only requires a dataframe with occurrence records:

* one row per occurrence record
* unique row names
* longitude and latitude columns
* one grouping column with two (or more) groups to compare

Example input:

```r
head(occurrence_data[, c("ID", "Longitude", "Latitude", "Species")])
rownames(occurrence_data) <- occurrence_data$ID
```

```text
     ID Longitude Latitude Species
1 ID_001  -121.50   38.40   Sp1
2 ID_002  -121.75   38.55   Sp1
3 ID_003  -117.25   34.15   Sp2
4 ID_004  -117.10   34.40   Sp2
```

The dataframe can include other columns as long as they are specified under `exclude_cols` (see below).
Groups can be species, populations, lineages or any other predefined groupings or clusters.
The dataframe can also include multiple species if you want to perform multiple pairwise comparisons (see section "How to include multiple pairwise comparisons" at the end).


## Recommended workflow

The NicheDiv workflow has several major steps. The code below describes the full workflow using recommended default parameters throughout. Parameters that may require tuning are discussed explicitly.

Below is a schematic overview of the niche divergence framework, using two theoretical taxon pairs and three environmental layers as example (figure 1 from Schönberger et al. preprint):

![NicheDiv workflow](man/figures/README-schoenberger-etal-figure-1.png)


## Set working environment and input parameters

Before starting, we need to define all directories, file names and parameters.

```r
#### Set working environment and input #########################################

## Set directories
base_dir <- "path/to/project"
results_dir <- file.path(base_dir, "Results")
figure_dir <- file.path(results_dir, "Figure_files")
intermediate_files_dir_name <- "Intermediate_files"

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


## Set input occurrence file
occurrence_data_file <- file.path(base_dir, "Data/occurrences.csv")


## Set output occurrence files with environmental data
csv_occurrence_out_file <- "Occurrences_env.csv"
csv_background_out_file <- "Background_env.csv"


## Set parameters and column names
Sp1_name <- "Group_1"
Sp2_name <- "Group_2"
Sp1_label <- "Group 1"
Sp2_label <- "Group 2"

Species_col <- "Species"
ID_col <- "ID"

Longitude_col <- "Longitude"
Latitude_col <- "Latitude"
CRS_all <- "EPSG:4326"

buffer_km <- 5
base_colors <- c("#A331A3", "#6CB3A5")
exclude_cols <- c(ID_col, "Locality", "CollectionDate")
```

Use `Sp1_name` and `Sp2_name` for the group names exactly as they appear in the grouping column of your input dataframe (e.g., `"Hemileuca_nevadensis"`), and use `Sp1_label` and `Sp2_label` for the labels displayed in plots (e.g., `"H. nevadensis"`).

`buffer_km` should be chosen to reflect the estimated approximate dispersal distance of the species group.

Use `exclude_cols` to list columns that should be excluded from environmental predictor variables throughout the workflow, such as IDs, locality names, or collection dates.

`CRS_all` defines the coordinate reference system of the occurrence coordinates. Use `"EPSG:4326"` when your longitude and latitude columns are in decimal degrees, which is the most common format for occurrence data. If your coordinates are already projected, provide the corresponding projected CRS instead.


## 1. Extract environmental data and generate background points

We start by extracting environmental values for occurrence records and generating background points within the accessible area.
This step usually takes the most time because the environmental layers need to be downloaded. Fortunately, our approach uses minimal GIS layer processing/projection, saving hours of time and a lot of memory.
In the example below, all implemented environmental datasets are used, which can take several hours. Using all datasets is typically a good approach to describe the niche as comprehensively as possible, but your study system may require excluding datasets that are biologically less relevant.
Furthermore, some datasets are only available for North America (namely, `"ClimateNA"`, `"daylength"`, and `"snow_water_equivalent"`). If your study system is outside North America, remove these datasets from `env.datasets`.

For terrestrial taxa, we recommend setting `remove.hydrolakes.background = TRUE`, which prevents background points from being sampled from lakes and other major inland water bodies.

```r
#### Extract environmental data ################################################

## Import occurrences
occurrence_data <- read.csv(occurrence_data_file, check.names = FALSE)
rownames(occurrence_data) <- occurrence_data[[ID_col]]
if (anyDuplicated(rownames(occurrence_data)) > 0) stop("Occurrence IDs must be unique")

## Extract environmental data and background points
extract.env.and.background(occurrence.data = occurrence_data,
                           longitude.col = Longitude_col,
                           latitude.col = Latitude_col,
                           generate.background.data = TRUE,
                           N.background.points = 300000,
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

Optional custom rasters can also be supplied by the user as one or more GeoTIFF files:

```r
custom_raster_path <- file.path(base_dir, "Data/custom_environmental_layers.tif")
custom_raster_variable_names <- names(terra::rast(custom_raster_path))

NicheDiv::extract.env.and.background(occurrence.data = occurrence_data,
                                     longitude.col = Longitude_col,
                                     latitude.col = Latitude_col,
                                     generate.background.data = TRUE,
                                     N.background.points = 300000,
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

## 2. Import and prepare extracted data

Next, we import and process the extracted environmental data. In this section, no changes in parameters are needed.

```r
#### Import and prepare extracted data #########################################

## Import extracted occurrence and background data
Env_data_occurrences <- read.csv(file.path(results_dir, csv_occurrence_out_file))
dim(Env_data_occurrences)
Env_data_background <- read.csv(file.path(results_dir, csv_background_out_file), check.names = FALSE)
dim(Env_data_background)


## Remove metadata columns not used as environmental predictors
Env_data_occurrences <- Env_data_occurrences[, setdiff(colnames(Env_data_occurrences), exclude_cols), drop = FALSE]
Env_data_background <- Env_data_background[, setdiff(colnames(Env_data_background), exclude_cols), drop = FALSE]


## Convert integer columns to numeric
Env_data_occurrences <- convert.integer.to.numeric(Env_data_occurrences)
Env_data_background <- convert.integer.to.numeric(Env_data_background)


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

## 3. Crop and downsample background data

The next step crops the shared background to a buffered convex hull around each group’s occurrence records and then downsamples each background to the same target size. Only the `buffer.method` argument needs to be considered in this section. This argument defines how the accessible area is buffered around each group’s occurrence records, with larger or more inclusive buffers retaining more background environments and smaller or stricter buffers focusing the comparison on environments closer to the observed occurrences.
Available background geometries are `"hull"`, `"points"`, `"alpha"`, and `"bbox"`. Below we use the convex hull which is usually a robust default. Point buffers or alpha hulls may be useful for fragmented or spatially complex distributions.

```r
#### Prepare background data ###################################################

## Crop background to each group-specific accessible area
Sp1_background_data <- crop.background.buffered(occurrence.data = Sp1_occurrence_data,
                                                background.data = Env_data_background,
                                                latitude.col = Latitude_col,
                                                longitude.col = Longitude_col,
                                                CRS = CRS_all,
                                                buffer.method = "hull",
                                                buffer.dist.meters = buffer_km * 1000)

Sp2_background_data <- crop.background.buffered(occurrence.data = Sp2_occurrence_data,
                                                background.data = Env_data_background,
                                                latitude.col = Latitude_col,
                                                longitude.col = Longitude_col,
                                                CRS = CRS_all,
                                                buffer.method = "hull",
                                                buffer.dist.meters = buffer_km * 1000)


## Downsample background data
Sp1_background_data <- sample.down(Sp1_background_data, N.rows = 10000)
Sp2_background_data <- sample.down(Sp2_background_data, N.rows = 10000)
```

## 4. Spatially thin and balance occurrence records

To reduce spatial autocorrelation, we thin our occurrence records. A thinning distance (`thinning.dist.km`) of one kilometer is usually an appropriate value, as set below. If hundreds of occurrence records remain after thinning, you can consider increasing the thinning distance threshold to reduce spatial autocorrelation. 
We also downsample both groups to the same number of occurrences (to avoid bias in the discriminant analysis caused by unequal sample sizes).

```r
#### Spatial thinning and sample-size balancing ################################

## Thin occurrence records
Sp1_occurrence_thinned <- thin.occurrence(Sp1_occurrence_data,
                                          latitude.col = Latitude_col,
                                          longitude.col = Longitude_col,
                                          thinning.dist.km = 1)
Sp2_occurrence_thinned <- thin.occurrence(Sp2_occurrence_data,
                                          latitude.col = Latitude_col,
                                          longitude.col = Longitude_col,
                                          thinning.dist.km = 1)


## Downsample to equal sample size
n_min_occurrence_thinned <- min(nrow(Sp1_occurrence_thinned), nrow(Sp2_occurrence_thinned))

Sp1_occurrence_thinned <- sample.down(Sp1_occurrence_thinned,
                                      N.rows = n_min_occurrence_thinned)
Sp2_occurrence_thinned <- sample.down(Sp2_occurrence_thinned,
                                      N.rows = n_min_occurrence_thinned)
```

## 5. Transform and filter environmental variables

Next, we transform skewed variables and remove variables with low variation. This reduces the influence of extreme values and removes predictors that contribute little information to group separation. In this section, no changes in parameters are needed.

```r
#### Transform skewed environmental variables ##################################

## Combine occurrence and background datasets
Sp1_Sp2_occurrence_thinned <- rbind(Sp1_occurrence_thinned, Sp2_occurrence_thinned)
Sp1_background_data[[Species_col]] <- Sp1_label
Sp2_background_data[[Species_col]] <- Sp2_label
Sp1_Sp2_background_data <- rbind(Sp1_background_data, Sp2_background_data)


## Transform skewed variables
transformation_results <- transform.skewed.variables(data.frame = Sp1_Sp2_occurrence_thinned,
                                                     exclude.cols = c(Latitude_col, Longitude_col, Species_col, ID_col),
                                                     background.dataframe = Sp1_Sp2_background_data)
Sp1_Sp2_occurrence_transformed <- transformation_results$transformed
Sp1_Sp2_background_transformed <- transformation_results$background.transformed


## Split transformed data by group
Sp1_occurrence_transformed <- Sp1_Sp2_occurrence_transformed[Sp1_Sp2_occurrence_transformed[[Species_col]] == Sp1_label, , drop = FALSE]
Sp2_occurrence_transformed <- Sp1_Sp2_occurrence_transformed[Sp1_Sp2_occurrence_transformed[[Species_col]] == Sp2_label, , drop = FALSE]

Sp1_background_transformed <- Sp1_Sp2_background_transformed[Sp1_Sp2_background_transformed[[Species_col]] == Sp1_label, , drop = FALSE]
Sp2_background_transformed <- Sp1_Sp2_background_transformed[Sp1_Sp2_background_transformed[[Species_col]] == Sp2_label, , drop = FALSE]
```

```r
#### Remove low-information variables ##########################################
CV_removal_results <- remove.low.CV.vars(Sp1.occurrence.data = Sp1_occurrence_transformed,
                                         Sp2.occurrence.data = Sp2_occurrence_transformed,
                                         Sp1.background.data = Sp1_background_transformed,
                                         Sp2.background.data = Sp2_background_transformed,
                                         exclude.cols = c(Latitude_col, Longitude_col, Species_col),
                                         CV.threshold = 0.01)

Sp1_occurrence_filtered <- CV_removal_results$occurrence_Sp1
Sp2_occurrence_filtered <- CV_removal_results$occurrence_Sp2
Sp1_background_filtered <- CV_removal_results$background.Sp1
Sp2_background_filtered <- CV_removal_results$background.Sp2

Sp1_Sp2_occurrence_filtered <- rbind(Sp1_occurrence_filtered, Sp2_occurrence_filtered)
```

## 6. Filter to analogous environmental variables

This step reduces bias from non-analogous environments by filtering out variables that show insufficient overlap between the accessible background environments of the two groups. Environmental analogy is assessed using univariate kernel-density overlap and bivariate histogram overlap. Variables with little overlap between background spaces are removed so that the DAPC comparison is restricted to environmental conditions that are comparably available to both groups. This helps prevent apparent niche divergence from being driven by environmental gradients that one group could access but the other group could not.

```r
#### Filter to analogous environmental variables ###############################
Sp1_Sp2_analogous <- filter.analogous.variables(Sp1.Sp2.occurrence.data = Sp1_Sp2_occurrence_filtered,
                                                Sp1.background.data = Sp1_background_filtered,
                                                Sp2.background.data = Sp2_background_filtered,
                                                exclude.cols = c(Latitude_col, Longitude_col, Species_col),
                                                CV.threshold = 0.01,
                                                overlap.threshold = 0.7)
```


## 7. Run DAPC

Finally, we run the main niche divergence analysis by applying DAPC to the filtered environmental data. The function first performs a PCA to reduce dimensionality and collinearity, and then uses discriminant analysis to identify the axis that best separates the two groups in multivariate environmental space.

Cross-validation is used to select the number of PCs retained for DAPC. This helps retain enough environmental variation to separate the groups while avoiding overfitting caused by retaining too many PCs.

To assess significance, a permutation test is also performed that compares the observed DAPC assignment accuracy to a null distribution generated by randomly permuting group labels emulating a single shared niche (k = 1). A significant result indicates that group separation along the discriminant axis is stronger than expected under random group membership.
Based on simulation testing, the permutation test is highly sensitive and can become significant already at low to moderate levels of niche divergence. Therefore, statistical significance should be interpreted together with the divergence metrics and discriminant density plots

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
DAPC_results <- run.DAPC.crossval.permutation(data.input = Sp1_Sp2_analogous,
                                              species.col = Species_col,
                                              exclude.cols = c(Latitude_col, Longitude_col),
                                              N.permutations = 1000,
                                              N.crossval.replicates = 300)
```

Based on the DAPC results, we can calculate the following five niche divergence metrics:

* `Schoener_D (D)`: niche overlap between the two groups along the discriminant axis. Values range from 0 to 1, where 1 indicates complete overlap and 0 indicates no overlap (Schoener 1968).
* `Niche_dissimilarity (NDS)`: density-based niche divergence along the discriminant axis. Values range from 0 to 1, where 0 indicates identical occurrence-density distributions and 1 indicates completely non-overlapping densities (Ascanio et al. 2024).
* `Niche_breadth_exclusivity (NE)`: range-based niche divergence along the discriminant axis. Values range from 0 to 1, where 0 indicates completely shared occupied ranges and 1 indicates completely exclusive occupied ranges (Ascanio et al. 2024).
* `Niche_divergence_magnitude (ND)`: combined divergence magnitude in the niche divergence plane. Values range from 0 to 1.41, where 0 indicates no divergence and 1.41 indicates maximum combined density-based and range-based divergence (Ascanio et al. 2024).
* `Niche_divergence_angle (θ)`: relative contribution of density-based versus range-based divergence. Values range from 0° to 90°, where values near 0° indicate divergence mainly driven by range exclusivity, values near 90° indicate divergence mainly driven by density differences within shared space, and intermediate values indicate mixed contributions (Ascanio et al. 2024).

The most important summary metrics are `D` and `ND`. Stronger niche divergence is indicated by lower `D` values and higher `ND` values. As a general rule of thumb: `D` values below 0.4 and `ND` values above 0.9 indicate strong divergence in the current framework.

```r
#### Calculate niche divergence metrics ########################################
Niche_divergence_metrics <- calc.niche.divergence.metrics(DAPC_results,
                                                          group.assignment = Sp1_Sp2_species_assignment)
```

Optionally, we can calculate background-corrected metrics (following Brown and Carnaval 2019) by up-weighting rare and down-weighting common environments along the discriminant axis to account for unequal environmental availability.

```r
Niche_divergence_metrics_weighted <- calc.niche.divergence.metrics(DAPC_results,
                                                                   weight.background = TRUE,
                                                                   Sp1.background.data = Sp1_background_filtered,
                                                                   Sp2.background.data = Sp2_background_filtered,
                                                                   group.assignment = Sp1_Sp2_species_assignment)
```

## 8. Plot results
In general, all plot functions include built-in saving options. Set save = TRUE to export figures directly as SVG, PNG, or JPEG files using the type argument. Figure dimensions can be adjusted with width and height. Set save = FALSE if you do not want to save the figures. The overwrite argument controls whether existing plot files are overwritten. Many plot functions also include additional arguments for adjusting font sizes and other plotting parameters.

We start by plotting the discriminant-axis density distributions, followed visualizing the permutation null distribution of classification accuracy (observed value shown as red line)

```r
#### Plot DAPC niche divergence ################################################
plot.DAPC.niche.divergence(DAPC_results,
                           group.colors = Sp1_Sp2_species_colors,
                           save = TRUE,
                           overwrite = TRUE,
                           type = "svg",
                           output.dir = figure_dir,
                           filename = "DAPC_niche_divergence",
                           width = 16,
                           height = 12)


#### Plot permutation test #####################################################
plot.DAPC.permutation(DAPC_results,
                      save = TRUE,
                      overwrite = TRUE,
                      type = "svg",
                      output.dir = figure_dir,
                      filename = "DAPC_permutation_test",
                      width = 16,
                      height = 9)
```

Here an example output from the two functions above showing strong multivariate niche divergence in this *Hemileuca maia* buck moth group (figure 4 from Schönberger et al. preprint):

![NicheDiv example result](man/figures/README-schoenberger-etal-figure-4.png)

Plot: environmental variable contributions to the discriminant axis

These values show which original environmental variables contribute most to the DAPC separation between the two groups. Contributions are calculated by back-transforming the discriminant axis from retained PCs to the original environmental variables. Higher values indicate variables that contribute more strongly to group separation, but they should not be interpreted as independent causal effects because correlated predictors can share the same signal.


```r
#### Plot variable contributions ##############################################
DAPC_results_short_names <- DAPC_results
DAPC_results_short_names$dapc_results$var.contr <- map.env.variable.names(DAPC_results_short_names$dapc_results$var.contr, "short")
DAPC_results_short_names$dapc_results$var.load <- map.env.variable.names(DAPC_results_short_names$dapc_results$var.load, "short")

DAPC_var_contr <- plot.DAPC.var.contributions(DAPC_results_short_names,
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

Plot: raw distributions of the top contributing predictors

```r
#### Plot top predictors #######################################################
Sp1_Sp2_analogous_short_names <- map.env.variable.names(Sp1_Sp2_analogous, "short")

plot.top.DAPC.predictors(dapc.results = DAPC_results_short_names,
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
Here is an example output from the two variable-contribution plotting functions above (figure 5 from Schönberger et al. preprint). The figure summarizes which environmental variables contribute most strongly to the DAPC-based separation of the two taxa in multivariate niche space. Panel A shows the relative contribution of each predictor to the discriminant axis and indicates which species has higher values for each variable. Panel B shows the distributions of the strongest contributing predictors, illustrating how univariate differences in these variables drive the estimated niche divergence. 

![NicheDiv example result](man/figures/README-schoenberger-etal-figure-5.png)


Plot: occurrences and background points

This map is useful for checking the geographic distribution of the two groups, the sampled background environments, and whether the accessible areas are biologically reasonable. 
Some map elements may need to be adjusted depending on the study area, map extent, and figure size. In particular, the arguments `north.arrow.length`, `north.arrow.N.position`, `north.arrow.position`, `scale.position`, `longitude.buffer.range`, `latitude.buffer.range`, and `north.arrow.lwd` may need manual tuning to avoid overlap with points or map boundaries.

```r
#### Plot occurrence and background map ########################################
background_labels <- factor(c(rep(levels(Sp1_Sp2_species_assignment)[1], nrow(Sp1_background_data)),
                              rep(levels(Sp1_Sp2_species_assignment)[2], nrow(Sp2_background_data))),
                            levels = levels(Sp1_Sp2_species_assignment))

background_data_combined <- rbind(Sp1_background_data, Sp2_background_data)

plot.occurrences.map(coordinates = Sp1_Sp2_analogous,
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
Here an example map (figure 3 from Schönberger et al. preprint): the large points represent occurrence records and the small points background records

![NicheDiv example result](man/figures/README-schoenberger-etal-figure-3.png)


## DAPC in full-environment

We can also run the DAPC test using the full-environment (without non-analogous filtering). 
This helps to evaluate how much non-analogous environmental space affects the final result.
This analysis is especially useful if many variables are removed during analogous-variable filtering, because some excluded variables may be biologically relevant and potentially contribute to divergence. However, in the full-environment analysis, we cannot determine whether these variables reflect niche divergence within shared accessible environmental space or differences in environmental availability between taxa.

```r
#### Optional DAPC without analogous-variable filtering ########################
Sp1_Sp2_species_assignment_no_analogy <- factor(Sp1_Sp2_occurrence_filtered[[Species_col]])

Sp1_Sp2_species_colors_no_analogy <- setNames(base_colors[seq_along(levels(Sp1_Sp2_species_assignment_no_analogy))],
                                              levels(Sp1_Sp2_species_assignment_no_analogy))

Sp1_Sp2_species_assignment_no_analogy <- factor(Sp1_Sp2_species_assignment_no_analogy,
                                                levels = names(Sp1_Sp2_species_colors_no_analogy))

DAPC_results_no_analogy <- run.DAPC.crossval.permutation(data.input = Sp1_Sp2_occurrence_filtered,
                                                         species.col = Species_col,
                                                         exclude.cols = c(Latitude_col, Longitude_col),
                                                         N.permutations = 1000,
                                                         N.crossval.replicates = 300)

Niche_divergence_metrics_no_analogy <- calc.niche.divergence.metrics(DAPC_results_no_analogy,
                                                                     group.assignment = Sp1_Sp2_species_assignment_no_analogy)
```

## Optional: Brown and Carnaval-style analogous trimming

In addition to variable-level analogy filtering (`filter.analogous.variables()`), *NicheDiv* includes `trim.to.analogous.environments()` to remove occurrence records from non-analogous environmental conditions following the logic of Brown and Carnaval-style environmental analogy correction (Brown & Carnaval 2019).

This occurrence record-based analogy filtering (`trim.to.analogous.environments()`) may be preferred over variable-based analogy filtering if the latter removes many variables and therefore limits inference about the environmental variables contributing most to niche separation. However, if occurrence record-based trimming removes many records, DAPC inference may become less stable because of reduced sample size. In that case, variable-based analogy filtering may be preferred.


```r
#### Optional Brown and Carnaval-style correction ##############################
Sp1_Sp2_analogous_trimmed <- trim.to.analogous.environments(Sp1.occurrence.data = Sp1_occurrence_filtered,
                                                            Sp2.occurrence.data = Sp2_occurrence_filtered,
                                                            Sp1.background.data = Sp1_background_filtered,
                                                            Sp2.background.data = Sp2_background_filtered,
                                                            exclude.cols = c(Latitude_col, Longitude_col, Species_col),
                                                            keep.occurrence.cols = c(Latitude_col, Longitude_col, Species_col))
```

The trimmed dataset can then be passed to `run.DAPC.crossval.permutation()` using the same DAPC workflow shown above.


## How to include multiple pairwise comparisons
If you have multiple taxa (e.g., all members of a species group), you can compare them by running *NicheDiv* in a pairwise fashion:
1) Use input dataframe with coordinates for all taxa of interest.
2) Extract environmental data and generate background points once for all taxa together (Step 1 in workflow: `extract.env.and.background()` function). We recommend increasing `N.background.points` to ensure enough background points for all comparisons (1 million worked well in our species group for North America). 
3) For each pairwise comparison, first set `Sp1_name`, `Sp2_name`, `Sp1_label`, `Sp2_label`, `base_colors`, and `filename` in each plotting call (or change `figure_dir` for each), and then run steps 2-7 of the workflow for each pair.


## Further recommendations

* We recommend first running the DAPC niche divergence test using only analogous environmental variables. Strong and significant divergence in this analysis suggests that the groups differ within shared accessible environmental space (Brown & Carnaval, 2019). If no analogous variables remain after filtering, this also provides evidence that the groups occupy strongly different accessible environments.

* Interpret permutation-test significance together with the divergence metrics and discriminant density plots. The permutation test can be highly sensitive and may become significant at low to moderate levels of divergence. Conversely, a non-significant result can reflect either true niche similarity or limited statistical power given the available environmental conditions and sample sizes.

* Interpret environmental variable contributions as hypothesis-generating rather than causal. Variable contributions identify predictors that contribute to multivariate separation along the discriminant axis, but they do not prove that these variables independently drive divergence. High contributions may reflect correlated sets of predictors, while low contributions do not rule out biological importance if the signal is absorbed by correlated variables.

* Ecological divergence should not by itself be interpreted as evidence of ecological speciation. NicheDiv tests realized niche divergence under current environmental and distributional conditions. Inferring the timing, mechanism, or evolutionary cause of divergence requires additional evidence, such as natural-history data, experiments, demographic analyses, or phylogeographic analyses.

* NicheDiv currently only supports continuous environmental variables. Because DAPC is widely used with biallelic genetic markers (Jombart et al. 2010, Miller et al. 2020), the framework could potentially be extended to binary or categorical ecological predictors in the future. If you want to include binary or categorical data (e.g., host presence/absence, habitat classes, symbionts, or pollinator types), running a SOM (self-organizing map) model may be useful (Pyron et al. 2023; see https://github.com/rpyron/delim-SOM).


## Main functions

| Function                           | Description                                                                     |
| ---------------------------------- | ------------------------------------------------------------------------------- |
| `extract.env.and.background()`     | Extract environmental variables and generate background data                    |
| `convert.integer.to.numeric()`     | Convert integer columns to numeric                                              |
| `crop.background.buffered()`       | Crop background points to buffered accessible areas                             |
| `sample.down()`                    | Downsample occurrence or background records                                     |
| `thin.occurrence()`                | Spatially thin occurrence records and evaluate residual spatial autocorrelation |
| `transform.skewed.variables()`     | Transform skewed environmental variables                                        |
| `remove.low.CV.vars()`             | Remove variables with low coefficient of variation                              |
| `filter.analogous.variables()`     | Retain predictors with analogous background distributions                       |
| `trim.to.analogous.environments()` | Remove occurrence records from non-analogous environmental conditions           |
| `run.DAPC.crossval.permutation()`  | Run cross-validated DAPC and permutation testing                                |
| `calc.niche.divergence.metrics()`  | Calculate Schoener’s D and niche divergence plane metrics                       |
| `plot.DAPC.niche.divergence()`     | Plot density distributions along the DAPC discriminant axis                     |
| `plot.DAPC.permutation()`          | Plot the permutation null distribution                                          |
| `plot.DAPC.var.contributions()`    | Plot variable contributions to the discriminant axis                            |
| `plot.top.DAPC.predictors()`       | Plot raw distributions of top contributing predictors                           |
| `plot.occurrences.map()`           | Plot occurrence and background records on a map                                 |
| `map.env.variable.names()`         | Convert environmental variable names to shorter or more readable labels         |


## References

* Ascanio, A., Bracken, J. T., Stevens, M. H. H., & Jezkova, T. (2024). New theoretical and analytical framework for quantifying and classifying ecological niche differentiation. *Ecological Monographs, 94*(4). https://doi.org/10.1002/ecm.1622

* Brown, J. L., & Carnaval, A. C. (2019). A tale of two niches: Methods, concepts, and evolution. *Frontiers of Biogeography, 11*(4). https://doi.org/10.21425/F5FBG44158

* Dormann, C. F., Elith, J., Bacher, S., Buchmann, C., Carl, G., Carré, G., Marquéz, J. R. G., Gruber, B., Lafourcade, B., Leitão, P. J., Münkemüller, T., McClean, C., Osborne, P. E., Reineking, B., Schröder, B., Skidmore, A. K., Zurell, D., & Lautenbach, S. (2013). Collinearity: A review of methods to deal with it and a simulation study evaluating their performance. *Ecography, 36*(1), 27–46. https://doi.org/10.1111/j.1600-0587.2012.07348.x

* Dormann, C. F., McPherson, J. M., Araújo, M. B., Bivand, R., Bolliger, J., Carl, G., Davies, R. G., Hirzel, A., Jetz, W., Kissling, W. D., Kühn, I., Ohlemüller, R., Peres-Neto, P. R., Reineking, B., Schröder, B., Schurr, F. M., & Wilson, R. (2007). Methods to account for spatial autocorrelation in the analysis of species distributional data: A review. *Ecography, 30*(5), 609–628. https://doi.org/10.1111/j.2007.0906-7590.05171.x

* Elith, J., & Leathwick, J. R. (2009). Species distribution models: Ecological explanation and prediction across space and time. *Annual Review of Ecology, Evolution, and Systematics, 40*, 677–697. https://doi.org/10.1146/annurev.ecolsys.110308.120159

* Fisher, R. A. (1936). The use of multiple measurements in taxonomic problems. *Annals of Eugenics, 7*(2), 179–188. https://doi.org/10.1111/j.1469-1809.1936.tb02137.x

* Hijmans, R. J., Cameron, S. E., Parra, J. L., Jones, P. G., & Jarvis, A. (2005). Very high resolution interpolated climate surfaces for global land areas. *International Journal of Climatology, 25*(15), 1965–1978. https://doi.org/10.1002/joc.1276

* Hutchinson, G. E. (1957). Concluding remarks. *Cold Spring Harbor Symposia on Quantitative Biology, 22*, 415–427. https://doi.org/10.1101/SQB.1957.022.01.039

* Jombart, T., Devillard, S., & Balloux, F. (2010). Discriminant analysis of principal components: A new method for the analysis of genetically structured populations. *BMC Genetics, 11*, 94. https://doi.org/10.1186/1471-2156-11-94

* Kearney, M., & Porter, W. (2009). Mechanistic niche modelling: Combining physiological and spatial data to predict species’ ranges. *Ecology Letters, 12*(4), 334–350. https://doi.org/10.1111/j.1461-0248.2008.01277.x

* Lachenbruch, P. A., & Goldstein, M. (1979). Discriminant analysis. *Biometrics, 35*(1), 69–85. https://doi.org/10.2307/2529937

* Miller, J. M., Cullingham, C. I., & Peery, R. M. (2020). The influence of a priori grouping on inference of genetic clusters: Simulation study and literature review of the DAPC method. *Heredity, 125*(5), 269–280. https://doi.org/10.1038/s41437-020-0348-2

* Prajzlerová, D., Barták, V., Balej, P., Moudrý, V., & Šímová, P. (2025). The time of acquisition of multispectral predictors matters: The role of seasonality in bird species distribution models. *Ecography*. https://doi.org/10.1002/ecog.07935

* Pyron, R. A., O’Connell, K. A., Duncan, S. C., Burbrink, F. T., & Beamer, D. A. (2023). Speciation hypotheses from phylogeographic delimitation yield an integrative taxonomy for Seal Salamanders (*Desmognathus monticola*). *Systematic Biology, 72*(1), 179–197. https://doi.org/10.1093/sysbio/syac065

* Schoener, T. W. (1968). The Anolis lizards of Bimini: Resource partitioning in a complex fauna. *Ecology, 49*(4), 704–726. https://doi.org/10.2307/1935534

* Soberón, J. (2007). Grinnellian and Eltonian niches and geographic distributions of species. *Ecology Letters, 10*(12), 1115–1123. https://doi.org/10.1111/j.1461-0248.2007.01107.x

* Zimmermann, N. E., Yoccoz, N. G., Edwards, T. C., Meier, E. S., Thuiller, W., Guisan, A., Schmatz, D. R., & Pearman, P. B. (2009). Climatic extremes improve predictions of spatial patterns of tree species. *Proceedings of the National Academy of Sciences, 106*, 19723–19728. https://doi.org/10.1073/pnas.0901643106


## Citation
Please cite the *NicheDiv* framework as follows:

Schönberger, D., MacDonald, Z. G., Schmidt, B. C., & Dupuis, J. R. *NicheDiv*: A DAPC framework to quantify niche divergence across highly multivariate environmental space. bioRxiv. https://doi.org/10.64898/2026.06.19.733388 


## License
*NicheDiv* is released under the MIT License. See the `LICENSE` file for details.
