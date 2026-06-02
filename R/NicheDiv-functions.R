#' @importFrom graphics abline hist legend par
#' @importFrom methods as
#' @importFrom stats complete.cases median na.omit ppois predict sd setNames
#' @importFrom utils installed.packages capture.output head read.csv tail unzip write.csv
#' @importFrom fields rdist.earth
#' @importFrom spdep knearneigh nb2listw moran.test
#' @importFrom magrittr %>%
#' @importFrom tidyselect where
#' @importFrom terra nlyr
#' @importFrom ggplot2 ggplot aes geom_density scale_fill_manual guides guide_legend labs theme_classic theme element_line element_text geom_rug scale_colour_manual ggsave geom_rect annotate coord_cartesian scale_y_continuous expansion geom_col scale_x_continuous facet_wrap element_rect
NULL

## Function to convert all integer columns to numeric
#' Convert integer columns to numeric
#'
#' Convert all integer columns in a data frame-like object to numeric columns
#' while leaving all non-integer columns unchanged. This is useful when imported
#' occurrence or background data contain integer storage modes that should be
#' treated as numeric in downstream calculations.
#'
#' @param dataframe A `data.frame` or tibble.
#'
#' @return A `data.frame` with the same dimensions, row order, and column names
#'   as `dataframe`, except that any columns of type integer are converted to
#'   numeric.
#'
#' @export
convert.integer.to.numeric <- function(dataframe) {
  if (!is.data.frame(dataframe)) stop("dataframe must be a data.frame or tibble")
  integer_cols <- vapply(dataframe, is.integer, FUN.VALUE = logical(1)) #identify integer columns
  if (any(integer_cols)) dataframe[integer_cols] <- lapply(dataframe[integer_cols], as.numeric) #convert to numeric
  return(dataframe)
}


## Function to crop background to buffered extent
#' Crop background points to a buffered occurrence extent
#'
#' Restrict a background dataset to the buffered accessible area defined from the
#' occurrence coordinates. The function converts occurrence and background
#' coordinates to `sf`, chooses or validates a metric CRS for buffering, builds a
#' clipping geometry around the occurrences, and returns only background rows
#' that intersect that buffered region.
#'
#' @param occurrence.data A `data.frame` or tibble containing occurrence records.
#'   It must contain the coordinate columns specified by `longitude.col` and
#'   `latitude.col`.
#' @param background.data A `data.frame` or tibble containing candidate
#'   background points to be clipped to the buffered extent. It must contain the
#'   same coordinate columns as `occurrence.data`.
#' @param latitude.col A single character string giving the latitude column name
#'   in both input tables. Default: `"Latitude"`.
#' @param longitude.col A single character string giving the longitude column
#'   name in both input tables. Default: `"Longitude"`.
#' @param CRS Coordinate reference system of the input coordinates, supplied as
#'   an EPSG string (for example `"EPSG:4326"`) or an `sf` CRS object. If the
#'   CRS is geographic, longitude/latitude ranges are checked and buffering is
#'   performed after internally projecting to a metric CRS. If the CRS is already
#'   projected in meters, geometric operations are performed directly in that CRS.
#' @param buffer.dist.meters A single positive numeric value giving the buffer
#'   distance, in meters, used to expand the occurrence-based geometry.
#' @param buffer.method Character string specifying how the base geometry is
#'   built before buffering. One of `"hull"` (default; convex hull of all
#'   occurrences), `"points"` (union of per-point buffers), `"alpha"` (concave
#'   hull built with `concaveman`), or `"bbox"` (bounding box around all
#'   occurrences).
#' @param alpha A single positive numeric value controlling concavity when
#'   `buffer.method = "alpha"`. Ignored for the other methods.
#' @param verbose Logical; if `TRUE`, progress messages about dropped rows and
#'   retained background points are printed.
#'
#' @details
#' The convex hull (`buffer.method = "hull"`) is the default geometry because it
#' provides a robust, simple, and biologically meaningful approximation of the
#' accessible area (`M`; Soberón & Peterson, 2005; Barve et al., 2011; Owens et
#' al., 2013). The convex hull traces the outer extent of occurrence records and
#' can approximate the area that has likely been accessible to the species,
#' especially when expanded by an ecologically informed dispersal buffer. Compared
#' with bounding boxes, convex hulls usually avoid strongly overinflated areas,
#' and compared with concave hulls, they avoid dependence on an additional shape
#' parameter.
#'
#' Three alternative geometries are also available. Concave hulls
#' (`buffer.method = "alpha"`) can follow the detailed shape of occurrence
#' distributions and may better approximate complex ranges, but they are
#' sensitive to the concavity parameter (`alpha`; Edelsbrunner & Mücke, 1994;
#' Pateiro-López & Rodríguez-Casal, 2010; Fourcade, 2016). Lower `alpha` values
#' produce more concave outlines that fit the points more closely but may create
#' holes, fragmentation, or unstable geometries. Higher `alpha` values produce
#' smoother outlines that increasingly resemble the convex hull. The default
#' `alpha = 3` is intended as a robust compromise that captures some range
#' complexity while reducing over-fragmentation (Fourcade, 2016;
#' Pateiro-López & Rodríguez-Casal, 2010).
#'
#' Per-point buffers (`buffer.method = "points"`) generate local circular buffers
#' around each occurrence record and unite those buffers into a single clipping
#' geometry. This can be useful for fragmented or disjunct distributions because
#' it does not force distant occurrence clusters into one continuous polygon
#' (Anderson & Raza, 2010; Fourcade et al., 2014). However, per-point buffers can
#' also produce patchy outlines and are most useful when background sampling is
#' dense enough to retain sufficient points within the buffered areas.
#'
#' Bounding boxes (`buffer.method = "bbox"`) provide the most inclusive geometry
#' by enclosing all occurrence records within a rectangle defined by the minimum
#' and maximum coordinate values. Bounding boxes are simple and can retain larger
#' background sample sizes, but they are usually the least realistic approximation
#' of the accessible area and can substantially overestimate the environmental
#' conditions available to a species (Phillips et al., 2006; Peterson et al.,
#' 2007, 2011; VanDerWal et al., 2009). They should generally be avoided unless
#' occurrence sample sizes or very small species distributions make the other
#' methods too restrictive.
#'
#' @return A filtered `data.frame` containing only the rows of
#'   `background.data` whose coordinates fall within the buffered clipping
#'   geometry. Column structure and row order from the retained background rows
#'   are preserved.
#'
#' @references
#' Anderson, R. P., & Raza, A. (2010). The effect of the extent of the study
#'   region on GIS models of species geographic distributions and estimates of
#'   niche evolution: preliminary tests with montane rodents (genus Nephelomys)
#'   in Venezuela. \emph{Journal of Biogeography}, 37(7), 1378-1393.
#'   https://doi.org/10.1111/j.1365-2699.2010.02290.x
#'
#' Barve, N., Barve, V., Jiménez-Valverde, A., Lira-Noriega, A., Maher, S. P.,
#'   Peterson, A. T., Soberón, J., & Villalobos, F. (2011). The crucial role of
#'   the accessible area in ecological niche modeling and species distribution
#'   modeling. \emph{Ecological Modelling}, 222(11), 1810-1819.
#'   https://doi.org/10.1016/j.ecolmodel.2011.02.011
#'
#' Edelsbrunner, H., & Mücke, E. P. (1994). Three-dimensional alpha shapes.
#'   \emph{ACM Transactions on Graphics}, 13(1), 43-72.
#'   https://doi.org/10.1145/174462.156635
#'
#' Fourcade, Y. (2016). Comparing species distributions modelled from occurrence
#'   data and from expert-based range maps. Implication for predicting range
#'   shifts with climate change. \emph{Ecological Informatics}, 36, 8-14.
#'   https://doi.org/10.1016/j.ecoinf.2016.09.002
#'
#' Fourcade, Y., Engler, J. O., Rödder, D., & Secondi, J. (2014). Mapping species
#'   distributions with MAXENT using a geographically biased sample of presence
#'   data: A performance assessment of methods for correcting sampling bias.
#'   \emph{PLOS ONE}, 9(5), e97122.
#'   https://doi.org/10.1371/journal.pone.0097122
#'
#' Owens, H. L., Campbell, L. P., Dornak, L. L., Saupe, E. E., Barve, N.,
#'   Soberón, J., Ingenloff, K., Lira-Noriega, A., Hensz, C. M., Myers, C. E.,
#'   & Peterson, A. T. (2013). Constraints on interpretation of ecological niche
#'   models by limited environmental ranges on calibration areas.
#'   \emph{Ecological Modelling}, 263, 10-18.
#'   https://doi.org/10.1016/j.ecolmodel.2013.04.011
#'
#' Pateiro-López, B., & Rodríguez-Casal, A. (2010). Generalizing the convex hull
#'   of a sample: The R package alphahull. \emph{Journal of Statistical
#'   Software}, 34(5).
#'   https://doi.org/10.18637/jss.v034.i05
#'
#' Peterson, A. T., Papeş, T., & Eaton, M. (2007). Transferability and model
#'   evaluation in ecological niche modeling: A comparison of GARP and Maxent.
#'   \emph{Ecography}, 30(4), 550-560.
#'   https://doi.org/10.1111/j.0906-7590.2007.05102.x
#'
#' Peterson, A. T., Soberón, J., Pearson, R. G., Anderson, R. P.,
#'   Martínez-Meyer, E., Nakamura, M., & Araújo, M. B. (2011). \emph{Ecological
#'   niches and geographic distributions}. Princeton University Press.
#'   https://doi.org/10.1515/9781400840670
#'
#' Phillips, S. J., Anderson, R. P., & Schapire, R. E. (2006). Maximum entropy
#'   modeling of species geographic distributions. \emph{Ecological Modelling},
#'   190(3-4), 231-259.
#'   https://doi.org/10.1016/j.ecolmodel.2005.03.026
#'
#' Soberón, J., & Peterson, A. T. (2005). Interpretation of models of fundamental
#'   ecological niches and species' distributional areas. \emph{Biodiversity
#'   Informatics}, 2(0).
#'   https://doi.org/10.17161/bi.v2i0.4
#'
#' VanDerWal, J., Shoo, L. P., Graham, C., & Williams, S. E. (2009). Selecting
#'   pseudo-absence data for presence-only distribution modeling: How far should
#'   you stray from what you know? \emph{Ecological Modelling}, 220(4), 589-594.
#'   https://doi.org/10.1016/j.ecolmodel.2008.11.010
#'
#' @export
crop.background.buffered <- function(occurrence.data, #input data.frame or tibble with occurrence records (must include latitude and longitude cols)
                                     background.data, #background dataframe
                                     latitude.col = "Latitude", #name of latitude column
                                     longitude.col = "Longitude", #name of longitude column
                                     CRS = "EPSG:4326", #coordinate reference system (WKT or EPSG code)
                                     buffer.dist.meters, #buffer distance in meters
                                     buffer.method = c("hull", "points", "alpha", "bbox"), #buffer.method to build clipping geometry
                                     alpha = 3, #concavity parameter (only if buffer.method = "alpha")
                                     verbose = TRUE #show output messages
) {

  # Validate inputs
  if (!is.data.frame(occurrence.data)) stop("occurrence.data must be a data.frame or tibble")
  if (!is.data.frame(background.data)) stop("background.data must be a data.frame or tibble")
  if (nrow(occurrence.data) == 0) stop("occurrence.data has no rows")
  if (nrow(background.data) == 0) stop("background.data has no rows")
  if (!is.character(latitude.col) || length(latitude.col) != 1L) stop("latitude.col must be a single character string")
  if (!is.character(longitude.col) || length(longitude.col) != 1L) stop("longitude.col must be a single character string")
  if (identical(latitude.col, longitude.col)) stop("latitude.col and longitude.col must refer to different columns")
  if (!(latitude.col %in% colnames(occurrence.data) && longitude.col %in% colnames(occurrence.data))) stop("latitude.col and/or longitude.col not found in occurrence.data")
  if (!(latitude.col %in% colnames(background.data) && longitude.col %in% colnames(background.data))) stop("latitude.col and/or longitude.col not found in background.data")
  if (!is.numeric(occurrence.data[[latitude.col]]) || !is.numeric(occurrence.data[[longitude.col]])) stop("Latitude and longitude columns must be numeric")
  if (!is.numeric(background.data[[latitude.col]]) || !is.numeric(background.data[[longitude.col]])) stop("Latitude and longitude columns in background.data must be numeric")
  CRS <- tryCatch(sf::st_crs(CRS), error = function(e) sf::st_crs(NA))
  if (is.na(CRS)) stop("CRS must be a valid CRS accepted by sf::st_crs, for example 'EPSG:4326' or sf::st_crs(4326)")
  if (!is.numeric(buffer.dist.meters) || length(buffer.dist.meters) != 1L || !is.finite(buffer.dist.meters) || buffer.dist.meters <= 0) stop("buffer.dist.meters must be a single positive number (meters)")
  buffer.method <- match.arg(buffer.method)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0) stop("alpha must be a single positive number")
  if (!is.logical(verbose) || length(verbose) != 1L) stop("verbose must be TRUE or FALSE")

  # Drop rows with no coordinates
  occ_before <- nrow(occurrence.data)
  bg_before <- nrow(background.data)
  occurrence.data <- occurrence.data[stats::complete.cases(occurrence.data[, c(longitude.col, latitude.col), drop = FALSE]), , drop = FALSE]
  background.data <- background.data[stats::complete.cases(background.data[, c(longitude.col, latitude.col), drop = FALSE]), , drop = FALSE]
  occ_dropped <- occ_before - nrow(occurrence.data)
  bg_dropped <- bg_before - nrow(background.data)
  if (occ_dropped > 0 && verbose) message("Dropped ", occ_dropped, " of ", occ_before, " occurrence rows with missing coordinates")
  if (bg_dropped > 0 && verbose) message("Dropped ", bg_dropped, " of ", bg_before, " background rows with missing coordinates")
  if (nrow(occurrence.data) == 0) stop("No valid occurrence rows after removing NA coords")
  if (nrow(background.data) == 0) stop("No valid background rows after removing NA coords")

  # Warn on large buffer distances (>10,000 km)
  if (buffer.dist.meters > 1e7) warning("buffer.dist.meters is very large (>10,000 km) - ensure units and CRS are correct")

  # Convert to sf (x = Longitude, y = Latitude)
  occurrence_points_ll <- sf::st_as_sf(occurrence.data, coords = c(longitude.col, latitude.col), crs = CRS)
  background_points_ll <- sf::st_as_sf(background.data, coords = c(longitude.col, latitude.col), crs = CRS)

  # Validate longitude and latitude inputs if CRS is geographic
  if (sf::st_is_longlat(occurrence_points_ll)) {
    lon_occ <- occurrence.data[[longitude.col]]
    lat_occ <- occurrence.data[[latitude.col]]
    lon_bg <- background.data[[longitude.col]]
    lat_bg <- background.data[[latitude.col]]
    if (any(abs(lon_occ) > 180) || any(abs(lat_occ) > 90) || any(abs(lon_bg) > 180) || any(abs(lat_bg) > 90)) stop("Coordinate validation failed: values exceed valid lon/lat ranges (|lon| <= 180, |lat| <= 90) for a geographic CRS")
  }

  # Set CRS (always work in meters: UTM for modest extents; LAEA for large/multi-zone/polar)
  choose_metric_crs <- function(points_ll) {
    bounding_box <- sf::st_bbox(points_ll)
    lon_span <- bounding_box["xmax"] - bounding_box["xmin"]
    lat_span <- bounding_box["ymax"] - bounding_box["ymin"]
    if (lon_span > 180 && sf::st_is_longlat(points_ll)) warning("Occurrence extent spans more than 180 degrees longitude; antimeridian-crossing datasets may need manual CRS handling")
    ctr <- sf::st_coordinates(sf::st_centroid(sf::st_as_sfc(bounding_box)))
    ctr_lon <- ctr[1]
    ctr_lat <- ctr[2]
    near_pole <- (bounding_box["ymax"] > 83.5) || (bounding_box["ymin"] < -80) #UTM validity limits
    large_extent <- (lon_span > 12) || (lat_span > 12) #>~ two UTM zones or very tall
    if (!large_extent && !near_pole) {
      zone <- max(1, min(60, floor((ctr_lon + 180) / 6) + 1)) #UTM zone
      epsg <- if (ctr_lat >= 0) 32600 + zone else 32700 + zone #north/south
      return(sf::st_crs(epsg))
    } else {
      laea_wkt <- paste0("+proj=laea +lat_0=", ctr_lat, " +lon_0=", ctr_lon,
                         " +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
      return(sf::st_crs(laea_wkt))
    }
  }
  if (sf::st_is_longlat(occurrence_points_ll)) {
    metric_crs <- choose_metric_crs(occurrence_points_ll) #pick metric CRS from lon/lat
  } else {
    units_gdal <- tryCatch(sf::st_crs(occurrence_points_ll)$units_gdal, error = function(e) NA_character_) #units of provided CRS
    if (is.na(units_gdal) || !grepl("metre|meter", tolower(units_gdal))) {
      metric_crs <- choose_metric_crs(sf::st_transform(occurrence_points_ll, "EPSG:4326")) #fallback based on lon/lat
    } else {
      metric_crs <- sf::st_crs(occurrence_points_ll)
    }
  }
  occurrence_points_utm <- sf::st_transform(occurrence_points_ll, metric_crs) #to metric CRS
  background_points_utm <- sf::st_transform(background_points_ll, metric_crs) #to metric CRS

  # Build clipping geometry by buffer.method
  if (buffer.method == "hull") {
    occurrence_union <- sf::st_union(occurrence_points_utm)
    base_geometry <- sf::st_convex_hull(occurrence_union)
    buffer_geometry <- sf::st_buffer(base_geometry, dist = buffer.dist.meters)

  } else if (buffer.method == "points") {
    buffers_each_point <- sf::st_buffer(occurrence_points_utm, dist = buffer.dist.meters)
    buffer_geometry <- sf::st_union(buffers_each_point)

  } else if (buffer.method == "alpha") {
    coords_mat <- sf::st_coordinates(occurrence_points_utm)
    if (nrow(unique(coords_mat[, c("X", "Y"), drop = FALSE])) < 3) stop("Too few unique occurrence coordinates for alpha hull")
    alpha_geometry <- concaveman::concaveman(occurrence_points_utm, concavity = alpha)
    if (inherits(alpha_geometry, "sf") && all(sf::st_is_empty(alpha_geometry))) stop("alpha buffer.method failed - choose other buffer.method")
    buffer_geometry <- sf::st_buffer(alpha_geometry, dist = buffer.dist.meters)

  } else if (buffer.method == "bbox") {
    bbox_geometry <- sf::st_as_sfc(sf::st_bbox(occurrence_points_utm))
    buffer_geometry <- sf::st_buffer(bbox_geometry, dist = buffer.dist.meters)
  }

  # Repair invalid geometries
  suppressWarnings(buffer_geometry <- sf::st_make_valid(buffer_geometry))
  if (all(sf::st_is_empty(buffer_geometry))) stop("Buffer geometry is empty - check inputs and buffer distance")

  # Keep background points that intersect buffer
  background_intersects <- sf::st_intersects(background_points_utm, buffer_geometry) #sparse list
  keep_rows <- lengths(background_intersects) > 0
  if (sum(keep_rows) == 0) stop("No background points fell within buffered region")
  if (isTRUE(verbose)) message("Retained ", sum(keep_rows), " of ", nrow(background.data), " background points after buffering (buffer.method = ", buffer.method, ")")

  # Return results
  return(background.data[keep_rows, , drop = FALSE])
}


## Function to sample down background or occurrence data to N rows
#' Down-sample rows from a data frame
#'
#' Randomly sample rows without replacement from a `data.frame`-like object. By
#' default, rows with fewer missing values are favored probabilistically using
#' Poisson-tail weights; alternatively, uniform random sampling can be used.
#'
#' @param dataframe A `data.frame`, tibble, or `sf` object to sample from.
#' @param N.rows A positive integer giving the number of rows to retain.
#'   If `N.rows` is greater than or equal to the number of available rows, the
#'   original object is returned unchanged.
#' @param prioritize.NA.poisson Logical; if `TRUE` (default), rows with fewer
#'   missing values are given higher sampling probability using Poisson-tail
#'   weights. If `FALSE`, rows are sampled uniformly without replacement.
#' @param poisson.lambda Optional positive numeric value controlling the
#'   steepness of the Poisson-tail weighting applied to row-wise NA counts. If
#'   `NULL` (default), the median row-wise NA count is used.
#' @param seed A single numeric value used to set the random seed for
#'   reproducible sampling (default: 1).
#'
#' @return An object containing the sampled rows. The returned object has the
#'   same class and columns as `dataframe` whenever possible (including `sf`
#'   objects), but includes only the sampled subset of rows.
#'
#' @export
sample.down <- function(dataframe, #input dataframe to sample from
                        N.rows, #number of rows to sample
                        prioritize.NA.poisson = TRUE, #if TRUE, prefer rows with fewer NAs (probabilistically based on Poisson distribution)
                        poisson.lambda = NULL, #set Poisson lambda (if NULL, median NA count is used)
                        seed = 1 #set seed for reproducibility
) {

  # Validate inputs
  if (!(is.data.frame(dataframe) || inherits(dataframe, "sf"))) stop("dataframe must be a data.frame or tibble or sf object")
  if (nrow(dataframe) == 0) stop("dataframe has no rows")
  if (!is.numeric(N.rows) || length(N.rows) != 1L || !is.finite(N.rows) || N.rows < 1) stop("N.rows must be a single positive integer (>= 1)")
  if (abs(N.rows - round(N.rows)) > .Machine$double.eps^0.5) stop("N.rows must be an integer")
  N.rows <- as.integer(N.rows) #ensure integer
  if (!is.logical(prioritize.NA.poisson) || length(prioritize.NA.poisson) != 1L) stop("prioritize.NA.poisson must be TRUE or FALSE")
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) stop("seed must be a single finite number")
  if (!is.null(poisson.lambda) && (!is.numeric(poisson.lambda) || length(poisson.lambda) != 1L || !is.finite(poisson.lambda) || poisson.lambda <= 0)) {
    stop("poisson.lambda must be NULL (lambda determined automatically) or a single positive number (recommended: 1.0)")
  }

  # Set seed
  set.seed(seed)

  # Check N.rows
  if (N.rows >= nrow(dataframe)) {
    if (N.rows > nrow(dataframe)) message("Requested N.rows (", N.rows, ") exceeds available rows (", nrow(dataframe), ") - using all available rows instead")
    return(dataframe)
  }

  # If not prioritizing NA, do uniform sampling
  if (!isTRUE(prioritize.NA.poisson)) {
    keep_idx <- sample.int(n = nrow(dataframe), size = N.rows, replace = FALSE)
    return(dataframe[keep_idx, , drop = FALSE])
  }

  # Count NAs per row (drop geometry if present)
  df_no_geom <- if (inherits(dataframe, "sf")) sf::st_drop_geometry(dataframe) else dataframe
  na_counts <- rowSums(is.na(df_no_geom))
  if (all(is.na(na_counts))) stop("All rows contain only NA values")

  # Choose lambda: user-specified or adaptive (median NA count)
  lambda_value <- if (is.null(poisson.lambda)) median(na_counts) else poisson.lambda

  # Set weights: Poisson tail P(X >= n | lambda) = 1 - F(n - 1; lambda) - higher NA means smaller weight
  weights <- 1 - ppois(pmax(na_counts, 0) - 1, lambda = lambda_value)

  # Guard against degenerate weights
  if (!any(is.finite(weights))) {
    keep_idx <- sample.int(n = nrow(dataframe), size = N.rows, replace = FALSE) #fallback to uniform
    return(dataframe[keep_idx, , drop = FALSE])
  }

  # If too few positive-weight rows to fill sample, take all positives then fill uniformly
  positive_idx <- which(is.finite(weights) & weights > 0)
  if (length(positive_idx) >= N.rows) {
    keep_idx <- sample.int(n = nrow(dataframe), size = N.rows, replace = FALSE, prob = weights)
    return(dataframe[keep_idx, , drop = FALSE])
  } else if (length(positive_idx) > 0) {
    must_take <- sample(positive_idx, length(positive_idx)) #shuffle to avoid order bias
    remaining_idx <- setdiff(seq_len(nrow(dataframe)), must_take)
    need <- N.rows - length(must_take)
    fill_idx <- if (need > 0) sample(remaining_idx, size = need, replace = FALSE) else integer(0)
    keep_idx <- c(must_take, fill_idx)
    return(dataframe[keep_idx, , drop = FALSE])
  } else { #if all weights are zero but finite, use uniform sampling
    keep_idx <- sample.int(n = nrow(dataframe), size = N.rows, replace = FALSE)
    return(dataframe[keep_idx, , drop = FALSE])
  }
}


## Function to perform spatial thinning and check remaining spatial autocorrelation
#' Spatially thin occurrence records and evaluate remaining spatial autocorrelation
#'
#' Enforce a minimum nearest-neighbour distance among occurrence records by
#' iteratively removing spatially redundant points. The function runs multiple
#' thinning replicates, prioritizes retention of rows with fewer missing values,
#' and calculates and plots Moran's I values for retained numeric non-coordinate
#' variables as a diagnostic of residual spatial autocorrelation.
#'
#' @param occurrence.data A `data.frame` or tibble containing occurrence records
#'   with longitude and latitude coordinates in decimal degrees.
#' @param latitude.col A single character string giving the latitude column
#'   name (default: `"Latitude"`).
#' @param longitude.col A single character string giving the longitude column
#'   name (default: `"Longitude"`).
#' @param thinning.dist.km A single positive numeric value giving the minimum
#'   allowed distance, in kilometers, between retained occurrence records
#'   (default: `1`).
#' @param exclude.cols Optional character vector of non-environmental columns to
#'   exclude from Moran's I calculations. Coordinate columns are excluded
#'   automatically (default: `NULL`).
#' @param N.thinning.replicates A single positive integer-like numeric value
#'   giving the number of thinning replicates to run. The replicate retaining the
#'   most points is kept; ties are broken by total missingness when
#'   `prioritize.NA = TRUE` (default: `50`).
#' @param calc.Morans.I Logical; if `TRUE`, Moran's I is calculated for retained
#'   numeric non-coordinate variables as a diagnostic of residual spatial
#'   autocorrelation (default: `TRUE`).
#' @param plot.Morans.I Logical; if `TRUE` and `calc.Morans.I = TRUE`, a
#'   histogram of Moran's I values is plotted (default: `TRUE`).
#' @param prioritize.NA Logical; if `TRUE`, rows with more missing values are
#'   preferentially removed during thinning conflicts and, when necessary, the
#'   final replicate is chosen based on the lowest total missingness
#'   (default: `TRUE`).
#' @param seed A single numeric value used to set the random seed for
#'   reproducibility (default: `1`).
#' @param verbose Logical; if `TRUE`, progress messages about thinning and the
#'   Moran's I summary are printed (default: `TRUE`).
#'
#' @details
#' Spatial thinning is used to reduce occurrence clustering, sampling bias, and
#' residual spatial autocorrelation, all of which can influence ecological niche
#' models and downstream environmental comparisons (Dormann et al., 2007; Veloz,
#' 2009; Boria et al., 2014; Fourcade et al., 2014; Aiello-Lammens et al., 2015;
#' Kramer-Schadt et al., 2013; Inman et al., 2021; Lamboley & Fourcade, 2024).
#' Enforcing a minimum nearest-neighbour distance among occurrence records helps
#' reduce the disproportionate influence of spatially clustered samples, which
#' often reflect accessibility, collector effort, or database biases rather than
#' biological density.
#'
#' The default thinning threshold (`thinning.dist.km = 1`) is intended to
#' approximate the spatial resolution of many fine-scale environmental GIS layers
#' commonly used in ecological niche modeling. However, no single thinning
#' distance is universally appropriate. Suitable values depend on predictor
#' resolution, occurrence density, the spatial structure of sampling bias, the
#' biology of the study organism, and the geographic extent of the study system
#' (Aiello-Lammens et al., 2015; Boria et al., 2014; Fourcade et al., 2014;
#' Kramer-Schadt et al., 2013; Inman et al., 2021; Lamboley & Fourcade, 2024).
#' Larger thinning distances may be appropriate when occurrence clustering
#' reflects broad-scale sampling bias, whereas smaller distances may be preferable
#' for narrowly distributed taxa, sparse datasets, fine-resolution predictors, or
#' cases where stronger thinning would remove too many records (Veloz, 2009;
#' Aiello-Lammens et al., 2015; Boria et al., 2014; Kramer-Schadt et al., 2013).
#'
#' Because thinning can reduce but does not necessarily eliminate spatial
#' structure, residual spatial autocorrelation is evaluated with Moran's I
#' (Moran, 1950). Moran's I provides a standardized measure of whether similar
#' environmental values remain spatially clustered after thinning, which is
#' important because spatial autocorrelation can inflate model performance, bias
#' statistical inference, and reduce the effective independence of occurrence
#' records (Legendre, 1993; Fortin & Dale, 2005; Dormann et al., 2007; Veloz,
#' 2009).
#'
#' When sample sizes allow, users may consider increasing the thinning distance
#' until Moran's I values fall below approximately 0.5, corresponding to moderate
#' spatial autocorrelation (Legendre, 1993; Fortin & Dale, 2005). This threshold
#' is intentionally conservative and should be interpreted as a practical
#' diagnostic rather than a universal rule. The goal is to reduce residual spatial
#' dependence to biologically acceptable levels while retaining enough occurrence
#' records for downstream analyses (Dormann et al., 2007).
#'
#' For datasets with more than 5000 valid coordinate records, the function uses
#' a sparse distance-neighbour search instead of constructing a full pairwise
#' distance matrix. This reduces memory use and speeds up thinning for large
#' datasets by around 30%. Nonetheless, runtime can still be substantial when
#' records are highly clustered, the thinning distance is large, or many
#' thinning replicates are requested. For datasets with 5000 or fewer valid
#' coordinate records, the original full pairwise distance matrix thinning
#' algorithm is used.
#'
#' @return A `data.frame` containing the retained occurrence rows after spatial
#'   thinning. Moran's I statistics are produced as printed and/or graphical side
#'   effects and are not returned as part of the output object.
#'
#' @references
#' Aiello-Lammens, M. E., Boria, R. A., Radosavljevic, A., Vilela, B., &
#'   Anderson, R. P. (2015). spThin: An R package for spatial thinning of species
#'   occurrence records for use in ecological niche models. \emph{Ecography},
#'   38(5), 541-545. https://doi.org/10.1111/ecog.01132
#'
#' Boria, R. A., Olson, L. E., Goodman, S. M., & Anderson, R. P. (2014). Spatial
#'   filtering to reduce sampling bias can improve the performance of ecological
#'   niche models. \emph{Ecological Modelling}, 275, 73-77.
#'   https://doi.org/10.1016/j.ecolmodel.2013.12.012
#'
#' Dormann, C. F., M. McPherson, J., B. Araújo, M., Bivand, R., Bolliger, J.,
#'   Carl, G., G. Davies, R., Hirzel, A., Jetz, W., Daniel Kissling, W., Kühn, I.,
#'   Ohlemüller, R., R. Peres-Neto, P., Reineking, B., Schröder, B., M. Schurr,
#'   F., & Wilson, R. (2007). Methods to account for spatial autocorrelation in
#'   the analysis of species distributional data: A review. \emph{Ecography},
#'   30(5), 609-628. https://doi.org/10.1111/j.2007.0906-7590.05171.x
#'
#' Fortin, M.-J., & Dale, M. R. T. (2005). \emph{Spatial analysis: A guide for
#'   ecologists}. Cambridge University Press.
#'
#' Fourcade, Y., Engler, J. O., Rödder, D., & Secondi, J. (2014). Mapping species
#'   distributions with MAXENT using a geographically biased sample of presence
#'   data: A performance assessment of methods for correcting sampling bias.
#'   \emph{PLOS ONE}, 9(5), e97122.
#'   https://doi.org/10.1371/journal.pone.0097122
#'
#' Inman, R., Franklin, J., Esque, T., & Nussear, K. (2021). Comparing sample
#'   bias correction methods for species distribution modeling using virtual
#'   species. \emph{Ecosphere}, 12(3). https://doi.org/10.1002/ecs2.3422
#'
#' Kramer-Schadt, S., Niedballa, J., Pilgrim, J. D., Schröder, B., Lindenborn,
#'   J., Reinfelder, V., Stillfried, M., Heckmann, I., Scharf, A. K., Augeri,
#'   D. M., Cheyne, S. M., Hearn, A. J., Ross, J., Macdonald, D. W., Mathai, J.,
#'   Eaton, J., Marshall, A. J., Semiadi, G., Rustam, R., et al. (2013). The
#'   importance of correcting for sampling bias in MaxEnt species distribution
#'   models. \emph{Diversity and Distributions}, 19(11), 1366-1379.
#'   https://doi.org/10.1111/ddi.12096
#'
#' Lamboley, Q., & Fourcade, Y. (2024). No optimal spatial filtering distance for
#'   mitigating sampling bias in ecological niche models. \emph{Journal of
#'   Biogeography}, 51(9), 1783-1794. https://doi.org/10.1111/jbi.14854
#'
#' Legendre, P. (1993). Spatial autocorrelation: Trouble or new paradigm?
#'   \emph{Ecology}, 74(6), 1659-1673. https://doi.org/10.2307/1939924
#'
#' Moran, P. A. P. (1950). Notes on continuous stochastic phenomena.
#'   \emph{Biometrika}, 37, 17-23.
#'
#' Veloz, S. D. (2009). Spatially autocorrelated sampling falsely inflates
#'   measures of accuracy for presence-only niche models. \emph{Journal of
#'   Biogeography}, 36(12), 2290-2299.
#'   https://doi.org/10.1111/j.1365-2699.2009.02174.x
#'
#' @export
thin.occurrence <- function(occurrence.data, #input data.frame with occurrence records
                            latitude.col = "Latitude", #name of latitude column
                            longitude.col = "Longitude", #name of longitude column
                            thinning.dist.km = 1, #minimum distance (in km) between retained points
                            exclude.cols = NULL, #columns to exclude from Moran's I
                            N.thinning.replicates = 50, #number of thinning replicates
                            calc.Morans.I = TRUE, #whether to calculate Moran's I statistics
                            plot.Morans.I = TRUE, #whether to plot Moran's I distribution (only if calc.Morans.I = TRUE)
                            prioritize.NA = TRUE, #prioritize keeping rows with fewer NAs
                            seed = 1, #set seed for reproducibility
                            verbose = TRUE #show output messages
) {

  # Validate inputs
  if (!is.data.frame(occurrence.data)) stop("occurrence.data must be a data.frame or tibble")
  if (nrow(occurrence.data) == 0) stop("occurrence.data has no rows")
  if (!is.character(latitude.col) || length(latitude.col) != 1L) stop("latitude.col must be a single character string")
  if (!is.character(longitude.col) || length(longitude.col) != 1L) stop("longitude.col must be a single character string")
  if (!(latitude.col %in% colnames(occurrence.data) && longitude.col %in% colnames(occurrence.data))) stop("latitude.col and/or longitude.col not found in occurrence.data")
  if (identical(latitude.col, longitude.col)) stop("latitude.col and longitude.col must refer to different columns")
  if (!is.numeric(occurrence.data[[latitude.col]]) || !is.numeric(occurrence.data[[longitude.col]])) stop("latitude.col and longitude.col must be numeric columns")
  if (!is.numeric(thinning.dist.km) || length(thinning.dist.km) != 1L || !is.finite(thinning.dist.km) || thinning.dist.km <= 0) stop("thinning.dist.km must be a single positive number (km)")
  if (!is.null(exclude.cols) && !is.character(exclude.cols)) stop("exclude.cols must be NULL or a character vector")
  if (!is.numeric(N.thinning.replicates) || length(N.thinning.replicates) != 1L || !is.finite(N.thinning.replicates) || N.thinning.replicates < 1 || N.thinning.replicates %% 1 != 0) stop("N.thinning.replicates must be a single positive integer")
  if (!is.logical(calc.Morans.I) || length(calc.Morans.I) != 1L) stop("calc.Morans.I must be TRUE or FALSE")
  if (!is.logical(plot.Morans.I) || length(plot.Morans.I) != 1L) stop("plot.Morans.I must be TRUE or FALSE")
  if (!is.logical(prioritize.NA) || length(prioritize.NA) != 1L) stop("prioritize.NA must be TRUE or FALSE")
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) stop("seed must be a single finite numeric value")
  if (!is.logical(verbose) || length(verbose) != 1L) stop("verbose must be TRUE or FALSE")

  # Set seed for reproducibility
  set.seed(as.integer(seed))

  # Drop rows with non-finite coordinates
  occurrence.data <- occurrence.data[is.finite(occurrence.data[[latitude.col]]) & is.finite(occurrence.data[[longitude.col]]), , drop = FALSE]
  if (nrow(occurrence.data) < 2) stop("Not enough valid coordinates after filtering non-finite values")

  # Check coordinates
  if (any(abs(occurrence.data[[latitude.col]]) > 90, na.rm = TRUE) || any(abs(occurrence.data[[longitude.col]]) > 180, na.rm = TRUE)) {
    warning("Coordinates appear outside valid longitude or latitude ranges - function assumes degrees",
            " - verify units and columns (", longitude.col, ", ", latitude.col, ")")
  }

  # Extract coordinates
  coords <- occurrence.data[, c(longitude.col, latitude.col)]

  # Calculate NA counts for each row if prioritizing
  na.counts <- if (prioritize.NA) apply(occurrence.data, 1, function(x) sum(is.na(x))) else NULL

  # Create function to run modified thinning algorithm that prioritizes NA
  thinning.algorithm.NA <- function(rec.df.orig, thin.par, reps, na.counts = NULL, prioritize.NA = FALSE) {
    reduced.rec.dfs <- vector("list", reps)
    if (nrow(rec.df.orig) > 5000) {
      if(verbose) message("Using sparse thinning algorithm instead of full pairwise distance matrix for >5000 records (", nrow(rec.df.orig), " records) to speed up computations")
      if(verbose) message("This may still take a while")
      sf_pts <- sf::st_as_sf(rec.df.orig, coords = c(longitude.col, latitude.col), crs = 4326, remove = FALSE)
      old_s2 <- sf::sf_use_s2()
      on.exit(sf::sf_use_s2(old_s2), add = TRUE)
      sf::sf_use_s2(TRUE)
      neighbors <- sf::st_is_within_distance(sf_pts, dist = thin.par * 1000)
      neighbors <- lapply(seq_along(neighbors), function(i) {
        setdiff(as.integer(neighbors[[i]]), i)
      })
      degree.save <- lengths(neighbors)
      for (Rep in seq_len(reps)) {
        active <- rep(TRUE, length(neighbors))
        degree <- degree.save
        while (any(degree[active] > 0L) && sum(active) > 1) {
          active.idx <- which(active)
          RemoveRec <- active.idx[degree[active.idx] == max(degree[active.idx])]
          if (length(RemoveRec) > 1) {
            if (prioritize.NA) {
              max_na <- max(na.counts[RemoveRec]) #among ties, choose row with most NAs
              RemoveRec <- RemoveRec[which(na.counts[RemoveRec] == max_na)]
            }
            if (length(RemoveRec) > 1) RemoveRec <- sample(RemoveRec, 1) #if still tied, random among those
          }
          affected <- neighbors[[RemoveRec]]
          affected <- affected[active[affected]]
          degree[affected] <- degree[affected] - 1L
          degree[RemoveRec] <- 0L
          active[RemoveRec] <- FALSE
        }
        rec.df <- rec.df.orig[active, , drop = FALSE]
        colnames(rec.df) <- colnames(rec.df.orig)
        reduced.rec.dfs[[Rep]] <- rec.df
      }
    } else {
      DistMat.save <- rdist.earth(x1 = rec.df.orig, miles = FALSE) < thin.par
      diag(DistMat.save) <- FALSE
      DistMat.save[is.na(DistMat.save)] <- FALSE
      SumVec.save <- rowSums(DistMat.save)
      df.keep.save <- rep(TRUE, length(SumVec.save))
      for (Rep in seq_len(reps)) {
        DistMat <- DistMat.save
        SumVec <- SumVec.save
        df.keep <- df.keep.save
        while (any(DistMat) && sum(df.keep) > 1) {
          RemoveRec <- which(SumVec == max(SumVec))
          if (length(RemoveRec) > 1) {
            if (prioritize.NA) {
              max_na <- max(na.counts[RemoveRec]) #among ties, choose row with most NAs
              RemoveRec <- RemoveRec[which(na.counts[RemoveRec] == max_na)]
            }
            if (length(RemoveRec) > 1) RemoveRec <- sample(RemoveRec, 1) #if still tied, random among those
          }
          SumVec <- SumVec - DistMat[, RemoveRec]
          SumVec[RemoveRec] <- 0L
          DistMat[RemoveRec, ] <- FALSE
          DistMat[, RemoveRec] <- FALSE
          df.keep[RemoveRec] <- FALSE
        }
        rec.df <- rec.df.orig[df.keep, , drop = FALSE]
        colnames(rec.df) <- colnames(rec.df.orig)
        reduced.rec.dfs[[Rep]] <- rec.df
      }
    }
    reduced.rec.order <- unlist(lapply(reduced.rec.dfs, nrow))
    reduced.rec.order <- order(reduced.rec.order, decreasing = TRUE)
    reduced.rec.dfs <- reduced.rec.dfs[reduced.rec.order]
    return(reduced.rec.dfs)
  }

  # Run thinning algorithm (with NA prioritization)
  thinned.list <- thinning.algorithm.NA(rec.df.orig = coords,
                                        thin.par = thinning.dist.km,
                                        reps = N.thinning.replicates,
                                        na.counts = na.counts,
                                        prioritize.NA = prioritize.NA)
  if (length(thinned.list) == 0) stop("No thinning replicates found - check data")

  # If ties in nrow, choose replicate with lowest total NAs (excluding coords)
  rep.sizes <- sapply(thinned.list, nrow)
  max_size <- max(rep.sizes)
  best_reps <- which(rep.sizes == max_size)
  if (length(best_reps) == 1 || !prioritize.NA) {
    best.rep.idx <- best_reps[1]
  } else {
    na.counts.rep <- sapply(best_reps, function(i) {
      row.matches <- paste(occurrence.data[[latitude.col]], occurrence.data[[longitude.col]]) %in%
        paste(thinned.list[[i]][[latitude.col]], thinned.list[[i]][[longitude.col]])
      sum(rowSums(is.na(occurrence.data[row.matches, setdiff(names(occurrence.data), c(latitude.col, longitude.col)), drop = FALSE])))
    })
    best.rep.idx <- best_reps[which.min(na.counts.rep)]
  }
  best.rep <- thinned.list[[best.rep.idx]]
  colnames(best.rep) <- c(longitude.col, latitude.col)

  # Match thinned coordinates to full data
  original_df <- occurrence.data
  original_df$.row_id <- seq_len(nrow(original_df))
  original_df$.key <- paste(round(original_df[[latitude.col]], 6), round(original_df[[longitude.col]], 6))
  thin_df <- best.rep
  thin_df$.key <- paste(round(thin_df[[latitude.col]], 6), round(thin_df[[longitude.col]], 6))
  sel_ids <- integer(0)
  used <- rep(FALSE, nrow(original_df))
  for (kk in thin_df$.key) {
    idx <- which(original_df$.key == kk & !used)
    if (length(idx)) {
      sel_ids <- c(sel_ids, idx[1])
      used[idx[1]] <- TRUE
    }
  }
  thinned.full <- original_df[sel_ids, , drop = FALSE]
  rownames(thinned.full) <- rownames(occurrence.data)[sel_ids]

  # Report thinning result
  if(verbose) message("Kept ", nrow(thinned.full), " of ", nrow(occurrence.data), " rows after thinning")

  # Moran's I calculation and plotting
  if (calc.Morans.I) {
    numeric.cols <- names(thinned.full)[sapply(thinned.full, is.numeric)]
    exclude.cols <- unique(c(exclude.cols, latitude.col, longitude.col))
    test.vars <- setdiff(numeric.cols, exclude.cols)
    results <- data.frame(
      Variable = character(0),
      Morans_I = numeric(0),
      Expected_I = numeric(0),
      P_value = numeric(0),
      stringsAsFactors = FALSE
    )
    ignored.vars <- c()

    # Build metric coordinates (UTM)
    sf_pts <- sf::st_as_sf(thinned.full, coords = c(longitude.col, latitude.col), crs = 4326)
    cen <- sf::st_coordinates(sf::st_centroid(sf::st_union(sf_pts)))
    zone <- floor((cen[1] + 180) / 6) + 1
    epsg <- if (cen[2] >= 0) 32600 + zone else 32700 + zone
    pts_utm <- sf::st_transform(sf_pts, epsg)
    coords_utm_full <- sf::st_coordinates(pts_utm)
    for (var in test.vars) {
      var.values <- thinned.full[[var]]
      valid.rows <- which(is.finite(var.values) & is.finite(thinned.full[[latitude.col]]) & is.finite(thinned.full[[longitude.col]]))
      if (length(valid.rows) < 20) {
        ignored.vars <- c(ignored.vars, var)
        next
      }
      coords_utm <- coords_utm_full[valid.rows, , drop = FALSE]
      values.sub <- var.values[valid.rows]
      n_pts <- nrow(coords_utm)
      k <- max(1L, min(4L, n_pts - 1L))
      k_max <- max(1L, min(floor(n_pts / 3), 30L, n_pts - 1L))
      repeat {
        knn <- suppressWarnings(knearneigh(coords_utm, k = k))
        nb <- suppressWarnings(spdep::knn2nb(knn, sym = TRUE))
        n.comp <- spdep::n.comp.nb(nb)
        if (n.comp$nc == 1 || k >= k_max) break
        k <- k + 1L
      }
      listw <- nb2listw(nb, style = "W", zero.policy = TRUE)
      test <- tryCatch({
        moran.test(values.sub, listw, zero.policy = TRUE)
      }, error = function(e) NULL)
      if (!is.null(test)) {
        results <- rbind(results, data.frame(
          Variable = var,
          Morans_I = round(test$estimate["Moran I statistic"], 3),
          Expected_I = round(test$estimate["Expectation"], 3),
          P_value = signif(test$p.value, 3)
        ))
      }
    }
    if (plot.Morans.I && nrow(results) > 1) {
      par(mfrow = c(1, 1))
      hist(results$Morans_I,
           main = "Distribution of Moran's I",
           xlab = "Moran's I",
           col = "darkgray",
           border = "white")
      abline(v = median(results$Morans_I), col = "red", lwd = 2, lty = 2)
      legend("topright", legend = "Median", col = "red", lty = 2, lwd = 2)
    }
    if (nrow(results) > 0) {
      if(verbose) message("Median Moran's I: ", round(median(results$Morans_I), 2))
    } else {
      warning("No variables could successfully tested for Moran's I")
    }
  }

  # Remove helper columns
  thinned.full <- thinned.full[, !colnames(thinned.full) %in% c(".row_id", ".key"), drop = FALSE]

  # Return results
  return(thinned.full)
}


## Function to identify and transform skewed variables
#' Identify and transform skewed environmental variables
#'
#' Evaluate numeric variables for skewness and (when appropriate) apply a
#' transformation that reduces absolute skewness. The same transformation rules
#' can optionally be applied to a matching background dataset so that occurrence
#' and background data remain on the same scale.
#'
#' @param data.frame A `data.frame` containing the variables to evaluate and
#'   transform.
#' @param background.dataframe Optional `data.frame` to which the same
#'   transformations selected for `data.frame` are applied. Shared column names
#'   are required when this argument is supplied (default: `NULL`).
#' @param skewness.threshold A single positive numeric value giving the absolute
#'   skewness threshold above which variables are considered for transformation
#'   (default: `1`).
#' @param exclude.cols Optional character vector of column names to leave
#'   unchanged and carry through to the output (default: `NULL`).
#' @param verbose Logical; if `TRUE`, a summary of the number and type of
#'   transformations applied is printed (default: `TRUE`).
#'
#' @details
#' Skewed environmental variables are transformed to stabilize variance, reduce
#' the influence of extreme values, and limit the disproportionate effect of
#' heavy-tailed predictors on multivariate analyses. Strongly right- or
#' left-skewed variables can dominate the first few principal component axes and
#' bias discriminant functions toward variables with long tails, compressed
#' ranges, or extreme observations (Bartlett, 1947; Box & Cox, 1964; Osborne,
#' 2010). Transforming skewed variables before ordination or discrimination can
#' therefore improve comparability among predictors and reduce artifacts caused
#' by differences in distributional shape rather than biological signal.
#'
#' Skewness quantifies the asymmetry of a distribution by comparing the relative
#' weight of its left and right tails (Fisher, 1930; Joanes & Gill, 1998). Values
#' near zero indicate approximate symmetry, positive values indicate a longer
#' right tail, and negative values indicate a longer left tail. The default
#' threshold of `skewness.threshold = 1` treats variables with absolute skewness
#' greater than or equal to one as substantially asymmetric and therefore
#' candidates for transformation (Bulmer, 1979; Doane & Seward, 2011).
#'
#' Different transformations are appropriate for different data ranges because
#' environmental predictors can be strictly positive, non-negative with zeros,
#' bounded proportions, or variables that include negative values. For strictly
#' positive variables, logarithmic and square-root transformations are commonly
#' used to reduce right skew, with logarithmic transformations providing stronger
#' compression of long right tails (Bartlett, 1947; Box & Cox, 1964; Emerson &
#' Stoto, 1983). For non-negative variables that include zeros, transformations
#' that safely accommodate zero values are used to avoid undefined logarithms
#' while still reducing the influence of large values (Cleveland, 1984; Emerson &
#' Stoto, 1983).
#'
#' Continuous proportion variables bounded between zero and one require special
#' consideration because their variance and skewness often depend strongly on
#' proximity to the boundaries. Logit-type transformations can spread values near
#' zero and one while preserving mid-range differences, whereas arcsine
#' square-root transformations have historically been used for bounded
#' proportional data, especially when values are concentrated near the
#' boundaries. Because these transformations can behave differently depending on
#' the distribution of the data, the selected transformation should reduce
#' skewness rather than be applied automatically to all proportional variables
#' (Warton & Hui, 2011).
#'
#' Variables containing negative values require transformations that preserve
#' ordering while allowing the full observed range to be retained. Shifted
#' transformations and power-type transformations are useful in this context
#' because they can reduce asymmetry without discarding negative observations or
#' forcing arbitrary truncation of the data range (Manly, 1976; Tukey, 1977; John
#' & Draper, 1980). Signed power transformations, such as signed cube-root
#' transformations, can be especially useful when variables span both negative
#' and positive values because they preserve sign while reducing the influence of
#' extreme magnitudes.
#'
#' Variables with too few unique finite values, binary variables, and other
#' two-level variables are not transformed because continuous transformations
#' provide little benefit for such predictors and can make their interpretation
#' less clear. When occurrence and background datasets are supplied together, the
#' same selected transformation is applied to both datasets so that environmental
#' values remain directly comparable across occurrence and background samples.
#'
#' @return A named list with three elements: `transformed`, the transformed
#'   occurrence table; `summary`, a `data.frame` describing the transformation
#'   chosen for each evaluated variable, including skewness before and after
#'   transformation; and `background.transformed`, the transformed background
#'   table if `background.dataframe` was supplied, otherwise `NULL`.
#'
#' @references
#' Bartlett, M. S. (1947). The use of transformations. \emph{Biometrics}, 3(1),
#'   39. https://doi.org/10.2307/3001536
#'
#' Box, G. E. P., & Cox, D. R. (1964). An analysis of transformations.
#'   \emph{Journal of the Royal Statistical Society Series B: Statistical
#'   Methodology}, 26(2), 211-243.
#'   https://doi.org/10.1111/j.2517-6161.1964.tb00553.x
#'
#' Bulmer, M. G. (1979). \emph{Principles of statistics}. Dover Publications.
#'
#' Cleveland, W. S. (1984). Graphical methods for data presentation: Full scale
#'   breaks, dot charts, and multibased logging. \emph{The American Statistician},
#'   38(4), 270-280. https://doi.org/10.1080/00031305.1984.10483224
#'
#' Doane, D. P., & Seward, L. E. (2011). Measuring skewness: A forgotten
#'   statistic? \emph{Journal of Statistics Education}, 19(2).
#'   https://doi.org/10.1080/10691898.2011.11889611
#'
#' Emerson, J. D., & Stoto, M. A. (1983). Transforming data. In
#'   \emph{Understanding robust and exploratory data analysis} (pp. 97-128).
#'   Wiley.
#'
#' Fisher, R. A. (1930). Moments and product moments of sampling distributions.
#'   \emph{Proceedings of the London Mathematical Society}, s2-30(1), 199-238.
#'   https://doi.org/10.1112/plms/s2-30.1.199
#'
#' Joanes, D. N., & Gill, C. A. (1998). Comparing measures of sample skewness
#'   and kurtosis. \emph{Journal of the Royal Statistical Society: Series D
#'   (The Statistician)}, 47(1), 183-189.
#'   https://doi.org/10.1111/1467-9884.00122
#'
#' John, J. A., & Draper, N. R. (1980). An alternative family of transformations.
#'   \emph{Applied Statistics}, 29(2), 190.
#'   https://doi.org/10.2307/2986305
#'
#' Manly, B. F. J. (1976). Exponential data transformations. \emph{The
#'   Statistician}, 25(1), 37. https://doi.org/10.2307/2988129
#'
#' Osborne, J. (2010). Improving your data transformations: Applying the Box-Cox
#'   transformation. \emph{Practical Assessment, Research, and Evaluation},
#'   15(1), 12. https://doi.org/10.7275/qbpc-gk17
#'
#' Tukey, J. W. (1977). \emph{Exploratory data analysis}. Addison-Wesley.
#'
#' Warton, D. I., & Hui, F. K. C. (2011). The arcsine is asinine: The analysis of
#'   proportions in ecology. \emph{Ecology}, 92(1), 3-10.
#'   https://doi.org/10.1890/10-0340.1
#'
#' @examples
#' env.data <- data.frame(
#'   species = rep(c("species_1", "species_2"), each = 20),
#'   bio1 = rlnorm(40, meanlog = 1, sdlog = 0.5),
#'   bio12 = rlnorm(40, meanlog = 3, sdlog = 0.4)
#' )
#'
#' result <- transform.skewed.variables(
#'   data.frame = env.data,
#'   exclude.cols = "species",
#'   verbose = FALSE
#' )
#'
#' names(result)
#' head(result$transformed)
#'
#' @rawNamespace export(transform.skewed.variables)
transform.skewed.variables <- function(data.frame, #input data frame containing numeric variables to assess and transform
                                       background.dataframe = NULL, #optional background data frame to apply identical transformations
                                       skewness.threshold = 1, #absolute skewness threshold (variables with |skew| >= this value are transformed)
                                       exclude.cols = NULL, #optional vector of columns to exclude from transformations if present
                                       verbose = TRUE) { #print summary


  # Validate input arguments
  if (missing(data.frame) || !is.data.frame(data.frame)) stop("data.frame must be a data frame containing numeric variables")
  if (!is.null(background.dataframe) && !is.data.frame(background.dataframe)) stop("background.dataframe must be NULL or a data frame")
  if (!is.numeric(skewness.threshold) || length(skewness.threshold) != 1 || skewness.threshold <= 0) stop("skewness.threshold must be a positive numeric value (recommended: 1)")
  if (!is.logical(verbose) || length(verbose) != 1) stop("verbose must be TRUE or FALSE")
  if (!is.null(exclude.cols) && !is.character(exclude.cols)) stop("exclude.cols must be NULL or a character vector of column names")
  if (is.null(exclude.cols))
    exclude.cols <- character(0)
  helper_cols <- c(".row_id", ".key")
  if (!is.null(background.dataframe)) {
    shared.cols <- intersect(colnames(data.frame), colnames(background.dataframe))
    if (length(shared.cols) == 0)
      stop("No shared column names between data.frame and background.dataframe")
    occurrence.only.cols <- setdiff(colnames(data.frame), colnames(background.dataframe))
    occurrence.only.cols.not.excluded <- setdiff(occurrence.only.cols, c(exclude.cols, helper_cols))
    if (length(occurrence.only.cols.not.excluded) > 0) {
      stop(paste0("These columns are present in data.frame but missing from background.dataframe: ",
                  paste(occurrence.only.cols.not.excluded, collapse = ", "),
                  "- add columns to exclude.cols or add matching columns to background.dataframe"))
    }
  }
  numeric.cols <- sapply(data.frame[, !(names(data.frame) %in% c(exclude.cols, helper_cols)), drop = FALSE], is.numeric)
  if (sum(numeric.cols) == 0)
    stop("No numeric columns found in data.frame after removing exclude.cols - nothing to transform")

  data.frame_original <- data.frame
  background.dataframe_original <- background.dataframe

  # Create function to compute unbiased sample skewness (Fisher-Pearson g1)
  detect.skewness <- function(numeric.vector) {
    numeric.vector <- numeric.vector[is.finite(numeric.vector)]
    sample.size <- length(numeric.vector)
    if (sample.size < 3) return(NA_real_)
    standard.deviation <- sd(numeric.vector)
    if (!is.finite(standard.deviation) || standard.deviation == 0) return(0)
    mean.value <- mean(numeric.vector)
    sum(((numeric.vector - mean.value) / standard.deviation)^3) * sample.size / ((sample.size - 1) * (sample.size - 2))
  }

  # Create function to detect proportion variables (0-1)
  is.proportion.variable <- function(numeric.vector, tolerance = 1e-9) {
    variable.range <- range(numeric.vector[is.finite(numeric.vector)], na.rm = TRUE)
    if (!is.finite(variable.range[1])) return(FALSE)
    (variable.range[1] >= -tolerance) && (variable.range[2] <= 1 + tolerance)
  }

  # Create function to logit transform (with shrinkage to avoid +/-Inf at 0 and 1)
  perform.logit.transformation <- function(proportion.vector, epsilon = 1e-3) {
    clipped.values <- pmin(pmax(proportion.vector, 0), 1)
    adjusted.values <- (clipped.values + epsilon) / (1 + 2 * epsilon)
    log(adjusted.values / (1 - adjusted.values))
  }

  # Create function to apply selected numeric transformation
  apply.numeric.transformation <- function(variable.values,
                                           transformation.name,
                                           shift.minimum = NA_real_) {
    variable.values <- ifelse(abs(variable.values) < 1e-10, 0, variable.values) #prevent numerical issues near zero
    transformed.values <- switch(transformation.name,
                                 none = variable.values,
                                 identity = variable.values,
                                 log = log(variable.values),
                                 log1p = log1p(variable.values),
                                 sqrt = sqrt(variable.values),
                                 cuberoot = sign(variable.values) * abs(variable.values)^(1/3),
                                 log1p_shifted = log1p(variable.values - shift.minimum + 1e-6),
                                 sqrt_shifted = sqrt(variable.values - shift.minimum + 1e-6),
                                 variable.values)
    return(transformed.values)
  }

  # Create function to choose best transformation for numeric variables
  choose.best.transformation <- function(variable.values,
                                         reference.values = NULL,
                                         skewness.threshold.local,
                                         prefer.log.for.right.tail = TRUE) {
    finite.values <- variable.values[is.finite(variable.values)]
    reference.finite.values <- c(variable.values, reference.values)
    reference.finite.values <- reference.finite.values[is.finite(reference.finite.values)]
    if (length(unique(finite.values)) < 3) return(list(name = "none", transformed = variable.values, shift.minimum = NA_real_, diagnostics = list(skew.before = NA_real_, skew.after = NA_real_, selection.reason = "none")))
    skew.before <- detect.skewness(variable.values)
    if (!is.finite(skew.before)) return(list(name = "none", transformed = variable.values, shift.minimum = NA_real_, diagnostics = list(skew.before = skew.before, skew.after = skew.before, selection.reason = "skewness_not_finite")))
    if (abs(skew.before) < skewness.threshold.local) return(list(name = "none", transformed = variable.values, shift.minimum = NA_real_, diagnostics = list(skew.before = skew.before, skew.after = skew.before, selection.reason = "skew_below_threshold")))
    if (!length(reference.finite.values)) return(list(name = "none", transformed = variable.values, shift.minimum = NA_real_, diagnostics = list(skew.before = skew.before, skew.after = skew.before, selection.reason = "no_finite_reference_values")))
    minimum.reference.value <- suppressWarnings(min(reference.finite.values, na.rm = TRUE))
    transformation.candidates <- list(identity = variable.values)
    shift.minimum <- NA_real_
    variable.values <- ifelse(abs(variable.values) < 1e-10, 0, variable.values) #prevent numerical issues near zero
    if (is.finite(minimum.reference.value) && minimum.reference.value > 0) {
      transformation.candidates$log <- log(variable.values)
      transformation.candidates$sqrt <- sqrt(variable.values)
    } else if (is.finite(minimum.reference.value) && minimum.reference.value >= 0) {
      transformation.candidates$log1p <- log1p(variable.values)
      transformation.candidates$sqrt <- sqrt(variable.values)
    } else {
      shift.minimum <- minimum.reference.value
      transformation.candidates$log1p_shifted <- log1p(variable.values - shift.minimum + 1e-6)
      transformation.candidates$sqrt_shifted <- sqrt(variable.values - shift.minimum + 1e-6)
      transformation.candidates$cuberoot <- sign(variable.values) * abs(variable.values)^(1/3)
    }
    candidate.skews <- vapply(transformation.candidates, detect.skewness, numeric(1))
    candidate.abs.skews <- abs(candidate.skews)
    candidate.abs.skews[!is.finite(candidate.abs.skews)] <- Inf
    if (prefer.log.for.right.tail && !is.na(skew.before) && skew.before > 3 && "log" %in% names(transformation.candidates)) {
      chosen.name <- "log"
      selection.reason <- "heavy_right_tail"
    } else {
      improvements <- abs(skew.before) - candidate.abs.skews
      if (!any(is.finite(improvements) & improvements > 0)) {
        chosen.name <- "none"
        selection.reason <- "transformation_does_not_improve_skewness"
      } else {
        chosen.name <- names(which.max(improvements))
        selection.reason <- "best_transformation_was_chosen"
      }
    }
    if (chosen.name == "none") {
      transformed.values <- variable.values
      skew.after <- skew.before
    } else {
      transformed.values <- transformation.candidates[[chosen.name]]
      skew.after <- candidate.skews[[chosen.name]]
    }
    list(name = chosen.name,
         transformed = transformed.values,
         shift.minimum = shift.minimum,
         diagnostics = list(skew.before = skew.before,
                            skew.after = skew.after,
                            selection.reason = selection.reason))
  }

  # Normalize input and initialize output lists
  excluded.original <- data.frame[, names(data.frame) %in% exclude.cols, drop = FALSE]
  data.frame <- data.frame[, !(names(data.frame) %in% c(exclude.cols, helper_cols)), drop = FALSE]
  transformed.variable.list <- list()
  summary.record.list <- list()
  transform.counter <- list()
  transform.name.map <- list()

  # Initialize identity mapping for all original columns (including excluded)
  for (nm in colnames(data.frame_original)) transform.name.map[[nm]] <- nm

  # If background provided, initialize parallel copy
  background.transformed.dataframe <- NULL
  if (!is.null(background.dataframe)) background.transformed.dataframe <- background.dataframe

  # Count numeric variables used for transformation summary denominator
  all_numeric_vars <- colnames(data.frame)
  total.used <- length(all_numeric_vars)

  # Iterate through each variable
  for (variable.name in colnames(data.frame)) {
    variable.values <- data.frame[[variable.name]]
    bg.values <- if (!is.null(background.transformed.dataframe) && variable.name %in% colnames(background.transformed.dataframe)) background.transformed.dataframe[[variable.name]] else NULL

    # Handle non-numeric or too-few-unique values
    if (!is.numeric(variable.values) || length(unique(variable.values[is.finite(variable.values)])) < 3) {
      transformed.variable.list[[variable.name]] <- variable.values
      transform.name.map[[variable.name]] <- variable.name
      summary.record.list[[length(summary.record.list) + 1]] <- data.frame(variable = variable.name,
                                                                           transform_chosen = NA_character_,
                                                                           transformed = FALSE,
                                                                           skew_before = NA_real_,
                                                                           skew_after = NA_real_,
                                                                           selection_reason = "too_few_unique_values",
                                                                           stringsAsFactors = FALSE)
      next
    }

    # Detect binary 0/1 variables (presence-absence)
    unique.values <- sort(unique(variable.values[is.finite(variable.values)]))
    if (length(unique.values) == 2 && all(unique.values %in% c(0, 1))) {
      transformed.variable.list[[variable.name]] <- variable.values
      transform.name.map[[variable.name]] <- variable.name
      summary.record.list[[length(summary.record.list) + 1]] <- data.frame(variable = variable.name,
                                                                           transform_chosen = NA_character_,
                                                                           transformed = FALSE,
                                                                           skew_before = NA_real_,
                                                                           skew_after = NA_real_,
                                                                           selection_reason = "binary_variable_no_transformation",
                                                                           stringsAsFactors = FALSE)
      next
    }

    # Handle non-binary variables with only two unique values
    if (length(unique.values) == 2 && !all(unique.values %in% c(0, 1))) {
      skew.tmp <- detect.skewness(variable.values)
      transformed.variable.list[[variable.name]] <- variable.values
      transform.name.map[[variable.name]] <- variable.name
      summary.record.list[[length(summary.record.list) + 1]] <- data.frame(variable = variable.name,
                                                                           transform_chosen = NA_character_,
                                                                           transformed = FALSE,
                                                                           skew_before = round(skew.tmp, 2),
                                                                           skew_after = round(skew.tmp, 2),
                                                                           selection_reason = "two_unique_values_no_transformation",
                                                                           stringsAsFactors = FALSE)
      next
    }

    # Detect and transform proportion (0-1 continuous) variables
    if (is.proportion.variable(c(variable.values, bg.values))) {
      finite.values <- variable.values[is.finite(variable.values)]
      skew.before <- detect.skewness(finite.values)
      if (!is.finite(skew.before) || abs(skew.before) < skewness.threshold) {
        transformed.variable.list[[variable.name]] <- variable.values
        transform.name.map[[variable.name]] <- variable.name
        summary.record.list[[length(summary.record.list) + 1]] <- data.frame(variable = variable.name,
                                                                             transform_chosen = "none",
                                                                             transformed = FALSE,
                                                                             skew_before = round(skew.before, 2),
                                                                             skew_after = round(skew.before, 2),
                                                                             selection_reason = ifelse(!is.finite(skew.before), "skewness_not_finite", "skew_below_threshold"),
                                                                             stringsAsFactors = FALSE)
        next
      }
      logit.values <- perform.logit.transformation(variable.values)
      arcsin.values <- asin(sqrt(pmin(pmax(variable.values, 0), 1)))
      skew.logit <- detect.skewness(logit.values)
      skew.arcsin <- detect.skewness(arcsin.values)
      improvements <- c(abs(skew.before) - abs(skew.logit), abs(skew.before) - abs(skew.arcsin))
      names(improvements) <- c("logit_shrunk", "arcsine_sqrt")
      if (all(improvements <= 0, na.rm = TRUE)) {
        transformed.values <- variable.values
        chosen.transform <- "none"
        skew.after <- skew.before
        selection.reason <- "transformation_does_not_improve_skewness"
      } else {
        chosen.transform <- names(which.max(improvements))
        transformed.values <- switch(chosen.transform,
                                     logit_shrunk = logit.values,
                                     arcsine_sqrt = arcsin.values)
        skew.after <- switch(chosen.transform,
                             logit_shrunk = skew.logit,
                             arcsine_sqrt = skew.arcsin)
        selection.reason <- "best_transformation_was_chosen"
      }
      transformed.variable.list[[variable.name]] <- transformed.values
      transform.name.map[[variable.name]] <- if (chosen.transform == "none") variable.name else paste0(variable.name, "_", chosen.transform)
      if (!is.null(background.transformed.dataframe) && variable.name %in% colnames(background.transformed.dataframe)) {
        bg.values <- background.transformed.dataframe[[variable.name]]
        background.transformed.dataframe[[variable.name]] <- switch(chosen.transform,
                                                                    logit_shrunk = perform.logit.transformation(bg.values),
                                                                    arcsine_sqrt = asin(sqrt(pmin(pmax(bg.values, 0), 1))),
                                                                    bg.values)
      }
      summary.record.list[[length(summary.record.list) + 1]] <- data.frame(variable = variable.name,
                                                                           transform_chosen = chosen.transform,
                                                                           transformed = chosen.transform != "none",
                                                                           skew_before = round(skew.before, 2),
                                                                           skew_after = round(skew.after, 2),
                                                                           selection_reason = selection.reason,
                                                                           stringsAsFactors = FALSE)
      if (chosen.transform != "none") transform.counter[[chosen.transform]] <- (if (is.null(transform.counter[[chosen.transform]])) 0 else transform.counter[[chosen.transform]]) + 1
      next
    }

    # Regular numeric variable transformation
    chosen.result <- choose.best.transformation(variable.values,
                                                reference.values = bg.values,
                                                skewness.threshold.local = skewness.threshold)
    chosen.name <- chosen.result$name
    transformed.variable.list[[variable.name]] <- chosen.result$transformed
    transform.name.map[[variable.name]] <- if (chosen.name == "none") variable.name else paste0(variable.name, "_", chosen.name)
    if (!is.null(background.transformed.dataframe) && variable.name %in% colnames(background.transformed.dataframe)) {
      background.transformed.dataframe[[variable.name]] <- apply.numeric.transformation(bg.values,
                                                                                        chosen.name,
                                                                                        chosen.result$shift.minimum)
    }
    summary.record.list[[length(summary.record.list) + 1]] <- data.frame(variable = variable.name,
                                                                         transform_chosen = chosen.name,
                                                                         transformed = chosen.name != "none",
                                                                         skew_before = round(chosen.result$diagnostics$skew.before, 2),
                                                                         skew_after = round(chosen.result$diagnostics$skew.after, 2),
                                                                         selection_reason = chosen.result$diagnostics$selection.reason,
                                                                         stringsAsFactors = FALSE)
    if (chosen.name != "none") transform.counter[[chosen.name]] <- (if (is.null(transform.counter[[chosen.name]])) 0 else transform.counter[[chosen.name]]) + 1
  }

  # Build transformation summary
  transformation.summary <- do.call(rbind, summary.record.list)

  # Build full data frames from originals
  occ_full <- data.frame_original
  bg_full <- background.dataframe_original

  # Overwrite transformed variables in both occurrence + background (with suffixed names)
  for (old in names(transformed.variable.list)) {
    new <- transform.name.map[[old]]
    occ_full[[new]] <- transformed.variable.list[[old]]
    if (new != old) occ_full[[old]] <- NULL
    if (!is.null(bg_full) && old %in% colnames(bg_full)) {
      if (!is.null(background.transformed.dataframe)) {
        bg_full[[new]] <- background.transformed.dataframe[[old]]
        if (new != old) bg_full[[old]] <- NULL
      }
    }
  }

  # Preserve original column order, but replace names by mapped (transformed) ones
  final_order <- unlist(lapply(colnames(data.frame_original), function(nm) transform.name.map[[nm]]))
  final_order <- unique(final_order[final_order %in% colnames(occ_full)])
  occ_full <- occ_full[, final_order, drop = FALSE]
  if (!is.null(bg_full)) {
    bg_final_order <- final_order[final_order %in% colnames(bg_full)]
    bg_full <- bg_full[, bg_final_order, drop = FALSE]
  }

  # Verbose summary
  if (verbose) {
    transformed.vars <- sum(transformation.summary$transformed)
    message("Transformations completed:")
    message(sprintf("%d of %d variables transformed", transformed.vars, total.used))
    if (length(transform.counter) > 0)
      for (t in names(transform.counter))
        message(sprintf("  %s: %d", t, transform.counter[[t]]))
  }

  # Return results
  list(transformed = occ_full,
       summary = transformation.summary,
       background.transformed = bg_full)

}


## Function to remove variables with low coefficient of variation from occurrence data
#' Remove low-information variables
#'
#' Remove environmental predictors that show negligible variation based on
#' coefficient of variation in either species' occurrence data, while keeping
#' occurrence and background datasets synchronized. Non-numeric variables and
#' variables with too few finite values are also removed before
#' coefficient-of-variation filtering.
#'
#' @param Sp1.occurrence.data A `data.frame` or `sf` object containing occurrence
#'   data for species 1.
#' @param Sp2.occurrence.data A `data.frame` or `sf` object containing occurrence
#'   data for species 2.
#' @param Sp1.background.data A `data.frame` or `sf` object containing background
#'   data for species 1.
#' @param Sp2.background.data A `data.frame` or `sf` object containing background
#'   data for species 2.
#' @param exclude.cols Optional character vector of column names to exclude from
#'   filtering and retain in the returned datasets (default: `NULL`).
#' @param CV.threshold A single non-negative numeric value giving the minimum
#'   coefficient of variation required for a variable to be retained. Variables
#'   falling below this threshold in either species are removed (default: `0.01`).
#' @param verbose Logical; if `TRUE`, messages describing removed and retained
#'   variables are printed (default: `TRUE`).
#'
#' @details
#' Predictors with no or very low variability provide little discriminatory
#' information because they contribute minimally to between-group separation in
#' multivariate analyses. Including near-constant predictors increases
#' dimensionality without adding meaningful signal, can worsen collinearity and
#' numerical conditioning, and may contribute to unstable covariance estimates or
#' singular matrices during principal component analysis or discriminant analysis
#' (Dormann et al., 2013; Greenacre & Primicerio, 2014).
#'
#' The coefficient of variation is used to identify low-information predictors
#' because it measures relative dispersion as the standard deviation scaled by
#' the magnitude of the mean (Pearson, 1896; Sokal & Rohlf, 2012; Zar, 2010).
#' This makes it unitless and scale-invariant, allowing environmental predictors
#' measured in different units and magnitudes to be compared directly before
#' standardization. In contrast, variance is unit-bearing and depends strongly on
#' measurement scale, so variables with larger numeric units can appear more
#' variable even when their relative dispersion is negligible.
#'
#' The default threshold (`CV.threshold = 0.01`) removes variables with less than
#' one percent relative variability compared with their mean. Such variables are
#' effectively constant across occurrences and are unlikely to provide useful
#' discriminatory information. Filtering is applied using both species'
#' occurrence datasets, so a predictor is removed if it has low relative
#' variability in either species. This conservative rule avoids retaining
#' predictors that may separate poorly because one species has little or no
#' usable variation.
#'
#' Variables with too few finite observations are removed because sparse data
#' cannot provide reliable estimates of variability. A coefficient of variation
#' estimated from very few values is highly sensitive to individual observations
#' and can give a misleading impression of predictor informativeness. Removing
#' such variables reduces the risk of retaining predictors whose apparent
#' variation reflects missingness or sampling artifacts rather than biological or
#' environmental signal.
#'
#' For variables with means close to zero, the coefficient of variation can become
#' unstable or undefined because very small denominators can spuriously inflate
#' relative variability. In these cases, a scale-free fallback based on dispersion
#' relative to robust absolute deviation is used to avoid treating near-zero means
#' as evidence of high informativeness (Maronna et al., 2006; Zar, 2010). The
#' same threshold is used for this fallback so that low-information filtering
#' remains comparable across variables.
#'
#' The same retained variable set is applied to both occurrence and background
#' datasets. This keeps species-specific occurrence and background inputs
#' synchronized and prevents downstream analyses from using different predictor
#' spaces for occurrence records and environmental backgrounds.
#'
#' @return A named list with filtered occurrence and background datasets
#'   (`occurrence_Sp1`, `occurrence_Sp2`, `background.Sp1`, `background.Sp2`),
#'   together with character vectors listing variables removed because they were
#'   non-numeric (`dropped_non_numeric`), had too few finite observations
#'   (`dropped_NA_only`), or had low variation (`dropped_lowCV`), plus
#'   `kept_variables` for the retained environmental predictors.
#'
#' @references
#' Dormann, C. F., Elith, J., Bacher, S., Buchmann, C., Carl, G., Carré, G.,
#'   Marquéz, J. R. G., Gruber, B., Lafourcade, B., Leitão, P. J.,
#'   Münkemüller, T., McClean, C., Osborne, P. E., Reineking, B., Schröder, B.,
#'   Skidmore, A. K., Zurell, D., & Lautenbach, S. (2013). Collinearity: A
#'   review of methods to deal with it and a simulation study evaluating their
#'   performance. \emph{Ecography}, 36(1), 27-46.
#'   https://doi.org/10.1111/j.1600-0587.2012.07348.x
#'
#' Greenacre, M., & Primicerio, R. (2014). \emph{Multivariate analysis of
#'   ecological data}. Fundación BBVA.
#'
#' Maronna, R. A., Martin, R. D., & Yohai, V. J. (2006). \emph{Robust
#'   statistics}. Wiley. https://doi.org/10.1002/0470010940
#'
#' Pearson, K. (1896). VII. Mathematical contributions to the theory of
#'   evolution.—III. Regression, heredity, and panmixia. \emph{Philosophical
#'   Transactions of the Royal Society of London. Series A, Containing Papers of
#'   a Mathematical or Physical Character}, 187, 253-318.
#'   https://doi.org/10.1098/rsta.1896.0007
#'
#' Sokal, R. R., & Rohlf, F. J. (2012). \emph{Biometry: The principles and
#'   practice of statistics in biological research} (4th ed.). W. H. Freeman.
#'
#' Zar, J. H. (2010). \emph{Biostatistical analysis}. Pearson.
#'
#' @export
remove.low.CV.vars <- function(Sp1.occurrence.data, #occurrence data for species 1
                               Sp2.occurrence.data, #occurrence data for species 2
                               Sp1.background.data, #background data for species 1
                               Sp2.background.data, #background data for species 2
                               exclude.cols = NULL, #columns to exclude from CV filtering
                               CV.threshold = 0.01, #minimum coefficient of variation threshold
                               verbose = TRUE #whether to show messages
) {

  # Validate inputs
  if (!("data.frame" %in% class(Sp1.occurrence.data) || "sf" %in% class(Sp1.occurrence.data))) stop("Sp1.occurrence.data must be a data.frame or sf object")
  if (!("data.frame" %in% class(Sp2.occurrence.data) || "sf" %in% class(Sp2.occurrence.data))) stop("Sp2.occurrence.data must be a data.frame or sf object")
  if (!("data.frame" %in% class(Sp1.background.data) || "sf" %in% class(Sp1.background.data))) stop("Sp1.background.data must be a data.frame or sf object")
  if (!("data.frame" %in% class(Sp2.background.data) || "sf" %in% class(Sp2.background.data))) stop("Sp2.background.data must be a data.frame or sf object")
  if (is.null(exclude.cols)) exclude.cols <- character(0)
  if (!is.character(exclude.cols)) stop("exclude.cols must be NULL or a character vector")
  if (!is.numeric(CV.threshold) || length(CV.threshold) != 1L || !is.finite(CV.threshold) || CV.threshold < 0) stop("CV.threshold must be a single non-negative finite numeric value (recommended: 0.01)")
  if (CV.threshold < 0.001) warning("CV.threshold = ", CV.threshold, " is extremely low and may fail to remove near-constant variables (recommended: 0.01)")
  if (CV.threshold > 0.2) warning("Warning: CV.threshold = ", CV.threshold, " is high and may remove biologically meaningful predictors with moderate variation (recommended: 0.01)")
  if (!is.logical(verbose) || length(verbose) != 1L) stop("verbose must be TRUE or FALSE")

  # Warn if occurrence datasets differ in column count
  if (verbose && ncol(Sp1.occurrence.data) != ncol(Sp2.occurrence.data)) {
    warning("Sp1.occurrence.data and Sp2.occurrence.data have different numbers of columns: ",
            ncol(Sp1.occurrence.data), " vs ", ncol(Sp2.occurrence.data))
  }

  # Warn if background and occurrence datasets differ in column count
  if (verbose && ncol(Sp1.occurrence.data) != ncol(Sp1.background.data)) {
    warning("Sp1.occurrence.data and Sp1.background.data have different numbers of columns: ",
            ncol(Sp1.occurrence.data), " vs ", ncol(Sp1.background.data))
  }
  if (verbose && ncol(Sp2.occurrence.data) != ncol(Sp2.background.data)) {
    warning("Sp2.occurrence.data and Sp2.background.data have different numbers of columns: ",
            ncol(Sp2.occurrence.data), " vs ", ncol(Sp2.background.data))
  }

  # Check availability of exclude.cols
  missing_exclude <- setdiff(exclude.cols, names(Sp1.occurrence.data))
  if (length(missing_exclude) > 0 && verbose) warning("The following exclude.cols were not found in input data and will be ignored: ", paste(missing_exclude, collapse = ", "))
  exclude.cols <- setdiff(exclude.cols, missing_exclude)

  # Compute number of total usable variables
  all_vars <- union(names(Sp1.occurrence.data), names(Sp2.occurrence.data))
  TOTAL_variables <- length(setdiff(all_vars, exclude.cols))

  # Drop geometry if present
  Sp1_occurrence_df <- if (inherits(Sp1.occurrence.data, "sf")) sf::st_drop_geometry(Sp1.occurrence.data) else Sp1.occurrence.data
  Sp2_occurrence_df <- if (inherits(Sp2.occurrence.data, "sf")) sf::st_drop_geometry(Sp2.occurrence.data) else Sp2.occurrence.data

  # Identify non-numeric variables (except excluded)
  non_numeric_vars <- setdiff(names(Sp1_occurrence_df)[!sapply(Sp1_occurrence_df, is.numeric)], exclude.cols)
  non_numeric_vars <- union(non_numeric_vars, setdiff(names(Sp2_occurrence_df)[!sapply(Sp2_occurrence_df, is.numeric)], exclude.cols))
  X_non_numeric <- length(non_numeric_vars)
  if (X_non_numeric > 0 && verbose) message("Dropped ", X_non_numeric, " of ", TOTAL_variables, " variables due to non-numeric type: ", paste(non_numeric_vars, collapse = ", "))

  # Remove non-numeric variables
  Sp1_occurrence_df <- Sp1_occurrence_df[, setdiff(names(Sp1_occurrence_df), non_numeric_vars), drop = FALSE]
  Sp2_occurrence_df <- Sp2_occurrence_df[, setdiff(names(Sp2_occurrence_df), non_numeric_vars), drop = FALSE]

  # Extract numeric columns
  Sp1.occurrence.numeric <- names(Sp1_occurrence_df)[sapply(Sp1_occurrence_df, is.numeric)]
  Sp2.occurrence.numeric <- names(Sp2_occurrence_df)[sapply(Sp2_occurrence_df, is.numeric)]

  # Determine variables for testing
  candidate_vars_Sp1 <- setdiff(Sp1.occurrence.numeric, exclude.cols)
  candidate_vars_Sp2 <- setdiff(Sp2.occurrence.numeric, exclude.cols)
  test_vars_initial <- union(candidate_vars_Sp1, candidate_vars_Sp2)
  test_vars_initial <- intersect(test_vars_initial, intersect(names(Sp1_occurrence_df), names(Sp2_occurrence_df)))
  if (length(test_vars_initial) == 0) stop("No numeric variables to test after excluding: ", paste(exclude.cols, collapse = ", "))

  # Identify NA-only variables
  is_bad_NA <- function(variable_name) {
    values_species1 <- Sp1_occurrence_df[[variable_name]]
    values_species2 <- Sp2_occurrence_df[[variable_name]]
    (sum(is.finite(values_species1)) < 5) || (sum(is.finite(values_species2)) < 5)
  }
  NA_only_vars <- test_vars_initial[sapply(test_vars_initial, is_bad_NA)]
  Y_NA_only <- length(NA_only_vars)
  if (Y_NA_only > 0 && verbose) message("Dropped ", Y_NA_only, " of ", TOTAL_variables - X_non_numeric, " variables due to only NA or <5 finite values: ", paste(NA_only_vars, collapse = ", "))

  # Remove NA-only variables
  vars_after_NA <- setdiff(test_vars_initial, NA_only_vars)
  if (length(vars_after_NA) == 0) stop("No variables remain after filtering for missing or insufficient finite values")

  # Compute CV (with fallback when mean approx 0)
  compute_CV <- function(variable_values, variable_name) {
    variable_values <- variable_values[is.finite(variable_values)]
    if (length(variable_values) < 5) return(0)
    variable_mean <- mean(variable_values)
    variable_sd <- sd(variable_values)
    if (!is.finite(variable_sd) || variable_sd == 0) return(0) #no variation
    if (!is.finite(variable_mean) || abs(variable_mean) < 1e-7) { #mean near zero: CV invalid
      if (verbose) message("Variable '", variable_name, "' has mean near zero (", format(variable_mean, digits = 4), ") - falling back to SD-MAD-based variability")
      mad_val <- stats::mad(variable_values, constant = 1, na.rm = TRUE)
      if (!is.finite(mad_val) || mad_val == 0) return(0) #MAD=0 so truly constant
      return(variable_sd / (mad_val + 1e-12)) #SD/MAD scale-free variability
    }
    abs(variable_sd / variable_mean)
  }

  # Compute CV for each species
  CV_Sp1 <- sapply(vars_after_NA, function(v) compute_CV(Sp1_occurrence_df[[v]], v))
  CV_Sp2 <- sapply(vars_after_NA, function(v) compute_CV(Sp2_occurrence_df[[v]], v))

  # Identify low-CV variables
  lowCV_vars <- vars_after_NA[CV_Sp1 <= CV.threshold | CV_Sp2 <= CV.threshold]
  Z_lowCV <- length(lowCV_vars)
  if (Z_lowCV > 0 && verbose) message("Dropped ", Z_lowCV, " of ", length(vars_after_NA), " variables due to low variation (CV = ", CV.threshold, "): ", paste(lowCV_vars, collapse = ", "))

  # Determine retained variables
  retained_vars <- setdiff(vars_after_NA, lowCV_vars)
  if (length(retained_vars) == 0) stop("No variables remain after low-CV filtering")
  N_retained <- length(retained_vars)
  if (verbose) message("")
  if (verbose) message("Retained ", N_retained, " of ", TOTAL_variables, " variables after filtering")

  # Filter original datasets
  vars_to_keep <- c(exclude.cols, retained_vars)
  occurrence_Sp1.filtered <- Sp1.occurrence.data[, vars_to_keep, drop = FALSE]
  occurrence_Sp2.filtered <- Sp2.occurrence.data[, vars_to_keep, drop = FALSE]
  if (inherits(Sp1.occurrence.data, "sf")) occurrence_Sp1.filtered <- sf::st_set_geometry(occurrence_Sp1.filtered, sf::st_geometry(Sp1.occurrence.data))
  if (inherits(Sp2.occurrence.data, "sf")) occurrence_Sp2.filtered <- sf::st_set_geometry(occurrence_Sp2.filtered, sf::st_geometry(Sp2.occurrence.data))
  missing_bg1 <- setdiff(retained_vars, names(Sp1.background.data))
  missing_bg2 <- setdiff(retained_vars, names(Sp2.background.data))
  if (length(missing_bg1) > 0) stop("Sp1.background.data is missing retained variables: ", paste(missing_bg1, collapse = ", "))
  if (length(missing_bg2) > 0) stop("Sp2.background.data is missing retained variables: ", paste(missing_bg2, collapse = ", "))
  common_bg_cols1 <- intersect(vars_to_keep, names(Sp1.background.data))
  Sp1.background.data <- Sp1.background.data[, common_bg_cols1, drop = FALSE]
  common_bg_cols2 <- intersect(vars_to_keep, names(Sp2.background.data))
  Sp2.background.data <- Sp2.background.data[, common_bg_cols2, drop = FALSE]

  # Return results
  return(list(
    occurrence_Sp1 = occurrence_Sp1.filtered,
    occurrence_Sp2 = occurrence_Sp2.filtered,
    background.Sp1 = Sp1.background.data,
    background.Sp2 = Sp2.background.data,
    dropped_non_numeric = non_numeric_vars,
    dropped_NA_only = NA_only_vars,
    dropped_lowCV = lowCV_vars,
    kept_variables = retained_vars
  ))
}


## Function to account for environmental analogy bias by removing variables with non-analogous distributions in background data
#' Filter environmental variables to analogous background environment
#'
#' Screen background environmental variables for comparability between two taxa
#' by sequentially removing predictors with too few observations, low variation,
#' low univariate overlap, and low bivariate overlap. The retained variable set
#' is then used to subset the supplied occurrence table.
#'
#' @param Sp1.background.data A `data.frame` or `sf` object containing
#'   background environmental data for species 1.
#' @param Sp2.background.data A `data.frame` or `sf` object containing
#'   background environmental data for species 2.
#' @param Sp1.Sp2.occurrence.data A `data.frame` containing the occurrence data
#'   to subset and return after analogous-variable filtering. This table should
#'   contain the same environmental columns as the background tables, plus any
#'   optional columns listed in `exclude.cols`.
#' @param exclude.cols Optional character vector of columns to retain in the
#'   output occurrence table but exclude from environmental filtering
#'   (default: `NULL`).
#' @param CV.threshold A single non-negative numeric value giving the minimum
#'   coefficient of variation required in both species' backgrounds
#'   (default: `0.01`).
#' @param overlap.threshold A single numeric value between 0 and 1 giving the
#'   minimum acceptable environmental overlap for both the univariate and
#'   bivariate analogy filters. Variables below this threshold are removed. Set
#'   to `0` to skip overlap-based filtering after the low-variation screen
#'   (default: `0.7`).
#' @param max.NA.prop A single numeric value between 0 and 1 giving the maximum
#'   allowed proportion of missing values per row before that row is discarded
#'   from the background tables (default: `0.2`).
#' @param min.rows A single positive integer-like numeric value giving the
#'   minimum number of non-missing observations required per species for a
#'   variable to be tested (default: `15`).
#' @param impute.NA.median Logical; if `TRUE`, eligible missing values are
#'   replaced by the variable median during overlap calculations. If `FALSE`,
#'   incomplete cases are removed instead (default: `TRUE`).
#' @param plot.1D.overlap Logical; if `TRUE`, a histogram of univariate overlap
#'   values is plotted (default: `TRUE`).
#' @param bin.n.2D Optional single numeric value giving the number of bins per
#'   axis used in the bivariate overlap calculations. If `NULL`, the function
#'   determines a value automatically from the effective background sample size
#'   (default: `NULL`).
#' @param max.pairs A single positive integer-like numeric value giving the
#'   maximum number of variable pairs to evaluate in the bivariate overlap step
#'   (default: `10000`).
#' @param use.parallel Logical; if `TRUE`, parallel processing is used for the
#'   bivariate overlap calculations (default: `FALSE`).
#' @param N.cores A single positive integer-like numeric value giving the number
#'   of CPU cores to use when `use.parallel = TRUE` (default: `3`).
#' @param seed A single numeric value used to set the random seed for
#'   reproducible pair subsampling (default: `1`).
#' @param verbose Logical; if `TRUE`, progress messages describing each filtering
#'   step are printed (default: `TRUE`).
#'
#' @details
#' Environmental analogy screening reduces bias caused by comparing species
#' across background environments that are poorly comparable or partly
#' non-analogous. Non-analogous environments can distort estimates of niche
#' similarity, inflate apparent niche divergence, and increase extrapolation risk
#' when species have access to different portions of environmental space (Barve
#' et al., 2011; Peterson et al., 2011; Guisan et al., 2014; Brown & Carnaval,
#' 2019).
#'
#' Low-information predictors are removed before analogy screening because
#' near-constant variables contribute little to environmental discrimination and
#' can make overlap estimates unstable. The coefficient of variation is used
#' because it measures relative variability on a unitless scale, allowing
#' predictors measured in different units and magnitudes to be compared before
#' standardization (Pearson, 1896; Sokal & Rohlf, 2012; Zar, 2010). The default
#' threshold (`CV.threshold = 0.01`) treats variables with less than one percent
#' relative variability as effectively constant.
#'
#' Variables with too few observations are excluded because overlap estimates are
#' unreliable when based on sparse data. The default requirement
#' (`min.rows = 15`) is intended to remove variables with insufficient
#' information for stable estimates while avoiding unnecessarily strict filtering
#' in datasets with limited background availability. Missing data are also
#' restricted because high missingness can make apparent overlap depend more on
#' data availability than on environmental similarity. The default missingness
#' threshold (`max.NA.prop = 0.2`) follows common guidance that moderate levels of
#' missingness can often be handled cautiously, whereas higher levels increase
#' the risk of imputation-driven bias (Harrell, 2015; van Buuren, 2018).
#'
#' Univariate overlap screening evaluates whether each predictor has sufficient
#' marginal environmental overlap between species' accessible background spaces.
#' This step removes variables that are individually non-analogous and therefore
#' likely to drive comparisons through extrapolation rather than shared
#' environmental conditions. Overlap values are interpreted on the same general
#' scale as Schoener's D, where larger values indicate greater similarity between
#' environmental distributions (Schoener, 1968; Warren et al., 2008; Rödder &
#' Engler, 2011).
#'
#' Bivariate overlap screening is included because variables can appear
#' comparable in isolation but become non-analogous when considered jointly.
#' Pairwise screening therefore helps detect non-analogous combinations of
#' predictors and differences in environmental covariance structure that are not
#' visible from univariate comparisons alone (Peterson et al., 2011; Mesgaran et
#' al., 2014). Requiring retained predictors to participate in sufficiently
#' overlapping bivariate combinations reduces the risk that downstream analyses
#' are driven by environmental combinations available to one species but absent
#' from the other.
#'
#' The default overlap threshold (`overlap.threshold = 0.7`) is intended as a
#' practical compromise between retaining informative predictors and excluding
#' variables with poor analogy. Lower values allow greater environmental
#' non-analogy and higher extrapolation risk, whereas stricter thresholds may
#' remove biologically meaningful predictors and leave too few variables for
#' downstream analyses (Barve et al., 2011; Warren et al., 2008; Rödder & Engler,
#' 2011; Peterson et al., 2011).
#'
#' Screening is restricted to univariate and bivariate projections because
#' estimating overlap in higher-dimensional environmental space requires rapidly
#' increasing sample sizes and becomes increasingly sensitive to sparse data. In
#' practice, many major extrapolation risks are detectable in marginal or pairwise
#' dimensions, making univariate and bivariate screening a tractable compromise
#' between computational feasibility and effective detection of environmental
#' non-analogy (Silverman, 1986; Scott, 2015; Mesgaran et al., 2014).
#'
#' @return A filtered `data.frame` containing `Sp1.Sp2.occurrence.data` subset
#'   to the retained analogous environmental variables, plus any columns named in
#'   `exclude.cols` that were present in the input occurrence table.
#'
#' @references
#' Barve, N., Barve, V., Jiménez-Valverde, A., Lira-Noriega, A., Maher, S. P.,
#'   Peterson, A. T., Soberón, J., & Villalobos, F. (2011). The crucial role of
#'   the accessible area in ecological niche modeling and species distribution
#'   modeling. \emph{Ecological Modelling}, 222(11), 1810-1819.
#'   https://doi.org/10.1016/j.ecolmodel.2011.02.011
#'
#' Brown, J., & Carnaval, A. C. (2019). A tale of two niches: Methods, concepts,
#'   and evolution. \emph{Frontiers of Biogeography}, 11(4).
#'   https://doi.org/10.21425/F5FBG44158
#'
#' Guisan, A., Petitpierre, B., Broennimann, O., Daehler, C., & Kueffer, C.
#'   (2014). Unifying niche shift studies. Insights from biological invasions.
#'   \emph{Trends in Ecology & Evolution}, 29(5), 260-269.
#'   https://doi.org/10.1016/j.tree.2014.02.009
#'
#' Harrell, F. E. (2015). \emph{Regression modeling strategies: With applications
#'   to linear models, logistic and ordinal regression, and survival analysis}
#'   (2nd ed.). Springer.
#'
#' Mesgaran, M. B., Cousens, R. D., & Webber, B. L. (2014). Here be dragons: A
#'   tool for quantifying novelty due to covariate range and correlation change
#'   when projecting species distribution models. \emph{Diversity and
#'   Distributions}, 20(10), 1147-1159. https://doi.org/10.1111/ddi.12209
#'
#' Pearson, K. (1896). VII. Mathematical contributions to the theory of
#'   evolution.—III. Regression, heredity, and panmixia. \emph{Philosophical
#'   Transactions of the Royal Society of London. Series A, Containing Papers of
#'   a Mathematical or Physical Character}, 187, 253-318.
#'   https://doi.org/10.1098/rsta.1896.0007
#'
#' Peterson, A. T., Soberón, J., Pearson, R. G., Anderson, R. P.,
#'   Martínez-Meyer, E., Nakamura, M., & Araújo, M. B. (2011). \emph{Ecological
#'   niches and geographic distributions}. Princeton University Press.
#'   https://doi.org/10.1515/9781400840670
#'
#' Rödder, D., & Engler, J. O. (2011). Quantitative metrics of overlaps in
#'   Grinnellian niches: Advances and possible drawbacks. \emph{Global Ecology
#'   and Biogeography}, 20(6), 915-927.
#'   https://doi.org/10.1111/j.1466-8238.2011.00659.x
#'
#' Schoener, T. W. (1968). The anolis lizards of Bimini: Resource partitioning
#'   in a complex fauna. \emph{Ecology}, 49(4), 704-726.
#'   https://doi.org/10.2307/1935534
#'
#' Scott, D. W. (2015). \emph{Multivariate density estimation: Theory, practice,
#'   and visualization} (2nd ed.). Wiley.
#'
#' Silverman, B. W. (1986). \emph{Density estimation for statistics and data
#'   analysis}. Chapman and Hall.
#'
#' Sokal, R. R., & Rohlf, F. J. (2012). \emph{Biometry: The principles and
#'   practice of statistics in biological research} (4th ed.). W. H. Freeman.
#'
#' van Buuren, S. (2018). \emph{Flexible imputation of missing data}. Chapman &
#'   Hall/CRC.
#'
#' Warren, D. L., Glor, R. E., & Turelli, M. (2008). Environmental niche
#'   equivalency versus conservatism: Quantitative approaches to niche evolution.
#'   \emph{Evolution}, 62(11), 2868-2883.
#'   https://doi.org/10.1111/j.1558-5646.2008.00482.x
#'
#' Zar, J. H. (2010). \emph{Biostatistical analysis}. Pearson.
#'
#' @export
filter.analogous.variables <- function(Sp1.background.data, #input data.frame or sf with background data for species 1
                                       Sp2.background.data, #input data.frame or sf with background data for species 2
                                       Sp1.Sp2.occurrence.data, #input data.frame with shared occurrence records to subset/return
                                       exclude.cols = NULL, #columns to exclude from environmental filtering but keep in output
                                       CV.threshold = 0.01, #minimum coefficient of variation required in both backgrounds
                                       overlap.threshold = 0.7, #minimum overlap for univariate and bivariate overlap filters
                                       max.NA.prop = 0.2, #max proportion of NA allowed for imputation
                                       min.rows = 15, #minimum complete cases per species for variable to be used
                                       impute.NA.median = TRUE, #impute median for NA if <= max.NA.prop
                                       plot.1D.overlap = TRUE, #whether to plot histogram of univariate overlap values
                                       bin.n.2D = NULL, #2D bin grid per axis (if NULL, determined automatically)
                                       max.pairs = 10000, #subsample pairs for 2D overlap
                                       use.parallel = FALSE, #use parallel computing for 2D overlap
                                       N.cores = 3, #number of cores
                                       seed = 1, #set seed for reproducibility
                                       verbose = TRUE #whether to print filtering step summaries
) {

  # Input validations
  if (!("sf" %in% class(Sp1.background.data) || is.data.frame(Sp1.background.data))) stop("Sp1.background.data must be sf or data.frame")
  if (!("sf" %in% class(Sp2.background.data) || is.data.frame(Sp2.background.data))) stop("Sp2.background.data must be sf or data.frame")
  if (!is.data.frame(Sp1.Sp2.occurrence.data)) stop("Sp1.Sp2.occurrence.data must be data.frame")
  if (!(is.null(exclude.cols) || is.character(exclude.cols))) stop("exclude.cols must be NULL or a character vector")
  if (!is.numeric(CV.threshold) || length(CV.threshold) != 1 || CV.threshold < 0) stop("CV.threshold must be single non-negative numeric value")
  if (!is.numeric(overlap.threshold) || length(overlap.threshold) != 1 || overlap.threshold < 0 || overlap.threshold > 1) stop("overlap.threshold must be numeric between 0 and 1")
  if (!is.numeric(max.NA.prop) || length(max.NA.prop) != 1 || max.NA.prop < 0 || max.NA.prop > 1) stop("max.NA.prop must be numeric between 0 and 1")
  if (!is.numeric(min.rows) || length(min.rows) != 1 || min.rows < 1 || min.rows != floor(min.rows)) stop("min.rows must be positive integer")
  if (!is.logical(impute.NA.median) || length(impute.NA.median) != 1) stop("impute.NA.median must be single logical value")
  if (!is.logical(plot.1D.overlap) || length(plot.1D.overlap) != 1) stop("plot.1D.overlap must be single logical value")
  if (!is.logical(verbose) || length(verbose) != 1) stop("verbose must be single logical value")
  if (!is.null(bin.n.2D) && (!is.numeric(bin.n.2D) || length(bin.n.2D) != 1 || !is.finite(bin.n.2D) || bin.n.2D < 3)) stop("bin.n.2D must be NULL or single finite numeric >= 3")
  if (!is.null(bin.n.2D) && bin.n.2D > 200) warning("Very large bin.n.2D may be slow")
  if (!is.numeric(max.pairs) || length(max.pairs) != 1 || !is.finite(max.pairs) || max.pairs < 1 || max.pairs %% 1 != 0) stop("max.pairs must be single positive integer")
  if (!is.logical(use.parallel) || length(use.parallel) != 1) stop("use.parallel must be TRUE or FALSE")
  if (!is.numeric(N.cores) || length(N.cores) != 1 || !is.finite(N.cores) || N.cores < 1 || N.cores %% 1 != 0) stop("N.cores must be single positive integer")
  if (!(is.null(seed) || (is.numeric(seed) && length(seed) == 1 && is.finite(seed)))) stop("seed must be NULL or single finite numeric")

  # Processing message
  if (verbose) message("Filtering analogous environmental variables ...")

  # Drop geometry if sf
  drop_geom <- function(df) if ("sf" %in% class(df)) sf::st_drop_geometry(df) else df

  # Remove excluded columns from background before any filtering
  bg1 <- drop_geom(Sp1.background.data)
  bg2 <- drop_geom(Sp2.background.data)
  if (!is.null(exclude.cols)) {
    bg1 <- bg1[, setdiff(names(bg1), exclude.cols), drop = FALSE]
    bg2 <- bg2[, setdiff(names(bg2), exclude.cols), drop = FALSE]
  }

  # Select numeric background columns
  background.numeric.species1 <- bg1 %>% dplyr::select(where(is.numeric))
  background.numeric.species2 <- bg2 %>% dplyr::select(where(is.numeric))

  # Drop rows with too many NA or Inf
  drop.na.rows <- function(input_dataframe, max_na_proportion, impute_missing_values) {
    input_dataframe <- as.data.frame(input_dataframe)
    input_dataframe <- input_dataframe[vapply(input_dataframe, is.numeric, logical(1L))]
    if (ncol(input_dataframe) == 0L) return(input_dataframe[0, , drop = FALSE])
    for (variable_name in names(input_dataframe)) {
      variable_values <- input_dataframe[[variable_name]]
      variable_values[!is.finite(variable_values)] <- NA
      input_dataframe[[variable_name]] <- variable_values
    }
    numeric_matrix <- suppressWarnings(as.matrix(input_dataframe))
    row_missing_proportion <- rowMeans(is.na(numeric_matrix))
    rows_to_keep <- row_missing_proportion <= max_na_proportion
    input_dataframe <- input_dataframe[rows_to_keep, , drop = FALSE]
    if (!impute_missing_values) input_dataframe <- input_dataframe[stats::complete.cases(input_dataframe), , drop = FALSE]
    return(input_dataframe)
  }
  background.numeric.species1 <- drop.na.rows(background.numeric.species1, max.NA.prop, impute.NA.median)
  background.numeric.species2 <- drop.na.rows(background.numeric.species2, max.NA.prop, impute.NA.median)

  # Check background
  if (ncol(background.numeric.species1) == 0 || ncol(background.numeric.species2) == 0) stop("No numeric background variables remain after excluding columns")

  # Shared variable filtering
  variable.names.shared <- intersect(names(background.numeric.species1), names(background.numeric.species2))
  n.start <- length(variable.names.shared)
  if (n.start == 0) stop("No shared numeric variables found between species backgrounds after excluding columns")

  background.numeric.species1 <- background.numeric.species1[variable.names.shared]
  background.numeric.species2 <- background.numeric.species2[variable.names.shared]

  # Remove variables with insufficient data rows
  variables.with.sufficient.rows <- variable.names.shared[sapply(variable.names.shared, function(variable_name) {
    sum(!is.na(background.numeric.species1[[variable_name]])) >= min.rows &&
      sum(!is.na(background.numeric.species2[[variable_name]])) >= min.rows
  })]
  n.removed.rows <- length(variable.names.shared) - length(variables.with.sufficient.rows)
  if (n.removed.rows > 0 && verbose) message("Removed ", n.removed.rows, " of ", n.start, " variables due to too few non-NA rows (min.rows = ", min.rows, ")")
  variable.names.shared <- variables.with.sufficient.rows
  if (length(variable.names.shared) == 0) stop("No variables remain after filtering for minimum non-NA rows")
  background.numeric.species1 <- background.numeric.species1[variable.names.shared]
  background.numeric.species2 <- background.numeric.species2[variable.names.shared]

  # Filter low-variation variables by coefficient of variation
  calculate.CV <- function(numeric_vector, variable_name = NA_character_) {
    finite_values <- numeric_vector[is.finite(numeric_vector)]
    if (!length(finite_values)) return(NA_real_)
    mean_value <- mean(finite_values)
    sd_value <- stats::sd(finite_values)
    if (!is.finite(sd_value) || sd_value == 0) return(0)
    if (!is.finite(mean_value) || abs(mean_value) < 1e-7) {
      if (verbose) message("Variable '", variable_name, "' has mean near zero (", format(mean_value, digits = 4), ") - falling back to SD/MAD variability")
      mad_val <- stats::mad(finite_values, constant = 1, na.rm = TRUE)
      if (!is.finite(mad_val) || mad_val == 0) return(0)
      return(sd_value / (mad_val + 1e-12))
    }
    abs(sd_value / mean_value)
  }
  cv_species1 <- sapply(variable.names.shared, function(var) calculate.CV(background.numeric.species1[[var]], var))
  cv_species2 <- sapply(variable.names.shared, function(var) calculate.CV(background.numeric.species2[[var]], var))
  variables.low.cv <- variable.names.shared[is.na(cv_species1) | cv_species1 <= CV.threshold | is.na(cv_species2) | cv_species2 <= CV.threshold]
  n.removed.cv <- length(variables.low.cv)
  n.prev <- length(variable.names.shared)
  if (n.removed.cv > 0 && verbose) message("Removed ", n.removed.cv, " of ", n.prev, " variables due to low variation (CV.threshold <= ", CV.threshold, ")")
  variable.names.shared <- setdiff(variable.names.shared, variables.low.cv)
  if (length(variable.names.shared) == 0) stop("No variables remain after filtering for coefficient of variation")
  background.numeric.species1 <- background.numeric.species1[variable.names.shared]
  background.numeric.species2 <- background.numeric.species2[variable.names.shared]

  # Skip overlap filtering if overlap.threshold = 0
  if (overlap.threshold == 0) {
    if (verbose) message("Skipping all overlap filtering (overlap.threshold = 0)")
    cols.keep <- union(exclude.cols, variable.names.shared)
    output.data <- Sp1.Sp2.occurrence.data %>% dplyr::select(dplyr::any_of(cols.keep))
    return(output.data)
  }

  # Univariate KDE overlap
  compute.univariate.kde.overlap <- function(values_species1, values_species2) {
    na_prop_1 <- mean(is.na(values_species1))
    na_prop_2 <- mean(is.na(values_species2))
    if (na_prop_1 > max.NA.prop || na_prop_2 > max.NA.prop) return(NA_real_)
    if (impute.NA.median) {
      values_species1[is.na(values_species1)] <- median(values_species1, na.rm = TRUE)
      values_species2[is.na(values_species2)] <- median(values_species2, na.rm = TRUE)
    } else {
      values_species1 <- values_species1[!is.na(values_species1)]
      values_species2 <- values_species2[!is.na(values_species2)]
    }
    range_values <- range(c(values_species1, values_species2), na.rm = TRUE)
    if (!all(is.finite(range_values))) return(NA_real_)
    if (diff(range_values) == 0) return(1)
    pooled_values <- c(values_species1, values_species2)
    shared_bandwidth <- tryCatch(stats::bw.SJ(pooled_values, method = "dpi"), error = function(e) NA_real_)
    if (!is.finite(shared_bandwidth) || shared_bandwidth <= 0) shared_bandwidth <- stats::bw.nrd0(pooled_values)
    density_species1 <- stats::density(values_species1, from = range_values[1], to = range_values[2], n = 768, bw = shared_bandwidth)
    density_species2 <- stats::density(values_species2, from = range_values[1], to = range_values[2], n = 768, bw = shared_bandwidth)
    step_size <- diff(density_species1$x[1:2])
    prob_species1 <- density_species1$y / sum(density_species1$y * step_size)
    prob_species2 <- density_species2$y / sum(density_species2$y * step_size)
    sum(pmin(prob_species1, prob_species2) * step_size)
  }
  univariate.kde.overlap.values <- sapply(variable.names.shared, function(var) compute.univariate.kde.overlap(background.numeric.species1[[var]], background.numeric.species2[[var]]))
  if (plot.1D.overlap) {
    hist_values <- univariate.kde.overlap.values[is.finite(univariate.kde.overlap.values)]
    if (length(hist_values)) {
      par(mfrow = c(1, 1))
      hist(hist_values, main = "Histogram of 1D KDE overlap values", xlab = "Overlap (Schoener's D)", breaks = 20, xlim = c(0, 1))
      abline(v = overlap.threshold, col = 2, lwd = 2)
    } else if (verbose) message("No finite 1D overlap values to plot")
  }
  variables.removed.univariate <- names(univariate.kde.overlap.values)[is.na(univariate.kde.overlap.values) | univariate.kde.overlap.values < overlap.threshold]
  n.removed.uni <- length(variables.removed.univariate)
  n.prev <- length(variable.names.shared)
  if (n.removed.uni > 0 && verbose) message("Removed ", n.removed.uni, " of ", n.prev, " variables due to univariate overlap (overlap.threshold < ", overlap.threshold, ")")
  variable.names.shared <- setdiff(variable.names.shared, variables.removed.univariate)
  if (length(variable.names.shared) == 0) stop("No variables remain after univariate overlap filtering")
  background.numeric.species1 <- background.numeric.species1[variable.names.shared]
  background.numeric.species2 <- background.numeric.species2[variable.names.shared]

  # Bivariate overlap time
  if (is.null(bin.n.2D)) {
    effective_sample_size <- min(nrow(background.numeric.species1), nrow(background.numeric.species2))
    if (!is.finite(effective_sample_size) || effective_sample_size <= 0) effective_sample_size <- min.rows * 10
    bin.n.2D <- max(5, min(20, floor(sqrt(effective_sample_size / 20))))
  }
  bin.n.2D <- as.integer(bin.n.2D[1])
  if (!is.finite(bin.n.2D) || bin.n.2D < 3) stop("bin.n.2D must be a single integer >= 3")
  p_shared <- length(variable.names.shared)
  feasible_min_per_var <- floor((2L * max.pairs) / p_shared)
  auto_min_pairs_per_var <- max(1L, min(feasible_min_per_var, p_shared - 1L))
  total_pairs <- choose(p_shared, 2L)
  stratify.sample.pairs <- function(variable_names, max_pairs_total, min_pairs_per_variable, rng_seed = 1L) {
    variable_names <- as.character(variable_names)
    variable_count <- length(variable_names)
    if (variable_count < 2) return(list(pairs = list(), subsampled = FALSE, min_pairs_achieved = 0L))
    pair_matrix <- utils::combn(variable_names, 2)
    total_pair_count <- ncol(pair_matrix)
    if (total_pair_count <= max_pairs_total) {
      all_pairs <- lapply(seq_len(total_pair_count), function(col_idx) pair_matrix[, col_idx])
      min_pairs_achieved <- min(min_pairs_per_variable, variable_count - 1L)
      return(list(pairs = all_pairs, subsampled = FALSE, min_pairs_achieved = min_pairs_achieved))
    }
    feasible_min_pairs_per_variable <- floor((2L * max_pairs_total) / variable_count)
    min_pairs_achieved <- max(0L, min(min_pairs_per_variable, feasible_min_pairs_per_variable, variable_count - 1L))
    if (!is.null(rng_seed)) set.seed(rng_seed)
    if (min_pairs_achieved == 0L) {
      sample_size <- min(max_pairs_total, total_pair_count)
      chosen_indices <- sample.int(total_pair_count, sample_size)
      sampled_pairs <- lapply(chosen_indices, function(col_idx) pair_matrix[, col_idx])
      all_vars <- variable_names
      vars_in_pairs <- unique(unlist(sampled_pairs))
      missing_vars <- setdiff(all_vars, vars_in_pairs)
      if (length(missing_vars) > 0) {
        supplemental_pairs <- lapply(missing_vars, function(v) {
          partner <- sample(setdiff(all_vars, v), 1)
          c(v, partner)
        })
        sampled_pairs <- c(sampled_pairs, supplemental_pairs)
      }
      return(list(pairs = sampled_pairs, subsampled = (sample_size < total_pair_count), min_pairs_achieved = 0L))
    }
    selected_pair_flags <- rep(FALSE, total_pair_count)
    pair_count_per_variable <- setNames(integer(variable_count), variable_names)
    adjacency_indices_by_variable <- lapply(variable_names, function(var_name) {
      which(pair_matrix[1, ] == var_name | pair_matrix[2, ] == var_name)
    })
    names(adjacency_indices_by_variable) <- variable_names
    selection_limit <- min(max_pairs_total, total_pair_count)
    while (sum(selected_pair_flags) < selection_limit) {
      underrepresented_variables <- names(pair_count_per_variable)[pair_count_per_variable < min_pairs_achieved]
      if (!length(underrepresented_variables)) break
      coverage_deficit <- pmax(0L, min_pairs_achieved - pair_count_per_variable[underrepresented_variables])
      most_undercovered_variable <- underrepresented_variables[which.max(coverage_deficit)]
      if (is.na(most_undercovered_variable) || !(most_undercovered_variable %in% names(adjacency_indices_by_variable))) break
      candidate_pair_indices <- adjacency_indices_by_variable[[most_undercovered_variable]]
      candidate_pair_indices <- candidate_pair_indices[!selected_pair_flags[candidate_pair_indices]]
      if (!length(candidate_pair_indices)) {
        pair_count_per_variable[most_undercovered_variable] <- min_pairs_achieved
        next
      }
      deficit_reduction <- vapply(candidate_pair_indices, function(col_idx) {
        left_var <- pair_matrix[1, col_idx]
        right_var <- pair_matrix[2, col_idx]
        max(0L, min_pairs_achieved - pair_count_per_variable[left_var]) + max(0L, min_pairs_achieved - pair_count_per_variable[right_var])
      }, numeric(1))
      best_pair_indices <- candidate_pair_indices[deficit_reduction == max(deficit_reduction)]
      chosen_col_idx <- best_pair_indices[sample.int(length(best_pair_indices), 1L)]
      selected_pair_flags[chosen_col_idx] <- TRUE
      chosen_pair <- pair_matrix[, chosen_col_idx]
      pair_count_per_variable[chosen_pair] <- pair_count_per_variable[chosen_pair] + 1L
    }
    if (sum(selected_pair_flags) < selection_limit) {
      remaining_indices <- which(!selected_pair_flags)
      fill_count <- min(length(remaining_indices), selection_limit - sum(selected_pair_flags))
      if (fill_count > 0) {
        add_indices <- sample(remaining_indices, fill_count)
        selected_pair_flags[add_indices] <- TRUE
      }
    }
    final_indices <- which(selected_pair_flags)
    sampled_pairs <- lapply(final_indices, function(col_idx) pair_matrix[, col_idx])
    list(pairs = sampled_pairs, subsampled = (length(final_indices) < total_pair_count), min_pairs_achieved = min_pairs_achieved)
  }
  pair_sample <- stratify.sample.pairs(variable.names.shared, max.pairs, auto_min_pairs_per_var, rng_seed = seed)
  sampled_pairs <- pair_sample$pairs
  if (!length(sampled_pairs)) stop("Fewer than two variables remain")

  # Compute 2D histogram-overlap (bivariate environmental similarity)
  compute.bivariate.fast.overlap <- function(pair_vars) {
    vec1a <- background.numeric.species1[[pair_vars[1]]]
    vec1b <- background.numeric.species1[[pair_vars[2]]]
    vec2a <- background.numeric.species2[[pair_vars[1]]]
    vec2b <- background.numeric.species2[[pair_vars[2]]]
    if (mean(is.na(vec1a)) > max.NA.prop || mean(is.na(vec1b)) > max.NA.prop || mean(is.na(vec2a)) > max.NA.prop || mean(is.na(vec2b)) > max.NA.prop) return(NA_real_)
    if (impute.NA.median) {
      vec1a[is.na(vec1a)] <- median(vec1a, na.rm = TRUE)
      vec1b[is.na(vec1b)] <- median(vec1b, na.rm = TRUE)
      vec2a[is.na(vec2a)] <- median(vec2a, na.rm = TRUE)
      vec2b[is.na(vec2b)] <- median(vec2b, na.rm = TRUE)
    } else {
      ok1 <- complete.cases(vec1a, vec1b)
      vec1a <- vec1a[ok1]
      vec1b <- vec1b[ok1]
      ok2 <- complete.cases(vec2a, vec2b)
      vec2a <- vec2a[ok2]
      vec2b <- vec2b[ok2]
    }
    if (length(vec1a) < min.rows || length(vec2a) < min.rows) return(NA_real_)
    total_bins <- bin.n.2D
    pooled_x <- c(vec1a, vec2a)
    pooled_y <- c(vec1b, vec2b)
    if (!any(is.finite(pooled_x)) || !any(is.finite(pooled_y))) return(NA_real_)
    bin_boundary_epsilon <- .Machine$double.eps^0.5
    quantiles_x <- stats::quantile(pooled_x, probs = seq(0, 1, length.out = total_bins + 1), na.rm = TRUE, names = FALSE)
    quantiles_y <- stats::quantile(pooled_y, probs = seq(0, 1, length.out = total_bins + 1), na.rm = TRUE, names = FALSE)
    quantiles_x[c(1, length(quantiles_x))] <- quantiles_x[c(1, length(quantiles_x))] + c(-bin_boundary_epsilon, +bin_boundary_epsilon)
    quantiles_y[c(1, length(quantiles_y))] <- quantiles_y[c(1, length(quantiles_y))] + c(-bin_boundary_epsilon, +bin_boundary_epsilon)
    breaks_x <- unique(quantiles_x)
    breaks_y <- unique(quantiles_y)
    n_bins_x <- length(breaks_x) - 1L
    n_bins_y <- length(breaks_y) - 1L
    if (n_bins_x < 1L || n_bins_y < 1L) return(NA_real_)
    if (length(breaks_x) < 3 || length(breaks_y) < 3) return(NA_real_)
    ix1 <- cut(vec1a, breaks = breaks_x, include.lowest = TRUE, labels = FALSE)
    iy1 <- cut(vec1b, breaks = breaks_y, include.lowest = TRUE, labels = FALSE)
    ix2 <- cut(vec2a, breaks = breaks_x, include.lowest = TRUE, labels = FALSE)
    iy2 <- cut(vec2b, breaks = breaks_y, include.lowest = TRUE, labels = FALSE)
    tabA <- table(factor(ix1, levels = seq_len(n_bins_x)), factor(iy1, levels = seq_len(n_bins_y)))
    tabB <- table(factor(ix2, levels = seq_len(n_bins_x)), factor(iy2, levels = seq_len(n_bins_y)))
    tabA <- as.matrix(tabA)
    tabB <- as.matrix(tabB)
    totalA <- sum(tabA)
    totalB <- sum(tabB)
    if (totalA == 0 || totalB == 0) return(NA_real_)
    probA <- tabA / totalA
    probB <- tabB / totalB
    sum(pmin(probA, probB))
  }
  if (use.parallel) {
    avail <- tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
    if (is.finite(avail) && N.cores > avail && verbose) warning("Requested N.cores (", N.cores, ") exceeds available cores (", avail, ").")
    cl <- parallel::makeCluster(N.cores)
    if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, varlist = c("background.numeric.species1", "background.numeric.species2", "impute.NA.median", "min.rows", "bin.n.2D", "max.NA.prop", "compute.bivariate.fast.overlap"), envir = environment())
    bivariate.values <- parallel::parSapply(cl, sampled_pairs, compute.bivariate.fast.overlap)
  } else {
    if (!is.null(seed)) set.seed(seed)
    bivariate.values <- sapply(sampled_pairs, compute.bivariate.fast.overlap)
  }
  keep.vars <- unique(unlist(sampled_pairs[!is.na(bivariate.values) & bivariate.values >= overlap.threshold]))
  n.removed.biv <- length(variable.names.shared) - length(keep.vars)
  if (n.removed.biv > 0 && verbose) message("Removed ", n.removed.biv, " of ", length(variable.names.shared), " variables due to bivariate overlap (overlap.threshold < ", overlap.threshold, ")")

  # Print summary
  n.retained <- length(keep.vars)
  if (verbose) {
    message("")
    message("Filtering finished:")
    message("Retained ", n.retained, " of ", n.start, " variables after all filtering steps")
  }
  if (n.retained == 0) stop("No variables remain after 2D bin overlap filtering")

  # Final output assembly
  cols.keep <- union(exclude.cols, keep.vars)
  output.data <- Sp1.Sp2.occurrence.data %>% dplyr::select(dplyr::any_of(cols.keep))
  if (nrow(output.data) == 0) stop("No occurrence rows remain after filtering")

  # Return results
  return(output.data)
}


## Function to trim occurrences and backgrounds to shared analogous environmental space
#' Trim occurrence and background data to analogous environmental space
#'
#' Restrict species occurrence and background datasets to the portion of
#' environmental space shared between both species, following the logic of
#' analogous-environment trimming described by Brown and Carnaval (2019) and the
#' g2e workflow implemented in the Humboldt R package. The function compares
#' the environmental backgrounds of two species, identifies the shared analogous
#' portion of environmental space, and retains only occurrence records falling
#' within that shared space.
#'
#' @param Sp1.occurrence.data A `data.frame` or `sf` object containing occurrence
#'   records for species 1.
#' @param Sp2.occurrence.data A `data.frame` or `sf` object containing occurrence
#'   records for species 2.
#' @param Sp1.background.data A `data.frame` or `sf` object containing background
#'   environmental data for species 1.
#' @param Sp2.background.data A `data.frame` or `sf` object containing background
#'   environmental data for species 2.
#' @param exclude.cols Optional character vector of column names to exclude from
#'   analogous-space estimation (default: `NULL`).
#' @param keep.occurrence.cols Optional character vector of occurrence-data
#'   columns that should always be retained and placed first in the output
#'   (default: `NULL`).
#' @param analogous.window.size A single non-negative numeric value controlling
#'   the width of the analogous-environment retention window in environmental
#'   PCA space (default: `5`).
#' @param grid.resolution A single integer-like numeric value giving the number
#'   of grid divisions per PCA axis used to estimate shared environmental space
#'   (default: `50`).
#' @param max.NA.prop A single numeric value between 0 and 1 giving the maximum
#'   allowed proportion of missing values per background row before that row is
#'   removed from analogous-space estimation (default: `0.2`).
#' @param impute.NA.median Logical; if `TRUE`, eligible missing values are
#'   imputed using the variable median before estimating analogous environmental
#'   space. If `FALSE`, incomplete rows are removed instead (default: `TRUE`).
#' @param downsample.equal.sizes Logical; if `TRUE`, the two trimmed occurrence
#'   datasets are downsampled to equal sample sizes after trimming
#'   (default: `TRUE`).
#' @param verbose Logical; if `TRUE`, progress messages and summary information
#'   are printed during trimming (default: `TRUE`).
#'
#' @details
#' Humboldt-style analogous trimming is used as a complementary correction to
#' predictor-based environmental analogy screening. Whereas variable filtering
#' removes environmental predictors that are poorly comparable between species,
#' occurrence trimming removes occurrence records that fall outside the
#' environmentally shared region of the two species' accessible background
#' spaces. This helps reduce bias caused by non-analogous environments, which can
#' inflate apparent niche divergence and make comparisons depend on extrapolated
#' portions of environmental space rather than on conditions available to both
#' species (Brown & Carnaval, 2019; Guisan et al., 2014; Peterson et al., 2011).
#'
#' The method summarizes the combined environmental backgrounds of both species
#' in a shared principal component space. Using a common ordination space ensures
#' that both species are compared along the same environmental axes rather than
#' in separate species-specific coordinate systems. Principal component analysis
#' is appropriate here because it provides a reduced environmental space that
#' captures dominant gradients of covariation among predictors while allowing
#' occurrence and background points to be evaluated on the same scale (Chessel et
#' al., 2004; Dray & Dufour, 2007; Bougeard & Dray, 2018; Thioulouse et al.,
#' 2018).
#'
#' Prior to trimming, low-information predictors and rows with excessive
#' missingness are removed because near-constant variables and heavily incomplete
#' observations can distort environmental-space estimates. The missingness
#' threshold (`max.NA.prop = 0.2`) is intended as a practical default that
#' permits moderate missingness while reducing the risk that the shared
#' environmental space is driven by imputation or incomplete background records
#' (Harrell, 2015; van Buuren, 2018).
#'
#' The analogous trimming window controls how strictly environmental space is
#' shared between species. Smaller windows impose stricter analogy requirements
#' and remove more records, whereas larger windows retain more records but may
#' allow a broader range of marginally comparable environments. The default
#' window size (`analogous.window.size = 5`) is intended as a moderate compromise
#' that reduces extrapolation into non-analogous environmental regions while
#' retaining enough occurrence records for downstream analyses.
#'
#' Reciprocal trimming is used so that neither species defines the analogous
#' space unilaterally. By restricting each species relative to the environmental
#' space occupied by the other, the retained region represents a shared subset of
#' background environmental space rather than the environmental availability of
#' only one species. This is important for niche-divergence analyses because
#' asymmetric trimming can bias comparisons toward the species with the broader
#' or more densely sampled background environment (Brown & Carnaval, 2019).
#'
#' After the shared analogous region is identified, occurrence records are
#' retained only when they fall inside that shared environmental space. This
#' produces occurrence datasets restricted to comparable environmental
#' conditions, reducing the influence of extrapolated or non-analogous
#' occurrences that could otherwise bias downstream comparisons of niche
#' divergence, niche overlap, or discriminant environmental separation.
#'
#' Equalizing sample sizes after trimming is useful when downstream
#' methods are sensitive to unequal occurrence counts. Downsampling to matched
#' sample sizes can reduce imbalance between species and make subsequent
#' comparisons less dependent on differences in retained occurrence density.
#'
#' @return A `data.frame` combining the retained occurrence rows for both species
#'   after Humboldt-style trimming to analogous environmental space. If
#'   `downsample.equal.sizes = TRUE`, the returned species subsets are
#'   downsampled to equal sample sizes.
#'
#' @references
#' Bougeard, S., & Dray, S. (2018). Supervised multiblock analysis in R with the
#'   ade4 package. \emph{Journal of Statistical Software}, 86(1).
#'   https://doi.org/10.18637/jss.v086.i01
#'
#' Brown, J., & Carnaval, A. C. (2019). A tale of two niches: Methods, concepts,
#'   and evolution. \emph{Frontiers of Biogeography}, 11(4).
#'   https://doi.org/10.21425/F5FBG44158
#'
#' Chessel, D., Dufour, A., & Thioulouse, J. (2004). The ade4 Package - I:
#'   One-table methods. \emph{R News}, 4(1), 5-10.
#'
#' Dray, S., & Dufour, A.-B. (2007). The ade4 package: Implementing the duality
#'   diagram for ecologists. \emph{Journal of Statistical Software}, 22(4).
#'   https://doi.org/10.18637/jss.v022.i04
#'
#' Guisan, A., Petitpierre, B., Broennimann, O., Daehler, C., & Kueffer, C.
#'   (2014). Unifying niche shift studies. Insights from biological invasions.
#'   \emph{Trends in Ecology & Evolution}, 29(5), 260-269.
#'   https://doi.org/10.1016/j.tree.2014.02.009
#'
#' Harrell, F. E. (2015). \emph{Regression modeling strategies: With applications
#'   to linear models, logistic and ordinal regression, and survival analysis}
#'   (2nd ed.). Springer.
#'
#' Peterson, A. T., Soberón, J., Pearson, R. G., Anderson, R. P.,
#'   Martínez-Meyer, E., Nakamura, M., & Araújo, M. B. (2011). \emph{Ecological
#'   niches and geographic distributions}. Princeton University Press.
#'   https://doi.org/10.1515/9781400840670
#'
#' Thioulouse, J., Dray, S., Dufour, A.-B., Siberchicot, A., Jombart, T., &
#'   Pavoine, S. (2018). \emph{Multivariate analysis of ecological data with
#'   ade4}. Springer New York. https://doi.org/10.1007/978-1-4939-8850-1
#'
#' van Buuren, S. (2018). \emph{Flexible imputation of missing data}. Chapman &
#'   Hall/CRC.
#'
#' @export
trim.to.analogous.environments <- function(Sp1.occurrence.data, #occurrence data for species 1
                                           Sp2.occurrence.data, #occurrence data for species 2
                                           Sp1.background.data, #background data for species 1
                                           Sp2.background.data, #background data for species 2
                                           exclude.cols = NULL, #columns to exclude
                                           keep.occurrence.cols = NULL, #columns to retain
                                           analogous.window.size = 5, #width (in PCA grid bands) of moving window used for Humboldt G2E trimming
                                           grid.resolution = 50, #set number of grid divisions per PCA axis used to build the moving-window sequence (Humboldt G2E R parameter)
                                           max.NA.prop = 0.2, #max proportion of NA allowed per row
                                           impute.NA.median = TRUE, #impute median for NA if <= max.NA.prop
                                           downsample.equal.sizes = TRUE, #downsample to equal sample sizes
                                           verbose = TRUE #print messages
) {

  # Validate input
  if (!is.data.frame(Sp1.occurrence.data) || !is.data.frame(Sp2.occurrence.data)) stop("Sp1.occurrence.data and Sp2.occurrence.data must be data frames")
  if (!is.data.frame(Sp1.background.data) || !is.data.frame(Sp2.background.data)) stop("Sp1.background.data and Sp2.background.data must be data frames")
  if (!is.null(exclude.cols) && !is.character(exclude.cols)) stop("exclude.cols must be a character vector or NULL")
  if (!is.null(keep.occurrence.cols) && !is.character(keep.occurrence.cols)) stop("keep.occurrence.cols must be a character vector or NULL")
  if (is.null(exclude.cols)) exclude.cols <- character(0)
  if (is.null(keep.occurrence.cols)) keep.occurrence.cols <- character(0)
  if (!is.numeric(grid.resolution) || length(grid.resolution) != 1 || grid.resolution < 10) stop("grid.resolution must be >= 10")
  if (!is.numeric(analogous.window.size) || length(analogous.window.size) != 1 || analogous.window.size < 0) stop("analogous.window.size must be non-negative")
  if (!is.numeric(max.NA.prop) || length(max.NA.prop) != 1 || max.NA.prop < 0 || max.NA.prop > 1) stop("max.NA.prop must be between 0 and 1")
  if (!is.logical(impute.NA.median) || length(impute.NA.median) != 1) stop("impute.NA.median must be TRUE or FALSE")
  if (!is.logical(downsample.equal.sizes) || length(downsample.equal.sizes) != 1) stop("downsample.equal.sizes must be TRUE or FALSE")
  if (!is.logical(verbose) || length(verbose) != 1) stop("verbose must be TRUE or FALSE")

  # Warn if initial occurrence sample sizes differ
  initial_species1_count <- nrow(Sp1.occurrence.data)
  initial_species2_count <- nrow(Sp2.occurrence.data)
  if (initial_species1_count != initial_species2_count) warning("Initial sample sizes differ between species (species 1 = ", initial_species1_count, ", species 2 = ", initial_species2_count, ") - this may affect trimming balance")

  # Print starting message
  if (verbose) message("Trimming to shared analogous environment (Humboldt-like correction)")

  # Prepare background data
  species1_background_numeric <- Sp1.background.data[, setdiff(names(Sp1.background.data), exclude.cols), drop = FALSE]
  species2_background_numeric <- Sp2.background.data[, setdiff(names(Sp2.background.data), exclude.cols), drop = FALSE]
  species1_background_numeric <- species1_background_numeric[, sapply(species1_background_numeric, is.numeric), drop = FALSE]
  species2_background_numeric <- species2_background_numeric[, sapply(species2_background_numeric, is.numeric), drop = FALSE]
  shared_environment_variable_names <- intersect(names(species1_background_numeric), names(species2_background_numeric))
  if (length(shared_environment_variable_names) < 2) stop("At least two shared numeric environmental variables are required")
  species1_background_numeric <- species1_background_numeric[, shared_environment_variable_names, drop = FALSE]
  species2_background_numeric <- species2_background_numeric[, shared_environment_variable_names, drop = FALSE]
  species1_background_numeric[!is.finite(as.matrix(species1_background_numeric))] <- NA
  species2_background_numeric[!is.finite(as.matrix(species2_background_numeric))] <- NA

  # Combine backgrounds and drop NA/low-variance variables
  combined_background_matrix <- rbind(species1_background_numeric, species2_background_numeric)
  background_variable_all_na_mask <- vapply(combined_background_matrix, function(x) all(is.na(x)), logical(1))
  background_variable_low_variance_mask <- vapply(combined_background_matrix, function(x) {
    sd_x <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(sd_x) || sd_x < 1e-9) return(TRUE)
    return(FALSE)
  }, logical(1))
  background_variable_drop_mask <- background_variable_all_na_mask | background_variable_low_variance_mask
  if (any(background_variable_drop_mask)) {
    if (verbose && any(background_variable_all_na_mask)) message("Dropping background variables with all NA: ", paste(names(combined_background_matrix)[background_variable_all_na_mask], collapse = ", "))
    if (verbose && any(background_variable_low_variance_mask & !background_variable_all_na_mask)) message("Dropping background variables with near-zero variance: ", paste(names(combined_background_matrix)[background_variable_low_variance_mask & !background_variable_all_na_mask], collapse = ", "))
    combined_background_matrix <- combined_background_matrix[, !background_variable_drop_mask, drop = FALSE]
    species1_background_numeric <- species1_background_numeric[, !background_variable_drop_mask, drop = FALSE]
    species2_background_numeric <- species2_background_numeric[, !background_variable_drop_mask, drop = FALSE]
  }
  if (ncol(combined_background_matrix) < 2) stop("Fewer than two valid predictors remain after filtering variables")

  # Filter background rows based on NA proportion
  compute_row_drop_mask <- function(df) {
    apply(df, 1, function(row_values) {
      na.proportion <- mean(is.na(row_values))
      if (na.proportion > max.NA.prop) return(TRUE)
      if (all(is.na(row_values))) return(TRUE)
      return(FALSE)
    })
  }
  species1_background_row_drop_mask <- compute_row_drop_mask(species1_background_numeric)
  species2_background_row_drop_mask <- compute_row_drop_mask(species2_background_numeric)
  if (verbose && any(species1_background_row_drop_mask)) message("Dropping ", sum(species1_background_row_drop_mask), " background rows for species 1 due to NA proportion")
  if (verbose && any(species2_background_row_drop_mask)) message("Dropping ", sum(species2_background_row_drop_mask), " background rows for species 2 due to NA proportion")
  species1_background_numeric <- species1_background_numeric[!species1_background_row_drop_mask, , drop = FALSE]
  species2_background_numeric <- species2_background_numeric[!species2_background_row_drop_mask, , drop = FALSE]
  if (nrow(species1_background_numeric) < 2 || nrow(species2_background_numeric) < 2) stop("Fewer than two valid background rows remain after filtering")
  if (!impute.NA.median) {
    complete1 <- stats::complete.cases(species1_background_numeric)
    complete2 <- stats::complete.cases(species2_background_numeric)
    if (verbose && any(!complete1)) message("Dropping ", sum(!complete1), " background rows for species 1 with remaining NA")
    if (verbose && any(!complete2)) message("Dropping ", sum(!complete2), " background rows for species 2 with remaining NA")
    species1_background_numeric <- species1_background_numeric[complete1, , drop = FALSE]
    species2_background_numeric <- species2_background_numeric[complete2, , drop = FALSE]
    if (nrow(species1_background_numeric) < 2 || nrow(species2_background_numeric) < 2) stop("Fewer than two valid background rows remain after removing NA rows for PCA")
  }

  # Recombine backgrounds after row filtering
  combined_background_matrix <- rbind(species1_background_numeric, species2_background_numeric)
  if (ncol(combined_background_matrix) < 2) stop("Fewer than two valid predictors remain for PCA")

  # Impute background NA
  if (impute.NA.median) {
    for (variable_index in seq_len(ncol(combined_background_matrix))) {
      na_row_mask <- is.na(combined_background_matrix[, variable_index])
      if (any(na_row_mask)) combined_background_matrix[na_row_mask, variable_index] <- median(combined_background_matrix[, variable_index], na.rm = TRUE)
    }
  }
  if (!all(is.finite(as.matrix(combined_background_matrix)))) stop("Non-finite background values remain after imputation")

  # Run PCA (E-space)
  pca_model <- ade4::dudi.pca(combined_background_matrix, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)
  species1_background_row_count <- nrow(species1_background_numeric)
  species2_background_row_count <- nrow(species2_background_numeric)
  species1_background_pca_scores <- pca_model$li[1:species1_background_row_count, 1:2, drop = FALSE]
  species2_background_pca_scores <- pca_model$li[(species1_background_row_count + 1):(species1_background_row_count + species2_background_row_count), 1:2, drop = FALSE]

  # Initialize current scores for REDUC = 5 trimming
  species1_scores_current <- species1_background_pca_scores
  species2_scores_current <- species2_background_pca_scores

  # Create function to check that both species still have environments
  check_env_nonempty <- function(scores1, scores2) {
    if (nrow(scores1) == 0 || nrow(scores2) == 0) stop("No analogous environments remain for at least one species after trimming (Humboldt-style REDUC loop)")
  }

  # Perform full REDUC = 5 loops: PC1 trim S2 -> PC2 trim S2 -> PC2 trim S1 -> PC1 trim S1
  reduction_iterations <- 5
  for (reduction_iteration in seq_len(reduction_iterations)) {

    # Trim species 2 along PC1 based on species 1
    if (nrow(species1_scores_current) < 2 || nrow(species2_scores_current) < 2) check_env_nonempty(species1_scores_current, species2_scores_current)
    pc1_max_species1 <- max(species1_scores_current[, 1])
    pc1_min_species1 <- min(species1_scores_current[, 1])
    pc1_range_length <- pc1_max_species1 - pc1_min_species1
    if (pc1_range_length <= .Machine$double.eps) stop("Species 1 PC1 range collapsed during trimming")
    pc1_bandwidth <- pc1_range_length / grid.resolution
    pc1_sequence <- seq(pc1_min_species1, pc1_max_species1, length.out = grid.resolution)
    species2_trimmed_pc1 <- NULL
    for (pc1_window_center in pc1_sequence) {
      species1_band_subset <- species1_scores_current[species1_scores_current[, 1] <= (pc1_window_center + (analogous.window.size + 1) * pc1_bandwidth) & species1_scores_current[, 1] >= (pc1_window_center - analogous.window.size * pc1_bandwidth), , drop = FALSE]
      if (nrow(species1_band_subset) == 0) next
      pc2_max_for_window <- max(species1_band_subset[, 2])
      pc2_min_for_window <- min(species1_band_subset[, 2])
      species2_band_subset <- species2_scores_current[species2_scores_current[, 1] <= (pc1_window_center + pc1_bandwidth) & species2_scores_current[, 1] >= pc1_window_center, , drop = FALSE]
      if (nrow(species2_band_subset) == 0) next
      species2_band_subset <- species2_band_subset[species2_band_subset[, 2] <= pc2_max_for_window & species2_band_subset[, 2] >= pc2_min_for_window, , drop = FALSE]
      if (nrow(species2_band_subset) == 0) next
      if (is.null(species2_trimmed_pc1)) species2_trimmed_pc1 <- species2_band_subset else species2_trimmed_pc1 <- rbind(species2_trimmed_pc1, species2_band_subset)
    }
    if (is.null(species2_trimmed_pc1) || nrow(species2_trimmed_pc1) == 0) stop("No analogous environments retained for species 2 after PC1 trimming")
    species2_scores_current <- species2_trimmed_pc1

    # Trim species 2 along PC2 based on species 1
    if (nrow(species1_scores_current) < 2 || nrow(species2_scores_current) < 2) check_env_nonempty(species1_scores_current, species2_scores_current)
    pc2_max_species1 <- max(species1_scores_current[, 2])
    pc2_min_species1 <- min(species1_scores_current[, 2])
    pc2_range_length <- pc2_max_species1 - pc2_min_species1
    if (pc2_range_length <= .Machine$double.eps) stop("Species 1 PC2 range collapsed during trimming")
    pc2_bandwidth <- pc2_range_length / grid.resolution
    pc2_sequence <- seq(pc2_min_species1, pc2_max_species1, length.out = grid.resolution)
    species2_trimmed_pc2 <- NULL
    for (pc2_window_center in pc2_sequence) {
      species1_band_subset <- species1_scores_current[species1_scores_current[, 2] <= (pc2_window_center + (analogous.window.size + 1) * pc2_bandwidth) & species1_scores_current[, 2] >= (pc2_window_center - analogous.window.size * pc2_bandwidth), , drop = FALSE]
      if (nrow(species1_band_subset) == 0) next
      pc1_max_for_window <- max(species1_band_subset[, 1])
      pc1_min_for_window <- min(species1_band_subset[, 1])
      species2_band_subset <- species2_scores_current[species2_scores_current[, 2] <= (pc2_window_center + pc2_bandwidth) & species2_scores_current[, 2] >= pc2_window_center, , drop = FALSE]
      if (nrow(species2_band_subset) == 0) next
      species2_band_subset <- species2_band_subset[species2_band_subset[, 1] <= pc1_max_for_window & species2_band_subset[, 1] >= pc1_min_for_window, , drop = FALSE]
      if (nrow(species2_band_subset) == 0) next
      if (is.null(species2_trimmed_pc2)) species2_trimmed_pc2 <- species2_band_subset else species2_trimmed_pc2 <- rbind(species2_trimmed_pc2, species2_band_subset)
    }
    if (is.null(species2_trimmed_pc2) || nrow(species2_trimmed_pc2) == 0) stop("No analogous environments retained for species 2 after PC2 trimming")
    species2_scores_current <- species2_trimmed_pc2

    # Trim species 1 along PC2 based on species 2 (reciprocal)
    if (nrow(species1_scores_current) < 2 || nrow(species2_scores_current) < 2) check_env_nonempty(species1_scores_current, species2_scores_current)
    pc2_max_species2 <- max(species2_scores_current[, 2])
    pc2_min_species2 <- min(species2_scores_current[, 2])
    pc2_range_length_2 <- pc2_max_species2 - pc2_min_species2
    if (pc2_range_length_2 <= .Machine$double.eps) stop("Species 2 PC2 range collapsed during trimming")
    pc2_bandwidth_2 <- pc2_range_length_2 / grid.resolution
    pc2_sequence_2 <- seq(pc2_min_species2, pc2_max_species2, length.out = grid.resolution)
    species1_trimmed_pc2 <- NULL
    for (pc2_window_center in pc2_sequence_2) {
      species2_band_subset <- species2_scores_current[species2_scores_current[, 2] <= (pc2_window_center + (analogous.window.size + 1) * pc2_bandwidth_2) & species2_scores_current[, 2] >= (pc2_window_center - analogous.window.size * pc2_bandwidth_2), , drop = FALSE]
      if (nrow(species2_band_subset) == 0) next
      pc1_max_for_window <- max(species2_band_subset[, 1])
      pc1_min_for_window <- min(species2_band_subset[, 1])
      species1_band_subset <- species1_scores_current[species1_scores_current[, 2] <= (pc2_window_center + pc2_bandwidth_2) & species1_scores_current[, 2] >= pc2_window_center, , drop = FALSE]
      if (nrow(species1_band_subset) == 0) next
      species1_band_subset <- species1_band_subset[species1_band_subset[, 1] <= pc1_max_for_window & species1_band_subset[, 1] >= pc1_min_for_window, , drop = FALSE]
      if (nrow(species1_band_subset) == 0) next
      if (is.null(species1_trimmed_pc2)) species1_trimmed_pc2 <- species1_band_subset else species1_trimmed_pc2 <- rbind(species1_trimmed_pc2, species1_band_subset)
    }
    if (is.null(species1_trimmed_pc2) || nrow(species1_trimmed_pc2) == 0) stop("No analogous environments retained for species 1 after reciprocal PC2 trimming")
    species1_scores_current <- species1_trimmed_pc2

    # Trim species 1 along PC1 based on species 2 (reciprocal)
    if (nrow(species1_scores_current) < 2 || nrow(species2_scores_current) < 2) check_env_nonempty(species1_scores_current, species2_scores_current)
    pc1_max_species2 <- max(species2_scores_current[, 1])
    pc1_min_species2 <- min(species2_scores_current[, 1])
    pc1_range_length_2 <- pc1_max_species2 - pc1_min_species2
    if (pc1_range_length_2 <= .Machine$double.eps) stop("Species 2 PC1 range collapsed during trimming")
    pc1_bandwidth_2 <- pc1_range_length_2 / grid.resolution
    pc1_sequence_2 <- seq(pc1_min_species2, pc1_max_species2, length.out = grid.resolution)
    species1_trimmed_pc1 <- NULL
    for (pc1_window_center in pc1_sequence_2) {
      species2_band_subset <- species2_scores_current[species2_scores_current[, 1] <= (pc1_window_center + (analogous.window.size + 1) * pc1_bandwidth_2) & species2_scores_current[, 1] >= (pc1_window_center - analogous.window.size * pc1_bandwidth_2), , drop = FALSE]
      if (nrow(species2_band_subset) == 0) next
      pc2_max_for_window <- max(species2_band_subset[, 2])
      pc2_min_for_window <- min(species2_band_subset[, 2])
      species1_band_subset <- species1_scores_current[species1_scores_current[, 1] <= (pc1_window_center + pc1_bandwidth_2) & species1_scores_current[, 1] >= pc1_window_center, , drop = FALSE]
      if (nrow(species1_band_subset) == 0) next
      species1_band_subset <- species1_band_subset[species1_band_subset[, 2] <= pc2_max_for_window & species1_band_subset[, 2] >= pc2_min_for_window, , drop = FALSE]
      if (nrow(species1_band_subset) == 0) next
      if (is.null(species1_trimmed_pc1)) species1_trimmed_pc1 <- species1_band_subset else species1_trimmed_pc1 <- rbind(species1_trimmed_pc1, species1_band_subset)
    }
    if (is.null(species1_trimmed_pc1) || nrow(species1_trimmed_pc1) == 0) stop("No analogous environments retained for species 1 after reciprocal PC1 trimming")
    species1_scores_current <- species1_trimmed_pc1
    check_env_nonempty(species1_scores_current, species2_scores_current)
  }

  # Compute final shared PCA bounds after REDUC = 5 loop
  species1_pc1_range <- range(species1_scores_current[, 1])
  species1_pc2_range <- range(species1_scores_current[, 2])
  species2_pc1_range <- range(species2_scores_current[, 1])
  species2_pc2_range <- range(species2_scores_current[, 2])
  if (any(is.na(c(species1_pc1_range, species1_pc2_range, species2_pc1_range, species2_pc2_range)))) stop("Analogous region collapsed: at least one species has no remaining PCA range after trimming and NA filtering")
  shared_pc1_min <- max(species1_pc1_range[1], species2_pc1_range[1])
  shared_pc1_max <- min(species1_pc1_range[2], species2_pc1_range[2])
  shared_pc2_min <- max(species1_pc2_range[1], species2_pc2_range[1])
  shared_pc2_max <- min(species1_pc2_range[2], species2_pc2_range[2])
  if ((shared_pc1_max - shared_pc1_min) < .Machine$double.eps || (shared_pc2_max - shared_pc2_min) < .Machine$double.eps) stop("No shared analogous PCA extent after REDUC = 5 trimming")

  # Prepare occurrences (environmental variables matching PCA)
  pca_environment_variable_names <- names(combined_background_matrix)
  species1_occurrence_environment <- Sp1.occurrence.data[, intersect(pca_environment_variable_names, names(Sp1.occurrence.data)), drop = FALSE]
  species2_occurrence_environment <- Sp2.occurrence.data[, intersect(pca_environment_variable_names, names(Sp2.occurrence.data)), drop = FALSE]
  if (!all(pca_environment_variable_names %in% names(species1_occurrence_environment))) stop("Species 1 missing PCA vars: ", paste(setdiff(pca_environment_variable_names, names(species1_occurrence_environment)), collapse = ", "))
  if (!all(pca_environment_variable_names %in% names(species2_occurrence_environment))) stop("Species 2 missing PCA vars: ", paste(setdiff(pca_environment_variable_names, names(species2_occurrence_environment)), collapse = ", "))
  species1_occurrence_environment <- species1_occurrence_environment[, pca_environment_variable_names, drop = FALSE]
  species2_occurrence_environment <- species2_occurrence_environment[, pca_environment_variable_names, drop = FALSE]
  species1_occurrence_environment[!is.finite(as.matrix(species1_occurrence_environment))] <- NA
  species2_occurrence_environment[!is.finite(as.matrix(species2_occurrence_environment))] <- NA
  if (!impute.NA.median) {
    complete1_occ <- stats::complete.cases(species1_occurrence_environment)
    complete2_occ <- stats::complete.cases(species2_occurrence_environment)
    if (verbose && any(!complete1_occ)) message("Dropping ", sum(!complete1_occ), " occurrence rows for species 1 with remaining NA in PCA variables")
    if (verbose && any(!complete2_occ)) message("Dropping ", sum(!complete2_occ), " occurrence rows for species 2 with remaining NA in PCA variables")
    species1_occurrence_environment <- species1_occurrence_environment[complete1_occ, , drop = FALSE]
    species2_occurrence_environment <- species2_occurrence_environment[complete2_occ, , drop = FALSE]
    Sp1.occurrence.data <- Sp1.occurrence.data[complete1_occ, , drop = FALSE]
    Sp2.occurrence.data <- Sp2.occurrence.data[complete2_occ, , drop = FALSE]
    if (nrow(species1_occurrence_environment) == 0 || nrow(species2_occurrence_environment) == 0) stop("No occurrence rows remain after removing NA rows for projection when impute.NA.median = FALSE")
  }
  species1_occurrence_environment_temp <- species1_occurrence_environment
  species2_occurrence_environment_temp <- species2_occurrence_environment
  if (impute.NA.median) {
    background_medians <- apply(combined_background_matrix[, pca_environment_variable_names, drop = FALSE], 2, median, na.rm = TRUE)
    for (environment_variable in pca_environment_variable_names) {
      missing_mask_1 <- is.na(species1_occurrence_environment_temp[[environment_variable]])
      if (any(missing_mask_1)) species1_occurrence_environment_temp[[environment_variable]][missing_mask_1] <- background_medians[[environment_variable]]
      missing_mask_2 <- is.na(species2_occurrence_environment_temp[[environment_variable]])
      if (any(missing_mask_2)) species2_occurrence_environment_temp[[environment_variable]][missing_mask_2] <- background_medians[[environment_variable]]
    }
  }

  # Project occurrences as supplementary points into PCA
  species1_occurrence_pca_scores <- ade4::suprow(pca_model, species1_occurrence_environment_temp)$li[, 1:2, drop = FALSE]
  species2_occurrence_pca_scores <- ade4::suprow(pca_model, species2_occurrence_environment_temp)$li[, 1:2, drop = FALSE]

  # Identify retained occurrences within shared analogous PCA extent
  species1_occurrence_keep_mask <- species1_occurrence_pca_scores[, 1] >= shared_pc1_min & species1_occurrence_pca_scores[, 1] <= shared_pc1_max & species1_occurrence_pca_scores[, 2] >= shared_pc2_min & species1_occurrence_pca_scores[, 2] <= shared_pc2_max
  species2_occurrence_keep_mask <- species2_occurrence_pca_scores[, 1] >= shared_pc1_min & species2_occurrence_pca_scores[, 1] <= shared_pc1_max & species2_occurrence_pca_scores[, 2] >= shared_pc2_min & species2_occurrence_pca_scores[, 2] <= shared_pc2_max

  species1_occurrence_trimmed <- Sp1.occurrence.data[species1_occurrence_keep_mask, , drop = FALSE]
  species2_occurrence_trimmed <- Sp2.occurrence.data[species2_occurrence_keep_mask, , drop = FALSE]
  species1_occurrence_trimmed_count <- nrow(species1_occurrence_trimmed)
  species2_occurrence_trimmed_count <- nrow(species2_occurrence_trimmed)
  if (species1_occurrence_trimmed_count == 0 || species2_occurrence_trimmed_count == 0) stop("No occurrences retained within shared analogous PCA extent after REDUC = 5 trimming")
  target_sample_size <- min(species1_occurrence_trimmed_count, species2_occurrence_trimmed_count)

  # Downsample occurrences
  if (downsample.equal.sizes) {
    if (verbose) message("Equalizing sample sizes to ", target_sample_size, " (species 1 = ", species1_occurrence_trimmed_count, ", species 2 = ", species2_occurrence_trimmed_count, " originally)")
    species1_occurrence_equalized <- if (species1_occurrence_trimmed_count > target_sample_size)
      sample.down(species1_occurrence_trimmed, N.rows = target_sample_size, poisson.lambda = 1)
    else species1_occurrence_trimmed
    species2_occurrence_equalized <- if (species2_occurrence_trimmed_count > target_sample_size)
      sample.down(species2_occurrence_trimmed, N.rows = target_sample_size, poisson.lambda = 1)
    else species2_occurrence_trimmed
  } else {
    species1_occurrence_equalized <- species1_occurrence_trimmed
    species2_occurrence_equalized <- species2_occurrence_trimmed
  }

  # Reorder occurrence columns
  if (length(keep.occurrence.cols) > 0) {
    front1 <- intersect(keep.occurrence.cols, names(species1_occurrence_equalized))
    rest1 <- setdiff(names(species1_occurrence_equalized), front1)
    species1_occurrence_equalized <- species1_occurrence_equalized[, c(front1, rest1), drop = FALSE]
    front2 <- intersect(keep.occurrence.cols, names(species2_occurrence_equalized))
    rest2 <- setdiff(names(species2_occurrence_equalized), front2)
    species2_occurrence_equalized <- species2_occurrence_equalized[, c(front2, rest2), drop = FALSE]
  }

  # Print summary
  if (verbose) message("Retained ", nrow(species1_occurrence_equalized), " of ", nrow(Sp1.occurrence.data), " samples for species 1 and ", nrow(species2_occurrence_equalized), " of ", nrow(Sp2.occurrence.data), " samples for species 2 after Humboldt-style trimming (REDUC = 5)")

  # Return results
  dplyr::bind_rows(species1_occurrence_equalized, species2_occurrence_equalized)
}


## Function to run DAPC niche divergence test with cross-validation and permutation test
#' Run cross-validated DAPC with permutation testing
#'
#' Fit a DAPC-based niche divergence analysis for two groups and assess
#' significance using permutation testing. The function preprocesses
#' environmental predictors, selects the number of retained PCs by
#' cross-validation (unless a fixed value is supplied), fits the final
#' discriminant model, and evaluates whether observed group separation exceeds
#' null expectations.
#'
#' @param data.input A matrix or `data.frame` containing environmental predictor
#'   variables and the grouping column specified by `species.col`. Non-numeric
#'   predictor columns are removed before analysis. Exactly two groups are
#'   required.
#' @param N.crossval.replicates A single positive integer-like numeric value
#'   giving the number of cross-validation replicates (recommended default: `100`).
#' @param N.permutations A single positive integer-like numeric value giving the
#'   number of permutation replicates (default: `1000`).
#' @param fixed.n.pcs Optional positive integer-like numeric value giving a fixed
#'   number of PCs to retain. If `NULL`, the number of PCs is determined by
#'   cross-validation (default: `NULL`; recommended: `NULL` unless a fixed PC
#'   number is justified by prior analyses or sensitivity tests).
#' @param exclude.cols Optional character vector of columns to remove before
#'   analysis (default: `NULL`).
#' @param species.col A single character string giving the species or group
#'   column name in `data.input` (default: `"Species"`).
#' @param max.NA.prop A single numeric value between `0` and `1` giving the
#'   maximum allowed proportion of missing values per row before that row is
#'   removed (default: `0.2`; recommended: `0.2` to allow moderate missingness
#'   while limiting imputation-driven bias).
#' @param impute.NA.median Logical; if `TRUE`, remaining missing values are
#'   imputed using the median of each variable. If `FALSE`, incomplete rows are
#'   removed after applying `max.NA.prop` (default: `TRUE`).
#' @param save Logical; if `TRUE`, save the result object to disk
#'   (default: `FALSE`).
#' @param overwrite Logical; if `TRUE`, overwrite an existing saved result file
#'   when `save = TRUE` (default: `FALSE`).
#' @param output.dir Character string giving the output directory used when
#'   `save = TRUE`. Required only when `save = TRUE`.
#' @param output.filename Character string giving the output file name used when
#'   `save = TRUE`. Required only when `save = TRUE`.
#' @param Sp1.background.data Optional `data.frame` of background environmental
#'   values for species 1. Required when `background.permutation.test = TRUE`
#'   (default: `NULL`).
#' @param Sp2.background.data Optional `data.frame` of background environmental
#'   values for species 2. Required when `background.permutation.test = TRUE`
#'   (default: `NULL`).
#' @param background.permutation.test Logical; if `TRUE`, run a background-based
#'   permutation test in which one species' occurrences are replaced by random
#'   samples from the other species' background, and vice versa
#'   (default: `FALSE`).
#' @param use.parallel Character string specifying whether and how to use
#'   parallel computing. One of `"auto"`, `"Windows"`, `"Unix"`, or `"none"`.
#'   `"auto"` uses `"Unix"` on Unix-like systems and `"Windows"` otherwise
#'   (recommended default: `"auto"`).
#' @param N.cores A single positive integer-like numeric value giving the number
#'   of CPU cores to use for supported parallel operations (recommended
#'   default: `3`).
#' @param seed A single numeric value used to set the random seed for
#'   reproducibility (default: `1`).
#' @param verbose Logical; if `TRUE`, progress messages are printed
#'   (default: `TRUE`).
#'
#' @details
#' Discriminant analysis of principal components (DAPC; Jombart et al., 2010)
#' is used here to test for niche divergence between two a priori groups,
#' such as species, lineages, or populations, using environmental predictors
#' associated with occurrence records.
#'
#' The method first reduces the predictor space with principal component
#' analysis (PCA) before fitting a discriminant model. This PCA step transforms
#' correlated environmental variables into orthogonal axes and reduces
#' multicollinearity and high-dimensionality problems that can make discriminant
#' analysis unstable. Discriminant analysis is then used to identify axes that
#' maximize between-group separation while minimizing within-group variability,
#' providing a supervised measure of environmental differentiation between the
#' two groups (Lachenbruch & Goldstein, 1979; Jombart et al., 2010).
#'
#' Choosing the number of retained principal components (PCs) is important because
#' retaining too few PCs can discard biologically meaningful environmental signal,
#' whereas retaining too many PCs can overfit noise and inflate apparent
#' classification accuracy. Cross-validation is therefore used to identify a PC
#' number that balances information retention and predictive stability.
#' The selected number of PCs is interpreted as a compromise between underfitting
#' and overfitting, where additional PCs no longer improve out-of-sample
#' group assignment (Jombart, 2022).
#'
#' The final model is summarized using mean assignment accuracy and the adjusted
#' Rand index. Mean assignment accuracy describes how often samples are assigned
#' to their supplied groups, providing an intuitive measure of group
#' distinctiveness in discriminant environmental space. The adjusted Rand index
#' compares predicted and supplied groupings while accounting for agreement
#' expected by chance, making it useful when evaluating whether discriminant
#' assignments recover the original taxon labels.
#'
#' Permutation testing is used to evaluate whether the observed environmental
#' separation between groups exceeds chance expectations. By randomly
#' reassigning group labels while preserving sample sizes, the test constructs a
#' null distribution representing a single shared niche with no fixed association
#' between environmental predictors and group identity. Observed assignment
#' accuracy is then compared with this null distribution to assess whether the
#' fitted discrimination is stronger than expected under random group membership.
#'
#' The optional background permutation test provides a complementary assessment
#' conditioned on the environments available to each species. Rather than
#' permuting labels alone, one species' occurrences are replaced with random
#' samples from the other species' background, and the reciprocal comparison is
#' also performed. This follows the logic of background-based tests for assessing
#' whether apparent niche similarity or divergence may be influenced by the
#' environmental conditions available within each species' accessible area
#' (Brown & Carnaval, 2019). In this framework, background-conditioned null
#' comparisons help distinguish weak divergence from limited power caused by
#' similar or constrained environmental availability.
#'
#' Missing values are handled before ordination and discrimination because these
#' methods require complete numeric inputs. Rows with excessive missingness are
#' removed, and remaining missing values can either be imputed by the median or
#' handled by complete-case filtering. The default missingness threshold
#' (`max.NA.prop = 0.2`) is intended as a practical compromise that permits
#' moderate missingness while reducing the risk that ordination and
#' discrimination are driven by imputation or incomplete records (Harrell, 2015;
#' van Buuren, 2018).
#'
#' @return A named list containing:
#'   \describe{
#'     \item{crossval_run1}{First-stage cross-validation object, or `NULL` if
#'       the first stage was skipped.}
#'     \item{optimal_pcs_crossval_run1}{Optimal number of PCs from the first
#'       cross-validation stage, or `NA` if unavailable.}
#'     \item{crossval_run2}{Second-stage cross-validation object.}
#'     \item{optimal_pcs_crossval_run2}{Final number of PCs retained for DAPC.}
#'     \item{dapc_results}{Fitted DAPC-like result object containing
#'       discriminant assignments, group labels, retained PC count, discriminant
#'       scores, variable contributions, variable loadings, group coordinates,
#'       and discriminant scaling.}
#'     \item{pca_object}{PCA object used to generate the retained PC scores.}
#'     \item{var_explained}{Cumulative proportion of variance explained by the
#'       retained PCs.}
#'     \item{observed_assign_prop}{Observed mean assignment accuracy.}
#'     \item{permutation_assign_props}{Permutation null distribution of mean
#'       assignment accuracy. If `background.permutation.test = TRUE`, this is a
#'       list with `forward` and `reverse` background permutations.}
#'     \item{p_val_assign}{Permutation p-value for assignment accuracy. If
#'       `background.permutation.test = TRUE`, this is a named vector with
#'       `forward` and `reverse` p-values.}
#'     \item{ARI}{Adjusted Rand index comparing DAPC-predicted groupings to the
#'       supplied group labels.}
#'   }
#'
#' @references
#' Brown, J., & Carnaval, A. C. (2019). A tale of two niches: Methods, concepts,
#'   and evolution. \emph{Frontiers of Biogeography}, 11(4).
#'   https://doi.org/10.21425/F5FBG44158
#'
#' Harrell, F. E. (2015). \emph{Regression modeling strategies: With applications
#'   to linear models, logistic and ordinal regression, and survival analysis}
#'   (2nd ed.). Springer.
#'
#' Jombart, T. (2022). \emph{adegenet tutorials}.
#'   https://github.com/thibautjombart/adegenet/wiki/Tutorials
#'
#' Jombart, T., Devillard, S., & Balloux, F. (2010). Discriminant analysis of
#'   principal components: A new method for the analysis of genetically
#'   structured populations. \emph{BMC Genetics}, 11, 94.
#'   https://doi.org/10.1186/1471-2156-11-94
#'
#' Lachenbruch, P. A., & Goldstein, M. (1979). Discriminant analysis.
#'   \emph{Biometrics}, 35(1), 69. https://doi.org/10.2307/2529937
#'
#' van Buuren, S. (2018). \emph{Flexible imputation of missing data}. Chapman &
#'   Hall/CRC.
#'
#' @export
run.DAPC.crossval.permutation <- function(data.input, #matrix or data.frame of predictors (rows = samples, cols = variables)
                                          N.crossval.replicates = 100, #cross-validation replicates
                                          N.permutations = 1000, #number of permutations for test
                                          fixed.n.pcs = NULL, #override n.pca used for DAPC/permutations
                                          exclude.cols = NULL, #extra columns to drop from predictors
                                          species.col = "Species", #name of column with species/group labels (character)
                                          max.NA.prop = 0.2, #max proportion of NA allowed for imputation
                                          impute.NA.median = TRUE, #impute median for NA if <= max.NA.prop
                                          save = FALSE, #save results to file
                                          overwrite = FALSE, #overwrite result file if exists (only if save = TRUE)
                                          output.dir, #directory for output (only if save = TRUE)
                                          output.filename, #Rdata filename for results (only if save = TRUE)
                                          Sp1.background.data = NULL, #background data for species 1 (only required for background.permutation.test)
                                          Sp2.background.data = NULL, #background data for species 2 (only required for background.permutation.test)
                                          background.permutation.test = FALSE, #run Humboldt-style background permutation instead of random label permutation
                                          use.parallel = c("auto", "Windows", "Unix", "none"), #whether to use parallel computing
                                          N.cores = 3, #number of cores used for parallel computing
                                          seed = 1, #set seed for reproducibility
                                          verbose = TRUE #print messages/warnings
) {

  # Clear memory
  invisible(gc())

  # Set numeric tolerance for discriminant analysis
  lda_tol <- 1e-8

  # Save and restore RNG kind
  oldRNG <- RNGkind()
  on.exit({RNGkind(kind = oldRNG[1L], normal.kind = oldRNG[2L], sample.kind = if (length(oldRNG) >= 3) oldRNG[3L] else NULL)}, add = TRUE)

  # Validate input
  if (!is.data.frame(data.input) && !is.matrix(data.input)) stop("data.input must be a data.frame or matrix")
  if (is.matrix(data.input) && is.null(colnames(data.input))) stop("data.input must have column names")
  if (!is.numeric(N.crossval.replicates) || length(N.crossval.replicates) != 1 || !is.finite(N.crossval.replicates) || N.crossval.replicates <= 0 || N.crossval.replicates %% 1 != 0) stop("N.crossval.replicates must be a positive integer (recommended: 100)")
  if (!is.numeric(N.permutations) || length(N.permutations) != 1 || !is.finite(N.permutations) || N.permutations < 1 || N.permutations %% 1 != 0) stop("N.permutations must be a positive integer (recommended: 1000)")
  if (!is.null(fixed.n.pcs) && (!is.numeric(fixed.n.pcs) || length(fixed.n.pcs) != 1L || fixed.n.pcs < 1)) stop("fixed.n.pcs must be a positive numeric value or NULL")
  if (is.null(exclude.cols)) exclude.cols <- character(0)
  else if (!is.character(exclude.cols)) stop("exclude.cols must be a character vector or NULL")
  if (!is.character(species.col) || length(species.col) != 1 || is.na(species.col) || !nzchar(species.col)) stop("species.col must be a single non-empty character")
  if (!is.numeric(max.NA.prop) || max.NA.prop < 0 || max.NA.prop > 1) stop("max.NA.prop must be between 0 and 1 (recommended: 0.2)")
  if (!is.logical(impute.NA.median) || length(impute.NA.median) != 1) stop("impute.NA.median must be TRUE or FALSE")
  if (!is.logical(save) || length(save) != 1 || is.na(save)) stop("save must be TRUE or FALSE")
  if (!is.logical(overwrite) || length(overwrite) != 1 || is.na(overwrite)) stop("overwrite must be TRUE or FALSE")
  if (save) {
    if (missing(output.dir) || missing(output.filename)) stop("Both output.dir and output.filename must be specified when save = TRUE")
    if (!is.character(output.dir) || length(output.dir) != 1 || is.na(output.dir) || !nzchar(output.dir)) stop("output.dir must be a single non-empty character")
    if (!is.character(output.filename) || length(output.filename) != 1 || is.na(output.filename) || !nzchar(output.filename)) stop("output.filename must be a single non-empty character")
    if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(output.dir)) stop("Failed to create output.dir: ", output.dir)
    output_path <- file.path(output.dir, output.filename)
  }
  if (!is.null(Sp1.background.data) && !is.data.frame(Sp1.background.data)) stop("Sp1.background.data must be a data frame or NULL")
  if (!is.null(Sp2.background.data) && !is.data.frame(Sp2.background.data)) stop("Sp2.background.data must be a data frame or NULL")
  if (!is.logical(background.permutation.test) || length(background.permutation.test) != 1 || is.na(background.permutation.test)) stop("background.permutation.test must be TRUE or FALSE")
  if (isTRUE(background.permutation.test)) {
    if (is.null(Sp1.background.data) || is.null(Sp2.background.data)) stop("background.permutation.test = TRUE but background datasets (Sp1.background.data and Sp2.background.data) are missing")
    if (!is.data.frame(Sp1.background.data) || !is.data.frame(Sp2.background.data)) stop("Both Sp1.background.data and Sp2.background.data must be data frames when background.permutation.test = TRUE")
  }
  if (!is.character(use.parallel) || !any(use.parallel %in% c("auto", "Windows", "Unix", "none"))) stop("use.parallel must be one of: 'auto', 'Windows', 'Unix', or 'none'")
  if (!is.numeric(N.cores) || length(N.cores) != 1L || N.cores < 1 || N.cores %% 1 != 0) stop("N.cores must be a positive integer")
  if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed)) stop("seed must be a single finite numeric value")
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) stop("verbose must be TRUE or FALSE")

  # Normalize parallel backend selection
  use.parallel <- match.arg(use.parallel)
  parallel <- use.parallel
  if (parallel == "auto") parallel <- if (.Platform$OS.type == "unix") "Unix" else "Windows"
  N.cores <- as.integer(max(1L, N.cores))
  .backend_for_boot <- switch(parallel, "Windows" = "snow", "Unix" = "multicore", "none" = NULL)
  .use_parallel <- !is.null(.backend_for_boot) && N.cores > 1L

  # Set seed for reproducibility
  set.seed(seed)

  # Extract species factor
  if (!species.col %in% colnames(data.input)) stop("species.col not found in data.input")
  group.assignment <- data.input[[species.col]]
  if (!is.factor(group.assignment)) {
    if (is.character(group.assignment) || is.numeric(group.assignment)) {
      group.assignment <- as.factor(group.assignment)
    } else {
      stop("species.col must contain factor/character/numeric values")
    }
  }

  # Starting message
  if (verbose) message("Processing input data ...")

  # Drop excluded columns and add species
  drop_cols <- union(exclude.cols, species.col)
  if (is.data.frame(data.input)) {
    keep_cols <- setdiff(colnames(data.input), drop_cols)
    data.input <- data.input[, keep_cols, drop = FALSE]
  } else if (is.matrix(data.input)) {
    keep_cols <- setdiff(colnames(data.input), drop_cols)
    data.input <- data.input[, keep_cols, drop = FALSE]
  }
  if (is.matrix(data.input) && is.null(colnames(data.input))) stop("data.input (matrix) must have column names to apply exclude/keep logic")

  # Remove non-numeric columns
  if (is.data.frame(data.input)) {
    non_numeric_flags <- !vapply(data.input, is.numeric, TRUE)
    n_removed <- sum(non_numeric_flags)
    n_total <- ncol(data.input)
    if (n_removed > 0 && verbose) message(n_removed, " of ", n_total, " columns were removed due to non-numeric type")
    data.input <- data.input[, !non_numeric_flags, drop = FALSE]
    if (ncol(data.input) == 0) stop("No columns left after removing non-numeric columns")
    data.input <- as.matrix(data.input)
  } else if (!is.matrix(data.input)) {
    stop("data.input must be a matrix or data.frame")
  }
  storage.mode(data.input) <- "double"

  # Check lengths and groups
  if (length(group.assignment) != nrow(data.input)) {
    stop("Length of species labels (", length(group.assignment), ") does not match number of rows in data.input (", nrow(data.input), ")")
  }
  if (nlevels(group.assignment) != 2) stop("Need two groups to run DAPC niche divergence test (found ", nlevels(group.assignment), ")")

  # Replace invalid numeric values (NaN, Inf, -Inf) with NA
  data.input[!is.finite(data.input)] <- NA

  # Remove columns with all NA or zero variance
  n_cols_before <- ncol(data.input)
  drop_col <- apply(data.input, 2, function(x) {
    if (all(is.na(x))) return("all_NA")
    sd_x <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(sd_x) || sd_x < 1e-9) return("zero_var")
    return(NA)
  })
  drop_col_reason <- drop_col[!is.na(drop_col)]
  if (length(drop_col_reason) > 0) {
    n_all_na <- sum(drop_col_reason == "all_NA")
    n_zero_var <- sum(drop_col_reason == "zero_var")
    if (verbose) {
      if (n_all_na > 0) message(n_all_na, " of ", n_cols_before, " columns dropped due to all NA")
      if (n_zero_var > 0) message(n_zero_var, " of ", n_cols_before, " columns dropped due to ~zero variance")
    }
    data.input <- data.input[, is.na(drop_col), drop = FALSE]
  }
  if (ncol(data.input) == 0) stop("No informative columns remain after variance filtering")

  # Handle missing values
  if (impute.NA.median) {
    n_before <- nrow(data.input)
    na_prop_row <- rowMeans(is.na(data.input))
    keep_rows <- na_prop_row <= max.NA.prop
    data.input <- data.input[keep_rows, , drop = FALSE]
    group.assignment <- droplevels(group.assignment[keep_rows])
    n_after <- nrow(data.input)
    n_dropped <- n_before - n_after
    if (verbose && n_dropped > 0) message(n_dropped, " of ", n_before, " rows dropped due to NA higher than max.NA.prop")
    if (nrow(data.input) == 0) stop("No rows remain after applying max.NA.prop filter")
    for (j in seq_len(ncol(data.input))) {
      if (anyNA(data.input[, j])) {
        med <- stats::median(data.input[, j], na.rm = TRUE)
        if (is.finite(med)) data.input[is.na(data.input[, j]), j] <- med
      }
    }
  } else {
    n_before <- nrow(data.input)
    na_prop_row <- rowMeans(is.na(data.input))
    keep_rows <- na_prop_row <= max.NA.prop
    data.input <- data.input[keep_rows, , drop = FALSE]
    group.assignment <- droplevels(group.assignment[keep_rows])
    n_after <- nrow(data.input)
    n_dropped <- n_before - n_after
    if (verbose && n_dropped > 0) message(n_dropped, " of ", n_before, " rows dropped due to NA higher than max.NA.prop")
    if (nrow(data.input) == 0) stop("No rows remain after applying max.NA.prop filter")
    keep_complete <- stats::complete.cases(data.input)
    data.input <- data.input[keep_complete, , drop = FALSE]
    group.assignment <- droplevels(group.assignment[keep_complete])
    if (nrow(data.input) == 0) stop("No rows remain after removing incomplete cases")
  }
  data.input <- data.input[, colMeans(is.na(data.input)) < 1, drop = FALSE]
  if (length(group.assignment) != nrow(data.input)) stop("Internal mismatch: group labels and data rows differ after filtering")

  # Re-check for zero-variance or all-NA columns after row removal
  n_cols_before2 <- ncol(data.input)
  drop_col <- apply(data.input, 2, function(x) {
    if (all(is.na(x))) return("all_NA")
    sd_x <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(sd_x) || sd_x < 1e-9) return("zero_var")
    return(NA)
  })
  drop_col_reason <- drop_col[!is.na(drop_col)]
  if (length(drop_col_reason) > 0) {
    n_all_na <- sum(drop_col_reason == "all_NA")
    n_zero_var <- sum(drop_col_reason == "zero_var")
    if (verbose) {
      if (n_all_na > 0) message(n_all_na, " of ", n_cols_before2, " columns dropped after row filtering due to all NA")
      if (n_zero_var > 0) message(n_zero_var, " of ", n_cols_before2, " columns dropped after row filtering due to ~zero variance")
    }
    data.input <- data.input[, is.na(drop_col), drop = FALSE]
  }
  if (ncol(data.input) == 0) stop("No informative columns remain after final variance filtering")

  # Check sample sizes
  group_sizes <- table(group.assignment)
  if (any(group_sizes < 3)) stop("Each group must have at least three samples (found: ", paste(group_sizes, collapse = ", "), ")")
  if (any(group_sizes < 25)) warning("Each group should have at least 25 samples (found: ", paste(group_sizes, collapse = ", "), ")")

  # Decide whether we need to compute or can load from file
  if (save) {
    need_compute <- overwrite || !file.exists(output_path)
  } else {
    need_compute <- TRUE
  }

  # Create function to calculate CV parallel with Windows
  xval.DAPC.batched <- function(x, grp, n.pca.grid, training.set, center, scale, n.rep, ncores, base_seed) {
    predictor_matrix <- x
    group_labels <- grp
    result_metric <- "groupMean"

    # Perform PCA once
    pca_for_cv <- ade4::dudi.pca(predictor_matrix,
                                 nf = min(max(n.pca.grid), ncol(predictor_matrix), nrow(predictor_matrix) - 1),
                                 scannf = FALSE,
                                 center = center,
                                 scale = scale)
    pc_scores_matrix <- as.matrix(pca_for_cv$li)
    n_pc_available <- ncol(pc_scores_matrix)
    if (n_pc_available == 0L) stop("PCA produced no components")
    colnames(pc_scores_matrix) <- sprintf("PC%03d", seq_len(n_pc_available))
    if (any(n.pca.grid > n_pc_available)) {
      warning("PCA produced only ", n_pc_available, " components - trimming CV grid")
      n.pca.grid <- n.pca.grid[n.pca.grid <= n_pc_available]
      if (length(n.pca.grid) == 0L) n.pca.grid <- 1L
    }

    # Create function to robustly replicate loop for single npc (serial inside "chunk")
    run.replicates.for.npc <- function(npc, reps) {
      success_rates <- rep(NA_real_, reps)
      n_pcs_use <- min(as.integer(npc), n_pc_available)
      if (!is.finite(n_pcs_use) || n_pcs_use < 1L) return(success_rates)

      # Replicate loop for single npc (serial inside "chunk")
      for (rep_idx in seq_len(reps)) {
        RNGkind("L'Ecuyer-CMRG")
        set.seed(base_seed + npc * 1000L + rep_idx)

        # Stratified split with >=2 per class in train, >=1 in test
        idx_train <- unlist(lapply(split(seq_along(group_labels), group_labels), function(ix) {
          n_train <- floor(length(ix) * training.set)
          n_train <- max(2L, min(n_train, length(ix) - 1L))
          sample(ix, n_train, replace = FALSE)
        }), use.names = FALSE)
        idx_train <- sort(unique(idx_train))
        idx_test <- setdiff(seq_along(group_labels), idx_train)
        if (length(idx_test) == 0L) next
        train_scores <- pc_scores_matrix[idx_train, seq_len(n_pcs_use), drop = FALSE]
        test_scores  <- pc_scores_matrix[idx_test, seq_len(n_pcs_use), drop = FALSE]
        train_priors <- prop.table(table(group_labels[idx_train]))
        if (length(train_priors) < nlevels(group_labels)) {
          missing_lvls <- setdiff(levels(group_labels), names(train_priors))
          for (lv in missing_lvls) {
            cand <- which(group_labels == lv)
            add <- setdiff(cand, idx_train)[1]
            if (!is.na(add)) idx_train <- sort(c(idx_train, add))
          }
          train_scores <- pc_scores_matrix[idx_train, seq_len(n_pcs_use), drop = FALSE]
          test_scores <- pc_scores_matrix[idx_test, seq_len(n_pcs_use), drop = FALSE]
          idx_test <- setdiff(seq_along(group_labels), idx_train)
          if (length(idx_test) == 0L) next
          test_scores <- pc_scores_matrix[idx_test, seq_len(n_pcs_use), drop = FALSE]
        }
        lda_fit <- try(MASS::lda(x = train_scores, grouping = group_labels[idx_train], prior = as.vector(train_priors), tol = lda_tol), silent = TRUE)
        if (inherits(lda_fit, "try-error")) next
        predicted_classes <- try(predict(lda_fit, test_scores)$class, silent = TRUE)
        if (inherits(predicted_classes, "try-error")) next
        true_labels <- group_labels[idx_test]
        success_rates[rep_idx] <- mean(tapply(predicted_classes == true_labels, true_labels, mean))
      }
      success_rates
    }

    # Dispatch across PC grid with chunky jobs
    if (.use_parallel && parallel == "Unix") {
      successes_list <- parallel::mclapply(n.pca.grid, run.replicates.for.npc, reps = n.rep, mc.cores = ncores)
    } else if (.use_parallel && parallel == "Windows") {
      cl <- parallel::makeCluster(ncores, type = "PSOCK")
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterSetRNGStream(cl, iseed = base_seed)
      parallel::clusterExport(cl,
                              c("pc_scores_matrix",
                                "group_labels",
                                "training.set",
                                "run.replicates.for.npc",
                                "base_seed",
                                "lda_tol"),
                              envir = environment())
      parallel::clusterEvalQ(cl, {
        NULL
      })
      successes_list <- parallel::parLapplyLB(cl, n.pca.grid, run.replicates.for.npc, reps = n.rep)
    } else {
      successes_list <- lapply(n.pca.grid, run.replicates.for.npc, reps = n.rep)
    }

    # Assemble to xval-like structure
    xval_df <- data.frame(n.pca = rep(n.pca.grid, each = n.rep), success = unlist(successes_list, use.names = FALSE))
    n.pcaF <- as.factor(xval_df$n.pca)
    successV <- xval_df$success
    mean_success_by_pc <- tapply(successV, n.pcaF, function(v) mean(v, na.rm = TRUE))
    RMSE <- tapply(successV, n.pcaF, function(v) {
      v <- v[is.finite(v)]
      if (!length(v)) return(NA_real_)
      sqrt(mean((v - 1)^2))
    })
    best_by_RMSE <- names(which(RMSE == min(RMSE, na.rm = TRUE)))
    if (length(best_by_RMSE) == 0L) {
      best_by_RMSE <- as.character(stats::median(n.pca.grid))
    }
    if (length(best_by_RMSE) > 1) best_by_RMSE <- tail(best_by_RMSE, 1)
    best_by_RMSE <- as.integer(best_by_RMSE)
    out <- list(
      `Cross-Validation Results` = xval_df,
      `Median and Confidence Interval for Random Chance` = {
        phen <- group_labels
        random <- replicate(300, mean(tapply(sample(phen) == phen, phen, mean)))
        stats::quantile(random, c(0.025, 0.5, 0.975))
      },
      `Mean Successful Assignment by Number of PCs of PCA` = mean_success_by_pc,
      `Number of PCs Achieving Highest Mean Success` = names(which.max(mean_success_by_pc)),
      `Root Mean Squared Error by Number of PCs of PCA` = RMSE,
      `Number of PCs Achieving Lowest MSE` = as.character(best_by_RMSE),
      DAPC = list(n.pca = best_by_RMSE)
    )
    out
  }

  # Create function to generate one permutation replicate (deterministic per-replicate seeding)
  permute.once <- function(rep_id, n_pcs_keep, pc_scores_full) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed + as.integer(rep_id))
    permuted_groups <- sample(group.assignment)
    permuted_groups <- factor(permuted_groups, levels = levels(group.assignment))
    if (nlevels(droplevels(permuted_groups)) < 2) return(c(NA_real_, NA_real_))
    n_pcs_use <- min(as.integer(n_pcs_keep), ncol(pc_scores_full))
    if (!is.finite(n_pcs_use) || n_pcs_use < 1L) return(c(NA_real_, NA_real_))
    pc_scores_selected <- pc_scores_full[, seq_len(n_pcs_use), drop = FALSE]
    priors_perm <- as.vector(prop.table(table(permuted_groups)))
    lda_fit_perm <- try(MASS::lda(pc_scores_selected, grouping = permuted_groups, prior = priors_perm, tol = lda_tol), silent = TRUE)
    if (inherits(lda_fit_perm, "try-error")) return(c(NA_real_, NA_real_))
    predicted_perm <- try(predict(lda_fit_perm, pc_scores_selected)$class, silent = TRUE)
    if (inherits(predicted_perm, "try-error")) return(c(NA_real_, NA_real_))
    c(mean(predicted_perm == permuted_groups), sum(predicted_perm != permuted_groups))
  }

  # Create function to facilitate PSOCK (Windows) to reduce per-job overhead and keep IDs unique
  permute.chunk <- function(start_id, chunk_size, n_pcs_keep, pc_scores_full) {
    out_mat <- matrix(NA_real_, nrow = chunk_size, ncol = 2)
    for (i in seq_len(chunk_size)) {
      rep_id <- start_id + i - 1L
      out_mat[i, ] <- permute.once(rep_id, n_pcs_keep, pc_scores_full)
    }
    out_mat
  }

  # Perform crossvalidation
  if (need_compute) {
    if (verbose) message("")
    if (verbose) message("Running cross-validation ...")
    max_possible_pcs <- min(nrow(data.input), ncol(data.input))
    smallest_group_size <- min(group_sizes)
    training_set_fraction <- min(0.9, (smallest_group_size - 1) / smallest_group_size)
    N.training <- max(2L, round(nrow(data.input) * training_set_fraction))
    if (any(apply(data.input, 2, function(x) all(is.na(x))))) stop("Some variables contain only NA values after preprocessing - cannot run PCA")
    pca_prelim <- ade4::dudi.pca(data.input,
                                 nf = max_possible_pcs,
                                 scannf = FALSE,
                                 center = TRUE,
                                 scale = TRUE)
    n.pca.max <- min(max_possible_pcs, pca_prelim$rank, N.training - 1L)

    # Report if PCA rank or LDA training-size constraint capped maximum PCs
    if (verbose) {
      if (n.pca.max == pca_prelim$rank && pca_prelim$rank < max_possible_pcs) {
        message("n.pca.max limited by PCA rank (max number of non-redundant PCs): ", pca_prelim$rank)
      }
      if (n.pca.max == (N.training - 1L) && (N.training - 1L) < max_possible_pcs) {
        message("n.pca.max limited by training-set size (DA requires PCs < number of training samples): N.training - 1 = ", N.training - 1L)
      }
    }

    # Stage 1 cross-validation
    skip_stage1 <- (n.pca.max <= 30L)
    if (!skip_stage1) {
      pc_grid_stage1 <- unique(sort(pmin(n.pca.max - 1L, seq(10L, n.pca.max, by = 10L))))
      pc_grid_stage1 <- pc_grid_stage1[pc_grid_stage1 > 0L]
      if (length(pc_grid_stage1) == 0L) pc_grid_stage1 <- 1L
      xval1 <- xval.DAPC.batched(x = data.input,
                                 grp = group.assignment,
                                 n.pca.grid = pc_grid_stage1,
                                 training.set = training_set_fraction,
                                 center = TRUE,
                                 scale = TRUE,
                                 n.rep = N.crossval.replicates,
                                 ncores = N.cores,
                                 base_seed = seed + 10000L)
      optimal_pcs_1 <- NA_integer_
      if (!is.null(xval1$DAPC) && !is.null(xval1$DAPC$n.pca)) {
        optimal_pcs_1 <- as.integer(xval1$DAPC$n.pca)
      } else if (!is.null(xval1$`Number of PCs Achieving Lowest MSE`)) {
        optimal_pcs_1 <- as.integer(xval1$`Number of PCs Achieving Lowest MSE`)
      } else if (!is.null(xval1$MSE) && !is.null(xval1$`Number of PCs tested`)) {
        pcs_tested1 <- xval1$`Number of PCs tested`
        metric1 <- -xval1$MSE
        optimal_pcs_1 <- as.integer(pcs_tested1[which.max(metric1)])
      } else if (!is.null(xval1$`Mean Successful Assignment`) && !is.null(xval1$`Number of PCs tested`)) {
        pcs_tested1 <- xval1$`Number of PCs tested`
        metric1 <- xval1$`Mean Successful Assignment`
        optimal_pcs_1 <- as.integer(pcs_tested1[which.max(metric1)])
      }
      if (!is.finite(optimal_pcs_1)) optimal_pcs_1 <- max(1L, min(10L, n.pca.max - 1L))
    } else {
      if (verbose) message("Skipping stage 1 cross-validation because n.pca.max is <= 30")
      xval1 <- NULL
      optimal_pcs_1 <- max(1L, min(10L, n.pca.max - 1L))
    }

    # Stage 2 cross-validation (fine-tuning)
    window_len <- 30L
    lower_bound <- max(1L, optimal_pcs_1 - 14L)
    upper_bound <- lower_bound + (window_len - 1L)
    upper_bound <- min(upper_bound, n.pca.max - 1L)
    lower_bound <- max(1L, upper_bound - (window_len - 1L))
    pc_grid_stage2 <- seq.int(lower_bound, upper_bound, by = 1L)
    if (length(pc_grid_stage2) == 0L) pc_grid_stage2 <- max(1L, min(n.pca.max - 1L, optimal_pcs_1))
    xval2 <- xval.DAPC.batched(x = data.input,
                               grp = group.assignment,
                               n.pca.grid = pc_grid_stage2,
                               training.set = training_set_fraction,
                               center = TRUE,
                               scale = TRUE,
                               n.rep = N.crossval.replicates,
                               ncores = N.cores,
                               base_seed = seed + 20000L)
    optimal_pcs_2 <- NA_integer_
    if (!is.null(xval2$DAPC) && !is.null(xval2$DAPC$n.pca)) {
      optimal_pcs_2 <- as.integer(xval2$DAPC$n.pca)
    } else if (!is.null(xval2$`Number of PCs Achieving Lowest MSE`)) {
      optimal_pcs_2 <- as.integer(xval2$`Number of PCs Achieving Lowest MSE`)
    } else if (!is.null(xval2$MSE) && !is.null(xval2$`Number of PCs tested`)) {
      pcs_tested2 <- xval2$`Number of PCs tested`
      metric2 <- -xval2$MSE
      optimal_pcs_2 <- as.integer(pcs_tested2[which.max(metric2)])
    } else if (!is.null(xval2$`Mean Successful Assignment`) && !is.null(xval2$`Number of PCs tested`)) {
      pcs_tested2 <- xval2$`Number of PCs tested`
      metric2 <- xval2$`Mean Successful Assignment`
      optimal_pcs_2 <- as.integer(pcs_tested2[which.max(metric2)])
    }
    if (!is.null(fixed.n.pcs)) {
      optimal_pcs_2 <- as.integer(fixed.n.pcs)
      if (verbose) message("Overriding n.pca for DAPC or permutations: ", optimal_pcs_2)
    }
    if (!is.finite(optimal_pcs_2) || optimal_pcs_2 < 1) optimal_pcs_2 <- 1L
    optimal_pcs_2 <- min(optimal_pcs_2, n.pca.max)
    if (verbose) message("Optimal number of PCs retained for DAPC: ", optimal_pcs_2)
    cumulative_variance_explained <- sum(pca_prelim$eig[seq_len(optimal_pcs_2)]) / sum(pca_prelim$eig)
    if (verbose) message("Cumulative variance explained by retained PCs: ", round(cumulative_variance_explained * 100, 1), "%")

    # Run DAPC (PCA + LDA on retained PCs)
    if (verbose) {
      message("")
      message("Running DAPC ...")
    }
    pca_full <- ade4::dudi.pca(data.input,
                               nf = max(optimal_pcs_2, 2L),
                               scannf = FALSE,
                               center = TRUE,
                               scale = TRUE)
    pc_scores_full <- as.matrix(pca_full$li)
    pc_names <- sprintf("PC%03d", seq_len(ncol(pc_scores_full)))
    colnames(pc_scores_full) <- pc_names
    n_pcs_use <- min(optimal_pcs_2, ncol(pc_scores_full))
    pc_scores_used <- pc_scores_full[, seq_len(n_pcs_use), drop = FALSE]
    class_priors_observed <- as.vector(prop.table(table(group.assignment)))
    lda_fit_observed <- MASS::lda(pc_scores_used, grouping = group.assignment, prior = class_priors_observed, tol = lda_tol)
    predicted_observed <- predict(lda_fit_observed, pc_scores_used)
    observed_assign_prop <- mean(predicted_observed$class == group.assignment)
    dapc_results_scaling <- lda_fit_observed$scaling
    pca_object_full <- pca_full

    # Extract discriminant scores
    discriminant_scores <- as.matrix(predicted_observed$x)

    # Extract PCA variable loadings and contributions
    pca_loading_matrix <- as.matrix(pca_full$c1) #variables x PCs
    lda_scaling_matrix <- as.matrix(lda_fit_observed$scaling) #PCs x LDs

    # Align by index
    var_load_matrix <- pca_loading_matrix[, seq_len(n_pcs_use), drop = FALSE] %*%
      lda_scaling_matrix[seq_len(n_pcs_use), , drop = FALSE]

    # Safe dimnames for var.load
    if (!is.null(colnames(lda_scaling_matrix)) && ncol(var_load_matrix) == ncol(lda_scaling_matrix)) {
      colnames(var_load_matrix) <- colnames(lda_scaling_matrix)
    } else {
      colnames(var_load_matrix) <- paste0("LD", seq_len(ncol(var_load_matrix)))
    }
    if (!is.null(rownames(pca_loading_matrix)) && nrow(var_load_matrix) == nrow(pca_loading_matrix)) {
      rownames(var_load_matrix) <- rownames(pca_loading_matrix)
    }

    # Normalized squared contributions per LD (sum to 1 per LD)
    sq <- var_load_matrix^2
    den <- colSums(sq)
    den[!is.finite(den) | den == 0] <- NA_real_
    var_contrib_matrix <- sweep(sq, 2, den, "/")
    var_contrib_matrix[!is.finite(var_contrib_matrix)] <- 0

    # Back-compat 1D contribution for LD1
    var_contrib <- var_contrib_matrix[, 1, drop = TRUE]

    # Group centroids in LD space (robust with 1 LD)
    grp_coord_matrix <- rowsum(discriminant_scores, group.assignment) /
      as.vector(table(group.assignment))

    # Ensure row order matches factor levels
    grp_coord_matrix <- grp_coord_matrix[levels(group.assignment), , drop = FALSE]
    colnames(grp_coord_matrix) <- colnames(discriminant_scores)

    # Assemble results
    dapc_results <- list(
      assign = predicted_observed$class,
      grp = group.assignment,
      n.pca = n_pcs_use,
      n.da = nlevels(group.assignment) - 1L,
      ld = discriminant_scores,
      var.contr = var_contrib,
      var.contr.mat = var_contrib_matrix,
      var.load = var_load_matrix,
      grp.coord = grp_coord_matrix,
      scaling = dapc_results_scaling
    )
    if (verbose) message(sprintf("Mean assignment accuracy = %.2f", observed_assign_prop))

    # Calculate and report ARI
    predicted_groups <- factor(dapc_results$assign, levels = levels(group.assignment))
    ARI <- mclust::adjustedRandIndex(predicted_groups, group.assignment)
    if (verbose) message(sprintf("Adjusted Rand Index (ARI) = %.2f", ARI))

    # Permutation null distributions
    permutation_assign_props <- numeric(N.permutations)

    # Run permutation test
    if (verbose) message("")
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed + 30000L)
    if (background.permutation.test) {
      if (verbose) message("Running Humboldt-style background permutations using shared.background ...")
      Sp1.background <- Sp1.background.data[, vapply(Sp1.background.data, is.numeric, TRUE), drop = FALSE]
      Sp2.background <- Sp2.background.data[, vapply(Sp2.background.data, is.numeric, TRUE), drop = FALSE]
      if (!identical(colnames(Sp1.background), colnames(Sp2.background))) stop("Background datasets must have identical numeric environmental variables for background permutations")
      if (ncol(Sp1.background) == 0 || ncol(Sp2.background) == 0) stop("Background datasets contain no numeric environmental columns")
      species.levels <- levels(group.assignment)
      Sp1.label <- species.levels[1]
      Sp2.label <- species.levels[2]
      if (ncol(Sp1.background) == 0 || ncol(Sp2.background) == 0) stop("Background datasets contain no numeric columns")
      Sp1.data <- as.matrix(data.input[group.assignment == Sp1.label, , drop = FALSE])
      Sp2.data <- as.matrix(data.input[group.assignment == Sp2.label, , drop = FALSE])
      observed_acc <- observed_assign_prop
      perm_assign_props <- numeric(N.permutations)
      for (perm in seq_len(N.permutations)) {
        pseudo_Sp1 <- Sp2.background[sample(seq_len(nrow(Sp2.background)), nrow(Sp1.data), replace = TRUE), , drop = FALSE]
        pseudo_combined <- rbind(pseudo_Sp1, Sp2.data)
        pseudo_labels <- factor(c(rep(Sp1.label, nrow(pseudo_Sp1)), rep(Sp2.label, nrow(Sp2.data))), levels = c(Sp1.label, Sp2.label))
        pseudo_pca <- ade4::dudi.pca(pseudo_combined, nf = n_pcs_use, scannf = FALSE, center = TRUE, scale = TRUE)
        pseudo_pc_scores <- as.matrix(pseudo_pca$li[, seq_len(min(n_pcs_use, ncol(pseudo_pca$li))), drop = FALSE])
        priors_perm <- as.vector(prop.table(table(pseudo_labels)))
        lda_perm <- try(MASS::lda(pseudo_pc_scores, grouping = pseudo_labels, prior = priors_perm, tol = lda_tol), silent = TRUE)
        if (inherits(lda_perm, "try-error")) next
        predicted_perm <- try(predict(lda_perm, pseudo_pc_scores)$class, silent = TRUE)
        if (inherits(predicted_perm, "try-error")) next
        perm_assign_props[perm] <- mean(predicted_perm == pseudo_labels)
      }
      perm_assign_props <- perm_assign_props[is.finite(perm_assign_props)]
      p_val_assign <- (sum(perm_assign_props >= observed_acc) + 1) / (length(perm_assign_props) + 1)
      if (verbose) message(sprintf("Background permutation test 1 for mean accuracy (Species 1's occurrences replaced by random samples from Species 2's background): p-value = %.3f, mean = %.2f", p_val_assign, mean(perm_assign_props)))
      permutation_assign_props <- perm_assign_props
      perm_assign_props_rev <- numeric(N.permutations)
      for (perm in seq_len(N.permutations)) {
        pseudo_Sp2 <- Sp1.background[sample(seq_len(nrow(Sp1.background)), nrow(Sp2.data), replace = TRUE), , drop = FALSE]
        pseudo_combined <- rbind(Sp1.data, pseudo_Sp2)
        pseudo_labels <- factor(c(rep(Sp1.label, nrow(Sp1.data)), rep(Sp2.label, nrow(pseudo_Sp2))), levels = c(Sp1.label, Sp2.label))
        pseudo_pca <- ade4::dudi.pca(pseudo_combined, nf = n_pcs_use, scannf = FALSE, center = TRUE, scale = TRUE)
        pseudo_pc_scores <- as.matrix(pseudo_pca$li[, seq_len(min(n_pcs_use, ncol(pseudo_pca$li))), drop = FALSE])
        priors_perm <- as.vector(prop.table(table(pseudo_labels)))
        lda_perm <- try(MASS::lda(pseudo_pc_scores, grouping = pseudo_labels, prior = priors_perm, tol = lda_tol), silent = TRUE)
        if (inherits(lda_perm, "try-error")) next
        predicted_perm <- try(predict(lda_perm, pseudo_pc_scores)$class, silent = TRUE)
        if (inherits(predicted_perm, "try-error")) next
        perm_assign_props_rev[perm] <- mean(predicted_perm == pseudo_labels)
      }
      perm_assign_props_rev <- perm_assign_props_rev[is.finite(perm_assign_props_rev)]
      p_val_assign_rev <- (sum(perm_assign_props_rev >= observed_acc) + 1) / (length(perm_assign_props_rev) + 1)
      if (verbose) message(sprintf("Background permutation test 2 for mean accuracy (Species 2's occurrences replaced by random samples from Species 1's background): p-value = %.3f, mean = %.2f", p_val_assign_rev, mean(perm_assign_props_rev)))
      permutation_assign_props <- list(forward = perm_assign_props, reverse = perm_assign_props_rev)
      p_val_assign <- c(forward = p_val_assign, reverse = p_val_assign_rev)
    } else {
      if (verbose) message("Running permutation test with ", N.permutations, " permutations ...")
      if (!.use_parallel) {
        out <- lapply(seq_len(N.permutations), permute.once, n_pcs_keep = n_pcs_use, pc_scores_full = pc_scores_full)
        out <- do.call(rbind, out)
      } else if (parallel == "Unix") {
        out <- parallel::mclapply(seq_len(N.permutations), permute.once, n_pcs_keep = n_pcs_use, pc_scores_full = pc_scores_full, mc.cores = N.cores)
        out <- do.call(rbind, out)
      } else {
        cl <- parallel::makeCluster(N.cores, type = "PSOCK")
        on.exit(parallel::stopCluster(cl), add = TRUE)
        parallel::clusterSetRNGStream(cl, seed + 30000L)
        parallel::clusterExport(cl, varlist = c("permute.chunk", "permute.once", "seed", "lda_tol", "n_pcs_use", "pc_scores_full", "group.assignment"), envir = environment())
        parallel::clusterEvalQ(cl, {
          NULL
        })
        chunk_sizes <- rep(N.permutations %/% N.cores, N.cores)
        remainder <- N.permutations %% N.cores
        if (remainder > 0) chunk_sizes[seq_len(remainder)] <- chunk_sizes[seq_len(remainder)] + 1L
        starts <- cumsum(c(1L, head(chunk_sizes, -1)))
        chunk_results <- parallel::parLapplyLB(cl, seq_along(chunk_sizes), function(i) {
          permute.chunk(start_id = starts[i], chunk_size = chunk_sizes[i], n_pcs_keep = n_pcs_use, pc_scores_full = pc_scores_full)
        })
        out <- do.call(rbind, chunk_results)
      }
      permutation_assign_props <- out[, 1]
      perm_assign_props_finite <- permutation_assign_props[is.finite(permutation_assign_props)]
      p_val_assign <- (sum(perm_assign_props_finite > observed_assign_prop) + 1) / (length(perm_assign_props_finite) + 1)
      if (verbose) message(sprintf("Permutation test for mean accuracy: p-value = %.3f, mean = %.2f", p_val_assign, mean(perm_assign_props_finite)))
    }

    # Save results
    if (save) {
      save(list = c("xval1",
                    "optimal_pcs_1",
                    "xval2",
                    "optimal_pcs_2",
                    "dapc_results",
                    "cumulative_variance_explained",
                    "observed_assign_prop",
                    "permutation_assign_props",
                    "p_val_assign",
                    "ARI"),
           file = output_path)

      if (verbose) message("DAPC finished - results saved to: ", output_path)
    }

    # Return results
    return(list(
      crossval_run1 = if (exists("xval1")) xval1 else NULL,
      optimal_pcs_crossval_run1 = if (exists("optimal_pcs_1")) optimal_pcs_1 else NA_integer_,
      crossval_run2 = xval2,
      optimal_pcs_crossval_run2 = optimal_pcs_2,
      dapc_results = dapc_results,
      pca_object = pca_object_full,
      var_explained = cumulative_variance_explained,
      observed_assign_prop = observed_assign_prop,
      permutation_assign_props = permutation_assign_props,
      p_val_assign = p_val_assign,
      ARI = ARI
    ))
  } else {

    # Load results from existing file
    if (verbose) message("Results file exists - loading results from ", output_path)
    load(output_path)

    # Return loaded results
    return(list(
      crossval_run1 = if (exists("xval1")) xval1 else NULL,
      optimal_pcs_crossval_run1 = if (exists("optimal_pcs_1")) optimal_pcs_1 else NA_integer_,
      crossval_run2 = xval2,
      optimal_pcs_crossval_run2 = optimal_pcs_2,
      dapc_results = dapc_results,
      pca_object = if (exists("pca_object_full")) pca_object_full else NULL,
      var_explained = cumulative_variance_explained,
      observed_assign_prop = observed_assign_prop,
      permutation_assign_props = permutation_assign_props,
      p_val_assign = p_val_assign,
      ARI = ARI
    ))
  }
}


#' Plot DAPC niche divergence results
#'
#' Plot smoothed LD1 density distributions for the two groups from DAPC result.
#' The discriminant axis summarizes multivariate environmental separation in one
#' dimension, allowing group overlap and divergence along the fitted DAPC axis to
#' be visualized.
#'
#' @param dapc.results DAPC result object returned by
#'   `run.DAPC.crossval.permutation()`, or a nested object containing
#'   `$dapc_results`.
#' @param alpha.density A single numeric value between `0` and `1` controlling
#'   transparency of the density fill (default: `0.75`).
#' @param group.colors Character vector of two colors used for the two groups.
#'   Can be named with group levels to enforce a specific group-color mapping
#'   (default: `c("#00005A", "darkgrey")`).
#' @param legend.label A single character string giving the legend title
#'   (default: `"Species"`).
#' @param legend.position Legend position. Either one of `"right"`, `"left"`,
#'   `"top"`, `"bottom"`, `"none"`, or a numeric vector `c(x, y)`
#'   (default: `"right"`).
#' @param legend.title.font.size A single positive numeric value giving the
#'   legend title font size (default: `9.1`).
#' @param legend.text.font.size A single positive numeric value giving the
#'   legend text font size (default: `9.1`).
#' @param legend.text.italics Logical; if `TRUE`, legend entries are italicized
#'   (default: `FALSE`).
#' @param legend.symbol.size A single positive numeric value giving the legend
#'   key size in points (default: `15`).
#' @param axis.labels.font.size A single positive numeric value giving the axis
#'   title font size (default: `9.1`).
#' @param axis.ticks.font.size A single positive numeric value giving the axis
#'   tick-label font size (default: `7`).
#' @param add.axis.lines Logical; if `TRUE`, rug lines showing individual LD1
#'   scores are added below the density curves (default: `TRUE`).
#' @param axis.lines.alpha A single numeric value between `0` and `1`
#'   controlling transparency of rug lines (default: `0.75`).
#' @param axis.lines.thickness A single positive numeric value controlling rug
#'   line thickness (default: `0.4`).
#' @param add.title Logical; if `TRUE`, a plot title is added (default: `TRUE`).
#' @param title.text A single character string giving the plot title when
#'   `add.title = TRUE` (default: `"Multivariate niche divergence based on
#'   DAPC"`).
#' @param show.plot Logical; if `TRUE`, the plot is returned visibly
#'   (default: `TRUE`).
#' @param save Logical; if `TRUE`, the plot is saved to disk
#'   (default: `FALSE`).
#' @param overwrite Logical; if `TRUE`, an existing file is overwritten when
#'   `save = TRUE` (default: `FALSE`).
#' @param filename A single character string giving the output filename without
#'   extension when `save = TRUE` (default: `"DAPC_LD1_density"`).
#' @param output.dir Optional character string giving the directory for saved
#'   plots when `save = TRUE` (default: `NULL`; if `NULL`, the current working
#'   directory is used).
#' @param type A single character string giving the output file type. One of
#'   `"png"`, `"svg"`, or `"jpg"` (default: `"svg"`).
#' @param width A single positive numeric value giving plot width in centimeters
#'   when `save = TRUE` (default: `20`).
#' @param height A single positive numeric value giving plot height in
#'   centimeters when `save = TRUE` (default: `15`).
#' @param resolution A single positive numeric value giving plot resolution in
#'   dpi when saving raster formats (default: `300`).
#' @param verbose Logical; if `TRUE`, messages are printed when saving
#'   (default: `TRUE`).
#'
#' @details
#' DAPC niche-divergence analysis reduces multivariate environmental
#' differentiation to one or more discriminant axes that maximize separation
#' between predefined groups after dimensionality reduction by PCA. For two-group
#' niche-divergence analyses, the first discriminant axis provides a univariate
#' summary of multivariate environmental separation, making it useful for
#' visualizing whether groups occupy similar, partially overlapping, or strongly
#' separated regions of discriminant environmental space (Lachenbruch & Goldstein,
#' 1979; Jombart et al., 2010).
#'
#' Smoothed density curves along LD1 provide an interpretable visualization of
#' multivariate niche divergence because they show both the location and spread of
#' group scores along the discriminant axis. Shifts in density peaks indicate
#' differences in central environmental tendency, differences in curve width
#' indicate differences in multivariate breadth along the discriminant axis, and
#' limited overlap indicates stronger environmental separation. Rug lines show
#' the underlying sample distribution and help reveal whether density features are
#' supported by many observations or by sparse tails.
#'
#' This plot is descriptive and should be interpreted alongside quantitative
#' model outputs such as assignment accuracy, permutation p-values, adjusted Rand
#' index, niche-overlap metrics, and variable contributions. Density overlap
#' along LD1 summarizes separation in the fitted discriminant space but does not
#' by itself identify independent effects of individual predictors. Variable
#' contributions from DAPC should also be interpreted cautiously when predictors
#' are correlated, because the discriminant axis reflects shared covariance
#' structure rather than fully independent predictor effects.
#'
#' @return A `ggplot` object. If `show.plot = TRUE`, the plot is returned
#'   visibly; otherwise it is returned invisibly.
#'
#' @references
#' Jombart, T., Devillard, S., & Balloux, F. (2010). Discriminant analysis of
#'   principal components: A new method for the analysis of genetically
#'   structured populations. \emph{BMC Genetics}, 11, 94.
#'   https://doi.org/10.1186/1471-2156-11-94
#'
#' Lachenbruch, P. A., & Goldstein, M. (1979). Discriminant analysis.
#'   \emph{Biometrics}, 35(1), 69. https://doi.org/10.2307/2529937
#'
#' @rawNamespace export(plot.DAPC.niche.divergence)
plot.DAPC.niche.divergence <- function(dapc.results, #DAPC result object
                                       alpha.density = 0.75, #transparency for density fill (0-1)
                                       group.colors = c("#00005A", "darkgrey"), #two colors for groups
                                       legend.label = "Species", #legend title
                                       legend.position = "right", #legend position
                                       legend.title.font.size = 9.1, #legend title font size
                                       legend.text.font.size = 9.1, #legend text font size
                                       legend.text.italics = FALSE, #italicize legend entries
                                       legend.symbol.size = 15, #legend key box size (pt)
                                       axis.labels.font.size = 9.1, #axis title font size
                                       axis.ticks.font.size = 7, #axis tick font size
                                       add.axis.lines = TRUE, #draw rug with individual LD1 scores
                                       axis.lines.alpha = 0.75, #rug transparency
                                       axis.lines.thickness = 0.4, #rug line thickness
                                       add.title = TRUE, #whether to include plot title
                                       title.text = "Multivariate niche divergence based on DAPC", #title name (only if add.title = TRUE)
                                       show.plot = TRUE, #whether to print plot to console
                                       save = FALSE, #save plot to disk
                                       overwrite = FALSE, #whether to overwrite existing file (only if save = TRUE)
                                       filename = "DAPC_LD1_density", #filename for saved plot (no extension) (only if save = TRUE)
                                       output.dir = NULL, #directory for output (only if save = TRUE)
                                       type = "svg", #plot type: "png", "svg", "jpg" (only if save = TRUE)
                                       width = 20, #plot width in cm (only if save = TRUE)
                                       height = 15, #plot height in cm (only if save = TRUE)
                                       resolution = 300, #plot resolution in dpi (only if save = TRUE)
                                       verbose = TRUE #print messages when saving
) {

  # Validate input
  if (!is.numeric(alpha.density) || length(alpha.density) != 1L || !is.finite(alpha.density) || alpha.density < 0 || alpha.density > 1) stop("alpha.density must be between 0 and 1 (recommended: 0.75)")
  if (!is.character(group.colors) || length(group.colors) != 2L || any(!nzchar(group.colors))) stop("group.colors must be a character vector of length 2")
  if (!is.character(legend.label) || length(legend.label) != 1L) stop("legend.label must be a single character string")
  if (!(is.character(legend.position) || is.numeric(legend.position))) stop("legend.position must be a character or numeric")
  if (is.character(legend.position) && !legend.position %in% c("right", "left", "top", "bottom", "none")) stop("legend.position must be 'right', 'left', 'top', 'bottom', 'none'")
  if (is.numeric(legend.position) && (length(legend.position) != 2L || any(!is.finite(legend.position)))) stop("legend.position numeric must be c(x, y) with finite values")
  if (!is.numeric(legend.title.font.size) || length(legend.title.font.size) != 1L || !is.finite(legend.title.font.size) || legend.title.font.size <= 0) stop("legend.title.font.size must be positive (recommended: 9.1)")
  if (!is.numeric(legend.text.font.size) || length(legend.text.font.size) != 1L || !is.finite(legend.text.font.size) || legend.text.font.size <= 0) stop("legend.text.font.size must be positive (recommended: 9.1)")
  if (!is.logical(legend.text.italics) || length(legend.text.italics) != 1L) stop("legend.text.italics must be TRUE or FALSE")
  if (!is.numeric(legend.symbol.size) || length(legend.symbol.size) != 1L || !is.finite(legend.symbol.size) || legend.symbol.size <= 0) stop("legend.symbol.size must be positive")
  if (!is.numeric(axis.labels.font.size) || length(axis.labels.font.size) != 1L || !is.finite(axis.labels.font.size) || axis.labels.font.size <= 0) stop("axis.labels.font.size must be a single positive number (recommended: 9.1)")
  if (!is.numeric(axis.ticks.font.size) || length(axis.ticks.font.size) != 1L || !is.finite(axis.ticks.font.size) || axis.ticks.font.size <= 0) stop("axis.ticks.font.size must be a single positive number (recommended: 7)")
  if (!is.logical(add.axis.lines) || length(add.axis.lines) != 1L) stop("add.axis.lines must be TRUE or FALSE (recommended: TRUE)")
  if (!is.numeric(axis.lines.alpha) || length(axis.lines.alpha) != 1L || !is.finite(axis.lines.alpha) || axis.lines.alpha < 0 || axis.lines.alpha > 1) stop("axis.lines.alpha must be between 0 and 1 (recommended: 0.75)")
  if (!is.numeric(axis.lines.thickness) || length(axis.lines.thickness) != 1L || !is.finite(axis.lines.thickness) || axis.lines.thickness <= 0) stop("axis.lines.thickness must be positive (recommended: 0.4)")
  if (!is.logical(add.title) || length(add.title) != 1L) stop("add.title must be TRUE or FALSE")
  if (!is.character(title.text) || length(title.text) != 1L || !nzchar(title.text)) stop("title.text must be non-empty character")
  if (!is.logical(show.plot) || length(show.plot) != 1L) stop("show.plot must be TRUE or FALSE")
  if (!is.logical(save) || length(save) != 1L) stop("save must be TRUE or FALSE")
  if (!is.logical(overwrite) || length(overwrite) != 1L) stop("overwrite must be TRUE or FALSE")
  if (!is.character(filename) || length(filename) != 1L || !nzchar(filename)) stop("filename must be non-empty")
  if (!is.null(output.dir) && !is.character(output.dir)) stop("output.dir must be NULL or character")
  type <- tolower(type)
  if (!type %in% c("png", "svg", "jpg")) stop("type must be one of: png, svg, jpg")
  if (!is.numeric(width) || length(width) != 1L || !is.finite(width) || width <= 0) stop("width must be positive (cm)")
  if (!is.numeric(height) || length(height) != 1L || !is.finite(height) || height <= 0) stop("height must be positive (cm)")
  if (!is.numeric(resolution) || length(resolution) != 1L || !is.finite(resolution) || resolution <= 0) stop("resolution must be positive (dpi)")
  if (!is.logical(verbose) || length(verbose) != 1L) stop("verbose must be TRUE or FALSE")

  # Validate DAPC result object and extract fields
  if (!is.null(dapc.results$dapc_results)) dapc.results <- dapc.results$dapc_results
  if (is.null(dapc.results) || !is.list(dapc.results)) stop("dapc.results must be DAPC list or list with $dapc_results")
  has_ld_scores <- !is.null(dapc.results$ld) || !is.null(dapc.results$ind.coord)
  if (!has_ld_scores) stop("LD scores not found: expected $ld or $ind.coord in dapc.results")
  group_labels <- dapc.results$grp
  if (is.null(group_labels)) stop("dapc.results$grp is missing")
  group_labels <- factor(group_labels)
  if (nlevels(group_labels) != 2) stop("group_labels must contain two groups")

  # Extract LD1
  if (!is.null(dapc.results$ld)) {
    ld_matrix <- as.matrix(dapc.results$ld)
  } else {
    ld_matrix <- as.matrix(dapc.results$ind.coord)
  }
  if (ncol(ld_matrix) < 1) stop("LD matrix has no columns")
  if (length(group_labels) != nrow(ld_matrix)) stop("Length of dapc.results$grp must match number of rows in LD matrix")
  ld1_scores <- as.numeric(ld_matrix[, 1])

  # Drop non-finite LD1
  keep_index <- is.finite(ld1_scores)
  if (!all(keep_index)) {
    if (verbose) warning("Dropping ", sum(!keep_index), " non-finite LD1 values")
    ld1_scores <- ld1_scores[keep_index]
    group_labels <- droplevels(group_labels[keep_index])
  }
  if (length(ld1_scores) < 4) stop("Not enough LD1 values after filtering")
  if (any(table(group_labels) < 2)) stop("Each group must have at least two finite LD1 values for density estimation")

  # Color handling
  if (!is.null(names(group.colors))) {
    if (!all(levels(group_labels) %in% names(group.colors))) stop("names(group.colors) must include all group levels: ", paste(levels(group_labels), collapse = ", "))
    group_labels <- factor(group_labels, levels = names(group.colors))
    fill_values <- group.colors[names(group.colors)]
    legend_breaks <- names(group.colors)
  } else {
    fill_values <- setNames(group.colors[1:2], levels(group_labels))
    legend_breaks <- levels(group_labels)
  }

  # Plot
  plot_data <- data.frame(LD1 = ld1_scores,
                          group = group_labels,
                          stringsAsFactors = FALSE)
  plot_object <- ggplot(plot_data,
                        aes(x = LD1, fill = group)) +
    geom_density(alpha = alpha.density) +
    scale_fill_manual(values = fill_values,
                      breaks = legend_breaks,
                      name = legend.label) +
    guides(fill = guide_legend(keyheight = grid::unit(legend.symbol.size, "pt"),
                               keywidth = grid::unit(legend.symbol.size, "pt"))) +
    labs(x = "Discriminant axis (LD1)",
         y = "Density",
         fill = legend.label) +
    theme_classic() +
    theme(axis.ticks = element_line(colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black",
                                   size = axis.ticks.font.size),
          axis.title = element_text(colour = "black",
                                    size = axis.labels.font.size,
                                    face = "bold"),
          legend.title = element_text(size = legend.title.font.size,
                                      colour = "black",
                                      face = "bold"),
          legend.text = element_text(size = legend.text.font.size,
                                     colour = "black",
                                     face = if (legend.text.italics) "italic" else "plain"),
          legend.position = legend.position,
          plot.title = element_text(hjust = 0.5,
                                    colour = "black",
                                    face = "bold",
                                    size = axis.labels.font.size))

  # Add title
  if (isTRUE(add.title)) plot_object <- plot_object + labs(title = title.text)

  # Add rug (bottom only)
  if (isTRUE(add.axis.lines)) {
    plot_object <- plot_object +
      geom_rug(data = plot_data,
               mapping = aes(x = LD1, colour = group),
               inherit.aes = FALSE,
               sides = "b",
               alpha = axis.lines.alpha,
               linewidth = axis.lines.thickness) +
      scale_colour_manual(values = fill_values,
                          breaks = legend_breaks,
                          guide = "none")
  }

  # Save
  if (isTRUE(save)) {
    if (!is.null(output.dir) && !dir.exists(output.dir)) {
      dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
      if (verbose) message("Created directory: ", output.dir)
    }
    file_out <- if (is.null(output.dir)) paste0(filename, ".", type) else file.path(output.dir, paste0(filename, ".", type))
    if (file.exists(file_out) && !overwrite) stop("File already exists: ", file_out, " (set overwrite = TRUE to replace)")
    device_for_ggsave <- if (type == "jpg") "jpeg" else type
    ggsave(filename = file_out,
           plot = plot_object,
           device = device_for_ggsave,
           width = width,
           height = height,
           units = "cm",
           dpi = resolution)
    if (verbose) message("Plot saved as ", file_out)
  }

  # Return plot
  if (show.plot) return(plot_object) else invisible(plot_object)
}


#' Plot permutation histogram for DAPC assignment accuracy
#'
#' Plot null distribution of permutation assignment accuracy together with
#' the observed assignment accuracy from DAPC niche-divergence analysis.
#'
#' @param dapc_result Result object returned by
#'   `run.DAPC.crossval.permutation()`. The object must contain
#'   `observed_assign_prop` and a numeric vector of `permutation_assign_props`.
#' @param bar.color A single character string giving the histogram bar fill color
#'   (default: `"lightgrey"`).
#' @param line.color A single character string giving the color of the observed
#'   assignment-accuracy marker (default: `"firebrick"`).
#' @param axis.labels.font.size A single positive numeric value giving the axis
#'   title font size (default: `9.1`).
#' @param axis.ticks.font.size A single positive numeric value giving the axis
#'   tick-label font size (default: `7`).
#' @param N.bar.breaks A single positive integer-like numeric value giving the
#'   number of histogram bins (default: `20`).
#' @param add.title Logical; if `TRUE`, a plot title is added (default: `TRUE`).
#' @param title.text A single character string giving the plot title when
#'   `add.title = TRUE` (default: `"Null distribution of DAPC mean assignment
#'   accuracy"`).
#' @param show.plot Logical; if `TRUE`, the plot is returned visibly
#'   (default: `TRUE`).
#' @param save Logical; if `TRUE`, the plot is saved to disk
#'   (default: `FALSE`).
#' @param overwrite Logical; if `TRUE`, an existing file is overwritten when
#'   `save = TRUE` (default: `FALSE`).
#' @param filename A single character string giving the output filename without
#'   extension when `save = TRUE` (default: `"DAPC_permutation"`).
#' @param output.dir Optional character string giving the directory for saved
#'   plots when `save = TRUE` (default: `NULL`; if `NULL`, the current working
#'   directory is used).
#' @param type A single character string giving the output file type. One of
#'   `"png"`, `"svg"`, or `"jpg"` (default: `"svg"`).
#' @param width A single positive numeric value giving plot width in centimeters
#'   when `save = TRUE` (default: `20`).
#' @param height A single positive numeric value giving plot height in
#'   centimeters when `save = TRUE` (default: `15`).
#' @param resolution A single positive numeric value giving plot resolution in
#'   dpi when saving raster formats (default: `300`).
#' @param verbose Logical; if `TRUE`, messages are printed when saving
#'   (default: `TRUE`).
#'
#' @details
#' The permutation histogram visualizes the null expectation for DAPC assignment
#' accuracy under the hypothesis of a single shared niche (k = 1).
#' In this null model, group identity is treated as exchangeable, so the
#' distribution of permuted assignment accuracies represents the level of
#' discrimination expected when environmental predictors are not consistently
#' associated with the given group labels.
#'
#' The observed assignment accuracy is shown relative to this null distribution.
#' If the observed value falls in the upper tail of the permutation distribution,
#' this indicates that the fitted discriminant separation is stronger than
#' expected under random group membership. The plot therefore provides a visual
#' complement to the permutation p-value returned by
#' `run.DAPC.crossval.permutation()`.
#'
#' This diagnostic is important because apparent group separation can arise by
#' chance, especially in high-dimensional environmental datasets. Comparing the
#' empirical assignment accuracy against a random-label null distribution helps
#' distinguish niche divergence from discrimination caused by model flexibility,
#' random label structure, or sampling imbalance.
#'
#' A significant permutation result indicates that observed group separation
#' exceeds the random-label null expectation. The plot and p-value should be
#' interpreted together with niche-overlap metrics because a significant p-value
#' can be obtained even when overall niche divergence is low.
#'
#' @return A `ggplot` object. If `show.plot = TRUE`, the plot is returned
#'   visibly; otherwise it is returned invisibly.
#'
#' @rawNamespace export(plot.DAPC.permutation)
plot.DAPC.permutation <- function(dapc_result, #DAPC result object
                                  bar.color = "lightgrey", #histogram bar fill color
                                  line.color = "firebrick", #vertical line color for observed value
                                  axis.labels.font.size = 9.1, #axis title font size
                                  axis.ticks.font.size = 7, #axis tick font size
                                  N.bar.breaks = 20, #number of histogram bins
                                  add.title = TRUE, #whether to include plot title
                                  title.text = "Null distribution of DAPC mean assignment accuracy", #title name (only if add.title = TRUE)
                                  show.plot = TRUE, #whether to print plot to console
                                  save = FALSE, #whether to save plot to disk
                                  overwrite = FALSE, #whether to overwrite existing file (only if save = TRUE)
                                  filename = "DAPC_permutation", #filename for saved plot (no extension) (only if save = TRUE)
                                  output.dir = NULL, #directory for output (only if save = TRUE)
                                  type = "svg", #plot type: "png", "svg", "jpg" (only if save = TRUE)
                                  width = 20, #plot width in cm (only if save = TRUE)
                                  height = 15, #plot height in cm (only if save = TRUE)
                                  resolution = 300, #plot resolution in dpi (only if save = TRUE)
                                  verbose = TRUE #print messages when saving
) {

  # Validate input
  if (!is.list(dapc_result)) stop("dapc_result must be a list returned by run.DAPC.crossval.permutation function - potentially rerun run.DAPC.crossval.permutation function")
  required_fields <- c("observed_assign_prop", "permutation_assign_props")
  missing_fields <- setdiff(required_fields, names(dapc_result))
  if (length(missing_fields)) stop("rerun run.DAPC.crossval.permutation function - dapc_result is missing: ", paste(missing_fields, collapse = ", "))
  if (!is.character(bar.color) || length(bar.color) != 1L || !nzchar(bar.color)) stop("bar.color must be a single color string (recommended: lightgrey)")
  if (!is.character(line.color) || length(line.color) != 1L || !nzchar(line.color)) stop("line.color must be a single color string (recommended: firebrick)")
  if (!is.numeric(axis.labels.font.size) || length(axis.labels.font.size) != 1L || !is.finite(axis.labels.font.size) || axis.labels.font.size <= 0) stop("axis.labels.font.size must be single positive number (recommended: 9.1)")
  if (!is.numeric(axis.ticks.font.size) || length(axis.ticks.font.size) != 1L || !is.finite(axis.ticks.font.size) || axis.ticks.font.size <= 0) stop("axis.ticks.font.size must be single positive number (recommended: 7)")
  if (!is.numeric(N.bar.breaks) || length(N.bar.breaks) != 1L || !is.finite(N.bar.breaks) || N.bar.breaks <= 0) stop("N.bar.breaks must be single positive number (recommended: 20)")
  if (!is.logical(add.title) || length(add.title) != 1L) stop("add.title must be TRUE or FALSE")
  if (!is.character(title.text) || length(title.text) != 1L || !nzchar(title.text)) stop("title.text must be a non-empty character")
  if (!is.logical(show.plot) || length(show.plot) != 1L) stop("show.plot must be TRUE or FALSE")
  if (!is.logical(save) || length(save) != 1L) stop("save must be TRUE or FALSE")
  if (!is.logical(overwrite) || length(overwrite) != 1L) stop("overwrite must be TRUE or FALSE")
  if (!is.character(filename) || length(filename) != 1L || !nzchar(filename)) stop("filename must be non-empty")
  if (!is.null(output.dir) && !is.character(output.dir)) stop("output.dir must be NULL or a character string")
  type <- tolower(type)
  if (!type %in% c("png", "svg", "jpg")) stop("type must be one of: png, svg, jpg")
  if (!is.numeric(width) || length(width) != 1L || !is.finite(width) || width <= 0) stop("width must be positive (cm)")
  if (!is.numeric(height) || length(height) != 1L || !is.finite(height) || height <= 0) stop("height must be positive (cm)")
  if (!is.numeric(resolution) || length(resolution) != 1L || !is.finite(resolution) || resolution <= 0) stop("resolution must be positive (dpi)")
  if (!is.logical(verbose) || length(verbose) != 1L) stop("verbose must be TRUE or FALSE")

  # Coerce bins to integer
  N.bar.breaks <- max(1L, as.integer(round(N.bar.breaks)))

  # Extract and validate data
  observed_accuracy <- dapc_result$observed_assign_prop
  permutation_accuracies <- dapc_result$permutation_assign_props
  if (!is.numeric(observed_accuracy) || length(observed_accuracy) != 1L || !is.finite(observed_accuracy)) stop("observed_assign_prop must be single finite numeric")
  if (!is.numeric(permutation_accuracies) || !length(permutation_accuracies)) stop("permutation_assign_props must be a non-empty numeric vector")
  permutation_accuracies <- permutation_accuracies[is.finite(permutation_accuracies)]
  if (!length(permutation_accuracies)) stop("No finite values in permutation_assign_props")

  # Compute empirical p-value (one-sided)
  pval.format <- function(p) {
    if (is.na(p)) list(op = "=", txt = "NA")
    else if (p < 0.001) list(op = "<", txt = "0.001")
    else if (p < 0.01) list(op = "<", txt = "0.01")
    else list(op = "=", txt = sprintf("%.2f", round(p, 2)))
  }
  p_acc_raw <- mean(permutation_accuracies >= observed_accuracy, na.rm = TRUE)
  p_acc <- pval.format(p_acc_raw)
  .pt <- 72.27 / 25.4

  # Build histogram panel
  build_histogram_panel <- function(null_distribution, observed_value, p_info, x_label) {
    min_x <- min(c(null_distribution, observed_value), na.rm = TRUE)
    max_x <- max(c(null_distribution, observed_value), na.rm = TRUE)
    if (!is.finite(min_x) || !is.finite(max_x)) stop("Non-finite values in null/observed data")
    if (max_x <= min_x) max_x <- min_x + 1e-9

    brks <- seq(min_x, max_x, length.out = N.bar.breaks + 1)
    bin_index <- cut(null_distribution,
                     breaks = brks,
                     include.lowest = TRUE,
                     right = TRUE)
    bin_counts <- tabulate(bin_index, nbins = length(brks) - 1L)
    rect_df <- data.frame(xmin = brks[-length(brks)],
                          xmax = brks[-1],
                          ymin = 0,
                          ymax = bin_counts)

    y_max <- if (length(bin_counts)) max(bin_counts) else 0
    y_max_plot <- if (is.finite(y_max) && y_max > 0) y_max else 1

    pad <- max((max_x - min_x) * 0.07, 1e-9)
    x_limits <- c(brks[1] - pad, brks[length(brks)] + pad)
    center_x <- mean(x_limits)
    center_y <- y_max_plot * 1.10

    ggplot(rect_df,
           aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
      geom_rect(fill = bar.color,
                colour = "black") +
      annotate("segment",
               x = observed_value,
               xend = observed_value,
               y = y_max_plot / 2,
               yend = 0,
               colour = line.color,
               linewidth = 1.1) +
      annotate("point",
               x = observed_value,
               y = y_max_plot / 2,
               colour = line.color,
               shape = 18,
               size = 3.5) +
      annotate("text",
               x = center_x,
               y = center_y,
               label = sprintf("p-value %s %s", p_info$op, p_info$txt),
               vjust = 0,
               hjust = 0.5,
               size = 9.1 / .pt,
               colour = "black") +
      coord_cartesian(xlim = x_limits,
                      expand = FALSE,
                      clip = "off") +
      scale_y_continuous(limits = c(0, y_max_plot * 1.18),
                         expand = expansion(mult = c(0, 0.02))) +
      labs(title = NULL,
           x = x_label,
           y = "Frequency") +
      theme_classic() +
      theme(axis.title.x = element_text(face = "bold",
                                        size = axis.labels.font.size,
                                        colour = "black"),
            axis.title.y = element_text(face = "bold",
                                        size = axis.labels.font.size,
                                        colour = "black"),
            axis.text = element_text(size = axis.ticks.font.size,
                                     colour = "black"),
            axis.line = element_line(colour = "black",
                                     linewidth = 0.5),
            axis.ticks = element_line(colour = "black",
                                      linewidth = 0.5))
  }

  # Build plot
  permutation_plot <- build_histogram_panel(permutation_accuracies, observed_accuracy, p_acc, "Assignment accuracy")

  # Add title
  if (isTRUE(add.title)) {
    permutation_plot <- permutation_plot +
      labs(title = title.text) +
      theme(plot.title = element_text(hjust = 0.5,
                                      size = axis.labels.font.size,
                                      face = "bold",
                                      colour = "black"))
  }

  # Save plot
  if (isTRUE(save)) {
    out_dir <- if (is.null(output.dir) || !nzchar(output.dir)) getwd() else output.dir
    if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      if (verbose) message("Created directory: ", out_dir) }
    file_out <- file.path(out_dir, paste0(filename, ".", type))
    if (file.exists(file_out) && !overwrite) stop("File already exists: ", file_out, " (set overwrite = TRUE to replace)")
    device_for_ggsave <- if (type == "jpg") "jpeg" else type
    ggsave(filename = file_out, plot = permutation_plot, width = width, height = height, units = "cm", dpi = resolution, device = device_for_ggsave)
    if (verbose) message("Plot saved as ", file_out)
  }

  # Return plot
  if (show.plot) return(permutation_plot) else invisible(permutation_plot)
}


## Function to calculate niche divergence metrics from DAPC densities
#' Calculate niche divergence metrics
#'
#' Calculate Schoener's D overlap and derived niche divergence metrics from
#' smoothed density distributions for two groups along the DAPC discriminant
#' axis, optionally using background-based weighting.
#'
#' @param dapc_out A DAPC result object returned by
#'   `run.DAPC.crossval.permutation()`. The object must contain `$dapc_results`
#'   with LD scores and, when `weight.background = TRUE`, the PCA object and
#'   discriminant scaling needed to project background data into LD1 space.
#' @param group.assignment Factor or character vector giving group identity for
#'   each row of the input data. Exactly two groups are required.
#' @param density.grid.resolution A single positive integer-like numeric value
#'   giving the number of grid points used to estimate smoothed density overlap
#'   (default: `1024`).
#' @param weight.background Logical; if `TRUE`, apply background-availability
#'   weighting along the discriminant axis by up-weighting rare available
#'   environments and down-weighting common available environments
#'   (default: `FALSE`).
#' @param Sp1.background.data Optional `data.frame` of background environmental
#'   values for species 1. Required when `weight.background = TRUE`; columns must
#'   overlap the environmental variables used to fit the DAPC/PCA object
#'   (default: `NULL`).
#' @param Sp2.background.data Optional `data.frame` of background environmental
#'   values for species 2. Required when `weight.background = TRUE`; columns must
#'   overlap the environmental variables used to fit the DAPC/PCA object
#'   (default: `NULL`).
#' @param verbose Logical; if `TRUE`, progress messages and metric values are
#'   printed (default: `TRUE`).
#'
#' @details
#' DAPC (discriminant analysis of principal components; Jombart et al., 2010)
#' summarizes pairwise multivariate niche separation along a single discriminant
#' axis for two-group comparisons. Calculating niche metrics from density
#' distributions along this axis allows multivariate environmental
#' differentiation to be interpreted in terms of overlap, dissimilarity,
#' exclusivity, and overall divergence magnitude. This links the discriminant
#' separation returned by DAPC to niche-divergence quantities that are easier to
#' compare across species pairs or analyses.
#'
#' Schoener's D is used as an overlap metric measuring the proportion of shared
#' density between two distributions and ranges from zero for no overlap to one
#' for complete overlap (Schoener, 1968; Warren et al., 2008). Schoener's D is
#' calculated along the DAPC discriminant axis, so it summarizes overlap in the
#' multivariate environmental space most associated with group separation.
#'
#' The additional niche-divergence metrics extend the niche divergence plane of
#' Ascanio et al. (2024) from single environmental variables to the multivariate
#' DAPC discriminant axis. Niche dissimilarity (NDS) quantifies density separation
#' along the axis, whereas niche breadth exclusivity (NE) quantifies how much of
#' the occupied discriminant-axis range is not shared between groups. Together,
#' these metrics distinguish divergence driven mainly by density differences
#' within shared environmental space from divergence driven mainly by exclusive
#' environmental ranges.
#'
#' Niche divergence magnitude (ND) combines density dissimilarity and breadth
#' exclusivity into a single composite measure of divergence strength. The niche
#' divergence angle (theta) describes the relative contribution of dissimilarity
#' versus exclusivity: angles near zero indicate divergence dominated by exclusive
#' range differences, angles near ninety degrees indicate divergence dominated by
#' density differences within shared range, and intermediate values indicate mixed
#' contributions. These quantities help classify divergence into interpretable
#' forms such as weighted, soft, nested, or hard divergence (Ascanio et al.,
#' 2024).
#'
#' When `weight.background = TRUE`, density estimates are adjusted for unequal
#' environmental availability along the discriminant axis. This follows the logic
#' of background corrections that reduce bias caused by common environments
#' dominating density estimates and rare available environments being
#' underrepresented. Such weighting can be useful when species differ in
#' accessible background environments and when niche overlap should be
#' interpreted relative to the environmental conditions available to both
#' species (Brown & Carnaval, 2019).
#'
#' @return A named list containing:
#'   \describe{
#'     \item{Schoener_D}{Schoener's D overlap between the two groups along LD1.}
#'     \item{Niche_dissimilarity}{Niche dissimilarity based on density
#'       separation.}
#'     \item{Niche_breadth_exclusivity}{Niche breadth exclusivity based on
#'       non-overlap of occupied LD1 ranges.}
#'     \item{Niche_divergence_magnitude}{Composite divergence magnitude
#'       calculated from niche dissimilarity and niche breadth exclusivity.}
#'     \item{Niche_divergence_angle_degrees}{Divergence angle in degrees,
#'       describing the relative contribution of density dissimilarity versus
#'       breadth exclusivity.}
#'     \item{niche_limits}{A nested list with lower and upper LD1 limits for
#'       `group1` and `group2`.}
#'   }
#'
#' @references
#' Ascanio, A. K., Owens, H. L., Sousa, M. C., & Peterson, A. T. (2024).
#'   Quantifying niche shifts in one-dimensional environmental space.
#'   \emph{Ecography}, 2024(5), e07127. https://doi.org/10.1111/ecog.07127
#'
#' Brown, J., & Carnaval, A. C. (2019). A tale of two niches: Methods, concepts,
#'   and evolution. \emph{Frontiers of Biogeography}, 11(4).
#'   https://doi.org/10.21425/F5FBG44158
#'
#' Jombart, T., Devillard, S., & Balloux, F. (2010). Discriminant analysis of
#'   principal components: A new method for the analysis of genetically
#'   structured populations. \emph{BMC Genetics}, 11, 94.
#'   https://doi.org/10.1186/1471-2156-11-94
#'
#' Schoener, T. W. (1968). The anolis lizards of Bimini: Resource partitioning
#'   in a complex fauna. \emph{Ecology}, 49(4), 704-726.
#'   https://doi.org/10.2307/1935534
#'
#' Warren, D. L., Glor, R. E., & Turelli, M. (2008). Environmental niche
#'   equivalency versus conservatism: Quantitative approaches to niche evolution.
#'   \emph{Evolution}, 62(11), 2868-2883.
#'   https://doi.org/10.1111/j.1558-5646.2008.00482.x
#'
#' @export
calc.niche.divergence.metrics <- function(dapc_out, #DAPC results object
                                          group.assignment, #factor/character with labels for two groups
                                          density.grid.resolution = 1024, #KDE grid size (>=100)
                                          weight.background = FALSE, #apply Humboldt-like environmental background weighting
                                          Sp1.background.data = NULL, #optional background environmental data for species 1
                                          Sp2.background.data = NULL, #optional background environmental data for species 2
                                          verbose = TRUE #print messages
) {

  # Validate input
  if (is.null(dapc_out$dapc_results) || !is.list(dapc_out$dapc_results)) stop("Input validation failed - dapc_out$dapc_results is missing or not a list - potentially rerun run.DAPC.crossval.permutation function")
  if (!is.numeric(density.grid.resolution) || length(density.grid.resolution) != 1L || !is.finite(density.grid.resolution) || density.grid.resolution < 100) stop("density.grid.resolution must be a single numeric number >= 100 (recommended: 1024)")
  if (!is.logical(weight.background) || length(weight.background) != 1L) stop("weight.background must be TRUE or FALSE")
  if (!is.null(Sp1.background.data) && !is.data.frame(Sp1.background.data)) stop("Sp1.background.data must be a data frame or NULL")
  if (!is.null(Sp2.background.data) && !is.data.frame(Sp2.background.data)) stop("Sp2.background.data must be a data frame or NULL")
  if (isTRUE(weight.background) && (is.null(Sp1.background.data) || is.null(Sp2.background.data))) stop("weight.background = TRUE but one or both background datasets are missing - provide both Sp1.background.data and Sp2.background.data")
  if (!is.logical(verbose) || length(verbose) != 1L) stop("verbose must be TRUE or FALSE")

  # Project background environmental data into LD1 space if weighting requested
  if (isTRUE(weight.background)) {
    if (!is.null(Sp1.background.data) && !is.null(Sp2.background.data)) {
      if (isTRUE(verbose)) message("Projecting background environmental data into LD1 space for Humboldt-like background weighting")
      project.background.to.LD1 <- function(background.environment.data, dapc_out) {
        if (is.null(dapc_out$dapc_results$scaling)) stop("DAPC object missing $scaling - rerun run.DAPC.crossval.permutation with save.model = TRUE")
        if (is.null(dapc_out$pca_object)) stop("PCA object missing in dapc_out - ensure DAPC input included PCA step")
        pca_object <- dapc_out$pca_object
        if (is.null(pca_object$co)) stop("Expected ade4::dudi.pca object with $co loadings")
        variable_names <- rownames(pca_object$co)
        if (is.null(variable_names)) stop("ade4 PCA lacks variable names in $co")
        common_environmental_variables <- intersect(colnames(background.environment.data), variable_names)
        if (length(common_environmental_variables) == 0) stop("No overlapping environmental variables between background data and PCA object")
        background_environment_matrix <- as.matrix(background.environment.data[, common_environmental_variables, drop = FALSE])
        background_pca_scores <- ade4::suprow(pca_object, background_environment_matrix[, colnames(pca_object$tab), drop = FALSE])$lisup
        scaling_vector <- as.matrix(as.numeric(dapc_out$dapc_results$scaling[, 1, drop = FALSE]))
        n_pcs_needed <- length(scaling_vector)
        if (ncol(background_pca_scores) < n_pcs_needed) stop("Background PCA scores have fewer PCs than used in LDA")
        background_pca_scores_used <- background_pca_scores[, seq_len(n_pcs_needed), drop = FALSE]
        background_LD1_scores <- as.numeric(as.matrix(background_pca_scores_used) %*% as.matrix(as.numeric(dapc_out$dapc_results$scaling[, 1, drop = FALSE])))
        background_LD1_scores[!is.finite(background_LD1_scores)] <- NA_real_
        background_LD1_scores
      }
      Sp1.background.LD1 <- project.background.to.LD1(Sp1.background.data, dapc_out)
      Sp2.background.LD1 <- project.background.to.LD1(Sp2.background.data, dapc_out)
      Sp1.background.LD1 <- Sp1.background.LD1[is.finite(Sp1.background.LD1)]
      Sp2.background.LD1 <- Sp2.background.LD1[is.finite(Sp2.background.LD1)]
      if (length(Sp1.background.LD1) < 5 || length(Sp2.background.LD1) < 5) stop("Too few finite background LD1 scores for KDE (need >=5 per background)")
    } else {
      stop("weight.background = TRUE but background datasets missing")
    }
  }

  # Extract LD coordinates
  dapc_results <- dapc_out$dapc_results
  if (!is.null(dapc_results$ind.coord)) {
    ld_matrix <- as.matrix(dapc_results$ind.coord)
  } else if (!is.null(dapc_results$ld)) {
    ld_matrix <- as.matrix(dapc_results$ld)
  } else {
    stop("LD scores not found - expected $ind.coord or $ld in dapc_out$dapc_results - potentially rerun run.DAPC.crossval.permutation function")
  }
  if (ncol(ld_matrix) < 1L) stop("LD matrix has no columns")

  # Align group labels to LD rownames if possible
  ld_rownames <- rownames(ld_matrix)
  group_input <- group.assignment
  if (!is.null(names(group_input)) && !is.null(ld_rownames) && setequal(names(group_input), ld_rownames)) {
    group_input <- group_input[ld_rownames]
  }

  # Validate group labels
  group_factor <- factor(group_input)
  if (nlevels(group_factor) != 2L) stop("Exactly two groups are required")
  if (any(table(group_factor) == 0L)) stop("One group has zero samples")
  if (length(group_factor) != nrow(ld_matrix)) {
    stop("Length mismatch: group.assignment (", length(group_factor), ") vs LD rows (", nrow(ld_matrix), ")")
  }

  # Extract LD1 and drop non-finite
  ld1_values <- as.numeric(ld_matrix[, 1])
  keep_idx <- is.finite(ld1_values)
  if (!all(keep_idx)) {
    warning("Dropping ", sum(!keep_idx), " non-finite LD1 values")
    ld1_values <- ld1_values[keep_idx]
    group_factor <- group_factor[keep_idx]
  }
  if (length(ld1_values) < 4) stop("Not enough LD1 values after filtering")
  if (any(table(group_factor) < 2)) stop("Each group must have >=2 finite LD1 values for KDE")

  # Split by group
  group_levels <- levels(group_factor)
  ld1_group1 <- ld1_values[group_factor == group_levels[1]]
  ld1_group2 <- ld1_values[group_factor == group_levels[2]]

  # KDE range and bandwidth
  ld1_all <- c(ld1_group1, ld1_group2)
  ld1_range <- range(ld1_all, na.rm = TRUE)
  if (!all(is.finite(ld1_range)) || ld1_range[2] <= ld1_range[1]) stop("LD1 has zero or invalid range")

  bandwidth <- tryCatch(stats::bw.SJ(ld1_all, method = "dpi"), error = function(e) NA_real_)
  if (!is.finite(bandwidth) || bandwidth <= 0) bandwidth <- stats::bw.nrd0(ld1_all)

  # Shared evaluation grid
  density_grid <- seq(ld1_range[1],
                      ld1_range[2],
                      length.out = density.grid.resolution)

  # KDEs on shared grid
  kde_group1 <- stats::density(ld1_group1, from = density_grid[1],
                               to = density_grid[length(density_grid)],
                               n = density.grid.resolution,
                               bw = bandwidth)
  kde_group2 <- stats::density(ld1_group2, from = density_grid[1],
                               to = density_grid[length(density_grid)],
                               n = density.grid.resolution,
                               bw = bandwidth)
  grid_dx <- diff(kde_group1$x[1:2])

  ## Apply Humboldt-like environmental weighting if requested
  if (isTRUE(weight.background)) {
    if (isTRUE(verbose)) message("Applying Humboldt-like background weighting: up-weighting rare environments, down-weighting common ones")
    background_kde_sp1 <- stats::density(Sp1.background.LD1, from = density_grid[1],
                                         to = density_grid[length(density_grid)],
                                         n = density.grid.resolution, bw = bandwidth)
    background_kde_sp2 <- stats::density(Sp2.background.LD1, from = density_grid[1],
                                         to = density_grid[length(density_grid)],
                                         n = density.grid.resolution, bw = bandwidth)
    combined_background_density <- (background_kde_sp1$y + background_kde_sp2$y) / 2
    combined_background_density <- combined_background_density / max(combined_background_density, na.rm = TRUE)
    epsilon <- 1e-6
    dens1 <- kde_group1$y / (combined_background_density + epsilon)
    dens2 <- kde_group2$y / (combined_background_density + epsilon)
    dens1 <- dens1 / sum(dens1 * grid_dx)
    dens2 <- dens2 / sum(dens2 * grid_dx)
  } else {
    dens1 <- kde_group1$y
    dens2 <- kde_group2$y
  }

  # Max-normalization (Ascanio et al. 2024, Eq. 1): f_tilde(x) = f(x) / max_x f(x)
  peak1 <- max(dens1, na.rm = TRUE)
  if (!is.finite(peak1) || peak1 <= 0) stop("Non-positive peak (group 1)")
  peak2 <- max(dens2, na.rm = TRUE)
  if (!is.finite(peak2) || peak2 <= 0) stop("Non-positive peak (group 2)")
  f_norm1 <- dens1 / peak1
  f_norm2 <- dens2 / peak2

  # Support-based niche limits for NE (Ascanio et al. 2024, Eq. 4)
  support_threshold <- 1e-6
  support_idx1 <- which(f_norm1 > support_threshold)
  if (length(support_idx1) < 2) stop("Insufficient support for group 1")
  support_idx2 <- which(f_norm2 > support_threshold)
  if (length(support_idx2) < 2) stop("Insufficient support for group 2")
  niche1_lower <- density_grid[min(support_idx1)]
  niche1_upper <- density_grid[max(support_idx1)]
  niche2_lower <- density_grid[min(support_idx2)]
  niche2_upper <- density_grid[max(support_idx2)]

  # Create additional functions
  trapz.integral <- function(x, y) sum((head(y, -1) + tail(y, -1)) * diff(x) / 2)
  grid.between <- function(lower, upper, grid) which(grid >= lower & grid <= upper)

  # One-sided dissimilarities DS (Ascanio et al. 2024, Eq. 3): DS(A|B) = 1 - integral_A min(f_tilde_A, f_tilde_B) dx / integral_A f_tilde_A x and DS(B|A) analogous on B
  idx_A <- grid.between(niche1_lower, niche1_upper, density_grid)
  if (length(idx_A) < 2) stop("Too few grid points in group 1 limits")
  idx_B <- grid.between(niche2_lower, niche2_upper, density_grid)
  if (length(idx_B) < 2) stop("Too few grid points in group 2 limits")
  area_A_total <- trapz.integral(density_grid[idx_A], f_norm1[idx_A])
  if (!is.finite(area_A_total) || area_A_total <= 0) stop("Zero integral on group 1")
  area_A_overlap <- trapz.integral(density_grid[idx_A], pmin(f_norm1[idx_A], f_norm2[idx_A]))
  ds_group1_given_group2 <- (area_A_total - area_A_overlap) / area_A_total
  area_B_total <- trapz.integral(density_grid[idx_B], f_norm2[idx_B])
  if (!is.finite(area_B_total) || area_B_total <= 0) stop("Zero integral on group 2")
  area_B_overlap <- trapz.integral(density_grid[idx_B], pmin(f_norm1[idx_B], f_norm2[idx_B]))
  ds_group2_given_group1 <- (area_B_total - area_B_overlap) / area_B_total

  # Niche dissimilarity (NDS) (Ascanio et al. 2024, Eq. 2)
  niche_dissimilarity <- (ds_group1_given_group2 + ds_group2_given_group1) / 2
  niche_dissimilarity <- max(0, min(1, niche_dissimilarity))

  # Niche breadth exclusivity (NE) (Ascanio et al. 2024, Eq. 4)
  union_span <- max(niche1_upper, niche2_upper) - min(niche1_lower, niche2_lower)
  overlap_span <- max(0, min(niche1_upper, niche2_upper) - max(niche1_lower, niche2_lower))
  niche_exclusivity <- if (union_span > .Machine$double.eps) 1 - (overlap_span / union_span) else 0
  niche_exclusivity <- max(0, min(1, niche_exclusivity))

  # Niche divergence magnitude (ND) and angle (theta) (Ascanio et al. 2024, Eqs. 5-6): ND = sqrt(NDS^2 + NE^2) and theta = atan2(NDS, NE) in degrees
  niche_divergence_magnitude <- sqrt(niche_dissimilarity^2 + niche_exclusivity^2)
  niche_divergence_angle_degrees <- atan2(niche_dissimilarity, niche_exclusivity) * (180 / pi)

  # Schoener's D overlap (Schoener 1968): D = integral min(q1(x), q2(x)) dx, with q normalized to area 1
  area1 <- sum(dens1 * grid_dx)
  if (!is.finite(area1) || area1 <= 0) stop("KDE normalization failed (group 1)")
  area2 <- sum(dens2 * grid_dx)
  if (!is.finite(area2) || area2 <= 0) stop("KDE normalization failed (group 2)")
  q1 <- dens1 / area1
  q2 <- dens2 / area2
  schoener_D <- sum(pmin(q1, q2)) * grid_dx
  schoener_D <- max(0, min(1, schoener_D))

  # Messages
  if (isTRUE(verbose)) {
    message("")
    message("Niche divergence metrics:")
    message(sprintf("Schoener's D overlap = %.2f", schoener_D))
    message(sprintf("Niche dissimilarity (NDS) = %.2f", niche_dissimilarity))
    message(sprintf("Niche breadth exclusivity (NE) = %.2f", niche_exclusivity))
    message(sprintf("Niche divergence magnitude (ND) = %.2f", niche_divergence_magnitude))
    message(sprintf("Niche divergence angle (theta) = %.1f degree", niche_divergence_angle_degrees))
  }

  # Return results
  list(
    Schoener_D = schoener_D,
    Niche_dissimilarity = niche_dissimilarity,
    Niche_breadth_exclusivity = niche_exclusivity,
    Niche_divergence_magnitude = niche_divergence_magnitude,
    Niche_divergence_angle_degrees = niche_divergence_angle_degrees,
    niche_limits = list(
      group1 = c(lower = niche1_lower, upper = niche1_upper),
      group2 = c(lower = niche2_lower, upper = niche2_upper)
    )
  )
}


## Function to plot DAPC variable contributions
#' Plot DAPC variable contributions
#'
#' Plot the relative contributions of environmental variables to the DAPC
#' discriminant axis and indicate group associated with higher values for
#' each variable.
#'
#' @param dapc.results DAPC result object returned by
#'   `run.DAPC.crossval.permutation()`, or a nested object containing
#'   `$dapc_results`. The object must contain variable contributions, variable
#'   loadings, and group coordinates.
#' @param group.colors Character vector of two colors used for the plotted groups
#'   or contribution directions. Can be named with group levels to enforce a
#'   specific group-color mapping (default: `c("#00005A", "darkgrey")`).
#' @param min.contribution.threshold Optional single non-negative numeric value
#'   giving the minimum contribution required for a variable to be plotted
#'   (default: `NULL`).
#' @param top.N Optional single positive integer-like numeric value limiting the
#'   plot to the top contributing variables after ordering by contribution
#'   (default: `NULL`).
#' @param axis.labels.font.size A single positive numeric value giving the axis
#'   title font size (default: `9.1`).
#' @param axis.ticks.font.size A single positive numeric value giving the axis
#'   tick-label font size (default: `7`).
#' @param title.font.size A single positive numeric value giving the plot title
#'   font size (default: `9.1`).
#' @param legend.title.font.size A single positive numeric value giving the
#'   legend title font size (default: `9.1`).
#' @param legend.text.font.size A single positive numeric value giving the legend
#'   text font size (default: `9.1`).
#' @param legend.text.italics Logical; if `TRUE`, legend entries are italicized
#'   (default: `FALSE`).
#' @param legend.symbol.size A single positive numeric value giving the legend
#'   key size in points (default: `15`).
#' @param legend.position Legend position. Either one of `"right"`, `"left"`,
#'   `"top"`, `"bottom"`, `"none"`, or a numeric vector `c(x, y)`
#'   (default: `"right"`).
#' @param add.title Logical; if `TRUE`, a plot title is added (default: `TRUE`).
#' @param title.text A single character string giving the plot title when
#'   `add.title = TRUE` (default: `"Variable contributions to niche
#'   divergence"`).
#' @param show.plot Logical; if `TRUE`, the plot is displayed
#'   (default: `TRUE`).
#' @param save Logical; if `TRUE`, the plot is saved to disk
#'   (default: `FALSE`).
#' @param overwrite Logical; if `TRUE`, an existing file is overwritten when
#'   `save = TRUE` (default: `FALSE`).
#' @param filename A single character string giving the output filename without
#'   extension when `save = TRUE` (default: `"DAPC_var_contributions"`).
#' @param output.dir Optional character string giving the directory for saved
#'   plots when `save = TRUE` (default: `NULL`; if `NULL`, the current working
#'   directory is used).
#' @param type A single character string giving the output file type. One of
#'   `"png"`, `"svg"`, or `"jpg"` (default: `"svg"`).
#' @param width A single positive numeric value giving plot width in centimeters
#'   when `save = TRUE` (default: `20`).
#' @param height A single positive numeric value giving plot height in
#'   centimeters when `save = TRUE` (default: `15`).
#' @param resolution A single positive numeric value giving plot resolution in
#'   dpi when saving raster formats (default: `300`).
#' @param verbose Logical; if `TRUE`, progress messages are printed
#'   (default: `TRUE`).
#'
#' @details
#' Variable contributions identify environmental predictors that contribute most
#' strongly to separation along the DAPC discriminant axis. The discriminant axis
#' is a linear combination of retained principal components, which are themselves
#' linear combinations of the original predictors. Back-transforming the
#' discriminant axis to the original predictor space gives each variable a
#' loading that reflects its weight in the fitted discrimination.
#'
#' Contributions are calculated from squared loadings normalized to sum to one,
#' so each value represents the relative proportion of discrimination associated
#' with a predictor. Larger values indicate stronger contribution to separation
#' along the discriminant axis, whereas values near zero indicate little
#' contribution to the fitted group separation.
#'
#' The direction of effect is determined from the sign of each variable loading
#' relative to the group centroids on the discriminant axis. This indicates which
#' group is associated with higher values of each predictor along the fitted
#' discriminant direction. The direction should be interpreted as an axis-based
#' association, not as evidence of an independent causal effect.
#'
#' Variable contributions should be interpreted cautiously when predictors are
#' strongly correlated. PCA reduces collinearity for model fitting, but the
#' back-transformed contributions still reflect shared covariance structure among
#' the original predictors. Therefore, high contributions indicate variables or
#' correlated variable sets associated with niche divergence, rather than fully
#' independent effects of individual predictors.
#'
#' @return A `data.frame` containing the plotted variable contributions,
#'   loadings, and inferred contribution directions. If `show.plot = TRUE`, the
#'   plot is displayed as a side effect. If `save = TRUE`, the plot is also saved
#'   to disk using `filename`, `output.dir`, and `type`.
#'
#' @rawNamespace export(plot.DAPC.var.contributions)
plot.DAPC.var.contributions <- function(dapc.results, #DAPC object
                                        group.colors = c("#00005A", "darkgrey"), #two colors for groups
                                        min.contribution.threshold = NULL, #drop variables below this contribution
                                        top.N = NULL, #number of top variables to show
                                        axis.labels.font.size = 9.1, #axis title font size
                                        axis.ticks.font.size = 7, #axis tick font size
                                        title.font.size = 9.1, #title font size
                                        legend.title.font.size = 9.1, #legend title font size
                                        legend.text.font.size = 9.1, #legend text font size
                                        legend.text.italics = FALSE, #italicize legend text
                                        legend.symbol.size = 15, #legend key box size (pt)
                                        legend.position = "right", #legend position ("right", "left", "top", "bottom", "none")
                                        add.title = TRUE, #whether to include plot title
                                        title.text = "Variable contributions to niche divergence", #title name (only if add.title = TRUE)
                                        show.plot = TRUE, #whether to print plot to console
                                        save = FALSE, #save plot to disk
                                        overwrite = FALSE, #overwrite existing file if TRUE (only if save = TRUE)
                                        filename = "DAPC_var_contributions", #filename for saved plot (no extension) (only if save = TRUE)
                                        output.dir = NULL, #directory for output (only if save = TRUE)
                                        type = "svg", #plot type: "png", "svg", "jpg" (only if save = TRUE)
                                        width = 20, #plot width in cm (only if save = TRUE)
                                        height = 15, #plot height in cm (only if save = TRUE)
                                        resolution = 300, #plot resolution in dpi (only if save = TRUE)
                                        verbose = TRUE #print messages
) {

  # Validate input
  if (!is.null(dapc.results$dapc_results)) dapc.results <- dapc.results$dapc_results
  if (is.null(dapc.results) || !is.list(dapc.results)) stop("dapc.results must be a DAPC list or a list with $dapc_results")
  if (is.null(dapc.results$var.contr)) stop("dapc.results$var.contr is missing")
  if (is.null(dapc.results$var.load)) stop("dapc.results$var.load is missing")
  if (is.null(dapc.results$grp.coord)) stop("dapc.results$grp.coord is missing")
  if (ncol(as.matrix(dapc.results$var.contr)) < 1) stop("var.contr must have >=1 column (LD1)")
  if (ncol(as.matrix(dapc.results$var.load)) < 1) stop("var.load must have >=1 column (LD1)")
  if (!is.character(group.colors) || length(group.colors) != 2L || any(!nzchar(group.colors))) stop("group.colors must be a character vector of length 2")
  if (!is.null(min.contribution.threshold) && (!is.numeric(min.contribution.threshold) || length(min.contribution.threshold) != 1L || !is.finite(min.contribution.threshold) || min.contribution.threshold < 0)) stop("min.contribution.threshold must be NULL or a single non-negative number (recommended: NULL)")
  if (!is.null(top.N)) {
    if (!is.numeric(top.N) || length(top.N) != 1L || !is.finite(top.N) || top.N < 1) stop("top.N must be NULL or a single positive integer (recommended: NULL)")
    top.N <- as.integer(top.N)
  }
  if (!is.numeric(axis.labels.font.size) || length(axis.labels.font.size) != 1L || !is.finite(axis.labels.font.size) || axis.labels.font.size <= 0) stop("axis.labels.font.size must be a single positive number (recommended: 9.1)")
  if (!is.numeric(axis.ticks.font.size) || length(axis.ticks.font.size) != 1L || !is.finite(axis.ticks.font.size) || axis.ticks.font.size <= 0) stop("axis.ticks.font.size must be a single positive number (recommended: 7)")
  if (!is.numeric(title.font.size) || length(title.font.size) != 1L || !is.finite(title.font.size) || title.font.size <= 0) stop("title.font.size must be a single positive number (recommended: 9.1)")
  if (!is.numeric(legend.title.font.size) || length(legend.title.font.size) != 1L || !is.finite(legend.title.font.size) || legend.title.font.size <= 0) stop("legend.title.font.size must be a single positive number (recommended: 9.1)")
  if (!is.numeric(legend.text.font.size) || length(legend.text.font.size) != 1L || !is.finite(legend.text.font.size) || legend.text.font.size <= 0) stop("legend.text.font.size must be a single positive number (recommended: 9.1)")
  if (!is.logical(legend.text.italics) || length(legend.text.italics) != 1L) stop("legend.text.italics must be TRUE or FALSE")
  if (!is.numeric(legend.symbol.size) || length(legend.symbol.size) != 1L || !is.finite(legend.symbol.size) || legend.symbol.size <= 0) stop("legend.symbol.size must be a single positive number (pt) (recommended: 15)")
  if (!(is.character(legend.position) || is.numeric(legend.position))) stop("legend.position must be character or numeric")
  if (is.character(legend.position) && !legend.position %in% c("right", "left", "top", "bottom", "none")) stop("legend.position must be 'right', 'left', 'top', 'bottom', or 'none'")
  if (is.numeric(legend.position) && (length(legend.position) != 2L || any(!is.finite(legend.position)))) stop("legend.position numeric must be c(x, y) with finite values")
  if (!is.logical(add.title) || length(add.title) != 1L) stop("add.title must be TRUE or FALSE")
  if (!is.character(title.text) || length(title.text) != 1L || !nzchar(title.text)) stop("title.text must be a non-empty character")
  if (!is.logical(show.plot) || length(show.plot) != 1L) stop("show.plot must be TRUE or FALSE")
  if (!is.logical(save) || length(save) != 1L) stop("save must be TRUE or FALSE")
  if (!is.logical(overwrite) || length(overwrite) != 1L) stop("overwrite must be TRUE or FALSE")
  if (!is.character(filename) || length(filename) != 1L || !nzchar(filename)) stop("filename must be non-empty")
  if (!is.null(output.dir) && !is.character(output.dir)) stop("output.dir must be NULL or character")
  type <- tolower(type)
  if (!type %in% c("png", "svg", "jpg")) stop("type must be one of: png, svg, jpg")
  if (!is.numeric(width) || length(width) != 1L || !is.finite(width) || width <= 0) stop("width must be positive (cm)")
  if (!is.numeric(height) || length(height) != 1L || !is.finite(height) || height <= 0) stop("height must be positive (cm)")
  if (!is.numeric(resolution) || length(resolution) != 1L || !is.finite(resolution) || resolution <= 0) stop("resolution must be positive (dpi)")
  if (!is.logical(verbose) || length(verbose) != 1L) stop("verbose must be TRUE or FALSE")

  # Extract contributions and loadings (LD1)
  vc_mat <- as.matrix(dapc.results$var.contr)
  vl_mat <- as.matrix(dapc.results$var.load)
  variable_contribution <- as.numeric(vc_mat[, 1, drop = TRUE])
  variable_loading <- as.numeric(vl_mat[, 1, drop = TRUE])
  names(variable_contribution) <- rownames(vc_mat)
  names(variable_loading) <- rownames(vl_mat)
  if (is.null(names(variable_contribution)) || any(names(variable_contribution) == "")) stop("var.contr must have rownames")
  if (is.null(names(variable_loading)) || any(names(variable_loading) == "")) stop("var.load must have rownames")
  if (any(!is.finite(variable_contribution))) stop("var.contr contains non-finite values")
  if (any(!is.finite(variable_loading))) stop("var.load contains non-finite values")
  variable_loading <- variable_loading[names(variable_contribution)]

  # Filter by threshold
  if (!is.null(min.contribution.threshold)) {
    keep_index <- which(variable_contribution >= min.contribution.threshold)
    variable_contribution <- variable_contribution[keep_index]
    variable_loading <- variable_loading[keep_index]
    if (!length(variable_contribution)) stop("No variables meet min.contribution.threshold")
  }

  # Order and subset
  order_index <- order(variable_contribution, decreasing = TRUE)
  variable_names <- names(variable_contribution)[order_index]
  variable_contribution <- variable_contribution[order_index]
  variable_loading <- variable_loading[order_index]
  if (!length(variable_names)) stop("No variables available to plot after filtering")
  if (!is.null(top.N) && length(variable_names) > top.N) {
    variable_names <- variable_names[1:top.N]
    variable_contribution <- variable_contribution[1:top.N]
    variable_loading <- variable_loading[1:top.N]
    if (verbose) message("Showing ", length(variable_names), " variables (top.N)")
  }

  # Determine direction mapping from group centroids along LD1
  group_coordinates <- as.matrix(dapc.results$grp.coord)
  if (ncol(group_coordinates) < 1 || nrow(group_coordinates) < 2) stop("grp.coord must contain LD1 for >=2 groups")
  group_means <- group_coordinates[, 1]
  if (is.null(names(group_means)) || length(group_means) < 2) stop("grp.coord must have rownames (group names)")
  group_positive <- names(group_means)[which.max(group_means)]
  group_negative <- names(group_means)[which.min(group_means)]
  if (identical(group_positive, group_negative)) {
    all_groups <- rownames(group_coordinates)
    group_negative <- setdiff(all_groups, group_positive)[1]
  }

  # Legend levels and color names
  if (is.null(names(group.colors))) names(group.colors) <- rownames(dapc.results$grp.coord)
  legend_levels <- names(group.colors)
  if (!all(c(group_positive, group_negative) %in% legend_levels)) stop("names(group.colors) must include all group names in grp.coord")

  # Build plotting data
  contribution_df <- data.frame(
    Variable = factor(variable_names, levels = variable_names),
    Contribution = as.numeric(variable_contribution),
    Loading = as.numeric(variable_loading),
    Direction = factor(ifelse(variable_loading > 0, group_positive, group_negative), levels = legend_levels),
    stringsAsFactors = FALSE
  )

  # Plot
  contribution_plot <- ggplot(contribution_df, aes(x = Contribution, y = Variable, fill = Direction)) +
    geom_col() +
    scale_x_continuous(expand = c(0, 0.002)) +
    scale_fill_manual(breaks = legend_levels, values = group.colors[legend_levels]) +
    labs(x = "Variable contribution", y = NULL, fill = "Higher values in") +
    theme_classic() +
    theme(
      axis.ticks = element_line(colour = "black"),
      axis.line = element_line(colour = "black"),
      axis.title = element_text(size = axis.labels.font.size, colour = "black", face = "bold"),
      axis.text = element_text(size = axis.ticks.font.size, colour = "black"),
      legend.title = element_text(size = legend.title.font.size, colour = "black", face = "bold"),
      legend.text = element_text(size = legend.text.font.size, colour = "black",
                                 face = if (legend.text.italics) "italic" else "plain"),
      legend.key.size = grid::unit(legend.symbol.size, "pt"),
      legend.position = legend.position,
      plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = title.font.size)
    )

  # Add title
  if (isTRUE(add.title)) {
    contribution_plot <- contribution_plot + labs(title = title.text)
  }

  # Save
  if (isTRUE(save)) {
    out_dir <- if (is.null(output.dir) || !nzchar(output.dir)) getwd() else output.dir
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      if (!dir.exists(out_dir)) stop("Failed to create output.dir: ", out_dir)
    }
    file_out <- file.path(out_dir, paste0(filename, ".", type))
    if (file.exists(file_out) && !overwrite) stop("File already exists: ", file_out, " (set overwrite = TRUE to replace)")
    device_for_ggsave <- if (type == "jpg") "jpeg" else type
    ggsave(filename = file_out,
           plot = contribution_plot,
           width = width,
           height = height,
           units = "cm",
           dpi = resolution,
           device = device_for_ggsave)
    if (verbose) message("Plot saved as ", file_out)
  }

  # Return plot
  if (show.plot) print(contribution_plot) else invisible(contribution_plot)

  # Return results
  return(data = contribution_df)
}


## Function to plot top DAPC predictors
#' Plot top DAPC predictors
#'
#' Plot density distributions of the top environmental predictors contributing to
#' DAPC separation.
#'
#' @param dapc.results DAPC result object returned by
#'   `run.DAPC.crossval.permutation()`, or a nested object containing
#'   `$dapc_results`. The object must contain variable contribution information.
#' @param group.colors Character vector of two colors used for the plotted groups.
#'   Can be named with group levels to enforce a specific group-color mapping
#'   (default: `c("#00005A", "darkgrey")`).
#' @param predictor.data A `data.frame` or matrix containing the raw predictor
#'   values to plot for the top contributing variables.
#' @param species.labels Factor or character vector giving group identity for
#'   rows in `predictor.data`. Exactly two groups are required.
#' @param N.top.variables A single positive integer-like numeric value giving the
#'   number of top contributing variables to plot (default: `6`).
#' @param alpha.density A single numeric value between `0` and `1` controlling
#'   transparency of density fills (default: `0.75`).
#' @param axis.ticks.font.size A single positive numeric value giving the axis
#'   tick-label font size (default: `7`).
#' @param axis.labels.font.size A single positive numeric value giving the axis
#'   title font size (default: `9.1`).
#' @param variable.font.size A single positive numeric value giving the facet
#'   label font size for variable names (default: `9.1`).
#' @param legend.label A single character string giving the legend title
#'   (default: `"Species"`).
#' @param legend.position Legend position. Either one of `"right"`, `"left"`,
#'   `"top"`, `"bottom"`, `"none"`, or a numeric vector `c(x, y)`
#'   (default: `"right"`).
#' @param legend.title.font.size A single positive numeric value giving the
#'   legend title font size (default: `9.1`).
#' @param legend.text.font.size A single positive numeric value giving the legend
#'   text font size (default: `9.1`).
#' @param legend.text.italics Logical; if `TRUE`, legend entries are italicized
#'   (default: `FALSE`).
#' @param legend.symbol.size A single positive numeric value giving the legend
#'   key size in points (default: `15`).
#' @param add.title Logical; if `TRUE`, a plot title is added (default: `TRUE`).
#' @param title.font.size A single positive numeric value giving the plot title
#'   font size (default: `9.1`).
#' @param title.text Optional character string giving the plot title when
#'   `add.title = TRUE`. If `NULL`, a title is generated automatically
#'   (default: `NULL`).
#' @param plot.nrow Optional single positive integer-like numeric value giving
#'   the number of rows in the faceted plot layout. If `NULL`, this is determined
#'   automatically (default: `NULL`).
#' @param plot.ncol Optional single positive integer-like numeric value giving
#'   the number of columns in the faceted plot layout. If `NULL`, this is
#'   determined automatically (default: `NULL`).
#' @param show.plot Logical; if `TRUE`, the plot is returned visibly
#'   (default: `TRUE`).
#' @param save Logical; if `TRUE`, the plot is saved to disk
#'   (default: `FALSE`).
#' @param overwrite Logical; if `TRUE`, an existing file is overwritten when
#'   `save = TRUE` (default: `FALSE`).
#' @param filename A single character string giving the output filename without
#'   extension when `save = TRUE` (default: `"top_DAPC_predictors"`).
#' @param output.dir Optional character string giving the directory for saved
#'   plots when `save = TRUE` (default: `NULL`; if `NULL`, the current working
#'   directory is used).
#' @param type A single character string giving the output file type. One of
#'   `"png"`, `"svg"`, or `"jpg"` (default: `"svg"`).
#' @param width A single positive numeric value giving plot width in centimeters
#'   when `save = TRUE` (default: `16`).
#' @param height A single positive numeric value giving plot height in
#'   centimeters when `save = TRUE` (default: `10`).
#' @param resolution A single positive numeric value giving plot resolution in
#'   dpi when saving raster formats (default: `300`).
#' @param verbose Logical; if `TRUE`, progress messages are printed
#'   (default: `TRUE`).
#'
#' @details
#' This plot links the multivariate DAPC discriminant axis back to the original
#' environmental predictors by showing the distributions of the variables with
#' the highest contributions to group separation. These top-predictor density
#' plots show how the most influential original variables differ between groups.
#' Variables are ranked by their DAPC contribution values.
#'
#' @return A `ggplot` object. If `show.plot = TRUE`, the plot is returned
#'   visibly; otherwise it is returned invisibly.
#'
#' @rawNamespace export(plot.top.DAPC.predictors)
plot.top.DAPC.predictors <- function(dapc.results, #DAPC result object
                                     group.colors = c("#00005A", "darkgrey"), #two colors for groups
                                     predictor.data, #data.frame or matrix with predictor variables
                                     species.labels, #factor or vector with group labels (two levels)
                                     N.top.variables = 6, #number of top variables to plot
                                     alpha.density = 0.75, #transparency for density fill (0-1)
                                     axis.ticks.font.size = 7, #axis ticks font size
                                     axis.labels.font.size = 9.1, #axis labels font size
                                     variable.font.size = 9.1, #facet strip text font size
                                     legend.label = "Species", #legend label text
                                     legend.position = "right", #legend position ("right", "left", "top", "bottom", "none")
                                     legend.title.font.size = 9.1, #legend title font size
                                     legend.text.font.size = 9.1, #legend labels font size
                                     legend.text.italics = FALSE, #italicize legend text
                                     legend.symbol.size = 15, #legend symbol size (pt)
                                     add.title = TRUE, #whether to include plot title
                                     title.font.size = 9.1, #plot title font size
                                     title.text = NULL, #title name (if NULL, default is generated; only if add.title = TRUE)
                                     plot.nrow = NULL, #facet rows (NULL = auto)
                                     plot.ncol = NULL, #facet cols (NULL = auto)
                                     show.plot = TRUE, #whether to print plot to console
                                     save = FALSE, #save plot to disk
                                     overwrite = FALSE, #overwrite existing file if TRUE (only if save = TRUE)
                                     filename = "top_DAPC_predictors", #filename for saved plot (no extension) (only if save = TRUE)
                                     output.dir = NULL, #directory for output (only if save = TRUE)
                                     type = "svg", #plot type: "png", "svg", "jpg" (only if save = TRUE)
                                     width = 16, #plot width in cm (only if save = TRUE)
                                     height = 10, #plot height in cm (only if save = TRUE)
                                     resolution = 300, #plot resolution in dpi (only if save = TRUE)
                                     verbose = TRUE #print messages
) {

  # Validate input
  if (!is.null(dapc.results$dapc_results)) dapc.results <- dapc.results$dapc_results
  if (is.null(dapc.results) || !is.list(dapc.results)) stop("dapc.results must be a DAPC list or a list with $dapc_results")
  if (is.null(dapc.results$var.contr)) stop("dapc.results$var.contr is missing")
  if (!is.character(group.colors) || length(group.colors) != 2L || any(!nzchar(group.colors))) stop("group.colors must be a character vector of length 2")
  if (!is.data.frame(predictor.data) && !is.matrix(predictor.data)) stop("predictor.data must be a data.frame or matrix")
  if (missing(species.labels)) stop("species.labels must be provided")
  if (!is.numeric(N.top.variables) || length(N.top.variables) != 1L || !is.finite(N.top.variables) || N.top.variables < 1) stop("N.top.variables must be a single positive integer (recommended: 6)")
  if (!is.numeric(alpha.density) || length(alpha.density) != 1L || !is.finite(alpha.density) || alpha.density < 0 || alpha.density > 1) stop("alpha.density must be between 0 and 1 (recommended: 0.75)")
  if (!is.numeric(axis.ticks.font.size) || length(axis.ticks.font.size) != 1L || !is.finite(axis.ticks.font.size) || axis.ticks.font.size <= 0) stop("axis.ticks.font.size must be a single positive number (recommended: 7)")
  if (!is.numeric(axis.labels.font.size) || length(axis.labels.font.size) != 1L || !is.finite(axis.labels.font.size) || axis.labels.font.size <= 0) stop("axis.labels.font.size must be a single positive number (recommended: 9.1)")
  if (!is.numeric(variable.font.size) || length(variable.font.size) != 1L || !is.finite(variable.font.size) || variable.font.size <= 0) stop("variable.font.size must be a single positive number (recommended: 9.1)")
  if (!is.null(legend.label) && (!is.character(legend.label) || length(legend.label) != 1L || !nzchar(legend.label))) stop("legend.label must be NULL or a single non-empty character string (recommended: 'Species')")
  if (!(is.character(legend.position) || is.numeric(legend.position))) stop("legend.position must be character or numeric")
  if (is.character(legend.position) && !legend.position %in% c("right", "left", "top", "bottom", "none")) stop("legend.position must be 'right', 'left', 'top', 'bottom', or 'none'")
  if (is.numeric(legend.position) && (length(legend.position) != 2L || any(!is.finite(legend.position)))) stop("legend.position numeric must be c(x, y) with finite values")
  if (!is.numeric(legend.title.font.size) || length(legend.title.font.size) != 1L || !is.finite(legend.title.font.size) || legend.title.font.size <= 0) stop("legend.title.font.size must be a single positive number (recommended: 9.1)")
  if (!is.numeric(legend.text.font.size) || length(legend.text.font.size) != 1L || !is.finite(legend.text.font.size) || legend.text.font.size <= 0) stop("legend.text.font.size must be a single positive number (recommended: 9.1)")
  if (!is.logical(legend.text.italics) || length(legend.text.italics) != 1L) stop("legend.text.italics must be TRUE or FALSE (recommended: FALSE)")
  if (!is.numeric(legend.symbol.size) || length(legend.symbol.size) != 1L || !is.finite(legend.symbol.size) || legend.symbol.size <= 0) stop("legend.symbol.size must be a single positive number (pt) (recommended: 15)")
  if (!is.logical(add.title) || length(add.title) != 1L) stop("add.title must be TRUE or FALSE")
  if (!is.numeric(title.font.size) || length(title.font.size) != 1L || !is.finite(title.font.size) || title.font.size <= 0) stop("title.font.size must be a single positive number (recommended: 9.1)")
  if (!is.null(title.text) && (!is.character(title.text) || length(title.text) != 1L)) stop("title.text must be NULL or a single character string")
  if (!is.logical(show.plot) || length(show.plot) != 1L) stop("show.plot must be TRUE or FALSE")
  if (!is.logical(save) || length(save) != 1L) stop("save must be TRUE or FALSE")
  if (!is.logical(overwrite) || length(overwrite) != 1L) stop("overwrite must be TRUE or FALSE")
  if (!is.character(filename) || length(filename) != 1L || !nzchar(filename)) stop("filename must be non-empty (recommended: 'top_DAPC_predictors')")
  if (!is.null(output.dir) && !is.character(output.dir)) stop("output.dir must be NULL or character")
  type <- tolower(type)
  if (!type %in% c("png", "svg", "jpg")) stop("type must be one of: png, svg, jpg")
  if (!is.numeric(width) || length(width) != 1L || !is.finite(width) || width <= 0) stop("width must be positive (cm) (recommended: 16)")
  if (!is.numeric(height) || length(height) != 1L || !is.finite(height) || height <= 0) stop("height must be positive (cm) (recommended: 10)")
  if (!is.numeric(resolution) || length(resolution) != 1L || !is.finite(resolution) || resolution <= 0) stop("resolution must be positive (dpi) (recommended: 300)")
  if (!is.logical(verbose) || length(verbose) != 1L) stop("verbose must be TRUE or FALSE (recommended: TRUE)")

  # Validate DAPC input
  if (!is.null(dapc.results$dapc_results)) dapc.results <- dapc.results$dapc_results
  if (is.null(dapc.results) || !is.list(dapc.results)) stop("dapc.results must be a DAPC list or a list with $dapc_results")
  if (is.null(dapc.results$var.contr)) stop("dapc.results$var.contr is missing")

  # Validate predictor data and species labels
  if (!is.data.frame(predictor.data)) predictor.data <- as.data.frame(predictor.data)
  if (nrow(predictor.data) < 2) stop("predictor.data must have at least 2 rows")
  if (nrow(predictor.data) != length(species.labels)) stop("Length of species.labels must match number of rows in predictor.data")
  species.labels <- factor(species.labels)
  if (nlevels(species.labels) != 2L) stop("species.labels must contain exactly two levels")
  if (!is.null(names(group.colors))) {
    if (!all(levels(species.labels) %in% names(group.colors))) stop("names(group.colors) must include all group levels: ", paste(levels(species.labels), collapse = ", "))
    species.labels <- factor(species.labels, levels = names(group.colors))
  } else {
    group.colors <- stats::setNames(group.colors[seq_len(nlevels(species.labels))], levels(species.labels))
  }

  # Rank variables by LD1 contribution
  var_contr <- dapc.results$var.contr
  if (is.matrix(var_contr) || is.data.frame(var_contr)) {
    if (ncol(var_contr) < 1) stop("var.contr has no columns")
    contrib_values <- as.numeric(var_contr[, 1])
    if (is.null(rownames(var_contr))) stop("var.contr must have rownames")
    names(contrib_values) <- rownames(var_contr)
  } else {
    if (is.null(names(var_contr))) stop("var.contr vector must have names")
    contrib_values <- as.numeric(var_contr)
    names(contrib_values) <- names(var_contr)
  }
  if (any(!is.finite(contrib_values))) stop("var.contr contains non-finite values")
  variable_contributions <- sort(contrib_values, decreasing = TRUE)

  # Select top N variables present in predictor.data
  candidate_names <- names(variable_contributions)[seq_len(min(N.top.variables, length(variable_contributions)))]
  candidate_names <- candidate_names[candidate_names %in% colnames(predictor.data)]
  if (!length(candidate_names)) stop("None of the top variables are present in predictor.data")

  # Keep numeric and finite columns
  is_numeric <- vapply(predictor.data[, candidate_names, drop = FALSE], is.numeric, logical(1))
  if (any(!is_numeric)) {
    dropped <- setdiff(candidate_names, candidate_names[is_numeric])
    candidate_names <- candidate_names[is_numeric]
    if (verbose && length(dropped)) warning("Dropping non-numeric variables: ", paste(dropped, collapse = ", "))
  }
  if (!length(candidate_names)) stop("No numeric top variables available for density plots")
  has_finite <- vapply(predictor.data[, candidate_names, drop = FALSE], function(z) any(is.finite(z)), logical(1))
  if (any(!has_finite)) {
    dropped <- setdiff(candidate_names, candidate_names[has_finite])
    candidate_names <- candidate_names[has_finite]
    if (verbose && length(dropped)) warning("Dropping all-NA/non-finite variables: ", paste(dropped, collapse = ", "))
  }
  if (!length(candidate_names)) stop("No usable variables left after filtering")

  # Determine facet layout
  n_panels <- length(candidate_names)
  if (is.null(plot.nrow) && is.null(plot.ncol)) {
    plot.nrow <- floor(sqrt(n_panels))
    if (plot.nrow < 1L) plot.nrow <- 1L
    plot.ncol <- ceiling(n_panels / plot.nrow)
  } else if (!is.null(plot.nrow) && is.null(plot.ncol)) {
    if (!is.numeric(plot.nrow) || length(plot.nrow) != 1L || !is.finite(plot.nrow) || plot.nrow < 1) stop("plot.nrow must be NULL or single integer >=1")
    plot.nrow <- as.integer(round(plot.nrow))
    plot.ncol <- ceiling(n_panels / plot.nrow)
  } else if (is.null(plot.nrow) && !is.null(plot.ncol)) {
    if (!is.numeric(plot.ncol) || length(plot.ncol) != 1L || !is.finite(plot.ncol) || plot.ncol < 1) stop("plot.ncol must be NULL or single integer >=1")
    plot.ncol <- as.integer(round(plot.ncol))
    plot.nrow <- ceiling(n_panels / plot.ncol)
  } else {
    if (!is.numeric(plot.nrow) || length(plot.nrow) != 1L || !is.finite(plot.nrow) || plot.nrow < 1) stop("plot.nrow must be single integer >=1")
    if (!is.numeric(plot.ncol) || length(plot.ncol) != 1L || !is.finite(plot.ncol) || plot.ncol < 1) stop("plot.ncol must be single integer >=1")
    plot.nrow <- as.integer(round(plot.nrow))
    plot.ncol <- as.integer(round(plot.ncol))
  }

  # Prepare tidy data
  predictor_data_long <- tidyr::pivot_longer(
    cbind(predictor.data[, candidate_names, drop = FALSE], group = factor(species.labels)),
    cols = tidyselect::any_of(candidate_names),
    names_to = "Variable",
    values_to = "Value"
  )
  predictor_data_long$Variable <- factor(predictor_data_long$Variable, levels = candidate_names)
  group_levels <- levels(species.labels)
  if (is.null(names(group.colors))) names(group.colors) <- group_levels

  # Build plot
  plot_object <- ggplot(predictor_data_long, aes(x = Value, fill = group)) +
    geom_density(alpha = alpha.density, colour = "black", show.legend = TRUE) +
    scale_fill_manual(values = group.colors[group_levels], breaks = group_levels, name = legend.label) +
    guides(fill = guide_legend(keyheight = grid::unit(legend.symbol.size, "pt"), keywidth = grid::unit(legend.symbol.size, "pt"))) +
    labs(x = "Environmental variable value", y = "Density", fill = legend.label) +
    facet_wrap(~ Variable, scales = "free", nrow = plot.nrow, ncol = plot.ncol) +
    theme_classic() +
    theme(axis.ticks = element_line(colour = "black"),
          axis.line = element_line(colour = "black"),
          strip.background = element_rect(fill = "white", colour = "black"),
          strip.text = element_text(size = variable.font.size, colour = "black", face = "bold"),
          axis.text = element_text(size = axis.ticks.font.size, colour = "black"),
          axis.title = element_text(size = axis.labels.font.size, colour = "black", face = "bold"),
          legend.title = element_text(size = legend.title.font.size, colour = "black", face = "bold"),
          legend.text = element_text(size = legend.text.font.size, colour = "black", face = if (legend.text.italics) "italic" else "plain"),
          legend.key.size = grid::unit(legend.symbol.size, "pt"),
          legend.position = legend.position,
          plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = title.font.size))

  # Add title
  if (isTRUE(add.title)) {
    title_final <- if (is.null(title.text)) paste0("Distributions of ", length(candidate_names), " most contributing variables") else title.text
    plot_object <- plot_object + labs(title = title_final)
  }

  # Save plot
  if (isTRUE(save)) {
    out_dir <- if (is.null(output.dir) || !nzchar(output.dir)) getwd() else output.dir
    if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      if (verbose) message("Created directory: ", out_dir) }
    file_out <- file.path(out_dir, paste0(filename, ".", type))
    if (file.exists(file_out) && !overwrite) stop("File already exists: ", file_out, " (set overwrite = TRUE to replace)")
    device_for_ggsave <- if (type == "jpg") "jpeg" else type
    ggsave(filename = file_out, plot = plot_object, width = width, height = height, units = "cm", dpi = resolution, device = device_for_ggsave)
    if (verbose) message("Plot saved as ", file_out)
  }

  # Return plot
  if (show.plot) return(plot_object) else invisible(plot_object)
}


## Function to plot occurrences on a map
#' Plot occurrence records on a map
#'
#' Plot occurrence coordinates on a base map, optionally with background points
#' and map annotations.
#'
#' @param coordinates A `data.frame` or matrix containing occurrence coordinates
#'   to plot.
#' @param latitude.col A single character string giving the latitude column name
#'   in `coordinates` and, when supplied, `background.coords`
#'   (default: `"Latitude"`).
#' @param longitude.col A single character string giving the longitude column
#'   name in `coordinates` and, when supplied, `background.coords`
#'   (default: `"Longitude"`).
#' @param latitude.buffer.range A single non-negative numeric value giving the
#'   buffer added to the plotted latitude extent (default: `2`).
#' @param longitude.buffer.range A single non-negative numeric value giving the
#'   buffer added to the plotted longitude extent (default: `2`).
#' @param group.colors Character vector of two colors used for the occurrence
#'   groups. Can be named with group levels to enforce a specific group-color
#'   mapping (default: `c("#00005A", "darkgrey")`).
#' @param group.labels Factor or character vector giving group membership for
#'   occurrence points. Exactly two groups are required.
#' @param CRS A single character string giving the coordinate reference system
#'   used for plotting (default: `"EPSG:4326"`).
#' @param point.size A single positive numeric value giving the occurrence point
#'   size (default: `1.22`).
#' @param point.alpha A single numeric value between `0` and `1` controlling
#'   occurrence point transparency (default: `0.9`).
#' @param point.border.color A single character string giving the occurrence
#'   point border color (default: `"black"`).
#' @param plot.background.points Logical; if `TRUE`, background points are drawn
#'   on the map (default: `TRUE`).
#' @param background.coords Optional `data.frame` or matrix containing background
#'   point coordinates. Required when `plot.background.points = TRUE`
#'   (default: `NULL`).
#' @param background.group.labels Optional factor or character vector giving
#'   group membership for background points. If supplied, its levels must match
#'   `group.labels` (default: `NULL`).
#' @param background.point.size A single positive numeric value giving the
#'   background point size (default: `0.22`).
#' @param background.point.col A single character string giving the background
#'   point color when `background.group.labels = NULL` (default: `"grey60"`).
#' @param background.point.alpha A single numeric value between `0` and `1`
#'   controlling background point transparency (default: `0.9`).
#' @param axis.numbers.size A single positive numeric value giving the coordinate
#'   tick-label font size (default: `9.1`).
#' @param add.USA.states Logical; if `TRUE`, USA state borders are added
#'   (default: `TRUE`).
#' @param add.USA.counties Logical; if `TRUE`, USA county borders are added
#'   (default: `FALSE`).
#' @param country.lwd A single non-negative numeric value giving the country
#'   border line width (default: `1`).
#' @param USA.state.lwd A single non-negative numeric value giving the USA state
#'   border line width (default: `1`).
#' @param USA.county.lwd A single non-negative numeric value giving the USA
#'   county border line width (default: `0.3`).
#' @param map.background A single character string giving the fill color used for
#'   the map background (default: `"lightgrey"`).
#' @param north.arrow.position Numeric vector of length two giving the relative
#'   x/y position of the north arrow within the plotted extent
#'   (default: `c(0.03, 0.88)`).
#' @param north.arrow.length A single positive numeric value giving the north
#'   arrow length in map units (default: `0.7`).
#' @param north.arrow.lwd A single non-negative numeric value giving the north
#'   arrow line width (default: `2`).
#' @param north.arrow.font.size A single positive numeric value giving the font
#'   size of the north-arrow label (default: `9.1`).
#' @param north.arrow.N.position A single numeric value giving the vertical offset
#'   of the `"N"` label above the arrow tip (default: `0.3`).
#' @param scale.position Numeric vector of length two giving the relative x/y
#'   position of the scale bar within the plotted extent
#'   (default: `c(0.03, 0.05)`).
#' @param scale.size A single positive numeric value controlling the scale-bar
#'   width (default: `0.16`).
#' @param scale.font.size A single positive numeric value giving the scale-bar
#'   text font size (default: `7`).
#' @param legend.position A single character string giving the legend position
#'   used by base R graphics (default: `"topright"`).
#' @param legend.font.size A single positive numeric value giving the legend text
#'   font size (default: `9.1`).
#' @param legend.group.names Optional character vector giving legend labels for
#'   the two occurrence groups. If `NULL`, levels of `group.labels` are used
#'   (default: `NULL`).
#' @param legend.box Logical; if `TRUE`, a legend box is drawn
#'   (default: `TRUE`).
#' @param legend.text.italics Logical; if `TRUE`, legend text is italicized
#'   (default: `FALSE`).
#' @param legend.symbol.size A single positive numeric value giving the legend
#'   symbol size (default: `1.5`).
#' @param show.plot Logical; if `TRUE`, the map is drawn on the active graphics
#'   device (default: `TRUE`).
#' @param save Logical; if `TRUE`, the map is saved to disk
#'   (default: `FALSE`).
#' @param overwrite Logical; if `TRUE`, an existing output file is overwritten
#'   when `save = TRUE` (default: `TRUE`).
#' @param filename A single character string giving the output filename without
#'   extension when `save = TRUE` (default: `"occurrence_map"`).
#' @param output.dir Optional character string giving the directory for saved
#'   plots when `save = TRUE` (default: `NULL`; if `NULL`, the current working
#'   directory is used).
#' @param type A single character string giving the output file type. One of
#'   `"png"`, `"svg"`, or `"jpg"` (default: `"svg"`).
#' @param width A single positive numeric value giving plot width in centimeters
#'   when `save = TRUE` (default: `15`).
#' @param height A single positive numeric value giving plot height in
#'   centimeters when `save = TRUE` (default: `20`).
#' @param resolution A single positive numeric value giving plot resolution in
#'   dpi when saving raster formats (default: `300`).
#' @param verbose Logical; if `TRUE`, messages are printed when saving
#'   (default: `TRUE`).
#'
#' @details
#' Occurrence maps provide the geographic context for niche-divergence analyses
#' by showing where occurrence records and, optionally, accessible background
#' points are located in geographic space. This is useful because environmental
#' niche divergence is tested in environmental space, but the interpretation of
#' sampling, accessible areas, spatial clustering, and background availability
#' depends on the geographic distribution of the data.
#'
#' The function is designed for occurrence and background coordinates in decimal
#' degrees. Map extents are expanded by the latitude and longitude buffer
#' arguments so that plotted points are not placed directly on the map boundary.
#' Optional state and county overlays can be used when the study area falls within
#' the United States. The plot may require repeated adjustment of mapping
#' parameters to produce an appropriate layout, including legend placement, north
#' arrow position, scale-bar position, point sizes, and map extent.
#'
#' @return No return value. The function is called for its side effects: drawing
#'   a base R map and, when `save = TRUE`, writing the plot to disk.
#'
#' @rawNamespace export(plot.occurrences.map)
plot.occurrences.map <- function(coordinates, #data.frame/matrix with coordinates
                                 latitude.col = "Latitude", #name of latitude column
                                 longitude.col = "Longitude", #name of longitude column
                                 latitude.buffer.range = 2, #buffer added to latitude extent
                                 longitude.buffer.range = 2, #buffer added to longitude extent
                                 group.colors = c("#00005A", "darkgrey"), #two colors for groups
                                 group.labels, #factor/character (two levels) for points
                                 CRS = "EPSG:4326", #coordinate reference system for plotting (optional projection)
                                 point.size = 1.22, #occurrence point size
                                 point.alpha = 0.9, #occurrence point transparency (0-1)
                                 point.border.color = "black", #occurrence point border color
                                 plot.background.points = TRUE, #draw background points
                                 background.coords = NULL, #data.frame/matrix for background coordinates (only if plot.background.points = TRUE)
                                 background.group.labels = NULL, #optional factor/character for background groups (only if plot.background.points = TRUE)
                                 background.point.size = 0.22, #background point size (only if plot.background.points = TRUE)
                                 background.point.col = "grey60", #background point color (only if plot.background.points = TRUE)
                                 background.point.alpha = 0.9, #background point transparency (only if plot.background.points = TRUE)
                                 axis.numbers.size = 9.1, #coordinate tick font size
                                 add.USA.states = TRUE, #add US states
                                 add.USA.counties = FALSE, #add US counties
                                 country.lwd = 1, #country border line width
                                 USA.state.lwd = 1, #state border line width
                                 USA.county.lwd = 0.3, #county border line width
                                 map.background = "lightgrey", #background color of map
                                 north.arrow.position = c(0.03, 0.88), #relative (x,y)
                                 north.arrow.length = 0.7, #arrow length
                                 north.arrow.lwd = 2, #arrow line width
                                 north.arrow.font.size = 9.1, #font size of north arrow "N"
                                 north.arrow.N.position = 0.3, #offset above arrow tip for "N"
                                 scale.position = c(0.03, 0.05), #relative (x, y)
                                 scale.size = 0.16, #scale width
                                 scale.font.size = 7, #scale font size
                                 legend.position = "topright", #legend position
                                 legend.font.size = 9.1, #legend font size
                                 legend.group.names = NULL, #labels for two groups
                                 legend.box = TRUE, #draw legend box
                                 legend.text.italics = FALSE, #italicize legend text
                                 legend.symbol.size = 1.5, #legend symbol size
                                 show.plot = TRUE, #whether to print plot to console
                                 save = FALSE, #save plot
                                 overwrite = TRUE, #overwrite existing file (only if save = TRUE)
                                 filename = "occurrence_map", #filename for saved plot (no extension)
                                 output.dir = NULL, #directory for output (only if save = TRUE)
                                 type = "svg", #plot type: "png","svg","jpg" (only if save = TRUE)
                                 width = 15, #plot width in cm (only if save = TRUE)
                                 height = 20, #plot height in cm (only if save = TRUE)
                                 resolution = 300, #plot resolution in dpi (only if save = TRUE)
                                 verbose = TRUE #print messages when saving
) {

  # Validate input
  if (!is.data.frame(coordinates) && !is.matrix(coordinates)) stop("coordinates must be a data.frame or matrix")
  if (is.matrix(coordinates)) coordinates <- as.data.frame(coordinates)
  if (!all(c(latitude.col, longitude.col) %in% names(coordinates))) stop("coordinates must contain specified latitude and longitude columns")
  lat_vec <- as.numeric(coordinates[[latitude.col]])
  lon_vec <- as.numeric(coordinates[[longitude.col]])
  if (any(!is.finite(lat_vec)) || any(!is.finite(lon_vec))) stop("Latitude or Longitude contains non-finite values")
  if (is.null(group.labels)) stop("group.labels must be provided")
  if (length(group.labels) != nrow(coordinates)) stop("Length of group.labels must match number of rows in coordinates")
  group.labels <- factor(group.labels)
  if (nlevels(group.labels) != 2L) stop("group.labels must have exactly two levels")
  if (!is.character(group.colors) || length(group.colors) != 2L || any(!nzchar(group.colors))) stop("group.colors must be a character vector of length 2 (recommended: c('#00005A','darkgrey'))")
  if (!is.logical(show.plot) || length(show.plot) != 1L) stop("show.plot must be TRUE or FALSE")
  if (!is.logical(plot.background.points) || length(plot.background.points) != 1L) stop("plot.background.points must be TRUE or FALSE (recommended: TRUE)")
  if (isTRUE(plot.background.points)) {
    if (is.null(background.coords)) stop("background.coords must be provided when plot.background.points = TRUE")
    if (!is.data.frame(background.coords) && !is.matrix(background.coords)) stop("background.coords must be a data.frame or matrix")
    if (is.matrix(background.coords)) background.coords <- as.data.frame(background.coords)
    if (!all(c(latitude.col, longitude.col) %in% names(background.coords))) stop("background.coords must contain specified latitude and longitude columns")
    bg_lat <- as.numeric(background.coords[[latitude.col]])
    bg_lon <- as.numeric(background.coords[[longitude.col]])
    if (any(!is.finite(bg_lat)) || any(!is.finite(bg_lon))) stop("background.coords contains non-finite coordinates")
    if (!is.numeric(background.point.size) || length(background.point.size) != 1L || !is.finite(background.point.size) || background.point.size <= 0) stop("background.point.size must be a single positive number (recommended: 0.22)")
    if (!is.numeric(background.point.alpha) || length(background.point.alpha) != 1L || !is.finite(background.point.alpha) || background.point.alpha < 0 || background.point.alpha > 1) stop("background.point.alpha must be in between 0 and 1 (recommended: 0.9)")
    if (!is.null(background.group.labels)) {
      if (length(background.group.labels) != nrow(background.coords)) stop("Length of background.group.labels must equal nrow(background.coords)")
      background.group.labels <- factor(background.group.labels, levels = levels(group.labels))
      if (any(is.na(background.group.labels))) stop("background.group.labels must match group.labels levels")
    }
  }
  if (!is.character(CRS) || length(CRS) != 1L || !nzchar(CRS)) stop("CRS must be a single non-empty character string (recommended: 'EPSG:4326')")
  if (!is.numeric(point.size) || length(point.size) != 1L || !is.finite(point.size) || point.size <= 0) stop("point.size must be a single positive number (recommended: 1.22)")
  if (!is.numeric(point.alpha) || length(point.alpha) != 1L || !is.finite(point.alpha) || point.alpha < 0 || point.alpha > 1) stop("point.alpha must be in between 0 and 1 (recommended: 0.9)")
  if (!is.numeric(latitude.buffer.range) || length(latitude.buffer.range) != 1L || !is.finite(latitude.buffer.range) || latitude.buffer.range < 0) stop("latitude.buffer.range must be a single non-negative number (recommended: 2)")
  if (!is.numeric(longitude.buffer.range) || length(longitude.buffer.range) != 1L || !is.finite(longitude.buffer.range) || longitude.buffer.range < 0) stop("longitude.buffer.range must be a single non-negative number (recommended: 2)")
  if (!is.numeric(axis.numbers.size) || length(axis.numbers.size) != 1L || !is.finite(axis.numbers.size) || axis.numbers.size <= 0) stop("axis.numbers.size must be a single positive number (recommended: 9.1)")
  if (!is.numeric(country.lwd) || length(country.lwd) != 1L || !is.finite(country.lwd) || country.lwd < 0) stop("country.lwd must be a single non-negative number (recommended: 1)")
  if (!is.logical(add.USA.states) || length(add.USA.states) != 1L) stop("add.USA.states must be TRUE or FALSE")
  if (!is.logical(add.USA.counties) || length(add.USA.counties) != 1L) stop("add.USA.counties must be TRUE or FALSE")
  if (!is.numeric(USA.state.lwd) || length(USA.state.lwd) != 1L || !is.finite(USA.state.lwd) || USA.state.lwd < 0) stop("USA.state.lwd must be a single non-negative number (recommended: 1)")
  if (!is.numeric(USA.county.lwd) || length(USA.county.lwd) != 1L || !is.finite(USA.county.lwd) || USA.county.lwd < 0) stop("USA.county.lwd must be a single non-negative number (recommended: 0.3)")
  if (!is.numeric(north.arrow.length) || length(north.arrow.length) != 1L || !is.finite(north.arrow.length) || north.arrow.length <= 0) stop("north.arrow.length must be a single positive number (recommended: 0.7)")
  if (!is.numeric(north.arrow.lwd) || length(north.arrow.lwd) != 1L || !is.finite(north.arrow.lwd) || north.arrow.lwd < 0) stop("north.arrow.lwd must be a single non-negative number (recommended: 2)")
  if (!is.numeric(north.arrow.font.size) || length(north.arrow.font.size) != 1L || !is.finite(north.arrow.font.size) || north.arrow.font.size <= 0) stop("north.arrow.font.size must be a single positive number (recommended: 9.1)")
  if (!is.numeric(scale.size) || length(scale.size) != 1L || !is.finite(scale.size) || scale.size <= 0) stop("scale.size must be a single positive number (recommended: 0.16)")
  if (!is.numeric(scale.font.size) || length(scale.font.size) != 1L || !is.finite(scale.font.size) || scale.font.size <= 0) stop("scale.font.size must be a single positive number (recommended: 7)")
  if (!is.character(legend.position) || length(legend.position) != 1L || !nzchar(legend.position)) stop("legend.position must be a non-empty character string (recommended: 'topright')")
  if (!is.numeric(legend.font.size) || length(legend.font.size) != 1L || !is.finite(legend.font.size) || legend.font.size <= 0) stop("legend.font.size must be a single positive number (recommended: 9.1)")
  if (!is.logical(legend.box) || length(legend.box) != 1L) stop("legend.box must be TRUE or FALSE")
  if (!is.logical(legend.text.italics) || length(legend.text.italics) != 1L) stop("legend.text.italics must be TRUE or FALSE")
  if (!is.numeric(legend.symbol.size) || length(legend.symbol.size) != 1L || !is.finite(legend.symbol.size) || legend.symbol.size <= 0) stop("legend.symbol.size must be a single positive number (pt) (recommended: 1.5)")
  if (!is.logical(save) || length(save) != 1L) stop("save must be TRUE or FALSE")
  if (!is.logical(overwrite) || length(overwrite) != 1L) stop("overwrite must be TRUE or FALSE")
  if (!is.character(filename) || length(filename) != 1L || !nzchar(filename)) stop("filename must be a non-empty string without file extension")
  type <- tolower(type)
  if (!type %in% c("png","svg","jpg")) stop("type must be one of: 'png', 'svg', 'jpg'")
  if (!is.numeric(width) || length(width) != 1L || !is.finite(width) || width <= 0) stop("width must be a single positive number (cm)")
  if (!is.numeric(height) || length(height) != 1L || !is.finite(height) || height <= 0) stop("height must be a single positive number (cm)")
  if (!is.numeric(resolution) || length(resolution) != 1L || !is.finite(resolution) || resolution <= 0) stop("resolution must be a single positive number (dpi)")
  if (!is.logical(verbose) || length(verbose) != 1L) stop("verbose must be TRUE or FALSE")

  # Create function to draw map on current device
  draw.map <- function() {
    graphics::par(mar = c(1, 1, 1, 1), oma = c(0, 0, 0, 0))
    graphics::plot.new()
    graphics::plot.window(xlim = c(lon_min, lon_max),
                          ylim = c(lat_min, lat_max))
    maps::map("world",
              fill = TRUE,
              col = map.background,
              border = "black",
              lwd = country.lwd,
              add = TRUE,
              xlim = c(lon_min, lon_max),
              ylim = c(lat_min, lat_max))
    if (isTRUE(add.USA.counties)) maps::map("county",
                                            add = TRUE,
                                            col = "grey",
                                            lwd = USA.county.lwd,
                                            xlim = c(lon_min, lon_max),
                                            ylim = c(lat_min, lat_max))
    if (isTRUE(add.USA.states)) maps::map("state",
                                          add = TRUE,
                                          col = "black",
                                          lwd = USA.state.lwd)
    graphics::box()
    if (isTRUE(plot.background.points)) {
      if (!is.null(background.group.labels)) {
        col_map_bg <- stats::setNames(group.colors, levels(group.labels))
        bg_cols <- grDevices::adjustcolor(col_map_bg[as.character(background.group.labels)],
                                          alpha.f = background.point.alpha)
        graphics::points(x = bg_lon,
                         y = bg_lat,
                         pch = 16,
                         col = bg_cols,
                         cex = background.point.size)
      } else {
        bg_col <- grDevices::adjustcolor(background.point.col,
                                         alpha.f = background.point.alpha)
        graphics::points(x = bg_lon,
                         y = bg_lat,
                         pch = 16,
                         col = bg_col,
                         cex = background.point.size)
      }
    }
    col_map <- stats::setNames(group.colors, levels(group.labels))
    fill_cols <- grDevices::adjustcolor(col_map[as.character(group.labels)],
                                        alpha.f = point.alpha)
    graphics::points(x = lon_vec,
                     y = lat_vec,
                     pch = 21,
                     bg = fill_cols,
                     col = point.border.color,
                     cex = point.size)
    legend_labels <- if (is.null(legend.group.names)) levels(group.labels) else legend.group.names
    legend_box_type <- if (isTRUE(legend.box)) "o" else "n"
    graphics::legend(x = legend.position,
                     legend = legend_labels,
                     pch = 21,
                     cex = legend.font.size / graphics::par("ps"),
                     pt.cex = legend.symbol.size,
                     pt.bg = group.colors,
                     text.font = if (legend.text.italics) 3 else 1,
                     bty = legend_box_type)
    scale_x <- scale.position[1] * (lon_max - lon_min) + lon_min
    scale_y <- scale.position[2] * (lat_max - lat_min) + lat_min
    maps::map.scale(x = scale_x,
                    y = scale_y,
                    cex = scale.font.size / graphics::par("ps"),
                    relwidth = scale.size,
                    ratio = FALSE)
    north_x <- north.arrow.position[1] * (lon_max - lon_min) + lon_min
    north_y <- north.arrow.position[2] * (lat_max - lat_min) + lat_min
    graphics::arrows(x0 = north_x,
                     y0 = north_y,
                     x1 = north_x,
                     y1 = north_y + north.arrow.length,
                     length = 0.13,
                     col = "black",
                     lwd = north.arrow.lwd)
    graphics::text(x = north_x,
                   y = north_y + north.arrow.length + north.arrow.N.position,
                   labels = "N",
                   cex = north.arrow.font.size / graphics::par("ps"),
                   col = "black")
  }

  # Coordinate system awareness
  if (CRS != "EPSG:4326") {
    sf_points <- sf::st_as_sf(coordinates, coords = c(longitude.col, latitude.col), crs = 4326)
    sf_points <- sf::st_transform(sf_points, crs = CRS)
    coords_mat <- sf::st_coordinates(sf_points)
    lon_vec <- coords_mat[, 1]
    lat_vec <- coords_mat[, 2]
    if (isTRUE(plot.background.points) && !is.null(background.coords)) {
      bg_sf <- sf::st_as_sf(background.coords, coords = c(longitude.col, latitude.col), crs = 4326)
      bg_sf <- sf::st_transform(bg_sf, crs = CRS)
      bg_coords <- sf::st_coordinates(bg_sf)
      bg_lon <- bg_coords[, 1]
      bg_lat <- bg_coords[, 2]
    }
  }

  # Set extent
  if (isTRUE(plot.background.points)) {
    lon_min <- min(c(lon_vec, bg_lon)) - longitude.buffer.range
    lon_max <- max(c(lon_vec, bg_lon)) + longitude.buffer.range
    lat_min <- min(c(lat_vec, bg_lat)) - latitude.buffer.range
    lat_max <- max(c(lat_vec, bg_lat)) + latitude.buffer.range
  } else {
    lon_min <- min(lon_vec) - longitude.buffer.range
    lon_max <- max(lon_vec) + longitude.buffer.range
    lat_min <- min(lat_vec) - latitude.buffer.range
    lat_max <- max(lat_vec) + latitude.buffer.range
  }
  if (lon_max <= lon_min) lon_max <- lon_min + 0.01
  if (lat_max <= lat_min) lat_max <- lat_min + 0.01

  # Open device if saving
  if (isTRUE(save)) {
    if (!is.null(output.dir) && !dir.exists(output.dir)) {
      dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
      if (verbose) message("Created directory: ", output.dir)
    }
    file_out <- if (is.null(output.dir)) paste0(filename, ".", type) else file.path(output.dir, paste0(filename, ".", type))
    if (file.exists(file_out) && !overwrite) stop("File already exists: ", file_out)
    if (type == "svg") svglite::svglite(file_out,
                                        width = width / 2.54,
                                        height = height / 2.54,
                                        bg = "white")
    if (type == "png") grDevices::png(file_out,
                                      width = width,
                                      height = height,
                                      units = "cm",
                                      res = resolution,
                                      bg = "white")
    if (type == "jpg") grDevices::jpeg(file_out,
                                       width = width,
                                       height = height,
                                       units = "cm",
                                       res = resolution,
                                       bg = "white")
  }

  # Draw base map on current device
  draw.map()

  # Save plot
  if (isTRUE(save)) {
    grDevices::dev.off()
    if (verbose) message("Plot saved as ", file_out)
  }

  # Return plot
  if (isTRUE(show.plot)) draw.map()
}


## Function to map variable names to descriptive names
#' Map environmental variable names to readable labels
#'
#' Replace internal abbreviated environmental variable names with readable
#' full or shortened labels. The function can map column names of a `data.frame`,
#' names of a named vector, or row names of a matrix.
#'
#' @param input.data A `data.frame`, named vector, or matrix with row names. For
#'   a `data.frame`, column names are mapped. For a named vector, vector names are
#'   mapped. For a matrix, row names are mapped.
#' @param name.length Character string specifying the label length to use. One of
#'   `"full"` for descriptive labels with units where available, or `"short"` for
#'   abbreviated labels (default: `"full"`).
#' @param recognize.transformations Logical; if `TRUE`, transformation suffixes
#'   such as `"_log"`, `"_sqrt"`, or `"_arcsine_sqrt"` are detected and preserved
#'   in the mapped label (default: `TRUE`).
#'
#' @details
#' This utility is used to convert compact environmental variable names into
#' labels that are easier to interpret in tables and plots. Full labels are
#' intended for descriptive output where units and temporal information are
#' useful, whereas short labels are intended for compact visualizations or
#' summaries where space is limited.
#'
#' When `recognize.transformations = TRUE`, the function first identifies known
#' transformation suffixes and maps the base environmental variable name. The
#' transformation suffix is then appended to the readable label, allowing
#' transformed variables to remain identifiable while still using descriptive
#' environmental names. Variable names that are not found in the mapping table are
#' returned unchanged.
#'
#' @return The same object type as `input.data`. For a `data.frame`, column names
#'   are replaced by mapped environmental variable labels. For a named vector,
#'   vector names are replaced. For a matrix, row names are replaced. Data values,
#'   row order, and column order are unchanged.
#'
#' @export
map.env.variable.names <- function(input.data, #data frame or named numeric/vector to map variable names
                                   name.length = c("full", "short"), #use "full" for long descriptive names or "short" for abbreviated versions
                                   recognize.transformations = TRUE) { #whether to detect transformed variable names (e.g., MAT_log, CMD_sqrt)

  # Validate input
  if (!is.data.frame(input.data) &&
      !(is.vector(input.data) && !is.null(names(input.data))) &&
      !(is.matrix(input.data) && !is.null(rownames(input.data)))) {
    stop("Input must be either a data frame, a named vector, or a matrix with rownames")
  }
  if (!is.logical(recognize.transformations) || length(recognize.transformations) != 1) stop("recognize.transformations must be TRUE or FALSE")

  # Convert matrix input to named vector (map rownames)
  if (is.matrix(input.data) && !is.null(rownames(input.data))) {
    mapped.rownames <- vapply(rownames(input.data),
                              function(x) map.env.variable.names(setNames(1, x), name.length, recognize.transformations) |> names(),
                              FUN.VALUE = character(1))
    rownames(input.data) <- mapped.rownames
    return(input.data)
  }

  # Full mapping
  variable.mapping.full <- c(

    # ClimateNA
    "MAT" = "Mean annual temperature ( degreeC)",
    "MWMT" = "Mean warmest month temperature ( degreeC)",
    "MCMT" = "Mean coldest month temperature ( degreeC)",
    "TD" = "Temperature difference ( degreeC)",
    "MAP" = "Mean annual precipitation (mm)",
    "MSP" = "May to September precipitation (mm)",
    "AHM" = "Annual heat-moisture index",
    "SHM" = "Summer heat-moisture index",
    "bFFP" = "Beginning of frost-free period (day of year)",
    "eFFP" = "End of frost-free period (day of year)",
    "FFP" = "Frost-free period (days)",
    "CMD" = "Climatic moisture deficit (mm)",
    "CMI" = "Climatic moisture index",
    "DD_0" = "Degree days below 0  degreeC",
    "DD5" = "Degree days above 5  degreeC",
    "DD_18" = "Degree days below 18  degreeC",
    "DD18" = "Degree days above 18  degreeC",
    "DD1040" = "Degree days 10-40  degreeC",
    "EMT" = "Extreme minimum temperature ( degreeC)",
    "EXT" = "Extreme maximum temperature ( degreeC)",
    "Eref" = "Reference evaporation (mm)",
    "rsds" = "Mean solar radiation (MJ m^-^2 day^-^1)",
    "NFFD" = "Number of frost-free days (days)",
    "PAS" = "Precipitation as snow (mm)",
    "RH"  = "Mean relative humidity (%)",
    "Tmax01" = "Maximum temperature (Jan;  degreeC)",
    "Tmax02" = "Maximum temperature (Feb;  degreeC)",
    "Tmax03" = "Maximum temperature (Mar;  degreeC)",
    "Tmax04" = "Maximum temperature (Apr;  degreeC)",
    "Tmax05" = "Maximum temperature (May;  degreeC)",
    "Tmax06" = "Maximum temperature (Jun;  degreeC)",
    "Tmax07" = "Maximum temperature (Jul;  degreeC)",
    "Tmax08" = "Maximum temperature (Aug;  degreeC)",
    "Tmax09" = "Maximum temperature (Sep;  degreeC)",
    "Tmax10" = "Maximum temperature (Oct;  degreeC)",
    "Tmax11" = "Maximum temperature (Nov;  degreeC)",
    "Tmax12" = "Maximum temperature (Dec;  degreeC)",
    "Tmin01" = "Minimum temperature (Jan;  degreeC)",
    "Tmin02" = "Minimum temperature (Feb;  degreeC)",
    "Tmin03" = "Minimum temperature (Mar;  degreeC)",
    "Tmin04" = "Minimum temperature (Apr;  degreeC)",
    "Tmin05" = "Minimum temperature (May;  degreeC)",
    "Tmin06" = "Minimum temperature (Jun;  degreeC)",
    "Tmin07" = "Minimum temperature (Jul;  degreeC)",
    "Tmin08" = "Minimum temperature (Aug;  degreeC)",
    "Tmin09" = "Minimum temperature (Sep;  degreeC)",
    "Tmin10" = "Minimum temperature (Oct;  degreeC)",
    "Tmin11" = "Minimum temperature (Nov;  degreeC)",
    "Tmin12" = "Minimum temperature (Dec;  degreeC)",
    "Tave01" = "Mean temperature (Jan;  degreeC)",
    "Tave02" = "Mean temperature (Feb;  degreeC)",
    "Tave03" = "Mean temperature (Mar;  degreeC)",
    "Tave04" = "Mean temperature (Apr;  degreeC)",
    "Tave05" = "Mean temperature (May;  degreeC)",
    "Tave06" = "Mean temperature (Jun;  degreeC)",
    "Tave07" = "Mean temperature (Jul;  degreeC)",
    "Tave08" = "Mean temperature (Aug;  degreeC)",
    "Tave09" = "Mean temperature (Sep;  degreeC)",
    "Tave10" = "Mean temperature (Oct;  degreeC)",
    "Tave11" = "Mean temperature (Nov;  degreeC)",
    "Tave12" = "Mean temperature (Dec;  degreeC)",
    "PAS01" = "Precipitation as snow (Jan; mm)",
    "PAS02" = "Precipitation as snow (Feb; mm)",
    "PAS03" = "Precipitation as snow (Mar; mm)",
    "PAS04" = "Precipitation as snow (Apr; mm)",
    "PAS05" = "Precipitation as snow (May; mm)",
    "PAS06" = "Precipitation as snow (Jun; mm)",
    "PAS07" = "Precipitation as snow (Jul; mm)",
    "PAS08" = "Precipitation as snow (Aug; mm)",
    "PAS09" = "Precipitation as snow (Sep; mm)",
    "PAS10" = "Precipitation as snow (Oct; mm)",
    "PAS11" = "Precipitation as snow (Nov; mm)",
    "PAS12" = "Precipitation as snow (Dec; mm)",
    "PPT01" = "Precipitation (Jan; mm)",
    "PPT02" = "Precipitation (Feb; mm)",
    "PPT03" = "Precipitation (Mar; mm)",
    "PPT04" = "Precipitation (Apr; mm)",
    "PPT05" = "Precipitation (May; mm)",
    "PPT06" = "Precipitation (Jun; mm)",
    "PPT07" = "Precipitation (Jul; mm)",
    "PPT08" = "Precipitation (Aug; mm)",
    "PPT09" = "Precipitation (Sep; mm)",
    "PPT10" = "Precipitation (Oct; mm)",
    "PPT11" = "Precipitation (Nov; mm)",
    "PPT12" = "Precipitation (Dec; mm)",
    "rsds01" = "Surface downwelling shortwave radiation (Jan; MJ m^-^2 day^-^1)",
    "rsds02" = "Surface downwelling shortwave radiation (Feb; MJ m^-^2 day^-^1)",
    "rsds03" = "Surface downwelling shortwave radiation (Mar; MJ m^-^2 day^-^1)",
    "rsds04" = "Surface downwelling shortwave radiation (Apr; MJ m^-^2 day^-^1)",
    "rsds05" = "Surface downwelling shortwave radiation (May; MJ m^-^2 day^-^1)",
    "rsds06" = "Surface downwelling shortwave radiation (Jun; MJ m^-^2 day^-^1)",
    "rsds07" = "Surface downwelling shortwave radiation (Jul; MJ m^-^2 day^-^1)",
    "rsds08" = "Surface downwelling shortwave radiation (Aug; MJ m^-^2 day^-^1)",
    "rsds09" = "Surface downwelling shortwave radiation (Sep; MJ m^-^2 day^-^1)",
    "rsds10" = "Surface downwelling shortwave radiation (Oct; MJ m^-^2 day^-^1)",
    "rsds11" = "Surface downwelling shortwave radiation (Nov; MJ m^-^2 day^-^1)",
    "rsds12" = "Surface downwelling shortwave radiation (Dec; MJ m^-^2 day^-^1)",
    "DD_0_01" = "Degree days < 0  degreeC (Jan)",
    "DD_0_02" = "Degree days < 0  degreeC (Feb)",
    "DD_0_03" = "Degree days < 0  degreeC (Mar)",
    "DD_0_04" = "Degree days < 0  degreeC (Apr)",
    "DD_0_05" = "Degree days < 0  degreeC (May)",
    "DD_0_06" = "Degree days < 0  degreeC (Jun)",
    "DD_0_07" = "Degree days < 0  degreeC (Jul)",
    "DD_0_08" = "Degree days < 0  degreeC (Aug)",
    "DD_0_09" = "Degree days < 0  degreeC (Sep)",
    "DD_0_10" = "Degree days < 0  degreeC (Oct)",
    "DD_0_11" = "Degree days < 0  degreeC (Nov)",
    "DD_0_12" = "Degree days < 0  degreeC (Dec)",
    "DD5_01" = "Degree days > 5  degreeC (Jan)",
    "DD5_02" = "Degree days > 5  degreeC (Feb)",
    "DD5_03" = "Degree days > 5  degreeC (Mar)",
    "DD5_04" = "Degree days > 5  degreeC (Apr)",
    "DD5_05" = "Degree days > 5  degreeC (May)",
    "DD5_06" = "Degree days > 5  degreeC (Jun)",
    "DD5_07" = "Degree days > 5  degreeC (Jul)",
    "DD5_08" = "Degree days > 5  degreeC (Aug)",
    "DD5_09" = "Degree days > 5  degreeC (Sep)",
    "DD5_10" = "Degree days > 5  degreeC (Oct)",
    "DD5_11" = "Degree days > 5  degreeC (Nov)",
    "DD5_12" = "Degree days > 5  degreeC (Dec)",
    "DD_18_01" = "Degree days < 18  degreeC (Jan)",
    "DD_18_02" = "Degree days < 18  degreeC (Feb)",
    "DD_18_03" = "Degree days < 18  degreeC (Mar)",
    "DD_18_04" = "Degree days < 18  degreeC (Apr)",
    "DD_18_05" = "Degree days < 18  degreeC (May)",
    "DD_18_06" = "Degree days < 18  degreeC (Jun)",
    "DD_18_07" = "Degree days < 18  degreeC (Jul)",
    "DD_18_08" = "Degree days < 18  degreeC (Aug)",
    "DD_18_09" = "Degree days < 18  degreeC (Sep)",
    "DD_18_10" = "Degree days < 18  degreeC (Oct)",
    "DD_18_11" = "Degree days < 18  degreeC (Nov)",
    "DD_18_12" = "Degree days < 18  degreeC (Dec)",
    "DD18_01" = "Degree days > 18  degreeC (Jan)",
    "DD18_02" = "Degree days > 18  degreeC (Feb)",
    "DD18_03" = "Degree days > 18  degreeC (Mar)",
    "DD18_04" = "Degree days > 18  degreeC (Apr)",
    "DD18_05" = "Degree days > 18  degreeC (May)",
    "DD18_06" = "Degree days > 18  degreeC (Jun)",
    "DD18_07" = "Degree days > 18  degreeC (Jul)",
    "DD18_08" = "Degree days > 18  degreeC (Aug)",
    "DD18_09" = "Degree days > 18  degreeC (Sep)",
    "DD18_10" = "Degree days > 18  degreeC (Oct)",
    "DD18_11" = "Degree days > 18  degreeC (Nov)",
    "DD18_12" = "Degree days > 18  degreeC (Dec)",
    "NFFD01" = "Frost-free days (Jan; days)",
    "NFFD02" = "Frost-free days (Feb; days)",
    "NFFD03" = "Frost-free days (Mar; days)",
    "NFFD04" = "Frost-free days (Apr; days)",
    "NFFD05" = "Frost-free days (May; days)",
    "NFFD06" = "Frost-free days (Jun; days)",
    "NFFD07" = "Frost-free days (Jul; days)",
    "NFFD08" = "Frost-free days (Aug; days)",
    "NFFD09" = "Frost-free days (Sep; days)",
    "NFFD10" = "Frost-free days (Oct; days)",
    "NFFD11" = "Frost-free days (Nov; days)",
    "NFFD12" = "Frost-free days (Dec; days)",
    "Eref01" = "Reference evaporation (Jan; mm)",
    "Eref02" = "Reference evaporation (Feb; mm)",
    "Eref03" = "Reference evaporation (Mar; mm)",
    "Eref04" = "Reference evaporation (Apr; mm)",
    "Eref05" = "Reference evaporation (May; mm)",
    "Eref06" = "Reference evaporation (Jun; mm)",
    "Eref07" = "Reference evaporation (Jul; mm)",
    "Eref08" = "Reference evaporation (Aug; mm)",
    "Eref09" = "Reference evaporation (Sep; mm)",
    "Eref10" = "Reference evaporation (Oct; mm)",
    "Eref11" = "Reference evaporation (Nov; mm)",
    "Eref12" = "Reference evaporation (Dec; mm)",
    "CMD01" = "Climatic moisture deficit (Jan; mm)",
    "CMD02" = "Climatic moisture deficit (Feb; mm)",
    "CMD03" = "Climatic moisture deficit (Mar; mm)",
    "CMD04" = "Climatic moisture deficit (Apr; mm)",
    "CMD05" = "Climatic moisture deficit (May; mm)",
    "CMD06" = "Climatic moisture deficit (Jun; mm)",
    "CMD07" = "Climatic moisture deficit (Jul; mm)",
    "CMD08" = "Climatic moisture deficit (Aug; mm)",
    "CMD09" = "Climatic moisture deficit (Sep; mm)",
    "CMD10" = "Climatic moisture deficit (Oct; mm)",
    "CMD11" = "Climatic moisture deficit (Nov; mm)",
    "CMD12" = "Climatic moisture deficit (Dec; mm)",
    "RH01" = "Relative humidity (Jan; %)",
    "RH02" = "Relative humidity (Feb; %)",
    "RH03" = "Relative humidity (Mar; %)",
    "RH04" = "Relative humidity (Apr; %)",
    "RH05" = "Relative humidity (May; %)",
    "RH06" = "Relative humidity (Jun; %)",
    "RH07" = "Relative humidity (Jul; %)",
    "RH08" = "Relative humidity (Aug; %)",
    "RH09" = "Relative humidity (Sep; %)",
    "RH10" = "Relative humidity (Oct; %)",
    "RH11" = "Relative humidity (Nov; %)",
    "RH12" = "Relative humidity (Dec; %)",
    "CMI01" = "Climatic moisture index (Jan)",
    "CMI02" = "Climatic moisture index (Feb)",
    "CMI03" = "Climatic moisture index (Mar)",
    "CMI04" = "Climatic moisture index (Apr)",
    "CMI05" = "Climatic moisture index (May)",
    "CMI06" = "Climatic moisture index (Jun)",
    "CMI07" = "Climatic moisture index (Jul)",
    "CMI08" = "Climatic moisture index (Aug)",
    "CMI09" = "Climatic moisture index (Sep)",
    "CMI10" = "Climatic moisture index (Oct)",
    "CMI11" = "Climatic moisture index (Nov)",
    "CMI12" = "Climatic moisture index (Dec)",

    # ENVIREM
    "annualPET" = "Annual potential evapotranspiration (mm)",
    "aridityIndexThornthwaite" = "Thornthwaite aridity index",
    "climaticMoistureIndex" = "Climatic moisture index",
    "continentality" = "Avg temperature of warmest month minus avg temperature of coldest ( degreeC)",
    "embergerQ" = "Emberger's pluviothermic quotient",
    "growingDegDays0" = "Growing degree days > 0  degreeC",
    "growingDegDays5" = "Growing degree days > 5  degreeC",
    "maxTempColdest" = "Maximum temperature of coldest month ( degreeC)",
    "minTempWarmest" = "Minimum temperature of warmest month ( degreeC)",
    "monthCountByTemp10" = "Months with mean temperature > 10  degreeC (months)",
    "PETColdestQuarter" = "PET coldest quarter (mm)",
    "PETDriestQuarter" = "PET driest quarter (mm)",
    "PETseasonality" = "PET monthly variability",
    "PETWarmestQuarter" = "PET warmest quarter (mm)",
    "PETWettestQuarter" = "PET wettest quarter (mm)",
    "thermicityIndex" = "Compensated thermicity index ( degreeC)",

    # Urbanization (Venter et al. 2016)
    "Footprint" = "Human footprint",

    # Landcover variables (ESA WorldCover)
    "trees" = "Tree cover (%)",
    "shrubs" = "Shrub cover (%)",
    "grassland" = "Grassland cover (%)",
    "cropland" = "Cropland cover (%)",
    "built" = "Built-up land (%)",
    "bare" = "Bare & sparsely vegetated ground (%)",
    "snow" = "Snow & ice cover (%)",
    "water" = "Permanent water bodies (%)",
    "wetland" = "Herbaceous wetland cover (%)",
    "mangroves" = "Mangrove cover (%)",
    "moss" = "Moss & lichen cover (%)",

    # Soil variables (SoilGRIDS)
    "Bulk_density" = "Bulk density (kg/dm^3)",
    "Coarse_fragments_volume" = "Coarse fragments (vol. %)",
    "Clay_fraction" = "Clay fraction (%)",
    "Nitrogen_content" = "Total nitrogen content (%)",
    "Organic_carb_density" = "Organic carbon density (kg/m^2)",
    "pH_H2O" = "Soil pH (H_2O)",
    "Sand_fraction" = "Sand fraction (%)",
    "Silt_fraction" = "Silt fraction (%)",
    "Soil_organic_carb" = "Soil organic carbon (%)",

    # Forest height (ETH)
    "Forest_height" = "Forest canopy height (m)",

    # Atmosphere (Worldclim)
    "srad_01" = "Shortwave solar radiation (Jan-Feb; kJ m^-^2 day^-^1)",
    "srad_02" = "Shortwave solar radiation (Mar-Apr; kJ m^-^2 day^-^1)",
    "srad_03" = "Shortwave solar radiation (May-Jun; kJ m^-^2 day^-^1)",
    "srad_04" = "Shortwave solar radiation (Jul-Aug; kJ m^-^2 day^-^1)",
    "srad_05" = "Shortwave solar radiation (Sep-Oct; kJ m^-^2 day^-^1)",
    "srad_06" = "Shortwave solar radiation (Nov-Dec; kJ m^-^2 day^-^1)",
    "srad_median" = "Shortwave solar radiation (median; kJ m^-^2 day^-^1)",
    "srad_min" = "Shortwave solar radiation (minimum; kJ m^-^2 day^-^1)",
    "srad_max" = "Shortwave solar radiation (maximum; kJ m^-^2 day^-^1)",
    "wind_01" = "Mean wind speed (Jan-Feb; m s^-^1)",
    "wind_02" = "Mean wind speed (Mar-Apr; m s^-^1)",
    "wind_03" = "Mean wind speed (May-Jun; m s^-^1)",
    "wind_04" = "Mean wind speed (Jul-Aug; m s^-^1)",
    "wind_05" = "Mean wind speed (Sep-Oct; m s^-^1)",
    "wind_06" = "Mean wind speed (Nov-Dec; m s^-^1)",
    "wind_median" = "Mean wind speed (median; m s^-^1)",
    "wind_min" = "Mean wind speed (minimum; m s^-^1)",
    "wind_max" = "Mean wind speed (maximum; m s^-^1)",
    "vapr_01" = "Vapor pressure (Jan-Feb; kPa)",
    "vapr_02" = "Vapor pressure (Mar-Apr; kPa)",
    "vapr_03" = "Vapor pressure (May-Jun; kPa)",
    "vapr_04" = "Vapor pressure (Jul-Aug; kPa)",
    "vapr_05" = "Vapor pressure (Sep-Oct; kPa)",
    "vapr_06" = "Vapor pressure (Nov-Dec; kPa)",
    "vapr_median" = "Vapor pressure (median; kPa)",
    "vapr_min" = "Vapor pressure (minimum; kPa)",
    "vapr_max" = "Vapor pressure (maximum; kPa)",

    # Enhanced Vegetation Index (Open Land Map)
    "EVI_1" = "Enhanced Vegetation Index (Jan-Feb)",
    "EVI_2" = "Enhanced Vegetation Index (Mar-Apr)",
    "EVI_3" = "Enhanced Vegetation Index (May-Jun)",
    "EVI_4" = "Enhanced Vegetation Index (Jul-Aug)",
    "EVI_5" = "Enhanced Vegetation Index (Sep-Oct)",
    "EVI_6" = "Enhanced Vegetation Index (Nov-Dec)",
    "EVI_median" = "Enhanced Vegetation Index (median)",
    "EVI_min" = "Enhanced Vegetation Index (min)",
    "EVI_max" = "Enhanced Vegetation Index (max)",

    # Terrain metrics (Wilson et al. 2007; Beven & Kirkby 1979)
    "TRI" = "Terrain ruggedness index",
    "TPI" = "Topographic position index",
    "roughness" = "Surface roughness",
    "HLI" = "Heat load index",
    "TWI" = "Topographic wetness index",

    # Nighttime light (DMSP-OLS)
    "Nighttime_light" = "Nighttime light intensity",

    # Burned area (ESA FireCCI)
    "Burned_area_01" = "Burned area (Jan; km^2 per pixel)",
    "Burned_area_02" = "Burned area (Feb; km^2 per pixel)",
    "Burned_area_03" = "Burned area (Mar; km^2 per pixel)",
    "Burned_area_04" = "Burned area (Apr; km^2 per pixel)",
    "Burned_area_05" = "Burned area (May; km^2 per pixel)",
    "Burned_area_06" = "Burned area (Jun; km^2 per pixel)",
    "Burned_area_07" = "Burned area (Jul; km^2 per pixel)",
    "Burned_area_08" = "Burned area (Aug; km^2 per pixel)",
    "Burned_area_09" = "Burned area (Sep; km^2 per pixel)",
    "Burned_area_10" = "Burned area (Oct; km^2 per pixel)",
    "Burned_area_11" = "Burned area (Nov; km^2 per pixel)",
    "Burned_area_12" = "Burned area (Dec; km^2 per pixel)",

    # Snow water equivalent (Daymet v4)
    "snow_water_equivalent_01" = "Snow water equivalent (Jan; mm)",
    "snow_water_equivalent_02" = "Snow water equivalent (Feb; mm)",
    "snow_water_equivalent_03" = "Snow water equivalent (Mar; mm)",
    "snow_water_equivalent_04" = "Snow water equivalent (Apr; mm)",
    "snow_water_equivalent_05" = "Snow water equivalent (May; mm)",
    "snow_water_equivalent_06" = "Snow water equivalent (Jun; mm)",
    "snow_water_equivalent_07" = "Snow water equivalent (Jul; mm)",
    "snow_water_equivalent_08" = "Snow water equivalent (Aug; mm)",
    "snow_water_equivalent_09" = "Snow water equivalent (Sep; mm)",
    "snow_water_equivalent_10" = "Snow water equivalent (Oct; mm)",
    "snow_water_equivalent_11" = "Snow water equivalent (Nov; mm)",
    "snow_water_equivalent_12" = "Snow water equivalent (Dec; mm)",

    # Daylength (Daymet v4)
    "Daylength_01" = "Daylength (Jan; s day^-^1)",
    "Daylength_02" = "Daylength (Feb; s day^-^1)",
    "Daylength_03" = "Daylength (Mar; s day^-^1)",
    "Daylength_04" = "Daylength (Apr; s day^-^1)",
    "Daylength_05" = "Daylength (May; s day^-^1)",
    "Daylength_06" = "Daylength (Jun; s day^-^1)",
    "Daylength_07" = "Daylength (Jul; s day^-^1)",
    "Daylength_08" = "Daylength (Aug; s day^-^1)",
    "Daylength_09" = "Daylength (Sep; s day^-^1)",
    "Daylength_10" = "Daylength (Oct; s day^-^1)",
    "Daylength_11" = "Daylength (Nov; s day^-^1)",
    "Daylength_12" = "Daylength (Dec; s day^-^1)",

    # Soil moisture (ESA CCI)
    "soil_moisture_01" = "Soil moisture (Jan; m^3 m^-^3)",
    "soil_moisture_02" = "Soil moisture (Feb; m^3 m^-^3)",
    "soil_moisture_03" = "Soil moisture (Mar; m^3 m^-^3)",
    "soil_moisture_04" = "Soil moisture (Apr; m^3 m^-^3)",
    "soil_moisture_05" = "Soil moisture (May; m^3 m^-^3)",
    "soil_moisture_06" = "Soil moisture (Jun; m^3 m^-^3)",

    # Bird richness (IUCN Red List)
    "bird_species_richness" = "Bird species richness"
  )

  # Short mapping
  variable.mapping.short <- c(

    # ClimateNA
    "MAT" = "Avg annual temp",
    "MWMT" = "Avg warmest month temp",
    "MCMT" = "Avg coldest month temp",
    "TD" = "Temp difference",
    "MAP" = "Avg annual precipitation",
    "MSP" = "May to Sep precipitation",
    "AHM" = "Annual heat-moisture",
    "SHM" = "Summer heat-moisture",
    "bFFP" = "Beginning of frost-free period",
    "eFFP" = "End of frost-free period",
    "FFP" = "Annual frost-free period days",
    "CMD" = "Annual climatic moisture deficit",
    "CMI" = "Annual climatic moisture",
    "DD_0" = "Annual degree days below 0  degreeC",
    "DD5" = "Annual degree days above 5  degreeC",
    "DD_18" = "Annual degree days below 18  degreeC",
    "DD18" = "Annual degree days above 18  degreeC",
    "DD1040" = "Annual degree days 10-40  degreeC",
    "EMT" = "Extreme min temp",
    "EXT" = "Extreme max temp",
    "Eref" = "Annual reference evaporation",
    "rsds" = "Avg annual solar radiation",
    "NFFD" = "Annual frost-free days",
    "PAS" = "Annual precipitation as snow",
    "RH"  = "Avg annual relative humidity",
    "Tmax01" = "Max temp Jan",
    "Tmax02" = "Max temp Feb",
    "Tmax03" = "Max temp Mar",
    "Tmax04" = "Max temp Apr",
    "Tmax05" = "Max temp May",
    "Tmax06" = "Max temp Jun",
    "Tmax07" = "Max temp Jul",
    "Tmax08" = "Max temp Aug",
    "Tmax09" = "Max temp Sep",
    "Tmax10" = "Max temp Oct",
    "Tmax11" = "Max temp Nov",
    "Tmax12" = "Max temp Dec",
    "Tmin01" = "Min temp Jan",
    "Tmin02" = "Min temp Feb",
    "Tmin03" = "Min temp Mar",
    "Tmin04" = "Min temp Apr",
    "Tmin05" = "Min temp May",
    "Tmin06" = "Min temp Jun",
    "Tmin07" = "Min temp Jul",
    "Tmin08" = "Min temp Aug",
    "Tmin09" = "Min temp Sep",
    "Tmin10" = "Min temp Oct",
    "Tmin11" = "Min temp Nov",
    "Tmin12" = "Min temp Dec",
    "Tave01" = "Avg temp Jan",
    "Tave02" = "Avg temp Feb",
    "Tave03" = "Avg temp Mar",
    "Tave04" = "Avg temp Apr",
    "Tave05" = "Avg temp May",
    "Tave06" = "Avg temp Jun",
    "Tave07" = "Avg temp Jul",
    "Tave08" = "Avg temp Aug",
    "Tave09" = "Avg temp Sep",
    "Tave10" = "Avg temp Oct",
    "Tave11" = "Avg temp Nov",
    "Tave12" = "Avg temp Dec",
    "PAS01" = "Snow precipitation Jan",
    "PAS02" = "Snow precipitation Feb",
    "PAS03" = "Snow precipitation Mar",
    "PAS04" = "Snow precipitation Apr",
    "PAS05" = "Snow precipitation May",
    "PAS06" = "Snow precipitation Jun",
    "PAS07" = "Snow precipitation Jul",
    "PAS08" = "Snow precipitation Aug",
    "PAS09" = "Snow precipitation Sep",
    "PAS10" = "Snow precipitation Oct",
    "PAS11" = "Snow precipitation Nov",
    "PAS12" = "Snow precipitation Dec",
    "PPT01" = "Precipitation Jan",
    "PPT02" = "Precipitation Feb",
    "PPT03" = "Precipitation Mar",
    "PPT04" = "Precipitation Apr",
    "PPT05" = "Precipitation May",
    "PPT06" = "Precipitation Jun",
    "PPT07" = "Precipitation Jul",
    "PPT08" = "Precipitation Aug",
    "PPT09" = "Precipitation Sep",
    "PPT10" = "Precipitation Oct",
    "PPT11" = "Precipitation Nov",
    "PPT12" = "Precipitation Dec",
    "rsds01" = "Shortwave radiation Jan",
    "rsds02" = "Shortwave radiation Feb",
    "rsds03" = "Shortwave radiation Mar",
    "rsds04" = "Shortwave radiation Apr",
    "rsds05" = "Shortwave radiation May",
    "rsds06" = "Shortwave radiation Jun",
    "rsds07" = "Shortwave radiation Jul",
    "rsds08" = "Shortwave radiation Aug",
    "rsds09" = "Shortwave radiation Sep",
    "rsds10" = "Shortwave radiation Oct",
    "rsds11" = "Shortwave radiation Nov",
    "rsds12" = "Shortwave radiation Dec",
    "DD_0_01" = "Degree days < 0  degreeC Jan",
    "DD_0_02" = "Degree days < 0  degreeC Feb",
    "DD_0_03" = "Degree days < 0  degreeC Mar",
    "DD_0_04" = "Degree days < 0  degreeC Apr",
    "DD_0_05" = "Degree days < 0  degreeC May",
    "DD_0_06" = "Degree days < 0  degreeC Jun",
    "DD_0_07" = "Degree days < 0  degreeC Jul",
    "DD_0_08" = "Degree days < 0  degreeC Aug",
    "DD_0_09" = "Degree days < 0  degreeC Sep",
    "DD_0_10" = "Degree days < 0  degreeC Oct",
    "DD_0_11" = "Degree days < 0  degreeC Nov",
    "DD_0_12" = "Degree days < 0  degreeC Dec",
    "DD5_01" = "Degree days > 5  degreeC Jan",
    "DD5_02" = "Degree days > 5  degreeC Feb",
    "DD5_03" = "Degree days > 5  degreeC Mar",
    "DD5_04" = "Degree days > 5  degreeC Apr",
    "DD5_05" = "Degree days > 5  degreeC May",
    "DD5_06" = "Degree days > 5  degreeC Jun",
    "DD5_07" = "Degree days > 5  degreeC Jul",
    "DD5_08" = "Degree days > 5  degreeC Aug",
    "DD5_09" = "Degree days > 5  degreeC Sep",
    "DD5_10" = "Degree days > 5  degreeC Oct",
    "DD5_11" = "Degree days > 5  degreeC Nov",
    "DD5_12" = "Degree days > 5  degreeC Dec",
    "DD_18_01" = "Degree days < 18  degreeC Jan",
    "DD_18_02" = "Degree days < 18  degreeC Feb",
    "DD_18_03" = "Degree days < 18  degreeC Mar",
    "DD_18_04" = "Degree days < 18  degreeC Apr",
    "DD_18_05" = "Degree days < 18  degreeC May",
    "DD_18_06" = "Degree days < 18  degreeC Jun",
    "DD_18_07" = "Degree days < 18  degreeC Jul",
    "DD_18_08" = "Degree days < 18  degreeC Aug",
    "DD_18_09" = "Degree days < 18  degreeC Sep",
    "DD_18_10" = "Degree days < 18  degreeC Oct",
    "DD_18_11" = "Degree days < 18  degreeC Nov",
    "DD_18_12" = "Degree days < 18  degreeC Dec",
    "DD18_01" = "Degree days > 18  degreeC Jan",
    "DD18_02" = "Degree days > 18  degreeC Feb",
    "DD18_03" = "Degree days > 18  degreeC Mar",
    "DD18_04" = "Degree days > 18  degreeC Apr",
    "DD18_05" = "Degree days > 18  degreeC May",
    "DD18_06" = "Degree days > 18  degreeC Jun",
    "DD18_07" = "Degree days > 18  degreeC Jul",
    "DD18_08" = "Degree days > 18  degreeC Aug",
    "DD18_09" = "Degree days > 18  degreeC Sep",
    "DD18_10" = "Degree days > 18  degreeC Oct",
    "DD18_11" = "Degree days > 18  degreeC Nov",
    "DD18_12" = "Degree days > 18  degreeC Dec",
    "NFFD01" = "Frost-free days Jan",
    "NFFD02" = "Frost-free days Feb",
    "NFFD03" = "Frost-free days Mar",
    "NFFD04" = "Frost-free days Apr",
    "NFFD05" = "Frost-free days May",
    "NFFD06" = "Frost-free days Jun",
    "NFFD07" = "Frost-free days Jul",
    "NFFD08" = "Frost-free days Aug",
    "NFFD09" = "Frost-free days Sep",
    "NFFD10" = "Frost-free days Oct",
    "NFFD11" = "Frost-free days Nov",
    "NFFD12" = "Frost-free days Dec",
    "Eref01" = "Reference evaporation Jan",
    "Eref02" = "Reference evaporation Feb",
    "Eref03" = "Reference evaporation Mar",
    "Eref04" = "Reference evaporation Apr",
    "Eref05" = "Reference evaporation May",
    "Eref06" = "Reference evaporation Jun",
    "Eref07" = "Reference evaporation Jul",
    "Eref08" = "Reference evaporation Aug",
    "Eref09" = "Reference evaporation Sep",
    "Eref10" = "Reference evaporation Oct",
    "Eref11" = "Reference evaporation Nov",
    "Eref12" = "Reference evaporation Dec",
    "CMD01" = "Climatic moisture deficit Jan",
    "CMD02" = "Climatic moisture deficit Feb",
    "CMD03" = "Climatic moisture deficit Mar",
    "CMD04" = "Climatic moisture deficit Apr",
    "CMD05" = "Climatic moisture deficit May",
    "CMD06" = "Climatic moisture deficit Jun",
    "CMD07" = "Climatic moisture deficit Jul",
    "CMD08" = "Climatic moisture deficit Aug",
    "CMD09" = "Climatic moisture deficit Sep",
    "CMD10" = "Climatic moisture deficit Oct",
    "CMD11" = "Climatic moisture deficit Nov",
    "CMD12" = "Climatic moisture deficit Dec",
    "RH01" = "Relative humidity Jan",
    "RH02" = "Relative humidity Feb",
    "RH03" = "Relative humidity Mar",
    "RH04" = "Relative humidity Apr",
    "RH05" = "Relative humidity May",
    "RH06" = "Relative humidity Jun",
    "RH07" = "Relative humidity Jul",
    "RH08" = "Relative humidity Aug",
    "RH09" = "Relative humidity Sep",
    "RH10" = "Relative humidity Oct",
    "RH11" = "Relative humidity Nov",
    "RH12" = "Relative humidity Dec",
    "CMI01" = "Climatic moisture Jan",
    "CMI02" = "Climatic moisture Feb",
    "CMI03" = "Climatic moisture Mar",
    "CMI04" = "Climatic moisture Apr",
    "CMI05" = "Climatic moisture May",
    "CMI06" = "Climatic moisture Jun",
    "CMI07" = "Climatic moisture Jul",
    "CMI08" = "Climatic moisture Aug",
    "CMI09" = "Climatic moisture Sep",
    "CMI10" = "Climatic moisture Oct",
    "CMI11" = "Climatic moisture Nov",
    "CMI12" = "Climatic moisture Dec",

    # ENVIREM
    "annualPET" = "Annual PET",
    "aridityIndexThornthwaite" = "Annual aridity",
    "climaticMoistureIndex" = "Climatic moisture index",
    "continentality" = "Temp range warmest-coldest month",
    "embergerQ" = "Emberger's pluviothermic quotient",
    "growingDegDays0" = "Growing degree days > 0  degreeC",
    "growingDegDays5" = "Growing degree days > 5  degreeC",
    "maxTempColdest" = "Max temp of coldest month",
    "minTempWarmest" = "Min temp of warmest month",
    "monthCountByTemp10" = "Months with avg temp > 10  degreeC",
    "PETColdestQuarter" = "PET coldest quarter",
    "PETDriestQuarter" = "PET driest quarter",
    "PETseasonality" = "PET monthly variability",
    "PETWarmestQuarter" = "PET warmest quarter",
    "PETWettestQuarter" = "PET wettest quarter",
    "thermicityIndex" = "Annual thermicity",

    # Urbanization (Venter et al. 2016)
    "Footprint" = "Human footprint",

    # Landcover variables (ESA WorldCover)
    "trees" = "Tree cover",
    "shrubs" = "Shrub cover",
    "grassland" = "Grassland cover",
    "cropland" = "Cropland cover",
    "built" = "Built-up land",
    "bare" = "Bare ground",
    "snow" = "Snow cover",
    "water" = "Water cover",
    "wetland" = "Wetland cover",
    "mangroves" = "Mangrove cover",
    "moss" = "Moss cover",

    # Soil variables (SoilGRIDS)
    "Bulk_density" = "Soil bulk density",
    "Coarse_fragments_volume" = "Soil coarse fragments",
    "Clay_fraction" = "Soil clay fraction",
    "Nitrogen_content" = "Soil nitrogen content",
    "Organic_carb_density" = "Soil organic carbon density",
    "pH_H2O" = "Soil pH H_2O",
    "Sand_fraction" = "Soil sand fraction",
    "Silt_fraction" = "Soil silt fraction",
    "Soil_organic_carb" = "Soil organic carbon",

    # Forest height (ETH)
    "Forest_height" = "Forest height",

    # Atmosphere (Worldclim)
    "srad_01" = "Solar radiation Jan-Feb",
    "srad_02" = "Solar radiation Mar-Apr",
    "srad_03" = "Solar radiation May-Jun",
    "srad_04" = "Solar radiation Jul-Aug",
    "srad_05" = "Solar radiation Sep-Oct",
    "srad_06" = "Solar radiation Nov-Dec",
    "srad_median" = "Solar radiation median",
    "srad_min" = "Solar radiation min",
    "srad_max" = "Solar radiation max",
    "wind_01" = "Avg wind speed Jan-Feb",
    "wind_02" = "Avg wind speed Mar-Apr",
    "wind_03" = "Avg wind speed May-Jun",
    "wind_04" = "Avg wind speed Jul-Aug",
    "wind_05" = "Avg wind speed Sep-Oct",
    "wind_06" = "Avg wind speed Nov-Dec",
    "wind_median" = "Avg wind speed median",
    "wind_min" = "Avg wind speed min",
    "wind_max" = "Avg wind speed max",
    "vapr_01" = "Vapor pressure Jan-Feb",
    "vapr_02" = "Vapor pressure Mar-Apr",
    "vapr_03" = "Vapor pressure May-Jun",
    "vapr_04" = "Vapor pressure Jul-Aug",
    "vapr_05" = "Vapor pressure Sep-Oct",
    "vapr_06" = "Vapor pressure Nov-Dec",
    "vapr_median" = "Vapor pressure median",
    "vapr_min" = "Vapor pressure min",
    "vapr_max" = "Vapor pressure max",

    # Enhanced Vegetation Index (Open Land Map)
    "EVI_1" = "Enhanced vegetation Jan-Feb",
    "EVI_2" = "Enhanced vegetation Mar-Apr",
    "EVI_3" = "Enhanced vegetation May-Jun",
    "EVI_4" = "Enhanced vegetation Jul-Aug",
    "EVI_5" = "Enhanced vegetation Sep-Oct",
    "EVI_6" = "Enhanced vegetation Nov-Dec",
    "EVI_median" = "Enhanced vegetation median",
    "EVI_min" = "Enhanced vegetation min",
    "EVI_max" = "Enhanced vegetation max",

    # Terrain metrics (Wilson et al. 2007; Beven & Kirkby 1979)
    "TRI" = "Terrain ruggedness",
    "TPI" = "Topographic position",
    "roughness" = "Surface roughness",
    "HLI" = "Heat load",
    "TWI" = "Topographic wetness",

    # Nighttime light (DMSP-OLS)
    "Nighttime_light" = "Night light",

    # Burned area (ESA FireCCI)
    "Burned_area_01" = "Burned area Jan",
    "Burned_area_02" = "Burned area Feb",
    "Burned_area_03" = "Burned area Mar",
    "Burned_area_04" = "Burned area Apr",
    "Burned_area_05" = "Burned area May",
    "Burned_area_06" = "Burned area Jun",
    "Burned_area_07" = "Burned area Jul",
    "Burned_area_08" = "Burned area Aug",
    "Burned_area_09" = "Burned area Sep",
    "Burned_area_10" = "Burned area Oct",
    "Burned_area_11" = "Burned area Nov",
    "Burned_area_12" = "Burned area Dec",

    # Snow water equivalent (Daymet v4)
    "snow_water_equivalent_01" = "Snow water Jan",
    "snow_water_equivalent_02" = "Snow water Feb",
    "snow_water_equivalent_03" = "Snow water Mar",
    "snow_water_equivalent_04" = "Snow water Apr",
    "snow_water_equivalent_05" = "Snow water May",
    "snow_water_equivalent_06" = "Snow water Jun",
    "snow_water_equivalent_07" = "Snow water Jul",
    "snow_water_equivalent_08" = "Snow water Aug",
    "snow_water_equivalent_09" = "Snow water Sep",
    "snow_water_equivalent_10" = "Snow water Oct",
    "snow_water_equivalent_11" = "Snow water Nov",
    "snow_water_equivalent_12" = "Snow water Dec",

    # Daylength (Daymet v4)
    "Daylength_01" = "Daylength Jan",
    "Daylength_02" = "Daylength Feb",
    "Daylength_03" = "Daylength Mar",
    "Daylength_04" = "Daylength Apr",
    "Daylength_05" = "Daylength May",
    "Daylength_06" = "Daylength Jun",
    "Daylength_07" = "Daylength Jul",
    "Daylength_08" = "Daylength Aug",
    "Daylength_09" = "Daylength Sep",
    "Daylength_10" = "Daylength Oct",
    "Daylength_11" = "Daylength Nov",
    "Daylength_12" = "Daylength Dec",

    # Soil moisture (ESA CCI)
    "soil_moisture_01" = "Soil moisture Jan",
    "soil_moisture_02" = "Soil moisture Feb",
    "soil_moisture_03" = "Soil moisture Mar",
    "soil_moisture_04" = "Soil moisture Apr",
    "soil_moisture_05" = "Soil moisture May",
    "soil_moisture_06" = "Soil moisture Jun",

    # Bird richness (IUCN Red List)
    "bird_species_richness" = "Bird richness"
  )

  # Select which mapping version to use
  name.length <- match.arg(name.length)
  variable.mapping <- if (name.length == "full") variable.mapping.full else variable.mapping.short

  # Handle transformed variables
  if (recognize.transformations) {
    transformation.suffixes <- c("_log",
                                 "_log1p",
                                 "_log1p_shifted",
                                 "_sqrt", "_sqrt_shifted",
                                 "_cuberoot",
                                 "_logit_shrunk",
                                 "_arcsine_sqrt")
  }

  # Create function to map single variable name
  map.single.variable <- function(variable.name) {
    base.variable.name <- variable.name
    if (recognize.transformations) {
      base.variable.name <- sub(paste0("(", paste0(transformation.suffixes, collapse = "|"), ")$"), "", variable.name, perl = TRUE)
    }
    if (base.variable.name %in% names(variable.mapping)) {
      mapped.name <- variable.mapping[[base.variable.name]]
      if (recognize.transformations && !identical(base.variable.name, variable.name)) {
        suffix.raw <- sub(base.variable.name, "", variable.name, fixed = TRUE)
        suffix.clean <- gsub("^_", "", suffix.raw)
        if (grepl("\\)$", mapped.name)) {
          mapped.name <- sub("\\)$", paste0("; ", suffix.clean, ")"), mapped.name)
        } else {
          mapped.name <- paste0(mapped.name, " (", suffix.clean, ")")
        }
      }
      mapped.name
    } else {
      variable.name
    }
  }

  # Apply mapping to data frame columns or vector names
  if (is.data.frame(input.data)) {
    colnames(input.data) <- vapply(colnames(input.data), map.single.variable, FUN.VALUE = character(1))
  } else if (is.vector(input.data) && !is.null(names(input.data))) {
    names(input.data) <- vapply(names(input.data), map.single.variable, FUN.VALUE = character(1))
  }

  # Return dataframe or vector
  return(input.data)
}


## Function to extract environmental variables and generate background data
#' Extract environmental variables and optional background data
#'
#' Extract environmental variables for occurrence records, optionally generate
#' random background points within an accessible area, and write occurrence and
#' background environmental tables to disk.
#'
#' @param occurrence.data A `data.frame` containing occurrence records with
#'   unique row names used as record identifiers, and coordinate columns named by
#'   `longitude.col` and `latitude.col`.
#' @param longitude.col A single character string giving the longitude column
#'   name in `occurrence.data` (default: `"Longitude"`).
#' @param latitude.col A single character string giving the latitude column name
#'   in `occurrence.data` (default: `"Latitude"`).
#' @param generate.background.data Logical; if `TRUE`, random background points
#'   are generated within the accessible area (default: `FALSE`).
#' @param N.background.points A single positive numeric value giving the number
#'   of background points to generate when `generate.background.data = TRUE`
#'   (default: `100000`).
#' @param remove.hydrolakes.background Logical; if `TRUE`, HydroLAKES water
#'   bodies are excluded from background sampling (default: `FALSE`).
#' @param buffer.km A single positive numeric value giving the buffer distance,
#'   in kilometers, used to expand the occurrence-based accessible area.
#' @param landmask.largest.N.pieces A single positive integer-like numeric value
#'   giving the number of largest landmask polygon pieces to retain
#'   (recommended default: `5`).
#' @param csv.occurrence.out.file A single character string giving the output CSV
#'   file name for the occurrence environmental table.
#' @param csv.background.out.file Optional character string giving the output CSV
#'   file name for the background environmental table. Required when
#'   `generate.background.data = TRUE` (default: `NULL`).
#' @param output.dir A single character string giving the output directory for
#'   final CSV files and intermediate folders.
#' @param intermediate.files.dir A single character string giving the subfolder
#'   used for cached intermediate files (default: `"Intermediate_files"`).
#' @param env.datasets Optional character vector naming the built-in
#'   environmental datasets to extract. Valid entries are `"elevation"`,
#'   `"ClimateNA"`, `"EVI"`, `"terrain"`, `"ENVIREM"`, `"footprint"`,
#'   `"landcover"`, `"soil"`, `"forest_height"`, `"atmosphere"`,
#'   `"nightlight"`, `"burned_area"`, `"snow_water_equivalent"`,
#'   `"daylength"`, and `"soil_moisture"` (default: `NULL`). `"ClimateNA"` and
#'   `"terrain"` require `"elevation"`. `"ClimateNA"`, `"daylength"`, and
#'   `"snow_water_equivalent"` are only available for North America.
#' @param CRS.occurrences A single character string giving the coordinate
#'   reference system of the occurrence coordinates (default: `"EPSG:4326"`).
#' @param overwrite Logical; if `TRUE`, existing outputs and cached intermediate
#'   files are overwritten where applicable (default: `FALSE`).
#' @param delete.intermediate.files.folders Logical; if `TRUE`, intermediate
#'   raster and auxiliary folders are deleted after processing
#'   (default: `FALSE`).
#' @param redownload.rasters Logical; if `TRUE`, raster and auxiliary input files
#'   are downloaded again even when cached files are present (default: `FALSE`).
#' @param custom.env.rasters Optional character vector of custom raster sources.
#'   Each entry can be a local `.tif` or `.tiff`, a directory containing GeoTIFFs,
#'   a `.zip` archive containing GeoTIFFs, or an HTTP/HTTPS URL to a `.tif`,
#'   `.tiff`, or `.zip` file (default: `NULL`).
#' @param custom.env.rasters.names Optional character vector naming each custom
#'   raster dataset. If supplied, it must have the same length as
#'   `custom.env.rasters` (default: `NULL`).
#' @param custom.env.rasters.variable.names Optional character vector or list of
#'   character vectors giving variable names assigned to custom raster layers
#'   (default: `NULL`).
#' @param custom.env.rasters.crs Optional character vector of CRS definitions for
#'   custom rasters. This is used only when a custom raster lacks embedded CRS
#'   metadata (default: `NULL`).
#' @param custom.raster.min.size.mb A single positive numeric value giving the
#'   minimum file size, in megabytes, required for downloaded custom raster files
#'   to be considered valid (default: `5`).
#' @param seed A single numeric value used to set the random seed for
#'   reproducible background sampling (default: `1`).
#' @param verbose Logical; if `TRUE`, progress messages are printed
#'   (default: `TRUE`).
#'
#' @details
#' This function is designed to create occurrence and optional background
#' environmental tables for downstream DAPC niche-divergence analyses.
#'
#' Background points represent the environmental conditions available within an
#' estimated accessible area. Background sampling can optionally exclude major
#' lakes using HydroLAKES. This is useful for terrestrial taxa because aquatic
#' surfaces may represent inaccessible or biologically irrelevant environments,
#' and including them can distort the environmental background used in downstream
#' comparisons. For aquatic or semi-aquatic organisms, or for very large
#' geographic extents where lake masking is computationally costly, retaining
#' water bodies may be more appropriate.
#'
#' The function supports a broad suite of abiotic and biotic predictors because
#' ecological niches are multidimensional and often shaped by climate,
#' topography, hydrology, vegetation, soils, disturbance, land use, and biotic
#' context rather than climate alone (Guisan & Zimmermann, 2000; Soberón, 2007;
#' Elith & Leathwick, 2009; Kearney & Porter, 2009; Dormann et al., 2013; Title
#' & Bemmels, 2018). Seasonal and monthly predictors are included because annual
#' summaries can obscure phenological constraints, short-term physiological
#' stress, and seasonal resource availability (Hufkens et al., 2018; Zimmermann
#' et al., 2009; Prajzlerová et al., 2025).
#'
#' Elevation (`env.datasets = "elevation"`; 1 variable) is based on the
#' Copernicus Global Digital Elevation Model GLO-90 v.1.1
#' (https://doi.org/10.5069/G9028PQB), aggregated to approximately 250 m
#' resolution. Elevation can be useful for deriving climatic and topographic
#' predictors, but should usually be interpreted cautiously as a direct
#' ecological predictor because it often acts as a proxy for temperature,
#' moisture, oxygen availability, vegetation, and other mechanistic drivers
#' (Guisan & Zimmermann, 2000; Dormann et al., 2013; Title & Bemmels, 2018).
#'
#' Terrain variables (`env.datasets = "terrain"`; requires
#' `env.datasets = "elevation"`; 5 variables; approximately 250 m resolution) are
#' calculated to summarize local landscape structure, topographic exposure, and
#' hydrological position. The terrain ruggedness index (TRI) measures the mean
#' absolute difference between the elevation of a focal cell and its surrounding
#' cells, capturing local terrain heterogeneity. The topographic position index
#' (TPI) measures whether a focal cell is higher or lower than its local
#' neighborhood, distinguishing ridges, slopes, flats, and depressions. Surface
#' roughness measures the local elevation range and captures abrupt changes in
#' relief. Together, TRI, TPI, and roughness describe complementary aspects of
#' fine-scale topographic complexity and habitat heterogeneity (Wilson et al.,
#' 2007; Hengl & Evans, 2009; Conrad et al., 2015). The heat load index (HLI)
#' summarizes potential thermal exposure as a function of slope, aspect, and
#' latitude, providing a topographic proxy for solar radiation, heat balance, and
#' microclimatic exposure (McCune & Keon, 2002; Bennie et al., 2008; Dobrowski,
#' 2011). The topographic wetness index (TWI) represents soil-moisture potential
#' based on upslope contributing area and slope, linking topography to drainage,
#' hydrological accumulation, and local moisture availability (Beven & Kirkby,
#' 1979; Sørensen et al., 2006). Together, these five terrain variables capture
#' microclimatic buffering, exposure, drainage, and habitat heterogeneity relevant
#' to species distributions (Bennie et al., 2008; Dobrowski, 2011).
#'
#' ClimateNA (Wang et al., 2016) (`env.datasets = "ClimateNA"`; requires
#' `env.datasets = "elevation"`; scale-free; 205 variables) provides 25 annual
#' and 180 monthly North American climate estimates adjusted for local topography
#' and elevation. These variables represent temperature, precipitation, frost,
#' snow, humidity, evaporative demand, degree days, and heat-moisture indices,
#' capturing both long-term climatic gradients and seasonal climatic constraints
#' based on multi-year mean values for 2011-2020. ClimateNA can provide finer
#' ecological resolution than coarser gridded climatic products in mountainous or
#' otherwise environmentally heterogeneous landscapes because it combines
#' interpolation with local elevation adjustment (Wang et al., 2016; MacDonald et
#' al., 2025).
#'
#' ENVIREM (Environmental Rasters for Ecological Modeling; Title & Bemmels, 2018)
#' (`env.datasets = "ENVIREM"`; 16 variables; 30-arcsecond resolution) variables
#' provide ecologically derived climatic and moisture indices, including
#' potential evapotranspiration, aridity, climatic moisture, continentality,
#' growing degree days, thermicity, and seasonal PET summaries. These variables
#' are designed to represent ecophysiologically relevant climate gradients and
#' can improve ecological niche modeling relative to relying only on standard
#' bioclimatic variables (Title & Bemmels, 2018).
#'
#' Atmospheric variables (`env.datasets = "atmosphere"`; 3 variables with 9
#' estimates each; 30-arcsecond resolution) are based on WorldClim v.2.1 climate
#' layers (Hijmans et al., 2005; Fick & Hijmans, 2017) and include solar
#' radiation (srad), wind speed (wind), and vapor pressure (vapr). For each
#' variable, six bimonthly layers are extracted, and the median, minimum, and
#' maximum values across bimonthly layers are calculated, resulting in 27
#' atmospheric variables. Solar radiation captures atmospheric energy input, wind
#' speed reflects mechanical mixing and desiccation potential, and vapor pressure
#' represents atmospheric moisture conditions influencing evapotranspiration and
#' physiological water stress (Allen et al., 1998; Novick et al., 2016).
#'
#' Daylength (`env.datasets = "daylength"`; available only for North America; 1
#' variable with 12 monthly estimates; 30-arcsecond resolution) is based on
#' Daymet v.4 monthly surfaces (daymet.ornl.gov) and captures 12 monthly
#' photoperiod conditions. Photoperiod is a fundamental seasonal cue affecting
#' phenology, diapause, development, reproduction, and other circannual
#' biological processes, making it especially relevant for taxa with strong
#' seasonal life cycles (Hufkens et al., 2018).
#'
#' Snow water equivalent (`env.datasets = "snow_water_equivalent"`; available
#' only for North America; 1 variable with 12 monthly estimates; 30-arcsecond
#' resolution) is based on Daymet v.4 monthly surfaces (daymet.ornl.gov)
#' and represents monthly snowpack water content. This layer captures winter snow
#' accumulation and cryospheric constraints that can affect overwintering
#' conditions, surface moisture, phenology, and seasonal habitat accessibility.
#'
#' The Enhanced Vegetation Index (EVI; `env.datasets = "EVI"`; 9 variables;
#' 250 m resolution) captures vegetation productivity, canopy greenness, and
#' phenological dynamics while reducing saturation in high-biomass areas relative
#' to NDVI, making it useful for species whose distributions are linked to plant
#' productivity, host-plant availability, or seasonal vegetation structure (Huete
#' et al., 2002). EVI variables are based on OpenLandMap bimonthly vegetation
#' composites derived from MODIS MOD13Q1 products and provide six bimonthly
#' vegetation-greenness layers (Jan-Feb, Mar-Apr, May-Jun, Jul-Aug, Sep-Oct,
#' and Nov-Dec), plus median, minimum, and maximum annual summaries.
#'
#' Land-cover variables (`env.datasets = "landcover"`; 11 variables; 10 m source
#' resolution) are based on the European Space Agency WorldCover 2020 product
#' (esa-worldcover.org) and summarize proportional coverage of 11 major
#' land-cover classes: tree cover, grassland, shrubland, cropland, built-up land,
#' bare ground, snow and ice, water, wetlands, mangroves, and moss or lichen.
#' These variables capture habitat composition and broad structural differences
#' in land use and vegetation cover, which can affect resource availability,
#' dispersal, exposure, and habitat suitability.
#'
#' Forest height (`env.datasets = "forest_height"`; 1 variable; approximately
#' 250 m resolution for regional and continental products, 500 m for the global
#' product) is based on the ETH Global Canopy Height 2020 product (Lang et al.,
#' 2023) and provides one canopy-height variable representing forest structural
#' complexity. The data are derived from the original product, namely a global
#' 10 m canopy top height model estimated from Sentinel-2 imagery with NASA
#' Global Ecosystem Dynamics Investigation waveform Light Detection and Ranging
#' supervision (Lang et al., 2023). Canopy height can indicate vertical habitat
#' complexity, shading, forest maturity, and vegetation structure, all of which
#' may affect microclimate, host resources, and habitat availability (Potapov et
#' al., 2021; Lang et al., 2023).
#'
#' Soil variables (`env.datasets = "soil"`; 9 variables; 30-arcsecond
#' resolution) are based on SoilGrids 2.0 (Poggio et al., 2021; Turek et al.,
#' 2023) and provide nine topsoil properties: bulk density (`Bulk_density`),
#' coarse fragments volume (`Coarse_fragments_volume`), clay fraction
#' (`Clay_fraction`), sand fraction (`Sand_fraction`), silt fraction
#' (`Silt_fraction`), total nitrogen content (`Nitrogen_content`), soil pH in
#' water (`pH_H2O`), soil organic carbon (`Soil_organic_carb`), and organic carbon
#' density (`Organic_carb_density`). These variables represent substrate texture,
#' chemistry, nutrient availability, and water-holding capacity, which can
#' influence vegetation composition, host-plant distributions, and habitat
#' suitability (Hillel, 1998; Jackson et al., 1996; Weil & Brady, 2016).
#'
#' Soil moisture (`env.datasets = "soil_moisture"`; 1 variable with 6 bimonthly
#' estimates; 900-arcsecond resolution) is based on the ESA CCI Soil Moisture
#' product and provides six bimonthly volumetric soil-water availability layers
#' (Dorigo et al., 2017; Gruber et al., 2019; Preimesberger et al., 2021). Soil
#' moisture integrates climatic water balance, surface hydrology, and
#' vegetation-water interactions, and can be important for species limited by
#' drought stress, larval host-plant condition, or moisture-dependent habitat
#' structure.
#'
#' Human footprint (`env.datasets = "footprint"`; 1 variable; 30-arcsecond
#' resolution) is based on the Human Footprint Index and provides one variable
#' (`Footprint`) representing cumulative human modification of terrestrial
#' environments (Venter et al., 2016). It integrates built environments,
#' population density, electrical infrastructure, croplands, pasturelands, roads,
#' railways, and navigable waterways, providing a broad proxy for anthropogenic
#' disturbance and habitat alteration.
#'
#' Nighttime light (`env.datasets = "nightlight"`; 1 variable; 30-arcsecond
#' resolution) is based on the DMSP-OLS stable lights product and provides one
#' variable (`Nighttime_light`) representing persistent artificial light at night
#' (Elvidge et al., 1997; Baugh et al., 2010). Nighttime lights can capture
#' infrastructure, urban intensity, and human activity patterns that may not be
#' fully represented by land-cover classifications or human-footprint indices.
#'
#' Burned area (`env.datasets = "burned_area"`; 1 variable with 12 monthly
#' estimates; 900-arcsecond resolution) is based on the ESA Fire Disturbance
#' Climate Change Initiative burned-area product and provides 12 monthly
#' fire-disturbance variables (`Burned_area_01` to `Burned_area_12`). Fire regimes
#' can shape habitat structure, vegetation turnover, resource availability, and
#' landscape heterogeneity, making burned area useful for taxa whose distributions
#' are influenced by disturbance history or post-fire vegetation dynamics (Andela
#' et al., 2017; Chuvieco et al., 2016).
#'
#' Custom rasters (`custom.env.rasters`) allow users to extend the environmental
#' dataset with study-specific predictors. This is useful when biologically
#' important variables are not represented by the built-in datasets, when local
#' products have better spatial or temporal resolution, or when the analysis
#' requires predictors tailored to a specific taxon, region, or hypothesis.
#'
#' @return `NULL`. The function is called for its side effects: writing the
#'   occurrence environmental table to `csv.occurrence.out.file` in `output.dir`
#'   and, when `generate.background.data = TRUE`, writing the background
#'   environmental table to `csv.background.out.file`.
#'
#' @references
#' Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
#'   \emph{Crop evapotranspiration: Guidelines for computing crop water
#'   requirements}. FAO Irrigation and Drainage Paper 56.
#'
#' Andela, N., Morton, D. C., Giglio, L., Chen, Y., van der Werf, G. R.,
#'   Kasibhatla, P. S., DeFries, R. S., Collatz, G. J., Hantson, S.,
#'   Kloster, S., Bachelet, D., Forrest, M., Lasslop, G., Li, F., Mangeon, S.,
#'   Melton, J. R., Yue, C., & Randerson, J. T. (2017). A human-driven decline in
#'   global burned area. \emph{Science}, 356(6345), 1356-1362.
#'   https://doi.org/10.1126/science.aal4108
#'
#' Baugh, K., Hsu, F.-C., Elvidge, C. D., & Zhizhin, M. (2010). Global
#'   inventory modeling and mapping studies. \emph{Proceedings of the Asia-Pacific
#'   Advanced Network}, 30, 78-88.
#'
#' Bennie, J., Huntley, B., Wiltshire, A., Hill, M. O., & Baxter, R. (2008).
#'   Slope, aspect and climate: Spatially explicit and implicit models of
#'   topographic microclimate in chalk grassland. \emph{Ecological Modelling},
#'   216(1), 47-59. https://doi.org/10.1016/j.ecolmodel.2008.04.010
#'
#' Beven, K. J., & Kirkby, M. J. (1979). A physically based, variable
#'   contributing area model of basin hydrology. \emph{Hydrological Sciences
#'   Bulletin}, 24(1), 43-69. https://doi.org/10.1080/02626667909491834
#'
#' Chuvieco, E., Pettinari, M. L., Lizundia-Loiola, J., Storm, T., & Padilla
#'   Parellada, M. (2016). ESA Fire Climate Change Initiative (Fire_cci):
#'   MODIS Fire_cci burned area pixel product, version 5.1.
#'
#' Conrad, O., Bechtel, B., Bock, M., Dietrich, H., Fischer, E., Gerlitz, L.,
#'   Wehberg, J., Wichmann, V., & Böhner, J. (2015). System for Automated
#'   Geoscientific Analyses (SAGA) v. 2.1.4. \emph{Geoscientific Model
#'   Development}, 8, 1991-2007. https://doi.org/10.5194/gmd-8-1991-2015
#'
#' Dobrowski, S. Z. (2011). A climatic basis for microrefugia: The influence of
#'   terrain on climate. \emph{Global Change Biology}, 17(2), 1022-1035.
#'   https://doi.org/10.1111/j.1365-2486.2010.02263.x
#'
#' Dorigo, W. A., Gruber, A., de Jeu, R. A. M., Wagner, W., Stacke, T.,
#'   Loew, A., Albergel, C., Brocca, L., Chung, D., Parinussa, R. M., &
#'   Kidd, R. (2017). ESA CCI Soil Moisture for improved Earth system
#'   understanding: State-of-the art and future directions. \emph{Remote Sensing
#'   of Environment}, 203, 185-215. https://doi.org/10.1016/j.rse.2017.07.001
#'
#' Dormann, C. F., Elith, J., Bacher, S., Buchmann, C., Carl, G., Carré, G.,
#'   Marquéz, J. R. G., Gruber, B., Lafourcade, B., Leitão, P. J.,
#'   Münkemüller, T., McClean, C., Osborne, P. E., Reineking, B., Schröder, B.,
#'   Skidmore, A. K., Zurell, D., & Lautenbach, S. (2013). Collinearity: A
#'   review of methods to deal with it and a simulation study evaluating their
#'   performance. \emph{Ecography}, 36(1), 27-46.
#'   https://doi.org/10.1111/j.1600-0587.2012.07348.x
#'
#' Elith, J., & Leathwick, J. R. (2009). Species distribution models: Ecological
#'   explanation and prediction across space and time. \emph{Annual Review of
#'   Ecology, Evolution, and Systematics}, 40, 677-697.
#'   https://doi.org/10.1146/annurev.ecolsys.110308.120159
#'
#' Elvidge, C. D., Baugh, K. E., Kihn, E. A., Kroehl, H. W., Davis, E. R., &
#'   Davis, C. W. (1997). Relation between satellite observed visible-near
#'   infrared emissions, population, economic activity and electric power
#'   consumption. \emph{International Journal of Remote Sensing}, 18(6),
#'   1373-1379. https://doi.org/10.1080/014311697218485
#'
#' Fick, S. E., & Hijmans, R. J. (2017). WorldClim 2: New 1-km spatial
#'   resolution climate surfaces for global land areas. \emph{International
#'   Journal of Climatology}, 37(12), 4302-4315.
#'   https://doi.org/10.1002/joc.5086
#'
#' Gruber, A., Scanlon, T., van der Schalie, R., Wagner, W., & Dorigo, W.
#'   (2019). Evolution of the ESA CCI Soil Moisture climate data records and
#'   their underlying merging methodology. \emph{Earth System Science Data}, 11,
#'   717-739. https://doi.org/10.5194/essd-11-717-2019
#'
#' Guisan, A., & Zimmermann, N. E. (2000). Predictive habitat distribution
#'   models in ecology. \emph{Ecological Modelling}, 135(2-3), 147-186.
#'   https://doi.org/10.1016/S0304-3800(00)00354-9
#'
#' Hengl, T., & Evans, I. S. (2009). Mathematical and digital models of the land
#'   surface. In T. Hengl & H. I. Reuter (Eds.), \emph{Geomorphometry: Concepts,
#'   software, applications} (pp. 31-63). Elsevier.
#'   https://doi.org/10.1016/S0166-2481(08)00002-0
#'
#' Hijmans, R. J., Cameron, S. E., Parra, J. L., Jones, P. G., & Jarvis, A.
#'   (2005). Very high resolution interpolated climate surfaces for global land
#'   areas. \emph{International Journal of Climatology}, 25(15), 1965-1978.
#'   https://doi.org/10.1002/joc.1276
#'
#' Hillel, D. (1998). \emph{Environmental soil physics}. Academic Press.
#'
#' Hufkens, K., Basler, D., Milliman, T., Melaas, E. K., & Richardson, A. D.
#'   (2018). An integrated phenology modelling framework in R. \emph{Methods in
#'   Ecology and Evolution}, 9(5), 1276-1285.
#'   https://doi.org/10.1111/2041-210X.12970
#'
#' Huete, A., Didan, K., Miura, T., Rodriguez, E. P., Gao, X., & Ferreira, L. G.
#'   (2002). Overview of the radiometric and biophysical performance of the MODIS
#'   vegetation indices. \emph{Remote Sensing of Environment}, 83(1-2), 195-213.
#'   https://doi.org/10.1016/S0034-4257(02)00096-2
#'
#' Jackson, R. B., Canadell, J., Ehleringer, J. R., Mooney, H. A., Sala, O. E.,
#'   & Schulze, E. D. (1996). A global analysis of root distributions for
#'   terrestrial biomes. \emph{Oecologia}, 108, 389-411.
#'   https://doi.org/10.1007/BF00333714
#'
#' Kearney, M., & Porter, W. (2009). Mechanistic niche modelling: Combining
#'   physiological and spatial data to predict species' ranges. \emph{Ecology
#'   Letters}, 12(4), 334-350.
#'   https://doi.org/10.1111/j.1461-0248.2008.01277.x
#'
#' Lang, N., Jetz, W., Schindler, K., & Wegner, J. D. (2023). A high-resolution
#'   canopy height model of the Earth. \emph{Nature Ecology & Evolution}, 7,
#'   1778-1789. https://doi.org/10.1038/s41559-023-02206-6
#'
#' MacDonald, Z. G., Beninde, J., Matsunaga, K., Zhou, B., Gillespie, T. W., &
#'   Shaffer, H. B. (2025). Species distribution modeling for conservation
#'   science: New predictor layers, reproducible code, and an evaluation of
#'   California protected areas. \emph{bioRxiv}.
#'   https://doi.org/10.1101/2025.01.23.634559
#'
#' McCune, B., & Keon, D. (2002). Equations for potential annual direct incident
#'   radiation and heat load. \emph{Journal of Vegetation Science}, 13(4),
#'   603-606. https://doi.org/10.1111/j.1654-1103.2002.tb02087.x
#'
#' Novick, K. A., Ficklin, D. L., Stoy, P. C., Williams, C. A., Bohrer, G.,
#'   Oishi, A. C., Papuga, S. A., Blanken, P. D., Noormets, A., Sulman, B. N.,
#'   Scott, R. L., Wang, L., & Phillips, R. P. (2016). The increasing importance
#'   of atmospheric demand for ecosystem water and carbon fluxes. \emph{Nature
#'   Climate Change}, 6, 1023-1027. https://doi.org/10.1038/nclimate3114
#'
#' Poggio, L., de Sousa, L. M., Batjes, N. H., Heuvelink, G. B. M., Kempen, B.,
#'   Ribeiro, E., & Rossiter, D. (2021). SoilGrids 2.0: Producing soil
#'   information for the globe with quantified spatial uncertainty.
#'   \emph{SOIL}, 7(1), 217-240. https://doi.org/10.5194/soil-7-217-2021
#'
#' Potapov, P., Li, X., Hernandez-Serna, A., Tyukavina, A., Hansen, M. C.,
#'   Kommareddy, A., Pickens, A., Turubanova, S., Tang, H., Silva, C. E.,
#'   Armston, J., Dubayah, R., Blair, J. B., & Hofton, M. (2021). Mapping global
#'   forest canopy height through integration of GEDI and Landsat data.
#'   \emph{Remote Sensing of Environment}, 253, 112165.
#'   https://doi.org/10.1016/j.rse.2020.112165
#'
#' Prajzlerová, D., Barták, V., Balej, P., & Šímová, P. (2025). The time of
#'   acquisition of multispectral predictors matters: The role of seasonality in
#'   bird species distribution models. \emph{Ecography}, 2025, e07935.
#'   https://doi.org/10.1002/ecog.07935
#'
#' Preimesberger, W., Scanlon, T., Su, C.-H., Gruber, A., & Dorigo, W. (2021).
#'   Homogenization of structural breaks in the global ESA CCI soil moisture
#'   multisatellite climate data record. \emph{IEEE Transactions on Geoscience
#'   and Remote Sensing}, 59(4), 2845-2862.
#'   https://doi.org/10.1109/TGRS.2020.3012896
#'
#' Soberón, J. (2007). Grinnellian and Eltonian niches and geographic
#'   distributions of species. \emph{Ecology Letters}, 10(12), 1115-1123.
#'   https://doi.org/10.1111/j.1461-0248.2007.01107.x
#'
#' Soberón, J., & Peterson, A. T. (2005). Interpretation of models of fundamental
#'   ecological niches and species' distributional areas. \emph{Biodiversity
#'   Informatics}, 2. https://doi.org/10.17161/bi.v2i0.4
#'
#' Sørensen, R., Zinko, U., & Seibert, J. (2006). On the calculation of the
#'   topographic wetness index: Evaluation of different methods based on field
#'   observations. \emph{Hydrology and Earth System Sciences}, 10, 101-112.
#'   https://doi.org/10.5194/hess-10-101-2006
#'
#' Title, P. O., & Bemmels, J. B. (2018). ENVIREM: An expanded set of bioclimatic
#'   and topographic variables increases flexibility and improves performance of
#'   ecological niche modeling. \emph{Ecography}, 41(2), 291-307.
#'   https://doi.org/10.1111/ecog.02880
#'
#' Turek, M. E., Poggio, L., Batjes, N. H., Armindo, R. A.,
#'   de Jong van Lier, Q., de Sousa, L. M., & Heuvelink, G. B. M. (2023).
#'   Global mapping of volumetric water retention at 100, 330 and 15 000 cm
#'   suction using the WoSIS database. \emph{International Soil and Water
#'   Conservation Research}, 11(2), 225-239.
#'   https://doi.org/10.1016/j.iswcr.2022.08.001
#'
#' Venter, O., Sanderson, E. W., Magrach, A., Allan, J. R., Beher, J.,
#'   Jones, K. R., Possingham, H. P., Laurance, W. F., Wood, P., Fekete, B. M.,
#'   Levy, M. A., & Watson, J. E. M. (2016). Sixteen years of change in the
#'   global terrestrial human footprint and implications for biodiversity
#'   conservation. \emph{Nature Communications}, 7, 12558.
#'   https://doi.org/10.1038/ncomms12558
#'
#' Wang, T., Hamann, A., Spittlehouse, D., & Carroll, C. (2016). Locally
#'   downscaled and spatially customizable climate data for historical and future
#'   periods for North America. \emph{PLOS ONE}, 11(6), e0156720.
#'   https://doi.org/10.1371/journal.pone.0156720
#'
#' Weil, R. R., & Brady, N. C. (2016). \emph{The nature and properties of soils}
#'   (15th ed.). Pearson.
#'
#' Wilson, M. F. J., O'Connell, B., Brown, C., Guinan, J. C., & Grehan, A. J.
#'   (2007). Multiscale terrain analysis of multibeam bathymetry data for habitat
#'   mapping on the continental slope. \emph{Marine Geodesy}, 30(1-2), 3-35.
#'   https://doi.org/10.1080/01490410701295962
#'
#' Zimmermann, N. E., Yoccoz, N. G., Edwards, T. C., Meier, E. S.,
#'   Thuiller, W., Guisan, A., Schmatz, D. R., & Pearman, P. B. (2009).
#'   Climatic extremes improve predictions of spatial patterns of tree species.
#'   \emph{Proceedings of the National Academy of Sciences}, 106(Supplement 2),
#'   19723-19728. https://doi.org/10.1073/pnas.0901643106
#'
#' @export
extract.env.and.background <- function(occurrence.data, #input data.frame with coords (rownames need to be unique IDs)
                                       longitude.col = "Longitude", #column name for longitude
                                       latitude.col = "Latitude", #column name for latitude
                                       generate.background.data = FALSE, #whether to sample background points
                                       N.background.points = 100000, #number of background points to draw (default: 100000)
                                       remove.hydrolakes.background = FALSE, #erase lakes from background mask using HydroLAKES
                                       buffer.km, #buffer in km around convex hull for accessible area
                                       landmask.largest.N.pieces = 5, #keep N largest land polygons to form mask (recommended: 5)
                                       csv.occurrence.out.file, #output CSV name for occurrence+env values
                                       csv.background.out.file = NULL, #output CSV name for background+env values
                                       output.dir, #main output directory
                                       intermediate.files.dir = "Intermediate_files", #subfolder for cached intermediate CSVs
                                       env.datasets = NULL, #optional vector of built-in environmental datasets to extract
                                       CRS.occurrences = "EPSG:4326", #CRS string for input coordinates (default: "EPSG:4326")
                                       overwrite = FALSE, #overwrite existing occurrence/background CSVs
                                       delete.intermediate.files.folders = FALSE, #delete intermediate raster folders after run
                                       redownload.rasters = FALSE, #whether to force re-download of rasters/archives
                                       custom.env.rasters = NULL, #vector of custom raster paths/URLs (.tif/.tiff, folder containing .tif/.tiff files, or .zip containing .tif/.tiff files)
                                       custom.env.rasters.names = NULL, #dataset names for custom rasters
                                       custom.env.rasters.variable.names = NULL, #layer names to assign to custom rasters
                                       custom.env.rasters.crs = NULL, #optional CRS string (or one per custom raster) used only if raster CRS is missing
                                       custom.raster.min.size.mb = 5, #minimum size check (MB) for downloaded rasters
                                       seed = 1, #random seed
                                       verbose = TRUE #whether to print progress messages
) {

  # Validate datasets
  datasets_requested <- unique(env.datasets)
  valid_datasets <- c("elevation",
                      "ClimateNA",
                      "EVI",
                      "terrain",
                      "ENVIREM",
                      "footprint",
                      "landcover",
                      "soil",
                      "forest_height",
                      "atmosphere",
                      "nightlight",
                      "burned_area",
                      "snow_water_equivalent",
                      "daylength",
                      "soil_moisture")
  invalid <- setdiff(datasets_requested, valid_datasets)
  if (length(invalid)) stop(paste0("Invalid env.datasets entries: ", paste(invalid, collapse = ", "), " - choose from: ", paste(valid_datasets, collapse = ", ")))

  # Validate input arguments
  if (!is.data.frame(occurrence.data)) stop("occurrence.data must be a data.frame")
  if (!all(c(longitude.col, latitude.col) %in% names(occurrence.data))) stop("occurrence.data must contain longitude.col and latitude.col")
  if (is.null(rownames(occurrence.data))) stop("occurrence.data requires rownames as unique IDs - set rownames(occurrence.data)")
  if (anyDuplicated(rownames(occurrence.data))) stop("occurrence.data requires unique rownames - ensure they are unique before running")
  if (any(!is.finite(occurrence.data[[longitude.col]])) || any(!is.finite(occurrence.data[[latitude.col]]))) stop("Longitude/Latitude contain non-finite values (NA/Inf)")
  crs_occ <- tryCatch(sf::st_crs(CRS.occurrences), error = function(e) NULL)
  if (is.null(crs_occ)) stop("CRS.occurrences must be a valid CRS string (default: EPSG:4326)")
  if (isTRUE(sf::st_is_longlat(crs_occ)) && (any(occurrence.data[[longitude.col]] < -180 | occurrence.data[[longitude.col]] > 180) || any(occurrence.data[[latitude.col]] < -90 | occurrence.data[[latitude.col]] > 90))) stop("Longitude/Latitude are out of valid ranges (-180..180, -90..90)")
  if (!is.logical(generate.background.data)) stop("generate.background.data must be TRUE or FALSE")
  if (generate.background.data && nrow(occurrence.data) < 3) stop("Need at least three unique occurrence points to generate convex hull for accessible area")
  if (!is.numeric(N.background.points) || N.background.points <= 0) stop("N.background.points must be a positive numeric value (recommended: 10000)")
  if (!is.logical(remove.hydrolakes.background)) stop("remove.hydrolakes.background must be TRUE or FALSE")
  if (!is.numeric(buffer.km) || buffer.km <= 0) stop("buffer.km (dispersal distance) must be a positive numeric value in km")
  if (!is.numeric(landmask.largest.N.pieces) || landmask.largest.N.pieces <= 0) stop("landmask.largest.N.pieces must be a positive integer (recommended: 5)")
  if (!is.character(csv.occurrence.out.file) || !nzchar(csv.occurrence.out.file)) stop("csv.occurrence.out.file must be a valid non-empty file name ending with .csv")
  if (generate.background.data && (!is.character(csv.background.out.file) || !nzchar(csv.background.out.file))) stop("csv.background.out.file must be a valid non-empty file name ending with .csv when generate.background.data = TRUE")
  if (!is.character(output.dir) || length(output.dir) != 1 || !nzchar(output.dir)) stop("output.dir must be a non-empty character string")
  if (!is.character(intermediate.files.dir) || length(intermediate.files.dir) != 1 || !nzchar(intermediate.files.dir)) stop("intermediate.files.dir must be a non-empty character string")
  if (!is.character(CRS.occurrences) || length(CRS.occurrences) != 1 || !nzchar(CRS.occurrences)) stop("CRS.occurrences must be a valid CRS string (default: EPSG:4326)")
  if (!is.logical(overwrite)) stop("overwrite must be TRUE or FALSE")
  if (!is.logical(delete.intermediate.files.folders)) stop("delete.intermediate.files.folders must be TRUE or FALSE")
  if (!is.logical(redownload.rasters)) stop("redownload.rasters must be TRUE or FALSE")
  if (!is.null(custom.env.rasters)) {
    if (!is.character(custom.env.rasters)) stop("custom.env.rasters must be a character vector of file paths or URLs")
    if (!is.null(custom.env.rasters.names) && length(custom.env.rasters.names) != length(custom.env.rasters)) stop("custom.env.rasters.names (dataset names) must be NULL or have same length as custom.env.rasters")
    if (!is.null(custom.env.rasters.variable.names)) {
      if (is.character(custom.env.rasters.variable.names)) {
        if (length(custom.env.rasters) == 1) {
          custom.env.rasters.variable.names <- list(custom.env.rasters.variable.names)
        } else if (length(custom.env.rasters.variable.names) == length(custom.env.rasters)) {
          custom.env.rasters.variable.names <- as.list(custom.env.rasters.variable.names)
        } else {
          stop("custom.env.rasters.variable.names must be either a list with one character vector per raster, or a character vector of length equal to custom.env.rasters for single-layer rasters")
        }
      }
      if (!is.list(custom.env.rasters.variable.names) || length(custom.env.rasters.variable.names) != length(custom.env.rasters)) stop("custom.env.rasters.variable.names must be NULL, or a list with one character vector per raster")
      if (any(vapply(custom.env.rasters.variable.names, function(x) !is.character(x) || length(x) < 1, logical(1)))) stop("Each custom.env.rasters.variable.names entry must be a non-empty character vector")
    }
    if (!is.null(custom.env.rasters.crs)) {
      if (!is.character(custom.env.rasters.crs)) stop("custom.env.rasters.crs must be a character vector of CRS strings")
      if (!(length(custom.env.rasters.crs) %in% c(1, length(custom.env.rasters)))) stop("custom.env.rasters.crs must have length 1 or the same length as custom.env.rasters")
    }
    for (i in seq_along(custom.env.rasters)) {
      src <- custom.env.rasters[[i]]
      if (grepl("^https?://", src)) {
        if (!grepl("\\.(tif|tiff|zip)($|\\?)", src, ignore.case = TRUE)) stop("custom.env.rasters[", i, "] must point to a .tif/.tiff file or a .zip containing .tif/.tiff files: ", src)
      } else {
        if (!file.exists(src)) stop("custom.env.rasters[", i, "] local file not found: ", src)
        if (dir.exists(src)) {
          tif_in_dir <- list.files(src, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE)
          if (length(tif_in_dir) == 0) stop("custom.env.rasters[", i, "] is a directory but contains no .tif/.tiff files: ", src)
        } else if (grepl("\\.zip$", src, ignore.case = TRUE)) {
          zip_listing <- tryCatch(utils::unzip(src, list = TRUE), error = function(e) NULL)
          if (is.null(zip_listing)) stop("custom.env.rasters[", i, "] is a zip file but could not be read: ", src)
          tif_in_zip <- zip_listing$Name[grepl("\\.(tif|tiff)$", zip_listing$Name, ignore.case = TRUE)]
          if (length(tif_in_zip) == 0) stop("custom.env.rasters[", i, "] is a zip file but contains no .tif/.tiff files: ", src)
        } else if (!grepl("\\.(tif|tiff)$", src, ignore.case = TRUE)) {
          stop("custom.env.rasters[", i, "] must be a .tif/.tiff file, a folder containing .tif/.tiff files, or a .zip containing .tif/.tiff files: ", src)
        }
      }
    }
  }
  if (!is.numeric(seed) || length(seed) != 1) stop("seed must be a single numeric value (default: 1)")
  if (!is.logical(verbose)) stop("verbose must be TRUE or FALSE")
  if ("ClimateNA" %in% datasets_requested && !"elevation" %in% datasets_requested) stop("ClimateNA requires elevation dataset - add elevation to env.datasets")

  # Specify continental subset extents
  europe_extent_general <- terra::ext(-25, 41, 35, 71)
  asia_extent_general <- terra::ext(39, 177, -9, 80)
  eurasia_extent_general <- terra::ext(-25, 177, -9, 80)
  northamerica_extent_general <- terra::ext(-170, -50, 10, 76)
  southamerica_extent_general <- terra::ext(-90, -30, -59, 12)
  africa_extent_general <- terra::ext(-20, 51, -34, 37)
  australia_extent_general <- terra::ext(110, 176, -49, -7)
  indopacific_extent_general <- terra::ext(39, 177, -49, 80)
  holarctic_extent_general <- terra::ext(-170, 177, -9, 80)
  newworld_extent_general <- terra::ext(-170, -30, -59, 76)
  oldworld_extent_general <- terra::ext(-25, 177, -34, 80)

  # Create directories if needed
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)
  intermediate_files_dir <- file.path(output.dir, intermediate.files.dir)
  if (!dir.exists(intermediate_files_dir)) dir.create(intermediate_files_dir, recursive = TRUE)
  rasters.dir <- file.path(output.dir, "Env_rasters")
  if (!dir.exists(rasters.dir)) dir.create(rasters.dir, recursive = TRUE)

  # Create function to robustly download rasters or generic files (ZIP, CSV, etc.)
  robust.download.raster <- function(url, dest, max_attempts = 5, min_size_mb = 10, ...) {
    dest_ext <- tolower(tools::file_ext(dest))
    is_raster_ext <- dest_ext %in% c("tif", "tiff")
    old_timeout <- getOption("timeout")
    options(timeout = max(old_timeout, 8 * 60 * 60))
    on.exit(options(timeout = old_timeout), add = TRUE)
    for (attempt in 1:max_attempts) {
      if (file.exists(dest) && file.size(dest) < min_size_mb * 1024^2) file.remove(dest)
      download_try <- try(utils::download.file(url, destfile = dest, mode = "wb", quiet = TRUE, method = "libcurl"), silent = TRUE)
      if (inherits(download_try, "try-error")) message("Download failed: ", attr(download_try, "condition")$message)
      valid <- FALSE
      if (file.exists(dest) && file.size(dest) > min_size_mb * 1024^2) {
        if (is_raster_ext) {
          valid <- tryCatch({!is.null(dim(terra::rast(dest)))}, error = function(e) FALSE)
        } else if (dest_ext == "zip") {
          valid <- tryCatch(nrow(utils::unzip(dest, list = TRUE)) > 0, error = function(e) FALSE)
        } else {
          valid <- TRUE
        }
      }
      if (valid) return(dest)
      if (file.exists(dest)) file.remove(dest)
      Sys.sleep(10)
    }
    stop("Failed to download a valid file from ", url, " - try again or download manually")
  }

  # Create function to download and load Copernicus GLO-90 elevation tiles
  download.and.load.elevation.tile <- function(continent_name, url_elevation_zip, elevation_dir, min_size_mb = 1800, redownload = FALSE, verbose = TRUE) {
    elevation_zip_file <- file.path(elevation_dir, paste0("Copernicus_GLO90_", continent_name, "_250m.zip"))
    elevation_tif_file <- file.path(elevation_dir, paste0("Copernicus_GLO90_", continent_name, "_250m.tif"))
    if (file.exists(elevation_tif_file) && file.size(elevation_tif_file) > min_size_mb * 1e6 && !redownload) {
      if (verbose) message("Elevation raster already present - skipping download")
    } else {
      if (file.exists(elevation_zip_file) && (file.size(elevation_zip_file) < (min_size_mb - 100) * 1e6 || redownload)) {
        try(unlink(elevation_zip_file), silent = TRUE)
      }
      if (!file.exists(elevation_zip_file) || file.size(elevation_zip_file) < (min_size_mb - 100) * 1e6 || redownload) {
        ok <- FALSE
        old_timeout <- getOption("timeout")
        options(timeout = max(old_timeout, 8 * 60 * 60))
        on.exit(options(timeout = old_timeout), add = TRUE)
        for (i in 1:5) {
          download_try <- try(utils::download.file(url_elevation_zip, destfile = elevation_zip_file, mode = "wb", quiet = TRUE, method = "libcurl"), silent = TRUE)
          if (inherits(download_try, "try-error")) message("Download failed: ", attr(download_try, "condition")$message)
          if (file.exists(elevation_zip_file) && file.info(elevation_zip_file)$size > min_size_mb * 1e6) {
            ok <- TRUE
            break
          }
          if (file.exists(elevation_zip_file)) try(unlink(elevation_zip_file), silent = TRUE)
          Sys.sleep(5)
        }
        if (!ok) stop("Failed to download ", continent_name, " elevation zip from Zenodo")
      }
      utils::unzip(elevation_zip_file, exdir = elevation_dir, overwrite = TRUE)
      if (file.exists(elevation_tif_file) && file.info(elevation_tif_file)$size > min_size_mb * 1e6) {
        try(unlink(elevation_zip_file), silent = TRUE)
      } else {
        stop("Extraction failed or .tif missing after unzip for ", continent_name)
      }
    }
    terra::rast(elevation_tif_file)
  }

  # Extract raster values, cache results as CSVs, and reuse if already present
  extract.and.cache.env.dataset <- function(dataset_name, #dataset name
                                            raster_object, #SpatRaster or SpatRaster stack to extract values from
                                            coord_env, #SpatVector of occurrence coordinates
                                            coord_bg, #SpatVector of background coordinates (optional)
                                            environmental_dataset, #data.frame of occurrence data (to append)
                                            background.data, #data.frame of background data (to append)
                                            output.dir, #output directory for intermediate CSVs
                                            overwrite = FALSE, #whether to overwrite existing CSVs
                                            generate.background.data = TRUE, #whether to extract background values
                                            verbose = TRUE) { #print progress messages
    occ_csv_file <- file.path(intermediate_files_dir, paste0(dataset_name, "_extracted_occurrence.csv"))
    bg_csv_file <- file.path(intermediate_files_dir, paste0(dataset_name, "_extracted_background.csv"))
    if (file.exists(occ_csv_file) && (!generate.background.data || file.exists(bg_csv_file)) && !overwrite) {
      if (verbose) message(dataset_name, " extraction already saved - loading from intermediate CSV files")
      extracted_occurrences <- read.csv(occ_csv_file, row.names = 1, check.names = FALSE)
      extracted_occurrences <- extracted_occurrences[rownames(environmental_dataset), , drop = FALSE]
      extracted_background <- if (generate.background.data) read.csv(bg_csv_file, row.names = 1, check.names = FALSE) else NULL
      if (!is.null(extracted_background) && !is.null(background.data)) extracted_background <- extracted_background[rownames(background.data), , drop = FALSE]
    } else {
      raster_names <- names(raster_object)
      if (is.null(raster_names) || length(raster_names) != terra::nlyr(raster_object) || any(!nzchar(raster_names))) {
        raster_names <- paste0(dataset_name, "_", seq_len(terra::nlyr(raster_object)))
        names(raster_object) <- raster_names
      }
      extracted_occurrences <- terra::extract(raster_object, coord_env, ID = FALSE)
      extracted_occurrences <- as.data.frame(extracted_occurrences, check.names = FALSE)
      colnames(extracted_occurrences) <- raster_names
      rownames(extracted_occurrences) <- rownames(environmental_dataset)
      write.csv(extracted_occurrences, occ_csv_file, row.names = TRUE)
      if (!is.null(background.data) && generate.background.data) {
        extracted_background <- terra::extract(raster_object, coord_bg, ID = FALSE)
        extracted_background <- as.data.frame(extracted_background, check.names = FALSE)
        colnames(extracted_background) <- raster_names
        rownames(extracted_background) <- rownames(background.data)
        write.csv(extracted_background, bg_csv_file, row.names = TRUE)
      } else {
        extracted_background <- NULL
      }
    }
    environmental_dataset <- cbind(environmental_dataset, extracted_occurrences)
    if (!is.null(background.data) && generate.background.data && !is.null(extracted_background)) {
      background.data <- cbind(background.data, extracted_background)
    }
    suppressWarnings(terra::tmpFiles(remove = TRUE))
    invisible(gc())
    return(list(environmental_dataset = environmental_dataset, background.data = background.data))
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Set parameters
  occurrence_output_file <- file.path(output.dir, csv.occurrence.out.file)
  environmental_dataset <- occurrence.data
  rownames(environmental_dataset) <- rownames(occurrence.data)
  idx <- match(c(longitude.col, latitude.col), names(environmental_dataset))
  names(environmental_dataset)[idx] <- c("Longitude", "Latitude")
  base_ids <- rownames(environmental_dataset)
  existing <- if (file.exists(occurrence_output_file) && !overwrite) {
    tryCatch(read.csv(occurrence_output_file, row.names = 1, check.names = FALSE), error = function(e) NULL)
  } else NULL
  CRS_all <- CRS.occurrences

  # Set dataset counter
  total_datasets <- length(datasets_requested) + if (!is.null(custom.env.rasters)) length(custom.env.rasters) else 0
  if (total_datasets == 0) stop("No environmental datasets specified - include at least one in env.datasets or custom.env.rasters", call. = FALSE)
  counter <- 0 #sequential counter for ordered progress messages

  # Validate terrain requirements (requires elevation and whitebox R package)
  if ("terrain" %in% datasets_requested) {
    if (!"elevation" %in% datasets_requested) stop("terrain requires elevation dataset - add elevation to env.datasets")
    if (!requireNamespace("whitebox", quietly = TRUE)) stop("Package 'whitebox' is required for this function.")
  }

  # Generate study area and accessible area (background area)
  if (verbose) message("")
  if (verbose) message("-- Setting study area and accessible area --")
  presence_sf <- sf::st_as_sf(occurrence.data[, c(longitude.col, latitude.col)],
                              coords = c(longitude.col, latitude.col),
                              crs = sf::st_crs(CRS_all)) #presence points as sf
  if (sf::st_is_longlat(presence_sf)) { #if geographic CRS
    equal_area_extent_5070_ConusAlbers <- terra::ext(-2500000, 2500000, -2000000, 3500000)
    study_area_5070_ConusAlbers <- try(terra::project(terra::as.polygons(terra::ext(terra::vect(presence_sf))), "EPSG:5070"), silent = TRUE)
    is_within_5070_ConusAlbers <- !inherits(study_area_5070_ConusAlbers, "try-error") && terra::relate(study_area_5070_ConusAlbers, equal_area_extent_5070_ConusAlbers, relation = "within")[1]
    if (is_within_5070_ConusAlbers) {
      equal_area_crs <- "EPSG:5070" #NAD83 / Conus Albers - optimized for North America
    } else { #determine centroid and use locally centered Lambert Azimuthal Equal-Area
      s2_state <- getOption("sf_use_s2", TRUE)
      options(sf_use_s2 = FALSE) #disable s2
      invisible(suppressMessages(suppressWarnings(sf::sf_use_s2(FALSE))))
      suppressMessages(suppressWarnings(centroid_coords <- sf::st_coordinates(sf::st_centroid(sf::st_union(presence_sf)))))
      options(sf_use_s2 = s2_state) #restore s2 setting
      centroid_lon <- centroid_coords[1]
      centroid_lat <- centroid_coords[2]
      if (centroid_lat > -60 && centroid_lat < 60) {
        equal_area_crs <- sprintf("+proj=laea +lat_0=%f +lon_0=%f +datum=WGS84 +units=m +no_defs", centroid_lat, centroid_lon) #LAEA centered on data
      } else {
        equal_area_crs <- "EPSG:6933" #Global equal-area fallback for polar/extreme regions
      }
    }
    presence_equal_area <- sf::st_transform(presence_sf, equal_area_crs) #project to equal-area
    hull_equal_area <- sf::st_convex_hull(sf::st_union(presence_equal_area)) #convex hull in equal-area
    accessible_area_equal_area <- sf::st_buffer(hull_equal_area, dist = buffer.km * 1000) #buffer in meters
    accessible_area <- sf::st_transform(accessible_area_equal_area, CRS_all) #back to original CRS
    raster_25km_resolution <- c("soil_moisture", "burned_area")
    raster_1km_resolution <- c("ENVIREM", "footprint", "landcover", "soil", "atmosphere", "nightlight", "daylength", "snow_water_equivalent")
    raster_250m_resolution <- c("elevation", "EVI", "forest_height", "snow_water_equivalent")
    raster_scalefree_resolution <- "ClimateNA"
    min_raster_resolution_km <- if (any(datasets_requested %in% raster_25km_resolution)) 25 else
      if (any(datasets_requested %in% raster_1km_resolution)) 1 else
        if (any(datasets_requested %in% raster_250m_resolution)) 0.25 else
          if (any(datasets_requested %in% raster_scalefree_resolution)) 0.1 else
            1
    study_area_equal_area <- sf::st_buffer(accessible_area_equal_area, dist = (min_raster_resolution_km + 5) * 1000) #buffer background polygon by (min raster res + 2 km)
    study_area <- sf::st_transform(study_area_equal_area, CRS_all) #transform back to geographic CRS for cropping rasters later
    study_area_vect <- terra::vect(study_area) #convert sf study area to SpatVector
  } else { #already projected CRS
    raster_25km_resolution <- c("soil_moisture", "burned_area")
    raster_1km_resolution <- c("ENVIREM", "footprint", "landcover", "soil", "atmosphere", "nightlight", "daylength", "snow_water_equivalent")
    raster_250m_resolution <- c("elevation", "EVI", "forest_height", "snow_water_equivalent")
    raster_scalefree_resolution <- "ClimateNA"
    min_raster_resolution_km <- if (any(datasets_requested %in% raster_25km_resolution)) 25 else
      if (any(datasets_requested %in% raster_1km_resolution)) 1 else
        if (any(datasets_requested %in% raster_250m_resolution)) 0.25 else
          if (any(datasets_requested %in% raster_scalefree_resolution)) 0.1 else
            1
    if (identical(sf::st_crs(presence_sf)$units_gdal, "metre")) {
      hull_projected <- sf::st_convex_hull(sf::st_union(presence_sf)) #convex hull in native CRS (meters)
      accessible_area <- sf::st_buffer(hull_projected, dist = buffer.km * 1000) #buffer in meters
      study_area <- sf::st_buffer(accessible_area, dist = (min_raster_resolution_km + 5) * 1000) #buffer study area in meters
    } else {
      extent_north_america_background <- terra::ext(-170, -50, 10, 76)
      is_north_america_background <- terra::relate(terra::as.polygons(terra::ext(terra::vect(presence_sf))), extent_north_america_background, relation = "within")[1]
      equal_area_crs <- if (is_north_america_background) "EPSG:5070" else "EPSG:6933" #project to meters first
      presence_m <- sf::st_transform(presence_sf, equal_area_crs)
      hull_m <- sf::st_convex_hull(sf::st_union(presence_m))
      accessible_area_m <- sf::st_buffer(hull_m, dist = buffer.km * 1000) #buffer in meters
      study_area_m <- sf::st_buffer(accessible_area_m, dist = (min_raster_resolution_km + 5) * 1000) #buffer study area in meters
      accessible_area <- sf::st_transform(accessible_area_m, CRS_all) #back to original CRS
      study_area <- sf::st_transform(study_area_m, CRS_all) #back to original CRS
    }
    study_area_vect <- terra::vect(study_area) #convert sf study area to SpatVector
  }

  # Validate extent-dependent dataset requirements
  if ("snow_water_equivalent" %in% datasets_requested) {
    north_america_extent_daymetr  <- terra::ext(-178.1250, -52.8750, 14.0625, 82.9375)
    is_north_america_SWE <- terra::relate(study_area_vect, north_america_extent_daymetr, relation = "within")[1]
    if (!is_north_america_SWE) stop("snow_water_equivalent is only available for North America - remove from env.datasets")
  }

  if ("daylength" %in% datasets_requested) {
    north_america_extent_daymetr <- terra::ext(-178.1250, -52.8750, 14.0625, 82.9375)
    is_north_america_daylength <- terra::relate(study_area_vect, north_america_extent_daymetr, relation = "within")[1]
    if (!is_north_america_daylength) stop("daylength is only available for North America - remove from env.datasets")
  }

  if ("ClimateNA" %in% datasets_requested) {
    north_america_climatena_extent <- terra::ext(-179.133, -53.067, 14.075, 82.914)
    is_north_america_climatena <- terra::relate(study_area_vect, north_america_climatena_extent, relation = "within")[1]
    if (!is_north_america_climatena) stop("ClimateNA is only available for North America - remove from env.datasets")
  }

  # Create land mask
  if (("terrain" %in% datasets_requested) || generate.background.data) {
    if (verbose) message("")
    if (verbose) message("-- Creating land mask --")
    landmask_dir <- file.path(rasters.dir, "landmask")
    if (!dir.exists(landmask_dir)) dir.create(landmask_dir, recursive = TRUE)
    gadm_file <- file.path(landmask_dir, "gadm36_adm0_r5_pk.rds")
    gadm_url <- "https://geodata.ucdavis.edu/gadm/gadm3.6/gadm36_adm0_r5_pk.rds"
    if (!file.exists(gadm_file) || redownload.rasters) {
      tryCatch(
        utils::download.file(gadm_url, gadm_file, mode = "wb", quiet = TRUE),
        error = function(e) {
          stop("Failed to download GADM land mask from geodata.ucdavis.edu - this is usually a temporary server or network issue.\n",
               "Retry later or download file manually and place it at:\n",
               gadm_file,
               call. = FALSE
          )
        }
      )
    }
    suppressMessages(suppressWarnings(invisible(capture.output(countries_full <- readRDS(gadm_file)))))
    if (inherits(countries_full, "PackedSpatVector")) countries_full <- terra::unwrap(countries_full)
    if (!inherits(countries_full, "SpatVector")) countries_full <- terra::vect(countries_full)
    overlap_idx <- terra::relate(terra::vect(study_area), countries_full, relation = "intersects")[1, ]
    geodata.countries <- sort(unique(na.omit(countries_full$GID_0[which(overlap_idx)])))
    countries_crop <- terra::crop(countries_full, terra::vect(study_area))
    na_countries <- countries_crop[countries_crop$GID_0 %in% geodata.countries, ]
    if (nrow(na_countries) == 0) stop("No countries overlap study area - check CRS and coordinates")
    na_polys <- terra::disagg(na_countries)
    if (nrow(na_polys) == 0) stop("No polygons available to build land mask")
    areas <- terra::expanse(na_polys, unit = "km")
    land_mask <- na_polys[order(areas, decreasing = TRUE)[seq_len(min(landmask.largest.N.pieces, nrow(na_polys)))],]
  }

  # Remove lakes from land mask (using HydroLAKES)
  if (("terrain" %in% datasets_requested) || generate.background.data) {
    invisible(gc())
    if (remove.hydrolakes.background) {
      if (verbose) message("")
      if (verbose) message("-- Removing HydroLAKES polygons from land mask --")
      is_europe_hydrolakes <- terra::relate(study_area_vect, europe_extent_general, relation = "within")[1]
      is_asia_hydrolakes <- terra::relate(study_area_vect, asia_extent_general, relation = "within")[1]
      is_eurasia_hydrolakes <- terra::relate(study_area_vect, eurasia_extent_general, relation = "within")[1]
      is_north_america_hydrolakes <- terra::relate(study_area_vect, northamerica_extent_general, relation = "within")[1]
      is_south_america_hydrolakes <- terra::relate(study_area_vect, southamerica_extent_general, relation = "within")[1]
      is_africa_hydrolakes <- terra::relate(study_area_vect, africa_extent_general, relation = "within")[1]
      is_australia_hydrolakes <- terra::relate(study_area_vect, australia_extent_general, relation = "within")[1]
      is_indo_pacific_hydrolakes <- terra::relate(study_area_vect, indopacific_extent_general, relation = "within")[1]
      is_holarctic_hydrolakes <- terra::relate(study_area_vect, holarctic_extent_general, relation = "within")[1]
      is_new_world_hydrolakes <- terra::relate(study_area_vect, newworld_extent_general, relation = "within")[1]
      is_old_world_hydrolakes <- terra::relate(study_area_vect, oldworld_extent_general, relation = "within")[1]
      landmask_dir <- file.path(rasters.dir, "landmask")
      if (!dir.exists(landmask_dir)) dir.create(landmask_dir, recursive = TRUE)
      hydrolakes_dir <- file.path(landmask_dir, "HydroLAKES")
      if (!dir.exists(hydrolakes_dir)) dir.create(hydrolakes_dir, recursive = TRUE)
      if (is_africa_hydrolakes) {
        if (verbose) message("Downloading Africa HydroLAKES file")
        hydrolakes_url <- "https://zenodo.org/records/17503891/files/HydroLAKES_Africa.zip?download=1"; min_size_mb <- 3; subname <- "HydroLAKES_Africa"
      } else if (is_asia_hydrolakes) {
        if (verbose) message("Downloading Asia HydroLAKES file")
        hydrolakes_url <- "https://zenodo.org/records/17503891/files/HydroLAKES_Asia.zip?download=1"; min_size_mb <- 55; subname <- "HydroLAKES_Asia"
      } else if (is_australia_hydrolakes) {
        if (verbose) message("Downloading Australia HydroLAKES file")
        hydrolakes_url <- "https://zenodo.org/records/17503891/files/HydroLAKES_Australia.zip?download=1"; min_size_mb <- 2; subname <- "HydroLAKES_Australia"
      } else if (is_europe_hydrolakes) {
        message("Downloading Europe HydroLAKES file")
        hydrolakes_url <- "https://zenodo.org/records/17503891/files/HydroLAKES_Europe.zip?download=1"; min_size_mb <- 30; subname <- "HydroLAKES_Europe"
      } else if (is_eurasia_hydrolakes) {
        if (verbose) message("Downloading Eurasia HydroLAKES file")
        hydrolakes_url <- "https://zenodo.org/records/17503891/files/HydroLAKES_Eurasia.zip?download=1"; min_size_mb <- 90; subname <- "HydroLAKES_Eurasia"
      } else if (is_north_america_hydrolakes) {
        if (verbose) message("Downloading North America HydroLAKES file")
        hydrolakes_url <- "https://zenodo.org/records/17503891/files/HydroLAKES_NorthAmerica.zip?download=1"; min_size_mb <- 600; subname <- "HydroLAKES_NorthAmerica"
      } else if (is_south_america_hydrolakes) {
        if (verbose) message("Downloading South America HydroLAKES file")
        hydrolakes_url <- "https://zenodo.org/records/17503891/files/HydroLAKES_SouthAmerica.zip?download=1"; min_size_mb <- 8; subname <- "HydroLAKES_SouthAmerica"
      } else if (is_indo_pacific_hydrolakes) {
        if (verbose) message("Downloading Indo Pacific HydroLAKES file")
        hydrolakes_url <- "https://zenodo.org/records/17503891/files/HydroLAKES_IndoPacific.zip?download=1"; min_size_mb <- 50; subname <- "HydroLAKES_IndoPacific"
      } else if (is_new_world_hydrolakes) {
        if (verbose) message("Downloading New World HydroLAKES file")
        hydrolakes_url <- "https://zenodo.org/records/17503891/files/HydroLAKES_NewWorld.zip?download=1"; min_size_mb <- 620; subname <- "HydroLAKES_NewWorld"
      } else if (is_old_world_hydrolakes) {
        if (verbose) message("Downloading Old World HydroLAKES file")
        hydrolakes_url <- "https://zenodo.org/records/17503891/files/HydroLAKES_OldWorld.zip?download=1"; min_size_mb <- 100; subname <- "HydroLAKES_OldWorld"
      } else if (is_holarctic_hydrolakes) {
        if (verbose) message("Downloading Holarctic HydroLAKES file")
        hydrolakes_url <- "https://zenodo.org/records/17503891/files/HydroLAKES_Holarctic.zip?download=1"; min_size_mb <- 730; subname <- "HydroLAKES_Holarctic"
      } else {
        if (verbose) message("Downloading global HydroLAKES file")
        hydrolakes_url <- "https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_polys_v10_shp.zip"; min_size_mb <- 400; subname <- "HydroLAKES_Global"
      }
      hydrolakes_zip <- file.path(hydrolakes_dir, paste0(subname, ".zip"))
      hydrolakes_unzip_dir <- file.path(hydrolakes_dir, subname)
      if (!dir.exists(hydrolakes_unzip_dir)) dir.create(hydrolakes_unzip_dir, recursive = TRUE)
      if (length(list.files(hydrolakes_unzip_dir, pattern = "\\.shp$", full.names = TRUE)) > 0 && !redownload.rasters) {
        if (verbose) message("HydroLAKES file already present - skipping download")
      } else {
        if (!file.exists(hydrolakes_zip) || file.size(hydrolakes_zip) < min_size_mb * 1e6 || redownload.rasters) {
          robust.download.raster(hydrolakes_url, hydrolakes_zip, min_size_mb = min_size_mb)
        }
        unzip(hydrolakes_zip, exdir = hydrolakes_unzip_dir, overwrite = TRUE)
        if (length(list.files(hydrolakes_unzip_dir, pattern = "\\.shp$", full.names = TRUE)) > 0) try(unlink(hydrolakes_zip, force = TRUE), silent = TRUE)
      }
      hydrolakes_path <- list.files(hydrolakes_unzip_dir, pattern = "\\.shp$", full.names = TRUE)[1]
      hydrolakes_vector_rds <- file.path(landmask_dir, "HydroLAKES_vector.rds")
      if (file.exists(hydrolakes_vector_rds) && file.info(hydrolakes_vector_rds)$size > 1e6 && !redownload.rasters) {
        hydrolakes_vector <- readRDS(hydrolakes_vector_rds) #use cached
      } else {
        hydrolakes_vector <- terra::vect(hydrolakes_path) #read shapefile
        saveRDS(hydrolakes_vector, hydrolakes_vector_rds) #cache for next run
      }
      hydrolakes_crop <- suppressWarnings(terra::crop(hydrolakes_vector, study_area_vect)) #spatial crop
      land_mask_erased_rds <- file.path(intermediate_files_dir, "land_mask_erased_hydrolakes.rds")
      if (file.exists(land_mask_erased_rds) && !overwrite) {
        if (verbose) message("Using already saved land mask with Hydrolakes erased")
        land_mask <- readRDS(land_mask_erased_rds)
      } else {
        if (verbose) message("Erasing Hydrolakes polygons from land mask ")
        land_mask_erased <- suppressWarnings(terra::erase(land_mask, hydrolakes_crop)) #remove lakes from land
        saveRDS(land_mask_erased, land_mask_erased_rds)
        land_mask <- land_mask_erased
      }
      invisible(gc())
      if (delete.intermediate.files.folders) unlink(hydrolakes_dir, recursive = TRUE, force = TRUE)
    }
  }

  # Delete landmask folders
  if (delete.intermediate.files.folders && exists("landmask_dir")) unlink(landmask_dir, recursive = TRUE, force = TRUE)

  # Generate background data
  background.data <- NULL
  if (generate.background.data) {
    invisible(gc())
    if (verbose) message("")
    if (verbose) message("-- Generating random background points within landmask --")
    background_file <- file.path(intermediate_files_dir, "background_intermediate_file.csv")
    if (file.exists(background_file) && file.info(background_file)$size > 0 && !overwrite && !redownload.rasters) {
      if (verbose) message("Background intermediate file already exists - skipping random spatial sampling and using saved file")
      background.data <- read.csv(background_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
    } else {
      landmask_sf <- sf::st_as_sf(land_mask)
      accessible_area_sf <- sf::st_as_sf(accessible_area)
      if (remove.hydrolakes.background) {
        hydrolakes_raster_resolution <- 0.002245
        hydrolakes_raster_resolution_km <- round(hydrolakes_raster_resolution * 111.32, 3) #convert degrees to km
        hydro_mask_file <- file.path(intermediate_files_dir, "HydroLakes_sampling_mask.tif")
        if (!file.exists(hydro_mask_file) || file.size(hydro_mask_file) < 1e5 || overwrite) {
          if (verbose) message("Rasterizing HydroLAKES polygons at ", hydrolakes_raster_resolution_km, " km (~", hydrolakes_raster_resolution, " degree) resolution")
          hydro_lakes_sf <- suppressMessages(sf::st_read(hydrolakes_path, quiet = TRUE))
          if (sf::st_crs(hydro_lakes_sf) != sf::st_crs(accessible_area_sf)) hydro_lakes_sf <- sf::st_transform(hydro_lakes_sf, sf::st_crs(accessible_area_sf))
          s2_state <- getOption("sf_use_s2", TRUE)
          options(sf_use_s2 = FALSE)
          hydro_crop <- suppressMessages(sf::st_crop(hydro_lakes_sf, sf::st_bbox(accessible_area_sf)))
          options(sf_use_s2 = s2_state)
          if (nrow(hydro_crop) == 0) stop("No HydroLAKES polygons overlap with the accessible area - check CRS or extent")
          extent_object <- terra::ext(sf::st_bbox(accessible_area_sf)[c("xmin", "xmax", "ymin", "ymax")])
          template_raster <- terra::rast(as(extent_object, "SpatExtent"), res = hydrolakes_raster_resolution, crs = "EPSG:4326")
          hydro_raster <- suppressWarnings(terra::rasterize(terra::vect(hydro_crop), template_raster, field = 1))
          template_raster[] <- 1
          mask_raster <- template_raster
          mask_raster[!is.na(hydro_raster)] <- NA
          terra::writeRaster(mask_raster, hydro_mask_file, overwrite = TRUE)
        } else {
          if (verbose) message("HydroLAKES raster already exists - loading from file")
          mask_raster <- terra::rast(hydro_mask_file)
        }
        if (verbose) message("Sampling background points from accessible area in batches")
        target_background_points <- N.background.points
        candidate_batch_size <- if (target_background_points < 50000) ceiling(target_background_points * 1.10) else 50000
        max_sampling_batches <- 40
        sampling_batch_counter <- 0
        background_points_accumulated <- data.frame()
        temporary_background_file <- file.path(intermediate_files_dir, "background_points_sampling_intermediate.csv")
        repeat {
          sampling_batch_counter <- sampling_batch_counter + 1
          if (verbose) message(paste("Sampling batch", sampling_batch_counter, "with", candidate_batch_size, "candidate points"))
          candidate_points_vect <- terra::spatSample(terra::vect(accessible_area_sf), size = candidate_batch_size, method = "random")
          candidate_mask_values <- terra::extract(mask_raster, candidate_points_vect)[, 2]
          valid_indices <- which(!is.na(candidate_mask_values))
          if (length(valid_indices) > 0) {
            valid_points_batch <- as.data.frame(terra::crds(candidate_points_vect[valid_indices, ]))
            background_points_accumulated <- rbind(background_points_accumulated, valid_points_batch)
          }
          if (nrow(background_points_accumulated) >= target_background_points) {
            background_points_accumulated <- background_points_accumulated[seq_len(target_background_points), ]
            colnames(background_points_accumulated) <- c("Longitude", "Latitude")
            rownames(background_points_accumulated) <- paste0("bg_", seq_len(nrow(background_points_accumulated)))
            write.csv(background_points_accumulated, temporary_background_file, row.names = TRUE)
            background.data <- read.csv(temporary_background_file, row.names = 1, check.names = FALSE)
            if (file.exists(temporary_background_file)) unlink(temporary_background_file, force = TRUE)
            break
          }
          if (sampling_batch_counter >= max_sampling_batches) {
            if (nrow(background_points_accumulated) == 0) stop(paste("Reached maximum batch limit (", max_sampling_batches, ") - no valid points obtained"))
            background_points_accumulated <- background_points_accumulated[seq_len(min(nrow(background_points_accumulated), target_background_points)), ]
            colnames(background_points_accumulated) <- c("Longitude", "Latitude")
            rownames(background_points_accumulated) <- paste0("bg_", seq_len(nrow(background_points_accumulated)))
            write.csv(background_points_accumulated, temporary_background_file, row.names = TRUE)
            background.data <- read.csv(temporary_background_file, row.names = 1, check.names = FALSE)
            if (file.exists(temporary_background_file)) unlink(temporary_background_file, force = TRUE)
            break
          }
        }
      } else {
        if (verbose) message("Sampling background points from accessible land area")
        accessible_land_sf <- try(suppressMessages(suppressWarnings(sf::st_intersection(landmask_sf, accessible_area_sf))), silent = TRUE)
        if (inherits(accessible_land_sf, "try-error")) {
          if (verbose) message("Invalid geometry detected - retrying intersection after st_make_valid function")
          landmask_sf <- suppressWarnings(sf::st_make_valid(landmask_sf))
          accessible_area_sf <- suppressWarnings(sf::st_make_valid(accessible_area_sf))
          accessible_land_sf <- try(suppressMessages(suppressWarnings(sf::st_intersection(landmask_sf, accessible_area_sf))), silent = TRUE)
        }
        if (inherits(accessible_land_sf, "try-error")) {
          if (verbose) message("Geometry still invalid after st_make_valid function - retrying intersection with planar geometry (s2 disabled)")
          s2_state <- getOption("sf_use_s2", TRUE)
          options(sf_use_s2 = FALSE)
          landmask_sf <- suppressWarnings(sf::st_make_valid(landmask_sf))
          accessible_area_sf <- suppressWarnings(sf::st_make_valid(accessible_area_sf))
          accessible_land_sf <- suppressMessages(suppressWarnings(sf::st_intersection(landmask_sf, accessible_area_sf)))
          options(sf_use_s2 = s2_state)
        }
        if (nrow(accessible_land_sf) == 0) stop("No intersection between landmask and accessible area - check coordinates or CRS")
        accessible_land_vect <- suppressWarnings(terra::vect(accessible_land_sf))
        background_points_vect <- terra::spatSample(accessible_land_vect, size = N.background.points, method = "random")
        if (!inherits(background_points_vect, "SpatVector")) background_points_vect <- terra::as.points(background_points_vect)
        background.data <- as.data.frame(terra::crds(background_points_vect))
        rownames(background.data) <- paste0("bg_", seq_len(nrow(background.data)))
        colnames(background.data) <- c("Longitude", "Latitude")
      }
      write.csv(background.data, background_file, row.names = TRUE)
    }
  }

  # Validate ClimateNA requirements
  if ("ClimateNA" %in% datasets_requested) {
    install_nichediv_dependencies(install_climatena = TRUE,
                                  install_whitebox = FALSE,
                                  install_data_table = TRUE,
                                  force = FALSE,
                                  verbose = verbose)
    system_info <- Sys.info()[["sysname"]]
    if (!grepl("Windows", system_info, ignore.case = TRUE)) stop("ClimateNAr package is only available for Windows!")
    if (!("ClimateNAr" %in% rownames(utils::installed.packages()))) stop("Package 'ClimateNAr' is required for ClimateNA extraction, but installation failed.")
    if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required for ClimateNA extraction, but installation failed")
  }

  # Create spatial vector of points
  coordinate_vector_env <- terra::vect(environmental_dataset[, c("Longitude", "Latitude")],
                                       geom = c("Longitude", "Latitude"),
                                       crs = CRS_all)
  coordinate_vector_bg <- if (!is.null(background.data)) {
    terra::vect(background.data[, c("Longitude", "Latitude")],
                geom = c("Longitude", "Latitude"),
                crs = CRS_all)
  } else NULL

  # Download and process elevation data (250m resolution: aggregated from Copernicus GLO-90 Digital Elevation Model: https://portal.opentopography.org/raster?opentopoID=OTSDEM.032021.4326.1; download from Zenodo mirror)
  if ("elevation" %in% datasets_requested) {
    invisible(gc())
    counter <- counter + 1
    if (verbose) message("")
    if (verbose) message(sprintf("-- Downloading and extracting elevation data (ca. 1.5h): env.dataset %d of %d --",
                                 counter, total_datasets))
    elevation_variable_name <- "Elevation"
    europe_extent_elevation <- terra::ext(-25, 41, 35, 71)
    is_europe_elevation <- terra::relate(study_area_vect, europe_extent_elevation, relation = "within")[1]
    asia_extent_elevation <- terra::ext(39, 177, -9, 80)
    is_asia_elevation <- terra::relate(study_area_vect, asia_extent_elevation, relation = "within")[1]
    eurasia_extent_elevation <- terra::ext(-25, 177, -9, 80)
    is_eurasia_elevation <- terra::relate(study_area_vect, eurasia_extent_elevation, relation = "within")[1]
    north_america_extent_elevation <- terra::ext(-170, -50, 10, 76)
    is_north_america_elevation <- terra::relate(study_area_vect, north_america_extent_elevation, relation = "within")[1]
    south_america_extent_elevation <- terra::ext(-90, -30, -59, 12)
    is_south_america_elevation <- terra::relate(study_area_vect, south_america_extent_elevation, relation = "within")[1]
    africa_extent_elevation <- terra::ext(-20, 51, -34, 37)
    is_africa_elevation <- terra::relate(study_area_vect, africa_extent_elevation, relation = "within")[1]
    australia_extent_elevation <- terra::ext(110, 176, -49, -7)
    is_australia_elevation <- terra::relate(study_area_vect, australia_extent_elevation, relation = "within")[1]
    indo_pacific_extent_elevation <- terra::ext(39, 177, -49, 80)
    is_indo_pacific_elevation <- terra::relate(study_area_vect, indo_pacific_extent_elevation, relation = "within")[1]
    holarctic_extent_elevation <- terra::ext(-170, 177, -9, 80)
    is_holarctic_elevation <- terra::relate(study_area_vect, holarctic_extent_elevation, relation = "within")[1]
    new_world_extent_elevation <- terra::ext(-170, -30, -59, 76)
    is_new_world_elevation <- terra::relate(study_area_vect, new_world_extent_elevation, relation = "within")[1]
    old_world_extent_elevation <- terra::ext(-25, 177, -34, 80)
    is_old_world_elevation <- terra::relate(study_area_vect, old_world_extent_elevation, relation = "within")[1]
    elevation_dir <- file.path(rasters.dir, "elevation")
    if (!dir.exists(elevation_dir)) dir.create(elevation_dir, recursive = TRUE)
    elevation_occ_csv_file <- file.path(intermediate_files_dir, paste0(elevation_variable_name, "_extracted_occurrence.csv"))
    elevation_bg_csv_file <- file.path(intermediate_files_dir, paste0(elevation_variable_name, "_extracted_background.csv"))
    elevation_occ_exists <- file.exists(elevation_occ_csv_file)
    elevation_bg_exists <- file.exists(elevation_bg_csv_file)
    if (elevation_occ_exists && file.info(elevation_occ_csv_file)$size > 0 && (!generate.background.data || (elevation_bg_exists && file.info(elevation_bg_csv_file)$size > 0)) && !overwrite) {
      if (verbose) message("Elevation data already exist - skipping download and extraction")
      extracted_occurrences <- read.csv(elevation_occ_csv_file, row.names = 1, check.names = FALSE)
      extracted_occurrences <- extracted_occurrences[base_ids, , drop = FALSE]
      environmental_dataset <- cbind(environmental_dataset, extracted_occurrences)
      if (generate.background.data && elevation_bg_exists) {
        extracted_background <- read.csv(elevation_bg_csv_file, row.names = 1, check.names = FALSE)
        if (!elevation_variable_name %in% names(extracted_background)) stop("Cached elevation background file is missing column: ", elevation_variable_name)
        if (!is.null(background.data)) background.data <- cbind(background.data, extracted_background[rownames(background.data), elevation_variable_name, drop = FALSE])
      }
    }
    if (!elevation_occ_exists || (generate.background.data && !elevation_bg_exists) || overwrite) {
      if (is_europe_elevation) {
        if (verbose) message("Using Europe elevation raster (250m resolution)")
        Europe_elevation_zip_url <- "https://zenodo.org/records/17487973/files/Copernicus_GLO90_Europe_250m.zip?download=1"
        elevation_raster <- download.and.load.elevation.tile("Europe", Europe_elevation_zip_url, elevation_dir, min_size_mb = 600, redownload = redownload.rasters, verbose = verbose)
      } else if (is_asia_elevation) {
        if (verbose) message("Using Asia elevation raster (250m resolution)")
        Asia_elevation_zip_url <- "https://zenodo.org/records/17450334/files/Copernicus_GLO90_Asia_250m.zip?download=1"
        elevation_raster <- download.and.load.elevation.tile("Asia", Asia_elevation_zip_url, elevation_dir, min_size_mb = 3200, redownload = redownload.rasters, verbose = verbose)
      } else if (is_north_america_elevation) {
        if (verbose) message("Using North America elevation raster (250m resolution)")
        NorthAmerica_elevation_zip_url <- "https://zenodo.org/records/17487973/files/Copernicus_GLO90_NorthAmerica_250m.zip?download=1"
        elevation_raster <- download.and.load.elevation.tile("NorthAmerica", NorthAmerica_elevation_zip_url, elevation_dir, min_size_mb = 1500, redownload = redownload.rasters, verbose = verbose)
      } else if (is_south_america_elevation) {
        if (verbose) message("Using South America elevation raster (250m resolution)")
        SouthAmerica_elevation_zip_url <- "https://zenodo.org/records/17450334/files/Copernicus_GLO90_SouthAmerica_250m.zip?download=1"
        elevation_raster <- download.and.load.elevation.tile("SouthAmerica", SouthAmerica_elevation_zip_url, elevation_dir, min_size_mb = 600, redownload = redownload.rasters, verbose = verbose)
      } else if (is_africa_elevation) {
        if (verbose) message("Using Africa elevation raster (250m resolution)")
        Africa_elevation_zip_url <- "https://zenodo.org/records/17458753/files/Copernicus_GLO90_Africa_250m.zip?download=1"
        elevation_raster <- download.and.load.elevation.tile("Africa", Africa_elevation_zip_url, elevation_dir, min_size_mb = 1500, redownload = redownload.rasters, verbose = verbose)
      } else if (is_australia_elevation) {
        if (verbose) message("Using Australia elevation raster (250m resolution)")
        Australia_elevation_zip_url <- "https://zenodo.org/records/17458753/files/Copernicus_GLO90_Australia_250m.zip?download=1"
        elevation_raster <- download.and.load.elevation.tile("Australia", Australia_elevation_zip_url, elevation_dir, min_size_mb = 400, redownload = redownload.rasters, verbose = verbose)
      } else if (is_new_world_elevation) {
        if (verbose) message("Using New World elevation raster (250m resolution)")
        NewWorld_elevation_zip_url <- "https://zenodo.org/records/17485342/files/Copernicus_GLO90_NewWorld_250m.zip?download=1"
        elevation_raster <- download.and.load.elevation.tile("NewWorld", NewWorld_elevation_zip_url, elevation_dir, min_size_mb = 2500, redownload = redownload.rasters, verbose = verbose)
      } else if (is_indo_pacific_elevation) {
        if (verbose) message("Using Indo-Pacific elevation raster (250m resolution)")
        IndoPacific_elevation_zip_url <- "https://zenodo.org/records/17485342/files/Copernicus_GLO90_IndoPacific_250m.zip?download=1"
        elevation_raster <- download.and.load.elevation.tile("IndoPacific", IndoPacific_elevation_zip_url, elevation_dir, min_size_mb = 3800, redownload = redownload.rasters, verbose = verbose)
      } else if (is_eurasia_elevation) {
        if (verbose) message("Using Eurasia elevation raster (250m resolution)")
        Eurasia_elevation_zip_url <- "https://zenodo.org/records/17485342/files/Copernicus_GLO90_Eurasia_250m.zip?download=1"
        elevation_raster <- download.and.load.elevation.tile("Eurasia", Eurasia_elevation_zip_url, elevation_dir, min_size_mb = 4000, redownload = redownload.rasters, verbose = verbose)
      } else if (is_holarctic_elevation) {
        if (verbose) message("Using Holarctic elevation raster (250m resolution)")
        Holarctic_elevation_zip_url <- "https://zenodo.org/records/17485342/files/Copernicus_GLO90_Holarctic_250m.zip?download=1"
        elevation_raster <- download.and.load.elevation.tile("Holarctic", Holarctic_elevation_zip_url, elevation_dir, min_size_mb = 5900, redownload = redownload.rasters, verbose = verbose)
      } else {
        if (verbose) message("Using Global elevation raster (250m resolution)")
        Global_elevation_zip_url <- "https://zenodo.org/records/17485342/files/Copernicus_GLO90_Global_250m.zip?download=1"
        elevation_raster <- download.and.load.elevation.tile("Global", Global_elevation_zip_url, elevation_dir, min_size_mb = 8900, redownload = redownload.rasters, verbose = verbose)
      }
      names(elevation_raster) <- elevation_variable_name
      raster_crs <- terra::crs(elevation_raster)
      coord_env <- coordinate_vector_env
      coord_bg <- coordinate_vector_bg
      if (!terra::same.crs(coordinate_vector_env, elevation_raster)) {
        if (verbose) message("Projecting coordinates to match elevation raster CRS")
        coord_env <- terra::project(coordinate_vector_env, raster_crs)
        if (!is.null(coordinate_vector_bg)) coord_bg <- terra::project(coordinate_vector_bg, raster_crs)
      }
      result <- extract.and.cache.env.dataset(dataset_name = elevation_variable_name,
                                              raster_object = elevation_raster,
                                              coord_env = coord_env,
                                              coord_bg = coord_bg,
                                              environmental_dataset = environmental_dataset,
                                              background.data = background.data,
                                              output.dir = intermediate_files_dir,
                                              overwrite = overwrite,
                                              generate.background.data = generate.background.data,
                                              verbose = verbose)
      environmental_dataset <- result$environmental_dataset
      background.data <- result$background.data
      if (any(!elevation_variable_name %in% names(environmental_dataset))) stop("Elevation column missing from occurrence data after extraction - extraction seems to not have worked")
      if (generate.background.data && any(!elevation_variable_name %in% names(background.data))) stop("Elevation column missing from background data after extraction - extraction seems to not have worked")
      if (all(is.na(environmental_dataset$Elevation))) warning("All Elevation values in occurrence data are NA - check CRS or study extent")
      if (generate.background.data && all(is.na(background.data$Elevation))) warning("All Elevation values in background data are NA - check CRS or study extent")
      if (exists("elevation_raster")) rm(elevation_raster)
      if (exists("result")) rm(result)
      invisible(gc())
      if (delete.intermediate.files.folders && ("elevation" %in% datasets_requested) && !("terrain" %in% datasets_requested) && !("ClimateNA" %in% datasets_requested)) unlink(file.path(rasters.dir, "elevation"), recursive = TRUE, force = TRUE)
    }
  }

  # Download and process ClimateNA data (only for North America; https://climatena.ca/; scale-free; using ClimateNAr package)
  if ("ClimateNA" %in% datasets_requested) {
    invisible(gc())
    counter <- counter + 1
    if (verbose) message("")
    if (verbose) message(sprintf("-- Downloading and extracting ClimateNA data (ca. 3 min): env.dataset %d of %d --",
                                 counter, total_datasets))
    ClimateNAr.original <- eval(parse(text = "ClimateNAr::ClimateNAr"), envir = parent.frame())
    ClimateNAr_source_lines <- deparse(ClimateNAr.original, width.cutoff = 500)
    ClimateNAr_clean_lines <- grep("(message\\s*\\()|(print\\s*\\(.*Completed for)", ClimateNAr_source_lines, value = TRUE, invert = TRUE)
    ClimateNAr_clean_lines <- gsub("^\\s*\\)\\s*$", "", ClimateNAr_clean_lines)
    ClimateNAr_clean_lines <- gsub("^\\s*}\\s*$", "}", ClimateNAr_clean_lines)
    ClimateNAr_clean_lines <- gsub("if\\s*\\([^)]*\\)\\s*\\{\\s*\\}", "", ClimateNAr_clean_lines)
    ClimateNAr_clean_lines <- gsub("else\\s*\\{\\s*\\}", "", ClimateNAr_clean_lines)
    ClimateNAr_clean_lines <- gsub("^(\\s*stk_el\\s*<-\\s*terra::rast\\(.*\\))$",
                                   "\\1\n        terra::crs(stk_el) <- \"EPSG:4326\"",
                                   ClimateNAr_clean_lines)
    ClimateNAr_clean_lines <- gsub("^(\\s*stk_tmx\\s*<-\\s*terra::rast\\(.*\\))$",
                                   "\\1\n        terra::crs(stk_tmx) <- \"EPSG:4326\"",
                                   ClimateNAr_clean_lines)
    ClimateNAr_clean_lines <- gsub("^(\\s*stk_tmn\\s*<-\\s*terra::rast\\(.*\\))$",
                                   "\\1\n        terra::crs(stk_tmn) <- \"EPSG:4326\"",
                                   ClimateNAr_clean_lines)
    ClimateNAr_clean_lines <- gsub("^(\\s*stk_ppt\\s*<-\\s*terra::rast\\(.*\\))$",
                                   "\\1\n        terra::crs(stk_ppt) <- \"EPSG:4326\"",
                                   ClimateNAr_clean_lines)
    ClimateNAr_clean_lines <- gsub("^(\\s*stk_tmxr\\s*<-\\s*terra::rast\\(.*\\))$",
                                   "\\1\n        terra::crs(stk_tmxr) <- \"EPSG:4326\"",
                                   ClimateNAr_clean_lines)
    ClimateNAr_clean_lines <- gsub("^(\\s*stk_tmnr\\s*<-\\s*terra::rast\\(.*\\))$",
                                   "\\1\n        terra::crs(stk_tmnr) <- \"EPSG:4326\"",
                                   ClimateNAr_clean_lines)
    ClimateNAr_clean_lines <- gsub("^(\\s*stk_pptr\\s*<-\\s*terra::rast\\(.*\\))$",
                                   "\\1\n        terra::crs(stk_pptr) <- \"EPSG:4326\"",
                                   ClimateNAr_clean_lines)
    ClimateNAr_clean_lines <- gsub("^(\\s*stk_h\\s*<-\\s*terra::rast\\(.*\\))$",
                                   "\\1\n                  terra::crs(stk_h) <- \"EPSG:4326\"",
                                   ClimateNAr_clean_lines)
    ClimateNAr_clean_lines <- gsub("^(\\s*stk_h_py\\s*<-\\s*terra::rast\\(.*\\))$",
                                   "\\1\n                    terra::crs(stk_h_py) <- \"EPSG:4326\"",
                                   ClimateNAr_clean_lines)
    ClimateNAr_clean_lines <- gsub("^(\\s*cruNrm_stk\\s*<-\\s*terra::rast\\(.*\\))$",
                                   "\\1\n                    terra::crs(cruNrm_stk) <- \"EPSG:4326\"",
                                   ClimateNAr_clean_lines)
    ClimateNAr_clean_lines <- gsub("^(\\s*stk_gcm\\s*<-\\s*terra::rast\\(.*\\))$",
                                   "\\1\n                  terra::crs(stk_gcm) <- \"EPSG:4326\"",
                                   ClimateNAr_clean_lines)
    ClimateNAr_clean_src <- paste(ClimateNAr_clean_lines, collapse = "\n")
    ClimateNAr_clean_src <- gsub("(^|[^:[:alnum:]_.])as\\.data\\.table\\s*\\(", "\\1data.table::as.data.table(", ClimateNAr_clean_src, perl = TRUE)
    ClimateNAr_clean_src <- gsub("(^|[^:[:alnum:]_.])fread\\s*\\(", "\\1data.table::fread(", ClimateNAr_clean_src, perl = TRUE)
    ClimateNAr_env <- new.env(parent = asNamespace("ClimateNAr"))
    ClimateNAr_env$as.data.table <- data.table::as.data.table
    ClimateNAr_env$fread <- data.table::fread
    ClimateNAr.replicate <- eval(parse(text = ClimateNAr_clean_src), envir = ClimateNAr_env)
    ClimateNA_dir <- file.path(rasters.dir, "ClimateNA")
    if (!dir.exists(ClimateNA_dir)) dir.create(ClimateNA_dir, recursive = TRUE)
    ClimateNA_dir_abs <- normalizePath(ClimateNA_dir, winslash = "\\", mustWork = FALSE)
    climatena_occ_file <- file.path(intermediate_files_dir, "ClimateNA_extracted_occurrence.csv")
    climatena_bg_file <- file.path(intermediate_files_dir, "ClimateNA_extracted_background.csv")
    climatena_variable_names <- c("MAT", "MWMT", "MCMT", "TD", "MAP", "MSP",
                                  "AHM", "SHM", "bFFP", "eFFP", "FFP", "PAS",
                                  "NFFD", "CMD", "CMI", "Eref", "DD_0", "DD5",
                                  "DD18", "DD_18", "DD1040", "RH",
                                  paste0("Tmax", sprintf("%02d", 1:12)),
                                  paste0("Tmin", sprintf("%02d", 1:12)),
                                  paste0("Tave", sprintf("%02d", 1:12)),
                                  paste0("PPT", sprintf("%02d", 1:12)),
                                  paste0("Eref", sprintf("%02d", 1:12)))
    if (file.exists(climatena_occ_file) && (!generate.background.data || file.exists(climatena_bg_file)) && !overwrite) {
      if (verbose) message("ClimateNA data already exist - skipping download and extraction - loading all variables from intermediate files")
      occ_climatena <- read.csv(climatena_occ_file, row.names = 1, check.names = FALSE)
      occ_climatena <- occ_climatena[base_ids, , drop = FALSE]
      environmental_dataset <- cbind(environmental_dataset, occ_climatena)
      if (generate.background.data && file.exists(climatena_bg_file)) {
        bg_climatena <- read.csv(climatena_bg_file, row.names = 1, check.names = FALSE)
        if (!is.null(background.data)) background.data <- cbind(background.data, bg_climatena[rownames(background.data), , drop = FALSE])
      }
      if (all(is.na(environmental_dataset[, setdiff(names(occ_climatena), "ID")]))) warning("All ClimateNA values in occurrence data are NA - check CRS or extent")
      if (generate.background.data && all(is.na(background.data[, setdiff(names(bg_climatena), "ID")]))) warning("All ClimateNA values in background data are NA - check CRS or extent")
    } else {
      in_csv <- file.path(ClimateNA_dir_abs, "Lat_Long_Elev.csv")
      write.csv(data.frame(ID1 = base_ids,
                           ID2 = seq_len(nrow(environmental_dataset)),
                           lat = as.numeric(environmental_dataset$Latitude),
                           long = as.numeric(environmental_dataset$Longitude),
                           el = as.numeric(environmental_dataset$Elevation)),
                in_csv,
                row.names = FALSE)
      pkg_path <- find.package("ClimateNAr")
      oldwd <- getwd()
      setwd(pkg_path)
      ClimateNAr.replicate(inputFile = normalizePath(in_csv, winslash = "\\"),
                           varList = "YM",
                           periodList = "Decade_2011_2020.dcd",
                           outDir = "")
      setwd(oldwd)
      occurrences_output_candidates <- c(file.path(ClimateNA_dir_abs, "ClimateNALat_Long_Elev_Decade_2011_2020.csv"),
                                         file.path(ClimateNA_dir_abs, "Lat_Long_Elev_Decade_2011_2020.csv"))
      occurrences_output <- occurrences_output_candidates[file.exists(occurrences_output_candidates)]
      if (!length(occurrences_output)) {
        all_csvs <- list.files(ClimateNA_dir_abs, pattern = "\\.csv$", full.names = TRUE)
        all_csvs <- setdiff(all_csvs, in_csv)
        if (length(all_csvs)) occurrences_output <- all_csvs[which.max(file.info(all_csvs)$mtime)]
      } else occurrences_output <- occurrences_output[1]
      if (!length(occurrences_output)) stop("ClimateNA extraction failed for occurrence data")
      occurrences_climatena <- data.table::fread(occurrences_output)
      occurrences_climatena <- dplyr::select(occurrences_climatena, -ID2, -lat, -long, -el)
      occurrences_climatena <- dplyr::rename(occurrences_climatena, ID = ID1)
      rownames(occurrences_climatena) <- occurrences_climatena$ID
      occurrences_climatena$ID <- NULL
      occurrences_climatena <- as.data.frame(occurrences_climatena)
      occurrences_climatena$ID <- as.character(rownames(occurrences_climatena))
      rownames(occurrences_climatena) <- occurrences_climatena$ID
      occurrences_climatena$ID <- NULL
      if (!any(base_ids %in% rownames(occurrences_climatena))) {
        occurrences_climatena <- occurrences_climatena[seq_len(nrow(environmental_dataset)), , drop = FALSE]
        rownames(occurrences_climatena) <- base_ids
      }
      environmental_dataset <- cbind(environmental_dataset, occurrences_climatena[base_ids, , drop = FALSE])
      write.csv(occurrences_climatena, climatena_occ_file)
      if (!is.null(background.data)) {
        extra_point <- NULL
        if (nrow(background.data) %% 50000 == 0) {
          extra_point <- background.data[nrow(background.data), , drop = FALSE]
          background.data <- background.data[-nrow(background.data), ]
        }
        background_intermediate_temp <- file.path(intermediate_files_dir, sprintf("background_intermediate_file_run%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")))
        write.csv(background.data, background_intermediate_temp, row.names = FALSE)
        bg_csv <- file.path(ClimateNA_dir_abs, "Lat_Long_Elev_background.csv")
        write.csv(data.frame(ID1 = seq_len(nrow(background.data)),
                             ID2 = seq_len(nrow(background.data)),
                             lat = as.numeric(background.data$Latitude),
                             long = as.numeric(background.data$Longitude),
                             el = as.numeric(background.data$Elevation)),
                  bg_csv,
                  row.names = FALSE)
        pkg_path <- find.package("ClimateNAr")
        oldwd <- getwd()
        setwd(pkg_path)
        ClimateNAr.replicate(inputFile = normalizePath(bg_csv, winslash = "\\"),
                             varList = "YM",
                             periodList = "Decade_2011_2020.dcd",
                             outDir = "")
        setwd(oldwd)
        bg_output_candidates <- c(file.path(ClimateNA_dir_abs, "ClimateNALat_Long_Elev_background_Decade_2011_2020.csv"),
                                  file.path(ClimateNA_dir_abs, "Lat_Long_Elev_background_Decade_2011_2020.csv"))
        bg_output <- bg_output_candidates[file.exists(bg_output_candidates)]
        if (!length(bg_output)) {
          all_csvs <- list.files(ClimateNA_dir_abs, pattern = "\\.csv$", full.names = TRUE)
          all_csvs <- setdiff(all_csvs, bg_csv)
          if (length(all_csvs)) bg_output <- all_csvs[which.max(file.info(all_csvs)$mtime)]
        } else bg_output <- bg_output[1]
        if (!length(bg_output)) stop("ClimateNA extraction failed for background data")
        background_climatena <- data.table::fread(bg_output)
        background_climatena <- dplyr::select(background_climatena, -ID2, -lat, -long, -el)
        background_climatena <- dplyr::rename(background_climatena, ID = ID1)
        background_climatena$ID <- NULL
        background_climatena <- as.data.frame(background_climatena)
        if (!is.null(extra_point)) {
          two_csv <- file.path(ClimateNA_dir_abs, "Lat_Long_Elev_two_points.csv")
          extra_points_df <- rbind(extra_point, background.data[1, , drop = FALSE])
          write.csv(data.frame(ID1 = 1:2,
                               ID2 = 1:2,
                               lat = as.numeric(extra_points_df$Latitude),
                               long = as.numeric(extra_points_df$Longitude),
                               el = as.numeric(extra_points_df$Elevation)),
                    two_csv,
                    row.names = FALSE)
          pkg_path <- find.package("ClimateNAr")
          oldwd <- getwd()
          setwd(pkg_path)
          ClimateNAr.replicate(inputFile = normalizePath(two_csv, winslash = "\\"),
                               varList = "YM",
                               periodList = "Decade_2011_2020.dcd",
                               outDir = "")
          setwd(oldwd)
          two_output_candidates <- c(file.path(ClimateNA_dir_abs, "ClimateNALat_Long_Elev_two_points_Decade_2011_2020.csv"),
                                     file.path(ClimateNA_dir_abs, "Lat_Long_Elev_two_points_Decade_2011_2020.csv"))
          two_output <- two_output_candidates[file.exists(two_output_candidates)]
          if (!length(two_output)) {
            all_csvs <- list.files(ClimateNA_dir_abs, pattern = "\\.csv$", full.names = TRUE)
            all_csvs <- setdiff(all_csvs, two_csv)
            if (length(all_csvs)) two_output <- all_csvs[which.max(file.info(all_csvs)$mtime)]
          } else two_output <- two_output[1]
          if (!length(two_output)) stop("ClimateNA extraction failed for two background points")
          two_climatena <- data.table::fread(two_output)
          two_climatena <- dplyr::select(two_climatena, -ID2, -lat, -long, -el)
          two_climatena <- dplyr::rename(two_climatena, ID = ID1)
          two_climatena$ID <- NULL
          two_climatena <- as.data.frame(two_climatena)
          single_climatena <- two_climatena[1, , drop = FALSE]
          missing_from_single <- setdiff(names(background_climatena), names(single_climatena))
          if (length(missing_from_single) > 0) {
            for (missing_name in missing_from_single) single_climatena[[missing_name]] <- NA
          }
          missing_from_background <- setdiff(names(single_climatena), names(background_climatena))
          if (length(missing_from_background) > 0) {
            for (missing_name in missing_from_background) background_climatena[[missing_name]] <- NA
          }
          background_climatena <- background_climatena[, union(names(background_climatena), names(single_climatena)), drop = FALSE]
          single_climatena <- single_climatena[, union(names(background_climatena), names(single_climatena)), drop = FALSE]
          background_climatena <- rbind(background_climatena, single_climatena)
          background.data <- rbind(background.data, extra_point)
        }
        rownames(background_climatena) <- rownames(background.data)
        background.data <- cbind(background.data, background_climatena[rownames(background.data), , drop = FALSE])
        write.csv(background_climatena, climatena_bg_file)
        dest_file <- file.path(intermediate_files_dir, "ClimateNA_extracted_background.csv")
        if (file.exists(climatena_bg_file) &&
            normalizePath(climatena_bg_file) != normalizePath(dest_file)) {
          file.copy(climatena_bg_file, dest_file, overwrite = TRUE)
          file.remove(climatena_bg_file)
        }
      }
      if (any(!climatena_variable_names %in% names(environmental_dataset))) stop("ClimateNA columns missing from occurrence data after extraction - extraction seems to not have worked")
      if (generate.background.data && any(!climatena_variable_names %in% names(background.data))) stop("ClimateNA columns missing from background data after extraction - extraction seems to not have worked")
      suppressWarnings(try(terra::tmpFiles(remove = TRUE), silent = TRUE))
      if (exists("occurrences_climatena")) rm(occurrences_climatena)
      if (exists("background_climatena")) rm(background_climatena)
      if (exists("in_csv")) rm(in_csv)
      if (exists("bg_csv")) rm(bg_csv)
      if (exists("occurrences_output")) rm(occurrences_output)
      if (exists("bg_output")) rm(bg_output)
      invisible(gc())
      if (delete.intermediate.files.folders && "elevation" %in% datasets_requested && !"terrain" %in% datasets_requested && "ClimateNA" %in% datasets_requested)
        unlink(file.path(rasters.dir, "elevation"), recursive = TRUE, force = TRUE)
      if (delete.intermediate.files.folders && "ClimateNA" %in% datasets_requested)
        unlink(file.path(rasters.dir, "ClimateNA"), recursive = TRUE, force = TRUE)
    }
  }

  # Download and process Enhanced Vegetation Index (EVI) data (https://zenodo.org/records/17449851/files/EVI_North_America_250m.zip; 250m resolution)
  if ("EVI" %in% datasets_requested) {
    invisible(gc())
    counter <- counter + 1
    if (verbose) message("")
    if (verbose) message(sprintf("-- Downloading and extracting EVI data (ca. 1.5h): env.dataset %d of %d --",
                                 counter, total_datasets))
    evi_variable_names <- c("EVI_median", "EVI_min", "EVI_max", "EVI_1", "EVI_2", "EVI_3", "EVI_4", "EVI_5", "EVI_6")
    evi_dir <- file.path(rasters.dir, "EVI")
    if (!dir.exists(evi_dir)) dir.create(evi_dir, recursive = TRUE)
    evi_occ_csv_file <- file.path(intermediate_files_dir, "EVI_extracted_occurrence.csv")
    evi_bg_csv_file <- file.path(intermediate_files_dir, "EVI_extracted_background.csv")
    evi_occ_exists <- file.exists(evi_occ_csv_file)
    evi_bg_exists <- file.exists(evi_bg_csv_file)
    if (evi_occ_exists && (!generate.background.data || evi_bg_exists) && !overwrite) {
      if (verbose) message("EVI data already exist - skipping download and extraction - loading from intermediate files")
      occ_evi <- read.csv(evi_occ_csv_file, row.names = 1, check.names = FALSE)
      occ_evi <- occ_evi[base_ids, , drop = FALSE]
      environmental_dataset <- cbind(environmental_dataset, occ_evi)
      if (generate.background.data && evi_bg_exists) {
        bg_evi <- read.csv(evi_bg_csv_file, row.names = 1, check.names = FALSE)
        if (!is.null(background.data)) background.data <- cbind(background.data, bg_evi[rownames(background.data), , drop = FALSE])
      }
      if (all(is.na(environmental_dataset[, setdiff(names(occ_evi), "ID")]))) warning("All EVI values in occurrence data are NA - check CRS or extent")
      if (generate.background.data && all(is.na(background.data[, setdiff(names(bg_evi), "ID")]))) warning("All EVI values in background data are NA - check CRS or extent")
    } else {
      is.study.area.within.extent <- function(study_area_vect, extent_object) {
        extent_union <- terra::aggregate(do.call(rbind, lapply(extent_object, terra::as.polygons)))
        terra::relate(study_area_vect, extent_union, relation = "within")[1]
      }
      europe_extent_evi <- list(terra::ext(-26, 65, 34, 83))
      asia_extent_evi <- list(terra::ext(21, 180, -13, 84), terra::ext(-180, -168, 50, 75))
      northamerica_extent_evi <- list(terra::ext(-180, -10, 5, 85))
      southamerica_extent_evi <- list(terra::ext(-92, -30, -56, 15))
      africa_extent_evi <- list(terra::ext(-26, 57, -36, 39))
      australia_extent_evi <- list(terra::ext(90, 180, -56, 25), terra::ext(-180, -170, -56, 25))
      indopacific_extent_evi <- c(asia_extent_evi, australia_extent_evi)
      eurasia_extent_evi <- c(europe_extent_evi, asia_extent_evi)
      holarctic_extent_evi <- c(europe_extent_evi, asia_extent_evi, northamerica_extent_evi)
      newworld_extent_evi <- c(northamerica_extent_evi, southamerica_extent_evi)
      oldworld_extent_evi <- c(europe_extent_evi, africa_extent_evi, asia_extent_evi, australia_extent_evi)
      is_europe_evi <- is.study.area.within.extent(study_area_vect, europe_extent_evi)
      is_asia_evi <- is.study.area.within.extent(study_area_vect, asia_extent_evi)
      is_north_america_evi <- is.study.area.within.extent(study_area_vect, northamerica_extent_evi)
      is_south_america_evi <- is.study.area.within.extent(study_area_vect, southamerica_extent_evi)
      is_africa_evi <- is.study.area.within.extent(study_area_vect, africa_extent_evi)
      is_australia_evi <- is.study.area.within.extent(study_area_vect, australia_extent_evi)
      is_indo_pacific_evi <- is.study.area.within.extent(study_area_vect, indopacific_extent_evi)
      is_eurasia_evi <- is.study.area.within.extent(study_area_vect, eurasia_extent_evi)
      is_holarctic_evi <- is.study.area.within.extent(study_area_vect, holarctic_extent_evi)
      is_new_world_evi <- is.study.area.within.extent(study_area_vect, newworld_extent_evi)
      is_old_world_evi <- is.study.area.within.extent(study_area_vect, oldworld_extent_evi)
      evi_region <- NULL
      evi_region_label <- NULL
      evi_min_size_mb <- NULL
      if (is_europe_evi) {
        evi_region <- "europe"
        evi_region_label <- "Europe"
        evi_min_size_mb <- 400
      } else if (is_asia_evi) {
        evi_region <- "asia"
        evi_region_label <- "Asia"
        evi_min_size_mb <- 1500
      } else if (is_north_america_evi) {
        evi_region <- "northamerica"
        evi_region_label <- "North America"
        evi_min_size_mb <- 700
      } else if (is_south_america_evi) {
        evi_region <- "southamerica"
        evi_region_label <- "South America"
        evi_min_size_mb <- 400
      } else if (is_africa_evi) {
        evi_region <- "africa"
        evi_region_label <- "Africa"
        evi_min_size_mb <- 800
      } else if (is_australia_evi) {
        evi_region <- "australia"
        evi_region_label <- "Australia"
        evi_min_size_mb <- 350
      } else if (is_indo_pacific_evi) {
        evi_region <- "indopacific"
        evi_region_label <- "Indo-Pacific"
        evi_min_size_mb <- 1800
      } else if (is_eurasia_evi) {
        evi_region <- "eurasia"
        evi_region_label <- "Eurasia"
        evi_min_size_mb <- 1800
      } else if (is_holarctic_evi) {
        evi_region <- "holarctic"
        evi_region_label <- "Holarctic"
        evi_min_size_mb <- 2500
      } else if (is_new_world_evi) {
        evi_region <- "newworld"
        evi_region_label <- "New World"
        evi_min_size_mb <- 1100
      } else if (is_old_world_evi) {
        evi_region <- "oldworld"
        evi_region_label <- "Old World"
        evi_min_size_mb <- 2500
      }
      if (!is.null(evi_region)) {
        if (verbose) message("Using ", evi_region_label, " EVI rasters (250m resolution)")
        evi_record_ids <- c("18077970", "19582709", "18841516")
        EVI_files <- file.path(evi_dir, paste0("EVI_", 1:6, "_", evi_region, ".tif"))
        existing_files <- file.exists(EVI_files)
        valid_files <- existing_files & (file.size(EVI_files) > evi_min_size_mb * 1e6)
        if (all(valid_files) && !redownload.rasters) {
          if (verbose) message("EVI rasters already present - skipping download")
        } else {
          missing_indices <- which(!valid_files | redownload.rasters)
          if (verbose) message("Downloading ", length(missing_indices), " of 6 missing EVI rasters")
          for (i in seq_along(missing_indices)) {
            idx <- missing_indices[i]
            success <- FALSE
            for (record_id in evi_record_ids) {
              EVI_url <- paste0("https://zenodo.org/records/", record_id, "/files/EVI_", idx, "_", evi_region, ".tif?download=1")
              try_result <- try(robust.download.raster(EVI_url, EVI_files[idx], max_attempts = 1, min_size_mb = evi_min_size_mb), silent = TRUE)
              if (!inherits(try_result, "try-error") && file.exists(EVI_files[idx]) && file.size(EVI_files[idx]) > evi_min_size_mb * 1e6) {
                success <- TRUE
                break
              }
            }
            if (!success) stop("Failed to download EVI raster for region '", evi_region, "' and layer EVI_", idx, " from all provided Zenodo records")
          }
        }
      } else {
        if (verbose) message("Using global EVI raster (250m resolution)")
        EVI_urls <- c(
          "https://s3.openlandmap.org/arco/evi_mod13q1.tmwm.inpaint_p.90_250m_s_20200101_20200228_go_epsg.4326_v20230608.tif",
          "https://s3.openlandmap.org/arco/evi_mod13q1.tmwm.inpaint_p.90_250m_s_20200301_20200430_go_epsg.4326_v20230608.tif",
          "https://s3.openlandmap.org/arco/evi_mod13q1.tmwm.inpaint_p.90_250m_s_20200501_20200630_go_epsg.4326_v20230608.tif",
          "https://s3.openlandmap.org/arco/evi_mod13q1.tmwm.inpaint_p.90_250m_s_20200701_20200831_go_epsg.4326_v20230608.tif",
          "https://s3.openlandmap.org/arco/evi_mod13q1.tmwm.inpaint_p.90_250m_s_20200901_20201031_go_epsg.4326_v20230608.tif",
          "https://s3.openlandmap.org/arco/evi_mod13q1.tmwm.inpaint_p.90_250m_s_20201101_20201231_go_epsg.4326_v20230608.tif"
        )
        EVI_files <- file.path(evi_dir, basename(EVI_urls))
        existing_files <- file.exists(EVI_files)
        valid_files <- existing_files & (file.size(EVI_files) > 4.6e9)
        if (all(valid_files) && !redownload.rasters) {
          if (verbose) message("EVI rasters already present - skipping download")
        } else {
          missing_indices <- which(!valid_files | redownload.rasters)
          if (verbose) message("Downloading ", length(missing_indices), " of six missing EVI rasters")
          for (i in seq_along(missing_indices)) {
            idx <- missing_indices[i]
            robust.download.raster(EVI_urls[idx], EVI_files[idx])
          }
        }
      }
      EVI_stack <- terra::rast(EVI_files)
      names(EVI_stack) <- paste0("EVI_", seq_len(terra::nlyr(EVI_stack)))
      raster_crs <- terra::crs(EVI_stack)
      coord_env <- coordinate_vector_env
      coord_bg <- coordinate_vector_bg
      if (!terra::same.crs(coordinate_vector_env, EVI_stack)) {
        if (verbose) message("Projecting coordinates to match EVI raster CRS")
        coord_env <- terra::project(coordinate_vector_env, raster_crs)
        if (!is.null(coordinate_vector_bg)) coord_bg <- terra::project(coordinate_vector_bg, raster_crs)
      }
      extracted_occurrences <- as.matrix(terra::extract(EVI_stack, coord_env, ID = FALSE))
      colnames(extracted_occurrences) <- names(EVI_stack)
      environmental_dataset <- cbind(environmental_dataset, extracted_occurrences)
      environmental_dataset$EVI_median <- apply(extracted_occurrences, 1, function(x) ifelse(all(is.na(x)), NA, median(x, na.rm = TRUE)))
      environmental_dataset$EVI_min <- apply(extracted_occurrences, 1, function(x) ifelse(all(is.na(x)), NA, min(x, na.rm = TRUE)))
      environmental_dataset$EVI_max <- apply(extracted_occurrences, 1, function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE)))
      write.csv(environmental_dataset[, c(names(EVI_stack), "EVI_median", "EVI_min", "EVI_max")], evi_occ_csv_file)
      if (!is.null(background.data) && generate.background.data) {
        extracted_background <- as.matrix(terra::extract(EVI_stack, coord_bg, ID = FALSE))
        colnames(extracted_background) <- names(EVI_stack)
        background.data <- cbind(background.data, extracted_background)
        background.data$EVI_median <- apply(extracted_background, 1, function(x) ifelse(all(is.na(x)), NA, median(x, na.rm = TRUE)))
        background.data$EVI_min <- apply(extracted_background, 1, function(x) ifelse(all(is.na(x)), NA, min(x, na.rm = TRUE)))
        background.data$EVI_max <- apply(extracted_background, 1, function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE)))
        write.csv(background.data[, c(names(EVI_stack), "EVI_median", "EVI_min", "EVI_max")], evi_bg_csv_file)
      }
      if (any(!evi_variable_names %in% names(environmental_dataset))) stop("EVI columns missing from occurrence data after extraction - extraction seems to not have worked")
      if (generate.background.data && any(!evi_variable_names %in% names(background.data))) stop("EVI columns missing from background data after extraction - extraction seems to not have worked")
      suppressWarnings(try(terra::tmpFiles(remove = TRUE), silent = TRUE))
      if (exists("EVI_stack")) rm(EVI_stack)
      if (exists("extracted_occurrences")) rm(extracted_occurrences)
      if (exists("extracted_background")) rm(extracted_background)
      if (exists("EVI_files")) rm(EVI_files)
      if (exists("EVI_zip_file")) rm(EVI_zip_file)
      invisible(gc())
      if (all(is.na(environmental_dataset$EVI_median))) warning("All EVI values in occurrence data are NA - check CRS or study extent")
      if (generate.background.data && all(is.na(background.data$EVI_median))) warning("All EVI values in background data are NA - check CRS or study extent")
      if (delete.intermediate.files.folders && "EVI" %in% datasets_requested) unlink(file.path(rasters.dir, "EVI"), recursive = TRUE, force = TRUE)
    }
  }

  # Calculate and process terrain metrics (Wilson et al. 2007)
  if ("terrain" %in% datasets_requested) {
    invisible(gc())
    counter <- counter + 1
    if (verbose) message("")
    if (verbose) message(sprintf("-- Downloading and extracting terrain metrics (ca. 30min): env.dataset %d of %d --",
                                 counter, total_datasets))
    terrain_variable_names <- c("TRI", "TPI", "roughness", "HLI", "TWI")
    terrain_dir <- file.path(rasters.dir, "terrain")
    if (!dir.exists(terrain_dir)) dir.create(terrain_dir, recursive = TRUE)
    terrain_occ_csv_file <- file.path(intermediate_files_dir, "Terrain_extracted_occurrence.csv")
    terrain_bg_csv_file <- file.path(intermediate_files_dir, "Terrain_extracted_background.csv")
    terrain_occ_exists <- file.exists(terrain_occ_csv_file)
    terrain_bg_exists <- file.exists(terrain_bg_csv_file)
    if (terrain_occ_exists && file.info(terrain_occ_csv_file)$size > 0 && (!generate.background.data || terrain_bg_exists) && !overwrite) {
      if (verbose) message("Terrain data already exist - skipping download and extraction")
      extracted_occurrences <- read.csv(terrain_occ_csv_file, row.names = 1, check.names = FALSE)
      extracted_occurrences <- extracted_occurrences[base_ids, , drop = FALSE]
      environmental_dataset <- cbind(environmental_dataset, extracted_occurrences)
      if (generate.background.data && terrain_bg_exists) {
        extracted_background <- read.csv(terrain_bg_csv_file, row.names = 1, check.names = FALSE)
        if (!is.null(background.data)) background.data <- cbind(background.data, extracted_background)
      } else if (generate.background.data && !terrain_bg_exists) {
        if (verbose) message("Terrain data for background is missing - downloading and extracting")
        terrain_occ_exists <- TRUE
        terrain_bg_exists <- FALSE
      }
    }
    if (!terrain_occ_exists || (generate.background.data && !terrain_bg_exists) || overwrite) {
      terra::terraOptions(tempdir = terrain_dir, todisk = TRUE, memfrac = 0.6, progress = 1)
      suppressWarnings(terra::tmpFiles(remove = TRUE))
      invisible(gc())
            if (is_europe_elevation) {
        elevation_raster_file <- file.path(rasters.dir, "elevation", "Copernicus_GLO90_Europe_250m.tif")
      } else if (is_asia_elevation) {
        elevation_raster_file <- file.path(rasters.dir, "elevation", "Copernicus_GLO90_Asia_250m.tif")
      } else if (is_north_america_elevation) {
        elevation_raster_file <- file.path(rasters.dir, "elevation", "Copernicus_GLO90_NorthAmerica_250m.tif")
      } else if (is_south_america_elevation) {
        elevation_raster_file <- file.path(rasters.dir, "elevation", "Copernicus_GLO90_SouthAmerica_250m.tif")
      } else if (is_africa_elevation) {
        elevation_raster_file <- file.path(rasters.dir, "elevation", "Copernicus_GLO90_Africa_250m.tif")
      } else if (is_australia_elevation) {
        elevation_raster_file <- file.path(rasters.dir, "elevation", "Copernicus_GLO90_Australia_250m.tif")
      } else if (is_new_world_elevation) {
        elevation_raster_file <- file.path(rasters.dir, "elevation", "Copernicus_GLO90_NewWorld_250m.tif")
      } else if (is_indo_pacific_elevation) {
        elevation_raster_file <- file.path(rasters.dir, "elevation", "Copernicus_GLO90_IndoPacific_250m.tif")
      } else if (is_eurasia_elevation) {
        elevation_raster_file <- file.path(rasters.dir, "elevation", "Copernicus_GLO90_Eurasia_250m.tif")
      } else if (is_holarctic_elevation) {
        elevation_raster_file <- file.path(rasters.dir, "elevation", "Copernicus_GLO90_Holarctic_250m.tif")
      } else {
        elevation_raster_file <- file.path(rasters.dir, "elevation", "Copernicus_GLO90_Global_250m.tif")
      }
      if (!file.exists(elevation_raster_file)) stop("Elevation raster required for terrain does not exist: ", elevation_raster_file)
      elevation_raster <- terra::rast(elevation_raster_file)
      elevation_raster <- terra::crop(elevation_raster, study_area_vect)
      elev_crop_file <- file.path(terrain_dir, "elevation_cropped.tif")
      terra::writeRaster(elevation_raster, elev_crop_file, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
      elevation_raster <- terra::rast(elev_crop_file)
      elevation_lonlat <- elevation_raster
      if (!terra::is.lonlat(elevation_lonlat)) elevation_lonlat <- terra::project(elevation_lonlat, "EPSG:4326", method = "bilinear")
      hydrology_crs <- if (terra::relate(study_area_vect, terra::ext(-170, -50, 10, 76), "intersects")[1]) {
        "EPSG:5070"
      } else {
        study_area_centroid <- terra::crds(terra::centroids(study_area_vect))
        sprintf("+proj=laea +lat_0=%f +lon_0=%f +datum=WGS84 +units=m +no_defs", study_area_centroid[2], study_area_centroid[1])
      }
      hydrology_buffer_cells <- 20
      hydrology_buffer_meters <- hydrology_buffer_cells * 250
      study_area_buffered <- terra::buffer(study_area_vect, width = hydrology_buffer_meters)
      elevation_projected_path <- file.path(terrain_dir, "elevation_projected_for_hydrology.tif")
      if (file.exists(elevation_projected_path) && file.info(elevation_projected_path)$size > 0 && !overwrite) {
        if (verbose) message("Using already saved projected elevation raster for hydrology")
        elevation_projected_for_hydrology <- terra::rast(elevation_projected_path)
      } else {
        if (verbose) message("Projecting elevation raster to meter-based CRS for hydrology")
        elevation_projected_for_hydrology <- terra::project(elevation_raster, hydrology_crs, method = "bilinear")
        terra::writeRaster(elevation_projected_for_hydrology, elevation_projected_path, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
      }
      elevation_filled_path <- file.path(terrain_dir, "elevation_filled.tif")
      elevation_nopits_path <- file.path(terrain_dir, "elevation_nopits.tif")
      elevation_breached_path <- file.path(terrain_dir, "elevation_breached.tif")
      slope_degrees_path <- file.path(terrain_dir, "slope_degrees.tif")
      flow_accumulation_path <- file.path(terrain_dir, "flow_accumulation_sca.tif")
      twi_projected_path <- file.path(terrain_dir, "topographic_wetness_index_projected.tif")
      terra::writeRaster(elevation_projected_for_hydrology, elevation_filled_path, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
      whitebox::wbt_init(workdir = terrain_dir)
      if (file.exists(elevation_nopits_path) && file.info(elevation_nopits_path)$size > 0 && !overwrite) {
        if (verbose) message("Using already saved filled single-cell pits raster")
      } else {
        invisible(capture.output(whitebox::wbt_fill_single_cell_pits(dem = elevation_filled_path, output = elevation_nopits_path, verbose_mode = FALSE)))
      }
      if (file.exists(elevation_breached_path) && file.info(elevation_breached_path)$size > 0 && !overwrite) {
        if (verbose) message("Using already saved breached depressions raster")
      } else {
        invisible(capture.output(whitebox::wbt_breach_depressions_least_cost(dem = elevation_nopits_path, output = elevation_breached_path, dist = 4, fill = FALSE, verbose_mode = FALSE)))
      }
      if (file.exists(slope_degrees_path) && file.info(slope_degrees_path)$size > 0 && !overwrite) {
        if (verbose) message("Using already saved slope raster")
      } else {
        invisible(capture.output(whitebox::wbt_slope(dem = elevation_breached_path, output = slope_degrees_path, units = "degrees", verbose_mode = FALSE)))
      }
      if (file.exists(flow_accumulation_path) && file.info(flow_accumulation_path)$size > 0 && !overwrite) {
        if (verbose) message("Using already saved flow accumulation raster")
      } else {
        invisible(capture.output(whitebox::wbt_d8_flow_accumulation(i = elevation_breached_path, output = flow_accumulation_path, out_type = "sca", verbose_mode = FALSE)))
      }
      if (file.exists(twi_projected_path) && file.info(twi_projected_path)$size > 0 && !overwrite) {
        if (verbose) message("Using already saved Topographic Wetness Index (TWI) raster")
      } else {
        invisible(capture.output(whitebox::wbt_wetness_index(sca = flow_accumulation_path, slope = slope_degrees_path, output = twi_projected_path, verbose_mode = FALSE)))
      }
      twi_projected_raster <- terra::rast(twi_projected_path)
      names(twi_projected_raster) <- "TWI"
      twi_lonlat_raster <- terra::project(twi_projected_raster, terra::crs(elevation_lonlat), method = "bilinear")
      twi_lonlat_path <- file.path(terrain_dir, "topographic_wetness_index_lonlat.tif")
      terra::writeRaster(twi_lonlat_raster, twi_lonlat_path, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
      TRI <- terra::terrain(elevation_lonlat, v = "TRI", neighbors = 8)
      TPI <- terra::terrain(elevation_lonlat, v = "TPI", neighbors = 8)
      roughness <- terra::terrain(elevation_lonlat, v = "roughness", neighbors = 8)
      slope <- terra::terrain(elevation_lonlat, v = "slope", unit = "radians")
      aspect <- terra::terrain(elevation_lonlat, v = "aspect", unit = "radians")
      latitude_degrees_raster <- terra::init(elevation_lonlat, fun = "y")
      latitude_radians_raster <- latitude_degrees_raster * pi / 180
      HLI <- (0.339 + 0.808 * cos(latitude_radians_raster) * cos(slope) - 0.196 * sin(latitude_radians_raster) * sin(slope) - 0.482 * cos(aspect - 225 * pi / 180) * sin(slope))
      HLI_file <- file.path(terrain_dir, "HLI.tif")
      terra::writeRaster(HLI, HLI_file, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
      terrain_raster_lonlat <- c(TRI, TPI, roughness, HLI)
      names(terrain_raster_lonlat) <- c("TRI", "TPI", "roughness", "HLI")
      result_lonlat <- extract.and.cache.env.dataset(dataset_name = "Terrain",
                                                     raster_object = terrain_raster_lonlat,
                                                     coord_env = coordinate_vector_env,
                                                     coord_bg = coordinate_vector_bg,
                                                     environmental_dataset = environmental_dataset,
                                                     background.data = background.data,
                                                     output.dir = intermediate_files_dir,
                                                     overwrite = overwrite,
                                                     generate.background.data = generate.background.data,
                                                     verbose = verbose)
      environmental_dataset <- result_lonlat$environmental_dataset
      background.data <- result_lonlat$background.data
      coord_env_proj <- if (terra::same.crs(coordinate_vector_env, twi_projected_raster)) coordinate_vector_env else terra::project(coordinate_vector_env, terra::crs(twi_projected_raster))
      coord_bg_proj <- if (is.null(coordinate_vector_bg)) NULL else if (terra::same.crs(coordinate_vector_bg, twi_projected_raster)) coordinate_vector_bg else terra::project(coordinate_vector_bg, terra::crs(twi_projected_raster))
      result_twi <- extract.and.cache.env.dataset(dataset_name = "Terrain_TWI",
                                                  raster_object = twi_projected_raster,
                                                  coord_env = coord_env_proj,
                                                  coord_bg = coord_bg_proj,
                                                  environmental_dataset = environmental_dataset,
                                                  background.data = background.data,
                                                  output.dir = intermediate_files_dir,
                                                  overwrite = overwrite,
                                                  generate.background.data = generate.background.data,
                                                  verbose = verbose)
      environmental_dataset <- result_twi$environmental_dataset
      background.data <- result_twi$background.data
      if (any(!terrain_variable_names %in% names(environmental_dataset))) stop("Terrain columns missing from occurrence data after extraction - extraction seems to not have worked")
      if (generate.background.data && any(!terrain_variable_names %in% names(background.data))) stop("Terrain columns missing from background data after extraction - extraction seems to not have worked")
      if (exists("terrain_raster")) rm(terrain_raster)
      if (exists("result")) rm(result)
      invisible(gc())
      suppressWarnings(try(terra::tmpFiles(remove = TRUE), silent = TRUE))
      if (delete.intermediate.files.folders) {
        if ("terrain" %in% datasets_requested) unlink(file.path(rasters.dir, "terrain"), recursive = TRUE, force = TRUE)
        if ("terrain" %in% datasets_requested && "elevation" %in% datasets_requested) unlink(file.path(rasters.dir, "elevation"), recursive = TRUE, force = TRUE)
      }
    }
  }

  # Download and process ENVIREM data (Title & Bemmels 2018; 30s resolution; download from Zenodo mirrors)
  if ("ENVIREM" %in% datasets_requested) {
    invisible(gc())
    counter <- counter + 1
    if (verbose) message("")
    if (verbose) message(sprintf("-- Downloading and extracting ENVIREM data (ca. 1h): env.dataset %d of %d --",
                                 counter, total_datasets))
    envirem_variable_names <- c("annualPET",
                                "climaticMoistureIndex",
                                "aridityIndexThornthwaite",
                                "continentality",
                                "embergerQ",
                                "growingDegDays0",
                                "growingDegDays5",
                                "maxTempColdest",
                                "minTempWarmest",
                                "monthCountByTemp10",
                                "PETColdestQuarter",
                                "PETDriestQuarter",
                                "PETseasonality",
                                "PETWarmestQuarter",
                                "PETWettestQuarter",
                                "thermicityIndex")
    envirem_dir <- file.path(rasters.dir, "ENVIREM")
    if (!dir.exists(envirem_dir)) dir.create(envirem_dir, recursive = TRUE)
    envirem_occ_csv_file <- file.path(intermediate_files_dir, "ENVIREM_extracted_occurrence.csv")
    envirem_bg_csv_file <- file.path(intermediate_files_dir, "ENVIREM_extracted_background.csv")
    envirem_occ_exists <- file.exists(envirem_occ_csv_file)
    envirem_bg_exists <- file.exists(envirem_bg_csv_file)
    if (envirem_occ_exists && (!generate.background.data || envirem_bg_exists) && !overwrite) {
      if (verbose) message("ENVIREM data already exist - skipping download and extraction")
      occ_envirem <- read.csv(envirem_occ_csv_file, row.names = 1, check.names = FALSE)
      occ_envirem <- occ_envirem[base_ids, , drop = FALSE]
      environmental_dataset <- cbind(environmental_dataset, occ_envirem)
      if (generate.background.data && envirem_bg_exists) {
        bg_envirem <- read.csv(envirem_bg_csv_file, row.names = 1, check.names = FALSE)
        if (!is.null(background.data)) background.data <- cbind(background.data, bg_envirem[rownames(background.data), , drop = FALSE])
      }
      if (all(is.na(environmental_dataset[, setdiff(names(occ_envirem), "ID")]))) warning("All ENVIREM values in occurrence data are NA - check CRS or extent")
      if (generate.background.data && all(is.na(background.data[, setdiff(names(bg_envirem), "ID")]))) warning("All ENVIREM values in background data are NA - check CRS or extent")
    } else {
      europe_extent_envirem <- terra::ext(-25, 41, 35, 71)
      is_europe_envirem <- terra::relate(study_area_vect, europe_extent_envirem, relation = "within")[1]
      asia_extent_envirem <- terra::ext(39, 177, -9, 80)
      is_asia_envirem <- terra::relate(study_area_vect, asia_extent_envirem, relation = "within")[1]
      eurasia_extent_envirem <- terra::ext(-25, 177, -9, 80)
      is_eurasia_envirem <- terra::relate(study_area_vect, eurasia_extent_envirem, relation = "within")[1]
      north_america_extent_envirem <- terra::ext(-170, -50, 10, 76)
      is_north_america_envirem <- terra::relate(study_area_vect, north_america_extent_envirem, relation = "within")[1]
      south_america_extent_envirem <- terra::ext(-90, -30, -59, 12)
      is_south_america_envirem <- terra::relate(study_area_vect, south_america_extent_envirem, relation = "within")[1]
      africa_extent_envirem <- terra::ext(-20, 51, -34, 37)
      is_africa_envirem <- terra::relate(study_area_vect, africa_extent_envirem, relation = "within")[1]
      australia_extent_envirem <- terra::ext(110, 176, -49, -7)
      is_australia_envirem <- terra::relate(study_area_vect, australia_extent_envirem, relation = "within")[1]
      indo_pacific_extent_envirem <- terra::ext(39, 177, -49, 80)
      is_indo_pacific_envirem <- terra::relate(study_area_vect, indo_pacific_extent_envirem, relation = "within")[1]
      holarctic_extent_envirem <- terra::ext(-170, 177, -9, 80)
      is_holarctic_envirem <- terra::relate(study_area_vect, holarctic_extent_envirem, relation = "within")[1]
      new_world_extent_envirem <- terra::ext(-170, -30, -59, 76)
      is_new_world_envirem <- terra::relate(study_area_vect, new_world_extent_envirem, relation = "within")[1]
      old_world_extent_envirem <- terra::ext(-25, 177, -34, 80)
      is_old_world_envirem <- terra::relate(study_area_vect, old_world_extent_envirem, relation = "within")[1]
      envirem_zip_file <- file.path(envirem_dir, "ENVIREM_selected.zip")
      existing_tifs <- list.files(envirem_dir, pattern = "\\.tif$", full.names = TRUE)
      if (length(existing_tifs) == 16 && all(file.size(existing_tifs) > 1e5) && !redownload.rasters) {
        if (verbose) message("All 16 ENVIREM rasters already present - skipping download and extraction")
      } else {
        if (is_north_america_envirem) {
          if (verbose) message("Using North America ENVIREM rasters (1km resolution)")
          envirem_zip_url <- "https://zenodo.org/records/17507777/files/ENVIREM_north_america.zip?download=1"; min_size_mb <- 1170
        } else if (is_south_america_envirem) {
          if (verbose) message("Using South America ENVIREM rasters (1km resolution)")
          envirem_zip_url <- "https://zenodo.org/records/17507777/files/ENVIREM_south_america.zip?download=1"; min_size_mb <- 680
        } else if (is_africa_envirem) {
          if (verbose) message("Using Africa ENVIREM rasters (1km resolution)")
          envirem_zip_url <- "https://zenodo.org/records/17507777/files/ENVIREM_africa.zip?download=1"; min_size_mb <- 1170
        } else if (is_europe_envirem) {
          if (verbose) message("Using Europe ENVIREM rasters (1km resolution)")
          envirem_zip_url <- "https://zenodo.org/records/17507777/files/ENVIREM_europe.zip?download=1"; min_size_mb <- 470
        } else if (is_asia_envirem) {
          if (verbose) message("Using Asia ENVIREM rasters (1km resolution)")
          envirem_zip_url <- "https://zenodo.org/records/17507777/files/ENVIREM_asia.zip?download=1"; min_size_mb <- 2340
        } else if (is_eurasia_envirem) {
          if (verbose) message("Using Eurasia ENVIREM rasters (1km resolution)")
          envirem_zip_url <- "https://zenodo.org/records/17507777/files/ENVIREM_eurasia.zip?download=1"; min_size_mb <- 3510
        } else if (is_australia_envirem) {
          if (verbose) message("Using Australia ENVIREM rasters (1km resolution)")
          envirem_zip_url <- "https://zenodo.org/records/17507777/files/ENVIREM_australia.zip?download=1"; min_size_mb <- 310
        } else if (is_indo_pacific_envirem) {
          if (verbose) message("Using Indo-Pacific ENVIREM rasters (1km resolution)")
          envirem_zip_url <- "https://zenodo.org/records/17507777/files/ENVIREM_indo_pacific.zip?download=1"; min_size_mb <- 2700
        } else if (is_old_world_envirem) {
          if (verbose) message("Using Old World ENVIREM rasters (1km resolution)")
          envirem_zip_url <- "https://zenodo.org/records/17507777/files/ENVIREM_old_world.zip?download=1"; min_size_mb <- 4050
        } else if (is_new_world_envirem) {
          if (verbose) message("Using New World ENVIREM rasters (1km resolution)")
          envirem_zip_url <- "https://zenodo.org/records/17507777/files/ENVIREM_new_world.zip?download=1"; min_size_mb <- 1890
        } else if (is_holarctic_envirem) {
          if (verbose) message("Using Holarctic ENVIREM rasters (1km resolution)")
          envirem_zip_url <- "https://zenodo.org/records/17507777/files/ENVIREM_holarctic.zip?download=1"; min_size_mb <- 5130
        } else {
          if (verbose) message("Using Global ENVIREM rasters (1km resolution)")
          envirem_zip_url <- "https://zenodo.org/records/17507777/files/ENVIREM_global.zip?download=1"; min_size_mb <- 6210
        }
        if (!file.exists(envirem_zip_file) || file.size(envirem_zip_file) < min_size_mb * 1e6 || redownload.rasters) {
          robust.download.raster(envirem_zip_url, envirem_zip_file, min_size_mb = min_size_mb)
        }
        if (file.size(envirem_zip_file) < min_size_mb * 1e6) stop("ENVIREM zip appears incomplete - download failed or truncated")
        utils::unzip(envirem_zip_file, exdir = envirem_dir, overwrite = TRUE)
        existing_tifs <- list.files(envirem_dir, pattern = "\\.tif$", full.names = TRUE)
        if (length(existing_tifs) != 16) {
          if (verbose) message("Extraction incomplete - retrying download once")
          file.remove(envirem_zip_file)
          robust.download.raster(envirem_zip_url, envirem_zip_file, min_size_mb = min_size_mb)
          utils::unzip(envirem_zip_file, exdir = envirem_dir, overwrite = TRUE)
          existing_tifs <- list.files(envirem_dir, pattern = "\\.tif$", full.names = TRUE)
          if (length(existing_tifs) != 16) stop("Extraction failed again - expected 16 ENVIREM .tif layers")
        }
        if (file.exists(envirem_zip_file)) file.remove(envirem_zip_file)
      }
      envirem_tifs <- list.files(envirem_dir, pattern = "\\.tif$", full.names = TRUE)
      envirem_raster <- terra::rast(envirem_tifs)
      names(envirem_raster) <- envirem_variable_names
      raster_crs <- terra::crs(envirem_raster)
      coord_env <- coordinate_vector_env
      coord_bg <- coordinate_vector_bg
      if (!terra::same.crs(coordinate_vector_env, envirem_raster)) {
        if (verbose) message("Projecting coordinates to match ENVIREM raster CRS")
        coord_env <- terra::project(coordinate_vector_env, raster_crs)
        if (!is.null(coordinate_vector_bg)) coord_bg <- terra::project(coordinate_vector_bg, raster_crs)
      }
      result <- extract.and.cache.env.dataset(dataset_name = "ENVIREM",
                                              raster_object = envirem_raster,
                                              coord_env = coord_env,
                                              coord_bg = coord_bg,
                                              environmental_dataset = environmental_dataset,
                                              background.data = background.data,
                                              output.dir = intermediate_files_dir,
                                              overwrite = overwrite,
                                              generate.background.data = generate.background.data,
                                              verbose = verbose)
      environmental_dataset <- result$environmental_dataset
      background.data <- result$background.data
      if (any(!envirem_variable_names %in% names(environmental_dataset))) stop("ENVIREM columns missing from occurrence data after extraction - extraction seems to not have worked")
      if (generate.background.data && any(!envirem_variable_names %in% names(background.data))) stop("ENVIREM columns missing from background data after extraction - extraction seems to not have worked")
      if (all(is.na(environmental_dataset$annualPET))) warning("All ENVIREM values in occurrence data are NA - check CRS or study extent")
      if (generate.background.data && all(is.na(background.data$annualPET))) warning("All ENVIREM values in background data are NA - check CRS or study extent")
      if (exists("envirem_raster")) rm(envirem_raster)
      if (exists("result")) rm(result)
      invisible(gc())
      if (delete.intermediate.files.folders && "ENVIREM" %in% datasets_requested)
        unlink(file.path(rasters.dir, "ENVIREM"), recursive = TRUE, force = TRUE)
    }
  }

  # Download and process Human Footprint (HFP) data (Venter et al. 2016; 30s resolution; for 2009; download using geodata package)
  if ("footprint" %in% datasets_requested) {
    invisible(gc())
    counter <- counter + 1
    if (verbose) message("")
    if (verbose) message(sprintf("-- Downloading and extracting Human Footprint data (ca. 1 min): env.dataset %d of %d --",
                                 counter, total_datasets))
    footprint_variable_name <- "Footprint"
    footprint_dir <- file.path(rasters.dir, "footprint")
    if (!dir.exists(footprint_dir)) dir.create(footprint_dir, recursive = TRUE)
    footprint_raster_file <- file.path(footprint_dir, "footprint_2009.tif")
    footprint_occ_csv_file <- file.path(intermediate_files_dir, "Footprint_extracted_occurrence.csv")
    footprint_bg_csv_file <- file.path(intermediate_files_dir, "Footprint_extracted_background.csv")
    footprint_occ_exists <- file.exists(footprint_occ_csv_file)
    footprint_bg_exists <- file.exists(footprint_bg_csv_file)
    footprint_geodata_cache_dir <- file.path(footprint_dir, "landuse")
    if (file.exists(footprint_raster_file) && file.size(footprint_raster_file) > 1e6 && dir.exists(footprint_geodata_cache_dir)) unlink(footprint_geodata_cache_dir, recursive = TRUE, force = TRUE)
    if (footprint_occ_exists && (!generate.background.data || footprint_bg_exists) && !overwrite) {
      if (verbose) message("Human Footprint data already exist - skipping download and extraction")
      footprint_occurrence_values <- read.csv(footprint_occ_csv_file, row.names = 1, check.names = FALSE)
      footprint_occurrence_values <- footprint_occurrence_values[base_ids, , drop = FALSE]
      environmental_dataset <- cbind(environmental_dataset, footprint_occurrence_values)
      if (generate.background.data && footprint_bg_exists) {
        footprint_background_values <- read.csv(footprint_bg_csv_file, row.names = 1, check.names = FALSE)
        if (!is.null(background.data)) background.data <- cbind(background.data, footprint_background_values[rownames(background.data), , drop = FALSE])
      }
      if (all(is.na(environmental_dataset[, setdiff(names(footprint_occurrence_values), "ID")])))
        warning("All Human Footprint values in occurrence data are NA - check CRS or extent")
      if (generate.background.data && all(is.na(background.data[, setdiff(names(footprint_background_values), "ID")])))
        warning("All Human Footprint values in background data are NA - check CRS or extent")
    } else {
      footprint_geodata_cache_dir <- file.path(footprint_dir, "landuse")
      if (!file.exists(footprint_raster_file) || file.size(footprint_raster_file) < 1e6 || redownload.rasters) {
        if (dir.exists(footprint_geodata_cache_dir)) unlink(footprint_geodata_cache_dir, recursive = TRUE, force = TRUE)
        if (file.exists(footprint_raster_file)) unlink(footprint_raster_file, force = TRUE)
        footprint_url <- "https://geodata.ucdavis.edu/geodata/footprint/wildareas-v3-2009-human-footprint_geo.tif"
        utils::download.file(footprint_url, destfile = footprint_raster_file, mode = "wb", quiet = TRUE)
        if (!file.exists(footprint_raster_file) || file.size(footprint_raster_file) < 1e6) stop("Human Footprint download failed or produced an invalid file")
        footprint_raster <- terra::rast(footprint_raster_file)
      } else {
        if (verbose) message("Human Footprint raster already present - skipping download")
        footprint_raster <- terra::rast(footprint_raster_file)
      }
      if (file.exists(footprint_raster_file) && file.size(footprint_raster_file) > 1e6 && dir.exists(footprint_geodata_cache_dir)) unlink(footprint_geodata_cache_dir, recursive = TRUE, force = TRUE)
      names(footprint_raster) <- footprint_variable_name
      raster_crs <- terra::crs(footprint_raster)
      coord_env <- coordinate_vector_env
      coord_bg <- coordinate_vector_bg
      if (!terra::same.crs(coordinate_vector_env, footprint_raster)) {
        if (verbose) message("Projecting coordinates to match Human Footprint raster CRS")
        coord_env <- terra::project(coordinate_vector_env, raster_crs)
        if (!is.null(coordinate_vector_bg)) coord_bg <- terra::project(coordinate_vector_bg, raster_crs)
      }
      extracted_occurrences <- as.matrix(terra::extract(footprint_raster, coord_env, ID = FALSE))
      colnames(extracted_occurrences) <- footprint_variable_name
      environmental_dataset <- cbind(environmental_dataset, extracted_occurrences)
      write.csv(environmental_dataset[, footprint_variable_name, drop = FALSE], footprint_occ_csv_file)
      if (!is.null(background.data) && generate.background.data) {
        extracted_background <- as.matrix(terra::extract(footprint_raster, coord_bg, ID = FALSE))
        colnames(extracted_background) <- footprint_variable_name
        background.data <- cbind(background.data, extracted_background)
        write.csv(background.data[, footprint_variable_name, drop = FALSE], footprint_bg_csv_file)
      }
      if (any(!footprint_variable_name %in% names(environmental_dataset))) stop("Human Footprint column missing from occurrence data after extraction - extraction seems to not have worked")
      if (generate.background.data && any(!footprint_variable_name %in% names(background.data))) stop("Human Footprint column missing from background data after extraction - extraction seems to not have worked")
      suppressWarnings(try(terra::tmpFiles(remove = TRUE), silent = TRUE))
      if (exists("footprint_raster")) rm(footprint_raster)
      if (exists("extracted_occurrences")) rm(extracted_occurrences)
      if (exists("extracted_background")) rm(extracted_background)
      invisible(gc())
      if (all(is.na(environmental_dataset[[footprint_variable_name]]))) warning("All Human Footprint values in occurrence data are NA - check CRS or study extent")
      if (generate.background.data && all(is.na(background.data[[footprint_variable_name]]))) warning("All Human Footprint values in background data are NA - check CRS or study extent")
      if (delete.intermediate.files.folders && "footprint" %in% datasets_requested) unlink(file.path(rasters.dir, "footprint"), recursive = TRUE, force = TRUE)
    }
  }

  # Download and process landcover data (ESA WorldCover; 30s resolution; regional Zenodo subsets)
  if ("landcover" %in% datasets_requested) {
    invisible(gc())
    counter <- counter + 1
    if (verbose) message("")
    if (verbose) message(sprintf("-- Downloading and extracting landcover data (ca 5 min): env.dataset %d of %d --",
                                 counter, total_datasets))
    landcover_variable_names <- c("trees", "grassland", "shrubs", "cropland", "built", "bare", "snow", "water", "wetland", "mangroves", "moss")
    landcover_dir <- file.path(rasters.dir, "landcover")
    if (!dir.exists(landcover_dir)) dir.create(landcover_dir, recursive = TRUE)
    landcover_occ_csv_file <- file.path(intermediate_files_dir, "Landcover_extracted_occurrence.csv")
    landcover_bg_csv_file <- file.path(intermediate_files_dir, "Landcover_extracted_background.csv")
    landcover_occ_exists <- file.exists(landcover_occ_csv_file) && file.info(landcover_occ_csv_file)$size > 0
    landcover_bg_exists <- file.exists(landcover_bg_csv_file) && file.info(landcover_bg_csv_file)$size > 0
    if (landcover_occ_exists && (!generate.background.data || landcover_bg_exists) && !overwrite) {
      if (verbose) message("Landcover data already exist - skipping download and extraction - loading from intermediate files")
      occ_landcover <- read.csv(landcover_occ_csv_file, row.names = 1, check.names = FALSE)
      occ_landcover <- occ_landcover[base_ids, , drop = FALSE]
      environmental_dataset <- cbind(environmental_dataset, occ_landcover)
      if (generate.background.data && landcover_bg_exists) {
        bg_landcover <- read.csv(landcover_bg_csv_file, row.names = 1, check.names = FALSE)
        if (!is.null(background.data)) background.data <- cbind(background.data, bg_landcover[rownames(background.data), , drop = FALSE])
      }
      if (all(is.na(environmental_dataset[, setdiff(names(occ_landcover), "ID")]))) warning("All landcover values in occurrence data are NA - check CRS or extent")
      if (generate.background.data && all(is.na(background.data[, setdiff(names(bg_landcover), "ID")]))) warning("All landcover values in background data are NA - check CRS or extent")
    } else {
      europe_extent_landcover <- list(terra::ext(-26, 65, 34, 83))
      asia_extent_landcover <- list(terra::ext(21, 180, -13, 84), terra::ext(-180, -168, 50, 75))
      northamerica_extent_landcover <- list(terra::ext(-180, -10, 5, 85))
      southamerica_extent_landcover <- list(terra::ext(-92, -30, -56, 15))
      africa_extent_landcover <- list(terra::ext(-26, 57, -36, 39))
      australia_extent_landcover <- list(terra::ext(90, 180, -56, 25), terra::ext(-180, -170, -56, 25))
      indopacific_extent_landcover <- c(asia_extent_landcover, australia_extent_landcover)
      eurasia_extent_landcover <- c(europe_extent_landcover, asia_extent_landcover)
      holarctic_extent_landcover <- c(europe_extent_landcover, asia_extent_landcover, northamerica_extent_landcover)
      newworld_extent_landcover <- c(northamerica_extent_landcover, southamerica_extent_landcover)
      oldworld_extent_landcover <- c(europe_extent_landcover, africa_extent_landcover, asia_extent_landcover, australia_extent_landcover)
      is.within.landcover.extent <- function(study_area, extent_list) {
        any(vapply(extent_list, function(ext) terra::relate(study_area, ext, relation = "within")[1], logical(1)))
      }
      is_europe_landcover <- is.within.landcover.extent(study_area_vect, europe_extent_landcover)
      is_asia_landcover <- is.within.landcover.extent(study_area_vect, asia_extent_landcover)
      is_northamerica_landcover <- is.within.landcover.extent(study_area_vect, northamerica_extent_landcover)
      is_southamerica_landcover <- is.within.landcover.extent(study_area_vect, southamerica_extent_landcover)
      is_africa_landcover <- is.within.landcover.extent(study_area_vect, africa_extent_landcover)
      is_australia_landcover <- is.within.landcover.extent(study_area_vect, australia_extent_landcover)
      is_indopacific_landcover <- is.within.landcover.extent(study_area_vect, indopacific_extent_landcover)
      is_eurasia_landcover <- is.within.landcover.extent(study_area_vect, eurasia_extent_landcover)
      is_holarctic_landcover <- is.within.landcover.extent(study_area_vect, holarctic_extent_landcover)
      is_newworld_landcover <- is.within.landcover.extent(study_area_vect, newworld_extent_landcover)
      is_oldworld_landcover <- is.within.landcover.extent(study_area_vect, oldworld_extent_landcover)
      if (is_europe_landcover) {
        landcover_subset_name <- "europe"
        landcover_zip_url <- "https://zenodo.org/records/19600289/files/landcover_europe.zip?download=1"
        landcover_zip_min_size_mb <- 300
      } else if (is_asia_landcover) {
        landcover_subset_name <- "asia"
        landcover_zip_url <- "https://zenodo.org/records/19600289/files/landcover_asia.zip?download=1"
        landcover_zip_min_size_mb <- 900
      } else if (is_northamerica_landcover) {
        landcover_subset_name <- "northamerica"
        landcover_zip_url <- "https://zenodo.org/records/19600289/files/landcover_northamerica.zip?download=1"
        landcover_zip_min_size_mb <- 550
      } else if (is_southamerica_landcover) {
        landcover_subset_name <- "southamerica"
        landcover_zip_url <- "https://zenodo.org/records/19600289/files/landcover_southamerica.zip?download=1"
        landcover_zip_min_size_mb <- 180
      } else if (is_africa_landcover) {
        landcover_subset_name <- "africa"
        landcover_zip_url <- "https://zenodo.org/records/19600289/files/landcover_africa.zip?download=1"
        landcover_zip_min_size_mb <- 280
      } else if (is_australia_landcover) {
        landcover_subset_name <- "australia"
        landcover_zip_url <- "https://zenodo.org/records/19600289/files/landcover_australia.zip?download=1"
        landcover_zip_min_size_mb <- 150
      } else if (is_indopacific_landcover) {
        landcover_subset_name <- "indopacific"
        landcover_zip_url <- "https://zenodo.org/records/19600289/files/landcover_indopacific.zip?download=1"
        landcover_zip_min_size_mb <- 1000
      } else if (is_eurasia_landcover) {
        landcover_subset_name <- "eurasia"
        landcover_zip_url <- "https://zenodo.org/records/19600289/files/landcover_eurasia.zip?download=1"
        landcover_zip_min_size_mb <- 1000
      } else if (is_holarctic_landcover) {
        landcover_subset_name <- "holarctic"
        landcover_zip_url <- "https://zenodo.org/records/19600289/files/landcover_holarctic.zip?download=1"
        landcover_zip_min_size_mb <- 1500
      } else if (is_newworld_landcover) {
        landcover_subset_name <- "newworld"
        landcover_zip_url <- "https://zenodo.org/records/19600289/files/landcover_newworld.zip?download=1"
        landcover_zip_min_size_mb <- 700
      } else if (is_oldworld_landcover) {
        landcover_subset_name <- "oldworld"
        landcover_zip_url <- "https://zenodo.org/records/19600289/files/landcover_oldworld.zip?download=1"
        landcover_zip_min_size_mb <- 1200
      } else {
        landcover_subset_name <- "global"
        output_files <- file.path(landcover_dir, paste0("landcover_", landcover_variable_names, ".tif"))
        valid_files <- sapply(output_files, function(f) file.exists(f) && file.size(f) > 3e6)
        all_valid <- all(valid_files)
        if (all_valid && !redownload.rasters) {
          if (verbose) message("Study area does not fit any regional landcover subset extent - using global landcover rasters already present")
        } else {
          for (landcover_variable in landcover_variable_names) {
            landcover_url <- sprintf("https://geodata.ucdavis.edu/geodata/landuse/WorldCover_%s_30s.tif", landcover_variable)
            output_file <- file.path(landcover_dir, paste0("landcover_", landcover_variable, ".tif"))
            if (!file.exists(output_file) || file.size(output_file) < 3e6 || redownload.rasters) {
              if (verbose) message("Study area does not fit any regional landcover subset extent - downloading global landcover raster: ", landcover_variable)
              robust.download.raster(landcover_url, output_file, min_size_mb = 30)
            } else if (verbose) {
              message("Global landcover raster already present: ", landcover_variable)
            }
          }
        }
      }
      if (landcover_subset_name != "global") {
        landcover_subset_dir <- file.path(landcover_dir, landcover_subset_name)
        if (!dir.exists(landcover_subset_dir)) dir.create(landcover_subset_dir, recursive = TRUE)
        landcover_zip_file <- file.path(landcover_dir, paste0("landcover_", landcover_subset_name, ".zip"))
        output_files <- file.path(landcover_subset_dir, paste0("landcover_", landcover_variable_names, "_", landcover_subset_name, ".tif"))
        valid_files <- sapply(output_files, function(f) file.exists(f) && file.size(f) > 1e5)
        all_valid <- all(valid_files)
        if (all_valid && !redownload.rasters && file.exists(landcover_zip_file)) file.remove(landcover_zip_file)
        if (all_valid && !redownload.rasters) {
          if (verbose) message("Landcover subset rasters already present - skipping download")
        } else {
          if (!file.exists(landcover_zip_file) || file.size(landcover_zip_file) < landcover_zip_min_size_mb * 1e6 || redownload.rasters) {
            if (verbose) message("Downloading landcover subset zip: ", landcover_subset_name)
            robust.download.raster(landcover_zip_url, landcover_zip_file, min_size_mb = landcover_zip_min_size_mb)
          } else if (verbose) {
            message("Landcover subset zip already present: ", landcover_subset_name)
          }
          utils::unzip(landcover_zip_file, exdir = landcover_subset_dir, overwrite = TRUE)
          if (any(!file.exists(output_files))) stop("Not all landcover subset raster files were found after unzip for subset: ", landcover_subset_name)
          if (file.exists(landcover_zip_file)) file.remove(landcover_zip_file)
        }
      }
      landcover_raster_list <- lapply(output_files, terra::rast)
      landcover_raster_stack <- terra::rast(landcover_raster_list)
      names(landcover_raster_stack) <- landcover_variable_names
      raster_crs <- terra::crs(landcover_raster_stack)
      coord_env <- coordinate_vector_env
      coord_bg <- coordinate_vector_bg
      if (!terra::same.crs(coordinate_vector_env, landcover_raster_stack)) {
        if (verbose) message("Projecting coordinates to match landcover raster CRS")
        coord_env <- terra::project(coordinate_vector_env, raster_crs)
        if (!is.null(coordinate_vector_bg)) coord_bg <- terra::project(coordinate_vector_bg, raster_crs)
      }
      extracted_occurrences <- as.matrix(terra::extract(landcover_raster_stack, coord_env, ID = FALSE))
      colnames(extracted_occurrences) <- landcover_variable_names
      environmental_dataset <- cbind(environmental_dataset, extracted_occurrences)
      write.csv(environmental_dataset[, landcover_variable_names, drop = FALSE], landcover_occ_csv_file)
      if (!is.null(background.data) && generate.background.data) {
        extracted_background <- as.matrix(terra::extract(landcover_raster_stack, coord_bg, ID = FALSE))
        colnames(extracted_background) <- landcover_variable_names
        background.data <- cbind(background.data, extracted_background)
        write.csv(background.data[, landcover_variable_names, drop = FALSE], landcover_bg_csv_file)
      }
      if (any(!landcover_variable_names %in% names(environmental_dataset))) stop("Landcover columns missing from occurrence data after extraction - extraction seems to not have worked")
      if (generate.background.data && any(!landcover_variable_names %in% names(background.data))) stop("Landcover columns missing from background data after extraction - extraction seems to not have worked")
      suppressWarnings(try(terra::tmpFiles(remove = TRUE), silent = TRUE))
      if (exists("landcover_raster_stack")) rm(landcover_raster_stack)
      if (exists("extracted_occurrences")) rm(extracted_occurrences)
      if (exists("extracted_background")) rm(extracted_background)
      invisible(gc())
      if (all(is.na(environmental_dataset[, landcover_variable_names]))) warning("All landcover values in occurrence data are NA - check CRS or study extent")
      if (generate.background.data && all(is.na(background.data[, landcover_variable_names]))) warning("All landcover values in background data are NA - check CRS or study extent")
      if (delete.intermediate.files.folders && "landcover" %in% datasets_requested)
        unlink(file.path(rasters.dir, "landcover"), recursive = TRUE, force = TRUE)
    }
  }

  # Download and process Soil variables (derived from SoilGrids; 30s resolution; mean 0-5 cm depth)
  if ("soil" %in% datasets_requested) {
    invisible(gc())
    counter <- counter + 1
    if (verbose) message("")
    if (verbose) message(sprintf("-- Downloading and extracting soil data (ca. 5 min): env.dataset %d of %d --",
                                 counter, total_datasets))
    soil_variable_names <- c("Bulk_density",
                             "Coarse_fragments_volume",
                             "Clay_fraction",
                             "Nitrogen_content",
                             "Organic_carb_density",
                             "pH_H2O",
                             "Sand_fraction",
                             "Silt_fraction",
                             "Soil_organic_carb")
    soil_source_variables <- c("bdod", "cfvo", "clay", "nitrogen", "ocd", "phh2o", "sand", "silt", "soc")
    soil_dir <- file.path(rasters.dir, "soil")
    if (!dir.exists(soil_dir)) dir.create(soil_dir, recursive = TRUE)
    soil_occ_csv_file <- file.path(intermediate_files_dir, "Soil_extracted_occurrence.csv")
    soil_bg_csv_file <- file.path(intermediate_files_dir, "Soil_extracted_background.csv")
    soil_occ_exists <- file.exists(soil_occ_csv_file) && file.info(soil_occ_csv_file)$size > 0
    soil_bg_exists <- file.exists(soil_bg_csv_file) && file.info(soil_bg_csv_file)$size > 0
    if (soil_occ_exists && (!generate.background.data || soil_bg_exists) && !overwrite) {
      if (verbose) message("Soil data already exist - skipping download and extraction")
      soil_occurrence_values <- read.csv(soil_occ_csv_file, row.names = 1, check.names = FALSE)
      soil_occurrence_values <- soil_occurrence_values[base_ids, , drop = FALSE]
      environmental_dataset <- cbind(environmental_dataset, soil_occurrence_values)
      if (generate.background.data && soil_bg_exists) {
        soil_background_values <- read.csv(soil_bg_csv_file, row.names = 1, check.names = FALSE)
        if (!is.null(background.data)) background.data <- cbind(background.data, soil_background_values[rownames(background.data), , drop = FALSE])
      }
      if (all(is.na(environmental_dataset[, setdiff(names(soil_occurrence_values), "ID")]))) warning("All Soil values in occurrence data are NA - check CRS or extent")
      if (generate.background.data && all(is.na(background.data[, setdiff(names(soil_background_values), "ID")]))) warning("All Soil values in background data are NA - check CRS or extent")
    } else {
      europe_extent_soil <- list(terra::ext(-26, 65, 34, 83))
      asia_extent_soil <- list(terra::ext(21, 180, -13, 84), terra::ext(-180, -168, 50, 75))
      northamerica_extent_soil <- list(terra::ext(-180, -10, 5, 85))
      southamerica_extent_soil <- list(terra::ext(-92, -30, -56, 15))
      africa_extent_soil <- list(terra::ext(-26, 57, -36, 39))
      australia_extent_soil <- list(terra::ext(90, 180, -56, 25), terra::ext(-180, -170, -56, 25))
      indopacific_extent_soil <- c(asia_extent_soil, australia_extent_soil)
      eurasia_extent_soil <- c(europe_extent_soil, asia_extent_soil)
      holarctic_extent_soil <- c(europe_extent_soil, asia_extent_soil, northamerica_extent_soil)
      newworld_extent_soil <- c(northamerica_extent_soil, southamerica_extent_soil)
      oldworld_extent_soil <- c(europe_extent_soil, africa_extent_soil, asia_extent_soil, australia_extent_soil)
      is.within.soil.extent <- function(study_area, extent_list) {
        any(vapply(extent_list, function(ext) terra::relate(study_area, ext, relation = "within")[1], logical(1)))
      }
      is_europe_soil <- is.within.soil.extent(study_area_vect, europe_extent_soil)
      is_asia_soil <- is.within.soil.extent(study_area_vect, asia_extent_soil)
      is_northamerica_soil <- is.within.soil.extent(study_area_vect, northamerica_extent_soil)
      is_southamerica_soil <- is.within.soil.extent(study_area_vect, southamerica_extent_soil)
      is_africa_soil <- is.within.soil.extent(study_area_vect, africa_extent_soil)
      is_australia_soil <- is.within.soil.extent(study_area_vect, australia_extent_soil)
      is_indopacific_soil <- is.within.soil.extent(study_area_vect, indopacific_extent_soil)
      is_eurasia_soil <- is.within.soil.extent(study_area_vect, eurasia_extent_soil)
      is_holarctic_soil <- is.within.soil.extent(study_area_vect, holarctic_extent_soil)
      is_newworld_soil <- is.within.soil.extent(study_area_vect, newworld_extent_soil)
      is_oldworld_soil <- is.within.soil.extent(study_area_vect, oldworld_extent_soil)
      if (is_europe_soil) {
        soil_subset_name <- "europe"
        soil_zip_url <- "https://zenodo.org/records/19614207/files/soil_europe.zip?download=1"
        soil_zip_min_size_mb <- 140
      } else if (is_asia_soil) {
        soil_subset_name <- "asia"
        soil_zip_url <- "https://zenodo.org/records/19614207/files/soil_asia.zip?download=1"
        soil_zip_min_size_mb <- 520
      } else if (is_northamerica_soil) {
        soil_subset_name <- "northamerica"
        soil_zip_url <- "https://zenodo.org/records/19614207/files/soil_northamerica.zip?download=1"
        soil_zip_min_size_mb <- 240
      } else if (is_southamerica_soil) {
        soil_subset_name <- "southamerica"
        soil_zip_url <- "https://zenodo.org/records/19614207/files/soil_southamerica.zip?download=1"
        soil_zip_min_size_mb <- 100
      } else if (is_africa_soil) {
        soil_subset_name <- "africa"
        soil_zip_url <- "https://zenodo.org/records/19614207/files/soil_africa.zip?download=1"
        soil_zip_min_size_mb <- 180
      } else if (is_australia_soil) {
        soil_subset_name <- "australia"
        soil_zip_url <- "https://zenodo.org/records/19614207/files/soil_australia.zip?download=1"
        soil_zip_min_size_mb <- 70
      } else if (is_indopacific_soil) {
        soil_subset_name <- "indopacific"
        soil_zip_url <- "https://zenodo.org/records/19614207/files/soil_indopacific.zip?download=1"
        soil_zip_min_size_mb <- 600
      } else if (is_eurasia_soil) {
        soil_subset_name <- "eurasia"
        soil_zip_url <- "https://zenodo.org/records/19614207/files/soil_eurasia.zip?download=1"
        soil_zip_min_size_mb <- 600
      } else if (is_holarctic_soil) {
        soil_subset_name <- "holarctic"
        soil_zip_url <- "https://zenodo.org/records/19614207/files/soil_holarctic.zip?download=1"
        soil_zip_min_size_mb <- 900
      } else if (is_newworld_soil) {
        soil_subset_name <- "newworld"
        soil_zip_url <- "https://zenodo.org/records/19614207/files/soil_newworld.zip?download=1"
        soil_zip_min_size_mb <- 400
      } else if (is_oldworld_soil) {
        soil_subset_name <- "oldworld"
        soil_zip_url <- "https://zenodo.org/records/19614207/files/soil_oldworld.zip?download=1"
        soil_zip_min_size_mb <- 700
      } else {
        soil_subset_name <- "global"
        if (verbose) message("Downloading global layer")
        output_files <- file.path(soil_dir, paste0("soil_", soil_source_variables, ".tif"))
        valid_files <- sapply(output_files, function(f) file.exists(f) && file.size(f) > 1e6)
        all_valid <- all(valid_files)
        if (all_valid && !redownload.rasters) {
          if (verbose) message("Study area does not fit any regional soil subset extent - using global soil rasters already present")
        } else {
          soil_raster_list <- list()
          for (i in seq_along(soil_source_variables)) {
            source_var <- soil_source_variables[i]
            output_file <- file.path(soil_dir, paste0("soil_", source_var, ".tif"))
            soil_layer <- geodata::soil_world(var = source_var, depth = 5, stat = "mean", path = soil_dir, download = TRUE)
            terra::writeRaster(soil_layer, filename = output_file, overwrite = TRUE)
            soil_raster_list[[soil_variable_names[i]]] <- soil_layer
          }
        }
      }
      if (verbose && soil_subset_name != "global") message("Downloading continental subset: ", paste0(toupper(substr(soil_subset_name, 1, 1)), substring(soil_subset_name, 2)))
      if (soil_subset_name != "global") {
        soil_subset_dir <- file.path(soil_dir, soil_subset_name)
        if (!dir.exists(soil_subset_dir)) dir.create(soil_subset_dir, recursive = TRUE)
        soil_zip_file <- file.path(soil_dir, paste0("soil_", soil_subset_name, ".zip"))
        output_files <- file.path(soil_subset_dir, paste0("soil_", soil_source_variables, "_", soil_subset_name, ".tif"))
        valid_files <- sapply(output_files, function(f) file.exists(f) && file.size(f) > 1e5)
        all_valid <- all(valid_files)
        if (all_valid && !redownload.rasters && file.exists(soil_zip_file)) file.remove(soil_zip_file)
        if (all_valid && !redownload.rasters) {
          if (verbose) message("Soil subset rasters already present - skipping download")
        } else {
          if (!file.exists(soil_zip_file) || file.size(soil_zip_file) < soil_zip_min_size_mb * 1e6 || redownload.rasters) {
            robust.download.raster(soil_zip_url, soil_zip_file, min_size_mb = soil_zip_min_size_mb)
          } else if (verbose) {
            message("Soil subset zip already present: ", soil_subset_name)
          }
          utils::unzip(soil_zip_file, exdir = soil_subset_dir, overwrite = TRUE)
          if (any(!file.exists(output_files))) stop("Not all soil subset raster files were found after unzip for subset: ", soil_subset_name)
          if (file.exists(soil_zip_file)) file.remove(soil_zip_file)
        }
      }
      soil_raster_list <- list()
      for (i in seq_along(soil_source_variables)) {
        output_file <- output_files[i]
        soil_raster_list[[soil_variable_names[i]]] <- terra::rast(output_file)
      }
      soil_raster_stack <- terra::rast(soil_raster_list)
      names(soil_raster_stack) <- soil_variable_names
      raster_crs <- terra::crs(soil_raster_stack)
      coord_env <- coordinate_vector_env
      coord_bg <- coordinate_vector_bg
      if (!terra::same.crs(coordinate_vector_env, soil_raster_stack)) {
        if (verbose) message("Projecting coordinates to match Soil raster CRS")
        coord_env <- terra::project(coordinate_vector_env, raster_crs)
        if (!is.null(coordinate_vector_bg)) coord_bg <- terra::project(coordinate_vector_bg, raster_crs)
      }
      extracted_occurrences <- as.matrix(terra::extract(soil_raster_stack, coord_env, ID = FALSE))
      colnames(extracted_occurrences) <- soil_variable_names
      environmental_dataset <- cbind(environmental_dataset, extracted_occurrences)
      write.csv(environmental_dataset[, soil_variable_names, drop = FALSE], soil_occ_csv_file)
      if (!is.null(background.data) && generate.background.data) {
        extracted_background <- as.matrix(terra::extract(soil_raster_stack, coord_bg, ID = FALSE))
        colnames(extracted_background) <- soil_variable_names
        background.data <- cbind(background.data, extracted_background)
        write.csv(background.data[, soil_variable_names, drop = FALSE], soil_bg_csv_file)
      }
      if (any(!soil_variable_names %in% names(environmental_dataset))) stop("Soil columns missing from occurrence data after extraction - extraction seems to not have worked")
      if (generate.background.data && any(!soil_variable_names %in% names(background.data))) stop("Soil columns missing from background data after extraction - extraction seems to not have worked")
      suppressWarnings(try(terra::tmpFiles(remove = TRUE), silent = TRUE))
      if (exists("soil_raster_stack")) rm(soil_raster_stack)
      if (exists("extracted_occurrences")) rm(extracted_occurrences)
      if (exists("extracted_background")) rm(extracted_background)
      invisible(gc())
      if (all(is.na(environmental_dataset[, soil_variable_names]))) warning("All Soil values in occurrence data are NA - check CRS or study extent")
      if (generate.background.data && all(is.na(background.data[, soil_variable_names]))) warning("All Soil values in background data are NA - check CRS or study extent")
      if (delete.intermediate.files.folders && "soil" %in% datasets_requested) unlink(file.path(rasters.dir, "soil"), recursive = TRUE, force = TRUE)
    }
  }

  # Download and process forest height data (ETH Global Canopy Height 2020 derived rasters; 250 m regional/realm subsets and 500 m global raster; downloaded from Zenodo mirror: https://zenodo.org/records/19686625)
  if ("forest_height" %in% datasets_requested) {
    invisible(gc())
    counter <- counter + 1
    if (verbose) message("")
    if (verbose) message(sprintf("--- Downloading and extracting forest height data (ca. 1-20 min depending on subset): env.dataset %d of %d ---",
                                 counter, total_datasets))
    forest_height_variable_name <- "Forestheight"
    forest_height_dir <- file.path(rasters.dir, "forest_height")
    if (!dir.exists(forest_height_dir)) dir.create(forest_height_dir, recursive = TRUE)
    forest_height_occ_csv_file <- file.path(intermediate_files_dir, "forest_height_extracted_occurrence.csv")
    forest_height_bg_csv_file <- file.path(intermediate_files_dir, "forest_height_extracted_background.csv")
    forest_height_occ_exists <- file.exists(forest_height_occ_csv_file)
    forest_height_bg_exists <- file.exists(forest_height_bg_csv_file)
    if (forest_height_occ_exists && (!generate.background.data || forest_height_bg_exists) && !overwrite) {
      if (verbose) message("Forest height data already exist - skipping download and extraction - loading from intermediate files")
      occ_forest_height <- read.csv(forest_height_occ_csv_file, row.names = 1, check.names = FALSE)
      occ_forest_height <- occ_forest_height[base_ids, , drop = FALSE]
      environmental_dataset <- cbind(environmental_dataset, occ_forest_height)
      if (generate.background.data && forest_height_bg_exists) {
        bg_forest_height <- read.csv(forest_height_bg_csv_file, row.names = 1, check.names = FALSE)
        if (!is.null(background.data)) background.data <- cbind(background.data, bg_forest_height)
      }
      if (all(is.na(environmental_dataset[, setdiff(names(occ_forest_height), "ID")]))) warning("All forest height values in occurrence data are NA - check CRS or extent")
      if (generate.background.data && all(is.na(background.data[, setdiff(names(bg_forest_height), "ID")]))) warning("All forest height values in background data are NA - check CRS or extent")
    } else {
      forest_height_record_id <- "19686625"
      forest_height_subsets <- data.frame(
        subset = c("europe",
                   "asia",
                   "north_america",
                   "south_america",
                   "africa",
                   "australia",
                   "new_world",
                   "indo_pacific",
                   "eurasia",
                   "holarctic",
                   "old_world",
                   "global"),
        file = c("EUR_canopyheight_250m.tif",
                 "ASIA_canopyheight_250m.tif",
                 "NAM_canopyheight_250m.tif",
                 "SAM_canopyheight_250m.tif",
                 "AFR_canopyheight_250m.tif",
                 "AUS_canopyheight_250m.tif",
                 "NEW_WORLD_canopyheight_250m.tif",
                 "INDO_PACIFIC_canopyheight_250m.tif",
                 "EURASIA_canopyheight_250m.tif",
                 "HOLARCTIC_canopyheight_250m.tif",
                 "OLD_WORLD_canopyheight_250m.tif",
                 "WORLD_canopyheight_500m.tif"),
        min_size_mb = c(450,
                        2200,
                        900,
                        900,
                        800,
                        200,
                        1800,
                        2200,
                        2400,
                        3500,
                        3500,
                        1300),
        stringsAsFactors = FALSE
      )
      forest_height_subsets$url <- paste0("https://zenodo.org/records/",
                                          forest_height_record_id,
                                          "/files/",
                                          forest_height_subsets$file,
                                          "?download=1")
      is_europe_forest_height <- terra::relate(study_area_vect, europe_extent_general, relation = "within")[1]
      is_asia_forest_height <- terra::relate(study_area_vect, asia_extent_general, relation = "within")[1]
      is_eurasia_forest_height <- terra::relate(study_area_vect, eurasia_extent_general, relation = "within")[1]
      is_north_america_forest_height <- terra::relate(study_area_vect, northamerica_extent_general, relation = "within")[1]
      is_south_america_forest_height <- terra::relate(study_area_vect, southamerica_extent_general, relation = "within")[1]
      is_africa_forest_height <- terra::relate(study_area_vect, africa_extent_general, relation = "within")[1]
      is_australia_forest_height <- terra::relate(study_area_vect, australia_extent_general, relation = "within")[1]
      is_indo_pacific_forest_height <- terra::relate(study_area_vect, indopacific_extent_general, relation = "within")[1]
      is_holarctic_forest_height <- terra::relate(study_area_vect, holarctic_extent_general, relation = "within")[1]
      is_new_world_forest_height <- terra::relate(study_area_vect, newworld_extent_general, relation = "within")[1]
      is_old_world_forest_height <- terra::relate(study_area_vect, oldworld_extent_general, relation = "within")[1]
      forest_height_subset_name <- if (is_europe_forest_height) {
        "europe"
      } else if (is_asia_forest_height) {
        "asia"
      } else if (is_north_america_forest_height) {
        "north_america"
      } else if (is_south_america_forest_height) {
        "south_america"
      } else if (is_africa_forest_height) {
        "africa"
      } else if (is_australia_forest_height) {
        "australia"
      } else if (is_new_world_forest_height) {
        "new_world"
      } else if (is_indo_pacific_forest_height) {
        "indo_pacific"
      } else if (is_eurasia_forest_height) {
        "eurasia"
      } else if (is_holarctic_forest_height) {
        "holarctic"
      } else if (is_old_world_forest_height) {
        "old_world"
      } else {
        "global"
      }
      subset_row <- forest_height_subsets[forest_height_subsets$subset == forest_height_subset_name, , drop = FALSE]
      forest_height_tif_file <- file.path(forest_height_dir, subset_row$file)
      forest_height_url <- subset_row$url
      forest_height_min_size_mb <- subset_row$min_size_mb
      if (verbose) message("Using forest height subset: ", forest_height_subset_name)
      if (!file.exists(forest_height_tif_file) || file.size(forest_height_tif_file) < forest_height_min_size_mb * 1e6 || redownload.rasters) {
        robust.download.raster(forest_height_url,
                               forest_height_tif_file,
                               min_size_mb = forest_height_min_size_mb)
      } else {
        if (verbose) message("Forest height raster already present - skipping download")
      }
      forest_height_raster <- terra::rast(forest_height_tif_file)
      names(forest_height_raster) <- forest_height_variable_name
      raster_crs <- terra::crs(forest_height_raster)
      coord_env <- coordinate_vector_env
      coord_bg <- coordinate_vector_bg
      if (!terra::same.crs(raster_crs, CRS_all)) {
        if (verbose) message("Projecting coordinates to match forest height raster CRS")
        coord_env <- terra::project(coordinate_vector_env, raster_crs)
        if (!is.null(coordinate_vector_bg)) coord_bg <- terra::project(coordinate_vector_bg, raster_crs)
      }
      result <- extract.and.cache.env.dataset(dataset_name = "forest_height",
                                              raster_object = forest_height_raster,
                                              coord_env = coord_env,
                                              coord_bg = coord_bg,
                                              environmental_dataset = environmental_dataset,
                                              background.data = background.data,
                                              output.dir = intermediate_files_dir,
                                              overwrite = overwrite,
                                              generate.background.data = generate.background.data,
                                              verbose = verbose)
      environmental_dataset <- result$environmental_dataset
      background.data <- result$background.data
      if (any(!forest_height_variable_name %in% names(environmental_dataset))) stop("Forest height column missing from occurrence data after extraction - extraction seems to not have worked")
      if (generate.background.data && any(!forest_height_variable_name %in% names(background.data))) stop("Forest height column missing from background data after extraction - extraction seems to not have worked")
      if (all(is.na(environmental_dataset[, forest_height_variable_name, drop = TRUE]))) warning("All forest height values in occurrence data are NA - check CRS or study extent")
      if (generate.background.data && all(is.na(background.data[, forest_height_variable_name, drop = TRUE]))) warning("All forest height values in background data are NA - check CRS or study extent")
      if (exists("forest_height_raster")) rm(forest_height_raster)
      if (exists("result")) rm(result)
      invisible(gc())
      if (delete.intermediate.files.folders && "forest_height" %in% datasets_requested) unlink(file.path(rasters.dir, "forest_height"), recursive = TRUE, force = TRUE)
    }
  }

  # Download and process atmosphere variables (WorldClim 2.1; 1km resolution; bimonthly composites: srad, wind, vapr; downloaded from Zenodo mirror: https://zenodo.org/records/17495257)
  if ("atmosphere" %in% datasets_requested) {
    invisible(gc())
    counter <- counter + 1
    if (verbose) message("")
    if (verbose) message(sprintf("-- Downloading and extracting atmosphere data (ca. 1 min): env.dataset %d of %d --",
                                 counter, total_datasets))
    atmosphere_variable_names <- c("srad_median", "srad_min", "srad_max",
                                   "wind_median", "wind_min", "wind_max",
                                   "vapr_median", "vapr_min", "vapr_max")
    atmosphere_dir <- file.path(rasters.dir, "atmosphere")
    if (!dir.exists(atmosphere_dir)) dir.create(atmosphere_dir, recursive = TRUE)
    atmosphere_occ_csv_file <- file.path(intermediate_files_dir, "Atmosphere_extracted_occurrence.csv")
    atmosphere_bg_csv_file <- file.path(intermediate_files_dir, "Atmosphere_extracted_background.csv")
    atmosphere_occ_exists <- file.exists(atmosphere_occ_csv_file)
    atmosphere_bg_exists <- file.exists(atmosphere_bg_csv_file)
    if (atmosphere_occ_exists && file.info(atmosphere_occ_csv_file)$size > 0 && (!generate.background.data || (atmosphere_bg_exists && file.info(atmosphere_bg_csv_file)$size > 0)) && !overwrite) {
      if (verbose) message("Atmosphere data already exist - skipping download and extraction")
      extracted_occurrences <- read.csv(atmosphere_occ_csv_file, row.names = 1, check.names = FALSE)
      missing_atmosphere_occ_cols <- setdiff(atmosphere_variable_names, names(extracted_occurrences))
      if (length(missing_atmosphere_occ_cols) > 0) stop("Cached atmosphere occurrence file is missing columns: ", paste(missing_atmosphere_occ_cols, collapse = ", "))
      extracted_occurrences <- extracted_occurrences[base_ids, atmosphere_variable_names, drop = FALSE]
      environmental_dataset <- cbind(environmental_dataset, extracted_occurrences)
      if (generate.background.data) {
        extracted_background <- read.csv(atmosphere_bg_csv_file, row.names = 1, check.names = FALSE)
        missing_atmosphere_bg_cols <- setdiff(atmosphere_variable_names, names(extracted_background))
        if (length(missing_atmosphere_bg_cols) > 0) stop("Cached atmosphere background file is missing columns: ", paste(missing_atmosphere_bg_cols, collapse = ", "))
        if (!is.null(background.data)) background.data <- cbind(background.data, extracted_background[rownames(background.data), atmosphere_variable_names, drop = FALSE])
      }
    }
    if (!atmosphere_occ_exists || (generate.background.data && (!atmosphere_bg_exists || file.info(atmosphere_bg_csv_file)$size == 0)) || overwrite) {
      africa_extent_atmosphere <- terra::ext(-20, 60, -35, 38)
      asia_extent_atmosphere <- terra::ext(35, 180, -10, 80)
      europe_extent_atmosphere <- terra::ext(-25, 45, 34, 72)
      eurasia_extent_atmosphere <- terra::ext(-25, 180, -10, 80)
      north_america_extent_atmosphere <- terra::ext(-170, -30, 5, 85)
      south_america_extent_atmosphere <- terra::ext(-90, -30, -60, 15)
      australia_extent_atmosphere <- terra::ext(110, 180, -50, 0)
      old_world_extent_atmosphere <- terra::ext(-25, 180, -50, 80)
      new_world_extent_atmosphere <- terra::ext(-170, -30, -60, 85)
      if (terra::relate(study_area_vect, europe_extent_atmosphere, relation = "within")[1]) {
        continent_name_atmosphere <- "Europe"
      } else if (terra::relate(study_area_vect, asia_extent_atmosphere, relation = "within")[1]) {
        continent_name_atmosphere <- "Asia"
      } else if (terra::relate(study_area_vect, eurasia_extent_atmosphere, relation = "within")[1]) {
        continent_name_atmosphere <- "Eurasia"
      } else if (terra::relate(study_area_vect, north_america_extent_atmosphere, relation = "within")[1]) {
        continent_name_atmosphere <- "NorthAmerica"
      } else if (terra::relate(study_area_vect, south_america_extent_atmosphere, relation = "within")[1]) {
        continent_name_atmosphere <- "SouthAmerica"
      } else if (terra::relate(study_area_vect, africa_extent_atmosphere, relation = "within")[1]) {
        continent_name_atmosphere <- "Africa"
      } else if (terra::relate(study_area_vect, australia_extent_atmosphere, relation = "within")[1]) {
        continent_name_atmosphere <- "Australia"
      } else if (terra::relate(study_area_vect, old_world_extent_atmosphere, relation = "within")[1]) {
        continent_name_atmosphere <- "OldWorld"
      } else if (terra::relate(study_area_vect, new_world_extent_atmosphere, relation = "within")[1]) {
        continent_name_atmosphere <- "NewWorld"
      } else {
        continent_name_atmosphere <- "Global"
      }
      if (verbose) message(sprintf("Using %s atmosphere rasters (1km resolution)", continent_name_atmosphere))
      base_url_atmosphere <- "https://zenodo.org/records/17495257/files"
      atmosphere_variables <- c("srad", "wind", "vapr")
      atmosphere_raster_files <- file.path(atmosphere_dir, sprintf("worldclim_%s_bimonthly_0.5arcmin_%s.tif", atmosphere_variables, continent_name_atmosphere))
      atmosphere_raster_urls <- sprintf("%s/worldclim_%s_bimonthly_0.5arcmin_%s.tif?download=1", base_url_atmosphere, atmosphere_variables, continent_name_atmosphere)
      expected_size_mb <- list(
        Africa = c(600, 200, 230),
        Asia = c(1300, 490, 650),
        Australia = c(104, 50, 70),
        Europe = c(200, 80, 100),
        Eurasia = c(1900, 760, 1000),
        Global = c(4200, 1900, 2200),
        NewWorld = c(1000, 420, 550),
        NorthAmerica = c(700, 280, 330),
        OldWorld = c(2200, 900, 1200),
        SouthAmerica = c(280, 70, 150)
      )[[continent_name_atmosphere]]
      for (i in seq_along(atmosphere_variables)) {
        min_size_mb <- if (!is.null(expected_size_mb)) expected_size_mb[i] else 50
        if (!file.exists(atmosphere_raster_files[i]) || file.size(atmosphere_raster_files[i]) < (min_size_mb * 1e6) || redownload.rasters) robust.download.raster(atmosphere_raster_urls[i], atmosphere_raster_files[i], min_size_mb = min_size_mb)
      }
      srad_atmosphere_raster <- terra::rast(atmosphere_raster_files[1])
      wind_atmosphere_raster <- terra::rast(atmosphere_raster_files[2])
      vapr_atmosphere_raster <- terra::rast(atmosphere_raster_files[3])
      names(srad_atmosphere_raster) <- sprintf("srad_%02d", 1:6)
      names(wind_atmosphere_raster) <- sprintf("wind_%02d", 1:6)
      names(vapr_atmosphere_raster) <- sprintf("vapr_%02d", 1:6)
      atmosphere_stack <- c(srad_atmosphere_raster, wind_atmosphere_raster, vapr_atmosphere_raster)
      raster_crs_atmosphere <- terra::crs(atmosphere_stack)
      coord_env <- coordinate_vector_env
      coord_bg <- coordinate_vector_bg
      if (!terra::same.crs(coordinate_vector_env, atmosphere_stack)) {
        if (verbose) message("Projecting coordinates to match atmosphere raster CRS")
        coord_env <- terra::project(coordinate_vector_env, raster_crs_atmosphere)
        if (!is.null(coordinate_vector_bg)) coord_bg <- terra::project(coordinate_vector_bg, raster_crs_atmosphere)
      }
      extracted_occurrences_atmosphere <- as.matrix(terra::extract(atmosphere_stack, coord_env, ID = FALSE))
      colnames(extracted_occurrences_atmosphere) <- names(atmosphere_stack)
      srad_idx_atmosphere <- grep("^srad_", colnames(extracted_occurrences_atmosphere))
      wind_idx_atmosphere <- grep("^wind_", colnames(extracted_occurrences_atmosphere))
      vapr_idx_atmosphere <- grep("^vapr_", colnames(extracted_occurrences_atmosphere))
      if (!length(srad_idx_atmosphere) || !length(wind_idx_atmosphere) || !length(vapr_idx_atmosphere)) stop("Atmosphere raster stack is missing srad, wind, or vapr layers")
      occ_atmosphere_summary <- data.frame(srad_median = apply(extracted_occurrences_atmosphere[, srad_idx_atmosphere, drop = FALSE], 1, function(x) ifelse(all(is.na(x)), NA, median(x, na.rm = TRUE))),
                                           srad_min = apply(extracted_occurrences_atmosphere[, srad_idx_atmosphere, drop = FALSE], 1, function(x) ifelse(all(is.na(x)), NA, min(x, na.rm = TRUE))),
                                           srad_max = apply(extracted_occurrences_atmosphere[, srad_idx_atmosphere, drop = FALSE], 1, function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE))),
                                           wind_median = apply(extracted_occurrences_atmosphere[, wind_idx_atmosphere, drop = FALSE], 1, function(x) ifelse(all(is.na(x)), NA, median(x, na.rm = TRUE))),
                                           wind_min = apply(extracted_occurrences_atmosphere[, wind_idx_atmosphere, drop = FALSE], 1, function(x) ifelse(all(is.na(x)), NA, min(x, na.rm = TRUE))),
                                           wind_max = apply(extracted_occurrences_atmosphere[, wind_idx_atmosphere, drop = FALSE], 1, function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE))),
                                           vapr_median = apply(extracted_occurrences_atmosphere[, vapr_idx_atmosphere, drop = FALSE], 1, function(x) ifelse(all(is.na(x)), NA, median(x, na.rm = TRUE))),
                                           vapr_min = apply(extracted_occurrences_atmosphere[, vapr_idx_atmosphere, drop = FALSE], 1, function(x) ifelse(all(is.na(x)), NA, min(x, na.rm = TRUE))),
                                           vapr_max = apply(extracted_occurrences_atmosphere[, vapr_idx_atmosphere, drop = FALSE], 1, function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE))),
                                           check.names = FALSE)
      rownames(occ_atmosphere_summary) <- rownames(environmental_dataset)
      occ_atmosphere_summary <- occ_atmosphere_summary[, atmosphere_variable_names, drop = FALSE]
      environmental_dataset <- cbind(environmental_dataset, occ_atmosphere_summary)
      write.csv(occ_atmosphere_summary, atmosphere_occ_csv_file, row.names = TRUE)
      if (!is.null(background.data) && generate.background.data) {
        extracted_background_atmosphere <- as.matrix(terra::extract(atmosphere_stack, coord_bg, ID = FALSE))
        colnames(extracted_background_atmosphere) <- names(atmosphere_stack)
        bg_atmosphere_summary <- data.frame(srad_median = apply(extracted_background_atmosphere[, srad_idx_atmosphere, drop = FALSE], 1, function(x) ifelse(all(is.na(x)), NA, median(x, na.rm = TRUE))),
                                            srad_min = apply(extracted_background_atmosphere[, srad_idx_atmosphere, drop = FALSE], 1, function(x) ifelse(all(is.na(x)), NA, min(x, na.rm = TRUE))),
                                            srad_max = apply(extracted_background_atmosphere[, srad_idx_atmosphere, drop = FALSE], 1, function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE))),
                                            wind_median = apply(extracted_background_atmosphere[, wind_idx_atmosphere, drop = FALSE], 1, function(x) ifelse(all(is.na(x)), NA, median(x, na.rm = TRUE))),
                                            wind_min = apply(extracted_background_atmosphere[, wind_idx_atmosphere, drop = FALSE], 1, function(x) ifelse(all(is.na(x)), NA, min(x, na.rm = TRUE))),
                                            wind_max = apply(extracted_background_atmosphere[, wind_idx_atmosphere, drop = FALSE], 1, function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE))),
                                            vapr_median = apply(extracted_background_atmosphere[, vapr_idx_atmosphere, drop = FALSE], 1, function(x) ifelse(all(is.na(x)), NA, median(x, na.rm = TRUE))),
                                            vapr_min = apply(extracted_background_atmosphere[, vapr_idx_atmosphere, drop = FALSE], 1, function(x) ifelse(all(is.na(x)), NA, min(x, na.rm = TRUE))),
                                            vapr_max = apply(extracted_background_atmosphere[, vapr_idx_atmosphere, drop = FALSE], 1, function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE))),
                                            check.names = FALSE)
        rownames(bg_atmosphere_summary) <- rownames(background.data)
        bg_atmosphere_summary <- bg_atmosphere_summary[, atmosphere_variable_names, drop = FALSE]
        background.data <- cbind(background.data, bg_atmosphere_summary)
        write.csv(bg_atmosphere_summary, atmosphere_bg_csv_file, row.names = TRUE)
      }
      if (any(!atmosphere_variable_names %in% names(environmental_dataset))) stop("Atmosphere columns missing from occurrence data after extraction - extraction seems to not have worked")
      if (generate.background.data && any(!atmosphere_variable_names %in% names(background.data))) stop("Atmosphere columns missing from background data after extraction - extraction seems to not have worked")
      if (all(is.na(environmental_dataset$srad_median))) warning("All atmosphere values in occurrence data are NA - check CRS or study extent")
      if (generate.background.data && all(is.na(background.data$srad_median))) warning("All atmosphere values in background data are NA - check CRS or study extent")
      if (exists("atmosphere_stack")) rm(atmosphere_stack)
      if (exists("extracted_occurrences_atmosphere")) rm(extracted_occurrences_atmosphere)
      if (exists("extracted_background_atmosphere")) rm(extracted_background_atmosphere)
      invisible(gc())
      if (delete.intermediate.files.folders && "atmosphere" %in% datasets_requested) {
        unlink(file.path(rasters.dir, "atmosphere"), recursive = TRUE, force = TRUE)
    }
    }
  }

  # Download and process DMSP Nighttime Light (NTL) data (https://eogdata.mines.edu/products/dmsp/#rad_cal; download from Zenodo mirror: https://zenodo.org/records/17416839/files/F182013.v4c.global.intercal.stable_lights.avg_vis.tif; 1km resolution)
  if ("nightlight" %in% datasets_requested) {
    invisible(gc())
    counter <- counter + 1
    if (verbose) message("")
    if (verbose) message(sprintf("-- Downloading and extracting Nighttime Light data (ca. 6 min): env.dataset %d of %d --",
                                 counter, total_datasets))
    nightlight_variable_name <- "Nighttime_light"
    nightlight_dir <- file.path(rasters.dir, "nightlight")
    if (!dir.exists(nightlight_dir)) dir.create(nightlight_dir, recursive = TRUE)
    nightlight_raster_file <- file.path(nightlight_dir, "F182013.v4c.global.intercal.stable_lights.avg_vis.tif")
    nightlight_occ_csv_file <- file.path(intermediate_files_dir, "Nightlight_extracted_occurrence.csv")
    nightlight_bg_csv_file <- file.path(intermediate_files_dir, "Nightlight_extracted_background.csv")
    nightlight_occ_exists <- file.exists(nightlight_occ_csv_file)
    nightlight_bg_exists <- file.exists(nightlight_bg_csv_file)
    if (nightlight_occ_exists && (!generate.background.data || nightlight_bg_exists) && !overwrite) {
      if (verbose) message("Nighttime Light data already exist - skipping download and extraction - loading from intermediate files")
      occ_ntl <- read.csv(nightlight_occ_csv_file, row.names = 1, check.names = FALSE)
      occ_ntl <- occ_ntl[base_ids, , drop = FALSE]
      environmental_dataset <- cbind(environmental_dataset, occ_ntl)
      if (generate.background.data && nightlight_bg_exists) {
        bg_ntl <- read.csv(nightlight_bg_csv_file, row.names = 1, check.names = FALSE)
        if (!is.null(background.data)) background.data <- cbind(background.data, bg_ntl[rownames(background.data), , drop = FALSE])
      }
      if (all(is.na(environmental_dataset[, setdiff(names(occ_ntl), "ID")])))
        warning("All Nighttime Light values in occurrence data are NA - check CRS or extent")
      if (generate.background.data && all(is.na(background.data[, setdiff(names(bg_ntl), "ID")])))
        warning("All Nighttime Light values in background data are NA - check CRS or extent")
    } else {
      nightlight_url <- "https://zenodo.org/records/17416839/files/F182013.v4c.global.intercal.stable_lights.avg_vis.tif?download=1"
      if (!file.exists(nightlight_raster_file) || file.size(nightlight_raster_file) < 8e7 || redownload.rasters) {
        robust.download.raster(nightlight_url, nightlight_raster_file, min_size_mb = 80)
      } else {
        if (verbose) message("Nighttime Light raster already present - skipping download")
      }
      nightlight_raster <- terra::rast(nightlight_raster_file)
      names(nightlight_raster) <- nightlight_variable_name
      raster_crs <- terra::crs(nightlight_raster)
      coord_env <- coordinate_vector_env
      coord_bg <- coordinate_vector_bg
      if (!terra::same.crs(coordinate_vector_env, nightlight_raster)) {
        if (verbose) message("Projecting coordinates to match Nighttime Light raster CRS")
        coord_env <- terra::project(coordinate_vector_env, raster_crs)
        if (!is.null(coordinate_vector_bg)) coord_bg <- terra::project(coordinate_vector_bg, raster_crs)
      }
      extracted_occurrences <- as.matrix(terra::extract(nightlight_raster, coord_env, ID = FALSE))
      colnames(extracted_occurrences) <- nightlight_variable_name
      environmental_dataset <- cbind(environmental_dataset, extracted_occurrences)
      write.csv(environmental_dataset[, nightlight_variable_name, drop = FALSE], nightlight_occ_csv_file)
      if (!is.null(background.data) && generate.background.data) {
        extracted_background <- as.matrix(terra::extract(nightlight_raster, coord_bg, ID = FALSE))
        colnames(extracted_background) <- nightlight_variable_name
        background.data <- cbind(background.data, extracted_background)
        write.csv(background.data[, nightlight_variable_name, drop = FALSE], nightlight_bg_csv_file)
      }
      if (any(!nightlight_variable_name %in% names(environmental_dataset))) stop("Nighttime Light column missing from occurrence data after extraction - extraction seems to not have worked")
      if (generate.background.data && any(!nightlight_variable_name %in% names(background.data))) stop("Nighttime Light column missing from background data after extraction - extraction seems to not have worked")
      suppressWarnings(try(terra::tmpFiles(remove = TRUE), silent = TRUE))
      if (exists("nightlight_raster")) rm(nightlight_raster)
      if (exists("extracted_occurrences")) rm(extracted_occurrences)
      if (exists("extracted_background")) rm(extracted_background)
      invisible(gc())
      if (all(is.na(environmental_dataset[[nightlight_variable_name]]))) warning("All Nighttime Light values in occurrence data are NA - check CRS or study extent")
      if (generate.background.data && all(is.na(background.data[[nightlight_variable_name]]))) warning("All Nighttime Light values in background data are NA - check CRS or study extent")
      if (delete.intermediate.files.folders && "nightlight" %in% datasets_requested) unlink(file.path(rasters.dir, "nightlight"), recursive = TRUE, force = TRUE)
    }
  }

    # Download and process Burned area data (monthly means across 2019-2024; ESA FireCCI Sentinel-3 SYN v1.1; 25km resolution; download via zenodo mirror)
    if ("burned_area" %in% datasets_requested) {
      invisible(gc())
      counter <- counter + 1
      if (verbose) message("")
      if (verbose) message(sprintf("-- Downloading and extracting Burned area data (ca. 2 min): env.dataset %d of %d --",
                                   counter, total_datasets))
      burned_area_variable_name <- "Burned_area"
      burned_area_dir <- file.path(rasters.dir, "burned_area")
      if (!dir.exists(burned_area_dir)) dir.create(burned_area_dir, recursive = TRUE)
      burned_area_raster_file <- file.path(burned_area_dir, "FireCCI_BurnedArea_MonthlyMean_Global_0.25deg.zip")
      burned_area_occ_csv_file <- file.path(intermediate_files_dir, "burned_area_extracted_occurrence.csv")
      burned_area_bg_csv_file <- file.path(intermediate_files_dir, "burned_area_extracted_background.csv")
      burned_area_occ_exists <- file.exists(burned_area_occ_csv_file)
      burned_area_bg_exists <- file.exists(burned_area_bg_csv_file)
      burned_area_existing_tifs <- list.files(burned_area_dir, pattern = "^FireCCI_BurnedArea_Mean_[0-9]{2}\\.tif$", full.names = TRUE)
      burned_area_existing_tifs <- burned_area_existing_tifs[file.info(burned_area_existing_tifs)$size > 180000]
      if (length(burned_area_existing_tifs) == 12 && file.exists(burned_area_raster_file)) file.remove(burned_area_raster_file)
      if (burned_area_occ_exists && (!generate.background.data || burned_area_bg_exists) && !overwrite) {
        if (verbose) message("Burned area data already exist - skipping download and extraction - loading from intermediate files")
        occ_burn <- read.csv(burned_area_occ_csv_file, row.names = 1, check.names = FALSE)
        occ_burn <- occ_burn[base_ids, , drop = FALSE]
        environmental_dataset <- cbind(environmental_dataset, occ_burn)
        if (generate.background.data && burned_area_bg_exists) {
          bg_burn <- read.csv(burned_area_bg_csv_file, row.names = 1, check.names = FALSE)
          if (!is.null(background.data)) background.data <- cbind(background.data, bg_burn[rownames(background.data), , drop = FALSE])
        }
        if (all(sapply(sprintf("Burned_area_%02d", 1:12), function(nm) all(is.na(environmental_dataset[[nm]]))))) warning("All Burned area values in occurrence data are NA - check CRS or extent")
        if (generate.background.data && all(sapply(sprintf("Burned_area_%02d", 1:12), function(nm) all(is.na(background.data[[nm]]))))) warning("All Burned area values in background data are NA - check CRS or extent")
      } else {
        burned_area_url <- "https://zenodo.org/records/17487469/files/FireCCI_BurnedArea_MonthlyMean_Global_0.25deg.zip?download=1"
        if (length(burned_area_existing_tifs) == 12 && !redownload.rasters) {
          burned_tifs <- burned_area_existing_tifs
        } else {
          if (!file.exists(burned_area_raster_file) || file.size(burned_area_raster_file) < 2.6e6 || redownload.rasters) {
            robust.download.raster(burned_area_url, burned_area_raster_file, min_size_mb = 2.6)
          } else {
            if (verbose) message("Burned area raster zip already present - skipping download")
          }
          utils::unzip(burned_area_raster_file, exdir = burned_area_dir, overwrite = TRUE)
          burned_tifs <- list.files(burned_area_dir, pattern = "^FireCCI_BurnedArea_Mean_[0-9]{2}\\.tif$", full.names = TRUE)
          burned_tifs <- burned_tifs[file.info(burned_tifs)$size > 180000]
          if (length(burned_tifs) != 12) {
            robust.download.raster(burned_area_url, burned_area_raster_file, min_size_mb = 2.6)
            utils::unzip(burned_area_raster_file, exdir = burned_area_dir, overwrite = TRUE)
            burned_tifs <- list.files(burned_area_dir, pattern = "^FireCCI_BurnedArea_Mean_[0-9]{2}\\.tif$", full.names = TRUE)
            burned_tifs <- burned_tifs[file.info(burned_tifs)$size > 180000]
          }
        }
        if (length(burned_tifs) != 12) stop("Not all burned-area raster files were found after unzip")
        if (file.exists(burned_area_raster_file)) file.remove(burned_area_raster_file)
        burned_area_raster <- terra::rast(burned_tifs)
        names(burned_area_raster) <- sprintf("Burned_area_%02d", 1:12)
        raster_crs <- terra::crs(burned_area_raster)
        coord_env <- coordinate_vector_env
        coord_bg <- coordinate_vector_bg
        if (!terra::same.crs(coordinate_vector_env, burned_area_raster)) {
          if (verbose) message("Projecting coordinates to match Burned area raster CRS")
          coord_env <- terra::project(coordinate_vector_env, raster_crs)
          if (!is.null(coordinate_vector_bg)) coord_bg <- terra::project(coordinate_vector_bg, raster_crs)
        }
        extracted_occurrences <- as.matrix(terra::extract(burned_area_raster, coord_env, ID = FALSE))
        colnames(extracted_occurrences) <- names(burned_area_raster)
        rownames(extracted_occurrences) <- rownames(environmental_dataset)
        environmental_dataset <- cbind(environmental_dataset, extracted_occurrences)
        write.csv(environmental_dataset[, names(burned_area_raster), drop = FALSE], burned_area_occ_csv_file)
        if (!is.null(background.data) && generate.background.data) {
          extracted_background <- as.matrix(terra::extract(burned_area_raster, coord_bg, ID = FALSE))
          colnames(extracted_background) <- names(burned_area_raster)
          rownames(extracted_background) <- rownames(background.data)
          background.data <- cbind(background.data, extracted_background)
          write.csv(background.data[, names(burned_area_raster), drop = FALSE], burned_area_bg_csv_file)
        }
        if (any(!names(burned_area_raster) %in% names(environmental_dataset))) stop("Burned area column missing from occurrence data after extraction - extraction seems to not have worked")
        if (generate.background.data && any(!names(burned_area_raster) %in% names(background.data))) stop("Burned area column missing from background data after extraction - extraction seems to not have worked")
        if (all(sapply(sprintf("Burned_area_%02d", 1:12), function(nm) all(is.na(environmental_dataset[[nm]]))))) warning("All Burned area values in occurrence data are NA - check CRS or extent")
        if (generate.background.data && all(sapply(sprintf("Burned_area_%02d", 1:12), function(nm) all(is.na(background.data[[nm]]))))) warning("All Burned area values in background data are NA - check CRS or extent")
        suppressWarnings(try(terra::tmpFiles(remove = TRUE), silent = TRUE))
        if (exists("burned_area_raster")) rm(burned_area_raster)
        if (exists("extracted_occurrences")) rm(extracted_occurrences)
        if (exists("extracted_background")) rm(extracted_background)
        invisible(gc())
        if (delete.intermediate.files.folders && "burned_area" %in% datasets_requested) unlink(file.path(rasters.dir, "burned_area"), recursive = TRUE, force = TRUE)
      }
    }

  # Download and process Snow Water Equivalent (SWE) data (monthly means 2023, only for North America; 1 km resolution; download via zenodo mirror; derived from Daymet: https://daymet.ornl.gov/; Daymet Lambert Conformal Conic)
  if ("snow_water_equivalent" %in% datasets_requested && is_north_america_SWE) {
    invisible(gc())
    counter <- counter + 1
    if (verbose) message("")
    if (verbose) message(sprintf("-- Downloading and extracting SWE data (ca. 12 min): env.dataset %d of %d --",
                                 counter, total_datasets))
    swe_variable_names <- sprintf("snow_water_equivalent_%02d", 1:12)
    swe_dir <- file.path(rasters.dir, "snow_water_equivalent")
    if (!dir.exists(swe_dir)) dir.create(swe_dir, recursive = TRUE)
    swe_raster_file <- file.path(swe_dir, "daymet_v4_swe_monavg_NorthAmerica_combined_2023.tif")
    swe_occ_csv_file <- file.path(intermediate_files_dir, "SWE_extracted_occurrence.csv")
    swe_bg_csv_file <- file.path(intermediate_files_dir, "SWE_extracted_background.csv")
    swe_occ_exists <- file.exists(swe_occ_csv_file)
    swe_bg_exists <- file.exists(swe_bg_csv_file)
    if (swe_occ_exists && (!generate.background.data || swe_bg_exists) && !overwrite) {
      if (verbose) message("Daymet SWE data already exist - skipping download and extraction - loading from intermediate files")
      occ_swe <- read.csv(swe_occ_csv_file, row.names = 1, check.names = FALSE)
      occ_swe <- occ_swe[base_ids, , drop = FALSE]
      environmental_dataset <- cbind(environmental_dataset, occ_swe)
      if (generate.background.data && swe_bg_exists) {
        bg_swe <- read.csv(swe_bg_csv_file, row.names = 1, check.names = FALSE)
        if (!is.null(background.data)) background.data <- cbind(background.data, bg_swe[rownames(background.data), , drop = FALSE])
      }
      if (all(is.na(environmental_dataset[, setdiff(names(occ_swe), "ID")]))) warning("All SWE values in occurrence data are NA - check CRS or extent")
      if (generate.background.data && all(is.na(background.data[, setdiff(names(bg_swe), "ID")]))) warning("All SWE values in background data are NA - check CRS or extent")
    } else {
      swe_url <- "https://zenodo.org/records/17495170/files/daymet_v4_swe_monavg_NorthAmerica_combined_2023.tif?download=1"
      if (!file.exists(swe_raster_file) || file.size(swe_raster_file) < 5e8 || redownload.rasters) {
        robust.download.raster(swe_url, swe_raster_file, min_size_mb = 500)
      } else {
        if (verbose) message("Daymet SWE raster already present - skipping download")
      }
      swe_raster <- terra::rast(swe_raster_file)
      names(swe_raster) <- swe_variable_names
      raster_crs <- terra::crs(swe_raster)
      coord_env <- coordinate_vector_env
      coord_bg <- coordinate_vector_bg
      if (!terra::same.crs(coordinate_vector_env, swe_raster)) {
        if (verbose) message("Projecting coordinates to match Daymet SWE raster CRS")
        coord_env <- terra::project(coordinate_vector_env, raster_crs)
        if (!is.null(coordinate_vector_bg)) coord_bg <- terra::project(coordinate_vector_bg, raster_crs)
      }
      extracted_occurrences <- as.matrix(terra::extract(swe_raster, coord_env, ID = FALSE))
      colnames(extracted_occurrences) <- names(swe_raster)
      environmental_dataset <- cbind(environmental_dataset, extracted_occurrences)
      write.csv(environmental_dataset[, names(swe_raster), drop = FALSE], swe_occ_csv_file)
      if (!is.null(background.data) && generate.background.data) {
        extracted_background <- as.matrix(terra::extract(swe_raster, coord_bg, ID = FALSE))
        colnames(extracted_background) <- names(swe_raster)
        background.data <- cbind(background.data, extracted_background)
        write.csv(background.data[, names(swe_raster), drop = FALSE], swe_bg_csv_file)
      }
      if (any(!names(swe_raster) %in% names(environmental_dataset))) stop("SWE columns missing from occurrence data - extraction failed")
      if (generate.background.data && any(!names(swe_raster) %in% names(background.data))) stop("SWE columns missing from background data - extraction failed")
      suppressWarnings(try(terra::tmpFiles(remove = TRUE), silent = TRUE))
      if (exists("swe_raster")) rm(swe_raster)
      if (exists("extracted_occurrences")) rm(extracted_occurrences)
      if (exists("extracted_background")) rm(extracted_background)
      invisible(gc())
      if (all(is.na(environmental_dataset[[swe_variable_names[1]]]))) warning("All SWE values in occurrence data are NA - check CRS or study extent")
      if (generate.background.data && all(is.na(background.data[[swe_variable_names[1]]]))) warning("All SWE values in background data are NA - check CRS or study extent")
      if (delete.intermediate.files.folders && "snow_water_equivalent" %in% datasets_requested) unlink(file.path(rasters.dir, "snow_water_equivalent"), recursive = TRUE, force = TRUE)
    }
  }

  # Download and process Daylength data (monthly means 2024, only for North America; 1 km resolution; derived from Daymet: https://daymet.ornl.gov/; download via zenodo mirror; Daymet Lambert Conformal Conic)
  if ("daylength" %in% datasets_requested && is_north_america_daylength) {
    invisible(gc())
    counter <- counter + 1
    if (verbose) message("")
    if (verbose) message(sprintf("-- Downloading and extracting Daylength data (ca. 10 min): env.dataset %d of %d --",
                                 counter, total_datasets))
    dayl_variable_names <- sprintf("Daylength_%02d", 1:12)
    dayl_dir <- file.path(rasters.dir, "daylength")
    if (!dir.exists(dayl_dir)) dir.create(dayl_dir, recursive = TRUE)
    dayl_raster_file <- file.path(dayl_dir, "Daymet_Monthly_dayl_NorthAmerica_2024.tif")
    dayl_occ_csv_file <- file.path(intermediate_files_dir, "Daylength_extracted_occurrence.csv")
    dayl_bg_csv_file <- file.path(intermediate_files_dir, "Daylength_extracted_background.csv")
    dayl_occ_exists <- file.exists(dayl_occ_csv_file)
    dayl_bg_exists <- file.exists(dayl_bg_csv_file)
    if (dayl_occ_exists && (!generate.background.data || dayl_bg_exists) && !overwrite) {
      if (verbose) message("Daymet Daylength data already exist - skipping download and extraction - loading from intermediate files")
      occ_dayl <- read.csv(dayl_occ_csv_file, row.names = 1, check.names = FALSE)
      occ_dayl <- occ_dayl[base_ids, , drop = FALSE]
      environmental_dataset <- cbind(environmental_dataset, occ_dayl)
      if (generate.background.data && dayl_bg_exists) {
        bg_dayl <- read.csv(dayl_bg_csv_file, row.names = 1, check.names = FALSE)
        if (!is.null(background.data)) background.data <- cbind(background.data, bg_dayl[rownames(background.data), , drop = FALSE])
      }
      if (all(is.na(environmental_dataset[, setdiff(names(occ_dayl), "ID")]))) warning("All Daylength values in occurrence data are NA - check CRS or extent")
      if (generate.background.data && all(is.na(background.data[, setdiff(names(bg_dayl), "ID")]))) warning("All Daylength values in background data are NA - check CRS or extent")
    } else {
      dayl_url <- "https://zenodo.org/records/17468682/files/Daymet_Monthly_dayl_NorthAmerica_2024.tif?download=1"
      if (!file.exists(dayl_raster_file) || file.size(dayl_raster_file) < 3e8 || redownload.rasters) {
        robust.download.raster(dayl_url, dayl_raster_file, min_size_mb = 300)
      } else {
        if (verbose) message("Daymet Daylength raster already present - skipping download")
      }
      dayl_raster <- terra::rast(dayl_raster_file)
      names(dayl_raster) <- dayl_variable_names
      raster_crs <- terra::crs(dayl_raster)
      coord_env <- coordinate_vector_env
      coord_bg <- coordinate_vector_bg
      if (!terra::same.crs(coordinate_vector_env, dayl_raster)) {
        if (verbose) message("Projecting coordinates to match Daymet Daylength raster CRS")
        coord_env <- terra::project(coordinate_vector_env, raster_crs)
        if (!is.null(coordinate_vector_bg)) coord_bg <- terra::project(coordinate_vector_bg, raster_crs)
      }
      extracted_occurrences <- as.matrix(terra::extract(dayl_raster, coord_env, ID = FALSE))
      colnames(extracted_occurrences) <- names(dayl_raster)
      environmental_dataset <- cbind(environmental_dataset, extracted_occurrences)
      write.csv(environmental_dataset[, names(dayl_raster), drop = FALSE], dayl_occ_csv_file)
      if (!is.null(background.data) && generate.background.data) {
        extracted_background <- as.matrix(terra::extract(dayl_raster, coord_bg, ID = FALSE))
        colnames(extracted_background) <- names(dayl_raster)
        background.data <- cbind(background.data, extracted_background)
        write.csv(background.data[, names(dayl_raster), drop = FALSE], dayl_bg_csv_file)
      }
      if (any(!names(dayl_raster) %in% names(environmental_dataset))) stop("Daylength columns missing from occurrence data - extraction failed")
      if (generate.background.data && any(!names(dayl_raster) %in% names(background.data))) stop("Daylength columns missing from background data - extraction failed")
      suppressWarnings(try(terra::tmpFiles(remove = TRUE), silent = TRUE))
      if (exists("dayl_raster")) rm(dayl_raster)
      if (exists("extracted_occurrences")) rm(extracted_occurrences)
      if (exists("extracted_background")) rm(extracted_background)
      invisible(gc())
      if (all(is.na(environmental_dataset[[dayl_variable_names[1]]]))) warning("All Daylength values in occurrence data are NA - check CRS or study extent")
      if (generate.background.data && all(is.na(background.data[[dayl_variable_names[1]]]))) warning("All Daylength values in background data are NA - check CRS or study extent")
      if (delete.intermediate.files.folders && "daylength" %in% datasets_requested) unlink(file.path(rasters.dir, "daylength"), recursive = TRUE, force = TRUE)
    }
  }

  # Download and process ESA CCI Soil Moisture (combined v09.1)
  if ("soil_moisture" %in% datasets_requested) {
    invisible(gc())
    counter <- counter + 1
    if (verbose) message("")
    if (verbose) message(sprintf("-- Downloading and extracting soil moisture data (ca. 5 min): env.dataset %d of %d --",
                                 counter, total_datasets))
    soil_moisture_variable_names <- paste0("soil_moisture_", sprintf("%02d", 1:6))
    soil_moisture_dir <- file.path(rasters.dir, "soil_moisture")
    if (!dir.exists(soil_moisture_dir)) dir.create(soil_moisture_dir, recursive = TRUE)
    soil_moisture_files <- file.path(soil_moisture_dir, paste0("ESA_soil_moisture_", sprintf("%02d", 1:6), ".tif"))
    soil_moisture_occ_csv_file <- file.path(intermediate_files_dir, "soil_moisture_extracted_occurrence.csv")
    soil_moisture_bg_csv_file <- file.path(intermediate_files_dir, "soil_moisture_extracted_background.csv")
    soil_moisture_occ_exists <- file.exists(soil_moisture_occ_csv_file)
    soil_moisture_bg_exists <- file.exists(soil_moisture_bg_csv_file)
    if (soil_moisture_occ_exists && (!generate.background.data || soil_moisture_bg_exists) && !overwrite) {
      if (verbose) message("Soil Moisture data already exist - skipping download and extraction - loading from intermediate files")
      occ_soil_moisture <- read.csv(soil_moisture_occ_csv_file, row.names = 1, check.names = FALSE)
      occ_soil_moisture <- occ_soil_moisture[base_ids, , drop = FALSE]
      environmental_dataset <- cbind(environmental_dataset, occ_soil_moisture)
      if (generate.background.data && soil_moisture_bg_exists) {
        bg_soil_moisture <- read.csv(soil_moisture_bg_csv_file, row.names = 1, check.names = FALSE)
        if (!is.null(background.data)) background.data <- cbind(background.data, bg_soil_moisture[rownames(background.data), , drop = FALSE])
      }
    } else {
      base_url <- "https://zenodo.org/records/17496280/files/"
      all_present <- all(file.exists(soil_moisture_files)) && all(file.info(soil_moisture_files)$size > 1e5)
      if (all_present && !redownload.rasters) {
        if (verbose) message("Soil Moisture rasters already present - skipping download")
      } else {
        for (i in 1:6) {
          file_url <- paste0(base_url, "ESA_SoilMoisture_", sprintf("%02d", i), ".tif?download=1")
          if (!file.exists(soil_moisture_files[i]) || file.size(soil_moisture_files[i]) < 1e5 || redownload.rasters) {
            robust.download.raster(file_url, soil_moisture_files[i], min_size_mb = 0.1)
          }
        }
      }
      soil_moisture_stack <- terra::rast(soil_moisture_files)
      names(soil_moisture_stack) <- soil_moisture_variable_names
      result_soilmoisture <- extract.and.cache.env.dataset(dataset_name = "SoilMoisture",
                                                           raster_object = soil_moisture_stack,
                                                           coord_env = coordinate_vector_env,
                                                           coord_bg = coordinate_vector_bg,
                                                           environmental_dataset = environmental_dataset,
                                                           background.data = background.data,
                                                           output.dir = intermediate_files_dir,
                                                           overwrite = overwrite,
                                                           generate.background.data = generate.background.data,
                                                           verbose = verbose)
      environmental_dataset <- result_soilmoisture$environmental_dataset
      background.data <- result_soilmoisture$background.data
      if (any(!soil_moisture_variable_names %in% names(environmental_dataset))) stop("Soil_moisture columns missing from occurrence data after extraction - extraction seems to not have worked")
      if (generate.background.data && any(!soil_moisture_variable_names %in% names(background.data))) stop("Soil_moisture columns missing from background data after extraction - extraction seems to not have worked")
      suppressWarnings(try(terra::tmpFiles(remove = TRUE), silent = TRUE))
      invisible(gc())
      if (delete.intermediate.files.folders && "soil_moisture" %in% datasets_requested) unlink(file.path(rasters.dir, "soil_moisture"), recursive = TRUE, force = TRUE)
    }
  }

  # Download and process custom environmental rasters (.tif/.tiff files, folders containing .tif/.tiff files, or .zip files containing .tif/.tiff files)
  if (!is.null(custom.env.rasters)) {
    invisible(gc())
    for (i in seq_along(custom.env.rasters)) {
      custom_source <- custom.env.rasters[[i]] #URL or local path
      is_url <- grepl("^https?://", custom_source, ignore.case = TRUE)
      is_local_tif_dir <- !is_url && dir.exists(custom_source)
      is_local_tif_zip <- !is_url && !is_local_tif_dir && grepl("\\.zip$", custom_source, ignore.case = TRUE)
      default_custom_name <- if (is_local_tif_dir) {
        basename(normalizePath(custom_source, winslash = "/", mustWork = TRUE))
      } else {
        custom_source_clean <- basename(sub("\\?.*$", "", custom_source))
        custom_source_clean <- tools::file_path_sans_ext(custom_source_clean)
        gsub("[^[:alnum:]_.-]", "_", custom_source_clean)
      }
      custom_name <- if (!is.null(custom.env.rasters.names) && length(custom.env.rasters.names) >= i) custom.env.rasters.names[[i]] else default_custom_name
      custom_var_name_provided <- !is.null(custom.env.rasters.variable.names) && length(custom.env.rasters.variable.names) >= i
      if (custom_var_name_provided) {
        custom_var_name <- custom.env.rasters.variable.names[[i]]
      } else {
        custom_var_name <- custom_name
      }
      counter <- counter + 1
      if (verbose) message("")
      if (verbose) message(sprintf("-- Downloading and extracting custom raster %s: env.dataset %d of %d --",
                                   custom_name, counter, total_datasets))
      if (is_local_tif_dir) {
        custom_dir <- normalizePath(custom_source, winslash = "/", mustWork = TRUE)
        custom_tif_files_all <- sort(list.files(custom_dir,
                                                pattern = "\\.(tif|tiff)$",
                                                full.names = TRUE,
                                                ignore.case = TRUE))
        if (length(custom_tif_files_all) == 0) stop("No .tif/.tiff files found in custom raster directory: ", custom_dir)
        custom_tif_files <- character(0)
        skipped_tif_files <- character(0)
        for (tif_file in custom_tif_files_all) {
          tif_ok <- FALSE
          tif_size <- suppressWarnings(file.info(tif_file)$size)
          if (!is.na(tif_size) && tif_size > 0) {
            tif_ok <- tryCatch({
              test_rast <- terra::rast(tif_file)
              ok <- inherits(test_rast, "SpatRaster") && terra::nlyr(test_rast) >= 1
              rm(test_rast)
              invisible(gc())
              ok
            }, error = function(e) FALSE)
          }
          if (tif_ok) {
            custom_tif_files <- c(custom_tif_files, tif_file)
          } else {
            skipped_tif_files <- c(skipped_tif_files, basename(tif_file))
            if (verbose) message("Skipping empty or unreadable .tif/.tiff file in directory ", custom_name, ": ", basename(tif_file))
          }
        }
        if (length(custom_tif_files) == 0) stop("No readable non-empty .tif/.tiff files found in custom raster directory: ", custom_dir)
        if (!custom_var_name_provided) custom_var_name <- tools::file_path_sans_ext(basename(custom_tif_files))
        if (verbose) message("Found ", length(custom_tif_files), " .tif/.tiff files in directory: ", custom_name)
        custom_tif_file <- custom_tif_files[[1]]
        custom_dir_is_temporary <- FALSE
      } else if (is_local_tif_zip) {
        custom_dir <- file.path(rasters.dir, custom_name)
        if (!dir.exists(custom_dir)) dir.create(custom_dir, recursive = TRUE)
        custom_tif_file <- file.path(custom_dir, paste0(custom_name, ".tif"))
        custom_zip_file <- file.path(custom_dir, paste0(custom_name, ".zip"))
        custom_dir_is_temporary <- TRUE
      } else if (!is_url) {
        custom_tif_file <- normalizePath(custom_source, winslash = "/", mustWork = TRUE)
        custom_dir <- dirname(custom_tif_file)
        custom_dir_is_temporary <- FALSE
      } else {
        custom_dir <- file.path(rasters.dir, custom_name)
        if (!dir.exists(custom_dir)) dir.create(custom_dir, recursive = TRUE)
        custom_tif_file <- file.path(custom_dir, paste0(custom_name, ".tif"))
        custom_zip_file <- file.path(custom_dir, paste0(custom_name, ".zip"))
        custom_dir_is_temporary <- TRUE
      }
      custom_occ_csv_file <- file.path(intermediate_files_dir, paste0(custom_name, "_extracted_occurrence.csv"))
      custom_bg_csv_file <- file.path(intermediate_files_dir, paste0(custom_name, "_extracted_background.csv"))
      custom_occ_exists <- file.exists(custom_occ_csv_file) && file.info(custom_occ_csv_file)$size > 0
      custom_bg_exists <- file.exists(custom_bg_csv_file) && file.info(custom_bg_csv_file)$size > 0
      if (custom_occ_exists && (!generate.background.data || custom_bg_exists) && !overwrite) {
        if (verbose) message("Custom raster data already exist - skipping download and extraction - loading from intermediate files: ", custom_name)
        occ_custom <- read.csv(custom_occ_csv_file, row.names = 1, check.names = FALSE)
        occ_custom <- occ_custom[base_ids, , drop = FALSE]
        environmental_dataset <- cbind(environmental_dataset, occ_custom)
        if (generate.background.data && custom_bg_exists) {
          bg_custom <- read.csv(custom_bg_csv_file, row.names = 1, check.names = FALSE)
          if (!is.null(background.data)) background.data <- cbind(background.data, bg_custom[rownames(background.data), , drop = FALSE])
        }
        if (all(is.na(environmental_dataset[, setdiff(names(occ_custom), "ID"), drop = FALSE])))
          warning("All values in occurrence data are NA for custom raster: ", custom_name, " - check CRS or extent")
        if (generate.background.data && custom_bg_exists && all(is.na(background.data[, setdiff(names(bg_custom), "ID"), drop = FALSE])))
          warning("All values in background data are NA for custom raster: ", custom_name, " - check CRS or extent")
      } else {
        if (is_local_tif_dir) {
          if (verbose) message("Using all local .tif/.tiff files in directory: ", custom_name)
          custom_raster <- terra::rast(custom_tif_files)
        } else if (is_local_tif_zip) {
          if (verbose) message("Using local zip containing .tif/.tiff files: ", custom_name)
          custom_dir <- file.path(rasters.dir, custom_name)
          if (!dir.exists(custom_dir)) dir.create(custom_dir, recursive = TRUE)
          custom_zip_file <- file.path(custom_dir, paste0(custom_name, ".zip"))
          custom_dir_is_temporary <- TRUE
          if (!identical(normalizePath(custom_source, winslash = "/", mustWork = TRUE),
                         normalizePath(custom_zip_file, winslash = "/", mustWork = FALSE))) {
            ok <- file.copy(custom_source, custom_zip_file, overwrite = TRUE)
            if (!ok || !file.exists(custom_zip_file)) stop("Failed to copy custom raster zip: ", custom_source)
          } else {
            if (verbose) message("Custom raster zip already in target location: ", custom_name)
          }
          utils::unzip(custom_zip_file, exdir = custom_dir, overwrite = TRUE)
          tif_files_all <- sort(list.files(custom_dir, pattern = "\\.(tif|tiff)$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE))
          if (length(tif_files_all) == 0) stop("No .tif/.tiff found in custom zip: ", custom_name)
          tif_files <- character(0)
          for (tif_file in tif_files_all) {
            tif_ok <- FALSE
            tif_size <- suppressWarnings(file.info(tif_file)$size)
            if (!is.na(tif_size) && tif_size > 0) {
              tif_ok <- tryCatch({
                test_rast <- terra::rast(tif_file)
                ok <- inherits(test_rast, "SpatRaster") && terra::nlyr(test_rast) >= 1
                rm(test_rast)
                invisible(gc())
                ok
              }, error = function(e) FALSE)
            }
            if (tif_ok) {
              tif_files <- c(tif_files, tif_file)
            } else {
              if (verbose) message("Skipping empty or unreadable .tif/.tiff file in zip ", custom_name, ": ", basename(tif_file))
            }
          }
          if (length(tif_files) == 0) stop("No readable non-empty .tif/.tiff files found in custom zip: ", custom_name)
          if (!custom_var_name_provided) custom_var_name <- tools::file_path_sans_ext(basename(tif_files))
          if (verbose) message("Found ", length(tif_files), " .tif/.tiff files in zip: ", custom_name)
          custom_raster <- terra::rast(tif_files)
        } else if (is_url && grepl("\\.zip($|\\?)", custom_source, ignore.case = TRUE)) {
          if (verbose) message("Downloading zip containing .tif/.tiff files: ", custom_name)
          custom_dir <- file.path(rasters.dir, custom_name)
          if (!dir.exists(custom_dir)) dir.create(custom_dir, recursive = TRUE)
          custom_zip_file <- file.path(custom_dir, paste0(custom_name, ".zip"))
          custom_dir_is_temporary <- TRUE
          robust.download.raster(custom_source,
                                 custom_zip_file,
                                 min_size_mb = custom.raster.min.size.mb)
          utils::unzip(custom_zip_file, exdir = custom_dir, overwrite = TRUE)
          tif_files_all <- sort(list.files(custom_dir, pattern = "\\.(tif|tiff)$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE))
          if (length(tif_files_all) == 0) stop("No .tif/.tiff found in downloaded custom zip: ", custom_name)
          tif_files <- character(0)
          for (tif_file in tif_files_all) {
            tif_ok <- FALSE
            tif_size <- suppressWarnings(file.info(tif_file)$size)
            if (!is.na(tif_size) && tif_size > 0) {
              tif_ok <- tryCatch({
                test_rast <- terra::rast(tif_file)
                ok <- inherits(test_rast, "SpatRaster") && terra::nlyr(test_rast) >= 1
                rm(test_rast)
                invisible(gc())
                ok
              }, error = function(e) FALSE)
            }
            if (tif_ok) {
              tif_files <- c(tif_files, tif_file)
            } else {
              if (verbose) message("Skipping empty or unreadable .tif/.tiff file in downloaded zip ", custom_name, ": ", basename(tif_file))
            }
          }
          if (length(tif_files) == 0) stop("No readable non-empty .tif/.tiff files found in downloaded custom zip: ", custom_name)
          if (!custom_var_name_provided) custom_var_name <- tools::file_path_sans_ext(basename(tif_files))
          if (verbose) message("Found ", length(tif_files), " .tif/.tiff files in downloaded zip: ", custom_name)
          custom_raster <- terra::rast(tif_files)
        } else {
          if (!(file.exists(custom_tif_file) && file.size(custom_tif_file) > 1e6 && !redownload.rasters)) {
            if (grepl("^https?://", custom_source)) {
              if (verbose) message("Downloading custom raster tif directly from URL: ", custom_name)
              robust.download.raster(custom_source,
                                     custom_tif_file,
                                     min_size_mb = custom.raster.min.size.mb)
            } else {
              if (!file.exists(custom_tif_file)) stop("Local custom raster tif not found: ", custom_tif_file)
              if (verbose) message("Using local custom raster tif directly: ", custom_name)
            }
          } else {
            if (verbose) message("Custom raster already present - skipping download: ", custom_name)
          }
          custom_raster <- terra::rast(custom_tif_file)
        }
        if (!inherits(custom_raster, "SpatRaster")) stop("Failed to read custom raster as SpatRaster from: ", custom_name)
        custom_raster_crs <- terra::crs(custom_raster)
        custom_crs_i <- if (is.null(custom.env.rasters.crs)) NA_character_ else if (length(custom.env.rasters.crs) == 1) custom.env.rasters.crs[[1]] else custom.env.rasters.crs[[i]]
        if (is.na(custom_raster_crs) || !nzchar(custom_raster_crs)) {
          if (!is.na(custom_crs_i) && nzchar(custom_crs_i)) {
            terra::crs(custom_raster) <- custom_crs_i
          } else {
            stop("No CRS detected for custom raster: ", basename(custom_source), " - provide custom.env.rasters.crs")
          }
        }
        if (terra::nlyr(custom_raster) == 1 && (is.null(names(custom_raster)) || names(custom_raster) == "" || names(custom_raster) == "lyr.1")) {
          names(custom_raster) <- if (length(custom_var_name) == 1) custom_var_name else custom_name
        }
        raster_crs <- terra::crs(custom_raster)
        coord_env <- coordinate_vector_env
        coord_bg <- coordinate_vector_bg
        if (!terra::same.crs(coordinate_vector_env, custom_raster)) {
          if (verbose) message("Projecting coordinates to match custom raster CRS: ", custom_name)
          coord_env <- terra::project(coordinate_vector_env, raster_crs)
          if (!is.null(coordinate_vector_bg)) coord_bg <- terra::project(coordinate_vector_bg, raster_crs)
        }
        if (length(custom_var_name) == terra::nlyr(custom_raster)) {
          names(custom_raster) <- custom_var_name
        } else if (length(custom_var_name) == 1) {
          names(custom_raster) <- if (terra::nlyr(custom_raster) > 1) paste0(custom_var_name, "_", seq_len(terra::nlyr(custom_raster))) else custom_var_name
        } else {
          stop("custom.env.rasters.variable.names for ", custom_name, " must have length 1 or match terra::nlyr(custom_raster)")
        }
        result <- extract.and.cache.env.dataset(dataset_name = custom_name,
                                                raster_object = custom_raster,
                                                coord_env = coord_env,
                                                coord_bg = coord_bg,
                                                environmental_dataset = environmental_dataset,
                                                background.data = background.data,
                                                output.dir = intermediate_files_dir,
                                                overwrite = overwrite,
                                                generate.background.data = generate.background.data,
                                                verbose = verbose)
        environmental_dataset <- result$environmental_dataset
        background.data <- result$background.data
        extracted_custom_names <- names(custom_raster)
        if (is.null(extracted_custom_names) || length(extracted_custom_names) == 0 || any(!nzchar(extracted_custom_names))) {
          stop("Custom raster layer names are missing after extraction for: ", custom_name)
        }
        cols_occ <- intersect(extracted_custom_names, names(environmental_dataset))
        if (length(cols_occ) != length(extracted_custom_names)) stop("Custom raster columns missing from occurrence data after extraction - extraction seems to not have worked")
        write.csv(environmental_dataset[, extracted_custom_names, drop = FALSE], custom_occ_csv_file)
        if (!is.null(background.data) && generate.background.data) {
          cols_bg <- intersect(extracted_custom_names, names(background.data))
          if (length(cols_bg) != length(extracted_custom_names)) stop("Custom raster columns missing from background data after extraction - extraction seems to not have worked")
          write.csv(background.data[, extracted_custom_names, drop = FALSE], custom_bg_csv_file)
        }
        suppressWarnings(try(terra::tmpFiles(remove = TRUE), silent = TRUE))
        if (exists("custom_raster")) rm(custom_raster)
        if (exists("result")) rm(result)
        invisible(gc())
        if (delete.intermediate.files.folders && custom_dir_is_temporary && file.exists(custom_occ_csv_file) && (!generate.background.data || file.exists(custom_bg_csv_file))) unlink(custom_dir, recursive = TRUE, force = TRUE)
      }
    }
  }

  # Write final CSV files
  if (!"ID" %in% colnames(environmental_dataset) && !is.null(rownames(environmental_dataset))) environmental_dataset$ID <- rownames(environmental_dataset)
  keep_front <- intersect(unique(c("ID", "Species", "Latitude", "Longitude", latitude.col, longitude.col)), colnames(environmental_dataset)) #front-load canonical + aliases
  ordered_cols <- c(keep_front, setdiff(colnames(environmental_dataset), keep_front))
  environmental_dataset <- environmental_dataset[, ordered_cols, drop = FALSE]
  write.csv(environmental_dataset, file = occurrence_output_file, row.names = FALSE)
  occurrence_output_file_display <- sub("^\\.([/\\\\])+", "./", occurrence_output_file)
  if (verbose) message("")
  if (verbose) message("-- Extraction finished --")
  if (verbose) message("Occurrence environmental data was saved to: ", occurrence_output_file_display)
  if (generate.background.data) {
    background_env_file <- file.path(output.dir, csv.background.out.file)
    background_env_file_display <- sub("^\\.([/\\\\])+", "./", background_env_file)
    write.csv(background.data, file = background_env_file, row.names = FALSE)
    if (verbose) message("Background environmental data saved to: ", background_env_file_display)
  }
}
