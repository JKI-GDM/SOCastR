# ==============================================================================
# SCRIPT HEADER
# ==============================================================================
# Title:         SOCastR - Soil Organic Carbon Prediction Workflow
# Description:   A comprehensive digital soil mapping workflow for predicting 
#                soil organic carbon (SOC) content using Random Forest and 
#                Quantile Regression Forest models, with spatial cross-validation,
#                forward feature selection, and uncertainty quantification.
#
# Author:        Markus Möller
# Institution:   Julius Kühn Institute
# Date Created:  2025-10-23
# Last Modified: 2025-11-01
# Version:       1.0
#
# Purpose:       This script implements a complete DSM workflow including:
#                - Spatial data partitioning with stratified sampling
#                - Forward feature selection (FFS) using CAST package
#                - Random Forest model training with spatial cross-validation
#                - Model validation on independent test set
#                - Wall-to-wall spatial prediction
#                - Model uncertainty via Quantile Regression Forest (QRF)
#                - Distance-based uncertainty via Dissimilarity Index (DI/AOA)
#                - Tile-based processing for large rasters
#
# Requirements:  - Sample point data (shapefile) with SOC measurements
#                - Raster covariate stack (GeoTIFF)
#                - Sufficient RAM for model training
#                - Multi-core CPU recommended for parallel processing
#
# Output:        - SOC prediction maps (RF and QRF quantiles)
#                - Uncertainty maps (prediction interval width, DI, AOA)
#                - Model performance metrics and validation plots
#                - Variable importance rankings
#                - Spatial CV diagnostics
#
# License:       MIT
# Citation:      [Add citation information if applicable]
#
# Contact:       markus.moeller@julius-kuehn.de
# ==============================================================================

# ==============================================================================
# WORKFLOW OVERVIEW
# ==============================================================================
# STEP 1:  Load and check required R packages
# STEP 2:  Define main SOCastR function with parameters
# STEP 3:  Load sample data (point shapefile) and covariate rasters
# STEP 4:  Extract covariate values using 3x3 neighborhood filter
# STEP 5:  Spatial data partitioning (train/test split) using spatial folds
# STEP 6:  Spatial cross-validation setup using KNNDM
# STEP 7:  Forward feature selection (FFS) with spatial CV
# STEP 8:  Train final Random Forest model on training data
# STEP 9:  External validation on independent test set
# STEP 10: Retrain model on all data for spatial prediction
# STEP 11: Create SOC prediction map using Random Forest
# STEP 12: Model uncertainty quantification using Quantile Regression Forest
# STEP 13: Distance-based uncertainty using Dissimilarity Index (AOA)
# STEP 14: Generate summary statistics and visualization outputs
# ==============================================================================

# ==============================================================================
# INITIAL SETTINGS
# ==============================================================================

# Print workflow header
cat("=" , rep("=", 79), "\n", sep = "")
cat("SOCastR: Soil Organic Carbon (SOC) prediction workflow with uncertainties\n")
cat("=" , rep("=", 79), "\n\n", sep = "")

# ------------------------------------------------------------------------------
# Load required packages with version checking
# ------------------------------------------------------------------------------
# This section automatically checks for required packages and installs them
# if they are missing. All packages are loaded with library() after installation.

required_packages <- c("CAST",         # Spatial prediction and validation
                       "caret",        # Machine learning framework
                       "classInt",     # Class intervals for mapping
                       "dplyr",        # Data manipulation
                       "doParallel",   # Parallelization
                       "terra",        # Raster data handling
                       "tidyterra",    # Visualise SpatRaster objects
                       "sf",           # Vector data handling
                       "randomForest", # Random Forest algorithm
                       "RColorBrewer", # Color schemes for maps
                       "quantregForest", # Quantile regression forests for uncertainty
                       "ggplot2",      # Visualization
                       "viridis",      # Color palettes
                       "gridExtra",    # Multiple plots
                       "grid")         # Graphics

# Check and install missing packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing required package:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ==============================================================================
# MAIN FUNCTION: SOCastR
# ==============================================================================
#' SOCastR: Soil Organic Carbon Prediction Workflow
#'
#' @description
#' Main function to perform digital soil mapping of soil organic carbon using
#' Random Forest and Quantile Regression Forest with spatial cross-validation.
#'
#' @param working_dir Character. Path to working directory.
#' @param input_dir Character. Name of input directory (relative to working_dir).
#' @param output_dir Character. Name of output directory (relative to working_dir).
#' @param samples Character. Filename of sample point data (shapefile).
#' @param covariates Character. Filename of covariate raster stack (GeoTIFF).
#' @param soc_column Character. Name of column containing SOC values in sample data.
#' @param n.tile Integer. Number of tiles per dimension for large raster processing (default: 4).
#' @param model_uncertainty Logical. Whether to compute QRF uncertainty (default: TRUE).
#' @param distance_uncertainty Logical. Whether to compute DI/AOA uncertainty (default: TRUE).
#'
#' @return Creates output files in output_dir including:
#'   - SOC prediction rasters
#'   - Uncertainty maps
#'   - Model performance metrics
#'   - Validation plots
#'   - Variable importance rankings
#'
#' @details
#' The function implements a complete DSM workflow with:
#' - Spatial stratified train/test split using CreateSpacetimeFolds
#' - Forward feature selection with spatial CV
#' - External validation on independent test set
#' - Quantile regression for model uncertainty
#' - Area of Applicability for distance-based uncertainty
#' - Automatic tile-based processing for large rasters (>5M cells)
#'
#' @examples
#' SOCastR(working_dir = "~/SOCastR/",
#'         input_dir = "input",
#'         output_dir = "output",
#'         samples = "SAMPLES_EPSG25832.shp",
#'         covariates = "COVARIATS_EPSG25832.tif",
#'         soc_column = "SOC",
#'         n.tile = 4,
#'         model_uncertainty = TRUE,
#'         distance_uncertainty = TRUE)
#'
#' @export
SOCastR <- function(working_dir,
                    input_dir,
                    output_dir,
                    samples,
                    covariates,
                    soc_column,
                    n.tile = 4,
                    model_uncertainty=TRUE,
                    distance_uncertainty=TRUE){
  
  # ----------------------------------------------------------------------------
  # Environment setup
  # ----------------------------------------------------------------------------
  # Set random seed for reproducibility of spatial folds and model training
  set.seed(42)
  
  # Set working directory
  setwd(working_dir)
  getwd()
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # ============================================================================
  # DATA LOADING AND PREPARATION
  # ============================================================================
  cat("=" , rep("=", 79), "\n", sep = "")
  cat("DATA LOADING AND PREPARATION\n")
  cat("=" , rep("=", 79), "\n\n", sep = "")
  
  # Load sample point data (shapefile with SOC measurements)
  # Expected format: sf object with geometry column and SOC values
  soc_samples <- st_read(file.path(getwd(),input_dir,samples))
  cat(paste("Loaded", nrow(soc_samples), "sample points\n"))
  
  # Load raster covariate stack (GeoTIFF with multiple bands)
  # All bands should have identical extent, resolution, and CRS
  covariates <- rast(file.path(getwd(),input_dir,covariates))
  cat(paste("Loaded", nlyr(covariates), "covariate layers\n"))
  cat(paste("Raster resolution:", res(covariates)[1], "x", res(covariates)[2], "m\n"))
  
  # Rename covariate layers to generic names (COV1, COV2, etc.)
  # This ensures consistent naming throughout the workflow
  names(covariates) <- paste0("COV",1:length(names(covariates)))
  print("Data loaded successfully")
  
  # ============================================================================
  # HELPER FUNCTIONS
  # ============================================================================
  
  # --------------------------------------------------------------------------
  #' PlotR: Plot raster data with class intervals
  #'
  #' @description
  #' Creates a classified raster plot using ggplot2 with customizable color
  #' schemes and class interval methods.
  #'
  #' @param raster_result SpatRaster. Raster to be plotted.
  #' @param color_schema Character. RColorBrewer palette name (e.g., "YlOrBr").
  #' @param legend_title Character. Title for the legend.
  #' @param revers Logical. Whether to reverse color palette (default: FALSE).
  #' @param accuracy Numeric. Decimal precision for break labels (default: 0.01).
  #' @param ClassIntervallMethod Character. Method for class intervals: 
  #'        "quantile", "equal", "jenks", "pretty", "sd" (default: "quantile").
  #' @param n_classes Integer. Number of classes for classification (default: 9).
  #'
  #' @return ggplot object
  #' 
  #' @details
  #' Uses classInt package to calculate breaks and classifies the raster
  #' accordingly. Supports all RColorBrewer palettes and can extend to more
  #' than 9 classes using colorRampPalette.
  # --------------------------------------------------------------------------
  PlotR <- function(raster_result,
                    color_schema,
                    legend_title,
                    revers = FALSE,
                    accuracy = 0.01,
                    ClassIntervallMethod = "quantile",
                    n_classes = 9){
    
    # Calculate class breaks using classInt
    raster_values <- values(raster_result, na.rm = TRUE)
    breaks_result <- classIntervals(raster_values, n = n_classes, style = ClassIntervallMethod)
    breaks_values <- breaks_result$brks
    
    # Create reclassification matrix: from, to, new_value
    rcl <- cbind(breaks_values[-length(breaks_values)],
                 breaks_values[-1],
                 1:n_classes)
    
    # Classify raster based on breaks
    raster_class <- classify(raster_result, rcl, include.lowest = TRUE)
    
    # Create meaningful labels for each class
    breaks_labels <- sprintf(paste0("%.", abs(log10(accuracy)), "f"), breaks_values)
    class_labels <- paste0(breaks_labels[-length(breaks_labels)],
                           " - ",
                           breaks_labels[-1])
    
    # Set factor levels with labels
    levels(raster_class) <- data.frame(value = 1:n_classes,
                                       class = class_labels)
    
    # Generate color palette
    if(n_classes <= 9){
      colors <- brewer.pal(n_classes, color_schema)
    } else {
      # Extend palette if more than 9 classes needed
      base <- brewer.pal(9, color_schema)
      colors <- colorRampPalette(base)(n_classes)
    }
    
    # Reverse colors if requested
    if(revers){
      colors <- rev(colors)
    }
    
    # Create ggplot
    ggplot() +
      geom_spatraster(data = raster_class) +
      scale_fill_manual(values = colors,
                        name = legend_title,
                        na.value = "lightgrey") +
      theme_minimal() +
      theme(panel.background = element_rect(fill = "white", colour = NA),
            plot.background = element_rect(fill = "white", colour = NA),
            panel.grid.major = element_line(colour = "black", linewidth = 0.2),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
            axis.text = element_text(size = 11, colour = "black"),
            legend.key.height = unit(0.6, "cm"),
            legend.text = element_text(size = 11))
  }
  
  # --------------------------------------------------------------------------
  #' SFPlotR: Plot point data with train/test split
  #'
  #' @description
  #' Creates a map showing spatial distribution of training and test samples.
  #'
  #' @param sf_points sf object. Point data with train/test labels.
  #' @param split_col Character. Column name containing split labels (default: "split").
  #' @param legend_title Character. Title for legend (default: "Dataset").
  #' @param colors Named vector. Colors for train/test sets.
  #'
  #' @return ggplot object
  #'
  #' @details
  #' Useful for visual assessment of spatial distribution of training and
  #' test samples. Helps identify potential spatial clustering issues.
  # --------------------------------------------------------------------------
  SFPlotR <- function(sf_points,
                      split_col = "split",
                      legend_title = "Dataset",
                      colors = c("Train" = "#0077BB", "Test" = "#CC3311")) {
    
    ggplot() +
      geom_sf(data = sf_points,
              aes_string(color = split_col),
              size = 2,
              alpha = 0.7) +
      scale_color_manual(values = colors,
                         name = legend_title,
                         na.value = "lightgrey") +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_line(colour = "black", linewidth = 0.2),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        axis.text = element_text(size = 11, colour = "black"),
        legend.key.height = unit(1, "cm"),
        legend.text = element_text(size = 11),
        legend.position = "right"
      ) +
      labs(title = "Training and Test Samples",
           x = "Longitude", y = "Latitude")
  }
  
  # --------------------------------------------------------------------------
  #' extract_8_neighbors: Extract raster values using 3x3 neighborhood
  #'
  #' @description
  #' Extracts median values from a 3x3 pixel window (focal cell + 8 neighbors)
  #' for each sample point. This reduces noise and improves representativeness
  #' of covariate values at point locations.
  #'
  #' @param raster_stack SpatRaster. Raster stack with covariate layers.
  #' @param points sf object. Point locations for extraction.
  #'
  #' @return data.frame with extracted median values for each point and layer.
  #'
  #' @details
  #' For each point:
  #' 1. Identify the focal cell containing the point
  #' 2. Extract values from focal cell and 8 adjacent cells
  #' 3. Calculate median across 9 cells for each layer
  #' 4. Handle NA values gracefully
  #'
  #' This approach is more robust than simple point extraction, especially
  #' when dealing with:
  #' - GPS positioning uncertainty
  #' - Covariate data with localized artifacts
  #' - Misalignment between sample locations and raster cells
  # --------------------------------------------------------------------------
  extract_8_neighbors <- function(raster_stack, points) {
    
    # Get cell numbers for point locations
    cell_numbers <- cellFromXY(raster_stack, st_coordinates(points))
    
    # Initialize results list
    result_list <- list()
    
    # Loop through each point
    for(i in seq_along(cell_numbers)) {
      if(!is.na(cell_numbers[i])) {
        # Get adjacent cells (8 neighbors in all directions)
        adj_cells <- adjacent(raster_stack, cell_numbers[i], directions = 8)
        
        # Include the focal cell
        all_cells <- c(cell_numbers[i], adj_cells[,2])
        
        # Extract values from all 9 cells for each layer
        extracted_values <- extract(raster_stack, all_cells)
        
        # Calculate median across the 9 cells (focal + 8 neighbors)
        # Median is more robust to outliers than mean
        median_values <- apply(extracted_values[, ], 2, median, na.rm = TRUE)
        result_list[[i]] <- median_values
      } else {
        # If point falls outside raster extent, assign NA
        result_list[[i]] <- rep(NA, nlyr(raster_stack))
      }
    }
    
    # Combine results into data frame
    result_matrix <- do.call(rbind, result_list)
    result_df <- data.frame(ID = 1:nrow(result_matrix), result_matrix)
    names(result_df) <- c("ID", names(raster_stack))
    
    return(result_df)
  }
  
  # ============================================================================
  # COVARIATE EXTRACTION
  # ============================================================================
  cat("=" , rep("=", 79), "\n", sep = "")
  cat("EXTRACT (FILTERED) VALUES\n")
  cat("=" , rep("=", 79), "\n\n", sep = "")
  
  # Extract covariate values using 3x3 pixel median filter
  # This reduces the impact of GPS positioning errors and local artifacts
  extracted_covariates <- extract_8_neighbors(covariates, soc_samples)
  
  # Combine SOC values with extracted covariates
  model_data <- cbind(
    SOC = soc_samples[[paste(soc_column)]],
    extracted_covariates[, -1]  # Remove ID column
  )
  
  print("Covariate extraction completed")
  print(paste("Final dataset dimensions:", paste(dim(model_data), collapse = " x ")))
  
  # Remove samples with NA values (points outside raster extent or with missing covariate data)
  complete_idx <- complete.cases(model_data)
  model_data_clean <- model_data[complete_idx, ]
  soc_samples_clean <- soc_samples[complete_idx, ]
  
  # Report data cleaning results
  cat(paste("\nOriginal samples:", nrow(model_data), "\n"))
  cat(paste("Clean samples (after NA removal):", nrow(model_data_clean), "\n"))
  cat(paste("Removed samples:", nrow(model_data) - nrow(model_data_clean), "\n\n"))
  
  # Save basic data summary
  data_summary <- data.frame(
    Metric = c("Total samples", "Clean samples", "Removed samples",
               "Number of covariates", "SOC mean", "SOC sd", "SOC min", "SOC max"),
    Value = c(nrow(model_data), nrow(model_data_clean),
              nrow(model_data) - nrow(model_data_clean),
              nlyr(covariates),
              round(mean(model_data_clean$SOC), 2),
              round(sd(model_data_clean$SOC), 2),
              round(min(model_data_clean$SOC), 2),
              round(max(model_data_clean$SOC), 2))
  )
  
  # Export summary statistics
  write.csv(data_summary, file.path(getwd(),output_dir, "ExtractValues_SampleDataSummary.csv"), row.names = FALSE)
  
  # ============================================================================
  # SPATIAL DATA PARTITION: STRATIFIED TRAIN/TEST SPLIT
  # ============================================================================
  # Use spatial blocking to create train/test split that accounts for spatial
  # autocorrelation. This ensures more realistic validation performance estimates.
  
  set.seed(42)  # Ensure reproducibility
  
  cat("=" , rep("=", 79), "\n", sep = "")
  cat("SPATIAL DATA PARTITION\n")
  cat("=" , rep("=", 79), "\n\n", sep = "")
  
  # Create spatial folds for stratified splitting (k=5 folds)
  # CreateSpacetimeFolds performs spatially stratified splitting to ensure
  # training and test sets are spatially separated, reducing optimistic bias
  spatial_folds_obj <- CreateSpacetimeFolds(
    soc_samples_clean,
    spacevar = "geometry",
    k = 5
  )
  
  # Use first fold as test set, remaining folds as training set
  # Alternative: iterate through all folds for full k-fold CV
  test_idx_spatial <- spatial_folds_obj$indexOut[[1]]
  train_idx_spatial <- spatial_folds_obj$index[[1]]
  
  # Prepare train and test datasets
  train_data <- model_data_clean[train_idx_spatial, ]
  test_data <- model_data_clean[test_idx_spatial, ]
  train_samples <- soc_samples_clean[train_idx_spatial, ]
  test_samples <- soc_samples_clean[test_idx_spatial, ]
  
  # Report split statistics
  cat(paste("Spatial stratified split using CreateSpacetimeFolds\n"))
  cat(paste("Training samples:", nrow(train_data),
            "(", round(100 * length(train_idx_spatial) / nrow(model_data_clean), 1), "%)\n"))
  cat(paste("Test samples:", nrow(test_data),
            "(", round(100 * length(test_idx_spatial) / nrow(model_data_clean), 1), "%)\n"))
  
  # Visualization of the spatial split
  # Transform to WGS84 (EPSG:4326) for geographic visualization
  soc_samples_clean_epsg4326 <- st_transform(soc_samples_clean, crs = 4326)
  coords_all <- st_coordinates(soc_samples_clean_epsg4326)
  
  # Create data frame for visualization
  visual_split <- data.frame(
    x = coords_all[, 1],
    y = coords_all[, 2],
    SOC = model_data_clean$SOC,
    split = "Test"
  )
  visual_split$split[train_idx_spatial] <- "Train"
  visual_split$split <- factor(visual_split$split, levels = c("Train", "Test"))
  
  # Convert to sf object for spatial plotting
  visual_split_sf <- st_as_sf(visual_split,
                               coords = c("x", "y"),
                               crs = 4326)
  
  # Create spatial map of train/test split
  p_spatial <- SFPlotR(visual_split_sf)
  
  # Export map
  ggsave(file.path(getwd(),output_dir,"SpatialDataPartition_MapTrainTestSplit.png"),
         p_spatial,
         width = 5,
         height = 6,
         dpi = 300)
  
  # Distribution assessment: compare SOC distributions in train vs test
  p_soc_dist <- ggplot(visual_split, aes(x = SOC, fill = split)) +
    geom_density(alpha = 0.6) +
    scale_fill_manual(values = c("Train" = "#2E86AB", "Test" = "#A23B72")) +
    labs(
      title = "SOC Distribution: Train vs Test",
      x = "SOC (g/kg)",
      y = "Density",
      fill = "Set"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Export distribution plot
  ggsave(file.path(getwd(),output_dir,"SpatialDataPartition_SocDistSpatialBlock.png"),
         p_soc_dist,
         width = 7,
         height = 6,
         dpi = 300)
  
  # Statistical comparison using Kolmogorov-Smirnov test
  # Tests whether train and test distributions differ significantly
  ks_result <- ks.test(train_data$SOC, test_data$SOC)
  
  # Calculate summary statistics for both sets
  split_stats <- data.frame(
    Metric = c(
      "N_train", "N_test",
      "Train_SOC_mean", "Train_SOC_sd", "Train_SOC_min", "Train_SOC_max",
      "Test_SOC_mean", "Test_SOC_sd", "Test_SOC_min", "Test_SOC_max",
      "SOC_mean_difference", "SOC_mean_diff_percent", "KS_statistic", "KS_p_value"
    ),
    Value = c(
      nrow(train_data), nrow(test_data),
      round(mean(train_data$SOC), 2), round(sd(train_data$SOC), 2),
      round(min(train_data$SOC), 2), round(max(train_data$SOC), 2),
      round(mean(test_data$SOC), 2), round(sd(test_data$SOC), 2),
      round(min(test_data$SOC), 2), round(max(test_data$SOC), 2),
      round(abs(mean(train_data$SOC) - mean(test_data$SOC)), 2),
      round(abs(mean(train_data$SOC) - mean(test_data$SOC)) / mean(train_data$SOC) * 100, 2),
      round(ks_result$statistic, 4), round(ks_result$p.value, 4)
    )
  )
  
  # Export statistics
  write.csv(split_stats, file.path(getwd(),output_dir,"SpatialDataPartition_TrainTestStatisticsSpatialBlock.csv"), row.names = FALSE)
  cat("Spatial stratified split visualization and statistics completed.\n\n")
  
  # ============================================================================
  # SPATIAL CROSS-VALIDATION SETUP
  # ============================================================================
  # Use KNNDM (k-fold Nearest Neighbour Distance Matching) for spatial CV
  # This ensures CV folds approximate the geographic distance distribution
  # between training and prediction areas, providing more realistic validation
  
  cat("=" , rep("=", 79), "\n", sep = "")
  cat("SPATIAL CROSS-VALIDATION SETUP (Training Data Only)\n")
  cat("=" , rep("=", 79), "\n\n", sep = "")
  
  # Set up parallel processing for KNNDM
  detectCores()
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
  
  # Create spatial CV folds using KNNDM on training data only
  # KNNDM ensures CV folds match geographic distance distribution
  spatial_folds_train <- knndm(train_samples, covariates, k = 5)
  
  # Stop parallel cluster
  stopCluster(cl)
  
  cat(paste("Created", length(spatial_folds_train$indx_train), "spatial CV folds\n"))
  
  # Geodist analysis: compare distance distributions
  # Helps assess whether CV setup is representative of prediction task
  geodist_result <- geodist(
    x = train_samples,
    modeldomain = covariates,
    cvfolds = spatial_folds_train$indx_test
  )
  
  # Plot geodist ECDF (Empirical Cumulative Distribution Function)
  p_geodist <- plot(geodist_result, unit = "km", stat = "ecdf") +
    ggtitle("Spatial CV Validation: Distance Distributions") +
    theme_minimal()
  
  # Plot geodist density
  p_geodens <- plot(geodist_result, unit = "km") +
    ggtitle("Spatial CV Validation: Distance densities") +
    theme_minimal()
  
  # Export geodist plots
  ggsave(file.path(getwd(),output_dir, "SpatialCrossValidation_GeodistEcdf.png"),
         p_geodist,
         width = 7,
         height = 5,
         dpi = 300)
  
  ggsave(file.path(getwd(),output_dir, "SpatialCrossValidation_GeodistDensity.png"),
         p_geodens,
         width = 7,
         height = 5,
         dpi = 300)
  
  cat("Geodist analysis completed and saved.\n\n")
  
  # Calculate and export geodist summary metrics
  geodist_summary <- geodist_result %>%
    group_by(what) %>%
    summarise(
      n = n(),
      min = min(dist, na.rm = TRUE),
      q25 = quantile(dist, 0.25, na.rm = TRUE),
      median = median(dist, na.rm = TRUE),
      mean = mean(dist, na.rm = TRUE),
      q75 = quantile(dist, 0.75, na.rm = TRUE),
      max = max(dist, na.rm = TRUE),
      sd = sd(dist, na.rm = TRUE))
  
  # Export geodist metrics
  write.csv(geodist_summary, file.path(getwd(),output_dir, "SpatialCrossValidation_GeodistMetrics.csv"),
            row.names = FALSE)
  
  # ============================================================================
  # FORWARD FEATURE SELECTION
  # ============================================================================
  # Use forward feature selection (FFS) to identify the most relevant
  # covariates for SOC prediction. FFS starts with an empty model and
  # iteratively adds variables that improve cross-validation performance.
  
  cat("=" , rep("=", 79), "\n", sep = "")
  cat("FORWARD FEATURE SELECTION\n")
  cat("=" , rep("=", 79), "\n\n", sep = "")
  
  # Prepare predictor list (all covariates except SOC)
  predictors <- names(train_data)[names(train_data) != "SOC"]
  cat(paste("Total available predictors:", length(predictors), "\n"))
  cat("Preparing data for FFS...\n")
  cat(paste("Number of predictors:", length(predictors), "\n"))
  cat(paste("Training samples:", nrow(train_data), "\n"))
  
  # Check for zero variance predictors and remove them
  # Zero variance predictors have no predictive power
  nzv <- nearZeroVar(train_data[, predictors], saveMetrics = TRUE)
  if(any(nzv$zeroVar)) {
    cat("Removing zero-variance predictors:",
        paste(rownames(nzv)[nzv$zeroVar], collapse = ", "), "\n")
    predictors <- predictors[!nzv$zeroVar]
  }
  
  # Check for NA values (should not occur after data cleaning)
  if(any(is.na(train_data[, predictors]))) {
    stop("NA values detected in predictors!")
  }
  
  # Ensure data frame format (remove list columns if present)
  train_data <- as.data.frame(train_data)
  train_data <- train_data[, sapply(train_data, function(x) !is.list(x))]
  
  # Set up training control with spatial CV folds
  # Use KNNDM folds created earlier for consistent validation
  train_control <- trainControl(
    method = "cv",
    index = spatial_folds_train$indx_train,
    indexOut = spatial_folds_train$indx_test,
    savePredictions = "final",
    returnResamp = "final",
    allowParallel = FALSE)  # Disable parallel for stability
  
  # Perform forward feature selection
  # FFS iteratively adds variables that reduce RMSE in spatial CV
  cat("Running forward feature selection\n")
  ffs_model <- ffs(
    predictors = train_data[, predictors, drop = FALSE],
    response = train_data$SOC,
    method = "rf",           # Random Forest
    metric = "RMSE",         # Optimize for root mean squared error
    tuneLength = 3,          # Test 3 different mtry values
    trControl = train_control,
    ntree = 100,             # 100 trees for faster FFS (will use 500 for final model)
    verbose = FALSE
  )
  
  # Report FFS results
  cat("\nForward feature selection completed!\n")
  cat(paste("Selected", length(ffs_model$selectedvars), "variables out of",
            length(predictors), "\n"))
  cat("Selected variables:\n")
  cat(paste("  -", ffs_model$selectedvars, collapse = "\n"))
  cat("\n\n")
  
  # Save FFS results
  ffs_results <- data.frame(
    Variable = ffs_model$selectedvars,
    Selection_order = 1:length(ffs_model$selectedvars)
  )
  
  # Export FFS results
  write.csv(ffs_results, file.path(getwd(),output_dir, "ForwardFeatureSelection_SelectedVariables.csv"),
            row.names = FALSE)
  
  # ============================================================================
  # FINAL MODEL TRAINING ON TRAINING DATA
  # ============================================================================
  # Train final Random Forest model using only the selected features from FFS
  # Use more trees (500) for better stability and accuracy
  
  cat("=" , rep("=", 79), "\n", sep = "")
  cat("FINAL MODEL TRAINING\n")
  cat("=" , rep("=", 79), "\n\n", sep = "")
  
  # Train final model with selected features only
  final_model <- train(
    x = train_data[, ffs_model$selectedvars, drop = FALSE],
    y = train_data$SOC,
    method = "rf",
    tuneLength = 3,          # Test 3 mtry values
    trControl = train_control,
    ntree = 500,             # 500 trees for final model
    importance = TRUE        # Calculate variable importance
  )
  
  cat("Final model training completed.\n")
  cat("Training performance (Spatial CV):\n")
  print(final_model$results)
  cat("\n")
  
  # Extract variable importance
  # Higher importance = more influential in prediction
  var_importance <- varImp(final_model)
  var_imp_df <- data.frame(
    Variable = rownames(var_importance$importance),
    Importance = var_importance$importance$Overall
  )
  var_imp_df <- var_imp_df[order(var_imp_df$Importance, decreasing = TRUE), ]
  
  # Export variable importance
  write.csv(var_imp_df, file.path(getwd(),output_dir, "FinalModel_VariableImportance.csv"),
            row.names = FALSE)
  write.csv(data.frame(final_model$results), file.path(getwd(),output_dir, "FinalModel_Accuracy.csv"),
            row.names = FALSE)
  
  # Plot variable importance
  p_varimp <- ggplot(var_imp_df, aes(x = reorder(Variable, Importance), y = Importance)) +
    geom_col(fill = "#2E86AB") +
    coord_flip() +
    labs(
      title = "Variable Importance (Random Forest)",
      x = "Variable",
      y = "Importance"
    ) +
    theme_minimal()
  
  # Export variable importance plot
  ggsave(file.path(getwd(),output_dir,"ForwardFeatureSelection_VariableImportance.png"),
         p_varimp,
         width = 5,
         height = 5,
         dpi = 300)
  
  # ============================================================================
  # EXTERNAL VALIDATION
  # ============================================================================
  # Validate model on independent test set that was spatially separated
  # during data partitioning. This provides realistic performance estimates.
  
  cat("=" , rep("=", 79), "\n", sep = "")
  cat("EXTERNAL VALIDATION\n")
  cat("=" , rep("=", 79), "\n\n", sep = "")
  
  # Predict on test set using selected features
  test_predictions <- predict(final_model,
                               newdata = test_data[, ffs_model$selectedvars])
  
  # Calculate test set performance metrics
  test_performance <- data.frame(
    Dataset = "Test",
    RMSE = RMSE(test_predictions, test_data$SOC),
    MAE = MAE(test_predictions, test_data$SOC),
    R2 = cor(test_predictions, test_data$SOC)^2,
    Bias = mean(test_predictions - test_data$SOC),
    RMSE_percent = RMSE(test_predictions, test_data$SOC) / mean(test_data$SOC) * 100
  )
  
  # Get training CV performance for comparison
  cv_performance <- data.frame(
    Dataset = "Training_CV",
    RMSE = final_model$results$RMSE[which.min(final_model$results$RMSE)],
    MAE = final_model$results$MAE[which.min(final_model$results$RMSE)],
    R2 = final_model$results$Rsquared[which.min(final_model$results$RMSE)],
    Bias = 0,  # Not directly available from CV
    RMSE_percent = final_model$results$RMSE[which.min(final_model$results$RMSE)] /
      mean(train_data$SOC) * 100
  )
  
  # Combine performance metrics
  performance_comparison <- rbind(cv_performance, test_performance)
  performance_comparison$Difference <- c(0, test_performance$RMSE - cv_performance$RMSE)
  
  # Export performance comparison
  write.csv(performance_comparison,
            file.path(getwd(),output_dir, "Validation_PerformanceComparison.csv"),
            row.names = FALSE)
  
  cat("=== PERFORMANCE COMPARISON ===\n")
  print(performance_comparison)
  cat("\n")
  
  # Create validation scatter plot
  validation_df <- data.frame(
    Observed = test_data$SOC,
    Predicted = test_predictions)
  
  p_validation <- ggplot(validation_df, aes(x = Observed, y = Predicted)) +
    geom_point(alpha = 0.6, color = "#2E86AB", size = 3) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1) +
    geom_smooth(method = "lm", color = "blue", se = TRUE) +
    labs(
      title = "",
      subtitle = sprintf("R² = %.3f | RMSE = %.2f",
                         test_performance$R2, test_performance$RMSE),
      x = "Observed SOC content [%]",
      y = "Predicted SOC content [%]"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 12, color = "gray40"),
      axis.text = element_text(size = 14, colour = "black"),
      legend.key.height = unit(1, "cm"),
      legend.text = element_text(size = 14),
      axis.title.x = element_text(face = "bold",size = 14),
      axis.title.y = element_text(face = "bold",size = 14)
    )
  
  # Export validation plot
  ggsave(file.path(getwd(),output_dir, "Validation_ScatterPlot.png"),
         p_validation,
         width = 7,
         height = 5,
         dpi = 300)
  
  # Create residual plot to assess prediction bias
  validation_df$Residual <- validation_df$Predicted - validation_df$Observed
  
  p_residuals <- ggplot(validation_df, aes(x = Predicted, y = Residual)) +
    geom_point(alpha = 0.6, color = "#A23B72", size = 3) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    geom_smooth(method = "loess", color = "blue", se = TRUE) +
    labs(
      title = "Residual Plot (Test Set)",
      x = "Predicted SOC (g/kg)",
      y = "Residuals (Predicted - Observed)"
    ) +
    theme_minimal()
  
  # Export residual plot
  ggsave(file.path(getwd(),output_dir, "Validation_ResidualPlot.png"),
         p_residuals,
         width = 6,
         height = 6,
         dpi = 300)
  
  cat("Test set validation completed and visualizations saved.\n\n")
  
  # ============================================================================
  # RETRAIN ON ALL DATA FOR FINAL PREDICTION
  # ============================================================================
  # After validation, retrain model on ALL available data (train + test)
  # This maximizes use of available information for spatial prediction
  
  cat("=" , rep("=", 79), "\n", sep = "")
  cat("RETRAIN ON ALL DATA FOR FINAL PREDICTION\n")
  cat("=" , rep("=", 79), "\n\n", sep = "")
  
  # Stop any existing parallel processing
  stopImplicitCluster()
  registerDoSEQ()  # Register sequential back end
  closeAllConnections()
  
  # Create spatial folds for full dataset
  detectCores()
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
  
  # KNNDM on complete dataset for final model CV
  spatial_folds_full <- knndm(soc_samples_clean, covariates, k = 5)
  
  stopCluster(cl)
  
  # Train final model on ALL data with PARALLEL DISABLED for stability
  final_model_full <- train(
    x = model_data_clean[, ffs_model$selectedvars, drop = FALSE],
    y = model_data_clean$SOC,
    method = "rf",
    tuneLength = 3,
    trControl = trainControl(
      method = "cv",
      index = spatial_folds_full$indx_train,
      indexOut = spatial_folds_full$indx_test,
      savePredictions = "final",
      allowParallel = FALSE  # Disable for stability
    ),
    ntree = 500,
    importance = TRUE
  )
  
  # ============================================================================
  # SPATIAL PREDICTION
  # ============================================================================
  # Create wall-to-wall SOC prediction map across entire study area
  
  cat("=" , rep("=", 79), "\n", sep = "")
  cat("SPATIAL PREDICTION\n")
  cat("=" , rep("=", 79), "\n\n", sep = "")
  
  # Predict SOC across entire raster extent
  soc_prediction <- terra::predict(covariates, final_model_full, na.rm = TRUE)
  names(soc_prediction) <- "SOC_predicted"
  
  cat("SOC prediction map created.\n")
  
  # Export prediction raster
  writeRaster(soc_prediction,
              file.path(getwd(),output_dir, "FinalPrediction_SocRaster.tif"),
              overwrite = TRUE)
  
  cat("Prediction map saved.\n\n")
  
  # Create and export prediction map visualization
  plot_prediction <- PlotR(soc_prediction,"YlOrBr","SOC [%]")
  ggsave(file.path(getwd(),output_dir,"FinalPrediction_MapSocRF.png"),
         plot_prediction,
         width = 6,
         height = 7, dpi = 300)
  
  # ============================================================================
  # MODEL UNCERTAINTY: QUANTILE REGRESSION FOREST
  # ============================================================================
  # Use Quantile Regression Forest to estimate prediction intervals
  # QRF provides model-based uncertainty by quantifying prediction distribution
  
  if(model_uncertainty==TRUE){
    cat("=", rep("=", 79), "\n", sep = "")
    cat("MODEL UNCERTAINTY: QUANTILE REGRESSION FOREST UNCERTAINTY LAYERS (TILE-BASED)\n")
    cat("=", rep("=", 79), "\n\n", sep = "")
    
    # Prepare training data for QRF
    train_x <- model_data_clean[, ffs_model$selectedvars, drop = FALSE]
    train_y <- model_data_clean$SOC
    
    # Fit Quantile Regression Forest model
    # QRF maintains full prediction distribution, not just point estimates
    cat("Training quantile regression forest model...\n")
    qrf_model <- quantregForest(
      x = train_x,
      y = train_y,
      ntree = 500,
      importance = TRUE
    )
    
    # Calculate accuracy metrics for median prediction
    qrf_model.lm <- lm(train_y ~ qrf_model$predicted)
    qrf_model.lm.R2 <- summary(qrf_model.lm)$r.squared
    qrf_model.lm.RMSE <- sqrt(mean((train_y - qrf_model$predicted)^2))
    
    # Export QRF training results
    qrf_model.accuracy <- data.frame(
      R_squared = qrf_model.lm.R2,
      RMSE = qrf_model.lm.RMSE,
      Method = "Quantile_Regression_Forest"
    )
    
    write.csv(qrf_model.accuracy,
              file.path(getwd(), output_dir, "QuantilePrediction_Accuracy.csv"),
              row.names = FALSE)
    
    cat("QRF model training completed.\n\n")
    
    # --------------------------------------------------------------------------
    # Ensure raster layer names match model variables
    # --------------------------------------------------------------------------
    cat("Checking raster layer names...\n")
    cat("Model expects:", paste(ffs_model$selectedvars, collapse = ", "), "\n")
    cat("Raster has:", paste(names(covariates), collapse = ", "), "\n\n")
    
    missing_vars <- setdiff(ffs_model$selectedvars, names(covariates))
    if (length(missing_vars) > 0) {
      stop("ERROR: Missing variables in covariates: ", paste(missing_vars, collapse = ", "))
    }
    
    # Subset and reorder covariates to match model expectations
    covariates <- covariates[[ffs_model$selectedvars]]
    cat("Layer names matched successfully.\n\n")
    
    # --------------------------------------------------------------------------
    # Tile-based processing setup
    # --------------------------------------------------------------------------
    # Define quantiles for prediction intervals
    # 0.05 = 5th percentile, 0.5 = median, 0.95 = 95th percentile
    # Together these define a 90% prediction interval
    quantiles <- c(0.05, 0.5, 0.95)
    
    # Check if tiling is needed (rasters >5 million cells)
    raster_cells <- ncell(covariates)
    cat(paste("Total raster cells:", format(raster_cells, big.mark = ",")), "\n")
    
    # Determine tiling strategy based on raster size
    use_tiling <- raster_cells > 5e6
    
    if (use_tiling) {
      cat("Large raster detected. Using tile-based processing...\n\n")
      
      # ------------------------------------------------------------------------
      #' create_tiles_manual: Create spatial tiles from raster
      #'
      #' @description
      #' Divides large raster into smaller tiles for memory-efficient processing.
      #' Tiles are saved as separate GeoTIFF files and processed individually.
      #'
      #' @param raster SpatRaster. Input raster to be tiled.
      #' @param ntiles_x Integer. Number of tiles in x direction.
      #' @param ntiles_y Integer. Number of tiles in y direction.
      #' @param output_dir Character. Directory to save tiles.
      #'
      #' @return Character vector of tile file paths.
      #'
      #' @details
      #' Tiles are created by dividing the raster extent into equal-sized
      #' rectangular regions. Each tile maintains the original CRS and resolution.
      # ------------------------------------------------------------------------
      create_tiles_manual <- function(raster, ntiles_x = 4, ntiles_y = 4, output_dir = "tiles") {
        
        # Create output directory
        dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
        
        # Get raster extent
        ext_full <- ext(raster)
        xmin <- ext_full[1]
        xmax <- ext_full[2]
        ymin <- ext_full[3]
        ymax <- ext_full[4]
        
        # Calculate tile dimensions
        tile_width <- (xmax - xmin) / ntiles_x
        tile_height <- (ymax - ymin) / ntiles_y
        
        cat(paste("Creating", ntiles_x * ntiles_y, "tiles\n"))
        
        tile_files <- character()
        
        # Loop through grid to create tiles
        for (i in 1:ntiles_y) {
          for (j in 1:ntiles_x) {
            # Define tile extent
            tile_xmin <- xmin + (j - 1) * tile_width
            tile_xmax <- xmin + j * tile_width
            tile_ymin <- ymin + (i - 1) * tile_height
            tile_ymax <- ymin + i * tile_height
            
            tile_extent <- ext(tile_xmin, tile_xmax, tile_ymin, tile_ymax)
            
            # Crop raster to tile extent
            tile_raster <- crop(raster, tile_extent)
            
            # Save tile
            tile_filename <- file.path(output_dir, sprintf("tile_%02d_%02d.tif", i, j))
            writeRaster(tile_raster, tile_filename, overwrite = TRUE)
            tile_files <- c(tile_files, tile_filename)
          }
        }
        
        cat(paste("Created", length(tile_files), "tiles\n"))
        return(tile_files)
      }
      
      # Create tiles for processing
      cat("Creating raster tiles...\n")
      tiles_dir <- file.path(getwd(), output_dir, "qrf_tiles")
      quantile_tiles_dir <- file.path(getwd(), output_dir, "qrf_quantile_tiles")
      
      tile_files <- create_tiles_manual(
        raster = covariates,
        ntiles_x = n.tile,
        ntiles_y = n.tile,
        output_dir = tiles_dir
      )
      
      dir.create(quantile_tiles_dir, showWarnings = FALSE, recursive = TRUE)
      
      # ------------------------------------------------------------------------
      # Parallel tile processing
      # ------------------------------------------------------------------------
      library(parallel)
      library(doParallel)
      library(foreach)
      
      n_cores <- detectCores() - 1
      start_time <- Sys.time()
      
      cat(paste("\nProcessing", length(tile_files), "tiles across", n_cores, "cores...\n\n"))
      
      # Platform-specific parallel processing
      if (.Platform$OS.type == "windows") {
        # Windows: PSOCK cluster (required for Windows)
        cl <- makeCluster(n_cores, type = "PSOCK")
        registerDoParallel(cl)
        
        # Export necessary objects to cluster
        clusterExport(cl, varlist = c("qrf_model", "quantiles", "quantile_tiles_dir"),
                      envir = environment())
        
        # Load required packages on each worker
        clusterEvalQ(cl, {
          library(terra)
          library(quantregForest)
        })
        
        # Process tiles in parallel using foreach
        tile_results <- foreach(
          i = seq_along(tile_files),
          .packages = c("terra", "quantregForest"),
          .combine = 'c',
          .inorder = TRUE
        ) %dopar% {
          
          tile_path <- tile_files[i]
          tile_raster <- rast(tile_path)
          tile_basename <- basename(tile_path)
          tile_name <- sub(".tif", "", tile_basename)
          
          # Predict all quantiles for this tile
          quantile_results <- list()
          for (q_idx in seq_along(quantiles)) {
            pred <- terra::predict(
              tile_raster,
              qrf_model,
              na.rm = TRUE,
              fun = function(model, data) {
                predict(model, data, what = quantiles[q_idx])
              }
            )
            
            # Save quantile tile
            q_file <- file.path(quantile_tiles_dir,
                                paste0("Q", quantiles[q_idx]*100, "_", tile_name, ".tif"))
            writeRaster(pred, q_file, overwrite = TRUE)
            quantile_results[[q_idx]] <- q_file
          }
          
          return(list(quantile_results))
        }
        
        stopCluster(cl)
        
      } else {
        # Linux/Unix: Fork-based parallelization (more efficient)
        cat("Unix/Linux detected: Using fork-based parallelization\n")
        
        # Limit cores to avoid memory issues on some systems
        n_cores <- min(n_cores, 4)
        cat(paste("Reducing cores to", n_cores, "to avoid memory issues\n"))
        
        # Process tiles using mclapply
        tile_results <- mclapply(seq_along(tile_files), function(i) {
          
          tile_path <- tile_files[i]
          tile_raster <- rast(tile_path)
          tile_basename <- basename(tile_path)
          tile_name <- sub(".tif", "", tile_basename)
          
          # Predict all quantiles for this tile
          quantile_results <- list()
          for (q_idx in seq_along(quantiles)) {
            pred <- terra::predict(
              tile_raster,
              qrf_model,
              na.rm = TRUE,
              fun = function(model, data) {
                predict(model, data, what = quantiles[q_idx])
              }
            )
            
            # Save quantile tile
            q_file <- file.path(quantile_tiles_dir,
                                paste0("Q", quantiles[q_idx]*100, "_", tile_name, ".tif"))
            writeRaster(pred, q_file, overwrite = TRUE)
            quantile_results[[q_idx]] <- q_file
          }
          
          return(quantile_results)
        },
        mc.cores = n_cores,
        mc.preschedule = FALSE,
        mc.cleanup = TRUE)
      }
      
      end_time <- Sys.time()
      processing_time <- difftime(end_time, start_time, units = "mins")
      cat(paste("\nTile processing completed in", round(processing_time, 1), "minutes\n\n"))
      
      # ------------------------------------------------------------------------
      # Merge tiles into final rasters
      # ------------------------------------------------------------------------
      cat("Merging tiles into final quantile rasters...\n")
      
      # For each quantile, merge its tiles
      quantile_layers <- list()
      for (q_idx in seq_along(quantiles)) {
        q_pattern <- paste0("Q", quantiles[q_idx]*100, "_tile_.*\\.tif$")
        q_files <- list.files(quantile_tiles_dir, pattern = q_pattern, full.names = TRUE)
        
        cat(paste("Merging", length(q_files), "tiles for quantile", quantiles[q_idx], "...\n"))
        
        # Load all tiles for this quantile
        q_rasters <- lapply(q_files, rast)
        
        # Merge using mosaic (more robust than merge for tiles)
        q_merged <- do.call(mosaic, q_rasters)
        names(q_merged) <- paste0("SOC_q", quantiles[q_idx]*100)
        
        quantile_layers[[q_idx]] <- q_merged
      }
      
      # Combine all quantiles into one stack
      soc_qrf_pred <- rast(quantile_layers)
      cat("Tile merging completed.\n\n")
      
      # Optional: Clean up tile directories to save disk space
      # unlink(tiles_dir, recursive = TRUE)
      # unlink(quantile_tiles_dir, recursive = TRUE)
      
    } else {
      # ------------------------------------------------------------------------
      # Small raster: process directly without tiling
      # ------------------------------------------------------------------------
      cat("Small raster detected. Processing directly (no tiling needed)...\n\n")
      
      # Create progress bar for user feedback
      pb <- txtProgressBar(min = 0, max = length(quantiles), style = 3)
      quantile_list <- list()
      
      # Loop through each quantile
      for (i in seq_along(quantiles)) {
        # Predict quantile across entire raster
        pred <- terra::predict(
          covariates,
          qrf_model,
          na.rm = TRUE,
          fun = function(model, data) {
            predict(model, data, what = quantiles[i])
          }
        )
        
        names(pred) <- paste0("SOC_q", quantiles[i]*100)
        quantile_list[[i]] <- pred
        setTxtProgressBar(pb, i)
      }
      
      close(pb)
      
      # Combine quantiles into single stack
      soc_qrf_pred <- rast(quantile_list)
      cat("\nQuantile predictions completed.\n\n")
    }
    
    # --------------------------------------------------------------------------
    # Post-processing (same for both tiled and non-tiled approaches)
    # --------------------------------------------------------------------------
    
    # Calculate prediction interval width (90% PI = Q95 - Q05)
    # Wider intervals indicate higher uncertainty
    soc_qrf_interval_width <- soc_qrf_pred[[3]] - soc_qrf_pred[[1]]
    names(soc_qrf_interval_width) <- "SOC_PI_width"
    
    # Combine all layers into single stack
    soc_qrf_stack <- c(soc_qrf_pred, soc_qrf_interval_width)
    
    # Export final quantile rasters
    writeRaster(soc_qrf_stack,
                file.path(getwd(), output_dir, "FinalPrediction_SocQuantileLayers.tif"),
                overwrite = TRUE)
    
    cat("Quantile prediction raster layers saved.\n")
    
    # Calculate and export summary statistics
    pi_summary <- data.frame(
      Metric = c("Mean_PI_width", "Median_PI_width", "95th_percentile_width"),
      Value = c(
        round(mean(values(soc_qrf_interval_width), na.rm = TRUE), 3),
        round(median(values(soc_qrf_interval_width), na.rm = TRUE), 3),
        round(quantile(values(soc_qrf_interval_width), 0.95, na.rm = TRUE), 3)
      )
    )
    
    write.csv(pi_summary,
              file.path(getwd(), output_dir, "QuantilePrediction_UncertaintySummary.csv"),
              row.names = FALSE)
    
    cat("Quantile-based uncertainty summary statistics saved.\n\n")
    
    # Export visualizations for each quantile
    plot_prediction <- PlotR(soc_qrf_pred[[1]],"YlOrBr","SOC [%]")
    ggsave(file.path(getwd(),output_dir,"FinalPrediction_MapQuantialSoc5Percentile.png"),
           plot_prediction, width = 6, height = 7, dpi = 300)
    
    plot_prediction <- PlotR(soc_qrf_pred[[2]],"YlOrBr","SOC [%]")
    ggsave(file.path(getwd(),output_dir,"FinalPrediction_MapQuantialSoc50Percentile.png"),
           plot_prediction, width = 6, height = 7, dpi = 300)
    
    plot_prediction <- PlotR(soc_qrf_pred[[3]],"YlOrBr","SOC [%]")
    ggsave(file.path(getwd(),output_dir,"FinalPrediction_QuantialMapSoc95Percentile.png"),
           plot_prediction, width = 6, height = 7, dpi = 300)
    
    # Plot prediction interval width (reversed color scheme: red = high uncertainty)
    plot_PIW <- PlotR(soc_qrf_interval_width,"Spectral","PIW (90%)",revers=TRUE,accuracy = 0.1)
    ggsave(file.path(getwd(),output_dir,"FinalPrediction_MapPIW90.png"),
           plot_PIW, width = 6, height = 7, dpi = 300)
    
    cat("Quantile prediction maps saved.\n")
    cat("\n")
    cat("=", rep("=", 79), "\n", sep = "")
    cat("QUANTILE REGRESSION FOREST PREDICTION COMPLETED\n")
    cat("=", rep("=", 79), "\n\n", sep = "")
  }
  
  # ============================================================================
  # DISTANCE-BASED UNCERTAINTY: DISSIMILARITY INDEX AND AREA OF APPLICABILITY
  # ============================================================================
  # Calculate Dissimilarity Index (DI) and Area of Applicability (AOA)
  # DI measures how different prediction locations are from training data
  # AOA identifies areas where predictions are reliable (within training data envelope)
  
  if(distance_uncertainty==TRUE){
    cat("=" , rep("=", 79), "\n", sep = "")
    cat("DISTANCE UNCERTAINTY QUANTIFICATION\n")
    cat("=" , rep("=", 79), "\n\n", sep = "")
    
    cat("Calculating Area of Applicability with Tile-based processing ...\n")
    
    # --------------------------------------------------------------------------
    # Step 1: Pre-compute trainDI threshold
    # --------------------------------------------------------------------------
    # trainDI calculates the dissimilarity threshold based on training data
    # This threshold defines the boundary of the Area of Applicability
    cat("Step 1: Computing trainDI threshold...\n")
    trainDI_result <- trainDI(final_model_full)
    
    # Save trainDI object for potential reuse
    saveRDS(trainDI_result, paste0(getwd(),output_dir, "trainDI_object.rds"))
    cat(paste("TrainDI threshold:", round(trainDI_result$threshold, 4), "\n\n"))
    
    # --------------------------------------------------------------------------
    # Step 2: Check if tiling is needed
    # --------------------------------------------------------------------------
    raster_cells <- ncell(covariates)
    cat(paste("Total raster cells:", format(raster_cells, big.mark = ","), "\n"))
    
    if(raster_cells > 5e6) {
      cat("\nLarge raster detected. Using tile-based processing...\n\n")
      
      # Create directories for tiles
      dir.create("tiles", showWarnings = FALSE)
      dir.create("aoa_tiles", showWarnings = FALSE)
      
      # Manual tile creation function (same as for QRF)
      create_tiles_manual <- function(raster, n_tiles_x = 4, n_tiles_y = 4,
                                      output_dir = "tiles") {
        
        dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
        ext_full <- ext(raster)
        x_min <- ext_full[1]
        x_max <- ext_full[2]
        y_min <- ext_full[3]
        y_max <- ext_full[4]
        
        tile_width <- (x_max - x_min) / n_tiles_x
        tile_height <- (y_max - y_min) / n_tiles_y
        
        cat(paste("Creating", n_tiles_x * n_tiles_y, "tiles\n"))
        
        tile_files <- character()
        for(i in 1:n_tiles_y) {
          for(j in 1:n_tiles_x) {
            tile_xmin <- x_min + (j - 1) * tile_width
            tile_xmax <- x_min + j * tile_width
            tile_ymin <- y_min + (i - 1) * tile_height
            tile_ymax <- y_min + i * tile_height
            
            tile_extent <- ext(tile_xmin, tile_xmax, tile_ymin, tile_ymax)
            tile_raster <- crop(raster, tile_extent)
            
            tile_filename <- file.path(output_dir,
                                       sprintf("tile_%02d_%02d.tif", i, j))
            writeRaster(tile_raster, tile_filename, overwrite = TRUE)
            tile_files <- c(tile_files, tile_filename)
          }
        }
        
        cat(paste("Created", length(tile_files), "tiles\n"))
        return(tile_files)
      }
      
      # Create tiles
      cat("Creating raster tiles...\n")
      tile_files <- create_tiles_manual(
        raster = covariates,
        n_tiles_x = n.tile,
        n_tiles_y = n.tile,
        output_dir = "tiles"
      )
      
      cat("\n")
      
      # Detect available cores
      ncores <- detectCores() - 1
      cat(paste("Using", ncores, "cores for parallel processing\n\n"))
      
      # Process tiles in parallel
      cat("Processing AOA for each tile...\n")
      start_time <- Sys.time()
      
      # Platform-specific parallel processing
      if(.Platform$OS.type == "unix") {
        # Unix/Linux/Mac: Use fork-based parallelization
        tile_results <- mclapply(seq_along(tile_files), function(i) {
          
          tile_path <- tile_files[i]
          tile_raster <- rast(tile_path)
          
          # Calculate AOA for tile
          aoa_tile <- aoa(
            newdata = tile_raster,
            trainDI = trainDI_result,
            LPD = FALSE,         # Don't calculate Local Point Density (faster)
            verbose = FALSE
          )
          
          # Save DI and AOA results
          tile_basename <- basename(tile_path)
          tile_name <- sub("\\.tif$", "", tile_basename)
          
          di_file <- file.path("aoa_tiles", paste0("DI_", tile_name, ".tif"))
          aoa_file <- file.path("aoa_tiles", paste0("AOA_", tile_name, ".tif"))
          
          writeRaster(aoa_tile$DI, di_file, overwrite = TRUE)
          writeRaster(aoa_tile$AOA, aoa_file, overwrite = TRUE)
          
          return(list(di = di_file, aoa = aoa_file))
        }, mc.cores = ncores)
        
      } else {
        # Windows: Use PSOCK cluster
        library(doParallel)
        library(foreach)
        
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
        
        # Export trainDI to cluster
        clusterExport(cl, c("trainDI_result"))
        
        # Load packages on workers
        clusterEvalQ(cl, {
          library(CAST)
          library(terra)
        })
        
        # Process tiles
        tile_results <- foreach(i = seq_along(tile_files),
                                .packages = c("CAST", "terra")) %dopar% {
                                  
                                  tile_path <- tile_files[i]
                                  tile_raster <- rast(tile_path)
                                  
                                  # Calculate AOA
                                  aoa_tile <- aoa(
                                    newdata = tile_raster,
                                    trainDI = trainDI_result,
                                    LPD = FALSE,
                                    verbose = FALSE
                                  )
                                  
                                  # Save results
                                  tile_basename <- basename(tile_path)
                                  tile_name <- sub("\\.tif$", "", tile_basename)
                                  
                                  di_file <- file.path("aoa_tiles", paste0("DI_", tile_name, ".tif"))
                                  aoa_file <- file.path("aoa_tiles", paste0("AOA_", tile_name, ".tif"))
                                  
                                  writeRaster(aoa_tile$DI, di_file, overwrite = TRUE)
                                  writeRaster(aoa_tile$AOA, aoa_file, overwrite = TRUE)
                                  
                                  list(di = di_file, aoa = aoa_file)
                                }
        
        stopCluster(cl)
      }
      
      end_time <- Sys.time()
      processing_time <- difftime(end_time, start_time, units = "mins")
      cat(paste("\nTile processing completed in",
                round(processing_time, 1), "minutes\n\n"))
      
      # Merge tiles into final rasters
      cat("Merging tiles into final rasters...\n")
      
      di_files <- list.files("aoa_tiles", pattern = "^DI_.*\\.tif$",
                             full.names = TRUE)
      aoa_files <- list.files("aoa_tiles", pattern = "^AOA_.*\\.tif$",
                              full.names = TRUE)
      
      # Load all tiles
      di_rasters <- lapply(di_files, rast)
      aoa_rasters <- lapply(aoa_files, rast)
      
      # Merge using mosaic (handles overlaps better than merge)
      di_merged <- do.call(mosaic, di_rasters)
      aoa_merged <- do.call(mosaic, aoa_rasters)
      
      # Export final AOA layers
      writeRaster(di_merged,
                  file.path(getwd(),output_dir, "FinalPrediction_DissimilarityIndex.tif"),
                  overwrite = TRUE)
      writeRaster(aoa_merged,
                  file.path(getwd(),output_dir, "FinalPrediction_AOAmask.tif"),
                  overwrite = TRUE)
      
      aoa_result <- list(DI = di_merged, AOA = aoa_merged)
      cat("Tile merging completed.\n\n")
      
      # Optional: Clean up temporary tile files
      # unlink("tiles", recursive = TRUE)
      # unlink("aoa_tiles", recursive = TRUE)
      
    } else {
      # Small raster - process directly without tiling
      cat("Processing AOA directly (no tiling needed)...\n")
      
      aoa_result <- aoa(
        newdata = covariates,
        trainDI = trainDI_result,
        LPD = FALSE,
        verbose = TRUE
      )
      
      # Export AOA layers
      writeRaster(aoa_result$DI,
                  file.path(getwd(),output_dir, "FinalPrediction_DissimilarityIndex.tif"),
                  overwrite = TRUE)
      writeRaster(aoa_result$AOA,
                  file.path(getwd(),output_dir, "FinalPrediction_AOAmask.tif"),
                  overwrite = TRUE)
      
      di_merged <- aoa_result$DI
      aoa_merged <- aoa_result$AOA
    }
    
    # Calculate AOA summary statistics
    cat("Calculating AOA summary statistics...\n")
    
    aoa_summary <- data.frame(
      Metric = c("Total_prediction_cells", "Cells_within_AOA",
                 "Percent_within_AOA", "Mean_DI", "Max_DI"),
      Value = c(
        sum(!is.na(values(aoa_merged))),
        sum(values(aoa_merged) == 1, na.rm = TRUE),
        round(100 * sum(values(aoa_merged) == 1, na.rm = TRUE) /
                sum(!is.na(values(aoa_merged))), 2),
        round(mean(values(di_merged), na.rm = TRUE), 4),
        round(max(values(di_merged), na.rm = TRUE), 4)
      )
    )
    
    # Export AOA summary
    write.csv(aoa_summary, file.path(getwd(),output_dir, "FinalPrediction_AoaSummaryStatistics.csv"),
              row.names = FALSE)
    
    cat("\n=== AOA COMPUTATION COMPLETED ===\n")
    print(aoa_summary)
    cat("\n")
    
    # Create and export DI visualization
    # Higher DI (red) = less similar to training data = higher uncertainty
    plot_prediction <- PlotR(di_merged,
                             "Spectral",
                             "DI",
                             revers = TRUE,
                             accuracy = 0.01,
                             ClassIntervallMethod = "quantile",
                             n_classes = 9)
    ggsave(file.path(getwd(),output_dir,"FinalPrediction_MapDI.png"),
           plot_prediction,
           width = 6,
           height = 7, dpi = 300)
  }
}

# ==============================================================================
# EXAMPLE USAGE
# ==============================================================================
# Define working directory (modify as needed)
working_dir <- "~/SOCastR/"  # Alternative: use hard-coded absolute path

# Run SOCastR workflow
SOCastR(working_dir,
        input_dir = "input",
        output_dir = "output",
        samples = "SAMPLES_LUCAS_EPSG25832.shp",
        covariates = "COVARIATS_EPSG25832.tif",
        soc_column = "SOC",
        n.tile = 4,
        model_uncertainty = TRUE,
        distance_uncertainty = TRUE)

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
