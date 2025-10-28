# =============================================================================
# TILE-BASED PARALLELIZED QUANTILE REGRESSION FOREST PREDICTION
# Cross-platform solution (Windows & Linux) for terra SpatRaster
# Similar approach to CAST AOA tile processing
# =============================================================================

if(model_uncertainty==TRUE){
  
  cat("=", rep("=", 79), "\n", sep = "")
  cat("MODEL UNCERTAINTY: QUANTILE REGRESSION FOREST UNCERTAINTY LAYERS (TILE-BASED)\n")
  cat("=", rep("=", 79), "\n\n", sep = "")
  
  # Prepare training data for QRF
  train_x <- model_data_clean[, ffs_model$selectedvars, drop = FALSE]
  train_y <- model_data_clean$SOC
  
  # Fit quantile regression forest model
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
  
  # Export training results
  qrf_model.accuracy <- data.frame(
    R_squared = qrf_model.lm.R2,
    RMSE = qrf_model.lm.RMSE,
    Method = "Quantile_Regression_Forest"
  )
  
  write.csv(qrf_model.accuracy,
            file.path(getwd(), output_dir, "QuantilePrediction_Accuracy.csv"),
            row.names = FALSE)
  cat("QRF model training completed.\n\n")
  
  # =============================================================================
  # ENSURE RASTER LAYER NAMES MATCH MODEL VARIABLES
  # =============================================================================
  
  cat("Checking raster layer names...\n")
  cat("Model expects:", paste(ffs_model$selectedvars, collapse = ", "), "\n")
  cat("Raster has:", paste(names(covariates), collapse = ", "), "\n\n")
  
  missing_vars <- setdiff(ffs_model$selectedvars, names(covariates))
  if (length(missing_vars) > 0) {
    stop("ERROR: Missing variables in covariates: ", paste(missing_vars, collapse = ", "))
  }
  
  # Subset and reorder covariates to match model
  covariates <- covariates[[ffs_model$selectedvars]]
  cat("Layer names matched successfully.\n\n")
  
  # =============================================================================
  # TILE-BASED PROCESSING SETUP
  # =============================================================================
  
  # Define quantiles for prediction intervals
  quantiles <- c(0.05, 0.5, 0.95)
  
  # Check if tiling is needed (similar to AOA approach in your script)
  raster_cells <- ncell(covariates)
  cat(paste("Total raster cells:", format(raster_cells, big.mark = ",")), "\n")
  
  # Determine if tiling is needed (>5 million cells)
  use_tiling <- raster_cells > 5e6
  
  if (use_tiling) {
    cat("Large raster detected. Using tile-based processing...\n\n")
    
    # =============================================================================
    # CREATE TILES FUNCTION (from your SOCastR.R)
    # =============================================================================
    
    create_tiles_manual <- function(raster, ntiles_x = 4, ntiles_y = 4, output_dir = "tiles") {
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      
      ext_full <- ext(raster)
      xmin <- ext_full[1]
      xmax <- ext_full[2]
      ymin <- ext_full[3]
      ymax <- ext_full[4]
      
      tile_width <- (xmax - xmin) / ntiles_x
      tile_height <- (ymax - ymin) / ntiles_y
      
      cat(paste("Creating", ntiles_x * ntiles_y, "tiles\n"))
      
      tile_files <- character()
      
      for (i in 1:ntiles_y) {
        for (j in 1:ntiles_x) {
          tile_xmin <- xmin + (j - 1) * tile_width
          tile_xmax <- xmin + j * tile_width
          tile_ymin <- ymin + (i - 1) * tile_height
          tile_ymax <- ymin + i * tile_height
          
          tile_extent <- ext(tile_xmin, tile_xmax, tile_ymin, tile_ymax)
          tile_raster <- crop(raster, tile_extent)
          
          tile_filename <- file.path(output_dir, sprintf("tile_%02d_%02d.tif", i, j))
          writeRaster(tile_raster, tile_filename, overwrite = TRUE)
          tile_files <- c(tile_files, tile_filename)
        }
      }
      
      cat(paste("Created", length(tile_files), "tiles\n"))
      return(tile_files)
    }
    
    # =============================================================================
    # CREATE TILES
    # =============================================================================
    
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
    
    # =============================================================================
    # PARALLEL TILE PROCESSING
    # =============================================================================
    
    library(parallel)
    library(doParallel)
    library(foreach)
    
    n_cores <- detectCores() - 1
    start_time <- Sys.time()
    
    cat(paste("\nProcessing", length(tile_files), "tiles across", n_cores, "cores...\n\n"))
    
    # Platform-specific parallel processing
    if (.Platform$OS.type == "windows") {
      # Windows: PSOCK cluster
      cl <- makeCluster(n_cores, type = "PSOCK")
      registerDoParallel(cl)
      
      # Export necessary objects
      clusterExport(cl, varlist = c("qrf_model", "quantiles", "quantile_tiles_dir"), 
                    envir = environment())
      clusterEvalQ(cl, {
        library(terra)
        library(quantregForest)
      })
      
      # Process tiles in parallel
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
      # Linux/Unix: Fork or mclapply
      cat("Unix/Linux detected: Using fork-based parallelization\n")
      
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
      }, mc.cores = n_cores)
    }
    
    end_time <- Sys.time()
    processing_time <- difftime(end_time, start_time, units = "mins")
    cat(paste("\nTile processing completed in", round(processing_time, 1), "minutes\n\n"))
    
    # =============================================================================
    # MERGE TILES INTO FINAL RASTERS
    # =============================================================================
    
    cat("Merging tiles into final quantile rasters...\n")
    
    # For each quantile, merge its tiles
    quantile_layers <- list()
    
    for (q_idx in seq_along(quantiles)) {
      q_pattern <- paste0("Q", quantiles[q_idx]*100, "_tile_.*\\.tif$")
      q_files <- list.files(quantile_tiles_dir, pattern = q_pattern, full.names = TRUE)
      
      cat(paste("Merging", length(q_files), "tiles for quantile", quantiles[q_idx], "...\n"))
      
      # Load all tiles for this quantile
      q_rasters <- lapply(q_files, rast)
      
      # Merge using mosaic (more robust than merge)
      q_merged <- do.call(mosaic, q_rasters)
      names(q_merged) <- paste0("SOC_q", quantiles[q_idx]*100)
      
      quantile_layers[[q_idx]] <- q_merged
    }
    
    # Combine all quantiles into one stack
    soc_qrf_pred <- rast(quantile_layers)
    
    cat("Tile merging completed.\n\n")
    
    # Clean up tile directories (optional)
    # unlink(tiles_dir, recursive = TRUE)
    # unlink(quantile_tiles_dir, recursive = TRUE)
    
  } else {
    # =============================================================================
    # SMALL RASTER: PROCESS DIRECTLY WITHOUT TILING
    # =============================================================================
    
    cat("Small raster detected. Processing directly (no tiling needed)...\n\n")
    
    pb <- txtProgressBar(min = 0, max = length(quantiles), style = 3)
    
    quantile_list <- list()
    for (i in seq_along(quantiles)) {
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
    
    soc_qrf_pred <- rast(quantile_list)
    cat("\nQuantile predictions completed.\n\n")
  }
  
  # =============================================================================
  # POST-PROCESSING (SAME FOR BOTH APPROACHES)
  # =============================================================================
  
  # Calculate prediction interval width
  soc_qrf_interval_width <- soc_qrf_pred[[3]] - soc_qrf_pred[[1]]  # Q95 - Q05
  names(soc_qrf_interval_width) <- "SOC_PI_width"
  
  # Combine all layers
  soc_qrf_stack <- c(soc_qrf_pred, soc_qrf_interval_width)
  
  # Export final rasters
  writeRaster(soc_qrf_stack,
              file.path(getwd(), output_dir, "FinalPrediction_SocQuantileLayers.tif"),
              overwrite = TRUE)
  cat("Quantile prediction raster layers saved.\n")
  
  # Summary statistics
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
  
  # Export plots
  plot_prediction <- PlotR(soc_qrf_pred[[1]],"YlOrBr","SOC [%]")
  ggsave(file.path(getwd(),output_dir,"FinalPrediction_MapQuantialSoc5Percentile.png"),
         plot_prediction, width = 6, height = 7, dpi = 300)
  
  plot_prediction <- PlotR(soc_qrf_pred[[2]],"YlOrBr","SOC [%]")
  ggsave(file.path(getwd(),output_dir,"FinalPrediction_MapQuantialSoc50Percentile.png"),
         plot_prediction, width = 6, height = 7, dpi = 300)
  
  plot_prediction <- PlotR(soc_qrf_pred[[3]],"YlOrBr","SOC [%]")
  ggsave(file.path(getwd(),output_dir,"FinalPrediction_QuantialMapSoc95Percentile.png"),
         plot_prediction, width = 6, height = 7, dpi = 300)
  
  plot_PIW <- PlotR(soc_qrf_interval_width,"Spectral","PIW (90%)",revers=TRUE)
  ggsave(file.path(getwd(),output_dir,"FinalPrediction_MapPIW90.png"),
         plot_PIW, width = 6, height = 7, dpi = 300)
  
  cat("Quantile prediction maps saved.\n")
  
  cat("\n")
  cat("=", rep("=", 79), "\n", sep = "")
  cat("QUANTILE REGRESSION FOREST PREDICTION COMPLETED\n")
  cat("=", rep("=", 79), "\n\n", sep = "")
}

# =============================================================================
# KEY FEATURES OF TILE-BASED APPROACH:
# =============================================================================
# 1. Automatic detection: Uses 5M cells threshold (like AOA processing)
# 2. Tile creation: Uses your existing create_tiles_manual function
# 3. Parallel processing: Each tile processed independently on separate core
# 4. Memory efficient: Only loads one tile at a time per worker
# 5. Cross-platform: PSOCK for Windows, mclapply for Linux/Unix
# 6. Robust merging: Uses terra::mosaic instead of merge
# 7. Progress tracking: Reports processing time
# 8. Tile cleanup: Optional - can delete intermediate tiles after merging
# =============================================================================