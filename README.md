# SOCopt
Spatial SOC prediction considering spatial dependencies and uncertainties 

  working_dir <- "D:/Dropbox/GIT/SOCastR/"
  
  source(file.path(working_dir,"SOCastR.R"))
  
  SOCastR(working_dir,
          input_dir = "input",
          output_dir = "output",
          samples = "SAMPLES_EPSG25832.shp",
          covariates = "COVARIATS_EPSG25832.tif",
          soc_column = "SOC",
          n.tile = 4)
  
