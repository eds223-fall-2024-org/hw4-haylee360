

# minimum and maximum sea surface temperature
# minimum and maximum depth
# species name

aqua_fun <- function(species, max_sst, min_sst, max_depth, min_depth){
  
  # Read in EEZ and Depth data
  eez <- read_sf(here("data", "wc_regions_clean.shp"))
  depth <- terra::rast(here("data", "depth.tif"))
  
  # Create list of sst files
  sst_files <- list.files(here("data"), 
                          pattern = "average", 
                          full.names = TRUE)
  # Read in sst rasters 
  sst <- c(rast(sst_files))
  
  # Reproject CRS to EPSG:4326
  depth <- project(depth, crs("EPSG:4326"))
  sst <- project(sst, crs("EPSG:4326"))
  eez <- project(eez, crs("EPSG:4326"))

  # Find the mean sst and convert to kelvin
  mean_sst <- app(sst, fun = mean)
  mean_sst <- mean_sst - 273.15
  
  # Crop the depth raster to match the extent of the sst
  depth_crop <- crop(depth, sst)
  
  # Resample depth to match sst using nearest neighbor method
  depth_resample <- resample(depth_crop, y = sst, method = "near")
  
  # Create sst reclassification matrix
  rcl_sst <- matrix(c(-Inf, min_sst, NA, 
                      min_sst, max_sst, 1, 
                      max_sst, Inf, NA),
                    ncol = 3, byrow = TRUE)
  
  # use reclassification matrix to reclassify sst raster
  reclass_sst <- classify(mean_sst, rcl = rcl_sst)
  
  # create depth reclassification matrix
  rcl_depth <- matrix(c(-Inf, min_depth, NA, 
                        min_depth, max_depth, 1, 
                        max_depth, Inf, NA),
                      ncol = 3, byrow = TRUE)
  
  # use reclassification matrix to reclassify depth raster
  reclass_depth <- classify(depth_resample, rcl = rcl_depth)
  
  # Find locations that satisfy both sst and depth conditions
  sst_depth <- lapp(c(reclass_sst, reclass_depth), fun = "*")
  
  # Rasterize eez
  eez_raster <- rasterize(eez, sst_depth, field = "rgn")
  
  # Clip by eez raster
  sst_depth_eez <- sst_depth[eez_raster, drop = FALSE]

  # Find area of grid cells in suitable locations
  area_cell <- cellSize(sst_depth_eez, unit = "km")
  
  eez_raster <- rasterize(west_eez,
                          area_cell,
                          field = 'rgn')

  
  # Use zonal algebra to aggregate a grouping variable
  eez_suitable <- zonal(x = area_cell, 
                        z = eez_raster, # Raster representing zones
                        fun = 'sum', # To add up total area
                        na.rm = TRUE)
  
  }


# process SST (mean, Cel)
# crop, resample
# reclassify
# sst depth multiplication
# rasterize eez
# crop and mask to eez
# find area of grid cells
# find total suitable area