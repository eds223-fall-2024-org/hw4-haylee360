
#' Find Aquatic Species Suitable Habitat
#' 
#' @description
#' This function takes an aquatic species' sea surface temperature and depth requirements and visualizes where its suitable habitat is along the west coast of the United States.
#'
#' @param species str, the name of the species of interest
#' @param min_sst float, the minimum sea surface temperature of the species
#' @param max_sst float, the maximum sea surface temperature of the species
#' @param min_depth int, the minimum depth of the species
#' @param max_depth int, the maximum depth of the species
#'
#' @return A map of the suitable living areas of the species of interest
#' @export
#'
#' @examples
#' aqua_fun("Pacific Sea Star", 24.7, 29.3, -30, 0)

aqua_fun <- function(species, min_sst, max_sst, min_depth, max_depth){
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
  eez <- st_transform(eez, crs = crs("EPSG:4326"))

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

  # Use reclassification matrix to reclassify sst raster
  reclass_sst <- classify(mean_sst, rcl = rcl_sst)
  
  # Create depth reclassification matrix
  rcl_depth <- matrix(c(-Inf, min_depth, NA,
                        min_depth, max_depth, 1,
                        max_depth, Inf, NA),
                      ncol = 3, byrow = TRUE)

  # Use reclassification matrix to reclassify depth raster
  reclass_depth <- classify(depth_resample, rcl = rcl_depth)

  # Find locations that satisfy both sst and depth conditions
  sst_depth <- lapp(c(reclass_sst, reclass_depth), fun = "*")

  # Find area of suitable cells
  suitable_area <- cellSize(x = sst_depth, # places where sst and depth match
                            mask = TRUE, # keeps NA values
                            unit = 'km') # Selecting km from data
  # Rasterize eez
  eez_raster <- rasterize(eez,
                          suitable_area,
                          field = 'rgn')

  # Use zonal algebra to sum the suitable area by region
  eez_suitable <- zonal(x = suitable_area, 
                        z = eez_raster, # Raster representing zones
                        fun = 'sum', # To add up total area
                        na.rm = TRUE)
  
  # Join area with geometry to map later on
  eez <- left_join(eez, eez_suitable, by = join_by(rgn))

  # Plot regions and suitable areas together
  tm_shape(depth) +
    tm_raster(palette = "-GnBu",
              title = "Bathymetry\n(m above and below sea level)",
              # alpha = 0.5,
              midpoint = 0,
              legend.show = TRUE) +
    tm_shape(eez, raster.downsample = TRUE) +
    tm_polygons(col = "area",
                palette = "Reds",
                alpha = 0.8,
                linewidth = 0.2,
                title = expression("Suitable habitat area (km"^2*")")) +
    tm_text("rgn", size = 0.45) +
    tm_compass(size = 1,
               position = c("left", "bottom")) +
    tm_scale_bar(position = c("left", "bottom")) +
    tm_layout(legend.outside = TRUE,
              frame = FALSE,
              main.title = paste0(species, " Suitable Habitat in\nWest Coast Exclusive Economic Zones"))
  }
