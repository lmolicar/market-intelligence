# createROIMask.R
#
# Purpose:
#   Creates a mask for the region of interest and saves it to a GeoTIFF file.
#
# Record of revisions:
#     Date         Programmer                  Description of change
# ===========  ==================  =============================================
# Mar-01-2024  Luis Molina         Command line version.
#                                  Reads in the following parameters:
#                                  - the full path to the vector file will be
#                                   used to define the AOI;
#                                  - the file that defines the output grid.
#                                  TODO: include a helper function.        
#
# Jul-08-2023  Luis Molina         Reads in GeoJson vector file.
#
# Dec-2022     H. Achicanoy        Original code.
#
#


# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(if(!require(pacman)){install.packages('pacman');library(pacman)} else {library(pacman)})
suppressMessages(pacman::p_load(tidyverse,raster,terra,sp,geodata,rnaturalearthdata,rnaturalearth))

args <- commandArgs(trailingOnly = TRUE)

# Path to the geojson file that will define the extent of the area of interest
path_to_geojson <- args[1]

# Name of the file template that will be used to define the output grid
grid_template <- args[2]

filename <- tail(strsplit(path_to_geojson, "/")[[1]], 1)
roiName <- strsplit(filename, ".", fixed = TRUE)[[1]][1]

## Raster
tif  <- paste0('../ROI/', roiName, '.tif')

if(!file.exists(tif)){
    roi <- sf::read_sf(path_to_geojson)
    
    # Use temperature data as template (0.1 x 0.1 deg)
    ref <- raster(grid_template, lyrs = 1)
    
    # Crop the template raster to the extent of the ROI, but increase the border
    # by the spatial resolution in both directions
    roi_extent <- raster::extent(roi)
    minx <- xmin(roi_extent) - xres(ref)
    maxx <- xmax(roi_extent) + xres(ref)
    miny <- ymin(roi_extent) - yres(ref)
    maxy <- ymax(roi_extent) + yres(ref)
    roi_extent_larger <- extent(c(minx, maxx, miny, maxy))
    ref <- ref %>% raster::crop(roi_extent_larger)
    ref <- ref * 0 + 1

    raster::writeRaster(ref * 0 + 1, tif)
}
