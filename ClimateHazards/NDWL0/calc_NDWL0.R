#
#
# Program:
#         calc_NDWL.R
#
#
# Record of changes
#
#    Date         Programmer                 Description of changes
# ==========  ==================  ==============================================
# 2024-06-21  Luis Molina         Each input climate dataset is now read from a
#                                 single NetCDF file.
#
#

## Number of waterlogging days (NDWL0)
## By: H. Achicanoy
## December, 2022

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))
source('C:/Users/luismolina/Documents/OneSF/ANALYZER/NDWL0/CODE/_main_functions.R')

root <- 'E:/luismolina/ImpactSF/ANALYZER/NDWL0/LAOS'

roiName <- 'gadm41_LAO_0'
ref <- terra::rast(paste0(root,'/ROI/', roiName,'.tif'))

# Soil variables
scp <- terra::rast(paste(root,'SOILS', paste0(roiName, '_scp.tif'), sep = '/', collapse = ""))
scp <- scp %>% terra::resample(ref) %>% terra::mask(ref)
sst <- terra::rast(paste(root,'SOILS', paste0(roiName, '_ssat.tif'), sep = '/', collapse = ""))
sst <- sst %>% terra::resample(ref) %>% terra::mask(ref)

# Calculate NDWL function
calc_ndwl0 <- function(yr, mn){
  outfile <- paste0(root, '/OUTPUT','/NDWL0-',yr,'-',mn,'.tif')
  cat(outfile, "\n")
  if(!file.exists(outfile)){
    dir.create(dirname(outfile),F,T)
    # Last day of the month
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    # Sequence of dates
    if(as.numeric(yr) > 2020 & mn == '02'){
      dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-28')), by = 'day')
    } else {
      dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    }

    if(as.numeric(yr) > 2050){
      dts_chr <- as.character(dts)
      dts_chr <- gsub(pattern = yr, replacement = yrs_mpg$Baseline[yrs_mpg$Future == yr], x = dts_chr)
      dts <- as.Date(dts_chr); rm(dts_chr)
    }

    # Select the layers
    sel_layers <- numLayer[lubridate::year(timestamps) == as.numeric(yr) &
                           lubridate::month(timestamps) == as.numeric(mn)]
    
    prc <- terra::rast(pr_fls, lyrs = sel_layers)
    prc <- prc %>% terra::resample(ref) %>% terra::mask(ref)
    prc[prc == -9999] <- NA
    
    tmx <- terra::rast(tx_fls, lyrs = sel_layers)
    tmx <- tmx %>% terra::resample(ref) %>% terra::mask(ref)
    tmx[tmx == -9999] <- NA
    
    tmn <- terra::rast(tm_fls, lyrs = sel_layers)
    tmn <- tmn %>% terra::resample(ref) %>% terra::mask(ref)
    tmn[tmn == -9999] <- NA
    
    tav <- (tmx + tmn)/2
    
    srd <- terra::rast(sr_fls, lyrs = sel_layers)
    srd <- srd %>% terra::resample(ref) %>% terra::mask(ref)
    srd <- srd/1000000

    
    # Maximum evapotranspiration
    ETMAX <- terra::lapp(x = terra::sds(srd,tmn,tav,tmx), fun = peest)
    rm(list=c("tmn", "tmx", "tav", "srd"))
    gc(verbose=F, full=T, reset=T)
    
    # Compute water balance model
    date <- paste0(yr,'-',mn)
    if(date %in% c('2023-01')){
      AVAIL <<- ref
      AVAIL[!is.na(AVAIL)] <- 0
    } else {
      AVAIL <<- terra::rast(paste0(dirname(outfile),'/AVAIL.tif'))
    }
    
    eabyep_calc <- function(soilcp = scp, soilsat = ssat, avail = AVAIL,rain = prc[[1]], evap = ETMAX[[1]]){
      
      avail   <- min(avail, soilcp)
      
      # ERATIO
      percwt <- min(avail/soilcp*100, 100)
      percwt <- max(percwt, 1)
      eratio <- min(percwt/(97-3.868*sqrt(soilcp)), 1)
      
      demand  <- eratio * evap
      result  <- avail + rain - demand
      logging <- result - soilcp
      logging <- max(logging, 0)
      logging <- min(logging, soilsat)
      # runoff  <- result - logging + soilcp
      avail   <- min(soilcp, result)
      avail   <- max(avail, 0)
      # runoff  <- max(runoff, 0)
      
      out     <- list(Availability = c(AVAIL, avail),
                      # Demand       = demand,
                      # Eratio       = eratio,
                      Logging      = logging
      )
      return(out)
    }
    
    watbal <- 1:terra::nlyr(ETMAX) %>%
      purrr::map(.f = function(i){
        water_balance <- eabyep_calc(soilcp  = scp,
                                     soilsat = sst,
                                     avail   = AVAIL[[terra::nlyr(AVAIL)]],
                                     rain    = prc[[i]],
                                     evap    = ETMAX[[i]])
        AVAIL <<- water_balance$Availability
        return(water_balance)
      })
    LOGGING <- watbal %>% purrr::map('Logging') %>% terra::rast()
    # Calculate number of soil waterlogging days (if logging is above 0)
    # Note NDWL50 uses ssat * 0.5 (so this means soil is at 50% toward saturation)
    # here it suffices if the soil is above field capacity
    NDWL0  <- sum(LOGGING > 0)
    terra::writeRaster(NDWL0, outfile)
    terra::writeRaster(AVAIL[[terra::nlyr(AVAIL)]], paste0(dirname(outfile),'/AVAIL.tif'), overwrite = T)
    
    #clean up
    rm(list=c("prc", "ETMAX", "AVAIL", "watbal", "ERATIO", "LOGGING", "NDWL0"))
    gc(verbose=F, full=T, reset=T)
  }
}

# # Historical setup
yrs <- 2023:2027
mns <- c(paste0('0',1:9),10:12)
stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
names(stp) <- c('yrs','mns')
stp <- stp %>%
  dplyr::arrange(yrs, mns) %>%
  base::as.data.frame()


climate_path <- "E:/luismolina/ImpactSF/ANALYZER/NDWL0/LAOS/CLIMATE"
pr_fls <- paste0(climate_path, '/', "pr_bias-corrected_2021-2050.nc")
tx_fls <- paste0(climate_path, '/', "tasmax_day_EC-Earth3P-HR_highres-future_r1i1p2f1_gr_20210101-20501231_Laos.nc")
tm_fls <- paste0(climate_path, '/', "tasmin_day_EC-Earth3P-HR_highres-future_r1i1p2f1_gr_20210101-20501231_Laos.nc")
sr_fls <- paste0(climate_path, '/', "srd_bias-corrected_2021-2050.nc")


# Normally input datasets cover the same time period. We will read precipitation
# and use the time stamp to select the layer

tmp_ds <- terra::rast(pr_fls)
timestamps <- terra::time(tmp_ds)
numLayer <- 1:length(timestamps)
rm(tmp_ds)


#out_dir <- paste0(root,'/NDWL0')

1:nrow(stp) %>%
  purrr::map(.f = function(i){
    calc_ndwl0(yr = stp$yrs[i], mn = stp$mns[i]); gc(verbose=F, full=T, reset=T)
    if (i%%5 == 0) {
      tmpfls <- list.files(tempdir(), full.names=TRUE)
      1:length(tmpfls) %>% purrr::map(.f = function(k) {system(paste0("rm -f ", tmpfls[k]))})
    }
  })


# Future setup

# gcm <- 'INM-CM5-0' #'ACCESS-ESM1-5' 'MPI-ESM1-2-HR' 'EC-Earth3' 'INM-CM5-0' 'MRI-ESM2-0'
# #for (gcm in c("ACCESS-ESM1-5", "MPI-ESM1-2-HR", "EC-Earth3", "INM-CM5-0", "MRI-ESM2-0")) {
#     for (ssp in c('ssp245', 'ssp585')) {
#         for (prd in c('2021_2040', '2041_2060')) {
#             #ssp <- 'ssp245'
#             #prd <- '2021_2040'
#             
#             cmb <- paste0(ssp,'_',gcm,'_',prd)
#             prd_num <- as.numeric(unlist(strsplit(x = prd, split = '_')))
#             yrs <- prd_num[1]:prd_num[2]
#             mns <- c(paste0('0',1:9),10:12)
#             stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
#             names(stp) <- c('yrs','mns')
#             stp <- stp %>%
#               dplyr::arrange(yrs, mns) %>%
#               base::as.data.frame()
#             pr_pth <- paste0(root,'/chirps_cmip6_africa/Prec_',gcm,'_',ssp,'_',prd) # Precipitation
#             tx_pth <- paste0(root,'/chirts_cmip6_africa/Tmax_',gcm,'_',ssp,'_',prd) # Maximum temperature
#             tm_pth <- paste0(root,'/chirts_cmip6_africa/Tmin_',gcm,'_',ssp,'_',prd) # Minimum temperature
#             sr_pth <- paste0(root,'/ecmwf_agera5/solar_radiation_flux')             # Solar radiation
#             out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/',cmb,'/NDWL0')
# 
#             yrs_mpg <- data.frame(Baseline = as.character(rep(1995:2014, 2)),
#                                   Future = as.character(c(2021:2040,2041:2060)))
# 
#             1:nrow(stp) %>%
#               purrr::map(.f = function(i){
#                 calc_ndwl0(yr = stp$yrs[i], mn = stp$mns[i]); gc(verbose=F, full=T, reset=T)
#                 if (i%%5 == 0) {
#                   tmpfls <- list.files(tempdir(), full.names=TRUE)
#                   1:length(tmpfls) %>% purrr::map(.f = function(k) {system(paste0("rm -f ", tmpfls[k]))})
#                 }
#               })
#         }
#     }
#}
