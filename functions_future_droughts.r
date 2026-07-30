#' Compute Potential Evapotranspiration (PET) using the Penman-Monteith equation and SPEI
#'
#' This script contains the functions to compute the Potential Evapotranspiration (PET) using the Penman-Monteith 
#' equation and the Standardized Precipitation-Evapotranspiration Index (SPEI) for historical and future climate 
#' conditions.
#' 
#' The script is divided into four main sections:
#' 
#' 1. Evapotranspiration functions: This section contains the functions to compute the reference 
#' evapotranspiration (ETo) using the Penman-Monteith equation `eto_pixel` and the potential evapotranspiration (PET) 
#' using the combined method of Priestley-Taylor and forest conditions `pet_pixel`. it contains the wrapper function
#' `compute_et` to compute the ETo or PET for a given set of input rasters (temperature, solar radiation, relative 
#' humidity, and wind speed). which outouts a SpatRaster with the computed ETo or PET values.
#' 
#' 2. SPEI functions: reference period parameters: This section contains the functions to compute the parameters of
#' the log-logistic distribution from the historical data (precipitation and potential evapotranspiration) for a given 
#' reference period. It contains the wrapper function `compute_params` to compute the parameters for all pixels in the 
#' input rasters and output a SpatRaster with the computed parameters.
#' 
#' 3. SPEI functions: computing SPEI: This section contains the functions to compute the Standardized Precipitation-
#' Evapotranspiration Index (SPEI) for a given set of input rasters (precipitation, potential evapotranspiration,
#' and parameters from the reference period). It contains the wrapper function `compute_spei` to compute the SPEI
#'  for all pixels in the input rasters and output a SpatRaster with the computed SPEI values.
#' 
#' 4. Example of usage: This section contains an example  of how to use the functions defined in the previous sections to 
#' compute PET and SPEI for historical and future climate data. Substituting the names of files with real data should
#' allow the user to compute PET and SPEI for their own data, as long as the input rasters have the same extent and temporal 
#' resolution.
#' 
#' Author: Pedro Alencar
#' Date: 27.07.26


#%% 0. Load libraries --------

# data manipulation
library(tidyr)
library(dplyr)
library(lubridate)
library(stringr)
library(runner)

library(sf)
library(terra)
terra::terraOptions(progress = 1)   # show progress bar/messages

library(parallel)

# visualization
library(ggplot2)

# spei
library(SPEI)

# get httpgd running for plots (only if you are using VSCODE and want to plot in your BROWSER)
httpgd::hgd()
httpgd::hgd_browse()


#%% 1. evapotranspiration functions --------

########################################################################################
#' NOTES:
#' - to compute potential evapotranspiration, 4 variables are needed: temperature, solar radiation,
#' relative humidity and wind speed. They should have teh same extent and temporal resolution
#' 
#' - The units must be: wind in m/s, temperature in degrees Celsius, relative humidity in absolute (0-1), and 
#' solar radiation in W/m2 or MJ/m2/day. If the solar radiation is in W/m2, the function will convert it to MJ/m2/day 
#' (attention to the input `rsds_in_w`` parameter, which is TRUE by default, meaning that the input is in W/m2, if it is
#'  in MJ/m2/day, set it to FALSE)
#' 
#' - The output will be in mm/day, which is the unit used by the SPEI package to compute the index. 
#' 
#' - I made two functions to compute the evapotranspiration: one for the reference evapotranspiration (ETo) using 
#' the Penman-Monteith equation, and another for the potential evapotranspiration (PET) using the combined method
#'  of Priestley-Taylor and forest conditions. The ETo is more widely accepted (and I'd go for it), but the PET is 
#' fitted to the problem of evaporation in forests. In case you want to read more on this check:
#'      * https://doi.org/10.5194/hess-17-1331-2013
#'      * https://doi.org/10.1016/j.agwat.2020.106043
########################################################################################

#---------------------------------------------------------
# some code to be ignored, but I kept bcs I used it to get some example data
#---------------------------------------------------------
# example_path <- '/Volumes/PH_SSD_2/Brandenburg/BB-Projections/rcp85/ICHEC-EC-EARTH_r1i1p1_KNMI-RACMO22E'

# wind_file <- list.files(paste0(example_path, '/', 'sfcWind'), full.names = TRUE)[31]
# tas_file <- list.files(paste0(example_path, '/', 'tas'), full.names = TRUE)[31]
# rsds_file <- list.files(paste0(example_path, '/', 'rsds'), full.names = TRUE)[31]
# hurs_file <- list.files(paste0(example_path, '/', 'hurs'), full.names = TRUE)[31]

# wind <- rast(wind_file) # m.s-1
# tas <- rast(tas_file) # °C
# rsds <- rast(rsds_file) # W.m-2
# hurs <- rast(hurs_file) # absolute (0-1)


#' function to compute reference evapotranspiration (ETo) using the Penman-Monteith equation 
#' FAO 56 (Allen et al., 1998). Ther reference evapotranspiration is computed for a hypothetical
#'  grass reference crop with a height of 0.12 m, a fixed surface resistance of 70 s m-1, and an albedo of 0.23.
eto_pixel <- function(tas, rsds, hurs, sfcw, rsds_in_w = TRUE){

    # tas = -3.6
    # rsds = 26.3
    # hurs = 0.726
    # wind = 4.0
    # rsds_in_w = TRUE

    e_s <- 0.6108 * exp((17.27 * tas) / (tas + 237.3))
    delta <- (4098 * e_s) / ((tas + 237.3)^2)
    gamma <- 0.0668 # psychometric constant in kPa °C-1 for 101.3 kPa (sea level)

    if (rsds_in_w) {
        rsds <- rsds * 0.0864 # convert W.m-2 to MJ.m-2.day-1
        # cat('a')
    }

   eto <- (0.408 * delta * rsds + gamma * (900 / (tas + 273)) * sfcw * e_s * (1 - hurs)) / (delta + gamma * (1 + 0.34 * sfcw))

   return(eto)
}

#' function to compute potential evapotranspiration (PET) using the combined method of Priestley-Taylor
#'  (Priestley and Taylor, 1972) and forest conditions (Chow, 1988, table 2.8.2)
#' 
#' Note: althogh this approach is fitted to forest conditions, the use of Penman-Monteith widely recommended and well
#' accepted for multiple usages. The objective is also not to compute the evapotranspiration, but to compute SPEI, 
#' which is a drought index that is based on the difference between precipitation and potential evapotranspiration. Therefore, 
#' the absolute values of PET are not as important as the relative anomaly values, which are used to compute the index.
pet_pixel <- function(tas, rsds, hurs, sfcw, rsds_in_w = TRUE){

    # tas = -3.6
    # rsds = 26.3
    # hurs = 0.726
    # wind = 4.0
    # rsds_in_w = TRUE

    e_s <- 0.6108 * exp((17.27 * tas) / (tas + 237.3))
    delta <- (4098 * e_s) / ((tas + 237.3)^2)
    gamma <- 0.0668 # psychometric constant in kPa °C-1 for 101.3 kPa (sea level)

    if (rsds_in_w) {
        rsds <- rsds * 0.0864 # convert W.m-2 to MJ.m-2.day-1
        # cat('a')
    }

   pet <- (0.04 / (delta + gamma)) * (delta * rsds + gamma * sfcw * e_s * (1 - hurs))
   return(pet)
}


#---------------------------------------------------------
# Wrapper for SpatRaster inputs
#---------------------------------------------------------
compute_et <- function(tas, rsds, wind, hurs,
                       rsds_in_w = TRUE,
                       cores = parallel::detectCores() - 2,
                       filename = "files/test.nc",
                       overwrite = TRUE) {

  # Check inputs
  stopifnot(
    compareGeom(tas, rsds),
    compareGeom(tas, wind),
    compareGeom(tas, hurs)
  )

  terra::lapp(
    x = sds(tas, rsds, hurs, wind),
    fun = eto_pixel,
    rsds_in_w = rsds_in_w,
    cores = cores
  )
}

#---------------------------------------------------------
# Test
#---------------------------------------------------------
# test <- compute_et(tas = tas[[1]], rsds = rsds[[1]], wind = wind[[1]], hurs = hurs[[1]])
# plot(test)

#%% 2. spei functions: reference period paramenters --------

########################################################################################
#' NOTES: 
#' - Precipitation and evapotranspiration must have the same extent
#' - Both variables must be in mm/month (only positive values!)
#' - Make sure to save the parameters in a raster, as it will be reused several times
#'########################################################################################


#---------------------------------------------------------
# some code to be ignored, but I kept bcs I used it to get some example data
#---------------------------------------------------------
# example_path <- '/Volumes/KIRecover_1/Datasets/ECMWF/era5land_month/5d5501409f6b61a2857c669abfd99a07.nc'
# test <-rast(example_path) 

# # select period from 1970-1999
# test <- test[[which(year(as.Date(time(test), format = "X%Y.%m.%d")) >= 1970 & year(as.Date(time(test), format = "X%Y.%m.%d")) <= 1999)]]
# tp  <- test['tp'] 
# pev <- test['pev']

# # mask with shapefile
# shape_germany <- sf::read_sf('data/GIS/NUTS/NUTS_RG_20M_2021_4326.shp') |> filter(NUTS_ID == 'DE') 

# tp <- terra::mask(tp, shape_germany)
# terra::writeCDF(tp, 'files/tp_1970_1999.nc', overwrite = TRUE, compression = 5)

# pev <- terra::mask(pev, shape_germany)
# terra::writeCDF(pev, 'files/pev_1970_1999.nc', overwrite = TRUE, compression = 5)


# tp_de <- terra::rast('files/tp_1970_1999.nc')
# pev_de <- terra::rast('files/pev_1970_1999.nc')

# # trim NA!
# tp_de <- trim(tp_de)

# berlin <- terra::cellFromXY(tp_de, cbind(13.404954, 52.520008)) # Berlin
# tp_de[berlin]


# valid <- app(tp_de[[1]], fun = function(x) as.integer(any(!is.na(x))))

# source('/Users/alencar/Library/CloudStorage/OneDrive-Personal/@_PostDoc/@R scripts/easy_progress_bar.R')
# pb <- easy_progress_bar(359)
# pev_de_trim <- terra::trim(pev_de[[1]])
# for (i in 2:nlyr(tp_de)) {
#     pb$tick()
#     aux <- terra::trim(pev_de[[i]])

#     pev_de_trim <- c(pev_de_trim, aux)
# }

# writeCDF(pev_de_trim, 'files/pev_1970_1999_trim.nc', overwrite = TRUE, compression = 5) 

# tp_de_trim <- terra::rast('files/tp_1970_1999_trim.nc')
# pev_de_trim <- terra::rast('files/pev_1970_1999_trim.nc')

# plot(pev_de_trim[[1]])

# Cell-wise function: vector in -> c(a,b,c)

#' FUNCTIONS TO GET PARAMETERS FROM HISTORICAL DATA -----------------
#' To compute the parameters of the log-logistic distribution, we can use the method of probability-weighted
#'  moments (PWM), you should take the reference period data in raster format, and then apply the function to 
#' each pixel. The output will be a raster with 3 bands (a, b, c) for each month (36 bands in total).
#' This raster of parameters will be used to compute the SPEI for future periods, using the same function
#'  but with the future data and the parameters from the reference period.

#' function to compute log-logistic parameters from a vector 
pwm_to_abc <- function(x, nyears = NYEARS) {

    x <- x[is.finite(x)]
    n <- length(x)

    if (n < nyears) {
      return(c(NA_real_, NA_real_, NA_real_))
    }

    x <- sort(x)
    p <- 1 - ((seq_len(n) - 0.35) / n)

    w0 <- mean(x)
    w1 <- mean(x * p)
    w2 <- mean(x * p^2)

    den <- 6 * w1 - w0 - 6 * w2
    if (!is.finite(den) || den == 0) {
      return(c(NA_real_, NA_real_, NA_real_))
    }

    b <- (2 * w1 - w0) / den
    if (!is.finite(b) || abs(b) < 1e-8) {
      return(c(NA_real_, NA_real_, NA_real_))
    }

    g1 <- gamma(1 + 1 / b)
    g2 <- gamma(1 - 1 / b)
    if (!is.finite(g1) || !is.finite(g2) || g1 == 0 || g2 == 0) {
      return(c(NA_real_, NA_real_, NA_real_))
    }

    a <- (w0 - 2 * w1) * b / (g1 * g2)
    c <- w0 - (w0 - 2 * w1) * b

    return(c(a, b, c))
}

#' function to get all parameters for all timesteps (months)
get_all_abc <- function(def, nyears = NYEARS) {

    pwm_to_abc_local <- function(x, nyears) {
        x <- x[is.finite(x)]
        n <- length(x)

        if (n < nyears) {
          return(c(NA_real_, NA_real_, NA_real_))
        }

        x <- sort(x)
        p <- 1 - ((seq_len(n) - 0.35) / n)

        w0 <- mean(x)
        w1 <- mean(x * p)
        w2 <- mean(x * p^2)

        den <- 6 * w1 - w0 - 6 * w2
        if (!is.finite(den) || den == 0) {
          return(c(NA_real_, NA_real_, NA_real_))
        }

        b <- (2 * w1 - w0) / den
        if (!is.finite(b) || abs(b) < 1e-8) {
          return(c(NA_real_, NA_real_, NA_real_))
        }

        g1 <- gamma(1 + 1 / b)
        g2 <- gamma(1 - 1 / b)
        if (!is.finite(g1) || !is.finite(g2) || g1 == 0 || g2 == 0) {
          return(c(NA_real_, NA_real_, NA_real_))
        }

        a <- (w0 - 2 * w1) * b / (g1 * g2)
        c <- w0 - (w0 - 2 * w1) * b

        c(a, b, c)
    }

    mat_x <- matrix(def, ncol = nyears, byrow = FALSE)
    params <- t(apply(mat_x, 1, pwm_to_abc_local, nyears = nyears))

  return(as.vector(t(params)))
}

#---------------------------------------------------------
# Wrapper for SpatRaster inputs
#---------------------------------------------------------
compute_params <- function(pr, pet,
                           cores = parallel::detectCores() - 2,
                           filename = "files/loglogist_abc.tif",
                           overwrite = TRUE) {
  # Check inputs
  stopifnot(
    compareGeom(pr, pet)
  )

  out <- terra::app(
    x = pr - pet,
    fun = get_all_abc,
    nyears = NYEARS,
    cores = cores
  )

  month_names <- sprintf("m%02d", rep(1:12, each = 3))
  param_names <- rep(c("a", "b", "c"), times = 12)
  names(out) <- paste(month_names, param_names, sep = "_")
  out
}


#---------------------------------------------------------
# Test
#---------------------------------------------------------
NYEARS <- 30 # glocal var, number of years in the reference period
params <- compute_params(pr = tp_de_trim, pet = -1*pev_de_trim, cores = 4)
writeCDF(params, 'files/loglogist_abc.nc', overwrite = TRUE, compression = 5)

#%% 3. spei functions: computing SPEI --------

# tp_de_trim <- terra::rast('files/tp_1970_1999_trim.nc')
# pev_de_trim <- terra::rast('files/pev_1970_1999_trim.nc')
# params <- terra::rast('files/loglogist_abc.nc')



#---------------------------------------------------------
# Wrapper for SpatRaster inputs
#---------------------------------------------------------
compute_spei <- function(tp, pev, param,
                     n = 6,
                     cores = parallel::detectCores() - 2,
                     filename = "",
                     overwrite = TRUE) {
             

  stopifnot(
    compareGeom(tp, pev),
    compareGeom(tp, param)
  )

  if (terra::nlyr(param) != 36) {
    stop("param must have 36 layers: 12 months x 3 parameters (a, b, c)")
  }

  n_time <- terra::nlyr(tp)
  if (terra::nlyr(pev) != n_time) {
    stop("tp and pev must have the same number of layers")
  }
  if (n_time %% 12 != 0) {
    stop("tp and pev must have a number of layers divisible by 12")
  }

  nyear <- n_time / 12

  spei_cell <- function(x, n, n_time, nyear) {


    tp_vals <- x[seq_len(n_time)]
    pev_vals <- x[n_time +seq_len(n_time)]
    param_vals <- x[2*n_time + seq_len(36)]

    vdeficit <- unlist(tp_vals - pev_vals)
    deficit_acc <- runner::runner(vdeficit, k = n, f = function(y) sum(y, na.rm = TRUE), na_pad = TRUE)


    deficit_matrix <- matrix(deficit_acc, ncol = nyear, byrow = FALSE)
    ll_parameters <- matrix(param_vals, ncol = 3, byrow = TRUE)

    a_mat <- matrix(unlist(ll_parameters[, 1]), nrow = 12, ncol = nyear)
    b_mat <- matrix(unlist(ll_parameters[, 2]), nrow = 12, ncol = nyear)
    c_mat <- matrix(unlist(ll_parameters[, 3]), nrow = 12, ncol = nyear)

    prob <- (1 + (a_mat / (deficit_matrix - c_mat)) ^ b_mat) ^ -1
    prob[is.nan(prob) | is.infinite(prob)] <- 1e-4
    prob <- pmin(pmax(prob, 1e-4), 1 - 1e-4)

    c(qnorm(prob))
  }

  out <- terra::app(
    x = c(tp, pev, param),
    fun = spei_cell,
    n = n,
    n_time = n_time,
    nyear = nyear,
    cores = cores,
    filename = filename,
    overwrite = overwrite
  )

  names(out) <- names(tp)
  terra::time(out) <- terra::time(tp)
  out
}

#---------------------------------------------------------
# Test
#---------------------------------------------------------
# spei <- compute_spei(tp = tp_de_trim[[1:12]], pev = pev_de_trim[[1:12]], param = params, cores = 4)
 
# plot(spei[[6]])


#%% 4. Example of usage --------

# reference period datasets
ref_pr <- rast('ref_precipitation.nc')
sfcw <- rast('ref_wind_speed.nc')
rsds <- rast('ref_downward_solar_radiation.nc')
rh <- rast('ref_relative_humidity.nc')
temp <- rast('ref_temperature.nc')

# compute evapotranspiration
ref_eto <- compute_et(tas = temp, rsds = rsds, wind = sfcw, hurs = rh)

# compute parameters of log-logistic distribution
ref_params <- compute_params(pr = ref_pr, pet = ref_eto)

# compute SPEI for future datasets
fut_pr <- rast('fut_precipitation.nc')
fut_sfcw <- rast('fut_wind_speed.nc')
fut_rsds <- rast('fut_downward_solar_radiation.nc')
fut_rh <- rast('fut_relative_humidity.nc')
fut_temp <- rast('fut_temperature.nc')

fut_eto <- compute_et(tas = fut_temp, rsds = fut_rsds, wind = fut_sfcw, hurs = fut_rh)

fut_spei <- compute_spei(tp = fut_pr, pev = fut_eto, param = ref_params)


