"
Compute Potential Evapotranspiration (PET) using the Penman-Monteith equation and SPEI

For historical and future climate conditions

Pedro Alencar
27.07.26
"


#%% 0. Load libraries --------

# data manipulation
library(tidyr)
library(dplyr)
library(terra)
library(lubridate)
library(stringr)

# visualization
library(ggplot2)

# spei
library(SPEI)

# get httpgd running for plots
httpgd::hgd()
httpgd::hgd_browse()

#%%1. Read data --------


test_hyras <- rast('data/raw_data/TN_min/TN_min_hyras_5_1951-2020_v5-0_19510101-20201231_year-timeseries.nc')
test_hyras

test_proj <- rast('data/data_preprocessed_normalized/01b_norm_and_log/2006-2025/rcp26/RSDS.tif')

plot(test_hyras[[1]])
plot(test_proj[[1]])

# in test_proj get pixel in north that is non NA
north_pixel <- head(which(!is.na(values(test_hyras[[1]])), arr.ind = TRUE), 1)
south_pixel <- tail(which(!is.na(values(test_hyras[[1]])), arr.ind = TRUE), 1)

xyFromCell(test_hyras, 3920)
xyFromCell(test_hyras, 43543)

#%% 2. Functions evapotranspiration --------
atm_press <- function(alt){
    #' return pressire in kPa given altitude in meters

    if (is.na(alt)) {
        return(101.3)
    } 
    pressure <- 101.3 * ((293 - 0.0065 * alt) / 293) ^ 5.26
    return(pressure)
}

psychometric_constant <- function(alt) {
    
   #' 0.000665 = c_p / (lambda * epsilon) where:
   #'  c_p = specific heat of air at constant pressure (0.001013 MJ kg-1 °C-1), 
   #' lambda = latent heat of vaporization (2.45 MJ kg-1), 
   #' epsilon = ratio of molecular weight of water vapor to dry air (0.622)

    gamma <- 0.000665 * atm_press(alt)
    return(gamma)
}

sat_vapour_pressure <- function(tmean) {
    #' return saturation vapour pressure in kPa given mean temperature in °C
    es <- 0.6108 * exp((17.27 * tmean) / (tmean + 237.3))
    return(es)
}

slope_vapour_pressure_curve <- function(tmean) {
    #' return slope of vapour pressure curve in kPa °C-1 given mean temperature in °C
    delta <- (4098 * sat_vapour_pressure(tmean)) / ((tmean + 237.3)^2)
    return(delta)
}

get_eto <- function(dates, tmean, hrus, sfcw, rsds, alt = NA, resolution = "daily") {

    #' tmean in degrees Celsius, 
    #' hrus in absolute (0-1)
    #' sfcw in m/s
    #' rsds in MJ m-2 day-1
    #' alt in meters
    #' 
    if (is.na(alt)) {
        gamma = 0.000665 * 101.3
    } else {
        gamma = psychometric_constant(alt)
    }

    delta = slope_vapour_pressure_curve(tmean)
    e_s = sat_vapour_pressure(tmean)

    df <- data.frame(dates,eto = NA)

    df$eto = 0.408 * delta * rsds + gamma * (900 / (tmean + 273)) * sfcw * e_s * (1 - hrus) / (delta + gamma * (1 + 0.34 * sfcw))

    if (resolution == "monthly") {
        # count number of days per month
        df$eto <- df$eto * lubridate::days_in_month(df$dates)
    }

    return(df)
}

#' ETo = $\frac{0.408 \Delta (R_n - G) + \gamma \frac{900}{T+273}u_2(e_s - e_a)}{\Delta + \gamma (1+0.34u_2)}$
#' ETo reference evapotranspiration [mm day-1],
#' R_n net radiation at the crop surface [MJ m-2 day-1],
#' G soil heat flux density [MJ m-2 day-1],
#' T air temperature at 2 m height [°C],
#' u_2 wind speed at 2 m height [m s-1],
#' e_s saturation vapour pressure [kPa],
#' e_a actual vapour pressure [kPa],
#' e_s - e_a saturation vapour pressure deficit [kPa],
#' \Delta slope vapour pressure curve [kPa °C-1],
#' \gamma psychrometric constant [kPa °C-1].

#%% 3. Functions SPEI --------

#' Calibrate parameters of multiple log-logist functions using Probability
#' Weighted Moments (PWM) according to Singh (1998) "Entropy-based parameter
#' estimation in Hydrology" (Ch. 18)
#'
#' @description
#' Auxiliar function or calculate SPEI and EDDI.
#'
#' @param var A matrix or data frame with the variable (hydrologic deficit for
#' SPEI, potential evapotranspitation for EDDI). The data has to be organized with
#' weeks or months in rows (i.e. 52 rows of weeks, 12 rows of months) and years in columns.
#' @param n_param Number of parameters of the Log-logist distribution (either 2
#' or 3, default = 3; Vicente-Serrano 2010).
param_loglogist <- function(var, n_param = 3){
  number_years <- ncol(var)
  param <- 1:number_years
  param <- 1- (param-0.35)/number_years
  param <- matrix(rep(param,nrow(var)), ncol = number_years,
                  nrow = nrow(var), byrow = T)

  w0 <- rowSums(var, na.rm = T)/number_years
  w1 <- rowSums(var*param,na.rm = T)/number_years
  w2 <- rowSums(var*param^2, na.rm = T)/number_years

  b <- (2*w1 - w0)/(6*w1 - w0 - 6*w2)
  a <- (w0 - 2*w1)*b/(gamma(1 + 1/b) * gamma(1 - 1/b))
  c <- w0 - (w0 - 2*w1)*b

  if (n_param == 2){
    parameters <- data.frame(a = a,b = b)
  } else if (n_param == 3){
    parameters <- data.frame(a = a,b = b, c = c)
  } else {
    parameters <- NA
  }

  return(parameters)
}

#' Standard Precipitation Evaporation Index calculation
#'
#' @param vdeficit a data.frame column or vector with daily hydrological deficit
#' obtained by the difference of precipitation and potential evapotranspitation (P - ET0)
#' @param nyear a natural number that indicates the number of years in the time series.
#'  It is used to organize the data in a matrix with years in columns and months or weeks in rows.
#' @param ll_parameters a data.frame with the parameters of the log-logistic distribution (a, b, c)
#' obtaied from the param_loglogist function. The number of rows in the data frame has to be equal
#' to the number of periods in a year (e.g. 12 months) 
#' @param n a natural number that indicates the accumulation time (pentad, week, month, etc)
#'
#' @return The function return a list with two elements. One data frame with time stamped pentad values and a matrix with years organized in columns.
#'
#' @description Internal function to calculate the SPEI
get_spei <- function(vdeficit, nyear = 20, ll_parameters, n=6){

#   colnames(ll_parameters) <- c("a", "b", "c") # correct names of the parameters data frame
  
  #get number of years in the time series
  # nyear <- max(lubridate::year(vtime)) - min(lubridate::year(vtime)) + 1

  #get accumulated deficit over n periods
  deficit_acc <- runner::runner(vdeficit, f = function(x) sum(x), k = n)
  deficit_acc[1:(n-1)] <- -1e7 #necessary to sort

  #get data as matrix with nyear columns.
  #The number of rows is equal to the number of periods in a year (12 months, 52 weeks, etc)
  deficit_matrix <- as.data.frame(matrix(deficit_acc, ncol = nyear, byrow = F))

  #sort the values to assess parameters
  deficit_matrix_sort <- t(apply(deficit_matrix, 1, sort))
  deficit_matrix_sort[deficit_matrix_sort == -1e7] <- NA

  #get probabilities to assess SPEI
  deficit_matrix[deficit_matrix == -1e7] <- NA
  prob <- (1 + (ll_parameters$a/(deficit_matrix - ll_parameters$c))^ll_parameters$b)^-1
  aux <- unlist(prob)
  aux[is.nan(aux)] <- 1e-4
  prob <-matrix(aux, dim(prob))

  spei_ts <- c(qnorm(prob))

  return(spei_ts)
}

#' how to run spei functions:
#' 
#' 1) for the reference period (1970-2000) get the parameters of the log-logistic distribution
#'  using param_loglogist function. It will produce one dataframe per grid-point. 
#' 2) for the periods of interest, run using the computed monthly deficit and the 
#' calibrated parameters from step 1 to get the SPEI values.
