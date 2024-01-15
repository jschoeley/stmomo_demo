# StMoMo Lee-Carter ASFR forecast
#
# Jonas Schöley
#
# This script illustrates how to use the StMoMo framework
# to forecast age specific fertility rates based on historical
# information on births and person-years exposure. We also illustrate
# how to implement a Lee-Carter variant with a parametric beta
# specification which can aid with convergence for difficult series

# Init ------------------------------------------------------------

library(StMoMo)
library(dplyr)
library(ggplot2)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  denmark = 'DEN_Births_Exposures.rds'
)
paths$output <- list(
  data = 'forecast.rds'
)

# constants specific to this analysis
cnst <- within(list(), {
  lee_carter_fitting_period = 2000:2019
  jumpoff = 2020
  h = 20
  forecast_horizon = (jumpoff):(jumpoff + h-1)
  nsim = 250
  seed = 1987
  age_start = 12
  age_end = 55
  year_start = 1993
  year_end = 2022
})

# list containers for analysis artifacts
dat <- list()

# Load data -------------------------------------------------------

dnk <- readRDS(paths$input$denmark)
head(dnk)

# Functions -------------------------------------------------------

#' Function to create object to use with StMoMo()
#'
#' @param df data frame with time series of age specific
#'           births and exposures
#' @param name_births name of birth count variable in df
#' @param name_exposure name of exposure variable in df
#' @param years vector of unique years in the data
#' @param age vector of unique ages in the data
#' @param sex "female"
#' @param label Name of data series
#'
#' @return
#' An StMoMo object
Pop2StMoMo <- function (
    df,
    name_births = 'Births', name_exposure = 'Exposure',
    years, sex, label, age = cnst$age_start:cnst$age_end
) {
  nage = length(age)
  nyear = length(years)

  Dxt_vec <- df[[name_births]]
  Dxt <- matrix(Dxt_vec, nrow = nage, ncol = nyear)
  colnames(Dxt) <- years
  rownames(Dxt) <- age

  Ext_vec <- df[[name_exposure]]
  Ext <- matrix(Ext_vec, nrow = nage, ncol = nyear)
  colnames(Ext) <- years
  rownames(Ext) <- age

  stmomodata <- structure(
    list(Dxt = Dxt, Ext = Ext, ages = age, years = years,
         type = 'central', series = sex, label = label),
    class = "StMoMoData"
  )

  return(stmomodata)
}

#' StMoMo Lee-Carter Forecast
#'
#' @param stmomo StMoMo object created with Pop2StMoMo
#' @param h forecast horizont
#' @param betas specification for LC beta_x parameters
#'              ('np' for nonparametric, 'p0' for constant,
#'               'p1' for linear over age, 'p2' for quadratic,
#'               'p3' for cubic)
#' @param nsim number of simulations
#' @param seed simulation seed
#' @param maxit maximum number of iterations
#'
#' @return
#' List with elements "fit" and "sim"
ForecastLeeCarter <- function (stmomo, h = 1, betas = 'np',
                               nsim = 500, seed = 1987, maxit = 1e2) {

  constLC <- function(ax, bx, kt, b0x, gc, wxt, ages) {
    c1 <- mean(kt[1, ], na.rm = TRUE)
    c2 <- sum(bx[, 1], na.rm = TRUE)
    list(ax = ax + c1 * bx, bx = bx / c2, kt = c2 * (kt - c1))
  }

  f1 <- function (x, ages) {x}
  f2 <- function (x, ages) {x^2}
  f3 <- function (x, ages) {x^3}
  paf <- switch (betas,
    np = 'NP',
    p0 = '1',
    p1 = c(1, f1),
    p2 = c(1, f1, f2),
    p3 = c(1, f1, f2, f3)
  )

  LC_def <- StMoMo(link = 'log', staticAgeFun = TRUE,
                   periodAgeFun = paf)
  LC_fit <- fit(LC_def, data = stmomo, iterMax = maxit, trace = TRUE)
  LC_avg <-
    forecast(
      LC_fit, h = h,
      jumpchoice = 'actual', kt.method = 'mrwd'
    )[['rates']]
  LC_sim <-
    simulate(LC_fit, nsim = nsim, h = h,
             jumpchoice = 'actual',
             seed = seed,
             kt.method = 'mrwd'
    )[['rates']]

  return(list(forecast_fit = LC_fit, forecast_avg = LC_avg,
              forecast_sim = LC_sim))
}

# Forecast --------------------------------------------------------

# create stmomo object from danish birth data
dnk_stmomo <- Pop2StMoMo(
  dnk |> filter(Year %in% cnst$lee_carter_fitting_period),
  name_births = 'Births',
  name_exposure = 'Exposure',
  age = cnst$age_start:cnst$age_end,
  years = cnst$lee_carter_fitting_period,
  sex = 'female',
  label = 'DK'
)
dnk_stmomo

# fit Lee-carter model,
# usually we would set betas = 'np' for the standard Lee-Carter
# but in cases where the model does not converge, as for Denmark,
# we restrict the beta's to follow a cubic function over age
fit <- ForecastLeeCarter(
  dnk_stmomo,
  h = cnst$h, betas = 'p3',
  nsim = cnst$nsim,
  seed = cnst$seed
)

# plot estimated LC terms
plot(fit$forecast_fit, type = 's')
# plot residual map
plot(residuals(fit$forecast_fit), type = 'colourmap')

# Plot TFR forecast -----------------------------------------------

# calculate observed tfr
dnk_tfr <-
  dnk |>
  group_by(Year) |>
  summarise(tfr = sum(Births/Exposure))

# calculate forecasted tfr over simulations
as.data.frame.table(apply(fit$forecast_sim, 2:3, function (asfr) sum(asfr))) |>
  as_tibble() |>
  mutate(Var1 = as.integer(as.character(Var1))) |>
  ggplot() +
  geom_line(aes(x = Var1, y = Freq, group = Var2), alpha = 0.1) +
  geom_point(data = dnk_tfr, aes(x = as.integer(Year), y = tfr)) +
  coord_cartesian(ylim = c(1, 3)) +
  labs(y = 'TFR', x = 'Year', title = 'Lee-carter TFR forecast for Denmark starting 2020') +
  theme_minimal()
