#' @name new.env
#' @aliases new.env
#' @title Future Environmental Scenarios (2071–2100)
#' @description This dataset contains environmental covariates predicted under the Global Climate Model IPSL_CM6A_LR and the SSP126 socio-economic and emissions pathway for the 2071–2100 period.
#'
#' @format A \code{\link[terra:rast]{SpatRaster}} object with seven layers:
#' \describe{
#'   \item{\code{bio1}}{Annual Mean Temperature (°C).}
#'   \item{\code{bio2}}{Mean Diurnal Range (Mean of monthly (max temp - min temp), °C).}
#'   \item{\code{bio3}}{Isothermality (BIO2/BIO7) × 100.}
#'   \item{\code{bio4}}{Temperature Seasonality (standard deviation × 100).}
#'   \item{\code{bio12}}{Annual Precipitation (mm).}
#'   \item{\code{sand}}{Sand content (%).}
#'   \item{\code{radiation}}{Annual Solar Radiation (kJ/m²).}
#' }
#'
#' @source \doi{10.1111/ecog.07328}
#' @usage data(new.env)
"new.env"


