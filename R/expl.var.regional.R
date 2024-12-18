#' @name expl.var.regional
#' @aliases expl.var.regional
#' @title Regional Environmental Covariates
#' @description This dataset contains environmental covariates, primarily climatic and edaphic variables, at a regional scale.
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
#'
#' @source \doi{10.1111/ecog.07328}
#' @usage data(expl.var.regional)
"expl.var.regional"