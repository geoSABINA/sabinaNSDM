#' @name expl.var.global
#' @aliases expl.var.global
#' @title Global Environmental Covariates
#' @description This dataset contains environmental covariates, primarily climatic variables, at a global scale.
#'
#' @format A \code{\link[terra:rast]{SpatRaster}} object with four layers:
#' \describe{
#'   \item{\code{bio1}}{Annual Mean Temperature (°C).}
#'   \item{\code{bio2}}{Mean Diurnal Range (Mean of monthly (max temp - min temp), °C).}
#'   \item{\code{bio4}}{Temperature Seasonality (standard deviation × 100).}
#'   \item{\code{bio12}}{Annual Precipitation (mm).}
#' }
#'
#' @source \doi{10.1111/ecog.07328}
#' @usage data(expl.var.global)
"expl.var.global"
