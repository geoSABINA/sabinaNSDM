#' Data new.env
#'
#' A dataset containing multiple environmental covariates predicted in new scenarios for projecting models to different spatial or temporal scales. Specifically, these covariates correspond to the 2071–2100 period under the Global Climatic Model IPSL_CM6A_LR and the socio-economic and emissions pathway SSP126.  
#'
#' @format ## `new.env`
#' A \code{SpatRaster} object with the following layers:
#' \describe{
#'   \item{bio1}{Annual Mean Temperature}
#'   \item{bio2}{Mean Diurnal Range (Mean of monthly (max temp - min temp))}
#'   \item{bio3}{Isothermality (BIO2/BIO7) (×100)}
#'   \item{bio4}{Temperature Seasonality (standard deviation ×100)}
#'   \item{bio12}{Annual Precipitation}
#'   \item{sand}{Sand content}
#'   \item{radiation}{Annual solar radiation}
#' }
#'
#'
#' @source <https://nsojournals.onlinelibrary.wiley.com/doi/full/10.1111/ecog.07328>
"new.env"
