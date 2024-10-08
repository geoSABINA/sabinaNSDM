#' @name NSDM.Regional
#'
#' @title Perform regional scale species distribution modeling for spatially-nested hierarchical species distribution modeling (NSDM) analysis.
#'
#' @description This function calibrates, evaluates, and projects species distribution models at the \bold{regional} scale for \bold{NSDM} analysis. It generates ensemble models (combining a range of single-algorithm models), evaluates models performance, and projects models to current and new environmental conditions.
#'
#'
#' @param nsdm_selvars An object of class \code{nsdm.vinput} containing selected covariates for NSDM generated using the \code{\link{NSDM.SelectCovariates}} function.
#' @param algorithms (\emph{optional, default} \code{'c("GLM", "GAM", "RF")'}) \cr
#' A \code{vector} containing the statistical algorithms to use for modeling. Options are \code{'GLM'}, \code{'GAM'}, \code{'GBM'}, \code{'MAXNET'}, \code{'MARS'}, and/or \code{'RF'}.
#' @param CV.nb.rep (\emph{optional, default} \code{10}) \cr
#' An \code{integer} corresponding to the number of cross-validation sets (repetitions).
#' @param CV.perc (\emph{optional, default} \code{0.8}) \cr
#' A \code{numeric} between \code{0} and \code{1} defining the percentage of data that will be kept for calibration in each cross-validation set.
#' @param CustomModelOptions (\emph{optional, default} \code{NULL}) \cr
#' A \code{\link{BIOMOD.models.options}} object returned by the \code{\link{bm_ModelingOptions}} to tune models options. If \code{NULL} (the default), biomod2's default parameters are used.
#' @param metric.select.thresh (\emph{optional, default} \code{0.8}) \cr
#' A \code{numeric} between \code{0} and \code{1} corresponding to the minimum scores of AUC below which single models will be excluded from the ensemble model building.
#' @param save.output (\emph{optional, default} \code{TRUE}) \cr
#' A \code{logical} value defining whether the outputs should be saved at local.
#' @param rm.biomod.folder (\emph{optional, default} \code{TRUE}) \cr
#' A \code{logical} value indicating whether the intermediate biomod2's folders should be removed after processing.
#'
#'
#' @return An object of class \code{nsdm.predict.g} containing model information, predictions and evaluation statistics:
#' - `$SpeciesName` Name of the species.
#' - `$SpeciesData.XY.Global` Species occurrences at the global level at \code{data.frame} format after applying spatial thinning.
#' - `$SpeciesData.XY.Regional` Species occurrences at the regional level at \code{data.frame} format after applying spatial thinning.
#' - `$Background.XY.Global` Background data at the global level at \code{data.frame} format.
#' - `$Background.XY.Regional` Background data at the regional level at \code{data.frame} format.
#' - `$Absences.XY.Global` Absence data at the global level at \code{data.frame} format.
#' - `$Absences.XY.Regional` Absence data at the regional level at \code{data.frame} format.
#' - `$Scenarios` A \code{list} containing future scenarios in \code{\link[terra:rast]{PackedSpatRaster}} format.
#' - `$Selected.Variables.Global` A \code{character} vector specifying the names of the selected covariates at the global scale.
#' - `$IndVar.Global.Selected` Selected covariates at the global level in \code{\link[terra:rast]{PackedSpatRaster}} format.
#' - `$Selected.Variables.Regional` A \code{character} vector specifying the names of the selected covariates at the regional scale.
#' - `$IndVar.Regional.Selected` Selected covariates at the regional level in \code{\link[terra:rast]{PackedSpatRaster}} format.
#' - `$IndVar.Global.Selected.reg` Selected covariates at the global level for regional projections in \code{\link[terra:rast]{PackedSpatRaster}} format.
#' - `$args` A \code{list} containing the arguments used during modelling, including: `algorithms`, `CV.nb.rep`, `CV.perc` and `metric.select.thresh`.
#' - `$nbestreplicates` A \code{data.frame} containing  the number of replicates meeting or exceeding the specified \code{metric.select.thresh} for each algorithm used in the modeling.
#' - `$current.projections` A \code{list} containing: \code{Pred}, a \code{\link[terra:rast]{PackedSpatRaster}} representing the current continuous (suitability) projection; \code{Pred.bin.ROC} a \code{\link[terra:rast]{PackedSpatRaster}} representing binary projections generated through the optimization of the AUC statistic as a threshold; and \code{Pred.bin.TSS} a \code{\link[terra:rast]{PackedSpatRaster}} representing binary projections generated through the optimization of the TSS statistic as a threshold.
#' - `$myEMeval.replicates` Evaluation statistics for each replicate model according to different evaluation metrics (ROC, TSS, KAPPA, ACCURACY, SR, and BOYCE).
#' - `$myEMeval.Ensemble` Evaluation statistics for the ensemble model according to different evaluation metrics (ROC, TSS, KAPPA).
#' - `$myModelsVarImport` Variable importance measures for individual models.
#' - `$new.projections` A \code{list} containing: \code{Pred.Scenario}, the continuous (suitability) projections onto new scenarios in a \code{\link[terra:rast]{PackedSpatRaster}} format; \code{Pred.bin.ROC.Scenario} a \code{\link[terra:rast]{PackedSpatRaster}} representing binary projections onto new scenarios generated through the optimization of the AUC statistic as a threshold; and \code{Pred.bin.TSS.Scenario} a \code{\link[terra:rast]{PackedSpatRaster}} representing binary projections onto new scenarios generated through the optimization of the TSS statistic as a threshold.
#' - `Summary` Summary information about the modeling process.
#'
#'
#' @details
#' This function uses the (\emph{biomod2} package to generate, evaluate, and project species distribution models at the \bold{regional} scale for \bold{NSDM} analysis.
#' If `save.output=TRUE`, modelling results are stored out of R in the \emph{Results/} folder created in the current working directory:
#' - the \emph{Results/Regional/Projections/} folder, containing the continuous and binary current and new projections. Current projections are named with the species name followed by \file{.Current.tif}, \file{.bin.ROC.tif} and \file{.bin.TSS.tif}. New projections are named with the species name followed by the scenario name, and \file{.bin.ROC.tif}, \file{.bin.TSS.tif} when binary.
#' - the \emph{Results/Regional/Values/} folder, containing replicates statistics, the consensus model statistics, the covariate importance, and the \code{nbestreplicates}, named with the species name and \file{.__replica.csv}, \file{._ensemble.csv}, \file{._indvar.csv} and \file{._nbestreplicates.csv} respectively.
#'
#'
#' @seealso \code{\link{NSDM.InputData}}, \code{\link{NSDM.FormattingData}}, \code{\link{NSDM.SelectCovariates}}
#'
#'
#' @examples
#' library(sabinaNSDM)
#'
#' # Load species occurrences
#' data(Fagus.sylvatica.xy.global, package = "sabinaNSDM")
#' data(Fagus.sylvatica.xy.regional, package = "sabinaNSDM")
#'
#' # Load covariates
#' data(expl.var.global, package = "sabinaNSDM")
#' data(expl.var.regional, package = "sabinaNSDM")
#' expl.var.global<-terra::unwrap(expl.var.global)
#' expl.var.regional<-terra::unwrap(expl.var.regional)
#'
#' # Load new scenarios
#' data(new.env, package = "sabinaNSDM")
#' new.env<-terra::unwrap(new.env)
#'
#' # Prepare input data
#' myInputData<-NSDM.InputData(SpeciesName = "Fagus.sylvatica",
#'				spp.data.global = Fagus.sylvatica.xy.global,
#'				spp.data.regional = Fagus.sylvatica.xy.regional,
#'				expl.var.global = expl.var.global,
#'				expl.var.regional = expl.var.regional,
#'				new.env = new.env,
#'				new.env.names = c("Scenario1"),
#'				Background.Global = NULL,
#'				Background.Regional = NULL,
#'				Absences.Global = NULL,
#'				Absences.Regional = NULL)
#'
#' # Format the input data
#' myFormattedData <- NSDM.FormattingData(myInputData,
#'					nPoints=1000)
#' # Select covariates
#' mySelectedCovs <- NSDM.SelectCovariates(myFormattedData)
#'
#' # Perform regional scale SDMs with default parameters.
#' myRegionalModel <- NSDM.Regional(mySelectedCovs)
#'
#' summary(myRegionalModel)
#'
#' 
#' ## Perform regional scale SDMs with custom parameters.
#' ## This line shows an example how to customize modeling options using `bm_ModelingOptions`
#' ## `bm_ModelingOptions` from the `biomod2` package 
#' # opt.b <- bm_ModelingOptions(data.type = 'binary', 
#' #  				models = c("GAM", "GBM", "RF"),
#' #				strategy = 'bigboss')
#'   
#' # myRegionalModel <- NSDM.Regional(
#' #				     # Selected covariates output used as input
#' #				     mySelectedCovs,
#' #				     # Statistical models used in the ensemble
#' #				     algorithms = c("GBM", "RF", "GLM"),
#' #				     # Number of cross-validation replicates
#' #				     CV.nb.rep = 10,
#' #				     # Percentage of data used for each replicate
#' #				     CV.perc = 0.8,
#' #				     # Optional custom options for statistical models
#' #				     CustomModelOptions = opt.b,
#' #				     # Threshold for selecting models for ensemble
#' #				     metric.select.thresh = 0.8,
#' #				     #  Save the output externally
#' #				     save.output = TRUE,
#' #				     # Remove the temporary biomod2 output folder
#' #				     rm.biomod.folder = TRUE)
#'
#'
#' @export
NSDM.Regional <- function(nsdm_selvars,
                          algorithms=c( "GLM", "GAM", "RF"),
                          CV.nb.rep=10,
                          CV.perc=0.8,
                          CustomModelOptions=NULL,
                          metric.select.thresh = 0.8,
                          save.output=TRUE,
                          rm.biomod.folder=TRUE){
  if(!inherits(nsdm_selvars, "nsdm.vinput")){
          stop("nsdm_selvars must be an object of nsdm.vinput class.",
               " Consider running NSDM.SelectCovariates() function")
  }

  mod_call <- as.list(match.call()[-1])
  mod_call$model.type = "Regional"
  names(mod_call) <- gsub("nsdm_selvars.*", "nsdm.obj",
                          names(mod_call))
  sabina <- do.call(general_nsdm_model, mod_call)
  attr(sabina, "class") <- "nsdm.predict.r"
  class(sabina) <- c("nsdm.predict.r", "nsdm.predict")

  return(sabina)

}
