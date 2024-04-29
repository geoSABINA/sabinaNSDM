#' @name NSDM.Multiply
#'
#' @title Perform spatially-nested hierarchical species distribution modeling (NSDM) analysis with the multiply strategy
#'
#' @description This function generates and projects a \bold{NSDM} with the \bold{multiply} strategy. It averages the prediction maps of regional scale and global scale species distribution models.
#'
#'
#' @param nsdm_global An object of class \code{nsdm.predict.g} containing a global model generated using the \code{\link{NSDM.Global}} function.
#' @param nsdm_regional An object of class \code{nsdm.predict.r} containing a regional model generated using the \code{\link{NSDM.Regional}} function.
#' @param method (\emph{optional, default} \code{'Arithmetic'}) \cr
#'  A \code{character} specifying the method for averaging the global and regional predictions. Options are \code{'Arithmetic'}, or \code{'Geometric'}. The \code{'Arithmetic'} method (the default) calculates the arithmetic mean of global and regional predictions, and \code{'Geometric'} calculates a geometric mean.
#' @param rescale (\emph{optional, default} \code{FALSE}) \cr
#' A \code{logical} controlling whether global and regional model predictions should be rescaled before averaging them.
#' @param save.output (\emph{optional, default} \code{TRUE}) \cr
#' A \code{logical} value defining whether the outputs should be saved at local.
#'
#'
#' @return An object of class \code{nsdm.predict} containing model information, predictions and evaluation statistics: 
#' - `$SpeciesName` Name of the species.
#' - `$args` A \code{list} containing the arguments used during modelling, including: `method`, and `rescale`.
#' - `$current.projections` A \code{list} containing: \code{Pred}, a \code{\link[terra:rast]{PackedSpatRaster}} representing the current projection.
#' - `$new.projections` A \code{list} containing: \code{Pred.Scenario}, the projections onto new scenarios in a \code{\link[terra:rast]{PackedSpatRaster}} format.
#'
#'
#' @details
#' This function generates a \bold{NSDM} with the \bold{multiply} strategy. It averages the predictions of species distribution models at regional and global scales.
#' If `save.output=TRUE`, modelling results are stored out of R in the \emph{Results/} folder created in the current working directory:
#' - the \emph{Results/Multiply/Projections/} folder, containing the continuous current and new projections. Current projections are named with the species name followed by \file{.Current.tif}.
#'
#'
#' @examples
#' myMultiplyModel <- NSDM.Multiply(myGlobalModel,
#'
#'
#' @seealso \code{\link{NSDM.InputData}}, \code{\link{NSDM.FormattingData}}, \code{\link{NSDM.SelectCovariates}}, \code{\link{NSDM.Global}}, \code{\link{NSDM.Regional}}
#'
#'
#' @export
NSDM.Multiply <- function(nsdm_global,
                          nsdm_regional,
                          method="Arithmetic",
                          rescale=FALSE,
                          save.output=TRUE) {

  if(!inherits(nsdm_regional, "nsdm.predict.r") || !inherits(nsdm_global, "nsdm.predict.g")) {
    stop("nsdm_regional and nsdm_global must be objects of class nsdm.predict.r and nsdm.predict.r, respectively")
  }

  if(!method %in% c("Arithmetic", "Geometric")) {
    stop("Please select 'Arithmetic' or 'Geometric' method")
  }

  if(!identical(nsdm_regional$Species.Name, nsdm_global$Species.Name)) {
    stop("Different species in global and regional levels")
  } else {
    SpeciesName <- nsdm_regional$Species.Name
  }

  sabina<-nsdm_global[names(nsdm_global) %in% c("Species.Name")]
  sabina<-list()
  sabina$args <- list()
  sabina$args$method <- method
  sabina$args$rescale <- rescale

  Scenarios <- names(nsdm_regional$Scenarios)
  Scenarios <- c("Current",Scenarios)

  for(i in seq_along(Scenarios)) {
    projmodel <- Scenarios[i]
    # Load raster data at global and regional scales
    if(projmodel =="Current") {
      Pred.global <- terra::unwrap(nsdm_global$current.projections$Pred) # Unwrap object
      Pred.regional <- terra::unwrap(nsdm_regional$current.projections$Pred)  # Unwrap object
    } else {
        Pred.global <- terra::unwrap(nsdm_global$new.projections$Pred.Scenario[[i-1]])  # Unwrap object
        Pred.regional <- terra::unwrap(nsdm_regional$new.projections$Pred.Scenario[[i-1]])  # Unwrap object
      }

    # Rescale global prediction to a common range
    if(rescale==TRUE) {
      # Calculate minimum and maximum values of global prediction
      min_val <- min(terra::values(Pred.global), na.rm = TRUE)
      max_val <- max(values(Pred.global), na.rm = TRUE)

      Pred.global <- terra::app(Pred.global, fun = function(x) {
      ((x - min_val) / (max_val - min_val) * 999) + 1})
      # Calculate minimum and maximum values of regional prediction
      min_val <- min(values(Pred.regional), na.rm = TRUE)
      max_val <- max(values(Pred.regional), na.rm = TRUE)

      Pred.regional <- terra::app(Pred.regional, fun = function(x) {
      ((x - min_val) / (max_val - min_val) * 999) + 1})
    }

    # Average global and regional
    Stack.rasters <- c(Pred.global, Pred.regional)
    if(method=="Geometric") {
      res.average <-  sqrt(Pred.global*Pred.regional)
    } else if(method=="Arithmetic")  {
      res.average <-  terra::mean(Stack.rasters)
    }

    res.average<-terra::rast(wrap(res.average))
    if(projmodel =="Current") {
      sabina$current.projections$Pred <- setNames(res.average, paste0(SpeciesName, ".Current"))
    } else {
      sabina$new.projections$Pred.Scenario[[i-1]]<- setNames(res.average, paste0(SpeciesName,".", projmodel))
    }

    if(save.output){
      fs::dir_create(c("Results/Multiply/Projections/"))
      file_path<-paste0("Results/Multiply/Projections/",SpeciesName,".",projmodel,".tif")
      terra::writeRaster(res.average, file_path, overwrite = TRUE)
    }
  } # end for


  ## Evaluation multiply model 	
  myResp.xy <- rbind(nsdm_regional$SpeciesData.XY.Regional, nsdm_regional$Background.XY.Regional)
  myResp <- data.frame(c(rep(1,nrow(nsdm_regional$SpeciesData.XY.Regional)),rep(0,nrow(nsdm_regional$Background.XY.Regional))))
  pred_df <- sabina$current.projections$Pred
  pred_df <- terra::extract(pred_df, myResp.xy)[-1]
  myExpl <- pred_df
  myResp <- as.vector(myResp)[[1]]
  pred_df <- as.vector(pred_df)[[1]]


  myRespNA <- replace(myResp, myResp == 0, NA)
  myBiomodData <- biomod2::BIOMOD_FormatingData(resp.var = myRespNA,
	                                     resp.xy = myResp.xy,
	                                     expl.var = myExpl,
	                                     resp.name = SpeciesName,
	                                     PA.nb.rep = 1,
	                                     PA.nb.absences = nrow(nsdm_regional$Background.XY.Regional),
	                                     PA.strategy = "random")

  calib.lines <- bm_CrossValidation(bm.format = myBiomodData,
                                    strategy = "random",
                                    nb.rep = nsdm_regional$args$CV.nb.rep,
                                    perc =  nsdm_regional$args$CV.perc,
                                    k = NULL,
                                    balance = NULL,
                                    env.var = NULL,
                                    strat = NULL,
                                    user.table = NULL,
                                    do.full.models = FALSE)

  metric.eval<- c("ROC","TSS","KAPPA")
  stat.validation <- data.frame()
  cross.validation <- data.frame()

  for(i in seq_len(ncol(calib.lines))) {
    calib.lines.rep <- calib.lines[, i, drop = FALSE]
    eval.lines.rep <- which(rowSums(!calib.lines.rep) == ncol(calib.lines.rep))

    for(xx in metric.eval) {
      stat <-bm_FindOptimStat(metric.eval = xx,
  			obs = myResp[-eval.lines.rep],
  			fit = pred_df[-eval.lines.rep])
      cross.validation <- rbind(cross.validation, stat)
    }

    for(xx in metric.eval) {
      stat <-bm_FindOptimStat(metric.eval = xx,
  			obs = myResp[eval.lines.rep],
  			fit = pred_df[eval.lines.rep],
			threshold = cross.validation["cutoff", xx])
      stat.validation <- rbind(stat.validation, stat)
    }
  } # end for i calib.lines

  colnames(cross.validation)[which(colnames(cross.validation) == "best.stat")] <- "calibration"
  cross.validation$validation <- stat.validation$best.stat

  metric.means <- aggregate(validation ~ metric.eval, data = cross.validation, FUN = mean)

  # Summary
  summary <- data.frame(Values = c(SpeciesName,
  				paste(nsdm_regional$args$algorithms,collapse = ", "),
  				sum(nsdm_regional$myEMeval.replicates$metric.eval == "ROC" & nsdm_regional$myEMeval.replicates$validation >= nsdm_regional$args$CV.perc),
  				nsdm_regional$myEMeval.Ensemble$calibration[which(nsdm_regional$myEMeval.Ensemble$metric.eval=="ROC")],
  				nsdm_regional$myEMeval.Ensemble$calibration[which(nsdm_regional$myEMeval.Ensemble$metric.eval=="TSS")],
  				nsdm_regional$myEMeval.Ensemble$calibration[which(nsdm_regional$myEMeval.Ensemble$metric.eval=="KAPPA")],
  				method,
  				round(metric.means$validation[metric.means$metric.eval == "ROC"], 2),
  				round(metric.means$validation[metric.means$metric.eval == "TSS"], 2),
  				round(metric.means$validation[metric.means$metric.eval == "KAPPA"], 2)))

  rownames(summary) <- c("Species name",
  				"Statistical algorithms at regional level",
  				paste0("Number of replicates with AUC > ",nsdm_regional$args$CV.perc," at regional level"),
  				"AUC of ensemble model at regional level",
  				"TSS of ensemble model at regional level",
  				"KAPPA of ensemble model at regional level",
  				"Multiply method",
  				"AUC of hierarchical multiply ensemble model",
  				"TSS of hierarchical multiply ensemble model",
 				"KAPPA of hierarchical multiply ensemble model")

  sabina$Summary <- summary

  # Wrap objects
  sabina$current.projections <- rapply(sabina$current.projections, terra::wrap, how = "list")
  if(!is.null(nsdm_global$Scenarios)) {
    sabina$new.projections <- rapply(sabina$new.projections, terra::wrap, how = "list")
  }

  attr(sabina, "class") <- "nsdm.predict"

  # save.out messages
  if(save.output){
  message("Results saved in the following local folder/s:")
  message(paste(
    "- Hierarchical Multiply Models: /Results/Multiply/Projections/\n"
  ))
  }

  return(sabina)

}

