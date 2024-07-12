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
#' - `$current.projections` A \code{list} containing: \code{Pred}, a \code{\link[terra:rast]{PackedSpatRaster}} representing the current continuous (suitability) projection; \code{Pred.bin.ROC} a \code{\link[terra:rast]{PackedSpatRaster}} representing binary projections generated through the optimization of the AUC statistic as a threshold; and \code{Pred.bin.TSS} a \code{\link[terra:rast]{PackedSpatRaster}} representing binary projections generated through the optimization of the TSS statistic as a threshold.
#' - `$new.projections` A \code{list} containing: \code{Pred.Scenario}, the continuous (suitability) projections onto new scenarios in a \code{\link[terra:rast]{PackedSpatRaster}} format; \code{Pred.bin.ROC.Scenario} a \code{\link[terra:rast]{PackedSpatRaster}} representing binary projections onto new scenarios generated through the optimization of the AUC statistic as a threshold; and \code{Pred.bin.TSS.Scenario} a \code{\link[terra:rast]{PackedSpatRaster}} representing binary projections onto new scenarios generated through the optimization of the TSS statistic as a threshold.
#' - `$myEMeval.means` Evaluation statistics according to different evaluation metrics (ROC, TSS, KAPPA, ACCURACY, SR, and BOYCE). Please note that this is neither an independent evaluation nor a cross-validation, as the final model is the average of the global and regional models. Original occurrences and background data were used to validate the generated multiply model and to calculate optimal threshold for binary models. 
#'
#'
#' @details
#' This function generates a \bold{NSDM} with the \bold{multiply} strategy. It averages the predictions of species distribution models at regional and global scales.
#' If `save.output=TRUE`, modelling results are stored out of R in the \emph{Results/} folder created in the current working directory:
#' - the \emph{Results/Multiply/Projections/} folder, containing the continuous and binary current and new projections. Current projections are named with the species name followed by \file{.Current.tif}.
#' - the \emph{Results/Multiply/Values/} folder, containing statistics, named with the species name and \file{.__replica.csv}, and \file{._ensemble.csv},  respectively.
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
#'
#' # Select covariates
#' mySelectedCovs <- NSDM.SelectCovariates(myFormattedData)
#'
#' # Perform global scale SDMs
#' myGlobalModel <- NSDM.Global(mySelectedCovs)
#'
#' # Perform regional scale SDMs
#' myRegionalModel <- NSDM.Regional(mySelectedCovs)
#' 
#' # Perform NSDM analysis using the multiply strategy
#' myMultiplyModel <- NSDM.Multiply(myGlobalModel, 	# Output from the global scale SDM
#'				    myRegionalModel, 	# Output from the regional scale SDM
#'				    method = "Arithmetic", # Method for combining model outputs. Options include "Arithmetic" or "Geometric"
#'				    rescale = FALSE, 	# Whether to rescale the model outputs before combining
#'				    save.output = TRUE) # Save the combined model output externally
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
    stop("nsdm_regional and nsdm_global must be objects of class nsdm.predict.r and nsdm.predict.g, respectively.")
  }

  if(!method %in% c("Arithmetic", "Geometric")) {
    stop("Please select 'Arithmetic' or 'Geometric' method.")
  }

  if(!identical(nsdm_regional$Species.Name, nsdm_global$Species.Name)) {
    stop("Different species in global and regional levels.")
  } else {
    SpeciesName <- nsdm_regional$Species.Name
  }

  sabina<-nsdm_global[names(nsdm_global) %in% c("Species.Name")]
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
      max_val <- max(terra::values(Pred.global), na.rm = TRUE)

      Pred.global <- terra::app(Pred.global, fun = function(x) {
      ((x - min_val) / (max_val - min_val) * 999) + 1})
      # Calculate minimum and maximum values of regional prediction
      min_val <- min(terra::values(Pred.regional), na.rm = TRUE)
      max_val <- max(terra::values(Pred.regional), na.rm = TRUE)

      Pred.regional <- terra::app(Pred.regional, fun = function(x) {
      ((x - min_val) / (max_val - min_val) * 999) + 1})
    }

    # Average global and regional
    Stack.rasters <- c(Pred.global, Pred.regional)
    if(method=="Geometric") {
      res.average<-sqrt(Pred.global*Pred.regional)
    } else if(method=="Arithmetic")  {
      res.average <-  terra::mean(Stack.rasters)
    }

    res.average <- terra::rast(terra::wrap(res.average))
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
  }


  ## Evaluation multiply model
  if(!is.null(nsdm_regional$Background.XY.Regional)) {
    # Format the response (presence/background) and covariates data for BIOMOD2
    myResp.xy <- rbind(nsdm_regional$SpeciesData.XY.Regional, nsdm_regional$Background.XY.Regional)
    myResp <- data.frame(c(rep(1,nrow(nsdm_regional$SpeciesData.XY.Regional)),rep(0,nrow(nsdm_regional$Background.XY.Regional))))
  } else {
    # Format the response (presence/absence) and covariate data for BIOMOD2
    myResp.xy <- rbind(nsdm_regional$SpeciesData.XY.Regional, nsdm_regional$Absences.XY.Regional)
    myResp <- data.frame(c(rep(1,nrow(nsdm_regional$SpeciesData.XY.Regional)), rep(0,nrow(nsdm_regional$Absences.XY.Regional))))
  }
  pred_df <- sabina$current.projections$Pred
  pred_df <- terra::extract(pred_df, myResp.xy)[-1]
  myExpl <- pred_df
  myResp <- as.vector(myResp)[[1]]
  pred_df <- as.vector(pred_df)[[1]]

  myRespNA <- replace(myResp, myResp == 0, NA)
  myBiomodData <- biomod2::BIOMOD_FormatingData(resp.var = if(!is.null(nsdm_regional$Background.XY.Regional)) myRespNA else myResp,
	                                     resp.xy = myResp.xy,
	                                     expl.var = myExpl,
	                                     resp.name = SpeciesName,
	                                     PA.nb.rep = ifelse(!is.null(nsdm_regional$Absences.XY.Regional), 0, 1),
	                                     PA.nb.absences = ifelse(!is.null(nsdm_regional$Absences.XY.Regional), 0, nrow(nsdm_regional$Background.XY.Regional)),
	                                     PA.strategy =if(!is.null(nsdm_regional$Absences.XY.Regional)) NULL else "random")

  calib.lines <- biomod2::bm_CrossValidation(bm.format = myBiomodData,
                                    strategy = "random",
                                    nb.rep = nsdm_regional$args$CV.nb.rep,
                                    perc =  nsdm_regional$args$CV.perc,
                                    k = 0,
                                    balance = 'presences',
                                    env.var = NULL,
                                    strat = 'both',
                                    user.table = NULL,
                                    do.full.models = FALSE)

  metric.eval<- c("ROC","TSS","KAPPA","ACCURACY", "SR", "BOYCE")
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
      stat <- biomod2::bm_FindOptimStat(metric.eval = xx,
  			obs = myResp[eval.lines.rep],
  			fit = pred_df[eval.lines.rep],
			threshold = cross.validation["cutoff", xx])
      stat.validation <- rbind(stat.validation, stat)
    }
  }

  colnames(cross.validation)[which(colnames(cross.validation) == "best.stat")] <- "calibration"
  cross.validation$validation <- stat.validation$best.stat
  myEMeval.replicates<-cross.validation

  cross.validation<-cross.validation[which(cross.validation$metric.eval %in% c("ROC","TSS","KAPPA")),]
  metric.means <- aggregate(. ~ metric.eval, data = cross.validation, FUN = mean)
  sabina$myEMeval.means<-metric.means

  # Save some results
   if(save.output){
     fs::dir_create("Results/Multiply/Values/")
     write.csv(metric.means,file=paste0("Results/Multiply/Values/",SpeciesName,"_ensemble.csv"))
   }


  # Binary models
  for(i in seq_along(Scenarios)) {
    projmodel <- Scenarios[i]

    if(projmodel =="Current") {
      continuous<-sabina$current.projections$Pred
    } else {
      continuous<-sabina$new.projections$Pred.Scenario[[i-1]]
    }
    Threshold.ROC <- metric.means$cutoff[which(metric.means$metric.eval=="ROC")]
    Threshold.TSS <- metric.means$cutoff[which(metric.means$metric.eval=="TSS")]
    Pred.bin.ROC <- terra::classify(res.average,rbind(c(0,Threshold.ROC,0),c(Threshold.ROC,2000,1)))
    Pred.bin.TSS <- terra::classify(res.average,rbind(c(0,Threshold.TSS,0),c(Threshold.TSS,2000,1)))
    Pred.bin.ROC <- terra::rast(terra::wrap(Pred.bin.ROC))
    Pred.bin.TSS <- terra::rast(terra::wrap(Pred.bin.TSS))

    if(projmodel =="Current") {
      sabina$current.projections$Pred.bin.ROC <- setNames(Pred.bin.ROC, paste0(SpeciesName, ".Current.bin.ROC"))
      sabina$current.projections$Pred.bin.TSS <- setNames(Pred.bin.TSS, paste0(SpeciesName,".Current.bin.TSS"))

      if(save.output){
        file_path<-paste0("Results/Multiply/Projections/",SpeciesName,".Current.bin.ROC.tif")
        terra::writeRaster(Pred.bin.ROC, file_path, overwrite=TRUE)
        file_path<-paste0("Results/Multiply/Projections/",SpeciesName,".Current.bin.TSS.tif")
        terra::writeRaster(Pred.bin.TSS, file_path, overwrite=TRUE)
      }

    } else {
      Scenario.name<-projmodel
      sabina$new.projections$Pred.bin.ROC.Scenario[[i]] <- setNames(Pred.bin.ROC, paste0(SpeciesName,".",Scenario.name,".bin.ROC"))
      sabina$new.projections$Pred.bin.TSS.Scenario[[i]] <- setNames(Pred.bin.TSS, paste0(SpeciesName,".",Scenario.name,".bin.TSS"))

      if(save.output){
        file_path<-paste0("Results/Multiply/Projections/",SpeciesName,".",Scenario.name,".bin.ROC.tif")
        terra::writeRaster(Pred.bin.ROC, file_path, overwrite = TRUE)

        file_path<-paste0("Results/Multiply/Projections/",SpeciesName,".",Scenario.name,".bin.TSS.tif")
        terra::writeRaster(Pred.bin.TSS, file_path, overwrite = TRUE)
      }
    }
  }

  # Summary
  summary <- data.frame(Values = c(SpeciesName,
  				method,
  				round(metric.means$validation[metric.means$metric.eval == "ROC"], 3),
  				round(metric.means$validation[metric.means$metric.eval == "TSS"], 3),
  				round(metric.means$validation[metric.means$metric.eval == "KAPPA"], 3)))

  rownames(summary) <- c("Species name",
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
    message("\nResults saved in the following local folder/s:")
    message(paste(
    "- Hierarchical Multiply Models: /Results/Multiply/Projections/\n"
    ))
  }

  return(sabina)

}

