#' @name NSDM.Multiply
#'
#' @title Perform ......
#'
#' @description This function ....
#'
#'
#' @param nsdm_global An object of class \code{nsdm.predict.g} containing global model generated using the \code{\link{NSDM.Global}} function.
#' @param nsdm_regional An object of class \code{nsdm.predict.r} containing regional model generated using the \code{\link{NSDM.Regional}} function.
#' @param method (\emph{optional, default} \code{'Arithmetic'}) \cr 
#' Options are \code{'Arithmetic'}, or \code{'Geometric'}. The \code{'Arithmetic'} method calculates the arithmetic mean of global and regional predictions, and \code{'Geometric'} calculates a geometric mean.
#' @param rescale (\emph{optional, default} \code{FALSE}) \cr 
#' An \code{logical} controls whether global and regional model predictions should be rescaled before combining them.
#' @param save.output (\emph{optional, default} \code{TRUE}) \cr 
#' A \code{logical} value defining whether the outputs should be saved at local.
#'
#'
#' @return An object of class \code{nsdm.predict} containing model information, predictions and evaluation statistics: #@@@JMB statistics pendiente de definir
#' - `$SpeciesName` Name of the species.
#' - `$args` A \code{list} containing the arguments used during modelling, including: `method`, and `rescale`.
#' - `$current.projections` A \code{list} containing: \code{Pred}, a \code{\link[terra:rast]{PackedSpatRaster}} representing the current projection.....; \code{Pred.bin.ROC}, a \code{\link[terra:rast]{PackedSpatRaster}} representing projections ..........; and \code{Pred.bin.TSS}, a \code{\link[terra:rast]{PackedSpatRaster}} representing......
#' - `$new.projections` A \code{list} containing: \code{Pred.Scenario}, the projections onto new scenarios in a \code{\link[terra:rast]{PackedSpatRaster}} format; \code{Pred.bin.ROC.Scenario}, the binary projections onto new scenarios in a \code{\link[terra:rast]{PackedSpatRaster}} format, derived from ROC scores; and \code{Pred.bin.TSS.Scenario}, the binary projections onto new scenarios in a \code{\link[terra:rast]{PackedSpatRaster}} format, derived from ROC scores.
#'
#'
#' @details
#' This function conducts .....
#' If `save.output=TRUE`, modelling results are stored out of R in the \emph{Results/} folder created in the current working directory:
#' - the \emph{Results/Multiply/Projections/} folder, containing the continious and binary current and new projections. Current projections are named with the species name followed by \file{.Current.tif}, \file{.bin.ROC.tif} and \file{.bin.TSS.tif}. New projections are named with the species name followed by the scenario name, and \file{.bin.ROC.tif}, \file{.bin.TSS.tif} when binary.
#'
#'
#' @examples 
#' #@@@JMB Ver cómo hacen otros cuando una función depende de objetos anteriores
#' myMultiplyModel <- NSDM.Multiply(myGlobalModel,
#'
#'
#' @seealso \code{\link{NSDM.InputData}}, \code{\link{NSDM.FormattingData}}, \code{\link{NSDM.SelectCovariables}}, \code{\link{NSDM.Global}}, \code{\link{NSDM.Regional}}
#'
#' 
#' @export
NSDM.Multiply <- function(nsdm_global,
                          nsdm_regional,
                          method="Arithmetic",
                          rescale=FALSE,
                          save.output=TRUE) {

  if(!inherits(nsdm_regional, "nsdm.predict.r") || !inherits(nsdm_global, "nsdm.predict.g")) {
    stop("nsdm_regional and nsdm_global must be objects of class nsdm.predict.r and nsdm.predict.r, respectively.")
  }

  if(!method %in% c("Arithmetic", "Geometric")) {
    stop("Please select 'Arithmetic' or 'Geometric' method.")
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
    if(rescale==FALSE) {
      Pred.global <- Pred.global
      Pred.regional <- Pred.regional
    } else if(rescale==TRUE) {
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
    if(method=="Geometric") {
      Stack.rasters <- c(Pred.global, Pred.regional)
      res.average <-  terra::app(Stack.rasters, function(...) exp(mean(log(c(...)))))
    } else if(method=="Arithmetic")  {
      Stack.rasters <- c(Pred.global, Pred.regional)
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
      #message(paste("Multiply projections under",Scenarios[i] ,"conditions saved in:",file_path))
    }
  } # end for


  ## Evaluation multiply model 	#@@@JMB Esto queda pendiente de discusión (qué 0s usamos?)
  #myResp.xy <- rbind(nsdm_regional$SpeciesData.XY.Regional, nsdm_regional$Background.XY.Regional)
  #myResp <- data.frame(c(rep(1,nrow(nsdm_regional$SpeciesData.XY.Regional)),rep(0,nrow(nsdm_regional$Background.XY.Regional))))
  #pred_df <- sabina$current.projections$Pred
  #pred_df <- terra::extract(pred_df, myResp.xy)[-1]
  #myResp <- as.vector(myResp)[[1]]
  #pred_df <- as.vector(pred_df)[[1]]
 
  #valROC<-bm_FindOptimStat(metric.eval = "ROC",
  #			obs = myResp,
  #			fit = pred_df)
  #valTSS<-bm_FindOptimStat(metric.eval = "TSS",
  #			obs = myResp,
  #			fit = pred_df)
  #valKAPPA<-bm_FindOptimStat(metric.eval = "KAPPA",
  #			obs = myResp,
  #			fit = pred_df)

  #@@@JMB Guardar Evaluation multiply model en return?

  # Summary
  #summary <- data.frame(Values = c(SpeciesName,
  #				paste(nsdm_regional$args$algorithms,collapse = ", "),
  #				sum(nsdm_regional$myEMeval.replicates$metric.eval == "ROC" & nsdm_regional$myEMeval.replicates$validation >= CV.perc), 
  #				nsdm_regional$myEMeval.Ensemble$calibration[which(nsdm_regional$myEMeval.Ensemble$metric.eval=="ROC")],
  #				nsdm_regional$myEMeval.Ensemble$calibration[which(nsdm_regional$myEMeval.Ensemble$metric.eval=="TSS")],
  #				nsdm_regional$myEMeval.Ensemble$calibration[which(nsdm_regional$myEMeval.Ensemble$metric.eval=="KAPPA")],
  #				method,
  #				round(valROC[,"best.stat"],3),
  #				round(valTSS[,"best.stat"],3),
  #				round(valKAPPA[,"best.stat"],3)))

  #rownames(summary) <- c("Species name",
  #				"Statistical algorithms at regional level", 
  #				paste0("Number of replicates wit AUC > ",nsdm_regional$args$CV.perc," at regional level"), 
  #				"AUC of ensemble model at regional level", 
  #				"TSS of ensemble model at regional level",
  #				"KAPPA of ensemble model at regional level",
  #				"Multiply method",
  #				"AUC of hierarchical multiply ensemble model",
  #				"TSS of hierarchical multiply ensemble model",
  #				"KAPPA of hierarchical multiply ensemble model")

  #sabina$Summary <- summary

  # Wrap objects
  sabina$current.projections <- rapply(sabina$current.projections, terra::wrap, how = "list")
  if(!is.null(nsdm_selvars$Scenarios)) {
    sabina$new.projections <- rapply(sabina$new.projections, terra::wrap, how = "list")
  }

  attr(sabina, "class") <- "nsdm.predict"

  # save.out messages
  if(save.output){
  message("Results saved in the following local folder/s:")
  message(paste(
    " - Hierarchical Multiply Models: /Results/Multiply/Projections/\n"
  ))
  }

  return(sabina)

}

