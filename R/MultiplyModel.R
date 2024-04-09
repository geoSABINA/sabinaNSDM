#' @export
NSDM.MultiplyModels <- function(nsdm_global,
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
      Pred.global <- nsdm_global$current.projections$Pred
      Pred.regional <- nsdm_regional$current.projections$Pred
    } else {
        Pred.global <- nsdm_global$new.projections$Pred.Scenario[[i-1]]
        Pred.regional <- nsdm_regional$new.projections$Pred.Scenario[[i-1]]
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


  # Evaluation multiply model 	#@@@JMB Esto queda pendiente de discusión (qué 0s usamos?)
  myResp.xy <- rbind(nsdm_regional$SpeciesData.XY.Regional, nsdm_regional$Background.XY.Regional)
  myResp <- data.frame(c(rep(1,nrow(nsdm_regional$SpeciesData.XY.Regional)),rep(0,nrow(nsdm_regional$Background.XY.Regional))))
  pred_df <- sabina$current.projections$Pred
  pred_df <- terra::extract(pred_df, myResp.xy)[-1]
  myResp <- as.vector(myResp)[[1]]
  pred_df <- as.vector(pred_df)[[1]]
 
  valROC<-bm_FindOptimStat(metric.eval = "ROC",
 			obs = myResp,
 			fit = pred_df)
  valTSS<-bm_FindOptimStat(metric.eval = "TSS",
 			obs = myResp,
 			fit = pred_df)
  valKAPPA<-bm_FindOptimStat(metric.eval = "KAPPA",
 			obs = myResp,
 			fit = pred_df)

  #@@@JMB Guardar Evaluation multiply model en return?

  # Summary
  summary <- data.frame(Values = c(SpeciesName,
				paste(nsdm_regional$args$algorithms,collapse = ", "),
				nrow(nsdm_regional$myEMeval.replicates), 
				nsdm_regional$myEMeval.Ensemble$calibration[which(nsdm_regional$myEMeval.Ensemble$metric.eval=="ROC")],
				nsdm_regional$myEMeval.Ensemble$calibration[which(nsdm_regional$myEMeval.Ensemble$metric.eval=="TSS")],
				nsdm_regional$myEMeval.Ensemble$calibration[which(nsdm_regional$myEMeval.Ensemble$metric.eval=="KAPPA")],
				method,
				round(valROC[,"best.stat"],3),
				round(valTSS[,"best.stat"],3),
				round(valKAPPA[,"best.stat"],3)))

  rownames(summary) <- c("Species name",
				"Statistical algorithms at regional level", 
				paste0("Number of replicates wit AUC > ",nsdm_regional$args$CV.perc," at regional level"), 
				"AUC of ensemble model at regional level", 
				"TSS of ensemble model at regional level",
				"KAPPA of ensemble model at regional level",
				"Multiply method",
				"AUC of hierarchical multiply ensemble model",
				"TSS of hierarchical multiply ensemble model",
				"KAPPA of hierarchical multiply ensemble model")

  sabina$Summary <- summary

  attr(sabina, "class") <- "nsdm.predict"

  # Logs success or error messages
  #message("\nNSDM.MultiplyModels executed successfully!\n")

  if(save.output){
  message("Results saved in the following locations:")
  message(paste(
    " - Hierarchical Multiply Models: /Results/Multiply/Projections/\n"
  ))
  }

  return(sabina)

}

