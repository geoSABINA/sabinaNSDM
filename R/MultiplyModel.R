#' @export
HSDM.MultiplyModels <- function(hsdm_global,
                                    hsdm_regional,
                                    method="Arithmetic",
                                    rescale=FALSE,
                                    save.output=TRUE) {

  if(!inherits(hsdm_regional, "hsdm.predict.r") || !inherits(hsdm_global, "hsdm.predict.g")) {
    stop("hsdm_regional and hsdm_global must be objects of class hsdm.predict.r and hsdm.predict.r, respectively.")
  }

  if(!method %in% c("Arithmetic", "Geometric")) {
    stop("Please select 'Arithmetic' or 'Geometric' method.")
  }

  if(!identical(hsdm_regional$Species.Name, hsdm_global$Species.Name)) {
    stop("Different species in global and regional levels")
  } else {
    SpeciesName <- hsdm_regional$Species.Name
  }

  sabina<-hsdm_global[names(hsdm_global) %in% c("Species.Name")]
  sabina<-list()
  sabina$args <- list()
  sabina$args$method <- method
  sabina$args$rescale <- rescale

  Scenarios <- names(hsdm_regional$Scenarios)
  Scenarios <- c("Current",Scenarios)

  for(i in seq_along(Scenarios)) {
    projmodel <- Scenarios[i]
    # Load raster data at global and regional scales
    if(projmodel =="Current") {
      Pred.global <- terra::unwrap(hsdm_global$current.projections$Pred) # Unwrap object
      Pred.regional <- terra::unwrap(hsdm_regional$current.projections$Pred)  # Unwrap object
    } else {
        Pred.global <- terra::unwrap(hsdm_global$new.projections$Pred.Scenario[[i-1]])  # Unwrap object
        Pred.regional <- terra::unwrap(hsdm_regional$new.projections$Pred.Scenario[[i-1]])  # Unwrap object
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
  myResp.xy <- rbind(hsdm_regional$SpeciesData.XY.Regional, hsdm_regional$Background.XY.Regional)
  myResp <- data.frame(c(rep(1,nrow(hsdm_regional$SpeciesData.XY.Regional)),rep(0,nrow(hsdm_regional$Background.XY.Regional))))
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
				paste(hsdm_regional$args$algorithms,collapse = ", "),
				nrow(hsdm_regional$myEMeval.replicates), 
				hsdm_regional$myEMeval.Ensemble$calibration[which(hsdm_regional$myEMeval.Ensemble$metric.eval=="ROC")],
				hsdm_regional$myEMeval.Ensemble$calibration[which(hsdm_regional$myEMeval.Ensemble$metric.eval=="TSS")],
				hsdm_regional$myEMeval.Ensemble$calibration[which(hsdm_regional$myEMeval.Ensemble$metric.eval=="KAPPA")],
				method,
				round(valROC[,"best.stat"],3),
				round(valTSS[,"best.stat"],3),
				round(valKAPPA[,"best.stat"],3)))

  rownames(summary) <- c("Species name",
				"Statistical algorithms at regional level", 
				paste0("Number of replicates wit AUC > ",hsdm_regional$args$CV.perc," at regional level"), 
				"AUC of ensemble model at regional level", 
				"TSS of ensemble model at regional level",
				"KAPPA of ensemble model at regional level",
				"Multiply method",
				"AUC of hierarchical multiply ensemble model",
				"TSS of hierarchical multiply ensemble model",
				"KAPPA of hierarchical multiply ensemble model")

  # Wrap objects
  sabina$current.projections <- rapply(sabina$current.projections, terra::wrap, how = "list")
  if(!is.null(hsdm_selvars$Scenarios)) {
    sabina$new.projections <- rapply(sabina$new.projections, terra::wrap, how = "list")
  }

  sabina$Summary <- summary

  attr(sabina, "class") <- "hsdm.predict"

  # Logs success or error messages
  #message("\nHSDM.MultiplyModels executed successfully!\n")

  if(save.output){
  message("Results saved in the following locations:")
  message(paste(
    " - Hierarchical Multiply Models: /Results/Multiply/Projections/\n"
  ))
  }

  return(sabina)

}

