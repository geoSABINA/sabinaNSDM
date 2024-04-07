#' @export
NSH.SDM.Multiply.Models <- function(nshsdm_global,
                                    nshsdm_regional,
                                    method="Arithmetic",
                                    rescale=FALSE,
                                    save.output=TRUE) {

  if(!inherits(nshsdm_regional, "nshsdm.predict.r") || !inherits(nshsdm_global, "nshsdm.predict.g")) {
    stop("nshsdm_regional and nshsdm_global must be objects of class nshsdm.predict.r and nshsdm.predict.r, respectively.")
  }

  if(!method %in% c("Arithmetic", "Geometric")) {
    stop("Please select 'Arithmetic' or 'Geometric' method.")
  }

  if(!identical(nshsdm_regional$Species.Name, nshsdm_global$Species.Name)) {
    stop("Different species in global and regional levels")
  } else {
    SpeciesName <- nshsdm_regional$Species.Name
  }

  nshsdm_data<-nshsdm_global[names(nshsdm_global) %in% c("Species.Name")]
  nshsdm_data<-list()
  nshsdm_data$args <- list()
  nshsdm_data$args$method <- method
  nshsdm_data$args$rescale <- rescale

  Scenarios <- names(nshsdm_regional$Scenarios)
  Scenarios <- c("Current",Scenarios)

 for(i in 1:length(Scenarios)) {
    projmodel <- Scenarios[i]
    # Load raster data at global and regional scales
    if(projmodel =="Current") {
      Pred.global <- nshsdm_global$current.projections$Pred
      Pred.regional <- nshsdm_regional$current.projections$Pred
    } else {
      Pred.global <- nshsdm_global$new.projections$Pred.Scenario[[i-1]]
      Pred.regional <- nshsdm_regional$new.projections$Pred.Scenario[[i-1]]
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
      nshsdm_data$current.projections$Pred <- setNames(res.average, paste0(SpeciesName, ".Current"))
    } else {
      nshsdm_data$new.projections$Pred.Scenario[[i]]<- setNames(res.average, paste0(SpeciesName,".", projmodel))
    }

    if(save.output){
      #Create outoput folder
      dir_create(c("Results/Multiply/Projections/"),recurse=TRUE)
      file_path<-paste0("Results/Multiply/Projections/",SpeciesName,".",projmodel,".tif")
      terra::writeRaster(res.average, file_path, overwrite = TRUE)
      #message(paste("Multiply projections under",Scenarios[i] ,"conditions saved in:",file_path))
    }
  } # end for


  # Evaluation multiply model 	#@@@JMB Esto queda pendiente de discusión (qué 0s usamos?)
  myResp.xy <- rbind(nshsdm_regional$SpeciesData.XY.Regional, nshsdm_regional$Background.XY.Regional)
  myResp <- data.frame(c(rep(1,nrow(nshsdm_regional$SpeciesData.XY.Regional)),rep(0,nrow(nshsdm_regional$Background.XY.Regional))))
  pred_df <- nshsdm_data$current.projections$Pred
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

  # Summary
  summary <- data.frame(Values = c("",
				#SpeciesName,
				#paste(nshsdm_regional$args$algorithms,collapse = ", "), #@@@JMB si hacemos summary acumulado, todo esto ya está
				#nrow(nshsdm_regional$myEMeval.replicates), 
				#nshsdm_regional$myEMeval.Ensemble$calibration[which(nshsdm_regional$myEMeval.Ensemble$metric.eval=="ROC")],
				#nshsdm_regional$myEMeval.Ensemble$calibration[which(nshsdm_regional$myEMeval.Ensemble$metric.eval=="TSS")],
				#nshsdm_regional$myEMeval.Ensemble$calibration[which(nshsdm_regional$myEMeval.Ensemble$metric.eval=="KAPPA")],
				method,
				round(valROC[,"best.stat"],3),
				round(valTSS[,"best.stat"],3),
				round(valKAPPA[,"best.stat"],3)))

  rownames(summary) <- c("  - Title 6:",
				#"Species name",
				#"Statistical algorithms at regional level", 
				#"Number of replicates wit AUC > 0.8 at regional level", 
				#"AUC of ensemble model at regional level", 
				#"TSS of ensemble model at regional level",
				#"KAPPA of ensemble model at regional level",
				"Multiply method",
				"AUC of hierarchical multiply ensemble model",
				"TSS of hierarchical multiply ensemble model",
				"KAPPA of hierarchical multiply ensemble model")

  summary <- rbind(nshsdm_global$Summary, tail(nshsdm_regional$Summary, 6), summary)

  nshsdm_data$Summary <- summary

  #if(save.output){  #@@@JMB yo quitaría este bloque. Esta en el summary()
  #  write.table(summary, paste0("Results/",SpeciesName,"_summary.csv"), sep=",",  row.names=F, col.names=T)
  #}

  attr(nshsdm_data, "class") <- "nshsdm.predict"

  # Logs success or error messages
  #message("\nNSH.SDM.Multiply.Models executed successfully!\n")

  if(save.output){
  message("Results saved in the following locations:")
  message(paste(
    " - Hierarchical Multiply Models: /Results/Multiply/Projections/\n"
  ))
  }

  return(nshsdm_data)

}

