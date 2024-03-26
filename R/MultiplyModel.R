#' @export
NSH.SDM.Multiply.Models <- function(nshsdm_global,
                                    nshsdm_regional,
                                    #SpeciesName, #@@@#TG charge from objects
                                    method="Arithmetic",
                                    rescale=FALSE,
                                    save.output=TRUE) {

  SpeciesName <- nshsdm_regional$Species.Name
  if(SpeciesName != nshsdm_regional$Species.Name){
    stop("Different species in global and regional levels")
  }

  if(!inherits(nshsdm_regional, "nshsdm.predict")){
    stop("nshsdm_regional must be an object of nshsdm.predict class.")
  }

  if(!inherits(nshsdm_global, "nshsdm.predict")){
    stop("nshsdm_global must be an object of nshsdm.predict class.")
  }


  nshsdm_data<-list()
  nshsdm_data$args <- list()
  nshsdm_data$args$method <- method
  nshsdm_data$args$rescale <- rescale

  nshsdm_data$link <- list()
  nshsdm_data$link$multiply.models <- "/Results/Multiply/Projections/"

  tryCatch({
  	Scenarios <- sapply(nshsdm_regional$new.projections$Pred.Scenario, names) #dir_ls(paste0(VariablesPath,"/Regional"))
  	Scenarios <- gsub(paste0(SpeciesName,"."),"",Scenarios) #path_file(Scenarios) |> path_ext_remove()
	  Scenarios <- c("Current",Scenarios)

	  for (i in 1:length(Scenarios)) {

  	#walk(Scenarios, function(projmodel) {
	  projmodel <- Scenarios[i]
	  # Load raster data at global and regional scales
	  if (projmodel =="Current") {
	    Pred.global <- nshsdm_global$current.projections$Pred #terra::rast(paste0("Results/Global/Projections/",SpeciesName,".",projmodel,".tif"))
	    Pred.regional <- nshsdm_regional$current.projections$Pred #terra::rast(paste0("Results/Regional/Projections/",SpeciesName,".",projmodel,".tif"))

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
	  if (projmodel =="Current") {
	    nshsdm_data$current.projections$Pred <- setNames(res.average, paste0(SpeciesName, ".Current"))
	  } else {
	    nshsdm_data$new.projections$Pred.Scenario[[i]]<- setNames(res.average, paste0(SpeciesName,".", projmodel))
	  }

	  if(save.output){
	    #Create outoput folder
	    dir_create(c("Results/Multiply/Projections/"),recurse=TRUE)
	    file_path<-paste0("Results/Multiply/Projections/",SpeciesName,".",projmodel,".tif")
	    terra::writeRaster(res.average, file_path, overwrite = TRUE)
	    message(paste("Multiply projections under",Scenarios[i] ,"conditions saved in:",file_path,cat("\033[1;34m")))
	    cat("\033[0m")
	  }
  #	}) # walk
} #for

  # Logs success or error messages
  message("\nNSH.SDM.Multiply.Models executed successfully!\n",cat("\033[32m"))
  cat("\033[0m")

  if(save.output){
  message("Results saved in the following locations:",cat("\033[1;34m"))
  message(paste(
    " - Hierarchical Multiply Models: /Results/Multiply/Projections/\n"
  ),cat("\033[1;34m"))
  cat("\033[0m")
  }

  }, error = function(err) {
  message("\nError in NSH.SDM.Multiply.Models:", conditionMessage(err))
  return(list(result = NULL, error = err))
  })

  results<-nshsdm_global$Summary #@@@# careful! This is not charging the summaries of global!
  results<-rbind(results,
                 c("Statistical algorithms at regional level",paste(nshsdm_regional$args$models,collapse = ", ")),
                 c("Number of replicates wit AUC > 0.8 at regional level",nrow(nshsdm_regional$myEMeval.replicates)),
                 c("AUC of ensemble modle at regional level",nshsdm_regional$myEMeval.Ensemble$calibration[which(nshsdm_regional$myEMeval.Ensemble$metric.eval=="ROC")]),
                 c("TSS of ensemble modle at regional level",nshsdm_regional$myEMeval.Ensemble$calibration[which(nshsdm_regional$myEMeval.Ensemble$metric.eval=="TSS")]),
                 c("KAPPA of ensemble modle at regional level",nshsdm_regional$myEMeval.Ensemble$calibration[which(nshsdm_regional$myEMeval.Ensemble$metric.eval=="KAPPA")]),
                 c("Multiply method",method), #@@@# we need to calculate the evaluation metrics for multiply
                 c("AUC of hierarchical multiply ensemble modle","@@@@"), #@@@# we need to calculate the evaluation metrics for multiply
                 c("TSS of hierarchical multiply ensemble modle","@@@@"),
                 c("KAPPA of hierarchical multiply ensemble modle","@@@"))

  if(save.output){
    write.table(results, paste0("Results/",SpeciesName,"_summary.csv"), sep=",",  row.names=F, col.names=T)
  }
  nshsdm_data$Summary<-results
  nshsdm_selvars$Summary<-results

  return(nshsdm_data)

}
