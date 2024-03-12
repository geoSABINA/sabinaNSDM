#' @export
NSH.SDM.Multiply.Models <- function(method="Arithmetic",
				    rescale=FALSE,
				    SpeciesName) {

  nshsdm_data<-list()

  tryCatch({ 
  	Scenarios <- dir_ls(paste0(VariablesPath,"/Regional"))
  	Scenarios <- path_file(Scenarios) |> path_ext_remove()
	Scenarios <- Scenarios[1:3] #para probar

  	walk(Scenarios, function(projmodel) {
	  #projmodel <- Scenarios[1]
	  # Load raster dataat global scale
	  Pred.global <- terra::rast(paste0("Results/Global/Projections/",SpeciesName,".",projmodel,".tif"))
    
	  # Load raster at regional scale
	  Pred.regional <- terra::rast(paste0("Results/Regional/Projections/",SpeciesName,".",projmodel,".tif"))
    
	  # Rescale global prediction to a common range
	  if(rescale==FALSE) {
	    Pred.global <- Pred.global
	    Pred.regional <- Pred.regional
	  } else if(rescale==TRUE) {
	    # Calculate minimum and maximum values of global prediction
	    min_val <- min(values(Pred.global), na.rm = TRUE)
	    max_val <- max(values(Pred.global), na.rm = TRUE)

	    Pred.global <- terra::app(Pred.global, fun = function(x) {
	    ((x - min_val) / (max_val - min_val) * 999) + 1})
	    # Calculate minimum and maximum values of regional prediction
	    min_val <- min(values(Pred.regional), na.rm = TRUE)
	    max_val <- max(values(Pred.regional), na.rm = TRUE)

	    Pred.regional <- terra::app(Pred.regional, fun = function(x) {
	    ((x - min_val) / (max_val - min_val) * 999) + 1})
      	  } 
    
	  #Create outoput folder
    	  dir_create(c("Results/Multiply/Projections/"))
    
    	  # Average global and regional
    	  if(method=="Geometric") { 
    	    Stack.rasters <- c(Pred.global, Pred.regional) 
    	    res.average <-  terra::app(Stack.rasters, function(...) exp(mean(log(c(...)))))
    	  } else if(method=="Arithmetic")  {
    	    Stack.rasters <- c(Pred.global, Pred.regional)
    	    res.average <-  terra::mean(Stack.rasters)
	  }  
  	    
          terra::writeRaster(res.average, paste0("Results/Multiply/Projections/",SpeciesName,".",projmodel,".tif"), overwrite = TRUE)

  	}) 

  nshsdm_data$args <- list()
  nshsdm_data$args$method <- method
  nshsdm_data$args$rescale <- rescale

  nshsdm_data$link <- list()
  nshsdm_data$link$multiply.model <- "/Results/Multiply/Projections/"

  return(nshsdm_data)

  # Logs success or error messages 
  message("\nNSH.SDM.Multiply.Models executed successfully.\n")
  message("\nHierarchical Multiply Models saved in /Results/Multiply/Images/ \n")

  }, error = function(err) {
  message("\nError in NSH.SDM.Multiply.Models:", conditionMessage(err))
  return(list(result = NULL, error = err))
  }) 
} 
