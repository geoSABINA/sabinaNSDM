#' @export
#' @import terra
#' @importFrom terra rast writeRaster app

NSH.SDM.Multiply.Models <- function(method="Arithmetic",rescale=FALSE,SpeciesName) #@@@##(CAREFUL! I added parameters speciesName, rescale, amd method, and removed  selectedModel
{
  tryCatch({ 
  #@@@## This function was highly modified
  Scenarios.temp <- list.files(paste0(VariablesPath,"/Regional"))
  Scenario <- gsub(".tif","", Scenarios.temp, fixed=TRUE)
  
 # f <- "Current"
 for(f in Scenario) {
    # Load raster dataat global scale
    Pred.global <- rast(paste("Results/Global/Geotif/", SpeciesName, ".", f, ".tif", sep=""))
    
    # Load raster at regional scale
    Pred.regional <- rast(paste("Results/Regional/Geotif/", SpeciesName, ".", f, ".tif", sep=""))
    
    # Rescale global prediction to a common range
    if (rescale==TRUE) { #@@@##( new)
      # Calculate minimum and maximum values of global prediction
      min_val <- min(values(Pred.global), na.rm = TRUE)
      max_val <- max(values(Pred.global), na.rm = TRUE)
      Pred.global <- app(Pred.global, fun = function(x) {
        ((x - min_val) / (max_val - min_val) * 999) + 1})
      # Calculate minimum and maximum values of regional prediction
      min_val <- min(values(Pred.regional), na.rm = TRUE)
      max_val <- max(values(Pred.regional), na.rm = TRUE)
      Pred.regional <- app(Pred.regional, fun = function(x) {
        ((x - min_val) / (max_val - min_val) * 999) + 1})
      
    } 
    
    #Create outoput folder
    dir_create(c("Results/Multiply/", "Results/Multiply/Geotif/"))
    
    # Average global and regional
    if (method=="Geometric") { 
      Stack.rasters <- c(Pred.global, Pred.regional) 
      Geometric.average <-  app(Stack.rasters, function(...) exp(mean(log(c(...)))))
      writeRaster(Geometric.average, paste("Results/Multiply/Geotif/",SpeciesName,".",f,".tif",sep=""), overwrite = TRUE)
      
    } else  { # method="Arithmetic"
      Stack.rasters <- c(Pred.global, Pred.regional)
      Arithmetic.average <-  mean(Stack.rasters)
      writeRaster(Arithmetic.average, paste("Results/Multiply/Geotif/",SpeciesName,".",f,".tif",sep=""), overwrite = TRUE)
    }
 }

  # Logs success or error messages 
message("NSH.SDM.Multiply.Models executed successfully")
message("Hierarchical Multiply Models saved in Results/Multiply/Images/")
}, error = function(err) {
  message("Error in NSH.SDM.Multiply.Models:", conditionMessage(err))
  return(list(result = NULL, error = err))
  }) 
} 
