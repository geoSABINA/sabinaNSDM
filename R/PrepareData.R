#' @export
NSH.SDM.PrepareData <- function(VariablesPath, SpeciesFilePath, SpeciesName, nPoints, Min.Dist.Global, Min.Dist.Regional) {
  #library(terra) ##@@@# specify the version
  #library(fs)
  #library(ecospat)
  #prueba github desktop

  nshsdm_data<-list()
  
  # Create directories to save results of the modeling process
  dir_create(c("Results/", "Results/Global/", "Results/Global/SpeciesXY/", "Results/Global/Values/", "Results/Global/Geotif/", "Results/Global/Images/", "Results/Global/Background/"))#@@@# (I dont like the name of the Geotif folder, would prefer projections or raster)
  dir_create(c("Results/Regional/", "Results/Regional/SpeciesXY/", "Results/Regional/Background/", "Results/Regional/Values/", "Results/Regional/Geotif/", "Results/Regional/Images/")) 
  #dir_create(c("Results/RegionalNoClimatic/", "Results/RegionalNoClimatic/Values/", "Results/RegionalNoClimatic/Geotif/", "Results/RegionalNoClimatic/Images/")) #@@@# REMOVED, now creating only either with or without climate data, all stored in the same Regional folder. Also,it does not create the same folders as in regional, it lacks SpeciesXY and Background
  dir_create(c("Results/Global/Models/","Results/Regional/Models/")) #Added this to store BIOMOD models
  
  
  
  # Global scale 
  # Generate random background points for model calibration
  Mask <- rast(paste0(VariablesPath,"/Global/Current.tif"))[[1]]
  Points.random <- spatSample(Mask, nPoints, cells = TRUE, na.rm = TRUE)
  Background.xy.temp <- as.data.frame(xyFromCell(Mask, Points.random[, 1]))

  # Global independent variables (environmental layers)  
  IndVar.Global <- rast(paste0(VariablesPath,"/Global/Current.tif")) 
  IndVar.Global <- IndVar.Global[[names(IndVar.Global)]] 
  
  # Remove "no data" values when extracting info from environmental layers with background points #@@@# (remove the no data points also from presences)
  Background.xy.1 <- terra::extract(IndVar.Global, Background.xy.temp)[, -1]
  Background.temp2 <- cbind(Background.xy.1, Background.xy.temp)
  Background.temp3 <- na.omit(Background.temp2)
  Background.XY.Global <- Background.temp3[, c("x", "y")]
  write.csv(Background.XY.Global,  paste0("Results/Global/Background/Background.csv"), row.names = F)
  
  # Cargar datos de presencia de especies #@@@#(English: Load species data at global scale)
  SpeciesData.XY.Global <- read.csv(paste0(SpeciesFilePath,"/Global/", SpeciesName, ".csv"))
  names(SpeciesData.XY.Global) <- c("x","y") #@@@# (remove from bkground points the ones that are occurrence points?)
  
  # Remove "no data" values when extracting info from environmental layers with species presence data
  Species.data.temp <- terra::extract(IndVar.Global, SpeciesData.XY.Global)[, -1] 
  Species.data.temp2 <- cbind(SpeciesData.XY.Global, Species.data.temp)
  Species.data.temp3  <- na.omit(Species.data.temp2)
  XY <- Species.data.temp3[, c(1, 2)]
  colnames(XY) <- c("x", "y")
  #	head(XY)
  
  # Spatial filtering of presence data to remove duplicates and apply minimum distance criteria #@@@#(thinning instead of filtering?)
  tryCatch({
    XY <- unique(XY)
    XY.final.Global <- ecospat.occ.desaggregation(XY, min.dist = Min.Dist.Global, by = NULL)
  }, error = function(e) {
    # If an error occurs, run the alternative block
    XY <- unique(XY)
    XY <- round(XY, digits = 4)
    XY.final.Global <- ecospat.occ.desaggregation(XY, min.dist = Min.Dist.Global, by = NULL)
  })
  
  
  # Save filtered presence data for each species #@@@#(thinning instead of filtering?)
  write.csv(XY.final.Global, paste0("Results/Global/SpeciesXY/",SpeciesName,".csv"), row.names = F)
  
  # Sample size 
  Sample.size.temp <- dim(XY.final.Global)
  Sample.size <- Sample.size.temp[1] # Número de parcelas con presencias de esa especie
  write.table(Sample.size, paste("Results/Global/Values/",SpeciesName,"_samplesize.csv", sep=""), sep=",",  row.names=F, col.names=T) 
  
  
  # Regional scale 
  # Generate random background points for model calibration
  Mask.regional <- rast(paste0(VariablesPath,"/Regional/Current.tif"))[[1]]
  Points.random.regional <- spatSample(Mask.regional, nPoints, cells = TRUE, na.rm = TRUE) 
  Background.xy.temp.regional <- as.data.frame(xyFromCell(Mask.regional, Points.random.regional[, 1]))
  
  # Regional independent variables (environmental layers)  
  IndVar.Regional <- rast(paste0(VariablesPath,"/Regional/Current.tif")) 
  IndVar.Regional <- IndVar.Regional[[names(IndVar.Regional)]] 
  
  # Remove "no data" values when extracting info from environmental layers with background points
  Background.xy.1.regional <- terra::extract(IndVar.Regional, Background.xy.temp.regional)[, -1]
  Background.temp2.regional <- cbind(Background.xy.1.regional, Background.xy.temp.regional)
  Background.temp3.regional <- na.omit(Background.temp2.regional)
  Background.XY.Regional <- Background.temp3.regional[, c("x", "y")]
  write.csv(Background.XY.Regional,  paste0("Results/Regional/Background/Background.csv"), row.names = F)
  
  # Load species presence data #@@@# (at regional scale)
  SpeciesData.XY.Regional <- read.csv(paste0(SpeciesFilePath,"/Regional/", SpeciesName, ".csv"))
  names(SpeciesData.XY.Regional) <- c("x","y")
  
  # Remove "no data" values when extracting info from environmental layers with species presence data
  Species.data.temp.Regional <- terra::extract(IndVar.Regional, SpeciesData.XY.Regional)[, -1]
  Species.data.temp2.Regional <- cbind(SpeciesData.XY.Regional, Species.data.temp.Regional)
  Species.data.temp3.Regional <- na.omit(Species.data.temp2.Regional)
  XY.Regional <- Species.data.temp3.Regional[, c(1, 2)]
  colnames(XY.Regional) <- c("x", "y")
  
  # Spatial filtering of presence data to remove duplicates and apply minimum distance criteria #@@@#(thinning instead of filtering?)
  tryCatch({
    XY.Regional <- unique(XY.Regional)
    XY.final.Regional <- ecospat.occ.desaggregation(XY.Regional, min.dist = Min.Dist.Global, by = NULL)
  }, error = function(e) {
    # If an error occurs, run the alternative block
    XY.Regional <- unique(XY.Regional)
    XY.Regional <- round(XY.Regional, digits = 4)
    XY.final.Regional <- ecospat.occ.desaggregation(XY.Regional, min.dist = Min.Dist.Global, by = NULL)
  })
  
  # Save filtered presence data for each species 
  write.csv(XY.final.Regional, paste0("Results/Regional/SpeciesXY/", SpeciesName, ".csv"), row.names = FALSE)
  
  # Save sample size 
  Sample.size.temp.Regional <- dim(XY.final.Regional)
  Sample.size.Regional <- Sample.size.temp.Regional[1]  #@@@#REMOVED: Número de parcelas con presencias de esa especie
  write.table(Sample.size.Regional, paste("Results/Regional/Values/",SpeciesName,"_samplesize.csv", sep=""), sep=",",  row.names=F, col.names=T) 

  nshsdm_data$SpeciesData.XY.Global <- XY.final.Global
  nshsdm_data$Background.XY.Global <- Background.XY.Global
  nshsdm_data$SpeciesData.XY.Regional <- XY.final.Regional
  nshsdm_data$Background.XY.Regional <- Background.XY.Regional

  attr(nshsdm_data, "class") <- "nshsdm.input"

  return(nshsdm_data)
}

