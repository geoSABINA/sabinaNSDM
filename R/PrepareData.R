#' @export
NSH.SDM.PrepareData <- function(VariablesPath, 
				SpeciesFilePath, 
				SpeciesName, nPoints, 
				Min.Dist.Global=1, 
				Min.Dist.Regional=1) {
nshsdm_data<-list()
  
  dir_create(c("Results/Global/SpeciesXY/", 
               "Results/Global/Values/", 
               "Results/Global/Geotif/",
               "Results/Global/Images/",
               "Results/Global/Background/",
               "Results/Global/Models/"))#@@@# (I dont like the name of the Geotif folder, would prefer projections or raster)
  dir_create(c("Results/Regional/SpeciesXY/",
               "Results/Regional/Background/",
               "Results/Regional/Values/",
               "Results/Regional/Geotif/",
               "Results/Regional/Images/",
               "Results/Regional/Models/"))
  
  # GLOBAL SCALE 
  # Generate random background points for model calibration
  Mask <- rast(paste0(VariablesPath,"/Global/Current.tif"))
  Mask <- prod(Mask)  # propaga los NAs

  # Global independent variables (environmental layers)  
  IndVar.Global <- rast(paste0(VariablesPath,"/Global/Current.tif")) 
  IndVar.Global <- IndVar.Global[[names(IndVar.Global)]] 
  
  # Load species data at global scale
  SpeciesData.XY.Global <- read.csv(paste0(SpeciesFilePath,"/Global/", SpeciesName, ".csv"))
  names(SpeciesData.XY.Global) <- c("x","y") #@@@# (remove from bkground points the ones that are occurrence points?)

  # Occurrences from sites with no NAs
  XY <- terra::extract(Mask, SpeciesData.XY.Global, xy=TRUE, na.rm=TRUE)
  XY <- na.omit(XY)[, -c(1:2)]

  # Spatial thinning of presence data to remove duplicates and apply minimum distance criteria
  resGlobal<-res(Mask)[1]
  invisible(capture.output({
    tryCatch({
      XY <- unique(XY)
      XY.final.Global <- ecospat.occ.desaggregation(XY, min.dist = Min.Dist.Global*resGlobal, by = NULL)
    }, error = function(e) {
      # If an error occurs, run the alternative block
      XY <- unique(XY)
      XY <- round(XY, digits = 4)
      XY.final.Global <- ecospat.occ.desaggregation(XY, min.dist = Min.Dist.Global*resGlobal, by = NULL)
    })
  }))
  message(paste("Global data thinning: from", dim(XY), "to", dim(XY.final.Global), "species presences"))

  # Save thinning presence data for each species
  write.csv(XY.final.Global, paste0("Results/Global/SpeciesXY/",SpeciesName,".csv"), row.names = F)
  
  # Sample size 
  Sample.size.temp <- dim(XY.final.Global)
  Sample.size <- Sample.size.temp[1] # Número de parcelas con presencias de esa especie
  write.table(Sample.size, paste("Results/Global/Values/",SpeciesName,"_samplesize.csv", sep=""), sep=",",  row.names=F, col.names=T) 

  # Generate random background points for model calibration #@@@JMB# Only as a suggestion
  if(bckg.excluding.occu == FALSE) {
    Background.XY.Global <- spatSample(Mask, nPoints, cell= TRUE, replace = FALSE, xy = TRUE, na.rm = TRUE, values = TRUE, as.df = TRUE)[,-1]
    Background.XY.Global <- na.omit(Background.XY.Global)[,c("x","y")]
    #Background.XY.Global <- as.data.frame(Background.XY.Global) #no es df si values=FALSE
    #Background.XY.Global <- Background.XY.Global[,2:3] #@@@JMB# comprobar que no hay duplicados si quito cell=T arriba
    write.csv(Background.XY.Global,  paste0("Results/Global/Background/Background.csv"), row.names = F)
  } else {
    # Random background selection excluding species presence cells
    Mask2 <- Mask
    sp_cells <- terra::extract(Mask2, XY.final.Global, cells=T)$cell
    values(Mask2)[sp_cells] <- NA
    Background.XY.Global <- spatSample(Mask2, nPoints, cell= TRUE, replace = FALSE, xy = TRUE, na.rm = TRUE, values = TRUE, as.df = TRUE)[,-1] 
    Background.XY.Global <- na.omit(Background.XY.Global)[,c("x","y")]
    #Background.XY.Global <- as.data.frame(Background.XY.Global)
    #Background.XY.Global <- Background.XY.Global[,2:3] #@@@JMB# comprobar que no hay duplicados si quito cell=T arriba
    write.csv(Background.XY.Global,  paste0("Results/Global/Background/Background.csv"), row.names = F)
  }

  # REGIONAL SCALE 
  # Generate random background points for model calibration
  Mask.regional <- rast(paste0(VariablesPath,"/Regional/Current.tif"))
  Mask.regional <- prod(Mask.regional)

  # Load species presence data at regional scale
  SpeciesData.XY.Regional <- read.csv(paste0(SpeciesFilePath,"/Regional/", SpeciesName, ".csv"))
  names(SpeciesData.XY.Regional) <- c("x","y")

  # Occurrences from sites with no NA
  XY.Regional <- terra::extract(Mask.regional, SpeciesData.XY.Regional, xy=TRUE, na.rm=TRUE)
  XY.Regional <- na.omit(XY.Regional)[, -c(1:2)]

  # Spatial thinning of presence data to remove duplicates and apply minimum distance criteria
  resRegional<-res(Mask.regional)[1] 
  invisible(capture.output({
    tryCatch({
      XY.Regional <- unique(XY.Regional)
      XY.final.Regional <- ecospat.occ.desaggregation(XY.Regional, min.dist = Min.Dist.Regional*resRegional, by = NULL)
    }, error = function(e) {
      # If an error occurs, run the alternative block
      XY.Regional <- unique(XY.Regional)
      XY.Regional <- round(XY.Regional, digits = 4)
      XY.final.Regional <- ecospat.occ.desaggregation(XY.Regional, min.dist = Min.Dist.Regional*resRegional, by = NULL)
    })
  }))
  message(paste("Regional data thinning: from", dim(XY.Regional), "to", dim(XY.final.Regional), "species presences"))

  # Save filtered presence data for each species 
  write.csv(XY.final.Regional, paste0("Results/Regional/SpeciesXY/", SpeciesName, ".csv"), row.names = FALSE)
  
  # Save sample size 
  Sample.size.temp.Regional <- dim(XY.final.Regional)
  Sample.size.Regional <- Sample.size.temp.Regional[1]  #@@@#REMOVED: Número de parcelas con presencias de esa especie
  write.table(Sample.size.Regional, paste("Results/Regional/Values/",SpeciesName,"_samplesize.csv", sep=""), sep=",",  row.names=F, col.names=T) 

  # Generate random background points for model calibration #@@@JMB# only as a suggestion
  if(bckg.excluding.occu == FALSE) {
    Background.XY.Regional <- spatSample(Mask.regional, nPoints, cell= TRUE, replace = FALSE, xy = TRUE, na.rm = TRUE, values = TRUE, as.df = TRUE)[,-1]
    Background.XY.Regional <- na.omit(Background.XY.Regional)[,c("x","y")]
    #Background.XY.Regional <- as.data.frame(Background.XY.Regional)
    write.csv(Background.XY.Regional,  paste0("Results/Regional/Background/Background.csv"), row.names = F)
  } else {
    # Random background selection excluding spcies presence locations
    Mask.regional2 <- Mask.regional
    sp_cells.regional <- terra::extract(Mask.regional2, XY.final.Global, cells=T)$cell
    values(Mask.regional2)[sp_cells.regional] <- NA
    Background.XY.Regional <- spatSample(Mask.regional2, nPoints, cell= TRUE, replace = FALSE, xy = TRUE, na.rm = TRUE, values = TRUE, as.df = TRUE)[,-1]
    Background.XY.Regional <- na.omit(Background.XY.Regional)[,c("x","y")]
    #Background.XY.Regional <- as.data.frame(Background.XY.Regional)
    write.csv(Background.XY.Regional,  paste0("Results/Regional/Background/Background.csv"), row.names = F)
  }

  nshsdm_data$SpeciesData.XY.Global <- XY.final.Global
  nshsdm_data$Background.XY.Global <- Background.XY.Global
  nshsdm_data$SpeciesData.XY.Regional <- XY.final.Regional
  nshsdm_data$Background.XY.Regional <- Background.XY.Regional

  attr(nshsdm_data, "class") <- "nshsdm.input"

  return(nshsdm_data)
 
}
