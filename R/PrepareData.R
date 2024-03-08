#' @export
NSH.SDM.PrepareData <- function(VariablesPath, 
				SpeciesFilePath, 
				SpeciesName, 
				nPoints=10000, 
				Min.Dist.Global=1, 
				Min.Dist.Regional=1,
				bckg.excluding.occu = FALSE) { #@@@JMB# Only as a suggestion. Remove from valid cells those cells with species occurrences
  nshsdm_data<-list()
  
  dir_create(c("Results/Global/SpeciesXY/", 
               "Results/Global/Values/", 
               "Results/Global/Geotif/",#@@@# (I dont like the name of the Geotif folder, would prefer projections or raster)
               "Results/Global/Images/",
               "Results/Global/Background/",
               "Results/Global/Models/"))
  dir_create(c("Results/Regional/SpeciesXY/",
               "Results/Regional/Background/",
               "Results/Regional/Values/",
               "Results/Regional/Geotif/",
               "Results/Regional/Images/",
               "Results/Regional/Models/"))
  
  # GLOBAL SCALE
  # Generate random background points for model calibration 
  # from Global independent variables (environmental layers)  
  Mask <- rast(paste0(VariablesPath,"/Global/Current.tif"))
  # expand NAs
  Mask <- prod(Mask)  
  
  # Load species data at global scale
  SpeciesData.XY.Global <- read.csv(paste0(SpeciesFilePath,"/Global/", SpeciesName, ".csv"))
  names(SpeciesData.XY.Global) <- c("x","y")

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
  #write.csv(XY.final.Global, paste0("Results/Global/SpeciesXY/",SpeciesName,".csv"), row.names = F)
  
  # Sample size 
  Sample.size.temp <- dim(XY.final.Global)
  Sample.size <- Sample.size.temp[1] # Número de parcelas con presencias de esa especie
  #write.table(Sample.size, paste("Results/Global/Values/",SpeciesName,"_samplesize.csv", sep=""), sep=",",  row.names=F, col.names=T) 

  # Generate random background points for model calibration 
  if(bckg.excluding.occu == FALSE) {
    valid_cells <- which(!is.na(values(Mask)))
    if(length(valid_cells) < nPoints) {
      stop("The requested number of bockground points exceeds the number of valid cells.")
    }
    sampled_indices <- sample(valid_cells, nPoints)
    coords <- xyFromCell(Mask, sampled_indices)
    Background.XY.Global <- as.data.frame(coords)
    #write.csv(Background.XY.Global,  paste0("Results/Global/Background/Background.csv"), row.names = F)
  } else {
    # Random background selection excluding species presence cells
    Mask2 <- Mask
    sp_cells <- terra::extract(Mask2, XY.final.Global, cells=T)$cell
    values(Mask2)[sp_cells] <- NA
    valid_cells <- which(!is.na(values(Mask2)))
    if(length(valid_cells) < nPoints) {
      stop("The requested number of bockground points exceeds the number of valid cells.")
    }
    sampled_indices <- sample(valid_cells, nPoints)
    coords <- xyFromCell(Mask2, sampled_indices)
    Background.XY.Global <- as.data.frame(coords)
    #write.csv(Background.XY.Global,  paste0("Results/Global/Background/Background.csv"), row.names = F)
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
  #write.csv(XY.final.Regional, paste0("Results/Regional/SpeciesXY/", SpeciesName, ".csv"), row.names = FALSE)
  
  # Save sample size 
  Sample.size.temp.Regional <- dim(XY.final.Regional)
  Sample.size.Regional <- Sample.size.temp.Regional[1]  #@@@#REMOVED: Número de parcelas con presencias de esa especie
  #write.table(Sample.size.Regional, paste("Results/Regional/Values/",SpeciesName,"_samplesize.csv", sep=""), sep=",",  row.names=F, col.names=T) 

  # Generate random background points for model calibration
  if(bckg.excluding.occu == FALSE) {
    valid_cells <- which(!is.na(values(Mask)))
    if(length(valid_cells) < nPoints) {
      stop("The requested number of bockground points exceeds the number of valid cells.")
    }
    sampled_indices <- sample(valid_cells, nPoints)
    coords <- xyFromCell(Mask, sampled_indices)
    Background.XY.Regional <- as.data.frame(coords)
    #write.csv(Background.XY.Regional,  paste0("Results/Regional/Background/Background.csv"), row.names = F)
  } else {
    # Random background selection excluding spcies presence locations
    Mask.regional2 <- Mask.regional
    sp_cells.regional <- terra::extract(Mask.regional2, XY.final.Global, cells=T)$cell
    values(Mask.regional2)[sp_cells.regional] <- NA
    valid_cells <- which(!is.na(values(Mask2)))
    if(length(valid_cells) < nPoints) {
      stop("The requested number of bockground points exceeds the number of valid cells.")
    }
    sampled_indices <- sample(valid_cells, nPoints)
    coords <- xyFromCell(Mask, sampled_indices)
    Background.XY.Regional <- as.data.frame(coords)
    #write.csv(Background.XY.Regional,  paste0("Results/Regional/Background/Background.csv"), row.names = F)
  }

  nshsdm_data$Species.Name <- SpeciesName
  nshsdm_data$SpeciesData.XY.Global <- XY.final.Global
  nshsdm_data$Background.XY.Global <- Background.XY.Global
  nshsdm_data$SpeciesData.XY.Regional <- XY.final.Regional
  nshsdm_data$Background.XY.Regional <- Background.XY.Regional
  nshsdm_data$Sample.size.Global <- Sample.size
  nshsdm_data$Sample.size.Regional <- Sample.size.Regional

  attr(nshsdm_data, "class") <- "nshsdm.input"

  return(nshsdm_data)
 
}
