#' @export
NSH.SDM.PrepareData <- function(VariablesPath,
				SpeciesFilePath,
				SpeciesName,
				nPoints=10000,
				Min.Dist.Global=NULL,
				Min.Dist.Regional=NULL,
				Background.Global=NULL,
				Background.Regional=NULL,
				save.output=TRUE) {

  if (!is.null(Background.Global) && !(is.data.frame(Background.Global) && ncol(Background.Global) == 2 && all(c("x", "y") %in% names(Background.Global)))) {
    stop("background must be a data.frame with two columns 'x' and 'y'.")
  }

  if (!is.null(Background.Regional) && !(is.data.frame(Background.Regional) && ncol(Background.Regional) == 2 && all(c("x", "y") %in% names(Background.Regional)))) {
    stop("background must be a data.frame with two columns 'x' and 'y'.")
  }

  sabina_data<-list()

  if(save.output){
  dir_create(c("Results/", #@@@RGM no quitar estas carpetas
               "Results/Global/", #@@@RGM no quitar estas carpetas
               "Results/Global/SpeciesXY/",
               "Results/Global/Values/",
               "Results/Global/Projections/",
               "Results/Global/Background/"))
  dir_create(c("Results/Regional/", #@@@RGM no quitar estas carpetas
              "Results/Regional/SpeciesXY/",
               "Results/Regional/Background/",
               "Results/Regional/Values/",
               "Results/Regional/Projections/")) #@@@#TG esta no la debería crear aquí
  }

  # GLOBAL SCALE
  # Generate random background points for model calibration
  # from Global independent variables (environmental layers)

  # Global independent variables (environmental layers)
  IndVar.Global <- rast(paste0(VariablesPath,"/Global/Current.tif"))
  IndVar.Global <- IndVar.Global[[names(IndVar.Global)]]
  Mask.Global <- prod(IndVar.Global, 1)
  IndVar.Global <- terra::mask(IndVar.Global, Mask.Global) #@RGM creo que esto es necesario para evitar NODATA
  # Generate random background points for model calibration
  if(is.null(Background.Global)) {
    Valid.Cells.Global <- which(!is.na(values(Mask.Global)))
    if(length(Valid.Cells.Global) < nPoints) {
      stop("The requested number of bockground points exceeds the number of valid cells.")
    }
    Sampled.indices.Global <- sample(Valid.Cells.Global, nPoints)
    Coords.Global <- terra::xyFromCell(Mask.Global, Sampled.indices.Global)
    Background.XY.Global <- as.data.frame(Coords.Global)
  } else {
    Background.XY.Global <- Background.Global
  }
  if(save.output){
  write.csv(Background.XY.Global,  paste0("Results/Global/Background/Background.csv"))
  }

# Load species data at global scale
  SpeciesData.XY.Global <- read.csv(paste0(SpeciesFilePath,"/Global/", SpeciesName, ".csv"))
  names(SpeciesData.XY.Global) <- c("x","y")

  # Occurrences from sites with no NAs
  XY.Global <- terra::extract(Mask.Global, SpeciesData.XY.Global, xy=TRUE, na.rm=TRUE)
  XY.Global <- na.omit(XY.Global)[, -c(1:2)]
  XY.Global <- unique(XY.Global)

  # Spatial thinning of presence data to remove duplicates and apply minimum distance criteria
  if(is.null(Min.Dist.Global)) {
    Min.Dist.Global<-res(Mask.Global)[1]
  }

  invisible(capture.output({
    tryCatch({
      XY.final.Global <- ecospat::ecospat.occ.desaggregation(XY.Global, min.dist = Min.Dist.Global, by = NULL)
    }, error = function(e) {
      # If an error occurs, run the alternative block
      XY.Global <- round(XY.Global, digits = 4)
      XY.final.Global <- ecospat::ecospat.occ.desaggregation(XY.Global, min.dist = Min.Dist.Global, by = NULL)
    })
  }))
  message(paste("Global data thinning: from", nrow(XY.Global), "to", nrow(XY.final.Global), "species presences"))


  # Save thinning presence data for each species
  if(save.output){
  write.csv(XY.final.Global, paste0("Results/Global/SpeciesXY/",SpeciesName,".csv"))
  }

  # Sample size
  Sample.size.Global <- nrow(XY.final.Global)
  if(save.output){
  write.table(Sample.size.Global, paste0("Results/Global/Values/",SpeciesName,"_samplesize.csv", sep=""), sep=",",  row.names=F, col.names=T)
  }

  # REGIONAL SCALE
  # Generate random background points for model calibration
  # Regional independent variables (environmental layers)
  IndVar.Regional <- terra::rast(paste0(VariablesPath,"/Regional/Current.tif"))
  IndVar.Regional <- IndVar.Regional[[names(IndVar.Regional)]]
  Mask.regional <- prod(IndVar.Regional)
  IndVar.Regional <- terra::mask(IndVar.Regional, Mask.regional) #@RGM creo que esto es necesario para evitar NODATA

  # Generate random background points for model calibration
  if(is.null(Background.Regional)) {
    Valid.Cells.Regional <- which(!is.na(values(Mask.regional)))
    if(length(Valid.Cells.Regional) < nPoints) {
      stop("The requested number of bockground points exceeds the number of valid cells.")
    }
    Sampled.indices.Regional <- sample(Valid.Cells.Regional, nPoints)
    Coords.Regional <- terra::xyFromCell(Mask.regional, Sampled.indices.Regional)
    Background.XY.Regional <- as.data.frame(Coords.Regional)
  } else {
    Background.XY.Regional <- Background.Regional
  }
  if(save.output){
  write.csv(Background.XY.Regional,  paste0("Results/Regional/Background/Background.csv"))
  }

# Load species presence data at regional scale
  SpeciesData.XY.Regional <- read.csv(paste0(SpeciesFilePath,"/Regional/", SpeciesName, ".csv"))
  names(SpeciesData.XY.Regional) <- c("x","y")

  # Occurrences from sites with no NA
  XY.Regional <- terra::extract(Mask.regional, SpeciesData.XY.Regional, xy=TRUE, na.rm=TRUE)
  XY.Regional <- na.omit(XY.Regional)[, -c(1:2)]
  XY.Regional <- unique(XY.Regional)

  # Spatial thinning of presence data to remove duplicates and apply minimum distance criteria

  if(is.null(Min.Dist.Regional)) {
  Min.Dist.Regional<-res(Mask.regional)[1]
  }
    invisible(capture.output({
      tryCatch({
        XY.final.Regional <- ecospat::ecospat.occ.desaggregation(XY.Regional, min.dist = Min.Dist.Regional, by = NULL)
      }, error = function(e) {
        # If an error occurs, run the alternative block
        XY.Regional <- unique(XY.Regional)
        XY.Regional <- round(XY.Regional, digits = 4)
        XY.final.Regional <- ecospat::ecospat.occ.desaggregation(XY.Regional, min.dist = Min.Dist.Regional, by = NULL)
      })
    }))
    message(paste("Regional data thinning: from", nrow(XY.Regional), "to", nrow(XY.final.Regional), "species presences"))


  # Save filtered presence data for each species
  if(save.output){
  write.csv(XY.final.Regional, paste0("Results/Regional/SpeciesXY/", SpeciesName, ".csv"))
  }

  # Save sample size
  Sample.size.Regional <- nrow(XY.final.Regional)
  if(save.output){
  write.table(Sample.size.Regional, paste0("Results/Regional/Values/",SpeciesName,"_samplesize.csv"), sep=",",  row.names=F, col.names=T)
  }

  sabina_data$Species.Name <- SpeciesName
  sabina_data$VariablesPath <- VariablesPath
  sabina_data$SpeciesData.XY.Global <- XY.final.Global
  sabina_data$Background.XY.Global <- Background.XY.Global
  sabina_data$SpeciesData.XY.Regional <- XY.final.Regional
  sabina_data$Background.XY.Regional <- Background.XY.Regional
  sabina_data$IndVar.Global <- IndVar.Global
  sabina_data$IndVar.Regional <- IndVar.Regional

  attr(sabina_data, "class") <- "nshsdm.input"

  return(sabina_data)

}
