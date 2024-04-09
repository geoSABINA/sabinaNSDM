#' @export
NSDM.FormatingData <- function(nsdm_inputdata,
				nPoints=10000,
				Min.Dist.Global="resolution",
				Min.Dist.Regional="resolution",
				save.output=TRUE) {

  if(!inherits(nsdm_inputdata, "nsdm.input")){
      stop("nsdm_inputdata must be an object of nsdm.input class. Consider running NSDM.InputData() function.")
  }

  SpeciesName <- nsdm_inputdata$Species.Name

  sabina<-list()
  sabina$Species.Name <- SpeciesName
  sabina$args <- list()
  sabina$args$nPoints <- nPoints
  sabina$args$Min.Dist.Global <- Min.Dist.Global
  sabina$args$Min.Dist.Regional <- Min.Dist.Regional

  # Create directories
  if(save.output){
    fs::dir_create(c("Results/",
		"Results/Global/",
		"Results/Global/SpeciesXY/",
		"Results/Global/Values/",
		"Results/Global/Background/"))
    fs::dir_create(c("Results/Regional/",
		"Results/Regional/SpeciesXY/",
		"Results/Regional/Values/",
		"Results/Regional/Background/"))
  }

  # GLOBAL SCALE
  # Generate random background points for model calibration
  # from Global independent variables (environmental layers)

  # Global independent variables (environmental layers)
  IndVar.Global <- nsdm_inputdata$IndVar.Global
  IndVar.Global <- IndVar.Global[[names(IndVar.Global)]]
  Mask.Global <- prod(IndVar.Global, 1)
  IndVar.Global <- terra::mask(IndVar.Global, Mask.Global) 

  # Generate random background points for model calibration
  if(is.null(nsdm_inputdata$Background.Global.0)) {
    Valid.Cells.Global <- which(!is.na(values(Mask.Global)))
    if(length(Valid.Cells.Global) < nPoints) {
      stop(paste("The requested number of background nPoints exceeds the number of valid/available cells.
	Maximum number of background points at global level:",length(Valid.Cells.Global)))
    }
    Sampled.indices.Global <- sample(Valid.Cells.Global, nPoints)
    Coords.Global <- terra::xyFromCell(Mask.Global, Sampled.indices.Global)
    Background.XY.Global <- as.data.frame(Coords.Global)
  } else {
    #remove NAs and duplicates
    XY.Global <- terra::extract(Mask.Global, nsdm_inputdata$Background.Global.0) #@@@JMB xy=TRE creo que modifica las coordenadas originales
    XY.Global <- cbind(XY.Global, nsdm_inputdata$Background.Global.0)
    XY.Global <- na.omit(XY.Global)[, -c(1:2)]
    XY.Global <- unique(XY.Global)
    # Spatial thinning of background data to remove duplicates and apply minimum distance criteria 
    if(Min.Dist.Global == "resolution" ) {
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
    Background.XY.Global<-XY.final.Global
    if(!is.null(nsdm_inputdata$Background.Global.0)) {
      message(paste("Global background data: from", nrow(nsdm_inputdata$Background.Global.0), "to", nrow(Background.XY.Global), "points after cleaning and thinning."))
    } else {
      message(paste("Global background data: from", nPoints, "to", nrow(Background.XY.Global), "points after cleaning and thinning."))
    }
  }

  if(save.output){
  write.csv(Background.XY.Global,  paste0("Results/Global/Background/Background.csv"))
  }


  # Load species data at global scale
  SpeciesData.XY.Global <- nsdm_inputdata$SpeciesData.XY.Global.0
  names(SpeciesData.XY.Global) <- c("x","y")

  # Occurrences from sites with no NAs
  XY.Global <- terra::extract(Mask.Global, SpeciesData.XY.Global)
  XY.Global <- cbind(XY.Global, SpeciesData.XY.Global)
  XY.Global <- na.omit(XY.Global)[, -c(1:2)]
  XY.Global <- unique(XY.Global)

  # Spatial thinning of presence data to remove duplicates and apply minimum distance criteria
  if(Min.Dist.Global == "resolution" ) {
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
  message(paste("Global species data: from", nrow(SpeciesData.XY.Global), "to", nrow(XY.final.Global), "species presences after cleaning and thinning."))

  # Save thinning presence data for each species
  if(save.output){
  write.csv(XY.final.Global, paste0("Results/Global/SpeciesXY/",SpeciesName,".csv"))
  }

  # Summary global
  summary <- data.frame(Values = c(SpeciesName,
				nrow(SpeciesData.XY.Global), 
				nrow(XY.final.Global), 
				ifelse(is.null(nsdm_inputdata$Background.Global), nPoints, nrow(nsdm_inputdata$Background.Global)),
				nrow(Background.XY.Global)))

  rownames(summary) <- c("Species name",
			"Original number of species presences at global level", 
			"Final number of species presences at global level", 
			"Original number of background points at global level", 
			"Final number of background points global level")


  # REGIONAL SCALE
  # Generate random background points for model calibration
  # Regional independent variables (environmental layers)
  IndVar.Regional <- nsdm_inputdata$IndVar.Regional
  IndVar.Regional <- IndVar.Regional[[names(IndVar.Regional)]]
  Mask.Regional <- prod(IndVar.Regional)
  IndVar.Regional <- terra::mask(IndVar.Regional, Mask.Regional)

  # Generate random background points for model calibration
  if(is.null(nsdm_inputdata$Background.Regional.0)) {
    Valid.Cells.Regional <- which(!is.na(values(Mask.Regional)))
    if(length(Valid.Cells.Regional) < nPoints) {
      stop(paste("The requested number of background nPoints exceeds the number of valid/available cells.
	Maximum number of background points at regional level:",length(Valid.Cells.Regional)))
    }
    Sampled.indices.Regional <- sample(Valid.Cells.Regional, nPoints)
    Coords.Regional <- terra::xyFromCell(Mask.Regional, Sampled.indices.Regional)
    Background.XY.Regional <- as.data.frame(Coords.Regional)
  } else {
    #remove NAs and duplicates
    XY.Regional <- terra::extract(Mask.Regional, nsdm_inputdata$Background.Regional.0)
    XY.Regional <- cbind(XY.Regional, nsdm_inputdata$Background.Regional.0)
    XY.Regional <- na.omit(XY.Regional)[, -c(1:2)]
    XY.Regional <- unique(XY.Regional)
    # Spatial thinning of background data to remove duplicates and apply minimum distance criteria
    if(Min.Dist.Regional == "resolution" ) {
      Min.Dist.Regional<-res(Mask.Regional)[1]
    }
    invisible(capture.output({
      tryCatch({
        XY.final.Regional <- ecospat::ecospat.occ.desaggregation(XY.Regional, min.dist = Min.Dist.Regional, by = NULL)
      }, error = function(e) {
        # If an error occurs, run the alternative block
        XY.Regional <- round(XY.Regional, digits = 4)
        XY.final.Regional <- ecospat::ecospat.occ.desaggregation(XY.Regional, min.dist = Min.Dist.Regional, by = NULL)
      })
    }))
    Background.XY.Regional<-XY.final.Regional
    if(!is.null(nsdm_inputdata$Background.Regional.0)) {
      message(paste("Regional background data: from", nrow(nsdm_inputdata$Background.Regional.0), "to", nrow(Background.XY.Regional), "points after cleaning and thinning."))
    } else {
      message(paste("Regional background data: from", nPoints, "to", nrow(Background.XY.Regional), "points after cleaning and thinning."))
    }
  }


  if(save.output){
  write.csv(Background.XY.Regional,  paste0("Results/Regional/Background/Background.csv"))
  }

  # Load species presence data at regional scale
  SpeciesData.XY.Regional <- nsdm_inputdata$SpeciesData.XY.Regional.0
  #names(SpeciesData.XY.Regional) <- c("x","y")

  # Occurrences from sites with no NA
  XY.Regional <- terra::extract(Mask.Regional, SpeciesData.XY.Regional)
  XY.Regional <- cbind(XY.Regional, SpeciesData.XY.Regional)
  XY.Regional <- na.omit(XY.Regional)[, -c(1:2)]
  XY.Regional <- unique(XY.Regional)

  # Spatial thinning of presence data to remove duplicates and apply minimum distance criteria
  if(Min.Dist.Regional=="resolution") {
  Min.Dist.Regional<-res(Mask.Regional)[1]
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
  message(paste("Regional species data: from", nrow(SpeciesData.XY.Regional), "to", nrow(XY.final.Regional), "species presences after cleaning and thinning.\n"))
  
  # Save filtered presence data for each species
  if(save.output){
  write.csv(XY.final.Regional, paste0("Results/Regional/SpeciesXY/", SpeciesName, ".csv"))
  }

  # Summary regional
  summary_regional <- data.frame(Values = c(nrow(SpeciesData.XY.Regional), 
				nrow(XY.final.Regional), 
				ifelse(is.null(nsdm_inputdata$Background.Regional), nPoints, nrows(nsdm_inputdata$Background.Regional)),
				nrow(Background.XY.Regional)))

  rownames(summary_regional) <- c("Original number of species presences at regional level", 
			"Final number of species presences at regional level", 
			"Original number of background points at regional level", 
			"Final number of background points regional level")

  summary <- rbind(summary, summary_regional)

  #nScenarios <- names(nsdm_inputdata$Scenarios) 

  #if(length(nScenarios) == 0) { #@@@JMB parte de esto está en la función nueva NSDM.InputData
  #  message("There are no new scenarios different from Current.tif")
  #} #else  {
    #message("Future scenarios: ")
    #print(path_ext_remove(path_file(Scenarios)))
  #}
  
  summary_regional <- data.frame(Values = c(length(nsdm_inputdata$Scenarios))) 
  rownames(summary_regional) <- c("Number of new scenarios")
  summary <- rbind(summary, summary_regional)
  
  sabina$SpeciesData.XY.Global <- XY.final.Global
  sabina$SpeciesData.XY.Regional <- XY.final.Regional
  sabina$Background.XY.Global <- Background.XY.Global
  sabina$Background.XY.Regional <- Background.XY.Regional
  sabina$IndVar.Global <- IndVar.Global
  sabina$IndVar.Regional <- IndVar.Regional
  sabina$Scenarios <- nsdm_inputdata$Scenarios
  sabina$Summary<-summary

  attr(sabina, "class") <- "nsdm.input"

  # Logs success or error messages
  #message("\nNSH.SDM.PrepareData() executed successfully!\n") #@@@JMB necesario?

  if(save.output) {
    message("Results saved in the following locations:")
    message(paste0(
	" - Global background points: /Results/Global/Background/Background.csv\n",
	" - Global species occurrences: /Results/Global/SpeciesXY/", SpeciesName, ".csv\n",
	" - Regional background points: /Results/Regional/Background/Background.csv\n",
	" - Regional species occurrences: /Results/Regional/SpeciesXY/", SpeciesName, ".csv\n"
    ))
  }

  return(sabina)

  }


