#' @export
NSH.SDM.PrepareData <- function(nshsdm_data,
				nPoints=10000,
				Min.Dist.Global="resolution",
				Min.Dist.Regional="resolution",
				save.output=TRUE) {

  if(!inherits(nshsdm_data, "nshsdm.data")){
      stop("nshsdm_data must be an object of nshsdm.data class. Consider running NSH.SDM.InputData() function.")
  }

  SpeciesName <- nshsdm_data$Species.Name

  sabina_data<-list()
  sabina_data$Species.Name <- SpeciesName
  sabina_data$args <- list()
  sabina_data$args$nPoints <- nPoints
  sabina_data$args$Min.Dist.Global <- Min.Dist.Global
  sabina_data$args$Min.Dist.Regional <- Min.Dist.Regional

  # Create directories
  if(save.output){
    dir_create(c("Results/",
		"Results/Global/",
		"Results/Global/SpeciesXY/",
		"Results/Global/Values/",
		"Results/Global/Background/"))
    dir_create(c("Results/Regional/",
		"Results/Regional/SpeciesXY/",
		"Results/Regional/Values/",
		"Results/Regional/Background/"))
  }

  # GLOBAL SCALE
  # Generate random background points for model calibration
  # from Global independent variables (environmental layers)

  # Global independent variables (environmental layers)
  IndVar.Global <- nshsdm_data$IndVar.Global
  IndVar.Global <- IndVar.Global[[names(IndVar.Global)]]
  Mask.Global <- prod(IndVar.Global, 1)
  IndVar.Global <- terra::mask(IndVar.Global, Mask.Global) 

  # Generate random background points for model calibration
  if(is.null(nshsdm_data$Background.Global.0)) {
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
    XY.Global <- terra::extract(Mask.Global, nshsdm_data$Background.Global.0) #@@@JMB xy=TRE creo que modifica las coordenadas originales
    XY.Global <- cbind(XY.Global, nshsdm_data$Background.Global.0)
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
    if(!is.null(nshsdm_data$Background.Global.0)) {
      message(paste("Global background data: from", nrow(nshsdm_data$Background.Global.0), "to", nrow(Background.XY.Global), "points after cleaning and thinning."))
    } else {
      message(paste("Global background data: from", nPoints, "to", nrow(Background.XY.Global), "points after cleaning and thinning."))
    }
  }

  if(save.output){
  write.csv(Background.XY.Global,  paste0("Results/Global/Background/Background.csv"))
  }


  # Load species data at global scale
  SpeciesData.XY.Global <- nshsdm_data$SpeciesData.XY.Global.0
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

  ## Sample size  #@@@JMB quitaría este bloque si vamos a guardar el summary
  #Sample.size.Global <- nrow(XY.final.Global)
  #if(save.output){
  #write.table(Sample.size.Global, paste0("Results/Global/Values/",SpeciesName,"_samplesize.csv", sep=""), sep=",",  row.names=F, col.names=T)
  #}

  # Summary global
  summary <- data.frame(Values = c("",SpeciesName,  #@@@JMB Terminar de definir los summary de todas las funciones. independientes? acumulado?
				nrow(SpeciesData.XY.Global), 
				nrow(XY.final.Global), 
				ifelse(is.null(nshsdm_data$Background.Global), nPoints, nrow(nshsdm_data$Background.Global)),
				nrow(Background.XY.Global)))

  rownames(summary) <- c("  - Title 1:","Species name", #@@@JMB acabar de definir partes de summary
			"Original number of species presences at global level", 
			"Final number of species presences at global level", 
			"Original number of background points at global level", 
			"Final number of background points global level")


  # REGIONAL SCALE
  # Generate random background points for model calibration
  # Regional independent variables (environmental layers)
  IndVar.Regional <- nshsdm_data$IndVar.Regional
  IndVar.Regional <- IndVar.Regional[[names(IndVar.Regional)]]
  Mask.Regional <- prod(IndVar.Regional)
  IndVar.Regional <- terra::mask(IndVar.Regional, Mask.Regional)

  # Generate random background points for model calibration
  if(is.null(nshsdm_data$Background.Regional.0)) {
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
    XY.Regional <- terra::extract(Mask.Regional, nshsdm_data$Background.Regional.0)
    XY.Regional <- cbind(XY.Regional, nshsdm_data$Background.Regional.0)
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
    if(!is.null(nshsdm_data$Background.Regional.0)) {
      message(paste("Regional background data: from", nrow(nshsdm_data$Background.Regional.0), "to", nrow(Background.XY.Regional), "points after cleaning and thinning."))
    } else {
      message(paste("Regional background data: from", nPoints, "to", nrow(Background.XY.Regional), "points after cleaning and thinning."))
    }
  }


  if(save.output){
  write.csv(Background.XY.Regional,  paste0("Results/Regional/Background/Background.csv"))
  }

  # Load species presence data at regional scale
  SpeciesData.XY.Regional <- nshsdm_data$SpeciesData.XY.Regional.0
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

  ## Sample size  #@@@JMB quitaría este bloque si vamos a guardar el summary
  #Sample.size.Regional <- nrow(XY.final.Regional)
  #if(save.output){
  #write.table(Sample.size.Regional, paste0("Results/Regional/Values/",SpeciesName,"_samplesize.csv"), sep=",",  row.names=F, col.names=T)
  #}

  # Summary regional
  summary_regional <- data.frame(Values = c(nrow(SpeciesData.XY.Regional), 
				nrow(XY.final.Regional), 
				ifelse(is.null(nshsdm_data$Background.Regional), nPoints, nrow(nshsdm_data$Background.Regional)),
				nrow(Background.XY.Regional)))

  rownames(summary_regional) <- c("Original number of species presences at regional level", 
			"Final number of species presences at regional level", 
			"Original number of background points at regional level", 
			"Final number of background points regional level")

  summary <- rbind(summary, summary_regional)

  #nScenarios <- names(nshsdm_data$Scenarios) 

  #if(length(nScenarios) == 0) { #@@@JMB parte de esto está en la función nueva NSH.SDM.InputData
  #  message("There are no new scenarios different from Current.tif")
  #} #else  {
    #message("Future scenarios: ")
    #print(path_ext_remove(path_file(Scenarios)))
  #}
  
  summary_regional <- data.frame(Values = c(length(nshsdm_data$Scenarios))) 
  rownames(summary_regional) <- c("Number of new scenarios")
  summary <- rbind(summary, summary_regional)

  #if(save.output){   #@@@JMB quitaría este bloque si vamos a guardar un summary final.
  #  write.table(summary, paste0("Results/",SpeciesName,"_summary.csv"))
  #}
  
  sabina_data$SpeciesData.XY.Global <- XY.final.Global
  sabina_data$SpeciesData.XY.Regional <- XY.final.Regional
  sabina_data$Background.XY.Global <- Background.XY.Global
  sabina_data$Background.XY.Regional <- Background.XY.Regional
  sabina_data$IndVar.Global <- IndVar.Global
  sabina_data$IndVar.Regional <- IndVar.Regional
  sabina_data$Scenarios <- nshsdm_data$Scenarios
  sabina_data$Summary<-summary

  attr(sabina_data, "class") <- "nshsdm.input"

  # Logs success or error messages
  #message("\nNSH.SDM.PrepareData() executed successfully!\n")

  if(save.output) {
    message("Results saved in the following locations:")
    message(paste0(
        #" - Global sample size: /Results/Global/Values/", SpeciesName, "_samplesize.csv\n",  #@@@JMB quitar esto si no guardamos sample size
	" - Global background points: /Results/Global/Background/Background.csv\n",
	" - Global species occurrences: /Results/Global/SpeciesXY/", SpeciesName, ".csv\n",
	#" - Regional sample size: /Results/Regional/Values/", SpeciesName, "_samplesize.csv\n",  #@@@JMB quitar esto si no guardamos sample size
	" - Regional background points: /Results/Regional/Background/Background.csv\n",
	" - Regional species occurrences: /Results/Regional/SpeciesXY/", SpeciesName, ".csv\n"
    ))
  }

  return(sabina_data)

  }


