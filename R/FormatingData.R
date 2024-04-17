#' @name NSDM.FormatingData
#'
#' @title Prepare input data for the Hierarchical Species Distribution Modeling (NSDM) analysis.
#'
#' @description Format input data and background data for usage in \NSDM.
#'
#' @param nsdm_input An object of class "nsdm.input" generated using the \code{\link{NSDM.InputData}} function. #@@@JMB ver como ponemos el hsbm.input class
#' @param nPoints (\emph{optional, default} \code{10000}) \cr
#' An \code{integer} corresponding to the number of background points used to generate background data if absence/pseudo-absences/background points is not provided at \code{\link{NSDM.InputData}}.
#' @param Min.Dist.Global (\emph{optional, default} \code{'resolution'}) \cr
#' A \code{numeric} corresponding to the minimum distance between background points at the global level. If `Min.Dist.Global="resolution"`, the minimum distance is calculated based on the resolution of the input raster
#' @param Min.Dist.Regional (\emph{optional, default} \code{'resolution'}) \cr
#' A \code{numeric} corresponding to the minimum distance between background points at the regional level. If `Min.Dist.Regional="resolution"`, the minimum distance is calculated based on the resolution of the input raster
#' @param Background.method  (\emph{optional, default} \code{'random'}) \cr
#' If no background data is provided in the \code{\link{NSDM.InputData}} function, the generation method can be either \code{'random'} or \code{'stratified'}.
#' @param save.output (\emph{optional, default} \code{TRUE}) \cr
#' A \code{logical} value defining whether the outputs should be saved at local.  
#'
#' @return An object of class "nsdm.finput" containing formatted input data for the NSDM:
#' - `Species.Name` The name of the species provided as input.
#' - `args` A \code{list} containing the arguments used for data formatting, including: `nPoints`, `Min.Dist.Global`, `Min.Dist.Regional`, and `Background.method`.
#' - `SpeciesData.XY.Global` Species presence data at the global level at \code{data.frame} format after applying spatial thinning.
#' - `SpeciesData.XY.Regional` Species presence data at the regional level at \code{data.frame} format after applying spatial thinning.
#' - `Background.XY.Global` Background points data at the global level at \code{data.frame} format after applying spatial thinning.
#' - `Background.XY.Regional` Background points data at the regional level at \code{data.frame} format after applying spatial thinning.
#' - `IndVar.Global` Independent variables at the global level in \code{\link[terra:rast]{PackedSpatRaster}} format.
#' - `IndVar.Regional` Independent variables at the regional level in \code{\link[terra:rast]{PackedSpatRaster}} format.
#' - `Scenarios` A \code{list} containing future scenarios in \code{\link[terra:rast]{PackedSpatRaster}} format.
#' - `Summary` Summary of formated input data in \code{data.frame} format.
#'
#' @details
#' This function formates the input data for NSDM, including generating background points, cleaning and thinning presence and background data, and saving the results to local if specified. If `save.output=TRUE`, outputs (i.e., species occurrences and background points after appling spatial thinning, at both global and regional level, are stored out of R in the \emph{Results/} folder created in the current working directory:
#' - the \emph{Results/Global/SpeciesXY/} folder, containing the ocurrences species xy of global level after appling spatial thinning, named with the \code{resp.name} argument.
#' - the \emph{Results/Global/Background/} folder, containing the background points xy of global level after appling spatial thinning.
#' - the \emph{Results/Regional/SpeciesXY/} folder, containing the ocurrences species xy of global level after appling spatial thinning, named with the \code{resp.name} argument.
#' - the \emph{Results/Regional/Background/} folder, containing the background points xy of global level after appling spatial thinning.
#'
#' @examples
#' # Load the required packages  #@@@JMB en el ejemplo hay que poner también el NSDM.InputData() para tener myInputData? Ver cómo hacen otros
#' library(terra)
#' library(ecospat)
#' 
#' # Format the input data
#' myFormatedData <- NSDM.FormatingData(myInputData,
#					nPoints=1000)
#'
#' @export
NSDM.FormatingData <- function(nsdm_input,
				nPoints=10000,
				Min.Dist.Global="resolution",
				Min.Dist.Regional="resolution",
				Background.method="random", # "stratified"
				save.output=TRUE) {

  if(!inherits(nsdm_input, "nsdm.input")){
      stop("nsdm_input must be an object of nsdm.input class. Consider running NSDM.InputData() function.")
  }

  SpeciesName <- nsdm_input$Species.Name

  sabina<-list()
  sabina$Species.Name <- SpeciesName
  sabina$args <- list()
  sabina$args$nPoints <- nPoints
  sabina$args$Min.Dist.Global <- Min.Dist.Global
  sabina$args$Min.Dist.Regional <- Min.Dist.Regional
  sabina$args$Backround <- ifelse(!is.null(nsdm_input$Background.Global.0), "manually added", Background.method)

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


  # Unwrap objects if necessary
  IndVar.Global <- terra::unwrap(nsdm_input$IndVar.Global)
  IndVar.Regional <- terra::unwrap(nsdm_input$IndVar.Regional)

  # GLOBAL SCALE
  # Generate random background points for model calibration
  # from Global independent variables (environmental layers)

  # Global independent variables (environmental layers)
  IndVar.Global <- IndVar.Global[[names(IndVar.Global)]]
  Mask.Global <- prod(IndVar.Global, 1)
  IndVar.Global <- terra::mask(IndVar.Global, Mask.Global) 

  # Generate random background points for model calibration
  if(is.null(nsdm_input$Background.Global.0) && Background.method == "random") { #@@@JMB new 
    Valid.Cells.Global <- which(!is.na(values(Mask.Global)))
    if(length(Valid.Cells.Global) < nPoints) {
      stop(paste("The requested number of background nPoints exceeds the number of valid/available cells.
	Maximum number of background points at global level:",length(Valid.Cells.Global)))
    }
    Sampled.indices.Global <- sample(Valid.Cells.Global, nPoints)
    Coords.Global <- terra::xyFromCell(Mask.Global, Sampled.indices.Global)
    XY.Global <- as.data.frame(Coords.Global)
  } else if(is.null(nsdm_input$Background.Global.0) && Background.method == "stratified") {
    XY.Global <- background_stratified(IndVar.Global, nPoints=nPoints)
  } else {
    #remove NAs and duplicates of Background.Global.0
    XY.Global <- terra::extract(Mask.Global, nsdm_input$Background.Global.0) #@@@JMB xy=TRE creo que modifica las coordenadas originales
    XY.Global <- cbind(XY.Global, nsdm_input$Background.Global.0)
    XY.Global <- na.omit(XY.Global)[, -c(1:2)]
    XY.Global <- unique(XY.Global)
  }

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
  if(!is.null(nsdm_input$Background.Global.0)) {
    message(paste0("Global background data: from ", nrow(nsdm_input$Background.Global.0), " to ", nrow(Background.XY.Global), " points after cleaning and thinning."))
  } else {
    message(paste0("Global background data: from ", nPoints, " to ", nrow(Background.XY.Global), " points after cleaning and thinning."))
  }

  if(save.output){
  write.csv(Background.XY.Global,  paste0("Results/Global/Background/Background.csv"))
  }


  # Load species data at global scale
  SpeciesData.XY.Global <- nsdm_input$SpeciesData.XY.Global.0
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
  message(paste0("Global species data (",SpeciesName,"): from ", nrow(SpeciesData.XY.Global), " to ", nrow(XY.final.Global), " species presences after cleaning and thinning."))

  # Save thinning presence data for each species
  if(save.output){
  write.csv(XY.final.Global, paste0("Results/Global/SpeciesXY/",SpeciesName,".csv"))
  }

  # Summary global
  summary <- data.frame(Values = c(SpeciesName,
				nrow(SpeciesData.XY.Global), 
				nrow(XY.final.Global), 
				ifelse(is.null(nsdm_input$Background.Global), nPoints, nrow(nsdm_input$Background.Global)),
				nrow(Background.XY.Global)))

  rownames(summary) <- c("Species name",
			"Original number of species presences at global level", 
			"Final number of species presences at global level", 
			"Original number of background points at global level", 
			"Final number of background points global level")


  # REGIONAL SCALE
  # Generate random background points for model calibration
  # Regional independent variables (environmental layers)
  IndVar.Regional <- IndVar.Regional[[names(IndVar.Regional)]]
  Mask.Regional <- prod(IndVar.Regional)
  IndVar.Regional <- terra::mask(IndVar.Regional, Mask.Regional)

  # Generate random background points for model calibration
  if(is.null(nsdm_input$Background.Regional.0) && Background.method == "random") {
    Valid.Cells.Regional <- which(!is.na(values(Mask.Regional)))
    if(length(Valid.Cells.Regional) < nPoints) {
      stop(paste("The requested number of background nPoints exceeds the number of valid/available cells.
	Maximum number of background points at regional level:",length(Valid.Cells.Regional)))
    }
    Sampled.indices.Regional <- sample(Valid.Cells.Regional, nPoints)
    Coords.Regional <- terra::xyFromCell(Mask.Regional, Sampled.indices.Regional)
    XY.Regional <- as.data.frame(Coords.Regional)
  } else if(is.null(nsdm_input$Background.regional.0) && Background.method == "stratified") {
    XY.Regional <- background_stratified(IndVar.Regional, nPoints=nPoints)
  } else {
    #remove NAs and duplicates
    XY.Regional <- terra::extract(Mask.Regional, nsdm_input$Background.Regional.0)
    XY.Regional <- cbind(XY.Regional, nsdm_input$Background.Regional.0)
    XY.Regional <- na.omit(XY.Regional)[, -c(1:2)]
    XY.Regional <- unique(XY.Regional)
  }

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
    if(!is.null(nsdm_input$Background.Regional.0)) {
      message(paste0("Regional background data: from ", nrow(nsdm_input$Background.Regional.0), " to ", nrow(Background.XY.Regional), " points after cleaning and thinning."))
    } else {
      message(paste0("Regional background data: from ", nPoints, " to ", nrow(Background.XY.Regional), " points after cleaning and thinning."))
    }

  if(save.output){
  write.csv(Background.XY.Regional,  paste0("Results/Regional/Background/Background.csv"))
  }

  # Load species presence data at regional scale
  SpeciesData.XY.Regional <- nsdm_input$SpeciesData.XY.Regional.0
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
  message(paste0("Regional species data (",SpeciesName,"): from ", nrow(SpeciesData.XY.Regional), " to ", nrow(XY.final.Regional), " species presences after cleaning and thinning.\n"))
  
  # Save filtered presence data for each species
  if(save.output){
  write.csv(XY.final.Regional, paste0("Results/Regional/SpeciesXY/", SpeciesName, ".csv"))
  }

  # Summary regional
  summary_regional <- data.frame(Values = c(nrow(SpeciesData.XY.Regional), 
				nrow(XY.final.Regional), 
				ifelse(is.null(nsdm_input$Background.Regional), nPoints, nrows(nsdm_input$Background.Regional)),
				nrow(Background.XY.Regional)))

  rownames(summary_regional) <- c("Original number of species presences at regional level", 
			"Final number of species presences at regional level", 
			"Original number of background points at regional level", 
			"Final number of background points regional level")

  summary <- rbind(summary, summary_regional)

  #nScenarios <- names(nsdm_input$Scenarios) 

  #if(length(nScenarios) == 0) { #@@@JMB parte de esto está en la función nueva NSDM.InputData
  #  message("There are no new scenarios different from Current.tif")
  #} #else  {
    #message("Future scenarios: ")
    #print(path_ext_remove(path_file(Scenarios)))
  #}
  
  summary_regional <- data.frame(Values = c(length(nsdm_input$Scenarios))) 
  rownames(summary_regional) <- c("Number of new scenarios")
  summary <- rbind(summary, summary_regional)

  # Wrap objects
  IndVar.Global <- terra::wrap(IndVar.Global)
  IndVar.Regional <- terra::wrap(IndVar.Regional)
  
  sabina$SpeciesData.XY.Global <- XY.final.Global
  sabina$SpeciesData.XY.Regional <- XY.final.Regional
  sabina$Background.XY.Global <- Background.XY.Global
  sabina$Background.XY.Regional <- Background.XY.Regional
  sabina$IndVar.Global <- IndVar.Global
  sabina$IndVar.Regional <- IndVar.Regional
  sabina$Scenarios <- nsdm_input$Scenarios
  sabina$Summary<-summary

  attr(sabina, "class") <- "nsdm.finput"

  # save.out messages
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



background_stratified <- function(expl.var, nPoints) {
  vars <- as.data.frame(expl.var)
  df <- na.omit(vars)

  if(nrow(df) < nPoints) {
    stop(paste("The requested number of background nPoints exceeds the number of valid/available.
    Maximum number of background points at global level:",nrow(df)))
  }

  pca <- princomp(df)
  #dput(pca, file = "pca.csv") #@@@JMB guardamos esto?
  PC1 <- predict(expl.var, pca, index = 1)
  PC2 <- predict(expl.var, pca, index = 2)

  # Reclassify raster into 4 classes based on quartiles
  quartiles1 <- global(PC1, fun = quantile, na.rm = TRUE)
  cat1 <- cut(values(PC1), breaks = quartiles1, labels = c(1, 2, 3, 4), include.lowest = TRUE)
  PC1_cat <- setValues(PC1, cat1)
  quartiles2 <- global(PC2, fun = quantile, na.rm = TRUE)
  cat2 <- cut(values(PC2), breaks = quartiles2, labels = c(1, 2, 3, 4), include.lowest = TRUE)
  PC2_cat <- setValues(PC2, cat2)

  # Combine the 4 new categories in both PCs to create a final Stratum raster with 16 categories
  Stratum <- PC1_cat * PC2_cat
  #plot(Stratum)

  #writeRaster(Stratum, file = paste(VariablesPath, "/PC_stratums.tif", sep = "")) #@@@JMB guardamos este tif?
  
  # Create background sample
  Background <- sgsR::sample_balanced(Stratum, nPoints)

  coords <- sf::st_coordinates(Background)
  Background_df <- as.data.frame(coords)
  colnames(Background_df) <- c("x", "y")
  
  return(Background_df)

}
