#' @name NSDM.FormattingData
#'
#' @title Formatting input data for spatially-nested hierarchical species distribution modeling (NSDM) analysis.
#'
#' @description Format input data and background data (if necessary) for usage in \bold{NSDM}.
#'
#'
#' @param nsdm_input An object of class \code{nsdm.input} generated using the \code{\link{NSDM.InputData}} function.
#' @param nPoints (\emph{optional, default} \code{10000}) \cr
#' An \code{integer} corresponding to the number of background points used to generate background data if absence/pseudo-absences/background points are not provided at \code{\link{NSDM.InputData}}.
#' @param Min.Dist.Global (\emph{optional, default} \code{'resolution'}) \cr
#' A \code{numeric} corresponding to the minimum distance between species occurrences at the global level. If `Min.Dist.Global="resolution"`, the minimum distance is calculated based on the resolution of the input environmental covariates at the global scale provided in \code{nsdm_input}.
#' @param Min.Dist.Regional (\emph{optional, default} \code{'resolution'}) \cr
#' A \code{numeric} corresponding to the minimum distance between species occurrences at the regional level. If `Min.Dist.Regional="resolution"`, the minimum distance is calculated based on the resolution of the input environmental covariates at the regional scale provided in \code{nsdm_input}.
#' @param Background.method  (\emph{optional, default} \code{'random'}) \cr
#' If no background data is provided in the \code{\link{NSDM.InputData}} function, the generation method can be either \code{'random'} or \code{'stratified'}. The "random" (the default) option generates random background points considering the extension of the input environmental covariates at the global and the regional scales provided in \code{nsdm_input}. The "stratified" method is based on a PCA analysis from all environmental covariates, where the two principal component values are divided into quartiles, and multiplied to generate a total stratified variable of 16 categories (stratum). Then, the background points are generated randomly according to the area occupied by each stratum.
#' @param save.output (\emph{optional, default} \code{TRUE}) \cr
#' A \code{logical} value defining whether the outputs should be saved locally.
#'
#'
#' @return
#' An object of class \code{nsdm.finput} containing formatted input data for the \bold{NSDM}:
#' - `$Species.Name` The name of the species provided as input.
#' - `$args` A \code{list} containing the arguments used for data formatting, including: `nPoints`, `Min.Dist.Global`, `Min.Dist.Regional`, and `Background.method`.
#' - `$SpeciesData.XY.Global` Species occurrences at the global level at \code{data.frame} format after applying spatial thinning.
#' - `$SpeciesData.XY.Regional` Species occurrences at the regional level at \code{data.frame} format after applying spatial thinning.
#' - `$Background.XY.Global` Background data at the global level at \code{data.frame} format.
#' - `$Background.XY.Regional` Background data at the regional level at \code{data.frame} format.
#' - `$IndVar.Global` Covariates at the global level in \code{\link[terra:rast]{PackedSpatRaster}} format.
#' - `$IndVar.Regional` Covariates at the regional level in \code{\link[terra:rast]{PackedSpatRaster}} format.
#' - `$Scenarios` A \code{list} containing future scenarios in \code{\link[terra:rast]{PackedSpatRaster}} format.
#' - `$Summary` Summary of formatted input data in \code{data.frame} format.
#'
#'
#' @details
#' This function formats the input data for \bold{NSDM}, including generating background points, cleaning and thinning of occurrences data, and saving the results to local if specified. If `save.output=TRUE`, outputs (i.e., species occurrences after applying spatial thinning and background points, at both global and regional level, are stored out of R in the \emph{Results/} folder created in the current working directory:
#' - the \emph{Results/Global/SpeciesXY/} folder, containing the occurrences data (x and y coordinates) at the global scale after applying spatial thinning, named with the \code{Species.Name} argument.
#' - the \emph{Results/Global/Background/} folder, containing the background points (x and y coordinates) at the global scale.
#' - the \emph{Results/Regional/SpeciesXY/} folder, containing the occurrences data (x and y coordinates) at the regional scale, named with the \code{Species.Name} argument.
#' - the \emph{Results/Regional/Background/} folder, containing the background points (x and y coordinates) at the global scale.
#'
#'
#' @seealso \code{\link{NSDM.InputData}}
#'
#'
#' @examples
#' library(sabinaNSDM)
#'
#' # Load species occurrences
#' data(Fagus.sylvatica.xy.global, package = "sabinaNSDM")
#' data(Fagus.sylvatica.xy.regional, package = "sabinaNSDM")
#'
#' # Load covariates
#' data(expl.var.global, package = "sabinaNSDM")
#' data(expl.var.regional, package = "sabinaNSDM")
#' expl.var.global<-terra::unwrap(expl.var.global)
#' expl.var.regional<-terra::unwrap(expl.var.regional)
#'
#' # Load new scenarios
#' data(new.env, package = "sabinaNSDM")
#' new.env<-terra::unwrap(new.env)
#'
#' # Prepare input data
#' myInputData<-NSDM.InputData(SpeciesName = "Fagus.sylvatica",
#'				spp.data.global = Fagus.sylvatica.xy.global,
#'				spp.data.regional = Fagus.sylvatica.xy.regional,
#'				expl.var.global = expl.var.global,
#'				expl.var.regional = expl.var.regional,
#'				new.env = new.env,
#'				new.env.names = c("Scenario1"),
#'				Background.Global = NULL,
#'				Background.Regional = NULL)
#'
#' # Format the input data using default parameters.
#' myFormattedData <- NSDM.FormattingData(myInputData)
#' 
#' ## Format the input data specifying custom parameters.
#' # myFormattedData <- NSDM.FormattingData(nsdm_input, 
#' #				      nPoints = 1000, # Number of background points to generate
#' #				      Min.Dist.Global = "resolution", # Minimum distance between points at the global scale, based on raster resolution
#' #				      Min.Dist.Regional = "resolution",# Minimum distance between points at the regional scale, based on raster resolution
#' #				      Background.method="random",  # Method used to generate background points, here set to 'random'
#' #				      save.output = TRUE)  	# save the formatted data externally
#' 
#'
#' @export
NSDM.FormattingData <- function(nsdm_input,
				nPoints=10000,
				Min.Dist.Global="resolution",
				Min.Dist.Regional="resolution",
				Background.method="random",
				save.output=TRUE) {

  if(!inherits(nsdm_input, "nsdm.input")){
      stop("nsdm_input must be an object of nsdm.input class. Consider running NSDM.InputData() function")
  }

  if(any(!Background.method %in% c( "random", "stratified"))) {
    stop("Please select one valid Background.method (\"random\" or \"stratified\").")
  }

  SpeciesName <- nsdm_input$Species.Name

  sabina<-list()
  sabina$Species.Name <- SpeciesName
  sabina$args <- list()
  sabina$args$nPoints <- nPoints
  sabina$args$Min.Dist.Global <- Min.Dist.Global
  sabina$args$Min.Dist.Regional <- Min.Dist.Regional
  sabina$args$Backround.method <- ifelse(!is.null(nsdm_input$Background.Global.0), "manually added", Background.method)

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

  # Unwrap objects
  IndVar.Global <- terra::unwrap(nsdm_input$IndVar.Global)
  IndVar.Regional <- terra::unwrap(nsdm_input$IndVar.Regional)

  # GLOBAL SCALE
  # Generate random background points for model calibration
  # from Global covariates (environmental layers)

  # Global covariates (environmental layers)
  IndVar.Global <- IndVar.Global[[names(IndVar.Global)]]
  Mask.Global <- prod(IndVar.Global, 1)
  IndVar.Global <- terra::mask(IndVar.Global, Mask.Global)

  # Generate random background points for model calibration
  if(is.null(nsdm_input$Background.Global.0) && Background.method == "random") {
    Valid.Cells.Global <- which(!is.na(terra::values(Mask.Global)))
    if(length(Valid.Cells.Global) < nPoints) {
      stop(paste("The requested number of background nPoints exceeds the number of available cells.
	Maximum number of background points at global level:",length(Valid.Cells.Global)))
    }
    Sampled.indices.Global <- sample(Valid.Cells.Global, nPoints)
    Coords.Global <- terra::xyFromCell(Mask.Global, Sampled.indices.Global)
    Background.XY.Global <- as.data.frame(Coords.Global)
  } else if(is.null(nsdm_input$Background.Global.0) && Background.method == "stratified") {
    Background.XY.Global <- background_stratified(IndVar.Global, nPoints=nPoints)
  } else {
    #remove NAs and duplicates of Background.Global.0
    XY.Global <- terra::extract(Mask.Global, nsdm_input$Background.Global.0)
    XY.Global <- cbind(XY.Global, nsdm_input$Background.Global.0)
    XY.Global <- na.omit(XY.Global)[, -c(1:2)]
    Background.XY.Global <- unique(XY.Global)
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

  # Spatial thinning of occurrence data to remove duplicates and apply minimum distance criteria
  if(Min.Dist.Global == "resolution" ) {
    Min.Dist.Global<-terra::res(Mask.Global)[1]
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
  message(paste0("Global species data (",SpeciesName,"): from ", nrow(SpeciesData.XY.Global), " to ", nrow(XY.final.Global), " species occurrences after cleaning and thinning"))

  # Save thinning occurrence data for each species
  if(save.output){
  write.csv(XY.final.Global, paste0("Results/Global/SpeciesXY/",SpeciesName,".csv"))
  }

  # Summary global
  summary <- data.frame(Values = c(SpeciesName,
				nrow(SpeciesData.XY.Global),
				nrow(XY.final.Global),
				ifelse(is.null(nsdm_input$Background.Global.0), nPoints, nrow(nsdm_input$Background.Global.0))))

  rownames(summary) <- c("Species name",
			"Original number of species occurrences at global level",
			"Final number of species occurrences at global level",
			"Number of background points at global level")


  # REGIONAL SCALE
  # Generate random background points for model calibration
  # Regional covariates (environmental layers)
  IndVar.Regional <- IndVar.Regional[[names(IndVar.Regional)]]
  Mask.Regional <- prod(IndVar.Regional)
  IndVar.Regional <- terra::mask(IndVar.Regional, Mask.Regional)

  # Generate random background points for model calibration
  if(is.null(nsdm_input$Background.Regional.0) && Background.method == "random") {
    Valid.Cells.Regional <- which(!is.na(terra::values(Mask.Regional)))
    if(length(Valid.Cells.Regional) < nPoints) {
      stop(paste("The requested number of background nPoints exceeds the number of available cells.
	Maximum number of background points at regional level:",length(Valid.Cells.Regional)))
    }
    Sampled.indices.Regional <- sample(Valid.Cells.Regional, nPoints)
    Coords.Regional <- terra::xyFromCell(Mask.Regional, Sampled.indices.Regional)
    Background.XY.Regional <- as.data.frame(Coords.Regional)
  } else if(is.null(nsdm_input$Background.regional.0) && Background.method == "stratified") {
    Background.XY.Regional <- background_stratified(IndVar.Regional, nPoints=nPoints)
  } else {
    #remove NAs and duplicates
    XY.Regional <- terra::extract(Mask.Regional, nsdm_input$Background.Regional.0)
    XY.Regional <- cbind(XY.Regional, nsdm_input$Background.Regional.0)
    XY.Regional <- na.omit(XY.Regional)[, -c(1:2)]
    Background.XY.Regional <- unique(XY.Regional)
  }

  if(save.output){
    write.csv(Background.XY.Regional,  paste0("Results/Regional/Background/Background.csv"))
  }

  # Load species ocurrence data at regional scale
  SpeciesData.XY.Regional <- nsdm_input$SpeciesData.XY.Regional.0

  # Occurrences from sites with no NA
  XY.Regional <- terra::extract(Mask.Regional, SpeciesData.XY.Regional)
  XY.Regional <- cbind(XY.Regional, SpeciesData.XY.Regional)
  XY.Regional <- na.omit(XY.Regional)[, -c(1:2)]
  XY.Regional <- unique(XY.Regional)

  # Spatial thinning of occurrence data to remove duplicates and apply minimum distance criteria
  if(Min.Dist.Regional=="resolution") {
  Min.Dist.Regional<-terra::res(Mask.Regional)[1]
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
  message(paste0("Regional species data (",SpeciesName,"): from ", nrow(SpeciesData.XY.Regional), " to ", nrow(XY.final.Regional), " species occurrences after cleaning and thinning.\n"))

  # Save filtered occurrences data for each species
  if(save.output){
  write.csv(XY.final.Regional, paste0("Results/Regional/SpeciesXY/", SpeciesName, ".csv"))
  }

  # Summary regional
  summary_regional <- data.frame(Values = c(nrow(SpeciesData.XY.Regional),
				nrow(XY.final.Regional),
				ifelse(is.null(nsdm_input$Background.Regional.0), nPoints, nrows(nsdm_input$Background.Regional.0))))

  rownames(summary_regional) <- c("Original number of species occurrences at regional level",
			"Final number of species ocurrences at regional level",
			"Number of background points at regional level")

  summary <- rbind(summary, summary_regional)

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
    message("Results saved in the following local folder/s:")
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
  if(nrow(expl.var)*ncol(expl.var) > 20000) {
    points<-terra::xyFromCell(expl.var[[1]],which(complete.cases(terra::values(expl.var[[1]]))))
    indices <- sample(1:nrow(points), 20000, replace = FALSE)
    points<-points[indices,]
    df<-terra::extract(expl.var,points)
  } else {
  df <- as.data.frame(expl.var)
  }
  df <- na.omit(df)
  if(nrow(df) < nPoints) {
    stop(paste("The requested number of background nPoints exceeds the number of available cells.
    Maximum number of background points:",nrow(df)))
  }
  pca <- princomp(df)
  rm(df)
  PC1 <- predict(expl.var, pca, index = 1)
  PC2 <- predict(expl.var, pca, index = 2)

  gc()
  # Reclassify raster into 4 classes based on quartiles
  quartiles1 <- global(PC1, fun = quantile, na.rm = TRUE)
  cat1 <- cut(terra::values(PC1), breaks = quartiles1, labels = c(1, 2, 3, 4), include.lowest = TRUE)
  PC1_cat <- setValues(PC1, cat1)
  rm(PC1)
  quartiles2 <- global(PC2, fun = quantile, na.rm = TRUE)
  cat2 <- cut(terra::values(PC2), breaks = quartiles2, labels = c(1, 2, 3, 4), include.lowest = TRUE)
  PC2_cat <- setValues(PC2, cat2)
  rm(PC2)

  # Combine the 4 new categories in both PCs to create a final Stratum raster with 16 categories
  Stratum <- PC1_cat * PC2_cat
  rm(PC1_cat,PC2_cat)
  # Create background sample
  Background <- sgsR::sample_balanced(Stratum, nPoints)

  coords <- sf::st_coordinates(Background)
  Background_df <- as.data.frame(coords)
  colnames(Background_df) <- c("x", "y")

  return(Background_df)

}
