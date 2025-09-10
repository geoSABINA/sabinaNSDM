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
#' If no absence or background data is provided in the \code{\link{NSDM.InputData}} function, the generation method can be either \code{'random'} or \code{'stratified'}. The "random" (the default) option generates random background points considering the extension of the input environmental covariates at the global and the regional scales provided in \code{nsdm_input}. The "stratified" method is based on a PCA analysis from all environmental covariates, where the two principal component values are divided into quartiles, and multiplied to generate a total stratified variable of 16 categories (stratum). Then, the background points are generated randomly according to the area occupied by each stratum.
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
#' - `$Absences.XY.Global` Absence data at the global level at \code{data.frame} format after applying spatial thinning.
#' - `$Absences.XY.Regional` Absence data at the regional level at \code{data.frame} format after applying spatial thinning.
#' - `$IndVar.Global` Covariates at the global level in \code{\link[terra:rast]{PackedSpatRaster}} format.
#' - `$IndVar.Regional` Covariates at the regional level in \code{\link[terra:rast]{PackedSpatRaster}} format.
#' - `$Scenarios` A \code{list} containing future scenarios in \code{\link[terra:rast]{PackedSpatRaster}} format.
#' - `$Summary` Summary of formatted input data in \code{data.frame} format.
#'
#'
#' @details
#' This function formats the input data for \bold{NSDM}, including generating background points, cleaning and thinning of occurrences data (and absence data if available), and saving the results to local if specified. If `save.output=TRUE`, outputs (i.e., species occurrences after applying spatial thinning and background/absences points, at both global and regional level, are stored out of R in the \emph{Results/} folder created in the current working directory:
#' - the \emph{Results/Global/SpeciesXY/} folder, containing the occurrences data (x and y coordinates) at the global scale after applying spatial thinning, named with the \code{Species.Name} argument.
#' - the \emph{Results/Global/AbsencesXY/} folder, containing the background points (x and y coordinates) at the global scale.
#' - the \emph{Results/Regional/SpeciesXY/} folder, containing the occurrences data (x and y coordinates) at the regional scale, named with the \code{Species.Name} argument.
#' - the \emph{Results/Regional/AbsencesXY/} folder, containing the absences after applying spatial thinning or background points (x and y coordinates) at the global scale.
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
#' myInputData <- NSDM.InputData(SpeciesName = "Fagus.sylvatica",
#'				spp.data.global = Fagus.sylvatica.xy.global,
#'				spp.data.regional = Fagus.sylvatica.xy.regional,
#'				expl.var.global = expl.var.global,
#'				expl.var.regional = expl.var.regional,
#'				new.env = new.env,
#'				new.env.names = c("Scenario1"),
#'				Background.Global = NULL,
#'				Background.Regional = NULL,
#'				Absences.Global = NULL,
#'				Absences.Regional = NULL)
#'
#' # Format the input data using default parameters.
#' myFormattedData <- NSDM.FormattingData(myInputData, 
#'                                        nPoints = 1000,
#'                                        save.output=FALSE)
#'
#' summary(myFormattedData)
#' 
#' ## Format the input data specifying custom parameters.
#' # myFormattedData <- NSDM.FormattingData(
#' #				# Input data object
#' #				myInputData,
#' #				# Number of background points to generate
#' #				nPoints = 1000,
#' #				# Minimum global point distance, based on raster resolution
#' #				Min.Dist.Global = "resolution",
#' #				# Minimum regional point distance, based on raster resolution
#' #				Min.Dist.Regional = "resolution",
#' #				# Method used to generate background points, here set to 'random'
#' #				Background.method="random",
#' #				# Save the formatted data externally
#' #				save.output = TRUE)
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
    fs::dir_create(c("Results/Global/SpeciesXY/",
		"Results/Global/Values/",
		"Results/Global/AbsencesXY/"), recurse = TRUE)
    fs::dir_create(c("Results/Regional/SpeciesXY/",
		"Results/Regional/Values/",
		"Results/Regional/AbsencesXY/"), recurse = TRUE)
  }

  format_global <- gen_background_pts(nsdm_input, "Global",
                                      Background.method, nPoints, Min.Dist.Global,
                                      save.output)
  format_regional <- gen_background_pts(nsdm_input, "Regional",
                                        Background.method, nPoints, Min.Dist.Regional,
                                        save.output)

  main_summary <- rbind(format_global$Summary,
                        format_regional$Summary[-1, , drop = FALSE])

  summary_new <- data.frame(Values = c(length(nsdm_input$Scenarios)))
  rownames(summary_new) <- c("Number of new scenarios")
  main_summary <- rbind(main_summary, summary_new)

  sabina$SpeciesData.XY.Global <- format_global$SpeciesData.XY
  sabina$SpeciesData.XY.Regional <- format_regional$SpeciesData.XY
  if(!is.null(format_global$Background.XY)) {
    sabina$Background.XY.Global <- format_global$Background.XY
  } else {
    sabina$Background.XY.Global <- NULL
  }
  if(!is.null(format_regional$Background.XY)) {
    sabina$Background.XY.Regional <- format_regional$Background.XY
  } else {
    sabina$Background.XY.Regional <- NULL
  }
  if(!is.null(format_global$Absences.XY)) {
    sabina$Absences.XY.Global <- format_global$Absences.XY
  } else {
    sabina$Absences.XY.Global <- NULL
  }
  if(!is.null(format_regional$Absences.XY)) {
    sabina$Absences.XY.Regional <- format_regional$Absences.XY
  } else {
    sabina$Absences.XY.Regional <- NULL
  }
  sabina$IndVar.Global <- format_global$IndVar
  sabina$IndVar.Regional <- format_regional$IndVar
  sabina$Scenarios <- nsdm_input$Scenarios
  if(is.null(sabina$Scenarios)) {
    sabina <- c(sabina,list(Scenarios = nsdm_input$Scenarios))
  } 
  sabina$Summary <- main_summary

  attr(sabina, "class") <- "nsdm.finput"
  class(sabina) <- c("nsdm.finput", "nsdm.input")

  # save.out messages
  if(save.output) {
    message("Results saved in the following local folder/s:")
    message(paste0(
      " - Global species occurrences: /Results/Global/SpeciesXY/", SpeciesName, ".csv\n",
      if(!is.null(format_global$Background.XY)) {
        paste0(" - Global background points: /Results/Global/AbsencesXY/", SpeciesName, "_Background.csv\n")
      } else {
        ""
      },
      if(!is.null(format_global$Absences.XY)) {
        paste0(" - Global true absence points: /Results/Global/AbsencesXY/", SpeciesName, "_TrueAbsences.csv\n")
      } else {
        ""
      },
      " - Regional species occurrences: /Results/Regional/SpeciesXY/", SpeciesName, ".csv\n",
      if(!is.null(format_regional$Background.XY)) {
        paste0(" - Regional background points: /Results/Regional/AbsencesXY/", SpeciesName, "_Background.csv\n")
      } else {
        ""
      },
      if(!is.null(format_regional$Absences.XY)) {
        paste0(" - Regional true absence points: /Results/Regional/AbsencesXY/", SpeciesName, "_TrueAbsences.csv\n")
      } else {
        ""
      }
    ))
  }

  return(sabina)

}



background_stratified <- function(expl.var, nPoints) {
  if(nrow(expl.var)*ncol(expl.var) > 20000) {
    points<-terra::xyFromCell(expl.var[[1]],which(stats::complete.cases(terra::values(expl.var[[1]]))))
    indices <- sample(1:nrow(points), 20000, replace = FALSE)
    points<-points[indices,]
    df<-terra::extract(expl.var,points)
  } else {
  df <- as.data.frame(expl.var)
  }
  df <- stats::na.omit(df)
  if(nrow(df) < nPoints) {
    stop(paste("The requested number of background nPoints exceeds the number of available cells.
    Maximum number of background points:",nrow(df)))
  }
  pca <- stats::princomp(df)
  rm(df)
  PC1 <- terra::predict(expl.var, pca, index = 1)
  PC2 <- terra::predict(expl.var, pca, index = 2)

  gc()
  # Reclassify raster into 4 classes based on quartiles
  quartiles1 <- terra::global(PC1, fun = stats::quantile, na.rm = TRUE)
  cat1 <- cut(terra::values(PC1), breaks = quartiles1, labels = c(1, 2, 3, 4), include.lowest = TRUE)
  PC1_cat <- terra::setValues(PC1, cat1)
  rm(PC1)
  quartiles2 <- terra::global(PC2, fun = stats::quantile, na.rm = TRUE)
  cat2 <- cut(terra::values(PC2), breaks = quartiles2, labels = c(1, 2, 3, 4), include.lowest = TRUE)
  PC2_cat <- terra::setValues(PC2, cat2)
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

gen_background_pts <- function(nsdm_input, scale,
                               Background.method, nPoints, Min.Dist,
                               save.output){

  if(!(scale %in% c("Global", "Regional"))){
    stop("scale must be either Global or Regional")
  }

  SpeciesName <- nsdm_input$Species.Name

  lowcase_scale <- tolower(scale)
  IndVar <- terra::unwrap(nsdm_input[[paste0("IndVar.", scale)]])

  # Generate random background points for model calibration
  # Covariates (environmental layers)
  IndVar <- IndVar[[names(IndVar)]]
  Mask <- prod(IndVar)
  IndVar <- terra::mask(IndVar, Mask)

  background_scale <- nsdm_input[[paste0("Background.", scale, ".0")]]
  absences_scale <- nsdm_input[[paste0("Absences.", scale)]]
  # Generate random background points for model calibration
  if(is.null(absences_scale)) {
    if(is.null(background_scale) &&
         Background.method == "random") {
      Valid.Cells <- which(!is.na(terra::values(Mask)))
      if(length(Valid.Cells) < nPoints) {
        stop(paste("The requested number of background nPoints exceeds the number of available cells.
	  Maximum number of background points at", lowcase_scale, "level:",length(Valid.Cells)))
      }
      Sampled.indices <- sample(Valid.Cells, nPoints)
      Coords <- terra::xyFromCell(Mask, Sampled.indices)
      Background.XY <- as.data.frame(Coords)
    } else if(is.null(background_scale) && Background.method == "stratified") {
      Background.XY <- background_stratified(IndVar, nPoints=nPoints)
    } else {
      #remove NAs and duplicates
      XY <- clean_data(Mask, background_scale)
      Background.XY <- XY
    }

    if(save.output){
      utils::write.csv(Background.XY,
                paste0("Results/", scale, "/AbsencesXY/", SpeciesName, "_Background.csv"))
    }
  }

  # remove NAs and duplicates of true absences
  if(!is.null(absences_scale)) {
    XY <- clean_data(Mask, absences_scale)
    Absences.XY <- XY
  }

  # Load species ocurrence data at regional scale
  species_scale <- paste0("SpeciesData.XY.", scale, ".0")
  SpeciesData.XY <- nsdm_input[[species_scale]]

  # Occurrences from sites with no NA
  XY <- clean_data(Mask, SpeciesData.XY)

  # Spatial thinning of occurrence data to remove duplicates and apply minimum distance criteria
  if(Min.Dist == "resolution") {
    Min.Dist <- terra::res(Mask)[1]
  }

  # Convert distance to km for GeoThinneR
  #IsLonLat <- tryCatch(terra::is.lonlat(Mask), error = function(e) FALSE)
  #Min.Dist.km <- if(IsLonLat) as.numeric(Min.Dist) * 111.32 else as.numeric(Min.Dist)/1000  #@@@JMB Podemos asumir esto?
  # Convert Min.Dist (Mask CRS units) to local geodesic kilometers
  Min.Dist.km <- km_equivalent_from_mask(Mask, XY, Min.Dist)

  XY$id <- seq_len(nrow(XY)) # row id

  XY_2 <- terra::vect(XY, geom = c("x","y"), crs = terra::crs(Mask))
  XY_2 <- terra::project(XY_2, "EPSG:4326")
  XY_2 <- terra::crds(XY_2, df = TRUE)
  XY_2$id <- XY$id

  XY.thin <- GeoThinneR::thin_points(
    data = XY_2,
    lon_col = "x",
    lat_col = "y",
    method = "distance",
    thin_dist = Min.Dist.km,
    trials = 10,
    all_trials = FALSE,
    seed = 123,
    verbose = FALSE
  )
  kept_ids <- GeoThinneR::largest(XY.thin)$id
  XY.final <- XY[match(kept_ids, XY$id), , drop = FALSE]
  XY.final$id <- NULL
  XY$id <- NULL
  message(paste0(scale, " species data (",SpeciesName,"): from ",
                 nrow(SpeciesData.XY), " to ", nrow(XY.final),
                 " species occurrences after cleaning and thinning.\n"))

  # Save filtered occurrences data for each species
  if(save.output){
    utils::write.csv(XY.final, paste0("Results/", scale,
                               "/SpeciesXY/", SpeciesName, ".csv"))
  }

  # Spatial thinning of true absence data
  if(!is.null(absences_scale)) {
    if(Min.Dist == "resolution") {
      Min.Dist <- terra::res(Mask)[1]
    }

    ## Convert distance to km for GeoThinneR (absences)
    #IsLonLat <- tryCatch(terra::is.lonlat(Mask), error = function(e) FALSE)
    #Min.Dist.km <- if(IsLonLat) as.numeric(Min.Dist) * 111.32 else as.numeric(Min.Dist)/1000
    # Convert Min.Dist (Mask CRS units) to local geodesic kilometers
    Min.Dist.km <- km_equivalent_from_mask(Mask, Absences.XY, Min.Dist)   
    
    Absences.XY$id <- seq_len(nrow(Absences.XY))

    Abs_2 <- terra::vect(Absences.XY, geom = c("x","y"), crs = terra::crs(Mask))
    Abs_2 <- terra::project(Abs_2, "EPSG:4326")
    Abs_2 <- terra::crds(Abs_2, df = TRUE)
    Abs_2$id <- Absences.XY$id

    Absences.XY.thin <- GeoThinneR::thin_points(
      data = Abs_2,
      lon_col = "x",
      lat_col = "y",
      method = "distance",
      thin_dist = Min.Dist.km,
      trials = 10,
      all_trials = FALSE,
      seed = 123,
      verbose = FALSE
    )
    kept_ids_abs <- GeoThinneR::largest(Absences.XY.thin)$id
    Absences.XY.final <- Absences.XY[match(kept_ids_abs, Absences.XY$id), , drop = FALSE]
    Absences.XY.final$id <- NULL
    Absences.XY$id <- NULL
    message(paste0(scale, " absence data (",SpeciesName,"): from ",
                   nrow(Absences.XY), " to ", nrow(Absences.XY.final),
                   " absence points after cleaning and thinning.\n"))

    # Save filtered absences data for each species
    if(save.output){
    utils::write.csv(Absences.XY.final, paste0("Results/", scale,
                                        "/AbsencesXY/", SpeciesName, "_TrueAbsences.csv"))
    }
  }

  scale_summary <- data.frame(Values = c(SpeciesName,
                                         nrow(SpeciesData.XY),
                                         nrow(XY.final),
                                         ifelse(is.null(absences_scale),
                                           ifelse(is.null(background_scale),
                                                nPoints,
                                                nrow(background_scale)),
                                           "NA"),
                                        ifelse(is.null(absences_scale), "NA", nrow(Absences.XY)),
                                        ifelse(is.null(absences_scale), "NA", nrow(Absences.XY.final))))

  rownames(scale_summary) <- c("Species name",
                               paste0("Original number of species occurrences at ", lowcase_scale, " level"),
                               paste0("Final number of species occurrences at ", lowcase_scale, " level"),
                               paste0("Number of background points at ", lowcase_scale, " level"),
                               paste0("Original number of true absence points at ", lowcase_scale, " level"),
                               paste0("Final number of true absence points at ", lowcase_scale, " level"))

  return(list(
          SpeciesData.XY = XY.final,
          Background.XY = if(!is.null(Background.XY)) Background.XY else NULL,
          Absences.XY = if(!is.null(absences_scale)) Absences.XY.final else NULL,
          IndVar = terra::wrap(IndVar),
          Summary = scale_summary
  ))

}

