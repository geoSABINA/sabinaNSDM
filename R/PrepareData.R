#' @export
NSH.SDM.PrepareData <- function(VariablesPath,
				SpeciesFilePath,
				SpeciesName,
				nPoints=10000,
				Min.Dist.Global="Resolution",
				Min.Dist.Regional="Resolution",
				Background.Global=NULL,
				Background.Regional=NULL,
				save.output=TRUE) {

  if (!is.null(Background.Global) && !(is.data.frame(Background.Global) && ncol(Background.Global) == 2 && all(c("x", "y") %in% names(Background.Global)))) {
    stop("Background.Global must be a data.frame with two columns 'x' and 'y' or NULL.")
  }

  if (!is.null(Background.Regional) && !(is.data.frame(Background.Regional) && ncol(Background.Regional) == 2 && all(c("x", "y") %in% names(Background.Regional)))) {
    stop("Background.Regional must be a data.frame with two columns 'x' and 'y' or NULL.")
  }
  results<-data.frame(Name="Species name",Value=SpeciesName)
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
      stop(paste("The requested number of bockground points exceeds the number of valid cells. Maximum number of background points:",length(Valid.Cells.Global)))
    }
    Sampled.indices.Global <- sample(Valid.Cells.Global, nPoints)
    Coords.Global <- terra::xyFromCell(Mask.Global, Sampled.indices.Global)
    Background.XY.Global <- as.data.frame(Coords.Global)
  } else {
    #remova NAs
    XY.Global <- terra::extract(Mask.Global, Background.Global, xy=TRUE, na.rm=TRUE)
    XY.Global <- na.omit(XY.Global)[, -c(1:2)]
    XY.Global <- unique(XY.Global)
    # Spatial thinning of background data to remove duplicates and apply minimum distance criteria
    if(Min.Dist.Global == "Resolution" ) {
      Min.Dist.Global<-res(Mask.Global)[1]}
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
    message(paste("Global background data thinning: from", nrow(Background.Global), "to", nrow(Background.XY.Global)))
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
  if(Min.Dist.Global == "Resolution" ) {
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
results<-rbind(results,
               c("Original number of species presences at global level",nrow(XY.Global)),
               c("Final number of species presences at global level",nrow(XY.final.Global)),
               c("Original number of background points at global level",nrow(Background.Global)),
               c("Final number of background points global level",nrow(Background.XY.Global)))
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
  Mask.Regional <- prod(IndVar.Regional)
  IndVar.Regional <- terra::mask(IndVar.Regional, Mask.Regional) #@RGM creo que esto es necesario para evitar NODATA

  # Generate random background points for model calibration
  if(is.null(Background.Regional)) {
    Valid.Cells.Regional <- which(!is.na(values(Mask.Regional)))
    if(length(Valid.Cells.Regional) < nPoints) {
      stop("The requested number of bockground points exceeds the number of valid cells.")
    }
    Sampled.indices.Regional <- sample(Valid.Cells.Regional, nPoints)
    Coords.Regional <- terra::xyFromCell(Mask.Regional, Sampled.indices.Regional)
    Background.XY.Regional <- as.data.frame(Coords.Regional)
  } else {
    #remova NAs
    XY.Regional <- terra::extract(Mask.Regional, Background.Regional, xy=TRUE, na.rm=TRUE)
    XY.Regional <- na.omit(XY.Regional)[, -c(1:2)]
    XY.Regional <- unique(XY.Regional)
    # Spatial thinning of background data to remove duplicates and apply minimum distance criteria
    if(Min.Dist.Regional == "Resolution" ) {
      Min.Dist.Regional<-res(Mask.Regional)[1]}
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
    message(paste("Regional background data thinning: from", nrow(Background.Regional), "to", nrow(Background.XY.Regional)))
  }

  if(save.output){
  write.csv(Background.XY.Regional,  paste0("Results/Regional/Background/Background.csv"))
  }

# Load species presence data at regional scale
  SpeciesData.XY.Regional <- read.csv(paste0(SpeciesFilePath,"/Regional/", SpeciesName, ".csv"))
  names(SpeciesData.XY.Regional) <- c("x","y")

  # Occurrences from sites with no NA
  XY.Regional <- terra::extract(Mask.Regional, SpeciesData.XY.Regional, xy=TRUE, na.rm=TRUE)
  XY.Regional <- na.omit(XY.Regional)[, -c(1:2)]
  XY.Regional <- unique(XY.Regional)

  # Spatial thinning of presence data to remove duplicates and apply minimum distance criteria

  if(Min.Dist.Regional=="Resolution") {
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
    message(paste("Regional data thinning: from", nrow(XY.Regional), "to", nrow(XY.final.Regional), "species presences"))

     results<-rbind(results,
                   c("Original number of species presences at regional level",nrow(XY.Regional)),
                   c("Final number of species presences at regional level",nrow(XY.final.Regional)),
                   c("Original number of background points at regional level",nrow(Background.Regional)),
                   c("Final number of background points regional level",nrow(Background.XY.Regional)))
  # Save filtered presence data for each species

  if(save.output){
  write.csv(XY.final.Regional, paste0("Results/Regional/SpeciesXY/", SpeciesName, ".csv"))
  }

  # Save sample size
  Sample.size.Regional <- nrow(XY.final.Regional)
  if(save.output){
  write.table(Sample.size.Regional, paste0("Results/Regional/Values/",SpeciesName,"_samplesize.csv"), sep=",",  row.names=F, col.names=T)
  }

  Scenarios <- dir_ls(paste0(VariablesPath,"/Regional"), pattern="tif") #@@@# change this so it comes from an object
  Scenarios <- Scenarios[!grepl("Current.tif", Scenarios)]
  if(length(Scenarios) == 0) {
    message("There are no future scenarios",	cat("\033[0m"))
  } else  {
    message("Future scenarios: ",cat("\033[1;34m"))
    print(path_ext_remove(path_file(Scenarios)),cat("\033[1;34m"))
    cat("\033[0m")
  }
  results<-rbind(results,
                 c("Number of new scenarios",length(Scenarios)),
                 c("Variables path",VariablesPath))
  if(save.output){
    write.table(results, paste0("Results/",SpeciesName,"_summary.csv"), sep=",",  row.names=F, col.names=T)
  }
  sabina_data$Summary<-results
  sabina_data$Species.Name <- SpeciesName
  sabina_data$VariablesPath <- VariablesPath
  sabina_data$SpeciesData.XY.Global <- XY.final.Global
  sabina_data$Background.XY.Global <- Background.XY.Global
  sabina_data$SpeciesData.XY.Regional <- XY.final.Regional
  sabina_data$Background.XY.Regional <- Background.XY.Regional
  sabina_data$IndVar.Global <- IndVar.Global
  sabina_data$IndVar.Regional <- IndVar.Regional
  sabina_data$Scenarios<-	Scenarios

  attr(sabina_data, "class") <- "nshsdm.input"

   # Logs success or error messages
  # message("\nPrepareData executed successfully!\n",cat("\033[32m"))
  # cat("\033[0m")
  # if(save.output){
  #   message("Results saved in the following locations:",cat("\033[1;34m"))
  #   message(paste0(
  #     " - Global sample size: Results/Global/Values/\n",
  #     " - Global background points: Results/Global/Background/Background.csv\n",
  #     " - Global species occurrences: Results/Global/SpeciesXY/\n",
  #     " - Regional sample size: Results/Regional/Values/\n",
  #     " - Regional background points: Results/Regional/Background/Background.csv\n",
  #     " - Regional species occurrences: Results/Regional/SpeciesXY/\n"
  #   ),cat("\033[1;34m"))
  #   cat("\033[0m")
  #}

  return(sabina_data)

  }

