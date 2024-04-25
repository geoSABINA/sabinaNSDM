#' @name NSDM.SelectCovariates
#'
#' @title Select covariates for spatially-nested hierarchical species distribution modeling (NSDM) analysis.
#'
#' @description This function selects the best 'non-colinear' environmental covariates for \bold{NSDM} based on specified criteria and algorithms.
#'
#'
#' @param nsdm_finput An object of class \code{nsdm.finput} generated using the \code{\link{NSDM.FormatingData}} function. #@@@JMB ver como ponemos el hsbm.finput class
#' @param maxncov.Global (\emph{optional, default} \code{'nocorr'}) \cr
#' Maximum \code{numeric} value indicating the maximum number of covariates to select at the global scale. If `"nocorr"`, selects all non-correlated covariates.
#' @param maxncov.Regional (\emph{optional, default} \code{'nocorr'}) \cr
#' Maximum \code{numeric} value indicating the maximum number of covariates to select at the regional scale. If `"nocorr"`, selects all non-correlated covariates.
#' @param corcut (\emph{optional, default} \code{0.7}) \cr
#' A \code{numeric} value for the correlation coefficient threshold used for identifying collinearity.
#' @param algorithms (\emph{optional, default} \code{'c("glm", "gam", "rf")'}) \cr
#' Algorithms to use for ranking the covariates. Options are \code{'glm'}, \code{'gam'}, and/or \code{'rf'}.
#' @param ClimaticVariablesBands (\emph{optional, default} \code{NULL}) \cr
#' Indices of the regional environmental covariate bands to exclude from the selection. If \code{NULL} (the default), all regional-level covariates are considered. The excluded covariates typically include climatic covariates already included in global-level analyses.
#' A \code{logical} value defining whether the outputs should be saved at local.
#'
#'
#' @return An object of class \code{nsdm.vinput} containing selected covariates for \bold{NSDM}: #@@@JMB ver como ponemos el hsbm.input class
#' - `$SpeciesName` Name of the species.
#' - `$args` A \code{list} containing the arguments used during covariates selection procedure, including: `maxncov.Global`, `maxncov.Regional`, `corcut` and `algorithms`.
#' - `$SpeciesData.XY.Global` Species presence data at the global level at \code{data.frame} format after applying spatial thinning.
#' - `$SpeciesData.XY.Regional` Species presence data at the regional level at \code{data.frame} format after applying spatial thinning.
#' - `$Background.XY.Global` Background data at the global level at \code{data.frame} format.
#' - `$Background.XY.Regional` Species presence data at the regional level at \code{data.frame} format.
#' - `$Scenarios` A \code{list} containing future scenarios in \code{\link[terra:rast]{PackedSpatRaster}} format.
#' - `$Selected.Variables.Global` A \code{character} vector specifying the names of the selected covariates at the global scale.
#' - `$IndVar.Global.Selected` Selected independent variables at the global level in \code{\link[terra:rast]{PackedSpatRaster}} format.
#' - `$Selected.Variables.Regional` A \code{character} vector specifying the names of the selected covariates at the regional scale.
#' - `$IndVar.Regional.Selected` Selected independent variables at the regional level in \code{\link[terra:rast]{PackedSpatRaster}} format.
#' - `$IndVar.Global.Selected.reg` Selected variables at the global level for regional projections in \code{\link[terra:rast]{PackedSpatRaster}} format.
#' - `$Summary` Summary information about the covariates selection procedure.
#'
#'
#' @details
#' This function selects covariates for species distribution modeling with \emph{covsel} R package by combining (Step A) a collinearity-filtering algorithm and (Step B) three model-specific embedded regularization techniques, including GLM with elastic net regularization, GAM with null-space penalization, and guided regularized RF. More details can be found in (\emph{covsel} R package \doi{https://doi.org/10.1016/j.ecoinf.2023.102080} #@@@JMB dentro de doi, doi del artículo covsel)
#' If `save.output=TRUE`, selected variables at both global and regional level, are stored out of R in the \emph{Results/} folder created in the current working directory:
#' - the \emph{Results/Global/Values/} folder, containing the selected 'non-colinear' covariates at the global scale, named with the species name and \code{.variables.csv}.
#' - the \emph{Results/Regional/Values/} folder, containing the selected 'non-colinear' covariates at the regional scale, named with the species name and \code{.variables.csv}.
#'
#'
#' @seealso \code{\link{NSDM.InputData}}, \code{\link{NSDM.FormattingData}}
#'
#'
#' @examples
#' # Load required packages #@@@JMB en el ejemplo hay que poner también el NSDM.FormatingData() para tener myFormatedData? #@@@ creo que si. Ver cómo hacen otros
#' library(terra)
#' library(ecospat) #@@@JMB esto fuera cuando dependencias listas
#' library(covsel)
#'
#' # Select covariates
#' mySelectedCovs <- NSDM.SelectCovariates(myFormattedData)
#'
#' @import covsel
#'
#' @export
NSDM.SelectCovariates <- function(nsdm_finput,
				maxncov.Global="nocorr", #@@@#TG si usuario no pone ninguna usa todas las no correlacionadas #@@@JMB como sugerencia he cambiado el nombre pq "all" podría confundirse con todas las variables. Pendiente comprobar qué pasa si pone más de las correlacionadas?
				maxncov.Regional="nocorr",
				corcut=0.7,
				algorithms=c('glm','gam','rf'),
				ClimaticVariablesBands=NULL,
				save.output=TRUE) {

  if(!inherits(nsdm_finput, "nsdm.finput")){
      stop("nsdm_finput must be an object of nsdm.finput class. Consider running NSDM.FormatingData() function.")
  }

  algorithms <- tolower(algorithms)
  if(any(!algorithms %in% c('glm','gam','rf'))) {
    stop("Please select a valid algorithm (\"glm\", \"gam\", or \"rf\").")
  }

  sabina<-nsdm_finput[!names(nsdm_finput) %in% "Summary"]
  sabina$args <- list()
  sabina$args$maxncov.Global <- maxncov.Global
  sabina$args$maxncov.Regional <- maxncov.Regional
  sabina$args$corcut <- corcut
  sabina$args$algorithms <- algorithms

  SpeciesName <- nsdm_finput$Species.Name

  # Unwrap objects
  # Global independent variables (environmental layers)
  IndVar.Global <- terra::unwrap(nsdm_finput$IndVar.Global)
  # Regional independent variables (environmental layers)
  IndVar.Regional <- terra::unwrap(nsdm_finput$IndVar.Regional)

  # GLOBAL SCALE
  # Select the best subset of independent variables for each species using covsel package
  myResp.xy.Global <- rbind(nsdm_finput$SpeciesData.XY.Global, nsdm_finput$Background.XY.Global)
  names(myResp.xy.Global)<-c("x","y")
  row.names(myResp.xy.Global)<-c(1:nrow(myResp.xy.Global))
  myResp.Global <- as.vector(c(rep(1,nrow(nsdm_finput$SpeciesData.XY.Global)),rep(0,nrow(nsdm_finput$Background.XY.Global))))
  myExpl.covsel.Global <- terra::extract(IndVar.Global, myResp.xy.Global, as.df=TRUE)[, -1]

  # Variable selection process
  # Collinearity filtering
  Covdata.filter.Global<-covsel::covsel.filteralgo(covdata=myExpl.covsel.Global, pa=myResp.Global, corcut=corcut)

  # Embedding selected variables
  if(maxncov.Global=="nocorr") {
    maxncov.Global <- ncol(Covdata.filter.Global)
    Selected.Variables.Global <- names(Covdata.filter.Global)
  } else{
    Covdata.embed.Global<-covsel::covsel.embed(covdata=Covdata.filter.Global,
                                               pa=myResp.Global,
                                               algorithms=algorithms,
                                               maxncov=maxncov.Global,
                                               nthreads=detectCores()/2)
    Selected.Variables.Global <- labels(Covdata.embed.Global$covdata)[[2]]
  }


  # Save selected variables for each species
  if(save.output){
    #suffix <- 0
    #file_path <- paste0("Results/Global/Values/", SpeciesName, ".variables.csv")
    #old_path <- paste0("Results/Global/Values/", SpeciesName, ".variables.csv")
    # while (file.exists(file_path)) {
    #   suffix <- suffix + 1
    #   file_path <- gsub("\\.csv$", paste0("_", suffix, ".csv"), old_path)
   #}

    write.csv(Selected.Variables.Global, paste0("Results/Global/Values/", SpeciesName, ".variables.csv"))
    #message(paste("Selected variables at global level saved in:",file_path))
    }

  IndVar.Global.Selected <- IndVar.Global[[Selected.Variables.Global]]

  # Summary
  summary <- data.frame(Values = c(SpeciesName,
				nlyr(IndVar.Global),
				length(Selected.Variables.Global)))

  rownames(summary) <- c("Species name",
			"Original number of variables at global scale",
			"Final number of selected variables at global scale")

  # REGIONAL SCALE
  # Subset the global independent variables for regional projections
  IndVar.Global.Selected.reg <- IndVar.Regional[[Selected.Variables.Global]]

  # Exclude climatic bands specified by the user.
  Number.bands <- nlyr(IndVar.Regional)
  if(!is.null(ClimaticVariablesBands) && length(ClimaticVariablesBands) > 0) {
    # Non eliminated variables
    Bands.climatic <- setdiff(1:Number.bands, ClimaticVariablesBands)
    IndVar.Regional <- IndVar.Regional[[Bands.climatic]]
  } else {
    # If no bands are specified for exclusion, use all bands
    IndVar.Regional <- IndVar.Regional
  }

  # Select the best subset of independent variables for each species using covsel package
  myResp.xy.Regional <- rbind(nsdm_finput$SpeciesData.XY.Regional, nsdm_finput$Background.XY.Regional)
  row.names(myResp.xy.Regional) <- c(1:nrow(myResp.xy.Regional))
  myResp.Regional <- as.vector(c(rep(1, nrow(nsdm_finput$SpeciesData.XY.Regional)), rep(0, nrow(nsdm_finput$Background.XY.Regional))))
  myExpl.covsel.Regional <- terra::extract(IndVar.Regional, myResp.xy.Regional, rm.na=TRUE, df=TRUE)[, -1]

  Covdata.filter.Regional <- covsel::covsel.filteralgo(covdata = myExpl.covsel.Regional, pa = myResp.Regional, corcut = corcut)

  # Embedding selected variables
  if(maxncov.Regional=="nocorr") {
    maxncov.Regional <- ncol(Covdata.filter.Regional)
    Selected.Variables.Regional <-names(Covdata.filter.Regional)
  }else{
    Covdata.embed.Regional <- covsel::covsel.embed(covdata = Covdata.filter.Regional,
                                                   pa = myResp.Regional,
                                                   algorithms = algorithms,
                                                   maxncov = maxncov.Regional,
                                                   nthreads = detectCores() / 2)
    Selected.Variables.Regional <- labels(Covdata.embed.Regional$covdata)[[2]]
  }
  # Save selected variables for each species
  if(save.output){
    #suffix <- 0
    #file_path <- paste0("Results/Regional/Values/", SpeciesName, ".variables.csv")
    #old_path<-paste0("Results/Regional/Values/", SpeciesName, ".variables.csv")
    # while (file.exists(file_path)) {
    #   suffix <- suffix + 1
    #   file_path <- gsub("\\.csv$", paste0("_", suffix, ".csv"), old_path)
    # }

  write.csv(Selected.Variables.Regional, paste0("Results/Regional/Values/", SpeciesName, ".variables.csv"))
  #message(paste("Selected variables at regional level saved in:",file_path))
  }

  # Subset the regional independent variables for regional projections
  IndVar.Regional.Selected <- IndVar.Regional[[Selected.Variables.Regional]]

  # Summary
  summary_regional <- data.frame(Values = c(Number.bands,
				length(Selected.Variables.Regional)))

  rownames(summary_regional) <- c("Original number of variables at regiona scale",
			"Final number of of variables at regional scale")

  summary <- rbind(summary, summary_regional)

  # Wrap objects
  IndVar.Global.Selected <- terra::wrap(IndVar.Global.Selected)
  IndVar.Regional.Selected <- terra::wrap(IndVar.Regional.Selected)
  IndVar.Global.Selected.reg <- terra::wrap(IndVar.Global.Selected.reg)

  sabina$Selected.Variables.Global <- Selected.Variables.Global
  sabina$IndVar.Global.Selected <- IndVar.Global.Selected
  sabina$Selected.Variables.Regional <- Selected.Variables.Regional
  sabina$IndVar.Regional.Selected <- IndVar.Regional.Selected
  sabina$IndVar.Global.Selected.reg <- IndVar.Global.Selected.reg
  sabina$Summary<-summary
  # Make the output lighter
  sabina <- sabina[!names(sabina) %in% c("IndVar.Regional", "IndVar.Global")]

  attr(sabina, "class") <- "nsdm.vinput"

  # save.out messages
  if(save.output){
    message("Results saved in the following locations:")
    message(paste0(
       " - Selected variables at global level: Results/Global/Values/", SpeciesName, ".variables.csv\n",
       " - Selected variables at regional level: Results/Regional/Values/", SpeciesName, ".variables.csv\n"
      ))
  }

  return(sabina)

}
