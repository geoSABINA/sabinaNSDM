#' @name NSDM.SelectCovariates
#'
#' @title Select covariates for spatially-nested hierarchical species distribution modeling (NSDM) analysis.
#'
#' @description This function selects the best 'non-colinear' environmental covariates for \bold{NSDM} based on specified criteria and algorithms.
#'
#'
#' @param nsdm_finput An object of class \code{nsdm.finput} generated using the \code{\link{NSDM.FormatingData}} function.
#' @param maxncov.Global (\emph{optional, default} \code{'nocorr'}) \cr
#' Maximum \code{numeric} value indicating the maximum number of covariates to select at the global scale. If `"nocorr"`, selects all non-correlated covariates.
#' @param maxncov.Regional (\emph{optional, default} \code{'nocorr'}) \cr
#' Maximum \code{numeric} value indicating the maximum number of covariates to select at the regional scale. If `"nocorr"`, selects all non-correlated covariates.
#' @param corcut (\emph{optional, default} \code{0.7}) \cr
#' A \code{numeric} value for the correlation coefficient threshold used for identifying collinearity.
#' @param algorithms (\emph{optional, default} \code{'c("glm", "gam", "rf")'}) \cr
#' Algorithms to use for ranking the covariates. Options are \code{'glm'}, \code{'gam'}, and/or \code{'rf'}.
#' @param ClimaticVariablesBands (\emph{optional, default} \code{NULL}) \cr
#' Indices of the regional environmental covariate bands to exclude from the selection. If \code{NULL} (the default), all regional-level covariates are considered. The excluded covariates typically include climatic covariates already included in global-level analyses. A list of band numbers should be included. For example, use ClimaticVariablesBands = c(2, 3, 4) to exclude bands 2, 3, and 4 in the regional analysis.
#' A \code{logical} value defining whether the outputs should be saved at local.
#'
#'
#' @return An object of class \code{nsdm.vinput} containing selected covariates for \bold{NSDM}:
#' - `$SpeciesName` Name of the species.
#' - `$args` A \code{list} containing the arguments used during covariates selection procedure, including: `maxncov.Global`, `maxncov.Regional`, `corcut` and `algorithms`.
#' - `$SpeciesData.XY.Global` Species occurrences at the global level at \code{data.frame} format after applying spatial thinning.
#' - `$SpeciesData.XY.Regional` Species occurrences at the regional level at \code{data.frame} format after applying spatial thinning.
#' - `$Background.XY.Global` Background data at the global level at \code{data.frame} format.
#' - `$Background.XY.Regional` Background data at the regional level at \code{data.frame} format.
#' - `$Scenarios` A \code{list} containing future scenarios in \code{\link[terra:rast]{PackedSpatRaster}} format.
#' - `$Selected.Variables.Global` A \code{character} vector specifying the names of the selected covariates at the global scale.
#' - `$IndVar.Global.Selected` Selected covariates at the global level in \code{\link[terra:rast]{PackedSpatRaster}} format.
#' - `$Selected.Variables.Regional` A \code{character} vector specifying the names of the selected covariates at the regional scale.
#' - `$IndVar.Regional.Selected` Selected covariates at the regional level in \code{\link[terra:rast]{PackedSpatRaster}} format.
#' - `$IndVar.Global.Selected.reg` Selected covariates at the global level for regional projections in \code{\link[terra:rast]{PackedSpatRaster}} format.
#' - `$Summary` Summary information about the covariates selection procedure.
#'
#'
#' @details
#' This function selects covariates for species distribution modeling with \emph{covsel} R package by combining (Step A) a collinearity-filtering algorithm and (Step B) three model-specific embedded regularization techniques, including GLM with elastic net regularization, GAM with null-space penalization, and guided regularized RF. More details can be found in (\emph{covsel} R package \doi{https://doi.org/10.1016/j.ecoinf.2023.102080}
#' If `save.output=TRUE`, selected covariates at both global and regional level, are stored out of R in the \emph{Results/} folder created in the current working directory:
#' - the \emph{Results/Global/Values/} folder, containing the selected 'non-colinear' covariates at the global scale, named with the species name and \code{.variables.csv}.
#' - the \emph{Results/Regional/Values/} folder, containing the selected 'non-colinear' covariates at the regional scale, named with the species name and \code{.variables.csv}.
#'
#'
#' @seealso \code{\link{NSDM.InputData}}, \code{\link{NSDM.FormattingData}}
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
#' myInputData<-NSDM.InputData(
#'		SpeciesName = "Fagus.sylvatica",
#'		spp.data.global = Fagus.sylvatica.xy.global,
#'		spp.data.regional = Fagus.sylvatica.xy.regional,
#'		expl.var.global = expl.var.global,
#'		expl.var.regional = expl.var.regional,
#'		new.env = new_env,
#'		new.env.names = c("Scenario1"),
#'		Background.Global = NULL,
#'		Background.Regional = NULL
#' )
#'
#' # Format the input data
#' myFormatedData <- NSDM.FormatingData(myInputData,
#'					nPoints=1000)

#' # Select covariates using default parameters
#' mySelectedCovs <- NSDM.SelectCovariates(myFormattedData)
#' 
#' # Select covariates using custom parameters.
#' nsdm_selvars <- NSDM.SelectCovariates(nsdm_finput,
#'     maxncov.Global = 5, # Maximum number of global covariates
#'     maxncov.Regional = 7, # Maximum number of regional covariates
#'     corcut = 0.7, # Maximum number of regional covariates
#'     algorithms = c("glm","gam","rf"),  # Algorithms to use for selection
#'     ClimaticVariablesBands = c(2,3,5), # Bands to exclude in the analysis
#'     save.output = TRUE  # Save the output externally
#' )
#'
#' @import covsel
#'
#' @export
NSDM.SelectCovariates <- function(nsdm_finput,
				maxncov.Global="nocorr",
				maxncov.Regional="nocorr",
				corcut=0.7,
				algorithms=c('glm','gam','rf'),
				ClimaticVariablesBands=NULL,
				save.output=TRUE) {

  if(!inherits(nsdm_finput, "nsdm.finput")){
      stop("nsdm_finput must be an object of nsdm.finput class. Consider running NSDM.FormatingData() function.")
  }
  if(!is.null(ClimaticVariablesBands)) {if(!inherits(ClimaticVariablesBands, "numeric") & !inherits(ClimaticVariablesBands, "integer") ){
    stop("ClimaticVariablesBands must be either NULL or a vector indicating the regional environmental covariate bands to exclude from the selection at regional scale.")
  }}

  algorithms <- tolower(algorithms)
  if(any(!algorithms %in% c('glm','gam','rf'))) {
    stop("Please select a valid algorithm (\"glm\", \"gam\", or \"rf\")")
  }

  sabina<-nsdm_finput[!names(nsdm_finput) %in% "Summary"]
  sabina$args <- list()
  sabina$args$maxncov.Global <- maxncov.Global
  sabina$args$maxncov.Regional <- maxncov.Regional
  sabina$args$corcut <- corcut
  sabina$args$algorithms <- algorithms

  SpeciesName <- nsdm_finput$Species.Name

  # Unwrap objects
  # Global independent covariates (environmental layers)
  IndVar.Global <- terra::unwrap(nsdm_finput$IndVar.Global)
  # Regional independent covariates (environmental layers)
  IndVar.Regional <- terra::unwrap(nsdm_finput$IndVar.Regional)

  # GLOBAL SCALE
  # Select the best subset of independent covariates for each species using covsel package
  myResp.xy.Global <- rbind(nsdm_finput$SpeciesData.XY.Global, nsdm_finput$Background.XY.Global)
  names(myResp.xy.Global)<-c("x","y")
  row.names(myResp.xy.Global)<-c(1:nrow(myResp.xy.Global))
  myResp.Global <- as.vector(c(rep(1,nrow(nsdm_finput$SpeciesData.XY.Global)),rep(0,nrow(nsdm_finput$Background.XY.Global))))
  myExpl.covsel.Global <- terra::extract(IndVar.Global, myResp.xy.Global, as.df=TRUE)[, -1]

  # Covariates selection process
  # Collinearity filtering
  Covdata.filter.Global<-covsel::covsel.filteralgo(covdata=myExpl.covsel.Global, pa=myResp.Global, corcut=corcut)

  # Embedding selected covariates
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


  # Save selected covariates for each species
  if(save.output){

    write.csv(Selected.Variables.Global, paste0("Results/Global/Values/", SpeciesName, ".variables.csv"))
    }

  IndVar.Global.Selected <- IndVar.Global[[Selected.Variables.Global]]

  # Summary
  summary <- data.frame(Values = c(SpeciesName,
				nlyr(IndVar.Global),
				length(Selected.Variables.Global)))

  rownames(summary) <- c("Species name",
			"Original number of covariates at the global scale",
			"Final number of selected covariates at the global scale")

  # REGIONAL SCALE
  # Subset the global covariates for regional projections
  IndVar.Global.Selected.reg <- IndVar.Regional[[Selected.Variables.Global]]

  # Exclude climatic bands specified by the user.
  Number.bands <- nlyr(IndVar.Regional)
  if(!is.null(ClimaticVariablesBands) && length(ClimaticVariablesBands) > 0) {
    removed<-IndVar.Regional[[ClimaticVariablesBands]]
    names(removed)
    message(paste("The bands", paste(names(removed),collapse=", "), "have been excluded from the selection of regional environmental covariates."))
    # Non eliminated variables
    Bands.climatic <- setdiff(1:Number.bands, ClimaticVariablesBands)
    IndVar.Regional <- IndVar.Regional[[Bands.climatic]]
  } else {
    # If no bands are specified for exclusion, use all bands
    IndVar.Regional <- IndVar.Regional
  }

  # Select the best subset of covariates for each species using covsel package
  myResp.xy.Regional <- rbind(nsdm_finput$SpeciesData.XY.Regional, nsdm_finput$Background.XY.Regional)
  row.names(myResp.xy.Regional) <- c(1:nrow(myResp.xy.Regional))
  myResp.Regional <- as.vector(c(rep(1, nrow(nsdm_finput$SpeciesData.XY.Regional)), rep(0, nrow(nsdm_finput$Background.XY.Regional))))
  myExpl.covsel.Regional <- terra::extract(IndVar.Regional, myResp.xy.Regional, rm.na=TRUE, df=TRUE)[, -1]

  Covdata.filter.Regional <- covsel::covsel.filteralgo(covdata = myExpl.covsel.Regional, pa = myResp.Regional, corcut = corcut)

  # Embedding selected covariates
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
  # Save selected covariates for each species
  if(save.output){

  write.csv(Selected.Variables.Regional, paste0("Results/Regional/Values/", SpeciesName, ".variables.csv"))
  }

  # Subset the regional independent variables for regional projections
  IndVar.Regional.Selected <- IndVar.Regional[[Selected.Variables.Regional]]

  # Summary
  summary_regional <- data.frame(Values = c(Number.bands,
				length(Selected.Variables.Regional)))

  rownames(summary_regional) <- c("Original number of covariates at regiona scale",
			"Final number of covariates at regional scale")

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
       " - Selected covariates at the global level: Results/Global/Values/", SpeciesName, ".variables.csv\n",
       " - Selected covariates at the regional level: Results/Regional/Values/", SpeciesName, ".variables.csv\n"
      ))
  }

  return(sabina)

}
