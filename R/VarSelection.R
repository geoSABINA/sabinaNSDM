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
#' @param save.output (\emph{optional, default} \code{TRUE}) \cr
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
#' myInputData<-NSDM.InputData(SpeciesName = "Fagus.sylvatica",
#'				spp.data.global = Fagus.sylvatica.xy.global,
#'				spp.data.regional = Fagus.sylvatica.xy.regional,
#'				expl.var.global = expl.var.global,
#'				expl.var.regional = expl.var.regional,
#'				new.env = new_env,
#'				new.env.names = c("Scenario1"),
#'				Background.Global = NULL,
#'				Background.Regional = NULL)
#'
#' # Format the input data
#' myFormattedData <- NSDM.FormatingData(myInputData,
#'					nPoints=1000)
#'
#' # Select covariates using default parameters
#' mySelectedCovs <- NSDM.SelectCovariates(myFormattedData)
#'
#' summary(mySelectedCovs)
#' 
#' ## Select covariates using custom parameters.
#' # mySelectedCovs <- NSDM.SelectCovariates(nsdm_finput,
#' #					maxncov.Global = 5, 	# Maximum number of global covariates
#' #					maxncov.Regional = 7, 	# Maximum number of regional covariates
#' #					corcut = 0.7, 		# Value of the correlation coefficient threshold used for identifying collinearity
#' #					algorithms = c("glm","gam","rf"),  # Algorithms to use for selection
#' #					ClimaticVariablesBands = c(2,3,5), # Bands to exclude in the analysis
#' #					save.output = TRUE)  	# Save the output externally
#'
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
  if(!is.null(ClimaticVariablesBands)){
      if(!(is.numeric(ClimaticVariablesBands)) & !(is.integer(ClimaticVariablesBands)) ){
    stop("ClimaticVariablesBands must be either NULL or a vector indicating the regional environmental covariate bands to exclude from the selection at regional scale.")
  }}

  algorithms <- tolower(algorithms)
  if(any(!algorithms %in% c('glm','gam','rf'))) {
    stop("Please select a valid algorithm (\"glm\", \"gam\", or \"rf\")")
  }

  sabina<-nsdm_finput[!(names(nsdm_finput) %in% "Summary")]
  sabina$args <- list()
  sabina$args$maxncov.Global <- maxncov.Global
  sabina$args$maxncov.Regional <- maxncov.Regional
  sabina$args$corcut <- corcut
  sabina$args$algorithms <- algorithms

  SpeciesName <- nsdm_finput$Species.Name

  selected_global <- select_cov(nsdm_finput, scale = "Global",
                                ClimaticVariablesBands,
                                maxncov.Global, corcut, algorithms,
                                save.output)

  selected_regional <- select_cov(nsdm_finput, scale = "Regional",
                                  ClimaticVariablesBands,
                                  maxncov.Regional, corcut, algorithms,
                                  save.output,
                                  selected_global$Selected.Variables)

  main_summary <- rbind(selected_global$Summary,
                        selected_regional$Summary)

  sabina$Selected.Variables.Global <- selected_global$Selected.Variables
  sabina$IndVar.Global.Selected <- selected_global$IndVar.Selected
  sabina$Selected.Variables.Regional <- selected_regional$Selected.Variables
  sabina$IndVar.Regional.Selected <- selected_regional$IndVar.Selected
  sabina$IndVar.Global.Selected.reg <- selected_regional$IndVar.Global.Selected.reg
  sabina$Summary <- main_summary
  # Make the output lighter
  sabina <- sabina[!names(sabina) %in% c("IndVar.Regional", "IndVar.Global")]

  attr(sabina, "class") <- "nsdm.vinput"
  class(sabina) <- c("nsdm.vinput", "nsdm.input")

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

select_cov <- function(nsdm_finput, scale, ClimaticVariablesBands,
                       maxncov, corcut, algorithms, save.output,
                       Selected.Variables.Global = NULL){

  if(!(scale %in% c("Global", "Regional"))){
    stop("scale must be either Global or Regional")
  }

  if(scale == "Regional" && is.null(Selected.Variables.Global)){
      stop("Regional scale variable selecion requires running select_cov on Global scale ",
           "and provide Selected.Variables.Global as argument.")
  }

  IndVar <- terra::unwrap(nsdm_finput[[paste0("IndVar.", scale)]])

  if(scale == "Regional"){
    # Subset the global covariates for regional projections
    IndVar.Global.Selected.reg <- IndVar[[Selected.Variables.Global]]

    # Exclude climatic bands specified by the user.
    Number.bands <- terra::nlyr(IndVar)
    if(!is.null(ClimaticVariablesBands) && length(ClimaticVariablesBands) > 0) {
      removed<-IndVar[[ClimaticVariablesBands]]
      message(paste("The bands", paste(names(removed),collapse=", "),
                    "have been excluded from the selection of regional environmental covariates."))
      # Non eliminated variables
      Bands.climatic <- setdiff(1:Number.bands, ClimaticVariablesBands)
      IndVar <- IndVar[[Bands.climatic]]
    }
  }

  species_name <- paste0("SpeciesData.XY.", scale)
  background_name <- paste0("Background.XY.", scale)
  absences_name <- paste0("Absences.XY.", scale)


  # Select the best subset of covariates for each species using covsel package
  if(!is.null(nsdm_finput[[background_name]])) {
    myResp.xy <- rbind(nsdm_finput[[species_name]], nsdm_finput[[background_name]])
  } else {
    myResp.xy <- rbind(nsdm_finput[[species_name]], nsdm_finput[[absences_name]])
  }
  row.names(myResp.xy)<-c(1:nrow(myResp.xy))
  if(!is.null(nsdm_finput[[background_name]])) {
    myResp <- as.vector(c(rep(1,nrow(nsdm_finput[[species_name]])),rep(0,nrow(nsdm_finput[[background_name]]))))
  } else {
    myResp <- as.vector(c(rep(1, nrow(nsdm_finput[[species_name]])), rep(0, nrow(nsdm_finput[[absences_name]]))))
  }
  myExpl.covsel <- terra::extract(IndVar, myResp.xy, as.df=TRUE)[, -1]

  Covdata.filter <- covsel::covsel.filteralgo(covdata = myExpl.covsel,
                                              pa = myResp, corcut = corcut)

  # Embedding selected covariates
  if(maxncov == "nocorr") {
    maxncov <- ncol(Covdata.filter)
    Selected.Variables <-names(Covdata.filter)
  } else {
    Covdata.embed <- covsel::covsel.embed(covdata = Covdata.filter,
                                          pa = myResp,
                                          algorithms = algorithms,
                                          maxncov = maxncov,
                                          nthreads = detectCores() / 2)
    Selected.Variables <- labels(Covdata.embed$covdata)[[2]]
  }

  # Save selected covariates for each species
  if(save.output){
    write.csv(Selected.Variables, paste0("Results/", scale, "/Values/", SpeciesName, ".variables.csv"))
  }

  # Subset the regional independent variables for regional projections
  IndVar.Selected <- IndVar[[Selected.Variables]]

  if(scale == "Global"){
    summary <- data.frame(Values = c(SpeciesName,
                                     terra::nlyr(IndVar),
                                     length(Selected.Variables)))

    rownames(summary) <- c("Species name",
                           "Original number of covariates at the global scale",
                           "Final number of selected covariates at the global scale")
  }else if(scale == "Regional"){
    summary <- data.frame(Values = c(Number.bands,
                                     length(Selected.Variables)))
    rownames(summary) <- c("Original number of covariates at regional scale",
                           "Final number of covariates at regional scale")

  }

  select_list <- list(IndVar.Selected = terra::wrap(IndVar.Selected),
                      Selected.Variables = Selected.Variables,
                      Summary = summary)
  if(scale == "Regional"){
      select_list$IndVar.Global.Selected.reg <- terra::wrap(IndVar.Global.Selected.reg)
  }

  return(select_list)

}
