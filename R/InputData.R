#' @name NSDM.InputData
#'
#' @title Prepare input data for Nested spatial hierarchical Species Distribution Modeling (NSDM) analysis.
#'
#' @description This function gathers together all input data (\emph{species presences xy at global and regional level, explanatory variables at global and regional level, new environments, and if available, absence/pseudo-absences/background data}) needed to run \NSDM.
#'
#' @param SpeciesName A \code{character} specifying the species name for which data is being prepared.
#' @param spp.data.global A \code{data.frame} with two columns, 'x' and 'y', representing the species presence data at global level.
#' @param spp.data.regional A \code{data.frame} with two columns, 'x' and 'y', representing the species presence data at regional level.
#' @param expl.var.global An object of \code{\link[terra:rast]{SpatRaster}} class representing the explanatory variables at global level.
#' @param expl.var.regional An object of \code{\link[terra:rast]{SpatRaster}} class representing the explanatory variables at regional level.
#' @param new.env (\emph{optional, default} \code{NULL}) \cr 
#' An object of \code{\link[terra:rast]{SpatRaster}} class or a \code{list} of \code{\link[terra:rast]{SpatRaster}} objects representing new environmental scenarios.
#' @param new.env.names (\emph{optional, default} \code{NULL}) \cr 
#' A character \code{vector} specifying names for the new environmental scenarios.
#' @param Background.Global (\emph{optional, default} \code{NULL}) \cr
#' An \code{data.frame} with two columns, 'x' and 'y', representing background points at global level.
#' @param Background.Regional (\emph{optional, default} \code{NULL}) \cr
#' An \code{data.frame} with two columns, 'x' and 'y', representing background points at regional level.
#'
#' @return An object of class "nsdm.input" containing organized input data for NSDM.
#'
#' @details
#' - The `SpeciesName` parameter specifies the name of the species for which data is being prepared.
#' - The `spp.data.global` and `spp.data.regional` parameters should be \code{data.frame} with columns 'x' and 'y' representing species presence data.
#' - The `expl.var.global` and `expl.var.regional` parameters should be \code{\link[terra:rast]{SpatRaster}} objects representing explanatory variables.  #@@@JMB Aquí explicar que tienen que ser un .tif con las variables.  
#' - The `new.env` parameter can be either a \code{\link[terra:rast]{SpatRaster}} object or a \code{list} of \code{\link[terra:rast]{SpatRaster}} objects representing new environmental scenarios.
#' - The `Background.Global` and `Background.Regional` parameters are optional and represent background points at the global and regional levels, respectively.
#'
#' @examples
#' library(terra)
#' 
#' # Load species occurrences
#' data(Fagus.sylvatica.xy.global, package = "sabina")
#' data(Fagus.sylvatica.xy.regional, package = "sabina")
#' 
#' # Load explanatory variables
#' data(expl.var.global, package = "sabina")
#' data(expl.var.regional, package = "sabina")
#' expl.var.global<-terra::unwrap(expl.var.global)
#' expl.var.regional<-terra::unwrap(expl.var.regional)
#' 
#' # Load new escenarios
#' data(new.env, package = "sabina")
#' new.env<-terra::unwrap(new.env)
#'
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
#'
#' @export
NSDM.InputData <- function(SpeciesName,
				spp.data.global,
				spp.data.regional,
				expl.var.global,
				expl.var.regional,
				new.env = NULL,
				new.env.names=NULL,
				Background.Global = NULL,
				Background.Regional = NULL) {


  if(!(is.data.frame(spp.data.global) && 
        ncol(spp.data.global) == 2 && 
        all(c('x', 'y') %in% colnames(spp.data.global))) ||
      !(is.data.frame(spp.data.regional) && 
        ncol(spp.data.regional) == 2 && 
        all(c('x', 'y') %in% colnames(spp.data.regional)))) {
    stop("spp.data.global and spp.data.regional must be data.frames with 'x' and 'y' columns")
  }

  if(!inherits(expl.var.global, "SpatRaster") || !inherits(expl.var.regional, "SpatRaster")) {
    stop("expl.var.global and expl.var.regional must be SpatRaster objects.")
  }
  
  if(inherits(new.env, "SpatRaster")) {
    new.env <- list(new.env)
  }

  if(!is.null(new.env)) {
    if(!(is.list(new.env) && all(sapply(new.env, function(x) class(x) == "SpatRaster")))) {
      stop("new.env must be either a SpatRaster object or a list of SpatRaster objects")
    }
  }

  if(!is.null(Background.Global) && !is.null(Background.Regional)) {
    if(!(is.data.frame(Background.Regional) && 
          ncol(Background.Global) == 2 && 
          all(c('x', 'y') %in% colnames(Background.Global))) ||
        !(is.data.frame(Background.Regional) && 
          ncol(Background.Regional) == 2 && 
          all(c('x', 'y') %in% colnames(Background.Regional)))) {
      stop("Background.Global and Background.Regional must be data.frames with 'x' and 'y' columns")
    }
  }

  # Match variables?
  # Regionsal vs new scenarios
  match_vars <- sapply(new.env, function(file) {
    all(names(expl.var.regional) %in% names(file))
  })
  if(!all(match_vars)) {
   message("Not all scenarios have the same variables.")
  }
  #Global vars in regional?
  match_vars2 <- sapply(names(expl.var.global), function(var_name) {
    var_name %in% names(expl.var.regional)
  })
  if(!all(match_vars2)) {
    stop("All variables present in expl.var.global must also be present in expl.var.regional.")
  }

  # Name and Rename new.env
  if(!is.null(new.env)) {
    if(is.null(names(new.env)) && is.null(new.env.names)) {
      source_names <- lapply(new.env, function(x) {
        s <- sources(x)
        s |> path_file() |> path_ext_remove()
      })
      names(new.env) <- source_names
      valid_names <- sapply(names(new.env), function(x) sum(nchar(x)))
      if(!all(valid_names > 0)) {
        stop("Please provide names of new scenarios in names.new.env")
      }
    } else if(!is.null(new.env.names)) {
      if(length(new.env) != length(new.env.names)) {
        stop("The number of provided new.env.names does not match the number of elements in new.env.")
      } else {
        names(new.env) <- new.env.names
      }
    } else if(is.null(names(new.env)) && is.null(new.env.names)) {
    stop("Names of new.env is NULL. Please provide names of new scenarios in names.new.env.")
    }
  }

  # Wrap objects if necessary
  expl.var.global <- terra::wrap(expl.var.global)
  expl.var.regional <- terra::wrap(expl.var.regional)

  if(!is.null(new.env)) {
    new.env <- lapply(new.env, terra::wrap)
  }
  
  # Summary
  summary <- data.frame(Values = c(SpeciesName,
				nrow(spp.data.global), 
				ifelse(is.null(Background.Global), "NULL", nrow(Background.Global)),
				nrow(spp.data.regional),
				ifelse(is.null(Background.Regional), "NULL", nrow(Background.Regional)),
				length(new.env)))
  
  rownames(summary) <- c("Species name",
                         "Original number of species presences at global level", 
                         "Original number of background points at global level",
                         "Original number of species presences at regional level", 
                         "Original number of background points at regional level",
                         "Number of new scenarios")

  #
  sabina <- list(
    Species.Name = SpeciesName,
    SpeciesData.XY.Global.0 = spp.data.global,
    SpeciesData.XY.Regional.0 = spp.data.regional,
    IndVar.Global = expl.var.global,
    IndVar.Regional = expl.var.regional,
    Scenarios = new.env,
    Background.Global.0 = Background.Global,
    Background.Regional.0 = Background.Regional,
    Summary = summary
  )
  
  attr(sabina, "class") <- "nsdm.input"

  return(sabina)
}



