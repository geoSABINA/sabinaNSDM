#' @name NSDM.InputData
#'
#' @title Prepare input data for spatially-nested hierarchical species distribution modeling (NSDM) analysis.
#'
#' @description This function gathers together all input data (\emph{species occurrences xy at global and regional level, enviromental covariates at global and regional level, new environments, and if available, absence/pseudo-absences/background data}) needed to run \bold{NSDM}.
#'
#' @param SpeciesName A \code{character} specifying the species name for which data is being prepared.
#' @param spp.data.global A \code{data.frame} with two columns, 'x' and 'y', representing the species presence data at global level.
#' @param spp.data.regional A \code{data.frame} with two columns, 'x' and 'y', representing the species presence data at regional level.
#' @param expl.var.global An object of \code{\link[terra:rast]{SpatRaster}} class representing the environmental covariates (usually climatic) at global level.
#' @param expl.var.regional An object of \code{\link[terra:rast]{SpatRaster}} class representing the environmental covariates at regional level.
#' @param new.env (\emph{optional, default} \code{NULL}) \cr
#' An object of \code{\link[terra:rast]{SpatRaster}} class or a \code{list} of \code{\link[terra:rast]{SpatRaster}} objects representing the environmental covariates in new scenarios for projecting the models to different spatial or temporal scales.
#' @param new.env.names (\emph{optional, default} \code{NULL}) \cr
#' A character \code{vector} specifying names for the new environmental scenarios.
#' @param Background.Global (\emph{optional, default} \code{NULL}) \cr
#' An optional \code{data.frame} with two columns, 'x' and 'y', representing background points at global level.
#' @param Background.Regional (\emph{optional, default} \code{NULL}) \cr
#' An optional \code{data.frame} with two columns, 'x' and 'y', representing background points at regional level.
#' @param Absences.Global (\emph{optional, default} \code{NULL}) \cr
#' An optional \code{data.frame} with two columns, 'x' and 'y', representing absence points at global level.
#' @param Absences.Regional (\emph{optional, default} \code{NULL}) \cr
#' An optional \code{data.frame} with two columns, 'x' and 'y', representing absence points at regional level.
#'
#' @return An object of class \code{nsdm.input} containing organized input data for \bold{NSDM}.
#'
#'
#' @details
#' - The `expl.var.global` and `expl.var.regional` parameters should be \code{\link[terra:rast]{SpatRaster}} objects representing environmental covariates used to train the global-scale and regional-scale models, respectively. Each band corresponds to a different covariate. The regional-scale object must include all the covariates included in the global-scale object (usually climatic), and it can additionally include other covariates only available at this level.
#' - The `new.env` parameter can be either a \code{\link[terra:rast]{SpatRaster}} object or a \code{list} of \code{\link[terra:rast]{SpatRaster}} object representing the environmental covariates in new spatial or temporal scenarios. All new scenarios must have the same covariates (same band names) than `expl.var.regional`.
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
#'				Background.Regional = NULL,
#'				Absences.Global = NULL,
#'				Absences.Regional = NULL)
#'
#' summary(myInputData)
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
                           Background.Regional = NULL,
                           Absences.Global = NULL,
                           Absences.Regional = NULL) {

  if(!(is.coord.df(spp.data.global)) ||
       !(is.coord.df(spp.data.regional))) {
    stop(paste0("spp.data.global and spp.data.regional for ", SpeciesName, " must be data.frames with 'x' and 'y' columns."))
  }

  if(!inherits(expl.var.global, "SpatRaster") || !inherits(expl.var.regional, "SpatRaster")) {
    stop("expl.var.global and expl.var.regional must be SpatRaster objects.")
  }

  if(!is.null(new.env)) {
    if(inherits(new.env, "SpatRaster")) {
      new.env <- list(new.env)
    }
    if(!(is.list(new.env) && all(sapply(new.env, function(x) class(x) == "SpatRaster")))) {
      stop("new.env must be either a SpatRaster object or a list of SpatRaster objects.")
    }
  }

  if(!is.null(Absences.Global) && !is.null(Background.Global)) {
    stop("Background and Absences cannot both be incorporated at the global level.")
  }

  if(!is.null(Absences.Regional) && !is.null(Background.Regional)) {
    stop("Background and Absences cannot both be incorporated at the regional level.")
  }

  if(!is.null(Background.Global) && !is.null(Background.Regional)) {
    if(!(is.coord.df(Background.Global)) ||
       !(is.coord.df(Background.Regional))) {
      stop("Background.Global and Background.Regional must be data.frames with 'x' and 'y' columns.")
    }
  }

  if(!is.null(Absences.Global) && !is.coord.df(Absences.Global)) {
    stop("Absences.Global must be a data.frame with 'x' and 'y' columns.")
  }

  if(!is.null(Absences.Regional) && !is.coord.df(Absences.Regional)) {
    stop("Absences.Regional must be a data.frame with 'x' and 'y' columns.")
  }

  # Match variables?
  # Global vars in regional?
  match_vars2 <- names(expl.var.global) %in% names(expl.var.regional)
  if(!all(match_vars2)) {
    stop("All variables present in expl.var.global must also be present in expl.var.regional.")
  }
  # new.env and regional
  if(!is.null(new.env)) {
    match_vars <- sapply(new.env, function(x) {
      all(names(expl.var.regional) %in% names(x))
    })
    if(!all(match_vars)) {
      stop("Not all new scenarios have the same environmental covariates as expl.var.regional.")
    }

    # Name and Rename new.env scenarios
    if(is.null(names(new.env)) && is.null(new.env.names)) {
      source_names <- lapply(new.env, function(x) {
        s <- terra::sources(x)
        s |> fs::path_file() |> fs::path_ext_remove()
      })
      names(new.env) <- source_names
      valid_names <- sapply(names(new.env),
                            function(x){
                                !is.null(x) && all(nchar(x) > 0)})
      if(!all(valid_names)) {
        stop("Names of new.env are NULL or empty strings. Please provide names of new scenarios in names.new.env.")
      }
    } else if(!is.null(new.env.names)) {
      if(length(new.env) != length(new.env.names)) {
        stop("The number of provided new.env.names does not match the number of elements in new.env.")
      } else {
        names(new.env) <- new.env.names
      }
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
                                   ifelse(is.null(Background.Global),
                                          "NA", nrow(Background.Global)),
                                   ifelse(is.null(Absences.Global), 
                                          "NA", nrow(Absences.Global)),
                                   nrow(spp.data.regional),
                                   ifelse(is.null(Background.Regional),
                                          "NA", nrow(Background.Regional)),
                                   ifelse(is.null(Absences.Regional),
                                          "NA", nrow(Absences.Regional)),
                                   length(new.env)))

  rownames(summary) <- c("Species name",
                         "Original number of species occurrences at global level",
			 "Number of background points at global level",
			 "Number of true absence points at global level",
                         "Original number of species occurrences at regional level",
			 "Number of background points at regional level",
			 "Number of true absence points at regional level",
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
    Absences.Global = Absences.Global,
    Absences.Regional = Absences.Regional,
    Summary = summary
  )

  attr(sabina, "class") <- "nsdm.input"

  return(sabina)
}
