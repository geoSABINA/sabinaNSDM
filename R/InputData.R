#' @export
NSH.SDM.InputData <- function(SpeciesName,
				spp.data.global,		#@@@JMB debe ser data.frame con x,y
				spp.data.regional,	#@@@JMB debe ser data.frame con x,y
				expl.var.global,		#@@@JMB debe ser SpatRasrer
				expl.var.regional,	#@@@JMB debe ser SpatRasrer
				new.env = NULL,		#@@@JMB cuando no null, debe ser SpatRasrer o list de SpatRaster. Pensar si es más simple una cadena de spatrasters y hacer dentro la lista
				new.env.names = NULL,	#@@@JMB optional
				Background.Global = NULL,	#@@@JMB debe ser data.frame con x,y
				Background.Regional = NULL) {	#@@@JMB debe ser data.frame con x,y


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

  if(!is.null(Background.Global) && !is.null(Background.Regional)) { #@@@JMB esto igual se puede separar en dos para que puedan poner uno u otro.
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
  match_vars <- sapply(new.env, function(file) {
    all(names(expl.var.regional) %in% names(file))
  })
  if(!all(match_vars)) {
   message("Not all scenarios have the same variables.")
  }

  # Rename new.env
  if(!is.null(new.env)) {
    if(!is.null(names(new.env))) {
      names(new.env) <- path_ext_remove(names(new.env))
    } else {
    stop("names of new.env is null. Please provide names of new scenarios in names.new.env.")
    }
  }

  if(!is.null(new.env) && !is.null(new.env.names)) {
    if(length(new.env) != length(new.env.names)) {
      stop("The number of provided new.env.names does not match the number of elements in new.env.")
    } else {
      names(new.env) <- new.env.names
    }
  }

  #
  nshsdm_data <- list()
  nshsdm_data$Species.Name <- SpeciesName
  nshsdm_data$SpeciesData.XY.Global.0 <- spp.data.global
  nshsdm_data$SpeciesData.XY.Regional.0 <- spp.data.regional
  nshsdm_data$IndVar.Global <- expl.var.global
  nshsdm_data$IndVar.Regional <- expl.var.regional
  nshsdm_data$Scenarios <- new.env
  nshsdm_data$Background.Global.0 <- Background.Global
  nshsdm_data$Background.Regional.0 <- Background.Regional
  
  attr(nshsdm_data, "class") <- "nshsdm.data"

  return(nshsdm_data)
}



