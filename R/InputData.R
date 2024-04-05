#' @export
NSH.SDM.InputData <- function(SpeciesName,
				spp.data.global,		#@@@JMB debe ser data.frame con x,y
				spp.data.regional,	#@@@JMB debe ser data.frame con x,y
				expl.var.global,		#@@@JMB debe ser SpatRasrer
				expl.var.regional,	#@@@JMB debe ser SpatRasrer
				new.env = NULL,		#@@@JMB cuando no null, debe ser SpatRasrer o list de SpatRaster
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
  
  if(!is.null(new.env)) {
    if(!(class(new.env) == "SpatRaster" || is.list(new.env) && all(sapply(new.env, function(x) class(x) == "SpatRaster")))) {
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

  # Match variables
  match_vars <- sapply(new.env, function(file) { #@@@JMB check que las variables seleccionadas están en todos los escenarios
    all(names(expl.var.regional) %in% names(file))
  })
  if(!all(match_vars)) {
   message("Not all scenarios have the same variables.")
  }

  # Rename SpatRaster
  if(!is.null(new.env)) {
    names(new.env) <- path_ext_remove(path_file(names(new.env)))
  }

  if(!is.null(new.env) && !is.null(new.env.names)) {
    if(is.list(new.env) && length(new.env) != length(new.env.names)) {
      stop("The number of provided new.env.names does not match the number of elements in new.env.")
    } else {
      names(new.env) <- new.env.names
    }
  }

  #
  nshsdm_data <- list(Species.Name = SpeciesName,
    SpeciesData.XY.Global = spp.data.global,
    SpeciesData.XY.Regional = spp.data.regional,
    IndVar.Global = expl.var.global,
    IndVar.Regional = expl.var.regional,
    Scenarios = new.env,
    Background.Global = Background.Global,
    Background.Regional = Background.Regional
  )
  
  attr(nshsdm_data, "class") <- "nshsdm.data"

  return(nshsdm_data)
}



