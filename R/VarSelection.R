#' @export
NSDM.SelectVariables <- function(nsdm_finput,
				maxncov.Global="nocorr", #@@@#TG si usuario no pone ninguna usa todas las no correlacionadas #@@@JMB como sugerencia he cambiado el nombre pq "all" podría confundirse con todas las variables. Pendiente comprobar qué pasa si pone más de las correlacionadas?
				maxncov.Regional="nocorr",
				corcut=0.7,
				algorithms=c('glm','gam','rf'),
				ClimaticVariablesBands=NULL,
				save.output=TRUE) { 

  #nsdm_name <- as.list(match.call())$nsdm_finput
  if(!inherits(nsdm_finput, "nsdm.input")){
      stop("nsdm_finput must be an object of nsdm.input class. Consider running NSDM.FormatingData() function.")
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

  # GLOBAL SCALE
  # Global independent variables (environmental layers)
  IndVar.Global <- nsdm_finput$IndVar.Global

  # Select the best subset of independent variables for each species using covsel package
  myResp.xy.Global <- rbind(nsdm_finput$SpeciesData.XY.Global, nsdm_finput$Background.XY.Global)
  names(myResp.xy.Global)<-c("x","y")
  row.names(myResp.xy.Global)<-c(1:nrow(myResp.xy.Global))
  myResp.Global <- as.vector(c(rep(1,nrow(nsdm_finput$SpeciesData.XY.Global)),rep(0,nrow(nsdm_finput$Background.XY.Global))))
  myExpl.covsel.Global <- terra::extract(IndVar.Global, myResp.xy.Global, as.df=TRUE)[, -1]

  # Variable selection process
  Covdata.filter.Global<-covsel::covsel.filteralgo(covdata=myExpl.covsel.Global, pa=myResp.Global, corcut=corcut)

  # Embedding selected variables
  if(maxncov.Global=="nocorr") {
    maxncov.Global <- ncol(Covdata.filter.Global)
  }

  Covdata.embed.Global<-covsel::covsel.embed(covdata=Covdata.filter.Global,
                                        pa=myResp.Global,
                                        algorithms=algorithms,
                                        maxncov=maxncov.Global,
                                        nthreads=detectCores()/2) #@@@#TG why only half of cores?

  # Save selected variables for each species
  Selected.Variables.Global <- labels(Covdata.embed.Global$covdata)[[2]]
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

  IndVar.Global.Selected <- IndVar.Global[[Selected.Variables.Global]]  #@@@TG lo he vuelto a poner en resolucion global para entrenar el modelo global

  # Summary
  summary <- data.frame(Values = c(SpeciesName,
				nlyr(sabina$IndVar.Global), 
				length(Selected.Variables.Global)))

  rownames(summary) <- c("Species name",
			"Original number of variables at global scale", 
			"Final number of selected variables at global scale")

  # REGIONAL SCALE
  # Regional independent variables (environmental layers)
  IndVar.Regional <- nsdm_finput$IndVar.Regional

  # Subset the global independent variables for regional projections
  IndVar.Global.Selected.reg <- IndVar.Regional[[Selected.Variables.Global]]  #@@@#TG he puesto que se guarde tanto en resolucion global como regional para entrenar y proyectar 

  # Exclude climatic bands specified by the user.
  if(!is.null(ClimaticVariablesBands) && length(ClimaticVariablesBands) > 0) {
    Number.bands <- nlyr(IndVar.Regional)
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
  }

  Covdata.embed.Regional <- covsel::covsel.embed(covdata = Covdata.filter.Regional,
                                                   pa = myResp.Regional,
                                                   algorithms = algorithms,
                                                   maxncov = maxncov.Regional,
                                                   nthreads = detectCores() / 2) #@@@# why?

  # Save selected variables for each species
  Selected.Variables.Regional <- labels(Covdata.embed.Regional$covdata)[[2]]
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
  summary_regional <- data.frame(Values = c(nlyr(sabina$IndVar.Regional),
				length(Selected.Variables.Regional)))

  rownames(summary_regional) <- c("Original number of variables at regiona scale", 
			"Final number of of variables at regional scale")

  summary <- rbind(summary, summary_regional)

  sabina$Selected.Variables.Global <- Selected.Variables.Global
  sabina$IndVar.Global.Selected <- IndVar.Global.Selected #@@@# at global resolution for training
  sabina$Selected.Variables.Regional <- Selected.Variables.Regional
  sabina$IndVar.Regional.Selected <- IndVar.Regional.Selected
  sabina$IndVar.Global.Selected.reg <- IndVar.Global.Selected.reg #@@@# at regional resolution for projecting
  sabina$Summary<-summary
  # Make the output lighter
  sabina <- sabina[!names(sabina) %in% c("IndVar.Regional", "IndVar.Global")]

  attr(sabina, "class") <- "nsdm.input"

  # Logs success or error messages
  #message("\nVarSelection executed successfully!\n")

  if(save.output){
    message("Results saved in the following locations:")
    message(paste0(
       " - Selected variables at global level: Results/Global/Values/", SpeciesName, ".variables.csv\n",
       " - Selected variables at regional level: Results/Regional/Values/", SpeciesName, ".variables.csv\n"
      ))
  }

  return(sabina)

}
