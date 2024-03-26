#' @export
NSH.SDM.SelectVariables <- function(nshsdm_input,
				maxncov="All", #@@@#TG si usuario no pone ninguna usa todas las no correlacionadas
				corcut=0.7, #@@@TG lo dejo por defecto en 0.7
				algorithms=c('glm','gam','rf'), #@@@TG He puesto los tres como default, pero el usuario puede poner lo que quiera.
				ClimaticVariablesBands=NULL,
				save.output=TRUE) {

  nshsdm_name <- as.list(match.call())$nshsdm_input
  if(!inherits(nshsdm_input, "nshsdm.input")){
      stop("nshsdm_input must be an object of nshsdm.input class. Consider running NSH.SDM.PrepareData() function.")
  }

  algorithms <- tolower(algorithms)
  if(any(!algorithms %in% c('glm','gam','rf'))) {   #@@@JMB# Son todos lo que hay? los dejamos por defecto? #@@@ sí, los he puesto por defecto, pero esto es necesario por si ponen un nombre mal, por ej "maxent", o "GLM" en lugar de "glm"
    stop("Please select a valid algorithm (\"glm\", \"gam\", or \"rf\").")
  }


  SpeciesName <- nshsdm_input$Species.Name

  nshsdm_data<-nshsdm_input
  nshsdm_data$args <- list()
  nshsdm_data$args$VariablesPath <- VariablesPath
  nshsdm_data$args$maxncov <- maxncov
  nshsdm_data$args$corcut <- corcut
  nshsdm_data$args$algorithms <- algorithms

  # GLOBAL SCALE
  # Global independent variables (environmental layers)
  IndVar.Global <- nshsdm_input$IndVar.Global

  # Select the best subset of independent variables for each species using covsel package
  myResp.xy.Global <- rbind(nshsdm_input$SpeciesData.XY.Global, nshsdm_input$Background.XY.Global)
  names(myResp.xy.Global)<-c("x","y")
  row.names(myResp.xy.Global)<-c(1:nrow(myResp.xy.Global))
  myResp.Global <- as.vector(c(rep(1,nrow(nshsdm_input$SpeciesData.XY.Global)),rep(0,nrow(nshsdm_input$Background.XY.Global))))
  #myResp.Global <- as.numeric(as.vector(myResp.Global)) #@@@ esto sobra creo
  myExpl.covsel.Global <- terra::extract(IndVar.Global, myResp.xy.Global, as.df=TRUE)[, -1]

  # Variable selection process
  Covdata.filter.Global<-covsel::covsel.filteralgo(covdata=myExpl.covsel.Global, pa=myResp.Global, corcut=corcut)


  # Embedding selected variables
  if(maxncov=="all") {
    maxncov <- ncol(Covdata.filter.Global)
  }

  Covdata.embed.Global<-covsel::covsel.embed(covdata=Covdata.filter.Global,
                                        pa=myResp.Global,
                                        algorithms=algorithms,
                                        maxncov=maxncov,
                                        nthreads=detectCores()/2) #@@@#TG why only half of cores?

  # Save selected variables for each species
  Selected.Variables.Global <- labels(Covdata.embed.Global$covdata)[[2]]
  if(save.output){
    #suffix <- 0
    file_path <- paste0("Results/Global/Values/", SpeciesName, ".variables.csv")
    #old_path <- paste0("Results/Global/Values/", SpeciesName, ".variables.csv")
    # while (file.exists(file_path)) {
    #   suffix <- suffix + 1
    #   file_path <- gsub("\\.csv$", paste0("_", suffix, ".csv"), old_path)
    # }

    write.csv(Selected.Variables.Global, file_path)
    message(paste("Selected variables at global level saved in:",file_path,cat("\033[1;34m")))
    cat("\033[0m")
    }

  IndVar.Global.Selected <- IndVar.Global[[Selected.Variables.Global]]  #@@@TG lo he vuelto a poner en resolucion global para entrenar el modelo global

  # REGIONAL SCALE
  # Regional independent variables (environmental layers)
  IndVar.Regional <- nshsdm_input$IndVar.Regional

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
  myResp.xy.Regional <- rbind(nshsdm_input$SpeciesData.XY.Regional, nshsdm_input$Background.XY.Regional)
  row.names(myResp.xy.Regional) <- c(1:nrow(myResp.xy.Regional))
  myResp.Regional <- as.vector(c(rep(1, nrow(nshsdm_input$SpeciesData.XY.Regional)), rep(0, nrow(nshsdm_input$Background.XY.Regional))))
  #myResp.Regional <- as.numeric(as.vector(myResp.Regional)) #@@@#TG esto sobra
  myExpl.covsel.Regional <- terra::extract(IndVar.Regional, myResp.xy.Regional, rm.na=TRUE, df=TRUE)[, -1]


  Covdata.filter.Regional <- covsel::covsel.filteralgo(covdata = myExpl.covsel.Regional, pa = myResp.Regional, corcut = corcut)


  # Embedding selected variables

  Covdata.embed.Regional <- covsel::covsel.embed(covdata = Covdata.filter.Regional,
                                                   pa = myResp.Regional,
                                                   algorithms = algorithms,
                                                   maxncov = maxncov,
                                                   nthreads = detectCores() / 2) #@@@# why?

  # Save selected variables for each species
  Selected.Variables.Regional <- labels(Covdata.embed.Regional$covdata)[[2]]
  if(save.output){
    #suffix <- 0
    file_path <- paste0("Results/Regional/Values/", SpeciesName, ".variables.csv")
    #old_path<-paste0("Results/Regional/Values/", SpeciesName, ".variables.csv")
    # while (file.exists(file_path)) {
    #   suffix <- suffix + 1
    #   file_path <- gsub("\\.csv$", paste0("_", suffix, ".csv"), old_path)
    # }

  write.csv(Selected.Variables.Regional, file_path)
  message(paste("Selected variables at regional level saved in:",file_path,cat("\033[1;34m")))
  cat("\033[0m")
  }

  # Subset the regional independent variables for regional projections
  IndVar.Regional.Selected <- IndVar.Regional[[Selected.Variables.Regional]]



  nshsdm_data$Selected.Variables.Global <- Selected.Variables.Global
  nshsdm_data$Selected.Variables.Regional <- Selected.Variables.Regional
  nshsdm_data$IndVar.Global.Selected <- IndVar.Global.Selected #@@@# at global resolution for training
  nshsdm_data$IndVar.Regional.Selected <- IndVar.Regional.Selected
  nshsdm_data$IndVar.Global.Selected.reg <- IndVar.Global.Selected.reg #@@@# at regional resolution for projecting
  results<-  nshsdm_data$Summary
  results<-rbind(results,
                 c("Original number variables at global level",nlyr(nshsdm_data$IndVar.Global)),
                 c("Final number variables at global level",length(nshsdm_data$Selected.Variables.Global)),
                 c("Original number variables at regiona level",nlyr(nshsdm_data$IndVar.Regional)),
                 c("Final number variables at regional level",length(nshsdm_data$Selected.Variables.Regional)))
  if(save.output){
    write.table(results, paste0("Results/",SpeciesName,"_summary.csv"), sep=",",  row.names=F, col.names=T)
  }
  nshsdm_data$Summary<-results

   attr(nshsdm_data, "class") <- "nshsdm.input"

    # Logs success or error messages
  # message("\nVarSelection executed successfully!\n",cat("\033[32m"))
  # cat("\033[0m")
  #
  # if(save.output){
  #   message("Results saved in the following locations:",cat("\033[1;34m"))
  #   message(paste0(
  #     " - Global level: Results/Global/Values/\n",
  #     " - Regional level: Results/Regional/Values/\n"
  #   ),cat("\033[1;34m"))
  #   cat("\033[0m")
  # }

  return(nshsdm_data)

}
