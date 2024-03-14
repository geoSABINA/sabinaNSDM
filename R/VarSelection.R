#' @export
NSH.SDM.SelectVariables <- function(nshsdm_input, # #@@@RGM he quitado 	VariablesPath, 
				maxncov, #@@@RGM No podría por defecto un número de variables, no? Que las selccione el usuario
				corcut=NULL, #@@@RGM esto lo podemos dejar por defecto en 0.7 o incluso 0.8 
				algorithms=NULL, #@@@RGM Podemos poner tres algoritmos frecuentes o que el usario considere poner los mismos que usará en la modelización después 
				ClimaticVariablesBands=NULL,
				save.output=TRUE) {
		#@@@JMB# Pendiente: Buscar solución a sobreescritura de Results si ClimaticVariablesBands=NULL/o no null
  	  		
  nshsdm_name <- as.list(match.call())$nshsdm_input
  if(!inherits(nshsdm_input, "nshsdm.input")){
      stop("nshsdm_input must be an object of nshsdm.input class. Consider running NSH.SDM.PrepareData() function.")
  }

  # Esto no se si es necesario
#  if(any(!algorithms %in% c('glm','gam','rf'))) {   #@@@JMB# Son todos lo que hay? los dejamos por defecto?
#    stop("Please select a valid algorithms (\"glm\", \"gam\", or \"rf\").") #@RGM son en minusculas
#  }
#  algorithms <- tolower(algorithms)

  SpeciesName <- nshsdm_input$Species.Name
  
  nshsdm_data<-nshsdm_input
  nshsdm_data$args <- list()
  nshsdm_data$args$VariablesPath <- VariablesPath
  nshsdm_data$args$maxncov <- maxncov
  nshsdm_data$args$corcut <- corcut
  nshsdm_data$args$algorithms <- algorithms

  # GLOBAL SCALE
  # Global independent variables (environmental layers)  
  IndVar.Global <- nshsdm_input$IndVar.Global #@Ahora el stack se genera en la función anterior
  
  # Select the best subset of independent variables for each species using covsel package  #@@@RGM he cambio nombre de objetos
  myResp.xy.Global <- rbind(nshsdm_input$SpeciesData.XY.Global, nshsdm_input$Background.XY.Global) 
  names(myResp.xy.Global)<-c("x","y")
  row.names(myResp.xy.Global)<-c(1:nrow(myResp.xy.Global))
  myResp.Global <- as.vector(c(rep(1,nrow(nshsdm_input$SpeciesData.XY.Global)),rep(0,nrow(nshsdm_input$Background.XY.Global))))
  myResp.Global <- as.numeric(as.vector(myResp.Global))
  myExpl.covsel.Global <- terra::extract(IndVar.Global, myResp.xy.Global, as.df=TRUE)[, -1]
  
  # Variable selection process
  
  if(is.null(corcut)) { #@@@RGM he cambiado esto para que funcione con "null" (corcut=0.7) o corcut selecionado por el usuario
    Covdata.filter.Global<-covsel::covsel.filteralgo(covdata=myExpl.covsel.Global, pa=myResp.Global, corcut=0.7)
  } else {
    Covdata.filter.Global<-covsel::covsel.filteralgo(covdata=myExpl.covsel.Global, pa=myResp.Global, corcut=corcut)
  }
    
  # Embedding selected variables @RGM esto es necesario? 
  if(is.null(maxncov)) { 		#@@@JMB# Las maxncov posibles son las que salen en Covdata.filter?? 
    maxncov <- ncol(Covdata.filter.Global)
  }

  if(is.null(algorithms)) { #@@@RGM he cambiado esto para que funcione con "null" o algoritmos selecionados por el usuario
  Covdata.embed.Global<-covsel::covsel.embed(covdata=Covdata.filter.Global,
                              pa=myResp.Global,
                              algorithms=c('glm','gam','rf'), 
                              maxncov=maxncov, 
                              nthreads=detectCores()/2)  
  } else {
    Covdata.embed.Global<-covsel::covsel.embed(covdata=Covdata.filter.Global,
                                        pa=myResp,
                                        algorithms=algorithms, 
                                        maxncov=maxncov, 
                                        nthreads=detectCores()/2)  
  } 
  
  # Save selected variables for each species
  Selected.Variables.Global <- labels(Covdata.embed.Global$covdata)[[2]] 
  if(save.output){
  write.csv(Selected.Variables.Global, paste0("Results/Global/Values/",SpeciesName,".variables.csv"))
  }
  
  # REGIONAL SCALE
  # Regional independent variables (environmental layers)
  IndVar.Regional <- nshsdm_input$IndVar.Regional
  
  # Subset the global independent variables for regional projections
  IndVar.Global.Selected <- IndVar.Regional[[Selected.Variables.Global]]  #@RGM he cambiado el nombre y he quitado la proyección global, solo vamos a proyectar a escala regional. Proyectar a escalar global es como hacer un modelo "normal" lo pueden hacer con cualquier paquete de modelos
  
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
  myResp.Regional <- as.numeric(as.vector(myResp.Regional))
  myExpl.covsel.Regional <- terra::extract(IndVar.Regional, myResp.xy.Regional, rm.na=TRUE, df=TRUE)[, -1] 

  # Variable selection process
  if(is.null(corcut)) { #@@@RGM he cambiado esto para que funcione con "null" (corcut=0.7) o corcut selecionado por el usuario
  Covdata.filter.Regional <- covsel::covsel.filteralgo(covdata = myExpl.covsel.Regional, pa = myResp.Regional, corcut = 0.7)
  } else {
    Covdata.filter.Regional <- covsel::covsel.filteralgo(covdata = myExpl.covsel.Regional, pa = myResp.Regional, corcut = corcut)
  }
  
  # Embedding selected variables
  if(is.null(algorithms)) { #@@@RGM he cambiado esto para que funcione con "null" o algoritmos selecionados por el usuario
  Covdata.embed.Regional <- covsel::covsel.embed(covdata = Covdata.filter.Regional,
                                         pa = myResp.Regional,
                                         algorithms = c('glm','gam','rf'),
                                         maxncov = maxncov,
                                         nthreads = detectCores() / 2)
  } else {
    Covdata.embed.Regional <- covsel::covsel.embed(covdata = Covdata.filter.Regional,
                                                   pa = myResp.Regional,
                                                   algorithms = algorithms,
                                                   maxncov = maxncov,
                                                   nthreads = detectCores() / 2)
  }
  
  # Save selected variables for each species
  Selected.Variables.Regional <- labels(Covdata.embed.Regional$covdata)[[2]]
  if(save.output){
  write.csv(Selected.Variables.Regional, paste0("Results/Regional/Values/", SpeciesName, ".variables.csv"))
  }
  
  # Subset the regional independent variables for regional projections
  IndVar.Regional.Selected <- IndVar.Regional[[Selected.Variables.Regional]] #@RGM he cambiado el nombre

  nshsdm_data$Selected.Variables.Global <- Selected.Variables.Global 
  nshsdm_data$Selected.Variables.Regional <- Selected.Variables.Regional 
  nshsdm_data$IndVar.Global.Selected <- IndVar.Global.Selected 
  nshsdm_data$IndVar.Regional.Selected <- IndVar.Regional.Selected 

  attr(nshsdm_data, "class") <- "nshsdm.input"

  return(nshsdm_data)

}
