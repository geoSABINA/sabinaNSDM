#' @export
NSH.SDM.SelectVariables <- function(nshsdm_input, 
				VariablesPath, 
				maxncov=7, 
				corcut=0.7,
				algorithms=c("GLM","GAM","RF"), #@@@JMB lo ponemos cómo argumento?
				ClimaticVariablesBands=NULL) {
		#@@@JMB# Pendiente: Buscar solución a sobreescritura de Results si ClimaticVariablesBands=NULL/o no null
  	  		
  nshsdm_name <- as.list(match.call())$nshsdm_input
  if(!inherits(nshsdm_input, "nshsdm.input")){
      stop("nshsdm_input must be an object of nshsdm.input class. Consider running NSH.SDM.PrepareData() function.")
  }

  if(any(!algorithms %in% c("GLM", "GAM", "RF"))) {   #@@@JMB# Son todos lo que hay? va por defecto? Igualar con nombres con otras funciones?
    stop("Please select a valid algorithms (\"GLM\", \"GAM\", or \"RF\").")
  }
  algorithms <- tolower(algorithms)

  SpeciesName <- nshsdm_input$Species.Name
  
  
  nshsdm_data<-nshsdm_input
  #nshsdm_data<-list()

  # GLOBAL SCALE
  # Global independent variables (environmental layers)  
  IndVar.Global <- terra::rast(paste0(VariablesPath,"/Global/Current.tif"))
  Mask <- prod(IndVar.Global)
  IndVar.Global <- terra::mask(IndVar.Global, Mask)
  
  # Select the best subset of independent variables for each species using covsel package 
  myResp.xy <- rbind(nshsdm_input$SpeciesData.XY.Global, nshsdm_input$Background.XY.Global) 
  names(myResp.xy)<-c("x","y")
  row.names(myResp.xy)<-c(1:nrow(myResp.xy))
  myResp <- as.vector(c(rep(1,nrow(nshsdm_input$SpeciesData.XY.Global)),rep(0,nrow(nshsdm_input$Background.XY.Global))))
  myResp <- as.numeric(as.vector(myResp))
  myExpl.covsel <- terra::extract(IndVar.Global, myResp.xy, as.df=TRUE)[, -1]
  
  # Variable selection process
  Covdata.filter<-covsel::covsel.filteralgo(covdata=myExpl.covsel, pa=myResp, corcut=corcut)
  
  # Embedding selected variables
  if(is.null(maxncov)) { 		#@@@JMB# Las maxncov posibles son las que salen en Covdata.filter?? 
    maxncov <- ncol(Covdata.filter)
  }

  Covdata.embed<-covsel::covsel.embed(covdata=Covdata.filter,
                              pa=myResp,
                              algorithms=c('glm','gam','rf'), 
                              maxncov=maxncov, 
                              nthreads=detectCores()/2)  
  
  
  # Save selected variables for each species
  Selected.Variables.Global <- labels(Covdata.embed$covdata)[[2]]  
  write.csv(Selected.Variables.Global, paste0("Results/Global/Values/",SpeciesName,".variables.csv"))
  
  # REGIONAL SCALE
  # Regional independent variables (environmental layers)
  IndVar.Regional <- terra::rast(paste0(VariablesPath,"/Regional/Current.tif")) 
  Mask.Regional <- prod(IndVar.Regional)
  IndVar.Regional <- terra::mask(IndVar.Regional, Mask.Regional)
  
  # Subset the global independent variables for regional projections
  #if(!all(Selected.Variables.Global %in% names(IndVar.Regional))) { #@@@JMB# Esto es immportante??
  #  stop("Global and Regional variables must have matching names.")
  #}
  IndVar.Global.2 <- IndVar.Global[[Selected.Variables.Global]] #@@@## at global level to train the global model #@@@JMB rename object and varnames()
  IndVar.Global.3 <- IndVar.Regional[[Selected.Variables.Global]] #@@@## selected for global, but charged at regional scale to project the global model
  
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
  Covdata.filter.Regional <- covsel::covsel.filteralgo(covdata = myExpl.covsel.Regional, pa = myResp.Regional, corcut = corcut)
  
  # Embedding selected variables
  Covdata.embed.Regional <- covsel::covsel.embed(covdata = Covdata.filter.Regional,
                                         pa = myResp.Regional,
                                         algorithms = algorithms,
                                         maxncov = maxncov,
                                         nthreads = detectCores() / 2)
  
  
  # Save selected variables for each species
  Selected.Variables.Regional <- labels(Covdata.embed.Regional$covdata)[[2]]
  write.csv(Selected.Variables.Regional, paste0("Results/Regional/Values/", SpeciesName, ".variables.csv"))
  
  # Subset the regional independent variables for regional projections
  IndVar.Regional.2 <- IndVar.Regional[[Selected.Variables.Regional]]

  nshsdm_data$Selected.Variables.Global <- Selected.Variables.Global 
  nshsdm_data$Selected.Variables.Regional <- Selected.Variables.Regional 
  nshsdm_data$IndVar.Global.2 <- IndVar.Global.2 
  nshsdm_data$IndVar.Global.3 <- IndVar.Global.3 
  nshsdm_data$IndVar.Regional.2 <- IndVar.Regional.2 #@@@# Changed this to give only one regional result and the global variables at both global (to train the model) and regional (to project the model) scale

  attr(nshsdm_data, "class") <- "nshsdm.input"

  return(nshsdm_data)

}
