#' @export
NSH.SDM.SelectVariables <- function(NSH.SDM.Data, VariablesPath, Max.nCov, Cor.Cutoff, ClimaticVariablesBands=NULL,SpeciesName) { # #@@@##(CAREFUL!!add parameter speciesName), also: set as default ClimaticVariablesBands = NULL,
  #library(covsel)
  #library(terra)
  #library(ecospat)
  #library(fs)
  
  # Global scale
  # Global independent variables (environmental layers)  
  IndVar.Global <- rast(paste0(VariablesPath,"/Global/Current.tif")) 
  IndVar.Global <- IndVar.Global[[names(IndVar.Global)]]  #@@@##  Esto para qué se hace?
  
  # Select the best subset of independent variables for each species using covsel package 
  myResp.xy <- rbind(NSH.SDM.Data$SpeciesData.XY.Global ,NSH.SDM.Data$Background.XY.Global) 
  names(myResp.xy)<-c("x","y")
  row.names(myResp.xy)<-c(1:nrow(myResp.xy))
  myResp <- as.vector(c(rep(1,nrow(NSH.SDM.Data$SpeciesData.XY.Global)),rep(0,nrow(NSH.SDM.Data$Background.XY.Global))))
  myResp <- as.numeric(as.vector(myResp))
  myExpl.covsel <- data.frame (terra::extract (IndVar.Global, myResp.xy))[, -1]
  
  # Variable selection process
  Covdata.filter<-covsel.filteralgo(covdata=myExpl.covsel,pa=myResp,corcut=Cor.Cutoff)
  
  # Embedding selected variables
  Covdata.embed<-covsel.embed(covdata=Covdata.filter,
                              pa=myResp,
                              algorithms=c('glm','gam','rf'), 
                              maxncov=Max.nCov, 
                              nthreads=detectCores()/2)  
  
  # Get selected variables #@@@##(remove this)
  Selected.Variables.Global <- labels(Covdata.embed$covdata)[[2]] 
  
  # Save selected variables for each species
  write.csv(Selected.Variables.Global, paste("Results/Global/Values/", SpeciesName, ".variables.csv", sep = "")) #@@@##(CAREFUL!!if its saved with the species name, it should be given as a parameter)
  
  # Subset the global independent variables for regional projections
  IndVar.Regional.1 <- rast(paste0(VariablesPath,"/Regional/Current.tif")) 
  IndVar.Global.2 <- IndVar.Global[[Selected.Variables.Global]] #@@@## at global level to train the global model
  IndVar.Global.3 <- IndVar.Regional.1[[Selected.Variables.Global]] #@@@## selected for global, but charged at regional scale to project the global model
  
  # Regional scale
  # Regional independent variables (environmental layers)
  IndVar.Regional <- rast(paste0(VariablesPath, "/Regional/Current.tif"))
  IndVar.Regional <- IndVar.Regional[[names(IndVar.Regional)]]  #@@@# what is this for??
  
  # Exclude climatic bands specified by the user. #@@@# I moved and change this so it only gives one result excluding or not the climatic vars
  if (!is.null(ClimaticVariablesBands) && length(ClimaticVariablesBands) > 0) {
    Number.bands <- nlyr(IndVar.Regional)
    Bands.climatic <- setdiff(1:Number.bands, ClimaticVariablesBands)#@@@# non eliminated variables
    IndVar.Regional <- IndVar.Regional[[Bands.climatic]] 
  } else {
    # If no bands are specified for exclusion, use all bands
    IndVar.Regional <- IndVar.Regional
  }
  
  # Select the best subset of independent variables for each species using covsel package
  myResp.xy.Regional <- rbind(NSH.SDM.Data$SpeciesData.XY.Regional, NSH.SDM.Data$Background.XY.Regional)
  names(myResp.xy.Regional) <- c("x", "y")
  row.names(myResp.xy.Regional) <- c(1:nrow(myResp.xy.Regional))
  myResp.Regional <- as.vector(c(rep(1, nrow(NSH.SDM.Data$SpeciesData.XY.Regional)), rep(0, nrow(NSH.SDM.Data$Background.XY.Regional))))
  myResp.Regional <- as.numeric(as.vector(myResp.Regional))
  myExpl.covsel.Regional <- data.frame(terra::extract(IndVar.Regional, myResp.xy.Regional))[, -1]
  
  # Variable selection process
  Covdata.filter.Regional <- covsel.filteralgo(covdata = myExpl.covsel.Regional, pa = myResp.Regional, corcut = Cor.Cutoff)
  
  # Embedding selected variables
  Covdata.embed.Regional <- covsel.embed(covdata = Covdata.filter.Regional,
                                         pa = myResp.Regional,
                                         algorithms = c('glm', 'gam', 'rf'),
                                         maxncov = Max.nCov,
                                         nthreads = detectCores() / 2)
  
  # Get selected variables #@@@##(Remove this)
  Selected.Variables.Regional <- labels(Covdata.embed.Regional$covdata)[[2]]
  
  # Save selected variables for each species
  write.csv(Selected.Variables.Regional, paste0("Results/Regional/Values/", SpeciesName, ".variables.csv"), row.names = FALSE) #@@@##(CAREFUL!!if its saved with the species name, it should be given as a parameter)
  
  # Subset the regional independent variables for regional projections
  IndVar.Regional.2 <- IndVar.Regional[[Selected.Variables.Regional]]

  return(list(Selected.Variables.Global = Selected.Variables.Global, Selected.Variables.Regional = Selected.Variables.Regional, IndVar.Global.2 = IndVar.Global.2, IndVar.Global.3 =IndVar.Global.3, IndVar.Regional.2= IndVar.Regional.2)) #@@@# Changed this to give only one regional result and the global variables at both global (to train the model) and regional (to project the model) scale
}
