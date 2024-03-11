#' @export
NSH.SDM.Covariate.Models <- function(nshsdm_selvars, 
				models=c("GAM","GBM", "RF", "MAXNET","GLM"),
				CV.nb.rep=1, 
				CV.perc=0.8) {

  nshsdm_global <- as.list(match.call())$nshsdm_selvars
  if(!inherits(nshsdm_selvars, "nshsdm.input")){
      stop("nshsdm_selvars must be an object of nshsdm.input class.")
  }

  SpeciesName <- nshsdm_selvars$Species.Name
	##prueba RGM

  tryCatch({ 
	# Covariate model calibrated with all the independent variables.
    	# Regional model excluding climatic variables.
	
	# Create directories to save results of the modeling process
      	dir_create(c("Results/Covariate/Values/", "Results/Covariate/Projections/"))
      
	# Add the global model as an additional variable for the regional model.
	SDM.global <- terra::rast(paste0("Results/Global/Projections/",SpeciesName,".Current.tif"))
	names(SDM.global) <- c("SDM.global") 
	IndVar.Regional.temp <- nshsdm_selvars$IndVar.Regional.2
	IndVar.Regional.2 <- c(IndVar.Regional.temp, SDM.global)
      
	myResp.xy <- rbind(nshsdm_selvars$SpeciesData.XY.Regional, nshsdm_selvars$Background.XY.Regional)
	names(myResp.xy)<-c("x","y")
	row.names(myResp.xy)<-c(1:nrow(myResp.xy))
	myResp <- data.frame(c(rep(1,nrow(nshsdm_selvars$SpeciesData.XY.Regional)),rep(NA,nrow(nshsdm_selvars$Background.XY.Regional ))))
	names(myResp)<-"pa"
	row.names(myResp)<-c(1:nrow(myResp.xy))
	myExpl <- terra::extract(IndVar.Regional.2, myResp.xy, as.df=TRUE)[, -1]
      
	# Data required for the biomod2 package.
	myBiomodData <- biomod2::BIOMOD_FormatingData(resp.var = myResp, 
					resp.xy = myResp.xy, 
					expl.var = myExpl, 
					resp.name = SpeciesName, 
					PA.nb.rep = 1, 
					PA.nb.absences = nrow(nshsdm_selvars$Background.XY.Regional), 
					PA.strategy = "random")
      
	# Calibrate and evaluate individual models with specified statistical algorithms.
	myBiomodModelOut <- biomod2::BIOMOD_Modeling(bm.format = myBiomodData, 
					modeling.id = "AllModels", 
					#  bm.options = myThecniquesOptions,
					models = models, 
					CV.strategy = "random", 
					CV.nb.rep = CV.nb.rep, 
					CV.perc = CV.perc, 
					weights = NULL, 
					var.import = 3, 
					metric.eval = c("ROC", "TSS", "KAPPA", "ACCURACY", "SR", "BOYCE", "MPA"), 
					scale.models = FALSE, 
					do.progress = TRUE, 
					prevalence = 0.5, 
					seed.val = 42, 
 					CV.do.full.models = FALSE)
      
	# Generate and evaluate a single ensemble (i.e.,consensus) model that averages the individual models weighted
	# by the value of the AUC statistic.
	myBiomodEM.ROC  <- biomod2::BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                                 models.chosen = 'all',
                                                 em.by = 'all',
                                                 em.algo = c("EMwmean"),
                                                 metric.select = c('ROC'),
                                                 metric.select.thresh = NULL,
                                                 var.import = 3,
                                                 metric.eval = c('ROC', "TSS", "KAPPA"),
                                                 seed.val = 42) 
      
   	# Project the individual models to the study area at regional scale under training conditions.
	myBiomodProj <- biomod2::BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                        new.env = IndVar.Regional.2,
                                        proj.name = "Current",
                                        models.chosen = 'all',
                                        build.clamping.mask = FALSE) 
      
	# Project the ensemble model to the study area at regional scale under training conditions.
	biomod2::BIOMOD_EnsembleForecasting(bm.em = myBiomodEM.ROC, 
					bm.proj = myBiomodProj,
					models.chosen = 'all',
					metric.binary = 'all',
					metric.filter = 'all',
					build.clamping.mask = FALSE)
      
	# Load the model stored by biomod2 and save it in geotif format
	sp.name<-myBiomodData@sp.name
	Pred.ROC <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble.tif")) 
	terra::writeRaster(Pred.ROC, paste0("Results/Covariate/Projections/",SpeciesName,".Current.tif"), overwrite=TRUE)
      
	# Binary models
	Pred.bin.ROC <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_ROCbin.tif")) 
	terra::writeRaster(Pred.bin.ROC, paste0("Results/Covariate/Projections/",SpeciesName,".Current.bin.ROC.tif"), overwrite=TRUE)
      
	Pred.bin.TSS <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_TSSbin.tif")) 
	terra::writeRaster(Pred.bin.TSS, paste0("Results/Covariate/Projections/",SpeciesName,".Current.bin.TSS.tif"), overwrite=TRUE)
      
       # Save some results
	# Values of the statistics for each of the replicates   
	myEMeval.replicas <- biomod2::get_evaluations(myBiomodModelOut)
	write.csv(myEMeval.replicas,file=paste0("Results/Covariate/Values/",SpeciesName,"_replica.csv"))
	# Values of the statistics of the consensus model   
	myEMeval.Consenso <- biomod2::get_evaluations(myBiomodEM.ROC)
	write.csv(myEMeval.Consenso,file=paste0("Results/Covariate/Values/",SpeciesName,"_ensemble.csv"))
	# Variable importance
	myModelsVarImport <- biomod2::get_variables_importance(myBiomodModelOut)
	write.table(myModelsVarImport, file = paste0("Results/Covariate/Values/",SpeciesName,"_indvar.csv"), row.names = T, col.names = T)
      
	# Model projections for future climate scenarios
	################################################
	Scenarios <- dir_ls(paste0(VariablesPath,"/Regional"), pattern="tif")
        Scenarios <- Scenarios[!grepl("Current.tif", Scenarios)][1:2] # para probar con 2

	if(length(Scenarios) == 0) {
	  stop("There are no future climate variables")
	}
        
	walk(Scenarios, function(projmodel) { 
	  NewClim.temp <- terra::rast(projmodel)[[nshsdm_selvars$Selected.Variables.Regional]]
	  Scenario.name <- path_file(projmodel) |> path_ext_remove()
     
	  SDM.global.future <- terra::rast(paste0("Results/Global/Projections/",SpeciesName,".",Scenario.name,".tif"))
	  names(SDM.global.future) <- c("SDM.global")
	  NewClim <- c(NewClim.temp, SDM.global.future)
        
	  myBiomomodProjScenario <- biomod2::BIOMOD_Projection(bm.mod = myBiomodModelOut,
						new.env = NewClim,
						proj.name = Scenario.name,
						models.chosen = "all")
        
	  biomod2::BIOMOD_EnsembleForecasting(bm.em = myBiomodEM.ROC,
					bm.proj = myBiomomodProjScenario,
					models.chosen = "all",
					metric.binary = "all",
					metric.filter = "all")
        
	  Pred.ROC.Scenario <- terra::rast(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble.tif"))
	  terra::writeRaster(Pred.ROC.Scenario, paste0("Results/Covariate/Projections/",SpeciesName,".",Scenario.name,".tif"), overwrite = TRUE)
        
        
	  # Binarized models
          Pred.bin.TSS.Scenario <- terra::rast(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_ROCbin.tif"))
          terra::writeRaster(Pred.bin.TSS.Scenario, paste0("Results/Covariate/Projections/",SpeciesName,".",Scenario.name,".bin.TSS.tif"), overwrite = TRUE)
        
          Pred.bin.ROC.Scenario <- terra::rast(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_ROCbin.tif"))
          terra::writeRaster(Pred.bin.ROC.Scenario, paste0("Results/Covariate/Projections/",SpeciesName,".",Scenario.name,".bin.ROC.tif"), overwrite = TRUE)

	  gc()	
	}) 
      
	#@@@@## Move biomod2 results to Results/Covariate/Models folder
	source_folder <- sp.name
	destination_folder <- "Results/Covariate/Models"
	if(file.exists(destination_folder)) {
          unlink(destination_folder, recursive = TRUE)
	}
	file.rename(from = source_folder, to = destination_folder)
    
  # Logs success or error messages 
  message("\nNSH.SDM.Covariate.Models executed successfully.\n")
  }, error = function(err) {
    message("Error in NSH.SDM.Covariate.Models:", conditionMessage(err))
    return(list(result = NULL, error = err))
  })
}

