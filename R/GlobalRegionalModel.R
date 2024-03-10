#' @export
####################
# 3. GLOBAL MODELS
####################

NSH.SDM.Global.Model <- function(nshsdm_selvars, 
				models=c("GAM","GBM", "RF", "MAXNET","GLM"),
				CV.nb.rep=1, 
				CV.perc=0.8) {
  
  nshsdm_global <- as.list(match.call())$nshsdm_selvars
  if(!inherits(nshsdm_selvars, "nshsdm.input")){
      stop("nshsdm_selvars must be an object of nshsdm.input class.")
  }

  SpeciesName <- nshsdm_selvars$Species.Name

  tryCatch({ 
	# Here we format the response (presence/background) and explanatory (environmental variables) data for BIOMOD2
	myResp.xy <- rbind(nshsdm_selvars$SpeciesData.XY.Global,nshsdm_selvars$Background.XY.Global)
	row.names(myResp.xy)<-c(1:nrow(myResp.xy))
	myResp <- data.frame(c(rep(1,nrow(nshsdm_selvars$SpeciesData.XY.Global)),rep(NA,nrow(nshsdm_selvars$Background.XY.Global))))
	names(myResp)<-"pa"
	row.names(myResp)<-c(1:nrow(myResp.xy)) 
	myExpl <- terra::extract(nshsdm_selvars$IndVar.Global.2, myResp.xy, as.df=TRUE)[, -1] 
	
	# Remaining script sections involve executing biomod2 modeling procedures, including data formatting, model training, 
	# projection, evaluation, and visualization. Please refer to biomod2 documentation for detailed explanation of these steps.
	
	# Prepare data required to calibrate the model with biomod2 package (Result: A BIOMOD.formated.data object)
	myBiomodData <- biomod2::BIOMOD_FormatingData(resp.var = myResp, 
	                                     resp.xy = myResp.xy, 
	                                     expl.var = myExpl, 
	                                     resp.name = SpeciesName, 
	                                     PA.nb.rep = 1, 
	                                     PA.nb.absences = nrow(nshsdm_selvars$Background.XY.Global), 
	                                     PA.strategy = "random")

	# Calibrate and evaluate individual models with specified statistical algorithms
	myBiomodModelOut <- biomod2::BIOMOD_Modeling(bm.format = myBiomodData,  #@@@## (This function saves outputs, check them and specify them in the description). Outputs: /NSHSDM/Larix.decidua/models/AllModels/Larix.decidua_PA1_RUN1_GBM
					modeling.id = "AllModels", 
					#bm.options = myThecniquesOptions,
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

	# Replicates with ROC > 0.8	
	df <- myBiomodModelOut@models.evaluation
	df_slot <- slot(df, "val")
	df_slot <- df_slot[df_slot$metric.eval == "ROC", ]
	percentage <- 100 * sum(df_slot$validation >= 0.8)/nrow(df_slot)
	warning(sprintf("\n%.2f%% of replicates with AUC values >= 0.8.\n", percentage))			
      
	# Generate and evaluate a single ensemble (i.e.,consensus) model that averages the individual models weighted by 
	# the value of the AUC statistic.
	myBiomodEM.ROC  <- biomod2::BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
					models.chosen = 'all',
					em.by = 'all',
					em.algo = c("EMwmean"),
					metric.select = c('ROC'),
					metric.select.thresh = NULL,
					var.import = 3,
					metric.eval = c('ROC', "TSS", "KAPPA"),
					seed.val = 42) 
	
	
	# Project the individual models to the study area at regional scale under training conditions)
	myBiomodProj <- biomod2::BIOMOD_Projection(bm.mod = myBiomodModelOut,
					new.env = nshsdm_selvars$IndVar.Global.3,
					proj.name = "Current",
					models.chosen = 'all',
					build.clamping.mask = FALSE) 
	
	# Project the ensemble model to the study area at regional scale under training conditions
	biomod2::BIOMOD_EnsembleForecasting(bm.em = myBiomodEM.ROC, 
				   bm.proj = myBiomodProj,
	                           models.chosen = 'all',
	                           metric.binary = 'all',
	                           metric.filter = 'all',
	                           build.clamping.mask = FALSE)
	
	# Load the model stored by biomod2 and save it in geotif format
	sp.name<-myBiomodData@sp.name
	Pred.ROC <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble.tif"))
	terra::writeRaster(Pred.ROC, paste0("Results/Global/Projections/",SpeciesName,".Current.tif"), overwrite=TRUE)
	
	# Binary models
	Pred.bin.ROC <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_ROCbin.tif"))
	terra::writeRaster(Pred.bin.ROC, paste0("Results/Global/Projections/",SpeciesName,".Current.bin.ROC.tif"), overwrite=TRUE)
		  
	Pred.bin.TSS <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_TSSbin.tif"))
	terra::writeRaster(Pred.bin.TSS, paste0("Results/Global/Projections/",SpeciesName,".Current.bin.TSS.tif"), overwrite=TRUE)
		  
	# Save some results
	# Values of the evaluation statistics for each replica   
	myEMeval.replicates <- biomod2::get_evaluations(myBiomodModelOut)
	write.csv(myEMeval.replicates,file=paste0("Results/Global/Values/",SpeciesName,"_replica.csv"))
	# Values of the evaluation statistics of the consensus model   
	myEMeval.Ensemble <- biomod2::get_evaluations(myBiomodEM.ROC)
	write.csv(myEMeval.Ensemble,file=paste0("Results/Global/Values/",SpeciesName,"_ensemble.csv"))
	# Variable importance
	myModelsVarImport <- biomod2::get_variables_importance(myBiomodModelOut)
	write.csv(myModelsVarImport, file = paste0("Results/Global/Values/",SpeciesName,"_indvar.csv"), row.names = T)
		      
	# Model projections for future climate scenarios
	################################################
	Scenarios <- dir_ls(paste0(VariablesPath,"/Regional"), pattern="tif")
        Scenarios <- Scenarios[!grepl("Current.tif", Scenarios)][1:2] # para probar con dos

	if(length(Scenarios) == 0) {
	  stop("There are no future climate variables")
	}
 
	walk(Scenarios, function(projmodel) { 
	  #projmodel <- Scenarios[1]
	  NewClim <- terra::rast(projmodel)[[nshsdm_selvars$Selected.Variables.Global]]
	  Scenario.name <- path_file(projmodel) |> path_ext_remove()
  
	  #Project the individual models
	  myBiomomodProjScenario <- biomod2::BIOMOD_Projection(bm.mod = myBiomodModelOut,
							 new.env = NewClim,
							 proj.name = Scenario.name,
							 models.chosen = "all")

	  #Project the ensemble model
	  biomod2::BIOMOD_EnsembleForecasting(bm.em = myBiomodEM.ROC,
					bm.proj = myBiomomodProjScenario,
					models.chosen = "all",
					metric.binary = "all",
					metric.filter = "all")
		  #@@@@## All these tifs are saved both in the biomod folder and the Projections folder, I would suggest deleting them from one of these places
	  Pred.ROC.Scenario <- terra::rast(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble.tif"))
	  terra::writeRaster(Pred.ROC.Scenario, paste0("Results/Global/Projections/",SpeciesName,".",Scenario.name,".tif"), overwrite = TRUE)

	  # Binarized models
	  Pred.bin.TSS.Scenario <- terra::rast(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_ROCbin.tif")) 
	  terra::writeRaster(Pred.bin.TSS.Scenario, paste0("Results/Global/Projections/",SpeciesName,".",Scenario.name,".bin.TSS.tif"), overwrite = TRUE)
		  
	  Pred.bin.ROC.Scenario <- terra::rast(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_ROCbin.tif"))
	  terra::writeRaster(Pred.bin.ROC.Scenario, paste0("Results/Global/Projections/",SpeciesName,".",Scenario.name,".bin.ROC.tif"), overwrite = TRUE)

	gc()	
	})
		
	#@@@@## Move biomod2 results to Results/Global/Models folder
	source_folder <- sp.name
	destination_folder <- "Results/Global/Models"
	if(file.exists(destination_folder)) {
	  unlink(destination_folder, recursive = TRUE)
	}

	file.rename(from = source_folder, to = destination_folder) #@@@@## The tif values are saved both in the biomod folder and the Projections folder, I would suggest deleting them from any of the places

  # Logs success or error messages
  message("NSH.SDM.Global.Model executed successfully")
  }, error = function(err) {
    message("Error in NSH.SDM.Global.Model:", conditionMessage(err))
    return(list(result = NULL, error = err))
  })
}


####################
# 4. REGIONAL MODELS
####################

NSH.SDM.Regional.Models <- function(nshsdm_selvars, 
				models=c("GAM","GBM", "RF", "MAXNET","GLM"),
				CV.nb.rep=1, 
				CV.perc=0.8) {
  
  nshsdm_global <- as.list(match.call())$nshsdm_selvars
  if(!inherits(nshsdm_selvars, "nshsdm.input")){
      stop("nshsdm_selvars must be an object of nshsdm.input class.")
  }

  SpeciesName <- nshsdm_selvars$Species.Name

  tryCatch({ 
	# Regional model calibrated with all the independent variables
	# Here we format the response (presence/background) and explanatory (environmental variables) data for BIOMOD2
	myResp.xy <- rbind(nshsdm_selvars$SpeciesData.XY.Regional ,nshsdm_selvars$Background.XY.Regional)
	row.names(myResp.xy)<-c(1:nrow(myResp.xy))
	myResp <- data.frame(c(rep(1,nrow(nshsdm_selvars$SpeciesData.XY.Regional)),rep(NA,nrow(nshsdm_selvars$Background.XY.Regional ))))
	names(myResp)<-"pa"
	row.names(myResp)<-c(1:nrow(myResp.xy))
	myExpl <- terra::extract(nshsdm_selvars$IndVar.Regional.2, myResp.xy, as.df=TRUE)[, -1]
      
	# Data required for the biomod2 package
	myBiomodData <- biomod2::BIOMOD_FormatingData(resp.var = myResp, 
                                           resp.xy = myResp.xy, 
                                           expl.var = myExpl, 
                                           resp.name = SpeciesName, 
                                           PA.nb.rep = 1, 
                                           PA.nb.absences = nrow(nshsdm_selvars$Background.XY.Regional), 
                                           PA.strategy = "random")
      
	# Calibrate and evaluate individual models with specified statistical algorithms
	myBiomodModelOut <- biomod2::BIOMOD_Modeling(bm.format = myBiomodData, 
						modeling.id = "AllModels", 
						# bm.options = myThecniquesOptions,
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

	# Replicates with ROC > 0.8	
	df <- myBiomodModelOut@models.evaluation
	df_slot <- slot(df, "val")
	df_slot <- df_slot[df_slot$metric.eval == "ROC", ]
	percentage <- 100 * sum(df_slot$validation >= 0.8)/nrow(df_slot)
	warning(sprintf("\n%.2f%% of replicates with AUC values >= 0.8.\n", percentage))			
      
	# Generate and evaluate a single ensemble (i.e.,consensus) model that averages the individual models weighted by
	# the value of the AUC statistic.
	myBiomodEM.ROC <- biomod2::BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
						models.chosen = 'all',
						em.by = 'all',
						em.algo = c("EMwmean"),
						metric.select = c('ROC'),
						metric.select.thresh = NULL,
						var.import = 3,
						metric.eval = c('ROC', "TSS", "KAPPA"),
						seed.val = 42) 
      
      
	# Project the individual models to the study area at regional scale under training conditions)
	myBiomodProj <- biomod2::BIOMOD_Projection(bm.mod = myBiomodModelOut,
					new.env = nshsdm_selvars$IndVar.Regional.2,
					proj.name = "Current",
					models.chosen = 'all',
					build.clamping.mask = FALSE) # Faster 
      
	# Project the ensemble model the study area at regional scale under training conditions
	biomod2::BIOMOD_EnsembleForecasting(bm.em = myBiomodEM.ROC, 
                                 bm.proj = myBiomodProj,
                                 models.chosen = 'all',
                                 metric.binary = 'all',
                                 metric.filter = 'all',
                                 build.clamping.mask = FALSE)
      
	# Load the model stored by biomod2 and save it in geotif format
	sp.name<-myBiomodData@sp.name
	Pred.ROC <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble.tif")) 
	terra::writeRaster(Pred.ROC, paste0("Results/Regional/Projections/",SpeciesName,".Current.tif"), overwrite=TRUE)
      
	# Binary models
	Pred.bin.ROC <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_ROCbin.tif"))
	terra::writeRaster(Pred.bin.ROC, paste0("Results/Regional/Projections/",SpeciesName,".Current.bin.ROC.tif"), overwrite=TRUE)
      
	Pred.bin.TSS <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_TSSbin.tif"))
	terra::writeRaster(Pred.bin.TSS, paste0("Results/Regional/Projections/",SpeciesName,".Current.bin.TSS.tif"), overwrite=TRUE)
      
   	# Save some results
	# Values of the statistics for each of the replicates   
	myEMeval.replicates <- biomod2::get_evaluations(myBiomodModelOut)
	write.csv(myEMeval.replicates,file=paste0("Results/Regional/Values/",SpeciesName,"_replica.csv"))
	# Values of the statistics of the consensus model   
	myEMeval.Consenso <- biomod2::get_evaluations(myBiomodEM.ROC)
	write.csv(myEMeval.Consenso,file=paste0("Results/Regional/Values/",SpeciesName,"_ensemble.csv"))
	# Variable importance
	myModelsVarImport <- biomod2::get_variables_importance(myBiomodModelOut)
	write.csv(myModelsVarImport, file = paste0("Results/Regional/Values/",SpeciesName,"_indvar.csv"), row.names = T)

	# Model projections for future climate scenarios
	################################################
	Scenarios <- dir_ls(paste0(VariablesPath,"/Regional"), pattern="tif")
        Scenarios <- Scenarios[!grepl("Current.tif", Scenarios)][1:2] # para probar con dos

	if(length(Scenarios) == 0) {
	  stop("There are no future climate variables")
	}
        
	walk(Scenarios, function(projmodel) { 
	  #projmodel <- Scenarios[1]
	  NewClim <- terra::rast(projmodel)[[nshsdm_selvars$Selected.Variables.Regional]]
	  Scenario.name <- path_file(projmodel) |> path_ext_remove()
        
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
	  terra::writeRaster(Pred.ROC.Scenario, paste0("Results/Regional/Projections/",SpeciesName,".",Scenario.name,".tif"), overwrite = TRUE)
        
	  # Binarized models
	  Pred.bin.TSS.Scenario <- terra::rast(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_ROCbin.tif"))
	  terra::writeRaster(Pred.bin.TSS.Scenario, paste0("Results/Regional/Projections/",SpeciesName,".",Scenario.name,".bin.TSS.tif"), overwrite = TRUE)
        
	  Pred.bin.ROC.Scenario <- terra::rast(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_ROCbin.tif"))
	  terra::writeRaster(Pred.bin.ROC.Scenario, paste0("Results/Regional/Projections/",SpeciesName,".",Scenario.name,".bin.ROC.tif"), overwrite = TRUE)

	  gc()	
	})

    #@@@@## Move biomod2 results to Results/Regional/Models folder
	source_folder <- sp.name
	destination_folder <- "Results/Regional/Models"
	if(file.exists(destination_folder)) {
	  unlink(destination_folder, recursive = TRUE)
	}
	file.rename(from = source_folder, to = destination_folder)
    
  # Logs success or error messages 
  message("NSH.SDM.Regional.Model executed successfully")
  }, error = function(err) {
    message("Error in NSH.SDM.Regional.Model:", conditionMessage(err))
    return(list(result = NULL, error = err))
  })
}
