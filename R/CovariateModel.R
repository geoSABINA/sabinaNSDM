#' @export
NSH.SDM.Covariate.Models <- function(nshsdm_selvars, 
				models,
				CV.nb.rep=10, #@@@ he cambiado varias cosas aqui
				CV.perc=0.8,
				CustomModelOptions=NULL, #@@@ he cambiado varias cosas aqui,
				rm.biomod.folder=TRUE,
				save.output = TRUE) {

  #nshsdm_global <- as.list(match.call())$nshsdm_selvars
  if(!inherits(nshsdm_selvars, "nshsdm.input")){
      stop("nshsdm_selvars must be an object of nshsdm.input class.")
  }

  SpeciesName <- nshsdm_selvars$Species.Name

  nshsdm_data<-list()
  nshsdm_data$args <- list()
  nshsdm_data$args$models <- models
  nshsdm_data$args$CV.nb.rep <- CV.nb.rep
  nshsdm_data$args$CV.perc <- CV.perc

  current.projections <- list()
  new.projections <- list()
  new.projections$Pred.ROC.Scenario <- list()
  new.projections$Pred.bin.TSS.Scenario <- list()
  new.projections$Pred.bin.ROC.Scenario <- list()

  #tryCatch({ 
	# Covariate model calibrated with all the independent variables.
    	# Regional model excluding climatic variables.
	
	# Create directories to save results of the modeling process
	if(save.output){
      	dir_create(c("Results/Covariate/Values/", "Results/Covariate/Projections/"))
	}
      
	# Add the global model as an additional variable for the regional model.
	SDM.global <- terra::rast(paste0("Results/Global/Projections/",SpeciesName,".Current.tif")) #@RGM Cuidado si no guardamos los resultados en el disco duro del ordenador esto no funcionaria. Quizá habría que poner un aviso. 
	names(SDM.global) <- c("SDM.global") 
	IndVar.Regional.temp <- nshsdm_selvars$IndVar.Regional.Selected #@RGM cambiado 
	IndVar.Regional.Covariate <- c(IndVar.Regional.temp, SDM.global)
      
	myResp.xy <- rbind(nshsdm_selvars$SpeciesData.XY.Regional, nshsdm_selvars$Background.XY.Regional)
	names(myResp.xy)<-c("x","y")
	row.names(myResp.xy)<-c(1:nrow(myResp.xy))
	myResp <- data.frame(c(rep(1,nrow(nshsdm_selvars$SpeciesData.XY.Regional)),rep(NA,nrow(nshsdm_selvars$Background.XY.Regional ))))
	names(myResp)<-"pa"
	row.names(myResp)<-c(1:nrow(myResp.xy))
	myExpl <- terra::extract(IndVar.Regional.Covariate, myResp.xy, as.df=TRUE)[, -1]
      
	# Data required for the biomod2 package.
	myBiomodData <- biomod2::BIOMOD_FormatingData(resp.var = myResp, 
					resp.xy = myResp.xy, 
					expl.var = myExpl, 
					resp.name = SpeciesName, 
					PA.nb.rep = 1, 
					PA.nb.absences = nrow(nshsdm_selvars$Background.XY.Regional), 
					PA.strategy = "random")
      
	# Calibrate and evaluate individual models with specified statistical algorithms.
	if (is.null(CustomModelOptions)) {
	  # Use biomod2 default modeling options
	  myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,  #@@@## (This function saves outputs, check them and specify them in the description). Outputs: /NSHSDM/Larix.decidua/models/AllModels/Larix.decidua_PA1_RUN1_GBM
	                                      modeling.id = "AllModels", 
	                                      # bm.options = myThecniquesOptions,
	                                      models = models, 
	                                      CV.strategy = "random", 
	                                      CV.nb.rep = CV.nb.rep, CV.perc = CV.perc, 
	                                      weights = NULL, var.import = 3, #@@@## (var.import 3 or NULL) 
	                                      metric.eval = c("ROC", "TSS", "KAPPA", "ACCURACY", "SR", "BOYCE", "MPA"), 
	                                      scale.models = FALSE, do.progress = TRUE, 
	                                      prevalence = 0.5, seed.val = 42, 
	                                      CV.do.full.models = FALSE) # "CV.do.full.models = FALSE" and "var.import=0" to make it faster
	} else {
	  # Use custom modeling options provided by the user
	  # Model training using BIOMOD_Modeling
	  myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData, 
	                                      modeling.id = "AllModels", 
	                                      models = models, 
	                                      bm.options = CustomModelOptions, # Use the specified or default modeling options
	                                      CV.strategy = "random", 
	                                      CV.nb.rep = CV.nb.rep, CV.perc = CV.perc, 
	                                      weights = NULL, var.import = 3, 
	                                      metric.eval = c("ROC", "TSS", "KAPPA", "ACCURACY", "SR", "BOYCE", "MPA"), 
	                                      scale.models = FALSE, do.progress = TRUE, 
	                                      prevalence = 0.5, seed.val = 42, 
	                                      CV.do.full.models = FALSE)
	}

	# Replicates with ROC > 0.8	
	df <- myBiomodModelOut@models.evaluation
	df_slot <- slot(df, "val")
	df_slot <- df_slot[df_slot$metric.eval == "ROC", ]
	percentage <- 100 * sum(df_slot$validation >= 0.8)/nrow(df_slot)
	warning(sprintf("\n%.2f%% of replicates with AUC values >= 0.8.\n", percentage))
      
	# Generate and evaluate a single ensemble (i.e.,consensus) model that averages the individual models weighted
	# by the value of the AUC statistic.
	myBiomodEM.ROC  <- biomod2::BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                                 models.chosen = 'all',
                                                 em.by = 'all',
                                                 em.algo = c("EMmean"),#@RGM cambiado 
                                                 metric.select = c('ROC'),
                                                 metric.select.thresh = 0.8,
                                                 var.import = 0, #@RGM esto lo he cambiado
                                                 metric.eval = c('ROC', "TSS", "KAPPA"),
                                                 seed.val = 42) 
      
   	# Project the individual models to the study area at regional scale under training conditions.
	myBiomodProj <- biomod2::BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                        new.env = IndVar.Regional.Covariate, #@RGM cambiado 
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
	if(save.output){
 	terra::writeRaster(Pred.ROC, paste0("Results/Covariate/Projections/",SpeciesName,".Current.tif"), overwrite=TRUE)
	}
	nshsdm_data$current.projections$Pred.ROC <- setNames(Pred.ROC, paste0(SpeciesName, ".Current"))
      
	# Binary models
	Pred.bin.ROC <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_ROCbin.tif")) 
	if(save.output){
	terra::writeRaster(Pred.bin.ROC, paste0("Results/Covariate/Projections/",SpeciesName,".Current.bin.ROC.tif"), overwrite=TRUE)
	}
	nshsdm_data$current.projections$Pred.bin.ROC <- setNames(Pred.bin.ROC, paste0(SpeciesName, ".Current.bin.ROC"))
      
	Pred.bin.TSS <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_TSSbin.tif")) 
	if(save.output){
	terra::writeRaster(Pred.bin.TSS, paste0("Results/Covariate/Projections/",SpeciesName,".Current.bin.TSS.tif"), overwrite=TRUE)
	}
	nshsdm_data$current.projections$Pred.bin.TSS <- setNames(Pred.bin.TSS, paste0(SpeciesName,".Current.bin.TSS"))

        # Save some results
	# Values of the statistics for each of the replicates   
	myEMeval.replicates <- biomod2::get_evaluations(myBiomodModelOut)
	if(save.output){
	write.csv(myEMeval.replicates,file=paste0("Results/Covariate/Values/",SpeciesName,"_replica.csv"))
	}
	nshsdm_data$myEMeval.replicates <- myEMeval.replicates

	# Values of the statistics of the consensus model   
	myEMeval.Consenso <- biomod2::get_evaluations(myBiomodEM.ROC)
	if(save.output){
	write.csv(myEMeval.Consenso,file=paste0("Results/Covariate/Values/",SpeciesName,"_ensemble.csv"))
	}
	nshsdm_data$myEMeval.Consenso <- myEMeval.Consenso

	# Variable importance
	myModelsVarImport <- biomod2::get_variables_importance(myBiomodModelOut)
	if(save.output){
	write.table(myModelsVarImport, file = paste0("Results/Covariate/Values/",SpeciesName,"_indvar.csv"), row.names = T, col.names = T)
	}
	nshsdm_data$myModelsVarImport <- myModelsVarImport

      
	# Model projections for future climate scenarios
	################################################
	Scenarios <- dir_ls(paste0(VariablesPath,"/Regional"), pattern="tif")
        Scenarios <- Scenarios[!grepl("Current.tif", Scenarios)]

	if(length(Scenarios) == 0) {
	  stop("There are no future climate variables")
	}

	for(i in 1:length(Scenarios)) {   
	  projmodel <- Scenarios[i]
        #walk(Scenarios, function(projmodel) { 
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
	  if(save.output){
	  terra::writeRaster(Pred.ROC.Scenario, paste0("Results/Covariate/Projections/",SpeciesName,".",Scenario.name,".tif"), overwrite = TRUE)
          }
	  nshsdm_data$new.projections$Pred.ROC.Scenario[[i]] <- setNames(Pred.ROC.Scenario, paste0(SpeciesName,".",Scenario.name))

	  # Binarized models
          Pred.bin.TSS.Scenario <- terra::rast(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_ROCbin.tif"))
	  if(save.output){
          terra::writeRaster(Pred.bin.TSS.Scenario, paste0("Results/Covariate/Projections/",SpeciesName,".",Scenario.name,".bin.TSS.tif"), overwrite = TRUE)
	  }
	  nshsdm_data$new.projections$Pred.bin.TSS.Scenario[[i]] <- setNames(Pred.bin.TSS.Scenario, paste0(SpeciesName,".",Scenario.name,".bin.TSS"))

          Pred.bin.ROC.Scenario <- terra::rast(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_ROCbin.tif"))
	  if(save.output){
          terra::writeRaster(Pred.bin.ROC.Scenario, paste0("Results/Covariate/Projections/",SpeciesName,".",Scenario.name,".bin.ROC.tif"), overwrite = TRUE)
	  }
	  nshsdm_data$new.projections$Pred.bin.ROC.Scenario[[i]] <- setNames(Pred.bin.ROC.Scenario, paste0(SpeciesName,".",Scenario.name,".bin.ROC"))

	#}) #walk
	} #for

	if(rm.biomod.folder){
	# Remove species folder create by biomod2
	unlink(paste0(SpeciesName))
	} else {
	# Move biomod2 results to Results/Covariate/Models folder
  	dir_create(paste0("Results/Covariate/Models/",sp.name))
	source_folder <- sp.name
	destination_folder <- paste0("Results/Covariate/Models/",sp.name)
	if (file.exists(destination_folder)) {
	  unlink(destination_folder, recursive = TRUE)}
	file.rename(from = source_folder, to = destination_folder)
	nshsdm_data$links$biomod.folder <- destination_folder
	} 

  	gc()

  # Logs success or error messages 
  message("\nNSH.SDM.Covariate.Models executed successfully!\n")
  if(save.output){
  message("Results saved in the following locations:")
  message(paste(
    " - Current projections: /Results/Covariate/Projections/\n",
    "- ReplicateS statistics: /Results/Covariate/Values/\n",
    "- Consensus model statistics: /Results/Covariate/Values/\n",
    "- Variable importance: /Results/Covariate/Values/\n",
    "- New projections: /Results/Covariate/Projections/\n"
  ))
  }  	

  #}, error = function(err) {
  #  message("Error in NSH.SDM.Covariate.Models:", conditionMessage(err))
  #  return(list(result = NULL, error = err))
  #})	

  return(nshsdm_data)

}

