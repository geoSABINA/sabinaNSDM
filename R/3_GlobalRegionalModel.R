#' @export
NSH.SDM.Global.Model <- function(NSH.SDM.Data, NSH.SDM.Variables, ModelTechniques, CV.n.Repetitions=1, CV.Percentage=0.8,SpeciesName) #@@@##(CAREFUL! I added parameter speciesName, I also set defaults for CV.n.Repetitions and CV.Percentage parameters)
  {
  library(biomod2)
  tryCatch({ 
	# Here we format the response (presence/background) and explanatory (environmental variables) data for BIOMOD2
	myResp.xy <- rbind(NSH.SDM.Data$SpeciesData.XY.Global,NSH.SDM.Data$Background.XY.Global)
	names(myResp.xy)<-c("x","y")
	row.names(myResp.xy)<-c(1:nrow(myResp.xy))
	myResp <- data.frame(c(rep(1,nrow(NSH.SDM.Data$SpeciesData.XY.Global)),rep(NA,nrow(NSH.SDM.Data$Background.XY.Global))));names(myResp)<-"pa"
	row.names(myResp)<-c(1:nrow(myResp.xy)) 
	myExpl <- data.frame (terra::extract (NSH.SDM.Variables$IndVar.Global.2, myResp.xy)[, -1]) #@@@## before it was taking regional variables  to train the model, now its training with global variables
	
	# Remaining script sections involve executing #@@@## (biomod2 instead of BIOMOD2) BIOMOD2 modeling procedures,
	# including data formatting, model training, projection, evaluation, and visualization.
	# Each step is carefully structured to process and analyze species distribution data, #@@@## (REMOVE THIS LINE, redundant with previous)
	# generate models, evaluate their performance, and visualize the results #@@@## (REMOVE THIS LINE, redundant with previous)
	# Please refer to #@@@## (biomod2 instead of BIOMOD2) BIOMOD2 documentation for detailed explanation of these steps
	
	# Prepare data required to calibrate the model with biomod2 package (Result: A BIOMOD.formated.data object) #@@@## I changed this from: Data required for the biomod2 package 
	myBiomodData <- BIOMOD_FormatingData(resp.var = myResp, 
	                                     resp.xy = myResp.xy, 
	                                     expl.var = myExpl, 
	                                     resp.name = SpeciesName, 
	                                     PA.nb.rep = 1, 
	                                     PA.nb.absences = nrow(NSH.SDM.Data$Background.XY.Global), 
	                                     PA.strategy = "random")
	

	# Calibrate and evaluate individual models with specified statistical algorithms #@@@## (I changed this from "Ensamble models" to "Train individual models", their are not ensambled yet)
	myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,  #@@@## (This function saves outputs, check them and specify them in the description). Outputs: /NSHSDM/Larix.decidua/models/AllModels/Larix.decidua_PA1_RUN1_GBM
	                                    modeling.id = "AllModels", 
	                                   # bm.options = myThecniquesOptions,
	                                    models = ModelTechniques, 
	                                    CV.strategy = "random", 
	                                    CV.nb.rep = CV.n.Repetitions, CV.perc = CV.Percentage, 
	                                    weights = NULL, var.import = 3, #@@@## (var.import 3 or NULL) 
	                                    metric.eval = c("ROC", "TSS", "KAPPA", "ACCURACY", "SR", "BOYCE", "MPA"), 
	                                    scale.models = FALSE, do.progress = TRUE, 
	                                    prevalence = 0.5, seed.val = 42, 
	                                    CV.do.full.models = FALSE) # "CV.do.full.models = FALSE" and "var.import=0" to make it faster
	
	# Generate and evaluate a single ensemble (i.e.,consensus) model that averages the individual models weighted by the value of the AUC statistic, eliminating the replicas whose AUC value is less than 0.8 #@@@## (Instead of: "Generate a single consensus model that is the average of the replicas generated in the previous step weighted by the value of the area under the ROC curve statistic, eliminating the replicas whose AUC value is less than 0.8"
	myBiomodEM.ROC  <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
	                                                          models.chosen = 'all',
	                                                          em.by = 'all',
	                                                          em.algo = c("EMwmean"),
	                                                          metric.select = c('ROC'),
	                                                          metric.select.thresh = c(0.8),
	                                                          var.import = 3, #@@@## change this to NULL?
	                                                          metric.eval = c('ROC', "TSS", "KAPPA"),
	                                                          seed.val = 42) 
	
	
	# Project the 30 replicas in the study area #@@@## (which 30?? I would say: Project the individual models to the study area at regional scale under training conditions)
	myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
	                                 new.env = NSH.SDM.Variables$IndVar.Global.3, #@@@## changed from NSH.SDM.Variables$IndVar.Global.2 to use regional instead of global data for projecting
	                                 proj.name = "Current",
	                                 models.chosen = 'all',
	                                  build.clamping.mask = FALSE) # Faster 
	
	# Project the ensemble model to the study area at regional scale under training conditions #@@@## (Instead of: "Final consensus model")
	BIOMOD_EnsembleForecasting(bm.em = myBiomodEM.ROC, 
	                           bm.proj = myBiomodProj,
	                           models.chosen = 'all',
	                           metric.binary = 'all',
	                           metric.filter = 'all',
	                           build.clamping.mask = FALSE) # Faster
	
	# Load the model stored by biomod2 and save it in geotif format
	sp.name<-myBiomodData@sp.name
	Pred.ROC <- terra::rast (paste(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble.tif",sep="")) #@@@## I changed SpeciesName to sp.name because it FAILED if the name was not separated with a dot (I separated it with a _)
	terra::writeRaster (Pred.ROC, paste("Results/Global/Geotif/",SpeciesName,".Current.tif",sep=""), overwrite=TRUE) #@@@## Change the name of geotif folder
	
	# Binary models
	Pred.bin.ROC <- rast (paste(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_ROCbin.tif",sep="")) # terra #@@@## I changed SpeciesName to sp.nañe because it FAILED if the na,e was not separated with a dot (I separated it with a _)
	terra::writeRaster(Pred.bin.ROC, paste("Results/Global/Geotif/",SpeciesName,".Current.bin.ROC.tif", sep = "" ), overwrite=TRUE)
		  
	Pred.bin.TSS <- rast (paste(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_TSSbin.tif",sep="")) # terra #@@@## I changed SpeciesName to sp.nañe because it FAILED if the na,e was not separated with a dot (I separated it with a _)
	terra::writeRaster(Pred.bin.TSS, paste("Results/Global/Geotif/",SpeciesName,".Current.bin.TSS.tif", sep = "" ), overwrite=TRUE)
		  
		 
	# Save some results
		# Values of the #@@@## (evaluation) statistics for each of the replicates   
		myEMeval.replicates <- get_evaluations(myBiomodModelOut) #@@@## I changed replicas to replicates
		write.csv(myEMeval.replicates,file=paste("Results/Global/Values/",SpeciesName,"_replica.csv",sep=""))
		# Values of the #@@@## (evaluation) statistics of the consensus model   
		myEMeval.Ensemble <- get_evaluations(myBiomodEM.ROC) #@@@## I changed Consenso to Ensemble
		write.csv(myEMeval.Ensemble,file=paste("Results/Global/Values/",SpeciesName,"_ensemble.csv",sep=""))
		# Variable importance
		myModelsVarImport <- get_variables_importance(myBiomodModelOut)
		write.csv(myModelsVarImport, file = paste("Results/Global/Values/", SpeciesName, "_indvar.csv", sep = ""), row.names = T) #@@@##  I changed this so is stored also as csv
		      
		# Model projections for future climate scenarios
		################################################
		Scenarios.temp <- list.files(paste0(VariablesPath,"/Regional"),pattern = "tif")  #@@@## pattern="tif"
		Scenanios.dim <- length(Scenarios.temp)
		Scenarios.temp <- Scenarios.temp[-which(Scenarios.temp =="Current.tif")] #@@@## I changed this from Scenarios.temp[c(-1)]
		Scenario <- gsub(".tif","", Scenarios.temp, fixed=TRUE)
		
		# i <- 1
		if (Scenanios.dim == 1) {
		  stop("Only current projection, there are no future climate variables")
		} else {
		  
		for  (i in 1:length(Scenario)) {  

		  Scenario.name <- Scenario[i]
		  NewClim <- c(rast(paste0(VariablesPath,"/Regional/",Scenario[i],".tif")))[[NSH.SDM.Variables$Selected.Variables.Global]]
		  
		  #Project the individual models #@@@## I added this
		  myBiomomodProjScenario <- BIOMOD_Projection(
		    bm.mod = myBiomodModelOut,
		    new.env = NewClim,
		    proj.name = Scenario.name,
		    models.chosen = "all"
		  )
		  #Project the ensemble model #@@@## I added this
		  BIOMOD_EnsembleForecasting(
		    bm.em = myBiomodEM.ROC,
		    bm.proj = myBiomomodProjScenario,
		    models.chosen = "all",
		    metric.binary = "all",
		    metric.filter = "all"
		  )
		  #@@@@## All these tifs are saved both in the biomod folder and the Geotif folder, I would suggest deleting them from one of these places
		  Pred.ROC.Scenario <- rast(paste(sp.name, "/proj_", Scenario.name, "/proj_", Scenario.name, "_", sp.name, "_ensemble.tif", sep = "")) #@@@## I changed from SpeciesName to sp.name
		  terra::writeRaster(Pred.ROC.Scenario, paste("Results/Global/Geotif/", SpeciesName, ".", Scenario.name, ".tif", sep = ""), overwrite = TRUE)
		  
		  
		  # Binarized models
		  Pred.bin.TSS.Scenario <- rast(paste(sp.name, "/proj_", Scenario.name, "/proj_", Scenario.name, "_", sp.name, "_ensemble_ROCbin.tif", sep = "")) #@@@## I changed from SpeciesName to sp.name
		  terra::writeRaster(Pred.bin.TSS.Scenario, paste("Results/Global/Geotif/", SpeciesName, ".", Scenario.name, ".bin.TSS.tif", sep = ""), overwrite = TRUE)
		  
		  Pred.bin.ROC.Scenario <- rast(paste(sp.name, "/proj_", Scenario.name, "/proj_", Scenario.name, "_", sp.name, "_ensemble_ROCbin.tif", sep = "")) #@@@## I changed from SpeciesName to sp.name
		  terra::writeRaster(Pred.bin.ROC.Scenario, paste("Results/Global/Geotif/", SpeciesName, ".", Scenario.name, ".bin.ROC.tif", sep = ""), overwrite = TRUE)
		
		}
		
		save.image(paste("Results/Global/Images/", SpeciesName, ".RData", sep = ""))
		gc()
		}		
		
		#@@@@## Move biomod2 results to Results/Global/Models folder
		source_folder <- sp.name
		destination_folder <- "Results/Global/Models"
		if (file.exists(destination_folder)) {
		  unlink(destination_folder, recursive = TRUE)}
		file.rename(from = source_folder, to = destination_folder) #@@@@## The tif values are saved both in the biomod folder and the Geotif folder, I would suggest deleting them from any of the places

  # Logs success or error messages
		message("NSH.SDM.Global.Model executed successfully")
  }, error = function(err) {
    message("Error in NSH.SDM.Global.Model:", conditionMessage(err))
    return(list(result = NULL, error = err))
  })
}

## EXAMPLE 

NSH.SDM.Global.Model.Results <- NSH.SDM.Global.Model(NSH.SDM.Data, NSH.SDM.Variables, c("GAM","GBM", "RF", "MAXNET","GLM"), 1, 0.8,SpeciesName= "Larix_decidua")   #@@@@#NSH.SDM.Global.Model.Results <- NSH.SDM.Global.Model(NSH.SDM.Data, NSH.SDM.Variables, c("GBM", "RF", "GLM"), 1, 0.8,SpeciesName= "Larix_decidua") 

####################
# 4. REGIONAL MODELS
####################
NSH.SDM.Regional.Models <- function(NSH.SDM.Data, NSH.SDM.Variables, ModelTechniques, CV.n.Repetitions=1, CV.Percentage=0.8,SpeciesName) #@@@##(CAREFUL! I added parameter speciesName, and I removed ClimaticVariables , I also set defaults for CV.n.Repe
{
  library(biomod2)
  tryCatch({ 
    
    # Regional model calibrated with all the independent variables #@@@##(CAREFUL! I  removed ClimaticVariables)

    # Regional model excluding climatic variables 
    #if("FALSE" %in% ClimaticVariables) { #@@@##(CAREFUL! I  removed ClimaticVariables)
      # Here we format the response (presence/background) and explanatory (environmental variables) data for BIOMOD2
      myResp.xy <- rbind(NSH.SDM.Data$SpeciesData.XY.Regional ,NSH.SDM.Data$Background.XY.Regional)
      names(myResp.xy)<-c("x","y")
      row.names(myResp.xy)<-c(1:nrow(myResp.xy))
      myResp <- data.frame(c(rep(1,nrow(NSH.SDM.Data$SpeciesData.XY.Regional)),rep(NA,nrow(NSH.SDM.Data$Background.XY.Regional ))));names(myResp)<-"pa"
      row.names(myResp)<-c(1:nrow(myResp.xy))
      myExpl <- data.frame (terra::extract (NSH.SDM.Variables$IndVar.Regional.2, myResp.xy)[, -1]) #@@@##(I removed .noClimatic, now all results of regional with or without climatic data are stored as .Regional)
      
      # Data required for the biomod2 package
      myBiomodData <- BIOMOD_FormatingData(resp.var = myResp, 
                                           resp.xy = myResp.xy, 
                                           expl.var = myExpl, 
                                           resp.name = SpeciesName, 
                                           PA.nb.rep = 1, 
                                           PA.nb.absences = nrow(NSH.SDM.Data$Background.XY.Regional), 
                                           PA.strategy = "random")
      
      
    	# Calibrate and evaluate individual models with specified statistical algorithms #@@@## (I changed this from "Ensamble models" to "Train individual models", their are not ensambled yet)
      myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData, 
                                          modeling.id = "AllModels", 
                                        #  bm.options = myThecniquesOptions,
                                          models = ModelTechniques, 
                                          CV.strategy = "random", 
                                          CV.nb.rep = CV.n.Repetitions, CV.perc = CV.Percentage, 
                                          weights = NULL, var.import = 3, 
                                          metric.eval = c("ROC", "TSS", "KAPPA", "ACCURACY", "SR", "BOYCE", "MPA"), 
                                          scale.models = FALSE, do.progress = TRUE, 
                                          prevalence = 0.5, seed.val = 42, 
                                          CV.do.full.models = FALSE) # "CV.do.full.models = FALSE" and "var.import=0" to make it faster
      
      # Generate and evaluate a single ensemble (i.e.,consensus) model that averages the individual models weighted by the value of the AUC statistic, eliminating the replicas whose AUC value is less than 0.8 #@@@## (Instead of: "Generate a single consensus model that is the average of the replicas generated in the previous step weighted by the value of the area under the ROC curve statistic, eliminating the replicas whose AUC value is less than 0.8"
      myBiomodEM.ROC  <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                                 models.chosen = 'all',
                                                 em.by = 'all',
                                                 em.algo = c("EMwmean"),
                                                 metric.select = c('ROC'),
                                                 metric.select.thresh = c(0.8), #@@@## it gives me error because none of my models have an AUC higher than 0.8, change to smaller value? put it as parameter?
                                                 var.import = 3,
                                                 metric.eval = c('ROC', "TSS", "KAPPA"),
                                                 seed.val = 42) 
      
      
      # Project the 30 replicas in the study area #@@@## (which 30?? I would say: Project the individual models to the study area at regional scale under training conditions)
      myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                        new.env = NSH.SDM.Variables$IndVar.Regional.2, #@@@##(I removed .noClimatic, now all results of regional with or without climatic data are stored as .Regional)
                                        proj.name = "Current",
                                        models.chosen = 'all',
                                        build.clamping.mask = FALSE) # Faster 
      
      # Project the ensemble model the study area at regional scale under training conditions #@@@## (Instead of: "Final consensus model")
      BIOMOD_EnsembleForecasting(bm.em = myBiomodEM.ROC, 
                                 bm.proj = myBiomodProj,
                                 models.chosen = 'all',
                                 metric.binary = 'all',
                                 metric.filter = 'all',
                                 build.clamping.mask = FALSE) # Faster
      
      # Load the model stored by biomod2 and save it in geotif format
      sp.name<-myBiomodData@sp.name
      Pred.ROC <- rast (paste(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble.tif",sep="")) 
      terra::writeRaster (Pred.ROC, paste("Results/Regional/Geotif/",SpeciesName,".Current.tif",sep=""), overwrite=TRUE)
      
      # Binary models
      Pred.bin.ROC <- rast (paste(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_ROCbin.tif",sep="")) # terra
      terra::writeRaster(Pred.bin.ROC, paste("Results/Regional/Geotif/",SpeciesName,".Current.bin.ROC.tif", sep = "" ), overwrite=TRUE)
      
      Pred.bin.TSS <- rast (paste(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_TSSbin.tif",sep="")) # terra
      terra::writeRaster(Pred.bin.TSS, paste("Results/Regional/Geotif/",SpeciesName,".Current.bin.TSS.tif", sep = "" ), overwrite=TRUE)
      
      
      # Save some results
      # Values of the statistics for each of the replicates   
      myEMeval.replicas <- get_evaluations(myBiomodModelOut)
      write.csv(myEMeval.replicas,file=paste("Results/Regional/Values/",SpeciesName,"_replica.csv",sep=""))
      # Values of the statistics of the consensus model   
      myEMeval.Consenso <- get_evaluations(myBiomodEM.ROC)
      write.csv(myEMeval.Consenso,file=paste("Results/Regional/Values/",SpeciesName,"_ensemble.csv",sep=""))
      # Variable importance
      myModelsVarImport <- get_variables_importance(myBiomodModelOut)
      write.csv(myModelsVarImport, file = paste("Results/Regional/Values/", SpeciesName, "_indvar.csv", sep = ""), row.names = T) #@@@##  I changed this so is stored also as csv

      # Model projections for future climate scenarios
      ################################################
      Scenarios.temp <- list.files(paste0(VariablesPath,"/Regional"))
      Scenanios.dim <- length(Scenarios.temp)
      Scenarios.temp <- Scenarios.temp[-which(Scenarios.temp =="Current.tif")] #@@@## I changed this from Scenarios.temp[c(-1)]
      Scenario <- gsub(".tif","", Scenarios.temp, fixed=TRUE)
      
      # i <- 1
      if (Scenanios.dim == 1) {
        stop("Only current projection, there are no future climate variables")
      } else {
        
      for  (i in 1:length(Scenario)) {  
        Scenario.name <- Scenario[i]
        NewClim <- c(rast(paste0(VariablesPath,"/Regional/",Scenario[i],".tif")))[[NSH.SDM.Variables$Selected.Variables.Regional ]] #@@@##(I removed .noClimatic, now all results of regional with or without climatic data are stored as .Regional)
        
        myBiomomodProjScenario <- BIOMOD_Projection(
          bm.mod = myBiomodModelOut,
          new.env = NewClim,
          proj.name = Scenario.name,
          models.chosen = "all"
        )
        
        BIOMOD_EnsembleForecasting(
          bm.em = myBiomodEM.ROC,
          bm.proj = myBiomomodProjScenario,
          models.chosen = "all",
          metric.binary = "all",
          metric.filter = "all"
        )
        
        Pred.ROC.Scenario <- rast(paste(sp.name, "/proj_", Scenario.name, "/proj_", Scenario.name, "_", sp.name, "_ensemble.tif", sep = ""))
        terra::writeRaster(Pred.ROC.Scenario, paste("Results/Regional/Geotif/", SpeciesName, ".", Scenario.name, ".tif", sep = ""), overwrite = TRUE)
        
        
        # Binarized models
        Pred.bin.TSS.Scenario <- rast(paste(sp.name, "/proj_", Scenario.name, "/proj_", Scenario.name, "_", sp.name, "_ensemble_ROCbin.tif", sep = ""))
        terra::writeRaster(Pred.bin.TSS.Scenario, paste("Results/Regional/Geotif/", SpeciesName, ".", Scenario.name, ".bin.TSS.tif", sep = ""), overwrite = TRUE)
        
        Pred.bin.ROC.Scenario <- rast(paste(sp.name, "/proj_", Scenario.name, "/proj_", Scenario.name, "_", sp.name, "_ensemble_ROCbin.tif", sep = "")) # terra
        terra::writeRaster(Pred.bin.ROC.Scenario, paste("Results/Regional/Geotif/", SpeciesName, ".", Scenario.name, ".bin.ROC.tif", sep = ""), overwrite = TRUE)
        
      }
      
      save.image(paste("Results/Regional/Images/", SpeciesName, ".RData", sep = "")) #@@@@#set the scenario
      gc()
      }
    #} 
    #@@@@## Move biomod2 results to Results/Regional/Models folder
    source_folder <- sp.name
    destination_folder <- "Results/Regional/Models"
    if (file.exists(destination_folder)) {
      unlink(destination_folder, recursive = TRUE)}
    file.rename(from = source_folder, to = destination_folder)
    
    # Logs success or error messages 
    message("NSH.SDM.Regional.Model executed successfully")
  }, error = function(err) {
    message("Error in NSH.SDM.Regional.Model:", conditionMessage(err))
    return(list(result = NULL, error = err))
  })

}
