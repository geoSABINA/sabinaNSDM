#' @export
NSH.SDM.Covariate.Models <- function(NSH.SDM.Data, NSH.SDM.Variables, ModelTechniques, CV.n.Repetitions, CV.Percentage, SpeciesName) #@@@##(CAREFUL! I added parameter speciesName, and removed ClimaticVariables)
{
  library(biomod2)
  tryCatch({ 
    
    # Covariate model calibrated with all the independent variables #@@@##(I removed .noClimatic, now all results of regional with or without climatic data are stored as .Regional)
    
    # Regional model excluding climatic variables #@@@##(I removed .noClimatic, now all results of regional with or without climatic data are stored as .Regional)
    #if("FALSE" %in% ClimaticVariables) { #@@@##(I removed .noClimatic, now all results of regional with or without climatic data are stored as .Regional)
      # Create directories to save results of the modeling process
      dir_create(c("Results/Covariate/", "Results/Covariate/Values/", "Results/Covariate/Geotif/", "Results/Covariate/Images/"))
      
      #  Format the response (presence/background) and explanatory (environmental variables) data for BIOMOD2 #@@@## (REMOVEthis line)
      # Add the global model as an additional variable for the regional model#@@@## INSTEAD OF:Load Global model
      SDM.global <- rast(paste("Results/Global/Geotif/",SpeciesName,".Current.tif",sep=""))
      names(SDM.global) <- c("SDM.global") 
      IndVar.Regional.temp <- NSH.SDM.Variables$IndVar.Regional.2 #@@@##(I removed .noClimatic, now all results of regional with or without climatic data are stored as .Regional)
      IndVar.Regional.2 <- c(IndVar.Regional.temp, SDM.global)
      
      myResp.xy <- rbind(NSH.SDM.Data$SpeciesData.XY.Regional ,NSH.SDM.Data$Background.XY.Regional)
      names(myResp.xy)<-c("x","y")
      row.names(myResp.xy)<-c(1:nrow(myResp.xy))
      myResp <- data.frame(c(rep(1,nrow(NSH.SDM.Data$SpeciesData.XY.Regional)),rep(NA,nrow(NSH.SDM.Data$Background.XY.Regional ))));names(myResp)<-"pa"
      row.names(myResp)<-c(1:nrow(myResp.xy))
      myExpl <- data.frame (terra::extract (IndVar.Regional.2 , myResp.xy)[, -1])
      
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
                                                 metric.select.thresh = c(0.8),
                                                 var.import = 3,
                                                 metric.eval = c('ROC', "TSS", "KAPPA"),
                                                 seed.val = 42) 
      
      
      # Project the 30 replicas in the study area  #@@@## (which 30?? I would say: Project the individual models to the study area at regional scale under training conditions)
      myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                        new.env = IndVar.Regional.2,
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
      Pred.ROC <- rast (paste(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble.tif",sep="")) 
      terra::writeRaster (Pred.ROC, paste("Results/Covariate/Geotif/",SpeciesName,".Current.tif",sep=""), overwrite=TRUE)
      
      # Binary models
      Pred.bin.ROC <- rast (paste(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_ROCbin.tif",sep="")) 
      terra::writeRaster(Pred.bin.ROC, paste("Results/Covariate/Geotif/",SpeciesName,".Current.bin.ROC.tif", sep = "" ), overwrite=TRUE)
      
      Pred.bin.TSS <- rast (paste(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_TSSbin.tif",sep="")) 
      terra::writeRaster(Pred.bin.TSS, paste("Results/Covariate/Geotif/",SpeciesName,".Current.bin.TSS.tif", sep = "" ), overwrite=TRUE)
      
      
      # Save some results
      # Values of the statistics for each of the replicates   
      myEMeval.replicas <- get_evaluations(myBiomodModelOut)
      write.csv(myEMeval.replicas,file=paste("Results/Covariate/Values/",SpeciesName,"_replica.csv",sep=""))
      # Values of the statistics of the consensus model   
      myEMeval.Consenso <- get_evaluations(myBiomodEM.ROC)
      write.csv(myEMeval.Consenso,file=paste("Results/Covariate/Values/",SpeciesName,"_ensemble.csv",sep=""))
      # Variable importance
      myModelsVarImport <- get_variables_importance(myBiomodModelOut)
      write.table(myModelsVarImport, file = paste("Results/Covariate/Values/", SpeciesName, "_indvar.csv", sep = ""), row.names = T, col.names = T, sep = ",") #@@@##  I changed this so is stored also as csv
      
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
      
      # i <- 1
      for  (i in 1:length(Scenario)) {  
        Scenario.name <- Scenario[i]
        NewClim.temp <- c(rast(paste0(VariablesPath,"/Regional/",Scenario[i],".tif")))[[NSH.SDM.Variables$Selected.Variables.Regional  ]]#@@@##(I removed .noClimatic, now all results of regional with or without climatic data are stored as .Regional)
        
        SDM.global.future <- rast(paste("Results/Global/Geotif/",SpeciesName,".",Scenario.name,".tif",sep=""))
        names(SDM.global.future) <- c("SDM.global")
        NewClim <- c(NewClim.temp, SDM.global.future)
        
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
        terra::writeRaster(Pred.ROC.Scenario, paste("Results/Covariate/Geotif/", SpeciesName, ".", Scenario.name, ".tif", sep = ""), overwrite = TRUE)
        
        
        # Binarized models
        Pred.bin.TSS.Scenario <- rast(paste(sp.name, "/proj_", Scenario.name, "/proj_", Scenario.name, "_", sp.name, "_ensemble_ROCbin.tif", sep = ""))
        terra::writeRaster(Pred.bin.TSS.Scenario, paste("Results/Covariate/Geotif/", SpeciesName, ".", Scenario.name, ".bin.TSS.tif", sep = ""), overwrite = TRUE)
        
        Pred.bin.ROC.Scenario <- rast(paste(sp.name, "/proj_", Scenario.name, "/proj_", Scenario.name, "_", sp.name, "_ensemble_ROCbin.tif", sep = "")) # terra
        terra::writeRaster(Pred.bin.ROC.Scenario, paste("Results/Covariate/Geotif/", SpeciesName, ".", Scenario.name, ".bin.ROC.tif", sep = ""), overwrite = TRUE)
        
      }
      
      save.image(paste("Results/Covariate/Images/", SpeciesName, ".RData", sep = ""))
      gc()
      }  
      
      #@@@@## Move biomod2 results to Results/Covariate/Models folder
      source_folder <- sp.name
      destination_folder <- "Results/Covariate/Models"
      if (file.exists(destination_folder)) {
        unlink(destination_folder, recursive = TRUE)}
      file.rename(from = source_folder, to = destination_folder)
    
    # Logs success or error messages 
    message("NSH.SDM.Covariate.Models executed successfully")
  }, error = function(err) {
    message("Error in NSH.SDM.Covariate.Models:", conditionMessage(err))
    return(list(result = NULL, error = err))
  })
}
