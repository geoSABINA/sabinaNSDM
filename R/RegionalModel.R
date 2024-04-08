#' @export
#'
####################
# 4. REGIONAL MODELS
####################

NSH.SDM.Regional.Models <- function(nshsdm_selvars,
                                    algorithms=c( "GLM", "GAM", "RF"),
                                    CV.nb.rep=10,
                                    CV.perc=0.8,
                                    CustomModelOptions=NULL,
                                    save.output=TRUE,
                                    rm.biomod.folder=TRUE){

  #nshsdm_global <- as.list(match.call())$nshsdm_selvars
  if(!inherits(nshsdm_selvars, "nshsdm.input")){
      stop("nshsdm_selvars must be an object of nshsdm.input class.")
  }

  models <- toupper(algorithms)
  if(any(!models %in% c( "GAM", "GBM", "GLM", "MARS", "MAXNET", "RF"))) {
    stop("Please select at least one valid algorithm (\"GLM\", \"GAM\",\"MARS\",\"GBM\",\"MAXNET\", or \"RF\").")
  }

  SpeciesName <- nshsdm_selvars$Species.Name
  Level="Regional"

  nshsdm_data<-nshsdm_selvars[!names(nshsdm_selvars) %in% c("Summary", "args")]
  nshsdm_data$args <- list()
  nshsdm_data$args$algorithms <- algorithms
  nshsdm_data$args$CV.nb.rep <- CV.nb.rep
  nshsdm_data$args$CV.perc <- CV.perc

  current.projections <- list()
  new.projections <- list()
  new.projections$Pred.Scenario <- list()
  new.projections$Pred.bin.TSS.Scenario <- list()
  new.projections$Pred.bin.ROC.Scenario <- list()

  # Regional model calibrated with all the independent variables
  # Format the response (presence/background) and explanatory (environmental variables) data for BIOMOD2
  myResp.xy <- rbind(nshsdm_selvars$SpeciesData.XY.Regional ,nshsdm_selvars$Background.XY.Regional)
  row.names(myResp.xy)<-c(1:nrow(myResp.xy))
  myResp <- data.frame(c(rep(1,nrow(nshsdm_selvars$SpeciesData.XY.Regional)),rep(NA,nrow(nshsdm_selvars$Background.XY.Regional ))))
  names(myResp)<-"pa"
  row.names(myResp)<-c(1:nrow(myResp.xy))
  myExpl <- terra::extract(nshsdm_selvars$IndVar.Regional.Selected , myResp.xy, as.df=TRUE)[, -1]

  # Data required for the biomod2 package
  myBiomodData <- biomod2::BIOMOD_FormatingData(resp.var = myResp,
                                           resp.xy = myResp.xy,
                                           expl.var = myExpl,
                                           resp.name = SpeciesName,
                                           PA.nb.rep = 1,
                                           PA.nb.absences = nrow(nshsdm_selvars$Background.XY.Regional),
                                           PA.strategy = "random")

  # Calibrate and evaluate individual models with specified statistical algorithms
  myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
	                                      modeling.id = "AllModels",
	                                      models = models,
	                                      bm.options = CustomModelOptions, # Use the specified or default modeling options
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
  nreplicates<-sum(df_slot$validation >= 0.8)
  percentage <- 100 * nreplicates/nrow(df_slot)
  message(sprintf("\n%.2f%% of replicates with AUC values >= 0.8.\n", percentage))
  nreplicates<-data.frame(Algorithm="All",'Number of replicates'=nreplicates)
  for (algorithm.i in models) {
    nreplicates<-rbind(nreplicates,c(algorithm.i,sum(df_slot$validation[which(df_slot$algo==algorithm.i)] >= 0.8)))
  }
  
  nshsdm_data$nreplicates <- nreplicates

  # Generate and evaluate a single ensemble (i.e.,consensus) model that averages the individual models
  myBiomodEM.ROC <- biomod2::BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
						models.chosen = 'all',
						em.by = 'all',
						em.algo = c("EMmean"),
						metric.select = c('ROC'),
						metric.select.thresh = 0.8,
						var.import = 0, #@RGM esto lo he cambiado, creo que solo es necesario en el paso anterior
						metric.eval = c('ROC', "TSS", "KAPPA"),
						seed.val = 42)


  # Project the individual models to the study area at regional scale under training conditions)
  myBiomodProj <- biomod2::BIOMOD_Projection(bm.mod = myBiomodModelOut,
					new.env = nshsdm_selvars$IndVar.Regional.Selected ,
					proj.name = "Current",
					models.chosen = 'all',
					build.clamping.mask = FALSE)

  # Project the ensemble model the study area at regional scale under training conditions
  biomod2::BIOMOD_EnsembleForecasting(bm.em = myBiomodEM.ROC,
                                 bm.proj = myBiomodProj,
                                 models.chosen = 'all',
                                 metric.binary = 'all',
                                 metric.filter = 'all',
                                 build.clamping.mask = FALSE)

  # Load the model stored by biomod2 and save it in geotif format
  sp.name<-myBiomodData@sp.name
  Pred <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble.tif"))
  Pred<-terra::rast(wrap(Pred))
  
  nshsdm_data$current.projections$Pred <- setNames(Pred, paste0(SpeciesName, ".Current"))

  if(save.output){
    dir_create(paste0("Results/",Level,"/Projections/"))
    file_path<-paste0("Results/",Level,"/Projections/",SpeciesName,".Current.tif")
    terra::writeRaster(Pred,file_path , overwrite=TRUE)
    #message(paste("Projections at regional level under training conditions saved in:",file_path))
    file.remove(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble.tif"))
  }

  # Binary models
  Pred.bin.ROC <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_ROCbin.tif"))
  Pred.bin.TSS <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_TSSbin.tif"))
  Pred.bin.ROC<-terra::rast(wrap(Pred.bin.ROC))
  Pred.bin.TSS<-terra::rast(wrap(Pred.bin.TSS))
  
  nshsdm_data$current.projections$Pred.bin.ROC <- setNames(Pred.bin.ROC, paste0(SpeciesName, ".Current.bin.ROC"))
  nshsdm_data$current.projections$Pred.bin.TSS <- setNames(Pred.bin.TSS, paste0(SpeciesName, ".Current.bin.TSS"))

  if(save.output){
    file_path<-paste0("Results/",Level,"/Projections/",SpeciesName,".Current.bin.ROC.tif")
    terra::writeRaster(Pred.bin.ROC, file_path, overwrite=TRUE)
    #message(paste("ROC binary projections at regional level under training conditions saved in:",file_path))
    file.remove(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_ROCbin.tif"))

    file_path<-paste0("Results/",Level,"/Projections/",SpeciesName,".Current.bin.TSS.tif")
    terra::writeRaster(Pred.bin.TSS,file_path , overwrite=TRUE)
    #message(paste("TSS binary projections at regional level under training conditions saved in:",file_path))
    file.remove(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_TSSbin.tif"))
  }

  # Values of the statistics for each of the replicates
  myEMeval.replicates <- biomod2::get_evaluations(myBiomodModelOut)
  
  nshsdm_data$myEMeval.replicates <- myEMeval.replicates

  if(save.output){
    dir_create(paste0("Results/",Level,"/Values/",recurse=T))
    write.csv(myEMeval.replicates,file=paste0("Results/",Level,"/Values/",SpeciesName,"_replica.csv"))
    write.csv(nreplicates,file=paste0("Results/",Level,"/Values/",SpeciesName,"_nreplicates.csv"))
  }

  # Values of the statistics of the consensus model
  myEMeval.Ensemble <- biomod2::get_evaluations(myBiomodEM.ROC)
  
  nshsdm_data$myEMeval.Ensemble <- myEMeval.Ensemble
  
  if(save.output){
    write.csv(myEMeval.Ensemble,file=paste0("Results/",Level,"/Values/",SpeciesName,"_ensemble.csv"))
  }
	
  nshsdm_data$myEMeval.Ensemble <- myEMeval.Ensemble

  # Variable importance
  myModelsVarImport <- biomod2::get_variables_importance(myBiomodModelOut)

  if(save.output){
    write.csv(myModelsVarImport, file = paste0("Results/",Level,"/Values/",SpeciesName,"_indvar.csv"), row.names = T)
  }

  nshsdm_data$myModelsVarImport <- myModelsVarImport

  # Model projections for future climate scenarios
  ################################################
  Scenarios <- nshsdm_selvars$Scenarios
  
  if(length(Scenarios) == 0) {
    message("There are no new scenarios different from Current.tif!\n")
  } else {
    for(i in 1:length(Scenarios)) {
      new.env <- Scenarios[[i]][[nshsdm_selvars$Selected.Variables.Regional]]
      #new.env <- terra::mask(new.env, nshsdm_selvars$IndVar.Global.Selected[[1]])
      Scenario.name <- names(Scenarios[i])

      myBiomomodProjScenario <- biomod2::BIOMOD_Projection(bm.mod = myBiomodModelOut,
						new.env = new.env,
						proj.name = Scenario.name,
						models.chosen = "all")

      biomod2::BIOMOD_EnsembleForecasting(bm.em = myBiomodEM.ROC,
					bm.proj = myBiomomodProjScenario,
					models.chosen = "all",
					metric.binary = "all",
					metric.filter = "all")

      Pred.Scenario <- terra::rast(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble.tif"))
      Pred.Scenario<-terra::rast(wrap(Pred.Scenario))
      nshsdm_data$new.projections$Pred.Scenario[[i]] <- setNames(Pred.Scenario, paste0(SpeciesName,".",Scenario.name))

      if(save.output){
       file_path<-paste0("Results/",Level,"/Projections/",SpeciesName,".",Scenario.name,".tif")
       terra::writeRaster(Pred.Scenario, file_path, overwrite = TRUE)
       #message(paste("Projections at regional level under", Scenario.name,"conditions saved in:",file_path))
      file.remove(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble.tif"))
      }

      # Binarized models
      Pred.bin.ROC.Scenario <- terra::rast(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_ROCbin.tif"))
      Pred.bin.TSS.Scenario <- terra::rast(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_TSSbin.tif"))
      Pred.bin.ROC.Scenario<-terra::rast(wrap(Pred.bin.ROC.Scenario))
      Pred.bin.TSS.Scenario<-terra::rast(wrap(Pred.bin.TSS.Scenario))
    
      nshsdm_data$new.projections$Pred.bin.ROC.Scenario[[i]] <- setNames(Pred.bin.ROC.Scenario, paste0(SpeciesName,".",Scenario.name,".bin.ROC"))
      nshsdm_data$new.projections$Pred.bin.TSS.Scenario [[i]]<- setNames(Pred.bin.TSS.Scenario, paste0(SpeciesName,".",Scenario.name,".bin.TSS"))

      if(save.output){           #@@@# create the folders if they are not created yet
	file_path<-paste0("Results/",Level,"/Projections/",SpeciesName,".",Scenario.name,".bin.ROC.tif")
	terra::writeRaster(Pred.bin.ROC.Scenario, file_path, overwrite = TRUE)
	#message(paste("ROC binary projections at regiona level under", Scenario.name,"conditions saved in:",file_path))
	file.remove(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_ROCbin.tif"))

	file_path<-paste0("Results/",Level,"/Projections/",SpeciesName,".",Scenario.name,".bin.TSS.tif")
	terra::writeRaster(Pred.bin.TSS.Scenario,file_path , overwrite = TRUE)
	#message(paste("TSS binary projections at regional level under", Scenario.name,"conditions saved in:",file_path))
	file.remove(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_TSSbin.tif"))
      }
    } # end for
  } # end if(length(Scenarios) > 0 & all(match_vars))

  if(rm.biomod.folder || !save.output){
    # Remove species folder create by biomod2
    unlink(paste(sp.name,sep=""), recursive = TRUE) #@RGM no se estaba borrando
    
  } else {
    # Move biomod2 results to Results/Regional/Models folder
    dir_create(paste0("Results/",Level,"/Models/",sp.name))
    source_folder <- sp.name
    destination_folder <- paste0("Results/",Level,"/Models/",sp.name)
    if(file.exists(destination_folder)) {
      unlink(destination_folder, recursive = TRUE)}
      file.rename(from = source_folder, to = destination_folder)
      unlink(sp.name)
  }


  # Summary
  summary <- data.frame(Values = c("",
				#SpeciesName,
				paste(toupper(algorithms),collapse = ", "), 
				nrow(nshsdm_data$myEMeval.replicates), 
				myEMeval.Ensemble$calibration[which(myEMeval.Ensemble$metric.eval=="ROC")],
				myEMeval.Ensemble$calibration[which(myEMeval.Ensemble$metric.eval=="TSS")],
				myEMeval.Ensemble$calibration[which(myEMeval.Ensemble$metric.eval=="KAPPA")]))

  rownames(summary) <- c("  - Title 4:",
				#"Species name",
				"Statistical algorithms at regional level", 
				"Number of replicates with AUC > 0.8 at regional level", 
				"AUC of ensemble model at regional level", 
				"TSS of ensemble model at regional level",
				"KAPPA of ensemble model at regional level")

  summary <- rbind(nshsdm_selvars$Summary, summary) #@@@JMB guarda el summary acumulado
  
  nshsdm_data$Summary <- summary 
	#@@@# careful! This is not charging the summaries of global! #@@@JMB Yo no haria rbind con summaries anteriores

  #if(save.output){
  #  write.table(summary, paste0("Results/",SpeciesName,"_summary.csv"), sep=",",  row.names=F, col.names=T)
  #}

  attr(nshsdm_data, "class") <- "nshsdm.predict.r" #@@@JMB cambio atributo como propuesta para uso covariate()

  # Logs success or error messages
  #message("\nNSH.SDM.Regional.Models executed successfully!\n")

  if(save.output){
    message("Results saved in the following locations:")
    message(paste(
    " - Current and new projections: /Results/Regional/Projections/\n",
    "- Replicates statistics: /Results/Regional/Values/\n",
    "- Consensus model statistics: /Results/Regional/Values/\n",
    "- Variable importance: /Results/Regional/Values/"
    ))
    if(!rm.biomod.folder) { 
      message("- BIOMOD results: /Results/Regional/Models/\n")
    }
  }

  return(nshsdm_data)

}
