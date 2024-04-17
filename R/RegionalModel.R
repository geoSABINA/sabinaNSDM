#' @export
NSDM.Regional <- function(nsdm_selvars,
			  algorithms=c( "GLM", "GAM", "RF"),
			  CV.nb.rep=10,
			  CV.perc=0.8,
			  CustomModelOptions=NULL,
			  metric.select.thresh = 0.8,
			  save.output=TRUE,
			  rm.biomod.folder=TRUE){

  if(!inherits(nsdm_selvars, "nsdm.vinput")){
      stop("nsdm_selvars must be an object of nsdm.vinput class. Consider running NSDM.SelectCovariables() function.")
  }

  models <- toupper(algorithms)
  if(any(!models %in% c( "GAM", "GBM", "GLM", "MARS", "MAXNET", "RF"))) {
    stop("Please select at least one valid algorithm (\"GLM\", \"GAM\",\"MARS\",\"GBM\",\"MAXNET\", or \"RF\").")
  }

  SpeciesName <- nsdm_selvars$Species.Name
  Level="Regional"

  sabina<-nsdm_selvars[!names(nsdm_selvars) %in% c("Summary", "args")]
  sabina$args <- list()
  sabina$args$algorithms <- algorithms
  sabina$args$CV.nb.rep <- CV.nb.rep
  sabina$args$CV.perc <- CV.perc
  sabina$args$metric.select.thresh <- metric.select.thresh

  current.projections <- list()
  new.projections <- list()
  new.projections$Pred.Scenario <- list()
  new.projections$Pred.bin.TSS.Scenario <- list()
  new.projections$Pred.bin.ROC.Scenario <- list()

  # Unwrap objects
  IndVar.Regional.Selected <- terra::unwrap(nsdm_selvars$IndVar.Regional.Selected)
  if(!is.null(nsdm_selvars$Scenarios)) {
  Scenarios <- lapply(nsdm_selvars$Scenarios, terra::unwrap)
  }

  # Regional model calibrated with all the independent variables
  # Format the response (presence/background) and explanatory (environmental variables) data for BIOMOD2
  myResp.xy <- rbind(nsdm_selvars$SpeciesData.XY.Regional ,nsdm_selvars$Background.XY.Regional)
  row.names(myResp.xy)<-c(1:nrow(myResp.xy))
  myResp <- data.frame(c(rep(1,nrow(nsdm_selvars$SpeciesData.XY.Regional)),rep(NA,nrow(nsdm_selvars$Background.XY.Regional ))))
  names(myResp)<-"pa"
  row.names(myResp)<-c(1:nrow(myResp.xy))
  myExpl <- terra::extract(IndVar.Regional.Selected , myResp.xy, as.df=TRUE)[, -1]

  # Data required for the biomod2 package
  myBiomodData <- biomod2::BIOMOD_FormatingData(resp.var = myResp,
                                           resp.xy = myResp.xy,
                                           expl.var = myExpl,
                                           resp.name = SpeciesName,
                                           PA.nb.rep = 1,
                                           PA.nb.absences = nrow(nsdm_selvars$Background.XY.Regional),
                                           PA.strategy = "random")

  # Calibrate and evaluate individual models with specified statistical algorithms
  myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
	                                      modeling.id = "AllModels",
	                                      models = models,
	                                      OPT.user = CustomModelOptions, # Use the specified or default modeling options
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
  nreplicates<-sum(df_slot$validation >= metric.select.thresh)
  if(nreplicates == 0) {
    stop(paste0("\nNo replica for ", SpeciesName, " has reached an AUC value >= ", metric.select.thresh, ".\n"))
  }
  percentage <- 100 * nreplicates/nrow(df_slot)
  nreplicates<-data.frame(Algorithm="All",'Number of replicates'=nreplicates)
  for(algorithm.i in models) {
    nreplicates<-rbind(nreplicates,c(algorithm.i,sum(df_slot$validation[which(df_slot$algo==algorithm.i)] >= metric.select.thresh)))
  }
  
  sabina$nbestreplicates <- nreplicates #@@@JMB sugiero poner nbestreplicates en lugar de nreplicates.

  # Generate and evaluate a single ensemble (i.e.,consensus) model that averages the individual models
  myBiomodEM.ROC <- biomod2::BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
						models.chosen = 'all',
						em.by = 'all',
						em.algo = c("EMmean"),
						metric.select = c('ROC'),
						metric.select.thresh = metric.select.thresh,
						var.import = 0, #@RGM esto lo he cambiado, creo que solo es necesario en el paso anterior
						metric.eval = c('ROC', "TSS", "KAPPA"),
						seed.val = 42)


  # Project the individual models to the study area at regional scale under training conditions)
  myBiomodProj <- biomod2::BIOMOD_Projection(bm.mod = myBiomodModelOut,
					new.env = IndVar.Regional.Selected,
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
  
  sabina$current.projections$Pred <- setNames(Pred, paste0(SpeciesName, ".Current"))

  if(save.output){
    fs::dir_create(paste0("Results/",Level,"/Projections/"))
    file_path<-paste0("Results/",Level,"/Projections/",SpeciesName,".Current.tif")
    terra::writeRaster(Pred,file_path , overwrite=TRUE)
    #message(paste("Projections at regional level under training conditions saved in:",file_path))
    fs::file_delete(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble.tif"))
  }

  # Binary models
  Pred.bin.ROC <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_ROCbin.tif"))
  Pred.bin.TSS <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_TSSbin.tif"))
  Pred.bin.ROC<-terra::rast(wrap(Pred.bin.ROC))
  Pred.bin.TSS<-terra::rast(wrap(Pred.bin.TSS))
  
  sabina$current.projections$Pred.bin.ROC <- setNames(Pred.bin.ROC, paste0(SpeciesName, ".Current.bin.ROC"))
  sabina$current.projections$Pred.bin.TSS <- setNames(Pred.bin.TSS, paste0(SpeciesName, ".Current.bin.TSS"))

  if(save.output){
    file_path<-paste0("Results/",Level,"/Projections/",SpeciesName,".Current.bin.ROC.tif")
    terra::writeRaster(Pred.bin.ROC, file_path, overwrite=TRUE)
    #message(paste("ROC binary projections at regional level under training conditions saved in:",file_path))
    fs::file_delete(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_ROCbin.tif"))

    file_path<-paste0("Results/",Level,"/Projections/",SpeciesName,".Current.bin.TSS.tif")
    terra::writeRaster(Pred.bin.TSS,file_path , overwrite=TRUE)
    #message(paste("TSS binary projections at regional level under training conditions saved in:",file_path))
    fs::file_delete(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_TSSbin.tif"))
  }

  # Values of the statistics for each of the replicates
  myEMeval.replicates <- biomod2::get_evaluations(myBiomodModelOut)
  
  sabina$myEMeval.replicates <- myEMeval.replicates

  if(save.output){
    fs::dir_create(paste0("Results/",Level,"/Values/"))
    write.csv(myEMeval.replicates,file=paste0("Results/",Level,"/Values/",SpeciesName,"_replica.csv"))
    write.csv(nreplicates,file=paste0("Results/",Level,"/Values/",SpeciesName,"_nreplicates.csv"))
  }

  # Values of the statistics of the consensus model
  myEMeval.Ensemble <- biomod2::get_evaluations(myBiomodEM.ROC)
  
  sabina$myEMeval.Ensemble <- myEMeval.Ensemble
  
  if(save.output){
    write.csv(myEMeval.Ensemble,file=paste0("Results/",Level,"/Values/",SpeciesName,"_ensemble.csv"))
  }
	
  sabina$myEMeval.Ensemble <- myEMeval.Ensemble

  # Variable importance
  myModelsVarImport <- biomod2::get_variables_importance(myBiomodModelOut)

  if(save.output){
    write.csv(myModelsVarImport, file = paste0("Results/",Level,"/Values/",SpeciesName,"_indvar.csv"), row.names = T)
  }

  sabina$myModelsVarImport <- myModelsVarImport

  # Model projections for future climate scenarios
  ################################################
  if(length(Scenarios) == 0) {
    warning("No new scenarios for further projections!\n") #Aquí pondría un warning en lugar de message
  } else {
    for(i in 1:length(Scenarios)) {
      new.env <- Scenarios[[i]][[nsdm_selvars$Selected.Variables.Regional]]
      #new.env <- terra::mask(new.env, IndVar.Global.Selected[[1]])
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
      Pred.Scenario <- terra::rast(wrap(Pred.Scenario))

      sabina$new.projections$Pred.Scenario[[i]] <- setNames(Pred.Scenario, paste0(SpeciesName,".",Scenario.name))

      if(save.output){
       fs::dir_create(paste0("Results/",Level,"/Projections/"))
       file_path<-paste0("Results/",Level,"/Projections/",SpeciesName,".",Scenario.name,".tif")
       terra::writeRaster(Pred.Scenario, file_path, overwrite = TRUE)
       #message(paste("Projections at regional level under", Scenario.name,"conditions saved in:",file_path))
       fs::file_delete(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble.tif"))
      }

      # Binarized models
      Pred.bin.ROC.Scenario <- terra::rast(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_ROCbin.tif"))
      Pred.bin.TSS.Scenario <- terra::rast(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_TSSbin.tif"))
      Pred.bin.ROC.Scenario<-terra::rast(wrap(Pred.bin.ROC.Scenario))
      Pred.bin.TSS.Scenario<-terra::rast(wrap(Pred.bin.TSS.Scenario))
    
      sabina$new.projections$Pred.bin.ROC.Scenario[[i]] <- setNames(Pred.bin.ROC.Scenario, paste0(SpeciesName,".",Scenario.name,".bin.ROC"))
      sabina$new.projections$Pred.bin.TSS.Scenario [[i]]<- setNames(Pred.bin.TSS.Scenario, paste0(SpeciesName,".",Scenario.name,".bin.TSS"))

      if(save.output){
	file_path<-paste0("Results/",Level,"/Projections/",SpeciesName,".",Scenario.name,".bin.ROC.tif")
	terra::writeRaster(Pred.bin.ROC.Scenario, file_path, overwrite = TRUE)
	#message(paste("ROC binary projections at regiona level under", Scenario.name,"conditions saved in:",file_path))
	fs::file_delete(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_ROCbin.tif"))

	file_path<-paste0("Results/",Level,"/Projections/",SpeciesName,".",Scenario.name,".bin.TSS.tif")
	terra::writeRaster(Pred.bin.TSS.Scenario,file_path , overwrite = TRUE)
	#message(paste("TSS binary projections at regional level under", Scenario.name,"conditions saved in:",file_path))
	fs::file_delete(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_TSSbin.tif"))
      }
    } # end for
  } # end if(length(Scenarios) > 0 & all(match_vars))

  source_folder <- sp.name
  destination_folder <- paste0("Results/",Level,"/Models/",sp.name)

  if(rm.biomod.folder){
    # Remove species folder created by biomod2
    fs::dir_delete(source_folder)
  } else { 
    if(save.output){
      # Move and remove biomod2 results from /sp.name/ to Results/Regional/Models/ folder
      fs::dir_create(paste0("Results/",Level,"/Models/",sp.name))
      fs::dir_copy(source_folder, destination_folder, overwrite = TRUE)
      fs::dir_delete(source_folder)
    }
  }


  # Summary
  summary <- data.frame(Values = c(SpeciesName,
				paste(toupper(algorithms),collapse = ", "), 
				sum(sabina$myEMeval.replicates$metric.eval == "ROC" & sabina$myEMeval.replicates$validation >= metric.select.thresh), 
				myEMeval.Ensemble$calibration[which(myEMeval.Ensemble$metric.eval=="ROC")],
				myEMeval.Ensemble$calibration[which(myEMeval.Ensemble$metric.eval=="TSS")],
				myEMeval.Ensemble$calibration[which(myEMeval.Ensemble$metric.eval=="KAPPA")]))

  rownames(summary) <- c("Species name",
				"Statistical algorithms at regional level", 
				paste0("Number of replicates with AUC > ",metric.select.thresh, " at regional level"), 
				"AUC of ensemble model at regional level", 
				"TSS of ensemble model at regional level",
				"KAPPA of ensemble model at regional level")

  # Wrap objects
  sabina$current.projections <- rapply(sabina$current.projections, terra::wrap, how = "list")
  if(!is.null(nsdm_selvars$Scenarios)) {
    sabina$new.projections <- rapply(sabina$new.projections, terra::wrap, how = "list")
  }
  
  sabina$Summary <- summary 

  attr(sabina, "class") <- "nsdm.predict.r"

  # % best replicates messages
  message(sprintf("\n%.2f%% of %s replicates with AUC values >= %.2f.\n", percentage, SpeciesName, metric.select.thresh))
  # save.out messages
  if(save.output){
    message("Results saved in the following local folder/s:")
    message(paste(
    " - Current and new projections: /Results/Regional/Projections/\n",
    "- Replicates statistics: /Results/Regional/Values/\n",
    "- Consensus model statistics: /Results/Regional/Values/\n",
    "- Variable importance: /Results/Regional/Values/"
    ))
    if(!rm.biomod.folder) { 
      message(" - BIOMOD results: /Results/Regional/Models/\n")
    }
  }

  return(sabina)

}

