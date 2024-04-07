#' @export
#'
####################
# 3. GLOBAL MODELS
####################

NSH.SDM.Global.Model <- function(nshsdm_selvars,
				algorithms=c("GLM", "GAM", "RF"),
				CV.nb.rep=10,
				CV.perc=0.8,
				CustomModelOptions=NULL,
				save.output=TRUE,
				rm.biomod.folder=TRUE) {

  #nshsdm_global <- as.list(match.call())$nshsdm_selvars
  if(!inherits(nshsdm_selvars, "nshsdm.input")){
      stop("nshsdm_selvars must be an object of nshsdm.input class. Consider running VarSelection() function.")
  }

  models <- toupper(algorithms)
  if(any(!models %in% c( "GAM", "GBM", "GLM", "MARS", "MAXNET", "RF"))) {
    stop("Please select at least one valid algorithm (\"GLM\", \"GAM\", \"MARS\", \"GBM\", \"MAXNET\" or \"RF\").")
  }

  SpeciesName <- nshsdm_selvars$Species.Name
  Level="Global"

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

  # GLOBAL SCALE
  # Format the response (presence/background) and explanatory (environmental variables) data for BIOMOD2
  myResp.xy <- rbind(nshsdm_selvars$SpeciesData.XY.Global,nshsdm_selvars$Background.XY.Global)
  row.names(myResp.xy)<-c(1:nrow(myResp.xy))
  myResp <- data.frame(c(rep(1,nrow(nshsdm_selvars$SpeciesData.XY.Global)),rep(NA,nrow(nshsdm_selvars$Background.XY.Global))))
  names(myResp)<-"pa"
  row.names(myResp)<-c(1:nrow(myResp.xy))
  myExpl <- terra::extract(nshsdm_selvars$IndVar.Global.Selected , myResp.xy, as.df=TRUE)[, -1]  #@@@#TG antes usaba misma resolucion para entrenar y proyectar el global

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
  # Train and evaluate individual models using BIOMOD_Modeling
  myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,   #@@@## (This function saves outputs, check them and specify them in the description). Outputs: /NSHSDM/Larix.decidua/models/AllModels/Larix.decidua_PA1_RUN1_GBM
	                                      modeling.id = "AllModels",
	                                      models = models,
	                                      bm.options = CustomModelOptions,
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
  message(sprintf("\n%.2f%% of replicates with AUC values >= 0.8.\n", percentage)) #@@@JMB con message se qeuda perdido por la consola. Pendiente Bajarlo?
  nreplicates<-data.frame(Algorithm="All",'Number of replicates'=nreplicates) #@@@JMB nrepg0.8 en lugar de nreplicates?
  for(algorithm.i in models) {
    nreplicates<-rbind(nreplicates,c(algorithm.i,sum(df_slot$validation[which(df_slot$algo==algorithm.i)] >= 0.8)))
  }

  nshsdm_data$nreplicates <- nreplicates

  # Generate and evaluate a single ensemble (i.e.,consensus) model that averages the individual models
  myBiomodEM.ROC  <- biomod2::BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
					models.chosen = 'all',
					em.by = 'all',
					em.algo = c("EMmean"),
					metric.select = c('ROC'),
					metric.select.thresh = 0.8,
					var.import = 0, 	#@RGM esto lo he cambiado, creo que solo es necesario en el paso anterior
					metric.eval = c('ROC', "TSS", "KAPPA"),
					seed.val = 42)


  # Project the individual models to the study area at regional scale under training conditions
  myBiomodProj <- biomod2::BIOMOD_Projection(bm.mod = myBiomodModelOut,
					new.env = nshsdm_selvars$IndVar.Global.Selected.reg,  #@@@TG Lo he cambiado porque estabamos usando la misma resolucion para entrenar y para proyectar
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
  Pred <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble.tif")) #remove duplicated files
  Pred<-terra::rast(wrap(Pred))

  nshsdm_data$current.projections$Pred <- c(setNames(Pred, paste0(SpeciesName, ".Current")))

  if(save.output){
    dir_create(paste0("Results/",Level,"/Projections/")) #@@@JMB Pendiente revisar. Crea una folder en Projections que se llama TRUE. Ver en todas f()
    file_path <- paste0("Results/",Level,"/Projections/",SpeciesName,".Current.tif")
    terra::writeRaster(Pred, file_path, overwrite=TRUE)
    #message(paste("Projections at global level under training conditions saved in:",file_path)) #@@@JMB todos estos message están resumidos también al final. Decidir unos u otros.
    file.remove(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble.tif"))
  }

  # Binary models
  Pred.bin.ROC <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_ROCbin.tif"))
  Pred.bin.TSS <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_TSSbin.tif"))
  Pred.bin.ROC<-terra::rast(wrap(Pred.bin.ROC))
  Pred.bin.TSS<-terra::rast(wrap(Pred.bin.TSS))

  nshsdm_data$current.projections$Pred.bin.ROC <- setNames(Pred.bin.ROC, paste0(SpeciesName, ".Current.bin.ROC"))
  nshsdm_data$current.projections$Pred.bin.TSS <- setNames(Pred.bin.TSS, paste0(SpeciesName,".Current.bin.TSS"))

  if(save.output){
    file_path <- paste0("Results/",Level,"/Projections/",SpeciesName,".Current.bin.ROC.tif")
    terra::writeRaster(Pred.bin.ROC, file_path, overwrite=TRUE)
    #message(paste("ROC binary projections at global level under training conditions saved in:",file_path))
    file.remove(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_ROCbin.tif"))

    file_path <- paste0("Results/",Level,"/Projections/",SpeciesName,".Current.bin.TSS.tif")
    terra::writeRaster(Pred.bin.TSS, file_path, overwrite=TRUE)
    #message(paste("TSS binary projections at global level under training conditions saved in:",file_path))
    file.remove(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_TSSbin.tif"))
  }

  # Values of the evaluation statistics for each replica
  myEMeval.replicates <- biomod2::get_evaluations(myBiomodModelOut)
  
  nshsdm_data$myEMeval.replicates <- myEMeval.replicates

  if(save.output){
    dir_create(paste0("Results/",Level,"/Values/",recurse=T))
    file_path <- paste0("Results/",Level,"/Values/",SpeciesName,"_replica.csv")
    write.csv(myEMeval.replicates,file=file_path)
    write.csv(nreplicates,file=paste0("Results/",Level,"/Values/",SpeciesName,"_nreplicates.csv"))
  }

  # Values of the evaluation statistics of the consensus model
  myEMeval.Ensemble <- biomod2::get_evaluations(myBiomodEM.ROC)

  nshsdm_data$myEMeval.Ensemble <- myEMeval.Ensemble

  if(save.output){
    file_path <- paste0("Results/",Level,"/Values/",SpeciesName,"_ensemble.csv")
    write.csv(myEMeval.Ensemble,file=file_path)
  }

  # Variable importance
  myModelsVarImport <- biomod2::get_variables_importance(myBiomodModelOut)

  nshsdm_data$myModelsVarImport <- myModelsVarImport

  if(save.output){
    file_path <- paste0("Results/",Level,"/Values/",SpeciesName,"_indvar.csv")
    write.csv(myModelsVarImport, file = file_path, row.names = T)
  }


  # Model projections for future climate scenarios
  ################################################
  Scenarios <- nshsdm_selvars$Scenarios
  
  if(length(Scenarios) == 0) {
    message("There are no new scenarios different from Current.tif!\n")
  } else {
    for(i in 1:length(Scenarios)) {
      new.env <- Scenarios[[i]][[nshsdm_selvars$Selected.Variables.Global]]
      #new.env <- terra::mask(new.env, nshsdm_selvars$IndVar.Global.Selected[[1]])
      Scenario.name <- names(Scenarios[i])

      #Project the individual models
      myBiomomodProjScenario <- biomod2::BIOMOD_Projection(bm.mod = myBiomodModelOut,
							 new.env = new.env,
							 proj.name = Scenario.name,
							 models.chosen = "all")

      #Project the ensemble model
      biomod2::BIOMOD_EnsembleForecasting(bm.em = myBiomodEM.ROC,
					bm.proj = myBiomomodProjScenario,
					models.chosen = "all",
					metric.binary = "all",
					metric.filter = "all")
		  #@@@@## All these tifs are saved both in the biomod folder and the Projections folder, I would suggest deleting them from one of these places
      Pred.Scenario <- terra::rast(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble.tif"))
      Pred.Scenario <- terra::rast(wrap(Pred.Scenario))
	    
      nshsdm_data$new.projections$Pred.Scenario[[i]] <- setNames(Pred.Scenario, paste0(SpeciesName,".",Scenario.name))

      if(save.output){
        dir_create(paste0("Results/",Level,"/Projections/",recurse=T))
        file_path <- paste0("Results/",Level,"/Projections/",SpeciesName,".",Scenario.name,".tif")
        terra::writeRaster(Pred.Scenario, file_path, overwrite = TRUE)
        file.remove(paste0(sp.name,"/proj_",Scenario.name,"/proj_", Scenario.name,"_",sp.name,"_ensemble.tif"))
        #message(paste("Projections at global level under", Scenario.name,"conditions saved in:",file_path))
      }

      # Binarized models
      Pred.bin.ROC.Scenario <- terra::rast(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_ROCbin.tif"))
      Pred.bin.TSS.Scenario <- terra::rast(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_TSSbin.tif"))
      Pred.bin.ROC.Scenario<-terra::rast(wrap(Pred.bin.ROC.Scenario))
      Pred.bin.TSS.Scenario<-terra::rast(wrap(Pred.bin.TSS.Scenario))
      
      nshsdm_data$new.projections$Pred.bin.ROC.Scenario[[i]] <- setNames(Pred.bin.ROC.Scenario, paste0(SpeciesName,".",Scenario.name,".bin.ROC"))
      nshsdm_data$new.projections$Pred.bin.TSS.Scenario[[i]] <- setNames(Pred.bin.TSS.Scenario, paste0(SpeciesName,".",Scenario.name,".bin.TSS"))

      if(save.output){
        file_path <- paste0("Results/",Level,"/Projections/",SpeciesName,".",Scenario.name,".bin.ROC.tif")
        terra::writeRaster(Pred.bin.ROC.Scenario, file_path, overwrite = TRUE)
        #message(paste("ROC binary projections at global level under", Scenario.name,"conditions saved in:",file_path))
        file.remove(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_ROCbin.tif"))

        file_path <- paste0("Results/",Level,"/Projections/",SpeciesName,".",Scenario.name,".bin.TSS.tif")
        terra::writeRaster(Pred.bin.TSS.Scenario, file_path, overwrite = TRUE)
        #message(paste("TSS binary projections at global level under", Scenario.name,"conditions saved in:",file_path))
        file.remove(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name, "_ensemble_TSSbin.tif"))
      }
    } # end for
  } # end if(length(Scenarios) == 0)

  if(rm.biomod.folder || !save.output){ #rm= T save=F
    # Remove species folder created by biomod2
    unlink(sp.name)
  } else {
    # Move biomod2 results to Results/Global/Models folder
    dir_create(paste0("Results/",Level,"/Models/",sp.name))
    source_folder <- sp.name
    destination_folder <- paste0("Results/",Level,"/Models/",sp.name)
    if (file.exists(destination_folder)) {
      unlink(source_folder, recursive = TRUE)}  #@@@JMB igual hay que quitar el recursive=T para que no se borren folders anteriores. Probar
      file.rename(from = source_folder, to = destination_folder)
      unlink(sp.name)
  }

  gc()

  # Summary
  summary <- data.frame(Values = c("",
				#SpeciesName,
				paste(toupper(algorithms),collapse = ", "), 
				nrow(nshsdm_data$myEMeval.replicates), 
				myEMeval.Ensemble$calibration[which(myEMeval.Ensemble$metric.eval=="ROC")],
				myEMeval.Ensemble$calibration[which(myEMeval.Ensemble$metric.eval=="TSS")],
				myEMeval.Ensemble$calibration[which(myEMeval.Ensemble$metric.eval=="KAPPA")]))

  rownames(summary) <- c("  - Title 3:",
				#"Species name",
				"Statistical algorithms at global level", 
				"Number of replicates with AUC > 0.8 at global level", 
				"AUC of ensemble model at global level", 
				"TSS of ensemble model at global level",
				"KAPPA of ensemble model at global level")

  summary <- rbind(nshsdm_selvars$Summary, summary) #@@@JMB guarda el summary acumulado
  
  nshsdm_data$Summary <- summary 

  #if(save.output){
  #  write.table(summary, paste0("Results/",SpeciesName,"_summary.csv"), sep=",",  row.names=F, col.names=T)
  #}

  attr(nshsdm_data, "class") <- "nshsdm.predict.g" #@@@JMB cambio atributo como propuesta para uso en covariate() y multiply()

  # Logs success or error messages
  #message("\nNSH.SDM.Global.Model executed successfully!\n")

  if(save.output){
    message("Results saved in the following locations:")
    message(paste(
    " - Current and new projections: /Results/Global/Projections/\n",
    "- Replicates statistics: /Results/Global/Values/\n",
    "- Consensus model statistics: /Results/Global/Values/\n",
    "- Variable importance: /Results/Global/Values/"
    ))
    if(!rm.biomod.folder) { 
      message("- BIOMOD results: /Results/Global/Models/\n")
    }
  }

  return(nshsdm_data)

}



