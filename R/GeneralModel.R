general_nsdm_model <- function(nsdm.obj,
                               model.type,
                               algorithms=c("GLM", "GAM", "RF"),
                               CV.nb.rep=10,
                               CV.perc=0.8,
                               CustomModelOptions=NULL,
                               metric.select.thresh = 0.8,
                               rm.corr = TRUE,
                               save.output=TRUE,
                               rm.biomod.folder=TRUE){

  if(model.type == "Global"){
      scale <- "Global"
  }else if(model.type == "Regional"){
      scale <- "Regional"
  }else if(model.type != "Covariate"){
      scale <- NULL
    stop("model.type must be either Global or Regional or Covariate")
  }

  models <- toupper(algorithms)
  if(any(!models %in% c( "GAM", "GBM", "GLM", "MARS", "MAXNET", "RF"))) {
    stop("Please select at least one valid algorithm (\"GLM\", \"GAM\", \"MARS\", \"GBM\", \"MAXNET\" or \"RF\").")
  }

  SpeciesName <- nsdm.obj$Species.Name

  if(model.type != "Covariate"){
    sabina<-nsdm.obj[!names(nsdm.obj) %in% c("Summary", "args")]
    sabina$args <- list()
    sabina$corcut <- nsdm.obj$args$corcut
  }else{
    sabina<-nsdm.obj[names(nsdm.obj) %in% c("Species.Name")]
    sabina$args <- list()
    sabina$args$rm.corr <- rm.corr
  }
  sabina$args$algorithms <- algorithms
  sabina$args$CV.nb.rep <- CV.nb.rep
  sabina$args$CV.perc <- CV.perc
  sabina$args$metric.select.thresh <- metric.select.thresh

  current.projections <- list()
  new.projections <- list()
  new.projections$Pred.Scenario <- list()
  new.projections$Pred.bin.TSS.Scenario <- list()
  new.projections$Pred.bin.ROC.Scenario <- list()

  if(!is.null(nsdm.obj$Scenarios)) {
    Scenarios <- lapply(nsdm.obj$Scenarios, terra::unwrap)
  } else {
    Scenarios <- NULL
  }

  if(is.null(Scenarios) || length(Scenarios) == 0) {
    warning("No new scenarios for further projections!\n")
  }

  if(model.type == "Covariate"){
    # Covariate model calibrated with all the covariates.
    # Regional model excluding climatic variables

    # Add the global model as an additional variable for the regional model.
    SDM.global <- terra::unwrap(nsdm.obj$current.projections$Pred)
    names(SDM.global) <- c("SDM.global")
    IndVar.Regional.temp <- terra::unwrap(nsdm.obj$IndVar.Regional.Selected)
    IndVar.Regional.Covariate <- c(IndVar.Regional.temp, SDM.global)

    sp_data <- nsdm.obj$SpeciesData.XY.Regional
    background_data <- nsdm.obj$Background.XY.Regional
    trueabsences_data <- nsdm.obj$Absences.XY.Regional
    biomod_inpt <- biomod_format(sp_data,
                                 background_data,
                                 trueabsences_data,
                                 IndVar.Regional.Covariate)
    myResp <- biomod_inpt$myResp
    myResp.xy <- biomod_inpt$myResp.xy
    myExpl <- biomod_inpt$myExpl

    if(rm.corr==TRUE) {
      # Remove correlated covariates with global model
      myResp.covsel <- replace(myResp,
                               is.na(myResp), 0)
      myResp.covsel <- as.vector(myResp.covsel)[[1]]
      myExpl.covsel <- myExpl
      myExpl <- covsel::covsel.filteralgo(covdata=myExpl.covsel, pa=myResp.covsel,
                                         force="SDM.global", corcut=nsdm.obj$corcut)
      matching<-intersect(names(IndVar.Regional.Covariate),names(myExpl))
      myExpl<-myExpl[,match( matching,names(myExpl))]
      IndVar.Regional.Covariate<-IndVar.Regional.Covariate[[which(names(IndVar.Regional.Covariate) %in% colnames(myExpl))]]

    }

    sabina$Selected.Variables.Covariate <- names(myExpl)

  }else{
    # Unwrap objects
    indvar_name <- paste0("IndVar.", scale, ".Selected")
    IndVar.Selected <- terra::unwrap(nsdm.obj[[indvar_name]])
    if(scale == "Global"){
      IndVar.Selected.reg <- terra::unwrap(nsdm.obj[[paste0(indvar_name, ".reg")]])
    }
    sp_data <- nsdm.obj[[paste0("SpeciesData.XY.", scale)]]
    background_data <- nsdm.obj[[paste0("Background.XY.", scale)]]
    trueabsences_data <- nsdm.obj[[paste0("Absences.XY.", scale)]]
    biomod_inpt <- biomod_format(sp_data,
                                 background_data,
                                 trueabsences_data,
                                 IndVar.Selected)
    myResp <- biomod_inpt$myResp
    myResp.xy <- biomod_inpt$myResp.xy
    myExpl <- biomod_inpt$myExpl

  }

  # Prepare data required to calibrate the model with biomod2 package (Result: A BIOMOD.formated.data object)
  myBiomodData <- biomod2::BIOMOD_FormatingData(resp.var = myResp,
                                                resp.xy = myResp.xy,
                                                expl.var = myExpl,
                                                resp.name = SpeciesName,
                                                PA.nb.rep = ifelse(!is.null(trueabsences_data), 0, 1),
                                                PA.nb.absences = ifelse(!is.null(trueabsences_data), 0, nrow(background_data)),
                                                PA.strategy = if(!is.null(trueabsences_data)) NULL else "random")

  # Calibrate and evaluate individual models with specified statistical algorithms
  # Train and evaluate individual models using BIOMOD_Modeling
  myBiomodModelOut <- biomod2::BIOMOD_Modeling(bm.format = myBiomodData,
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

  # Replicates with ROC > metric.select.thresh
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
    nreplicates<-rbind(nreplicates,
                       c(algorithm.i,
                         sum(df_slot$validation[which(df_slot$algo==algorithm.i)] >= metric.select.thresh)))
  }

  sabina$nbestreplicates <- nreplicates

  # Generate and evaluate a single ensemble (i.e.,consensus) model that averages the individual models
  myBiomodEM.ROC <- biomod2::BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                                     models.chosen = 'all',
                                                     em.by = 'all',
                                                     em.algo = c("EMmean", "EMcv"),
                                                     metric.select = c('ROC'),
                                                     metric.select.thresh = metric.select.thresh,
                                                     var.import = 0,
                                                     metric.eval = c('ROC', "TSS", "KAPPA"),
                                                     seed.val = 42)

  # Project the individual models to the study area at regional scale under training conditions
  if(model.type == "Global"){
    biomod_new.env <- IndVar.Selected.reg
  }else if(model.type == "Regional"){
    biomod_new.env <- IndVar.Selected
  }else if(model.type == "Covariate"){
    biomod_new.env <- IndVar.Regional.Covariate
  }

  myBiomodProj <- biomod2::BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                               new.env = biomod_new.env,
                                               proj.name = "Current",
                                               models.chosen = 'all',
                                               build.clamping.mask = FALSE)

  # Project the ensemble model to the study area at regional scale under training conditions
  myBiomodEMProj <- biomod2::BIOMOD_EnsembleForecasting(bm.em = myBiomodEM.ROC,
                                                        bm.proj = myBiomodProj,
                                                        models.chosen = 'all',
                                                        metric.binary = 'all',
                                                        metric.filter = 'all',
                                                        build.clamping.mask = FALSE)

  # Load the model stored by biomod2 and save it in geotif format
  sp.name<-myBiomodData@sp.name
  Pred <- terra::unwrap(myBiomodEMProj@proj.out@val)
  Pred <- terra::rast(terra::wrap(Pred))

  sabina$current.projections$Pred <- setNames(Pred[[1]], paste0(SpeciesName, ".Current"))

  # Uncertainty (coefficient of variation of the ensemble model projections).
  sabina$current.projections$EMcv <- setNames(Pred[[2]], paste0(SpeciesName, ".EMcv"))


  # Binary models
  proj_curr_prefix <- paste0(sp.name,"/proj_Current/proj_Current_",sp.name)
  Pred.bin.ROC <- terra::rast(paste0(proj_curr_prefix, "_ensemble_ROCbin.tif"))
  Pred.bin.TSS <- terra::rast(paste0(proj_curr_prefix, "_ensemble_TSSbin.tif"))
  Pred.bin.ROC <- terra::rast(terra::wrap(Pred.bin.ROC))
  Pred.bin.TSS <- terra::rast(terra::wrap(Pred.bin.TSS))

  sabina$current.projections$Pred.bin.ROC <- setNames(Pred.bin.ROC,
                                                      paste0(SpeciesName, ".Current.bin.ROC"))
  sabina$current.projections$Pred.bin.TSS <- setNames(Pred.bin.TSS,
                                                      paste0(SpeciesName, ".Current.bin.TSS"))

  # Values of the evaluation statistics for each replica
  myEMeval.replicates <- biomod2::get_evaluations(myBiomodModelOut)

  sabina$myEMeval.replicates <- myEMeval.replicates

  # Values of the evaluation statistics of the consensus model
  myEMeval.Ensemble <- biomod2::get_evaluations(myBiomodEM.ROC)

  sabina$myEMeval.Ensemble <- myEMeval.Ensemble

  # Variable importance
  myModelsVarImport <- biomod2::get_variables_importance(myBiomodModelOut)

  sabina$myModelsVarImport <- myModelsVarImport

  # Model projections for future climate scenarios
  if(!is.null(Scenarios)) {
    for(i in 1:length(Scenarios)) {
      if(model.type == "Covariate"){
        new.env.tmp <- Scenarios[[i]]
        SDM.global.future <- terra::unwrap(nsdm.obj$new.projections$Pred.Scenario[[i]])
        names(SDM.global.future) <- c("SDM.global")
        new.env <- c(new.env.tmp, SDM.global.future)
        new.env<-new.env[[which(names(new.env) %in% colnames(myExpl))]]

      }else{
        Selected.Variables <- nsdm.obj[[paste0("Selected.Variables.", model.type)]]
        new.env <- Scenarios[[i]][[Selected.Variables]]
      }
      Scenario.name <- names(Scenarios[i])

      #Project the individual models
      myBiomomodProjScenario <- biomod2::BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                                           new.env = new.env,
                                                           proj.name = Scenario.name,
                                                           models.chosen = "all")

      #Project the ensemble model
      myBiomodEMProjScenario<-biomod2::BIOMOD_EnsembleForecasting(bm.em = myBiomodEM.ROC,
                                                                  bm.proj = myBiomomodProjScenario,
                                                                  models.chosen = "all",
                                                                  metric.binary = "all",
                                                                  metric.filter = "all")

      Pred.Scenario <- terra::unwrap(myBiomodEMProjScenario@proj.out@val)
      Pred.Scenario <- terra::rast(terra::wrap(Pred.Scenario))

      sabina$new.projections$Pred.Scenario[[i]] <- setNames(Pred.Scenario,
                                                            paste0(SpeciesName,".",Scenario.name))

      # Binarized models
      scenario_prefix <- paste0(sp.name,"/proj_",Scenario.name,"/proj_", Scenario.name,"_",sp.name)
      Pred.bin.ROC.Scenario <- terra::rast(paste0(scenario_prefix,"_ensemble_ROCbin.tif"))
      Pred.bin.TSS.Scenario <- terra::rast(paste0(scenario_prefix,"_ensemble_TSSbin.tif"))
      Pred.bin.ROC.Scenario<-terra::rast(terra::wrap(Pred.bin.ROC.Scenario))
      Pred.bin.TSS.Scenario<-terra::rast(terra::wrap(Pred.bin.TSS.Scenario))

      sabina$new.projections$Pred.bin.ROC.Scenario[[i]] <- setNames(Pred.bin.ROC.Scenario,
                                                                    paste0(SpeciesName,".",Scenario.name,".bin.ROC"))
      sabina$new.projections$Pred.bin.TSS.Scenario[[i]] <- setNames(Pred.bin.TSS.Scenario,
                                                                    paste0(SpeciesName,".",Scenario.name,".bin.TSS"))

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
                         paste0("Statistical algorithms for ", model.type, " model"),
                         paste0("Number of replicates with AUC > ",metric.select.thresh, " for ", model.type, " model"),
                         paste0("AUC of ensemble for ", model.type, " model"),
                         paste0("TSS of ensemble for ", model.type, " model"),
                         paste0("KAPPA of ensemble for ", model.type, " model"))

  sabina$Summary <- summary

  # Wrap objects
  sabina$current.projections <- rapply(sabina$current.projections,
                                       terra::wrap, how = "list")
  if(!is.null(nsdm.obj$Scenarios)) {
    sabina$new.projections <- rapply(sabina$new.projections,
                                     terra::wrap, how = "list")
  }


  # % best replicates messages
  message(sprintf("\n%.2f%% of %s replicates with AUC values >= %.2f.\n", percentage, SpeciesName, metric.select.thresh))

  if(save.output){

    values_path <- paste0("Results/", model.type, "/Values/")
    fs::dir_create(values_path)

    projection_path <- paste0("Results/", model.type,"/Projections/")
    fs::dir_create(projection_path)

    # Variables
    write.csv(names(myExpl),
            paste0(values_path, SpeciesName, ".variables.csv"))

    #Ensemble
    file_path <- paste0(projection_path, SpeciesName, ".Current.tif")
    terra::writeRaster(Pred[[1]], file_path, overwrite=TRUE)
    fs::file_delete(paste0(proj_curr_prefix, "_ensemble.tif"))

    file_path <- paste0(projection_path, SpeciesName, ".EMcv.tif")
    terra::writeRaster(Pred[[2]], file_path, overwrite=TRUE)

    # Pred
    file_path <- paste0(projection_path, SpeciesName, ".Current.bin.ROC.tif")
    terra::writeRaster(Pred.bin.ROC, file_path, overwrite=TRUE)
    fs::file_delete(paste0(proj_curr_prefix, "_ensemble_ROCbin.tif"))

    file_path <- paste0(projection_path, SpeciesName,".Current.bin.TSS.tif")
    terra::writeRaster(Pred.bin.TSS, file_path, overwrite=TRUE)
    fs::file_delete(paste0(proj_curr_prefix, "_ensemble_TSSbin.tif"))

    # Evaluation replicates
    write.csv(myEMeval.replicates,file=paste0(values_path,SpeciesName,"_replica.csv"))
    write.csv(nreplicates,file=paste0(values_path,SpeciesName,"_nbestreplicates.csv"))

    # Evaluation consensus
    write.csv(myEMeval.Ensemble,file=paste0(values_path,SpeciesName,"_ensemble.csv"))

    # Variabale Importance
    file_path <- paste0(values_path,SpeciesName,"_indvar.csv")
    write.csv(myModelsVarImport, file = file_path, row.names = T)

    # New scenarios
    if(!is.null(Scenarios)){

      file_path <- paste0(projection_path, SpeciesName,".",Scenario.name,".tif")
      terra::writeRaster(Pred.Scenario, file_path, overwrite = TRUE)
      fs::file_delete(paste0(scenario_prefix,"_ensemble.tif"))

      file_path <- paste0(projection_path,SpeciesName,".",Scenario.name,".bin.ROC.tif")
      terra::writeRaster(Pred.bin.ROC.Scenario, file_path, overwrite = TRUE)
      fs::file_delete(paste0(scenario_prefix,"_ensemble_ROCbin.tif"))

      file_path <- paste0(projection_path,SpeciesName,".",Scenario.name,".bin.TSS.tif")
      terra::writeRaster(Pred.bin.TSS.Scenario, file_path, overwrite = TRUE)
      fs::file_delete(paste0(scenario_prefix, "_ensemble_TSSbin.tif"))

    }

    message("Results saved in the following local folder/s:")
    message(paste(
    " - Current and new projections: ", projection_path, "\n",
    "- Ensemble model CV: ", projection_path, "\n",
    "- Replicates statistics: ", values_path, "\n",
    "- Consensus model statistics: ", values_path, "\n",
    "- Covariate importance: ", values_path, "\n"
    ))
    if(!rm.biomod.folder) {
      message(" - BIOMOD results: /Results/", model.type, "/Models/\n")
    }
  }
  source_folder <- sp.name
  if(rm.biomod.folder){
    # Remove species folder created by biomod2
    fs::dir_delete(source_folder)
  } else {
    if(save.output){
      destination_folder <- paste0("Results/", model.type,"/Models/",SpeciesName)
      # Move and remove biomod2 results from /sp.name/ to Models/ folder
      fs::dir_create(destination_folder)
      fs::dir_copy(source_folder, destination_folder, overwrite = TRUE)
      fs::dir_delete(source_folder)
    }
  }

  return(sabina)

}
