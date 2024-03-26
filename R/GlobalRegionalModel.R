#' @export
#'
####################
# 3. GLOBAL MODELS
####################

NSH.SDM.Global.Model <- function(nshsdm_selvars,
				algorithms=c( "GLM", "GAM", "RF"), #@@@#TG i set as defaul the same as in varSelection, but user can change them
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
    stop("Please select at least one valid algorithm (\"GLM\", \"GAM\",\"MARS\",\"GBM\",\"MAXNET\", or \"RF\").")
  }

  SpeciesName <- nshsdm_selvars$Species.Name
  Level="Global"

  nshsdm_data<-nshsdm_selvars
  nshsdm_data$Species.Name<-SpeciesName
  nshsdm_data$args <- list()
  nshsdm_data$args$models <- models
  nshsdm_data$args$CV.nb.rep <- CV.nb.rep
  nshsdm_data$args$CV.perc <- CV.perc

  current.projections <- list()
  new.projections <- list()
  new.projections$Pred.Scenario <- list()
  new.projections$Pred.bin.TSS.Scenario <- list()
  new.projections$Pred.bin.ROC.Scenario <- list()

  #GLOBAL SCALE
  #tryCatch({
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
	if (is.null(CustomModelOptions)) {
	  # Use biomod2 default modeling options
	  myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,  #@@@## (This function saves outputs, check them and specify them in the description). Outputs: /NSHSDM/Larix.decidua/models/AllModels/Larix.decidua_PA1_RUN1_GBM
	                                      modeling.id = "AllModels",
	                                      # bm.options = myThecniquesOptions,
	                                      models = models,
	                                      CV.strategy = "random",
	                                      CV.nb.rep = CV.nb.rep, CV.perc = CV.perc,
	                                      weights = NULL, var.import = 3,
	                                      metric.eval = c("ROC", "TSS", "KAPPA", "ACCURACY", "SR", "BOYCE", "MPA"),
	                                      scale.models = FALSE, do.progress = TRUE,
	                                      prevalence = 0.5, seed.val = 42,
	                                      CV.do.full.models = FALSE) # "CV.do.full.models = FALSE" and "var.import=0" to make it faster
	} else {
	  # Use custom modeling options provided by the user
	  # Train and evaluate individual models using BIOMOD_Modeling
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
	nreplicates<-sum(df_slot$validation >= 0.8)
	percentage <- 100 * nreplicates/nrow(df_slot)
	message(sprintf("\n%.2f%% of replicates with AUC values >= 0.8.\n", percentage),	cat("\033[1;34m"))
	cat("\033[0m")
	nreplicates<-data.frame(Algorithm="All",'Number of replicates'=nreplicates)
	for (algorithm.i in models) {
	  nreplicates<-rbind(nreplicates,c(algorithm.i,sum(df_slot$validation[which(df_slot$algo==algorithm.i)] >= 0.8)))
	}
	nshsdm_data$args$nreplicates <- nreplicates
	# Generate and evaluate a single ensemble (i.e.,consensus) model that averages the individual models
	myBiomodEM.ROC  <- biomod2::BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
					models.chosen = 'all',
					em.by = 'all',
					em.algo = c("EMmean"),#@RGM cambiado, ya no es wmean (weighted)
					metric.select = c('ROC'),
					metric.select.thresh = 0.8,
					var.import = 0, #@RGM esto lo he cambiado, creo que solo es necesario en el paso anterior
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
	  dir_create(paste0("Results/",Level,"/Projections/",recurse=T))
	  file_path <- paste0("Results/",Level,"/Projections/",SpeciesName,".Current.tif")
	  terra::writeRaster(Pred, file_path, overwrite=TRUE)
	  message(paste("Projections at global level under training conditions saved in:",file_path,cat("\033[1;34m")))
	  cat("\033[0m")
	  file.remove(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble.tif"))
	}

	# Binary models
	Pred.bin.ROC <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_ROCbin.tif"))
	Pred.bin.TSS <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_TSSbin.tif"))
	Pred.bin.ROC<-terra::rast(wrap(Pred.bin.ROC))
	Pred.bin.TSS<-terra::rast(wrap(Pred.bin.TSS))
	nshsdm_data$current.projections$Pred.bin.ROC <- setNames(Pred.bin.ROC, paste0(SpeciesName, ".Current.bin.ROC"))
	nshsdm_data$current.projections$Pred.bin.TSS <- setNames(Pred.bin.TSS, paste0(SpeciesName,".Current.bin.TSS"))

	# Save some results
	if(save.output){
	  #suffix <- 0
	  file_path <- paste0("Results/",Level,"/Projections/",SpeciesName,".Current.bin.ROC.tif")
	  #old_path <- paste0("Results/",Level,"/Projections/",SpeciesName,".Current.bin.ROC.tif")
	  # while (file.exists(file_path)) {
	  #   suffix <- suffix + 1
	  #   file_path <- gsub("\\.tif$", paste0("_", suffix, ".tif"), old_path)
	  # }
	terra::writeRaster(Pred.bin.ROC, file_path, overwrite=TRUE)
	message(paste("ROC binary projections at global level under training conditions saved in:",file_path,cat("\033[1;34m")))
	cat("\033[0m")
	file.remove(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_ROCbin.tif"))
	#suffix <- 0
	file_path <- paste0("Results/",Level,"/Projections/",SpeciesName,".Current.bin.TSS.tif")
	#old_path <- paste0("Results/",Level,"/Projections/",SpeciesName,".Current.bin.TSS.tif")
	# while (file.exists(file_path)) {
	#   suffix <- suffix + 1
	#   file_path <- gsub("\\.tif$", paste0("_", suffix, ".tif"), old_path)
	# }
	terra::writeRaster(Pred.bin.TSS, paste0("Results/",Level,"/Projections/",SpeciesName,".Current.bin.TSS.tif"), overwrite=TRUE)
	message(paste("TSS binary projections at global level under training conditions saved in:",file_path,cat("\033[1;34m")))
	file.remove(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_TSSbin.tif"))
	cat("\033[0m")
	}

	# Values of the evaluation statistics for each replica
	myEMeval.replicates <- biomod2::get_evaluations(myBiomodModelOut)
	nshsdm_data$myEMeval.replicates <- myEMeval.replicates

	if(save.output){
	  dir_create(paste0("Results/",Level,"/Values/",recurse=T))
	  #suffix<-0
	  file_path <- paste0("Results/",Level,"/Values/",SpeciesName,"_replica.csv")
	  #old_path <- paste0("Results/",Level,"/Values/",SpeciesName,"_replica.csv")
	  # while (file.exists(file_path)) {
	  #   suffix <- suffix + 1
	  #   file_path <- gsub("\\.csv$", paste0("_", suffix, ".csv"), old_path)
	  # }
	write.csv(myEMeval.replicates,file=file_path)
	write.csv(nreplicates,file=paste0("Results/",Level,"/Values/",SpeciesName,"_nreplicates.csv"))
	}

	# Values of the evaluation statistics of the consensus model
	myEMeval.Ensemble <- biomod2::get_evaluations(myBiomodEM.ROC)
	nshsdm_data$myEMeval.Ensemble <- myEMeval.Ensemble

	if(save.output){

	  #suffix<-0
	  file_path <- paste0("Results/",Level,"/Values/",SpeciesName,"_ensemble.csv")
	  #old_path <- paste0("Results/",Level,"/Values/",SpeciesName,"_ensemble.csv")
	  # while (file.exists(file_path)) {
	  #   suffix <- suffix + 1
	  #   file_path <- gsub("\\.csv$", paste0("_", suffix, ".csv"), old_path)
	  # }
	write.csv(myEMeval.Ensemble,file=file_path)
	}

	# Variable importance
	myModelsVarImport <- biomod2::get_variables_importance(myBiomodModelOut)
	nshsdm_data$myModelsVarImport <- myModelsVarImport

	if(save.output){
	  #suffix<-0
	  file_path <- paste0("Results/",Level,"/Values/",SpeciesName,"_indvar.csv")
	  # old_path <- paste0("Results/",Level,"/Values/",SpeciesName,"_indvar.csv")
	  # while (file.exists(file_path)) {
	  #   suffix <- suffix + 1
	  #   file_path <- gsub("\\.csv$", paste0("_", suffix, ".csv"), old_path)
	  # }
	write.csv(myModelsVarImport, file = file_path, row.names = T)
	}


	# Model projections for future climate scenarios
	################################################
	Scenarios <- nshsdm_selvars$Scenarios

	if(length(Scenarios) == 0) {
	  message("There are no future scenarios",	cat("\033[0m"))
	} else{

	for(i in 1:length(Scenarios)) {
	  projmodel <- Scenarios[i]
	#walk(Scenarios, function(projmodel) {
	  new.env <- terra::rast(projmodel)[[nshsdm_selvars$Selected.Variables.Global]]
	  Scenario.name <- path_file(projmodel) |> path_ext_remove()

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
	  Pred.Scenario<-terra::rast(wrap(Pred.Scenario))
	  nshsdm_data$new.projections$Pred.Scenario[[i]] <- setNames(Pred.Scenario, paste0(SpeciesName,".",Scenario.name))

	  if(save.output){
	    dir_create(paste0("Results/",Level,"/Projections/",recurse=T))
	    #suffix <- 0
	    file_path <- paste0("Results/",Level,"/Projections/",SpeciesName,".",Scenario.name,".tif")
	    # old_path <- paste0("Results/",Level,"/Projections/",SpeciesName,".",Scenario.name,".tif")
	    # while (file.exists(file_path)) {
	    #   suffix <- suffix + 1
	    #   file_path <- gsub("\\.tif$", paste0("_", suffix, ".tif"), old_path)
	    # }
	  terra::writeRaster(Pred.Scenario, file_path, overwrite = TRUE)
	  file.remove(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble.tif"))
	  message(paste("Projections at global level under", Scenario.name,"conditions saved in:",file_path,cat("\033[1;34m")))
	  cat("\033[0m")
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
	    message(paste("ROC binary projections at global level under", Scenario.name,"conditions saved in:",file_path,cat("\033[1;34m")))
	    cat("\033[0m")
	    file.remove(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_ROCbin.tif"))

	    file_path <- paste0("Results/",Level,"/Projections/",SpeciesName,".",Scenario.name,".bin.TSS.tif")
	    terra::writeRaster(Pred.bin.TSS.Scenario, file_path, overwrite = TRUE)
	    message(paste("TSS binary projections at global level under", Scenario.name,"conditions saved in:",file_path,cat("\033[1;34m")))
	    cat("\033[0m")
	    file.remove(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_TSSbin.tif"))

	  }

	#}) #walk
	} #for
	  }

	if(rm.biomod.folder || !save.output){ #rm= T save=F
	# Remove species folder created by biomod2
	unlink(sp.name)
	} else {
	# Move biomod2 results to Results/Global/Models folder
  dir_create(paste0("Results/",Level,"/Models/",sp.name))
	source_folder <- sp.name
	destination_folder <- paste0("Results/",Level,"/Models/",sp.name)
	if (file.exists(destination_folder)) {
	  unlink(source_folder, recursive = TRUE)}
	file.rename(from = source_folder, to = destination_folder)
	nshsdm_data$links$biomod.folder <- destination_folder
	unlink(sp.name)
	}

  	gc()
  	results<-  nshsdm_selvars$Summary
  	results<-rbind(results,
  	               c("Statistical algorithms at global level",paste(toupper(algorithms),collapse = ", ")),
  	               c("Number of replicates with AUC > 0.8 at global level",nrow(nshsdm_data$myEMeval.replicates)),
  	               c("AUC of ensemble modle at global level",myEMeval.Ensemble$calibration[which(myEMeval.Ensemble$metric.eval=="ROC")]),
  	               c("TSS of ensemble modle at global level",myEMeval.Ensemble$calibration[which(myEMeval.Ensemble$metric.eval=="TSS")]),
  	               c("KAPPA of ensemble modle at global level",myEMeval.Ensemble$calibration[which(myEMeval.Ensemble$metric.eval=="KAPPA")]))
  	if(save.output){
  	  write.table(results, paste0("Results/",SpeciesName,"_summary.csv"), sep=",",  row.names=F, col.names=T)
  	}
  	nshsdm_selvars$Summary<-results
  	nshsdm_data$Summary<-results
	  attr(nshsdm_data, "class") <- "nshsdm.predict"

  # Logs success or error messages
  message("\nNSH.SDM.Global.Model executed successfully!\n",cat("\033[32m"))

  if(save.output){
  message("Results saved in the following locations:",cat("\033[1;34m"))
  message(paste(
    " - Current and new projections: /Results/Global/Projections/\n",
    "- Replicates statistics: /Results/Global/Values/\n",
    "- Consensus model statistics: /Results/Global/Values/\n",
    "- Variable importance: /Results/Global/Values/\n"
  ),cat("\033[1;34m"))
  if (!rm.biomod.folder) { message("- BIOMOD results: /Results/Global/Models/\n")}
  }
  cat("\033[0m")

  #}, error = function(err) {
  #  message("Error in NSH.SDM.Global.Model:", conditionMessage(err))
  #return(list(result = NULL, error = err))
  #})

  return(nshsdm_data)

}


####################
# 4. REGIONAL MODELS
####################

NSH.SDM.Regional.Models <- function(nshsdm_selvars,
                                    algorithms=c( "GLM", "GAM", "RF"),
                                    CV.nb.rep=10,  #@@@ he cambiado varias cosas aqui
                                    CV.perc=0.8,
                                    CustomModelOptions=NULL, #@@@ he cambiado varias cosas aqui,
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

  nshsdm_data<-nshsdm_selvars
  nshsdm_data$Species.Name<-SpeciesName
  nshsdm_data$args <- list()
  nshsdm_data$args$models <- models
  nshsdm_data$args$CV.nb.rep <- CV.nb.rep
  nshsdm_data$args$CV.perc <- CV.perc

  #option 1 link folder #@@@#TG que es esto? se usa o se puede quitar?
  # nshsdm_data$link <- list()
  # nshsdm_data$link$current.projections <- "/Results/Regional/Projections/"
  # nshsdm_data$link$statistics.replicates <- "/Results/Regional/Values/"
  # nshsdm_data$link$statistics.consensus.model <- "/Results/Regional/Values/"
  # nshsdm_data$link$variable.importance <- "/Results/Regional/Values/"
  # nshsdm_data$link$new.projections <- "/Results/Regional/Projections/"

  #option 2 save.output
  current.projections <- list()
  new.projections <- list()
  new.projections$Pred.Scenario <- list()
  new.projections$Pred.bin.TSS.Scenario <- list()
  new.projections$Pred.bin.ROC.Scenario <- list()

  #tryCatch({
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
	nreplicates<-sum(df_slot$validation >= 0.8)
	percentage <- 100 * nreplicates/nrow(df_slot)
	message(sprintf("\n%.2f%% of replicates with AUC values >= 0.8.\n", percentage),	cat("\033[1;34m"))
	cat("\033[0m")
	nreplicates<-data.frame(Algorithm="All",'Number of replicates'=nreplicates)
	for (algorithm.i in models) {
	  nreplicates<-rbind(nreplicates,c(algorithm.i,sum(df_slot$validation[which(df_slot$algo==algorithm.i)] >= 0.8)))
	}
	nshsdm_data$args$nreplicates <- nreplicates
	# Generate and evaluate a single ensemble (i.e.,consensus) model that averages the individual models
	myBiomodEM.ROC <- biomod2::BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
						models.chosen = 'all',
						em.by = 'all',
						em.algo = c("EMmean"),#@RGM cambiado
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

		# Save some results
	if(save.output){
	dir_create(paste0("Results/",Level,"/Projections/",recurse=T))
	file_path<-paste0("Results/",Level,"/Projections/",SpeciesName,".Current.tif")
	terra::writeRaster(Pred,file_path , overwrite=TRUE)
	message(paste("Projections at regional level under training conditions saved in:",file_path,cat("\033[1;34m")))
	cat("\033[0m")
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
	message(paste("ROC binary projections at regional level under training conditions saved in:",file_path,cat("\033[1;34m")))
	cat("\033[0m")
	file.remove(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_ROCbin.tif"))
  file_path<-paste0("Results/",Level,"/Projections/",SpeciesName,".Current.bin.TSS.tif")
	terra::writeRaster(Pred.bin.TSS,file_path , overwrite=TRUE)
  message(paste("TSS binary projections at regional level under training conditions saved in:",file_path,cat("\033[1;34m")))
  cat("\033[0m")
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
	#nshsdm_data$myEMeval.Ensemble <- myEMeval.Ensemble

	# Variable importance
	myModelsVarImport <- biomod2::get_variables_importance(myBiomodModelOut)
	if(save.output){
	write.csv(myModelsVarImport, file = paste0("Results/",Level,"/Values/",SpeciesName,"_indvar.csv"), row.names = T)
	}
	#nshsdm_data$myModelsVarImport <- myModelsVarImport

	# Model projections for future climate scenarios
	################################################
	Scenarios <- nshsdm_selvars$Scenarios

  if(length(Scenarios) == 0) {
    message("There are no future scenarios",	cat("\033[0m"))
  } else{

	for(i in 1:length(Scenarios)) {              #@@@# for or walk? decide
	 projmodel <- Scenarios[i]
	#walk(Scenarios, function(projmodel) {
	  new.env <- terra::rast(projmodel)[[nshsdm_selvars$Selected.Variables.Regional]]
	  Scenario.name <- path_file(projmodel) |> path_ext_remove()

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
	  message(paste("Projections at regional level under", Scenario.name,"conditions saved in:",file_path,cat("\033[1;34m")))
	  cat("\033[0m")
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
	    message(paste("ROC binary projections at regiona level under", Scenario.name,"conditions saved in:",file_path,cat("\033[1;34m")))
	    cat("\033[0m")
	    file.remove(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_ROCbin.tif"))
	    file_path<-paste0("Results/",Level,"/Projections/",SpeciesName,".",Scenario.name,".bin.TSS.tif")
	    terra::writeRaster(Pred.bin.TSS.Scenario,file_path , overwrite = TRUE)
	    message(paste("TSS binary projections at regional level under", Scenario.name,"conditions saved in:",file_path,cat("\033[1;34m")))
	    cat("\033[0m")
	    file.remove(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_TSSbin.tif"))
	    }

	#}) #walk
	} #for
}
	if(rm.biomod.folder || !save.output){
	# Remove species folder create by biomod2
	unlink(sp.name) #@@@@ sp.name not SpeciesName, to use the name of biomod that might be different thana the one you used
	} else {
	# Move biomod2 results to Results/Regional/Models folder
  dir_create(paste0("Results/",Level,"/Models/",sp.name))
	source_folder <- sp.name
	destination_folder <- paste0("Results/",Level,"/Models/",sp.name)
	if (file.exists(destination_folder)) {
	  unlink(destination_folder, recursive = TRUE)}
	file.rename(from = source_folder, to = destination_folder)
	nshsdm_data$links$biomod.folder <- destination_folder
	unlink(sp.name)
	}

  	gc()


  	results<-nshsdm_selvars$Summary #@@@# careful! This is not charging the summaries of global!
  	results<-rbind(results,
  	               c("Statistical algorithms at regional level",paste(toupper(algorithms),collapse = ", ")),
  	               c("Number of replicates with AUC > 0.8 at regional level",nrow(nshsdm_data$myEMeval.replicates)),
  	               c("AUC of ensemble modle at regional level",myEMeval.Ensemble$calibration[which(myEMeval.Ensemble$metric.eval=="ROC")]),
  	               c("TSS of ensemble modle at regional level",myEMeval.Ensemble$calibration[which(myEMeval.Ensemble$metric.eval=="TSS")]),
  	               c("KAPPA of ensemble modle at regional level",myEMeval.Ensemble$calibration[which(myEMeval.Ensemble$metric.eval=="KAPPA")]))

  	if(save.output){
  	  write.table(results, paste0("Results/",SpeciesName,"_summary.csv"), sep=",",  row.names=F, col.names=T)
  	}
  	nshsdm_data$Summary<-results
	attr(nshsdm_data, "class") <- "nshsdm.predict"

  # Logs success or error messages
  message("\nNSH.SDM.Regional.Models executed successfully!\n",cat("\033[32m"))

  if(save.output){
  message("Results saved in the following locations:",cat("\033[1;34m"))
  message(paste(
    " - Current and new projections: /Results/Regional/Projections/\n",
    "- ReplicateS statistics: /Results/Regional/Values/\n",
    "- Consensus model statistics: /Results/Regional/Values/\n",
    "- Variable importance: /Results/Regional/Values/\n"
  ),cat("\033[1;34m"))
  if (!rm.biomod.folder) { message("- BIOMOD results: /Results/Regional/Models/\n")}
  }
  cat("\033[0m")

  #}, error = function(err) {
  #  message("Error in NSH.SDM.Regional.Model:", conditionMessage(err))
  #  return(list(result = NULL, error = err))
  #})

  return(nshsdm_data)

}
