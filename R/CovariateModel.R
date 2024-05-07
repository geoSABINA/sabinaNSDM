#' @name NSDM.Covariate
#'
#' @title Perform spatially-nested hierarchical species distribution modeling (NSDM) analysis with the covariate strategy.
#'
#' @description This function calibrates, evaluates, and projects a \bold{NSDM} with the \bold{covariate} strategy. It uses a global scale species distribution model output  as an additional covariate to fit a regional scale species distribution model.
#'
#'
#' @param nsdm_global An object of class \code{nsdm.predict.g} containing a global model generated using the \code{\link{NSDM.Global}} function.
#' @param algorithms (\emph{optional, default} \code{'c("GLM", "GAM", "RF")'}) \cr
#' A \code{vector} containing the statistical algorithms to use for modeling. Options are \code{'GLM'}, \code{'GAM'}, \code{'GBM'}, \code{'MAXNET'}, \code{'MARS'}, and/or \code{'RF'}.
#' @param rm.corr (\emph{optional, default} \code{TRUE}) \cr
#' A \code{logical} controlling whether environmental covariates correlated with the global model should be removed. The threshold value used for identifying collinearity is the same used with \code{\link{NSDM.SelectCovariates}} function and stored in \code{nsdm_global}
#' @param CV.nb.rep (\emph{optional, default} \code{10}) \cr
#' An \code{integer} corresponding to the number of cross-validation sets (repetitions).
#' @param CV.perc (\emph{optional, default} \code{0.8}) \cr
#' A \code{numeric} between \code{0} and \code{1} defining the percentage of data that will be kept for calibration in each cross-validation set.
#' @param CustomModelOptions (\emph{optional, default} \code{NULL}) \cr
#' A \code{\link{BIOMOD.models.options}} object returned by the \code{\link{bm_ModelingOptions}} to tune models options. If \code{NULL} (the default), biomod2's default parameters are used.
#' @param metric.select.thresh (\emph{optional, default} \code{0.8}) \cr
#' A \code{numeric} between \code{0} and \code{1} corresponding to the minimum scores of AUC below which single models will be excluded from the ensemble model building.
#' @param save.output (\emph{optional, default} \code{TRUE}) \cr
#' A \code{logical} value defining whether the outputs should be saved at local.
#' @param rm.biomod.folder (\emph{optional, default} \code{TRUE}) \cr
#' A \code{logical} value indicating whether the intermediate biomod2's folders should be removed after processing.
#'
#'
#' @return An object of class \code{nsdm.predict.g} containing model information, predictions and evaluation statistics:
#' - `$SpeciesName` Name of the species.
#' - `$args` A \code{list} containing the arguments used during modeling, including: `algorithms`, `CV.nb.rep`, `CV.perc` and `metric.select.thresh`.
#' - `$Selected.Variables.Covariate` A \code{character} vector specifying the names of the selected covariates at the regional scale used for the covariate model.
#' - `$nbestreplicates` A \code{data.frame} containing  the number of replicates meeting or exceeding the specified \code{metric.select.thresh} for each algorithm used in the modeling.
#' - `$current.projections` A \code{list} containing: \code{Pred}, a \code{\link[terra:rast]{PackedSpatRaster}} representing the current continuous (suitability) projection; \code{Pred.bin.ROC} a \code{\link[terra:rast]{PackedSpatRaster}} representing binary projections generated through the optimization of the AUC statistic as a threshold; and \code{Pred.bin.TSS} a \code{\link[terra:rast]{PackedSpatRaster}} representing binary projections generated through the optimization of the TSS statistic as a threshold.
#' - `$myEMeval.replicates` Evaluation statistics for each replicate model according to different evaluation metrics (ROC, TSS, KAPPA, ACCURACY, SR, and BOYCE).
#' - `$myEMeval.Ensemble` Evaluation statistics for the ensemble model according to different evaluation metrics (ROC, TSS, KAPPA).
#' - `$myModelsVarImport` Covariate importance measures for individual models.
#' - `$new.projections` A \code{list} containing: \code{Pred.Scenario}, the continuous (suitability) projections onto new scenarios in a \code{\link[terra:rast]{PackedSpatRaster}} format; \code{Pred.bin.ROC.Scenario} a \code{\link[terra:rast]{PackedSpatRaster}} representing binary projections onto new scenarios generated through the optimization of the AUC statistic as a threshold; and \code{Pred.bin.TSS.Scenario} a \code{\link[terra:rast]{PackedSpatRaster}} representing binary projections onto new scenarios generated through the optimization of the TSS statistic as a threshold.
#' - `Summary` Summary information about the modeling process.
#'
#'
#' @details
#' This function generates a \bold{NSDM} with the \bold{covariate} strategy. It uses the (\emph{biomod2} package to generate, evaluate, and project species distribution models at the regional scale incorporating the global model prediction as an additional environmental covariate.
#' If `save.output=TRUE`, modeling results are stored out of R in the \emph{Results/} folder created in the current working directory:
#' - the \emph{Results/Covariate/Projections/} folder, containing the continuous (suitability) and binary current and new projections. Current projections are named with the species name followed by \file{.Current.tif}, \file{.bin.ROC.tif} and \file{.bin.TSS.tif}. New projections are named with the species name followed by the scenario name, and \file{.bin.ROC.tif}, \file{.bin.TSS.tif} when binary.
#' - the \emph{Results/Covariate/Values/} folder, containing replicates statistics, the consensus model statistics, the covariate importance, and the \code{nbestreplicates}, named with the species name and \file{.__replica.csv}, \file{._ensemble.csv}, \file{._indvar.csv} and \file{._nbestreplicates.csv} respectively.
#'
#'
#' @seealso \code{\link{NSDM.InputData}}, \code{\link{NSDM.FormattingData}}, \code{\link{NSDM.SelectCovariates}}, \code{\link{NSDM.Global}}
#'
#'
#' @examples
#' library(sabinaNSDM)
#'
#' # Load species occurrences
#' data(Fagus.sylvatica.xy.global, package = "sabinaNSDM")
#' data(Fagus.sylvatica.xy.regional, package = "sabinaNSDM")
#'
#' # Load covariates
#' data(expl.var.global, package = "sabinaNSDM")
#' data(expl.var.regional, package = "sabinaNSDM")
#' expl.var.global<-terra::unwrap(expl.var.global)
#' expl.var.regional<-terra::unwrap(expl.var.regional)
#'
#' # Load new scenarios
#' data(new.env, package = "sabinaNSDM")
#' new.env<-terra::unwrap(new.env)
#'
#' # Prepare input data
#' myInputData<-NSDM.InputData(SpeciesName = "Fagus.sylvatica",
#'				spp.data.global = Fagus.sylvatica.xy.global,
#'				spp.data.regional = Fagus.sylvatica.xy.regional,
#'				expl.var.global = expl.var.global,
#'				expl.var.regional = expl.var.regional,
#'				new.env = new_env,
#'				new.env.names = c("Scenario1"),
#'				Background.Global = NULL,
#'				Background.Regional = NULL)
#'
#' # Format the input data
#' myFormatedData <- NSDM.FormatingData(myInputData,
#'					nPoints=1000)
#'
#' # Select covariates
#' mySelectedCovs <- NSDM.SelectCovariates(myFormattedData)
#'
#' # Perform global scale SDMs
#' myGlobalModel <- NSDM.Global(mySelectedCovs)
#'
#' # Perform NSDM analysis using the covariate approach with default settings.
#' myCovariateModel <- NSDM.Covariate(nsdm_global)
#' 
#' ## Perform NSDM analysis using the covariate approach with custom settings.
#' # myGlobalModel <- NSDM.Covariate(nsdm_global,  	# Global model output used as input
#' #				   rm.corr=FALSE,  	# Do not remove correlated covariates
#' #				   algorithms = c("GBM", "RF", "GLM"), # Algorithms to use for modeling
#' #				   CV.nb.rep = 10,   	# Number of cross-validation replicates
#' #				   CV.perc = 0.8, 	# Percentage of data used in each cross-validation replicate
#' #				   CustomModelOptions = NULL, # Use default modeling options
#' #				   metric.select.thresh = 0.8, # Threshold for selecting models for ensemble
#' #				   rm.biomod.folder = TRUE, # Remove the temporary biomod2 output folder
#' #				   save.output = TRUE)	# Save the output externally
#'
#'
#' @export
NSDM.Covariate <- function(nsdm_global,
			   algorithms=c("GLM","GAM","RF"),
  			   rm.corr=TRUE,
			   CV.nb.rep=10,
			   CV.perc=0.8,
			   CustomModelOptions=NULL,
			   metric.select.thresh = 0.8,
			   rm.biomod.folder=TRUE,
			   save.output = TRUE) {

  if(!inherits(nsdm_global, "nsdm.predict.g")){
      stop("nsdm_global must be an object of nsdm.predict.g class.")
  }
  models <- toupper(algorithms)
  if(any(!models %in% c( "GAM", "GBM", "GLM", "MARS", "MAXNET", "RF"))) {
    stop("Please select at least one valid algorithm (\"GLM\", \"GAM\",\"MARS\",\"GBM\",\"MAXNET\", or \"RF\").")
  }

  SpeciesName <- nsdm_global$Species.Name

  sabina<-nsdm_global[names(nsdm_global) %in% c("Species.Name")]
  sabina$args <- list()
  sabina$args$rm.corr <- rm.corr
  sabina$args$algorithms <- algorithms
  sabina$args$CV.nb.rep <- CV.nb.rep
  sabina$args$CV.perc <- CV.perc
  sabina$args$metric.select.thresh <- metric.select.thresh

  current.projections <- list()
  new.projections <- list()
  new.projections$Pred.ROC.Scenario <- list()
  new.projections$Pred.bin.TSS.Scenario <- list()
  new.projections$Pred.bin.ROC.Scenario <- list()

  # Covariate model calibrated with all the covariates.
  # Regional model excluding climatic variables

  # Add the global model as an additional variable for the regional model.
  SDM.global <- terra::unwrap(nsdm_global$current.projections$Pred) # Unwrap objects
  names(SDM.global) <- c("SDM.global")
  IndVar.Regional.temp <- terra::unwrap(nsdm_global$IndVar.Regional.Selected)  # Unwrap objects
  IndVar.Regional.Covariate <- c(IndVar.Regional.temp, SDM.global)

  myResp.xy <- rbind(nsdm_global$SpeciesData.XY.Regional, nsdm_global$Background.XY.Regional)
  names(myResp.xy)<-c("x","y")
  row.names(myResp.xy)<-c(1:nrow(myResp.xy))
  myResp <- data.frame(c(rep(1,nrow(nsdm_global$SpeciesData.XY.Regional)),rep(NA,nrow(nsdm_global$Background.XY.Regional))))
  names(myResp)<-"pa"
  row.names(myResp)<-c(1:nrow(myResp.xy))
  myExpl <- terra::extract(IndVar.Regional.Covariate, myResp.xy, as.df=TRUE)[, -1]

  if(rm.corr==TRUE) {
    # Remove correlated covariates with global model
    myResp.covsel <- replace(myResp, is.na(myResp), 0)
    myResp.covsel <- as.vector(myResp.covsel)[[1]]
    myExpl.covsel <- myExpl
    myExpl <- covsel::covsel.filteralgo(covdata=myExpl.covsel, pa=myResp.covsel, force="SDM.global", corcut=nsdm_global$corcut)
    IndVar.Regional.Covariate<-IndVar.Regional.Covariate[[which(names(IndVar.Regional.Covariate) %in% colnames(myExpl))]]
  }

  sabina$Selected.Variables.Covariate <- names(myExpl)

  if(save.output){
    fs::dir_create("Results/Covariate/Values/")
    write.csv(names(myExpl), paste0("Results/Covariate/Values/", SpeciesName, ".variables.csv"))
  }

  # Data required for the biomod2 package.
  myBiomodData <- biomod2::BIOMOD_FormatingData(resp.var = myResp,
					resp.xy = myResp.xy,
					expl.var = myExpl,
					resp.name = SpeciesName,
					PA.nb.rep = 1,
					PA.nb.absences = nrow(nsdm_global$Background.XY.Regional),
					PA.strategy = "random")

  # Calibrate and evaluate individual models with specified statistical algorithms.
  # Model training using BIOMOD_Modeling
  myBiomodModelOut <- biomod2::BIOMOD_Modeling(bm.format = myBiomodData,
	                                      modeling.id = "AllModels",
	                                      models = models,
	                                      OPT.user = CustomModelOptions,
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
    nreplicates<-rbind(nreplicates,c(algorithm.i,sum(df_slot$validation[which(df_slot$algo==algorithm.i)] >= metric.select.thresh)))
  }

  sabina$nbestreplicates <- nreplicates

  # Generate and evaluate a single ensemble (i.e.,consensus) model that averages the individual models weighted
  # by the value of the AUC statistic.
  myBiomodEM.ROC  <- biomod2::BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                                 models.chosen = 'all',
                                                 em.by = 'all',
                                                 em.algo = c("EMmean"),
                                                 metric.select = c('ROC'),
                                                 metric.select.thresh = metric.select.thresh,
                                                 var.import = 0,
                                                 metric.eval = c('ROC', "TSS", "KAPPA"),
                                                 seed.val = 42)

  # Project the individual models to the study area at regional scale under training conditions.
  myBiomodProj <- biomod2::BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                        new.env = IndVar.Regional.Covariate,
                                        proj.name = "Current",
                                        models.chosen = 'all',
                                        build.clamping.mask = FALSE)

  # Project the ensemble model to the study area at regional scale under training conditions.
  myBiomodEMProj<-biomod2::BIOMOD_EnsembleForecasting(bm.em = myBiomodEM.ROC,
					bm.proj = myBiomodProj,
					models.chosen = 'all',
					metric.binary = 'all',
					metric.filter = 'all',
					build.clamping.mask = FALSE)

  # Load the model stored by biomod2 and save it in geotif format
  sp.name<-myBiomodData@sp.name
  Pred <- terra::unwrap(myBiomodEMProj@proj.out@val)
  Pred<-terra::rast(wrap(Pred))

  sabina$current.projections$Pred <- setNames(Pred, paste0(SpeciesName, ".Current"))

  if(save.output){
    fs::dir_create("Results/Covariate/Projections/")
    file_path<-paste0("Results/Covariate/Projections/",SpeciesName,".Current.tif")
    terra::writeRaster(Pred, file_path, overwrite=TRUE)
    fs::file_delete(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble.tif"))
  }

  # Binary models
  Pred.bin.ROC <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_ROCbin.tif"))
  Pred.bin.TSS <- terra::rast(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_TSSbin.tif"))
  Pred.bin.ROC<-terra::rast(wrap(Pred.bin.ROC))
  Pred.bin.TSS<-terra::rast(wrap(Pred.bin.TSS))

  sabina$current.projections$Pred.bin.ROC <- setNames(Pred.bin.ROC, paste0(SpeciesName, ".Current.bin.ROC"))
  sabina$current.projections$Pred.bin.TSS <- setNames(Pred.bin.TSS, paste0(SpeciesName,".Current.bin.TSS"))

  if(save.output){
    file_path<-paste0("Results/Covariate/Projections/",SpeciesName,".Current.bin.ROC.tif")
    terra::writeRaster(Pred.bin.ROC, file_path, overwrite=TRUE)
    fs::file_delete(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_ROCbin.tif"))

    file_path<-paste0("Results/Covariate/Projections/",SpeciesName,".Current.bin.TSS.tif")
    terra::writeRaster(Pred.bin.TSS, file_path, overwrite=TRUE)
    fs::file_delete(paste0(sp.name,"/proj_Current/proj_Current_",sp.name,"_ensemble_TSSbin.tif"))
  }

  # Save some results
  # Values of the statistics for each of the replicates
  myEMeval.replicates <- biomod2::get_evaluations(myBiomodModelOut)

  if(save.output){
    fs::dir_create("Results/Covariate/Values/")
    write.csv(myEMeval.replicates,file=paste0("Results/Covariate/Values/",SpeciesName,"_replica.csv"))
    write.csv(nreplicates,file=paste0("Results/Covariate/Values/",SpeciesName,"_nbestreplicates.csv"))
  }

  sabina$myEMeval.replicates <- myEMeval.replicates

  # Values of the statistics of the consensus model
  myEMeval.Ensemble <- biomod2::get_evaluations(myBiomodEM.ROC)

  if(save.output){
    write.csv(myEMeval.Ensemble,file=paste0("Results/Covariate/Values/",SpeciesName,"_ensemble.csv"))
  }

  sabina$myEMeval.Ensemble <- myEMeval.Ensemble

  # Covariate importance
  myModelsVarImport <- biomod2::get_variables_importance(myBiomodModelOut)
  if(save.output){
    write.csv(myModelsVarImport, file = paste0("Results/Covariate/Values/",SpeciesName,"_indvar.csv"), row.names = T)
  }

  sabina$myModelsVarImport <- myModelsVarImport

  # Model projections for future climate scenarios
  ################################################
  if(!is.null(nsdm_global$Scenarios)) {
    Scenarios <- lapply(nsdm_global$Scenarios, terra::unwrap) # Unwrap objects
  }else {Scenarios <-NULL}

  if(length(Scenarios) == 0) {
    warning("No new scenarios for further projections!\n")
  } else {
    for(i in 1:length(Scenarios)) {
      NewClim.temp <- Scenarios[[i]]
      Scenario.name <- names(Scenarios[i])

      SDM.global.future <- terra::unwrap(nsdm_global$new.projections$Pred.Scenario[[i]]) # Unwrap objects
      names(SDM.global.future) <- c("SDM.global")
      NewClim <- c(NewClim.temp, SDM.global.future)
      NewClim<-NewClim[[which(names(NewClim) %in% colnames(myExpl))]]

      myBiomomodProjScenario <- biomod2::BIOMOD_Projection(bm.mod = myBiomodModelOut,
						new.env = NewClim,
						proj.name = Scenario.name,
						models.chosen = "all")

      myBiomodEMProjScenario<-biomod2::BIOMOD_EnsembleForecasting(bm.em = myBiomodEM.ROC,
					bm.proj = myBiomomodProjScenario,
					models.chosen = "all",
					metric.binary = "all",
					metric.filter = "all")

      Pred.Scenario <- terra::unwrap(myBiomodEMProjScenario@proj.out@val)
      Pred.Scenario<-terra::rast(wrap(Pred.Scenario))
      sabina$new.projections$Pred.Scenario[[i]] <- setNames(Pred.Scenario, paste0(SpeciesName,".",Scenario.name))

      if(save.output){
        fs::dir_create(paste0("Results/Covariate/Projections/"))
        file_path <- paste0("Results/Covariate/Projections/",SpeciesName,".",Scenario.name,".tif")
        terra::writeRaster(Pred.Scenario, file_path, overwrite = TRUE)
        fs::file_delete(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble.tif"))

      }

      # Binarized models
      Pred.bin.ROC.Scenario <- terra::rast(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_ROCbin.tif"))
      Pred.bin.TSS.Scenario <- terra::rast(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_TSSbin.tif"))
      Pred.bin.ROC.Scenario<-terra::rast(wrap(Pred.bin.ROC.Scenario))
      Pred.bin.TSS.Scenario<-terra::rast(wrap(Pred.bin.TSS.Scenario))

      sabina$new.projections$Pred.bin.ROC.Scenario[[i]] <- setNames(Pred.bin.ROC.Scenario, paste0(SpeciesName,".",Scenario.name,".bin.ROC"))
      sabina$new.projections$Pred.bin.TSS.Scenario[[i]] <- setNames(Pred.bin.TSS.Scenario, paste0(SpeciesName,".",Scenario.name,".bin.TSS"))

      if(save.output){
        file_path<-paste0("Results/Covariate/Projections/",SpeciesName,".",Scenario.name,".bin.ROC.tif")
        terra::writeRaster(Pred.bin.ROC.Scenario, file_path, overwrite = TRUE)
        fs::file_delete(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_ROCbin.tif"))

        file_path<-paste0("Results/Covariate/Projections/",SpeciesName,".",Scenario.name,".bin.TSS.tif")
        terra::writeRaster(Pred.bin.TSS.Scenario, file_path, overwrite = TRUE)
        fs::file_delete(paste0(sp.name,"/proj_",Scenario.name,"/proj_",Scenario.name,"_",sp.name,"_ensemble_TSSbin.tif"))
      }
    }
  }

  source_folder <- sp.name

  if(rm.biomod.folder){
    # Remove species folder created by biomod2
    fs::dir_delete(source_folder)
  } else {
    if(save.output){
      destination_folder <- paste0("Results/Covariate/Models/",SpeciesName)
      # Move and remove biomod2 results from /sp.name/ to Results/Global/Models/ folder
      fs::dir_create(destination_folder)
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
			 "Statistical algorithms for covariate hierarchical model",
			 paste0("Number of replicates with AUC > ",metric.select.thresh, " for covariate hierarchical model"),
			 "AUC of hierarchical covariate ensemble model",
			 "TSS of hierarchical covariate ensemble model",
			 "KAPPA of hierarchical covariate ensemble model")

  # Wrap objects
  sabina$current.projections <- rapply(sabina$current.projections, terra::wrap, how = "list")

  if(!is.null(nsdm_global$Scenarios)) {
    sabina$new.projections <- rapply(sabina$new.projections, terra::wrap, how = "list")
  }

  sabina$Summary <- summary

  attr(sabina, "class") <- "nsdm.predict"

  # % best replicates messages
  message(sprintf("\n%.2f%% of %s replicates with AUC values >= %.2f.\n", percentage, SpeciesName, metric.select.thresh))
  # save.out messages
  if(save.output){
    message("Results saved in the following local folder/s:")
    message(paste(
    "- Current and new projections: /Results/Covariate/Projections/\n",
    "- Replicates statistics: /Results/Covariate/Values/\n",
    "- Consensus model statistics: /Results/Covariate/Values/\n",
    "- Covariate importance: /Results/Covariate/Values/"
    ))
    if(!rm.biomod.folder) {
    message("- BIOMOD results: /Results/Covariate/Models/\n")
    }
  }

  return(sabina)

}

