% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CovariateModel.R
\name{NSDM.Covariate}
\alias{NSDM.Covariate}
\title{Perform spatially-nested hierarchical species distribution modeling (NSDM) analysis with the covariate strategy.}
\usage{
NSDM.Covariate(
  nsdm_global,
  algorithms = c("GLM", "GAM", "RF"),
  rm.corr = TRUE,
  CV.nb.rep = 10,
  CV.perc = 0.8,
  CustomModelOptions = NULL,
  metric.select.thresh = 0.8,
  rm.biomod.folder = TRUE,
  save.output = TRUE
)
}
\arguments{
\item{nsdm_global}{An object of class \code{nsdm.predict.g} containing a global model generated using the \code{\link{NSDM.Global}} function.}

\item{algorithms}{(\emph{optional, default} \code{'c("GLM", "GAM", "RF")'}) \cr
A \code{vector} containing the statistical algorithms to use for modeling. Options are \code{'GLM'}, \code{'GAM'}, \code{'GBM'}, \code{'MAXNET'}, \code{'MARS'}, and/or \code{'RF'}.}

\item{rm.corr}{(\emph{optional, default} \code{TRUE}) \cr
A \code{logical} controlling whether environmental covariates correlated with the global model should be removed. The threshold value used for identifying collinearity is the same used with \code{\link{NSDM.SelectCovariates}} function and stored in \code{nsdm_global}}

\item{CV.nb.rep}{(\emph{optional, default} \code{10}) \cr
An \code{integer} corresponding to the number of cross-validation sets (repetitions).}

\item{CV.perc}{(\emph{optional, default} \code{0.8}) \cr
A \code{numeric} between \code{0} and \code{1} defining the percentage of data that will be kept for calibration in each cross-validation set.}

\item{CustomModelOptions}{(\emph{optional, default} \code{NULL}) \cr
A \code{\link{BIOMOD.models.options}} object returned by the \code{\link{bm_ModelingOptions}} to tune models options. If \code{NULL} (the default), biomod2's default parameters are used.}

\item{metric.select.thresh}{(\emph{optional, default} \code{0.8}) \cr
A \code{numeric} between \code{0} and \code{1} corresponding to the minimum scores of AUC below which single models will be excluded from the ensemble model building.}

\item{rm.biomod.folder}{(\emph{optional, default} \code{TRUE}) \cr
A \code{logical} value indicating whether the intermediate biomod2's folders should be removed after processing.}

\item{save.output}{(\emph{optional, default} \code{TRUE}) \cr
A \code{logical} value defining whether the outputs should be saved at local.}
}
\value{
An object of class \code{nsdm.predict.g} containing model information, predictions and evaluation statistics:
\itemize{
\item \verb{$SpeciesName} Name of the species.
\item \verb{$args} A \code{list} containing the arguments used during modeling, including: \code{algorithms}, \code{CV.nb.rep}, \code{CV.perc} and \code{metric.select.thresh}.
\item \verb{$Selected.Variables.Covariate} A \code{character} vector specifying the names of the selected covariates at the regional scale used for the covariate model.
\item \verb{$nbestreplicates} A \code{data.frame} containing  the number of replicates meeting or exceeding the specified \code{metric.select.thresh} for each algorithm used in the modeling.
\item \verb{$current.projections} A \code{list} containing: \code{Pred}, a \code{\link[terra:rast]{PackedSpatRaster}} representing the current continuous (suitability) projection; \code{Pred.bin.ROC} a \code{\link[terra:rast]{PackedSpatRaster}} representing binary projections generated through the optimization of the AUC statistic as a threshold; and \code{Pred.bin.TSS} a \code{\link[terra:rast]{PackedSpatRaster}} representing binary projections generated through the optimization of the TSS statistic as a threshold.
\item \verb{$myEMeval.replicates} Evaluation statistics for each replicate model according to different evaluation metrics (ROC, TSS, KAPPA, ACCURACY, SR, and BOYCE).
\item \verb{$myEMeval.Ensemble} Evaluation statistics for the ensemble model according to different evaluation metrics (ROC, TSS, KAPPA).
\item \verb{$myModelsVarImport} Covariate importance measures for individual models.
\item \verb{$new.projections} A \code{list} containing: \code{Pred.Scenario}, the continuous (suitability) projections onto new scenarios in a \code{\link[terra:rast]{PackedSpatRaster}} format; \code{Pred.bin.ROC.Scenario} a \code{\link[terra:rast]{PackedSpatRaster}} representing binary projections onto new scenarios generated through the optimization of the AUC statistic as a threshold; and \code{Pred.bin.TSS.Scenario} a \code{\link[terra:rast]{PackedSpatRaster}} representing binary projections onto new scenarios generated through the optimization of the TSS statistic as a threshold.
\item \code{Summary} Summary information about the modeling process.
}
}
\description{
This function calibrates, evaluates, and projects a \bold{NSDM} with the \bold{covariate} strategy. It uses a global scale species distribution model output  as an additional covariate to fit a regional scale species distribution model.
}
\details{
This function generates a \bold{NSDM} with the \bold{covariate} strategy. It uses the (\emph{biomod2} package to generate, evaluate, and project species distribution models at the regional scale incorporating the global model prediction as an additional environmental covariate.
If \code{save.output=TRUE}, modeling results are stored out of R in the \emph{Results/} folder created in the current working directory:
\itemize{
\item the \emph{Results/Covariate/Projections/} folder, containing the continuous (suitability) and binary current and new projections. Current projections are named with the species name followed by \file{.Current.tif}, \file{.bin.ROC.tif} and \file{.bin.TSS.tif}. New projections are named with the species name followed by the scenario name, and \file{.bin.ROC.tif}, \file{.bin.TSS.tif} when binary.
\item the \emph{Results/Covariate/Values/} folder, containing replicates statistics, the consensus model statistics, the covariate importance, and the \code{nbestreplicates}, named with the species name and \file{.__replica.csv}, \file{._ensemble.csv}, \file{._indvar.csv} and \file{._nbestreplicates.csv} respectively.
}
}
\examples{
library(sabinaNSDM)

# Load species occurrences
data(Fagus.sylvatica.xy.global, package = "sabinaNSDM")
data(Fagus.sylvatica.xy.regional, package = "sabinaNSDM")

# Load covariates
data(expl.var.global, package = "sabinaNSDM")
data(expl.var.regional, package = "sabinaNSDM")
expl.var.global<-terra::unwrap(expl.var.global)
expl.var.regional<-terra::unwrap(expl.var.regional)

# Load new scenarios
data(new.env, package = "sabinaNSDM")
new.env<-terra::unwrap(new.env)

# Prepare input data
myInputData <- NSDM.InputData(SpeciesName = "Fagus.sylvatica",
			spp.data.global = Fagus.sylvatica.xy.global,
			spp.data.regional = Fagus.sylvatica.xy.regional,
			expl.var.global = expl.var.global,
			expl.var.regional = expl.var.regional,
			new.env = new.env,
			new.env.names = c("Scenario1"),
			Background.Global = NULL,
			Background.Regional = NULL,
			Absences.Global = NULL,
			Absences.Regional = NULL)

# Format the input data
myFormattedData <- NSDM.FormattingData(myInputData,
                                       nPoints=1000)

# Select covariates
mySelectedCovs <- NSDM.SelectCovariates(myFormattedData)

# Perform global scale SDMs
myGlobalModel <- NSDM.Global(mySelectedCovs)

# Perform NSDM analysis using the covariate approach with default settings.
myCovariateModel <- NSDM.Covariate(myGlobalModel)

summary(myCovariateModel)

## Perform NSDM analysis using the covariate approach with custom settings.
# myCovariateModel <- NSDM.Covariate(
#					# Global model output used as input
#					myGlobalModel,
#					# Do not remove correlated covariates 	
#					rm.corr=FALSE,
#					# Algorithms to use for modeling
#					algorithms = c("GBM", "RF", "GLM"),
#					# Number of cross-validation replicates
#					CV.nb.rep = 10,
#					# Percentage of data used in each cross-validation replicate
#					CV.perc = 0.8,
#					# Use default modeling options
#					CustomModelOptions = NULL,
#					# Threshold for selecting models for ensemble
#					metric.select.thresh = 0.8,
#					# Remove the temporary biomod2 output folder
#					rm.biomod.folder = TRUE,
#					# Save the output externally
#					save.output = TRUE)


}
\seealso{
\code{\link{NSDM.InputData}}, \code{\link{NSDM.FormattingData}}, \code{\link{NSDM.SelectCovariates}}, \code{\link{NSDM.Global}}
}
