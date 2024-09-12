

<!-- Esto es para comentarios -->


[![DOI]( https://doi.org/10.1111/2041-210X.14417)



<img width="35%" align= "right" alt="logo_s-1" src="https://github.com/geoSABINA/sabinaNSDM/assets/168073517/d29288b9-c1a7-47aa-8753-918c931e4c53"/>




# sabinaNSDM: Spatially-nested Hierarchical Species Distribution Models (NSDM)



 <!-- <img width="252" alt="logo_s-1" src="https://github.com/geoSABINA/sabinaNSDM/assets/168073517/d29288b9-c1a7-47aa-8753-918c931e4c53">-->

 

   


## Overview

The <strong>sabinaNSDM</strong> R package generates <strong>spatially-nested hierarchical species distribution models (NSDMs)</strong> that integrates species distribution models (SDMs) at various spatial scales to address niche truncation and produce more reliable predictions than traditional non-hierarchical SDMs. <strong>sabinaNSDM</strong> combines two SDMs calibrated with species occurrences and environmental covariates at global and regional scales. The global-scale model allows capturing extensive ecological niches, while the regional-scale model features high-resolution drivers of species distributions. This toolkit is designed to facilitate the implementation of NSDMs for ecologists, conservationists, and researchers aiming to produce more reliable species distribution predictions.

<strong>sabinaNSDM</strong> streamlines the data preparation, calibration, integration, and projection of models across two scales. It automates (if necessary) the generation of background points, spatial thinning of species occurrences and absences (if available), covariate selection, single-scale modeling (global and regional), and the generation of NSDMs using two approaches (“covariate” and “multiply”). <strong>sabinaNSDM</strong> models use an ensemble modeling approach that combines multiple statistical techniques with the <strong>biomod2</strong> package (*Ecography*, <https://doi.org/10.1111/j.0906-7590.2004.03673.x>, Thuiller et al. 2009), thinning of species occurrences and absences with the <strong>ecospat</strong> package (*Ecography*, <https://doi.org/10.1111/ecog.02671.x>, Di Cola et al. 2017), and covariate selection of the <strong>covsel</strong>  package (*Ecological Informatics*,<https://doi.org/10.1016/j.ecoinf.2023.102080>, Adde et al. 2023b).


### Citing sabinaNSDM package

A research paper detailing the package details and functions is under review, but until it is published, please reference the package as following:

<code> <i> Mateo, R. G., Morales-Barbero, J., Zarzo-Arias, A., Lima, H., Gómez-Rubio, V., & Goicolea, T. (2024). sabinaNSDM: An R package for spatially nested hierarchical species distribution modelling. Methods in Ecology and Evolution, 00, 1–8
[![DOI]( https://doi.org/10.1111/2041-210X.14417)]( https://doi.org/10.1111/2041-210X.14417)
</code> </i>


## Installation

Depends:  R (\> 4.3.0)


You can install the released version of sabinaNSDM from
[GitHub](https://github.com) with:

```{r echo=TRUE, eval=FALSE}
library(remotes)
remotes::install_github("geoSABINA/sabinaNSDM")
```


## Summary of main sabinaNSDM functions


| Overall Step            | Function          | Objective                                       |
|-----------------|-------------------|:-----------------:|
| Data preparation| NSDM.InputData    | Provides the package with the species occurrences and environmental covariates at both global and regional scales|
|                 | NSDM.FormattingData| Background data generation and species occurrences (and absences if available) thinning|
|                 | NSDM.SelectCovariates| Selects uncorrelated and the most relevant environmental covariates|
|Single scale modeling | NSDM.Global| Calibrates, evaluates, and projects ensemble models at the global scale|
|     | NSDM.Regional    | Calibrates, evaluates, and projects ensemble models at the regional scale    |
|Nested modeling  | NSDM.Covariate    | Generate spatially-nested hierarchical species distribution models with the covariate approach. The covariate approach uses the output of the global model as an additional covariate for the regional scale model   |
|     | NSDM.Multiply   | Generate spatially-nested hierarchical species distribution models with the multiply approach. The multiply approach averages the global and regional models    |

## Example

This is an example on how to use the <strong>sabinaNSDM</strong> package for conducting spatially-nested hierarchical species distribution modeling. 
-   [Data preparation](#data_preparation)
-   [Single scale modeling](#single_scale_modeling)
-   [Nested modeling](#nested_modeling)

### Data preparation <a name="data_preparation">

First, set your working directory and load the required packages.

```{r eval = FALSE}
# setwd("/path/to/your/project")
#
if(!requireNamespace("sabinaNSDM", quietly = TRUE)) {
  install.packages("terra")
  install.packages("biomod2")
  install.packages("ecospat")
  install.packages("fs")
  install.packages("sgsR")
  install.packages("remotes")
  if(!requireNamespace("covsel", quietly = TRUE)) {
    remotes::install_github("N-SDM/covsel")
  }
  
  library(terra)
  library(covsel)
  library(biomod2)
  library(ecospat)
  library(fs)
  library(sgsR)
  remotes::install_github("geoSABINA/sabinaNSDM")
}

# Load the sabinaNSDM package
library(sabinaNSDM)
```
Define the species name
```{r eval = FALSE}
SpeciesName <- "Fagus.sylvativa"
```
Load species occurrence and environmental covariates data. Species occurrence *data.frame* must include only two columns: “x” and ”y” coordinates. No row names. The coordinate projection must match that used for the covariates

```{r eval = FALSE}
 # Species occurrences
 data(Fagus.sylvatica.xy.global, package = "sabinaNSDM")
 spp.data.global <- Fagus.sylvatica.xy.global
 data(Fagus.sylvatica.xy.regional, package = "sabinaNSDM")
 spp.data.regional <- Fagus.sylvatica.xy.regional
```
The covariates for each spatial scale (i.e., global and regional) should be provided as *SpatRaster*, with each band corresponding to a different covariate. The regional-scale *SpatRaster*  must include all the covariates included in the global-scale file, and it can additionally include covariates only available at this level.

 ```{r eval = FALSE}
 data(expl.var.global, package = "sabinaNSDM")
 data(expl.var.regional, package = "sabinaNSDM")
 expl.var.global <- terra::unwrap(expl.var.global)
 expl.var.regional <- terra::unwrap(expl.var.regional)
```
Additionally, regional-scale *SpatRaster*  corresponding to the covariates used to project the models at different scenarios (i.e., new scenarios) can be provided

 ```{r eval = FALSE}
# new scenarios
 data(new.env, package = "sabinaNSDM")
 new.env <- terra::unwrap(new.env)
```
Load the required data for the package with the *NSDM.InputData()* function
 ```{r eval = FALSE}
nsdm_input <- NSDM.InputData(SpeciesName = SpeciesName,
                    spp.data.global = Fagus.sylvatica.xy.global, 
                    spp.data.regional = Fagus.sylvatica.xy.regional, 
                    expl.var.global = expl.var.global, 
                    expl.var.regional = expl.var.regional,
                    new.env = new.env,
                    new.env.names = "scenario1",
                    Background.Global = NULL, 
                    Background.Regional = NULL,
                    Absences.Global = NULL,
                    Absences.Regional = NULL)

```
Format the data with the *NSDM.FormattingData()* function. This function generates random or stratified background points for model calibration when no specific background or true absence data was loaded in the *NSDM.InputData()* function. Additionally, it applies spatial thinning to species occurrences and absences (if available) to remove duplicates and enforce a minimum distance criterion (by default the resolution of the variables). 
 ```{r eval = FALSE}
nsdm_finput <- NSDM.FormattingData(nsdm_input,
                nPoints = 100, # number of background points
                Min.Dist.Global = "resolution",
                Min.Dist.Regional = "resolution",
                Background.method = "random", # method “random" or "stratified” to generate background points 
                save.output = TRUE) #save outputs locally

```
*NSDM.SelectCovariates()* function selects the most relevant and uncorrelated environmental covariates for both global and regional scales.
```{r eval = FALSE}
nsdm_selvars <- NSDM.SelectCovariates(nsdm_finput,
                maxncov.Global = 5,   # Max number of covariates to be selected at the global scale
                maxncov.Regional = 7, # Max number of covariates to be selected at the regional scale
                corcut = 0.7, #  correlation threshold
                algorithms = c("glm","gam","rf"),
                ClimaticVariablesBands = NULL, # covariate bands to be excluded in the covariate selection at the regional scale
                save.output = TRUE)

```


### Single scale modeling <a name="single_scale_modeling">

*NSDM.Global()* function generates the global component of the NSDM.

```{r eval = FALSE}
nsdm_global <- NSDM.Global(nsdm_selvars,
                algorithms = c("GAM","GBM", "RF", "MAXNET","GLM"),# Statistical algorithms used for modeling
                CV.nb.rep = 10, # number of cross-validation repetitions
                CV.perc = 0.8, # percentage of the data will be used for training in each cross-validation fold
                metric.select.thresh = 0.8, #  AUC threshold to include replicates in the final ensemble model
                CustomModelOptions = NULL, # Allows users to apply custom modeling options. 
                save.output = TRUE, 
                rm.biomod.folder = TRUE) # Remove the temporary folders created by `biomod2` 

```

*NSDM.Regional()* function generates the regional component of the NSDM.


```{r eval = FALSE}
nsdm_regional <- NSDM.Regional(nsdm_selvars,
                algorithms = c("GAM","GBM", "RF", "MAXNET","GLM"),
                CV.nb.rep = 10,
                CV.perc = 0.8,
                # metric.select.thresh = 0.8,
                CustomModelOptions = NULL, 
                save.output = TRUE,
                rm.biomod.folder = TRUE)

```



### Nested modeling <a name="nested_modeling">

*NSDM.Covariate()* function generates a NSDM with the covariate strategy. The covariate strategy incorporates the output of global models as an additional covariate in the regional model. 
```{r eval = FALSE}
nsdm_covariate <- NSDM.Covariate(nsdm_global,
                algorithms = c("GAM","GBM", "RF", "MAXNET","GLM"),
                rm.corr=TRUE,
                CV.nb.rep = 10,
                CV.perc = 0.8,
                # metric.select.thresh = 0.8,
                CustomModelOptions = NULL,
                save.output = TRUE,
                rm.biomod.folder = TRUE)
```
*NSDM.Multiply()* function generates a NSDM with the multiply strategy. The covariate averages the output of the global and the regional models. 

```{r eval = FALSE}
nsdm_multiply <- NSDM.Multiply(nsdm_global,
                nsdm_regional,
                method = "Arithmetic", # Method for averate model outputs: "Arithmetic" or "Geometric" mean
                rescale = FALSE,
                save.output=TRUE)
```



<!-- Así para cargar figuras: ![](man/figures/table_protconn.png)--> 



<!-- Multiple studies have shown that NSDMs outperform their non- hierarchical counterparts in various applications (*Journal of Vegetation Science*, <https://doi.org/10.1111/jvs.12726>,Mateo et al. 2019; *Frontiers in Ecology and Evolution*, <https://doi.org/10.3389/fevo.2022.944116>,Chevalier et al. 2022; *Ecography*, <https://doi.org/10.1016/j.ecoinf.2023.102080>; Goicolea et al. in press)--> 
