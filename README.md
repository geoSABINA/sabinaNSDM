

<!-- Esto es para comentarios -->

Quitar esto?
```{r, include = FALSE} knitr::opts_chunk$set(collapse = TRUE,comment = "#>",fig.path = "man/figures/README-",out.width = "60%"  message=FALSE, warning=FALSE)
```


# SabinaNSDM: Spaitally-nested Hierarchical Species Distribution Models (NSDM).

##### @Ruben can u upload the logo?
![](directory/logo.png)


## Overview

<strong>SabinaNSDM</strong> is an R package for ....

It provides a set of functions to ...

Include some references of NSDMs, like this: (*Landscape ecology*, <https://doi.org/10.1007/s10980-006-0013-z>; Saura
& Pascual-Hortal, 2007).


### Citing SabinaNSDM package

A research paper detailing the package details and functions is under review, but until it is
published, please reference the packe as following:

<code> <i> Rubén G. Mateo, Jennifer Morales-Barbero, Alejandra Zarzo-Arias,Herlander Lima, Teresa Goicolea 2024. sabinaNSDM: an R package for spatially-nested hierarchical species distribution modeling.
[![DOI](https://zenodo.org/...)
</code> </i>

## Installation

Depends: 
-   R (\> 4.0.0), ...


You can install the released version of SabinaNSDM from
[GitHub](https://github.com) with:

```{r echo=TRUE, eval=FALSE}
library(remotes)
remotes::install_github("geoSABINA/sabinaNSDM
```


## Summary of main *SabinaNSDM* functions

Put here the table with the functions


## Example

This is provides an example on how to use the sabinaNSDM package for conducting spatially-nested hierarchical species distribution modeling. This examples:
The following makes a hiper link to each of the parts:
-   [Load the input data](#load_input_data)
-   [Format data](#format_data)
-   ...

### Load input data {#load_input_data}

....

### Format data {#format_data}

....



#### Así para para poner codigo invisible ITS NOT WORKING
```{r echo=FALSE, eval=FALSE}
ecoregions <- read_sf("D:/Paper_Makurhini/Ejemplos/Ecoregiones_Colombia_amazonas.shp")
Protected_areas <- read_sf("D:/Paper_Makurhini/Ejemplos/PAs_Colombia_amazonas.shp")
```

#### Así para poner código visible:
```{r eval = FALSE}
test_protconn <- MK_ProtConnMult(nodes = Protected_areas, 
                                 region = ecoregions,
                                 area_unit = "ha",
                                 distance = list(type= "centroid"),
                                 distance_thresholds = 10000,
                                 probability = 0.5, 
                                 transboundary = 50000,
                                 plot = TRUE, 
                                 CI = NULL, 
                                 parallel = 4, 
                                 intern = FALSE)
test_protconn[[1]][[1]]
```



#### Así para cargar figuras:

![](man/figures/table_protconn.png)

#### Table:

| Edge depth distance (m) | Core Area (%) |
|-------------------------|:-------------:|
| 100                     |     83.5%     |
| 500                     |    34.14%     |
| 1000                    |     9.78%     |

or with R:

```{r echo=FALSE}
library(formattable)
functions_MK <- data.frame(Function = c("MK_Fragmentation", "distancefile", "MK_RMCentrality", "MK_BCentrality",  "MK_dPCIIC", "MK_dECA", "MK_ProtConn", "MK_ProtConnMult", "MK_ProtConn_raster", "MK_Connect_grid", "test_metric_distance"), Purpose = c("Calculate patch and landscape statistics (e.g., mean size patches, edge density, core area percent, shape index, fractal dimension index, effective mesh size).", "Get a table or matrix with the distances between pairs of nodes. Two Euclidean distances ('centroid' and 'edge') and two cost distances that consider the landscape heterogeneity ('least-cost' and 'commute-time, this last is analogous to the resistance distance of circuitscape, see ’gdistance’ package).", "Estimate centrality measures under one or several dispersal distances (e.g., betweenness centrality, node memberships, modularity). It uses the 'distancefile ()' to calculate the distances of the nodes so they can be calculated using Euclidean or cost distances that consider the landscape heterogeneity.", "Calculate the BC, BCIIC and BCPC indexes under one or several distance thresholds using the command line of CONEFOR. It uses the 'distancefile ()' to calculate the distances of the nodes so they can be calculated using Euclidean or cost distances that consider the landscape heterogeneity", "Calculate the integral index of connectivity (IIC) and probability of connectivity (PC) indices under one or several dispersal distances. It computes overall and index fractions (dPC or dIIC, intra, flux and connector) and the effect of restauration in the landscape connectivity when adding new nodes (restoration scenarios). It uses the 'distancefile()'.", "Estimate the Equivalent Connected Area (ECA) and compare the relative change in ECA (dECA) between time periods using one or several dispersal distances. It uses the 'distancefile()'.", "Estimate the Protected Connected (ProtConn) indicator and fractions for one region using one or several dispersal distances and transboundary buffer areas (e.g., ProtConn, ProtUnconn, RelConn, ProtConn[design], ProtConn[bound], ProtConn[Prot], ProtConn[Within], ProtConn[Contig], ProtConn[Trans], ProtConn[Unprot]). It uses the 'distancefile(). This function estimates what we call the ProtConn delta (dProtConn) which estimates the contribution of each protected area to connectivity in the region (ProtConn value)", "Estimate the ProtConn indicator and fractions for multiple regions. It uses the 'distancefile()'.", "Estimate Protected Connected (ProtConn) indicator and fractions for one region using raster inputs (nodes and region). It uses the 'distancefile()'.", "Compute the ProtConn indicator and fractions, PC or IIC overall connectivity metrics (ECA) in a regular grid. It uses the 'distancefile()'.", "Compare ECA or ProtConn connectivity metrics using one or up to four types of distances, computed in the 'distancefile()' function, and multiple dispersion distances."))

formattable(functions_MK,  align =c("l","l"), list(`Function` = formatter(
              "span", style = ~ style(font.style = "italic"))))

```
