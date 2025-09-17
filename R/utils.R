is.coord.df <- function(dat){
    coords <- c("x", "y")
    return(is.data.frame(dat) &&
           identical(colnames(dat), c("x", "y")))
}


clean_data <- function(mask, data){
    XY <- terra::extract(mask, data)
    XY <- cbind(XY, data)
    XY <- stats::na.omit(XY)[, -c(1:2)]
    XY <- unique(XY)

    return(XY)
}

km_equivalent_from_mask <- function(Mask, XY, Min.Dist) {
  if (as.numeric(Min.Dist) <= 0) return(0)

  # Representative anchor in Mask CRS
  anchor <- data.frame(
    x = stats::median(XY$x, na.rm = TRUE),
    y = stats::median(XY$y, na.rm = TRUE)
  )
  p0 <- terra::vect(anchor, geom = c("x","y"), crs = terra::crs(Mask))
  p1x <- terra::vect(data.frame(x = anchor$x + as.numeric(Min.Dist), y = anchor$y),
                     geom = c("x","y"), crs = terra::crs(Mask))
  p1y <- terra::vect(data.frame(x = anchor$x, y = anchor$y + as.numeric(Min.Dist)),
                     geom = c("x","y"), crs = terra::crs(Mask))

  # Project to lon/lat for geodesic distances
  p0_ll  <- terra::project(p0,  "EPSG:4326")
  p1x_ll <- terra::project(p1x, "EPSG:4326")
  p1y_ll <- terra::project(p1y, "EPSG:4326")

  # Geodesic distances (meters) -> km
  dx_km <- as.numeric(sf::st_distance(sf::st_as_sf(p0_ll), sf::st_as_sf(p1x_ll)))/1000
  dy_km <- as.numeric(sf::st_distance(sf::st_as_sf(p0_ll), sf::st_as_sf(p1y_ll)))/1000

  # Area-preserving equivalent circle radius
  sqrt(dx_km * dy_km)
}

# Format the response (presence/background) or (presence/true absences) and covariates data for BIOMOD2
biomod_format <- function(sp_data, background_data, trueabsences_data, indvar){

    if(!is.null(background_data)) {
      myResp.xy <- rbind(sp_data, background_data)
      row.names(myResp.xy)<-c(1:nrow(myResp.xy))
      myResp <- data.frame(c(rep(1,nrow(sp_data)),
                             rep(NA,nrow(background_data))))
    } else {
      myResp.xy <- rbind(sp_data, trueabsences_data)
      row.names(myResp.xy)<-c(1:nrow(myResp.xy))
      myResp <- data.frame(c(rep(1,nrow(sp_data)),
                             rep(0,nrow(trueabsences_data))))
    }
    names(myResp)<-"pa"
    row.names(myResp)<-c(1:nrow(myResp.xy))
    myExpl <- terra::extract(indvar, myResp.xy, as.df=TRUE)[, -1]

    return(list(myResp = myResp,
                myResp.xy = myResp.xy,
                myExpl = myExpl))
}


