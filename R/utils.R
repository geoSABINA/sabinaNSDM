is.coord.df <- function(dat){
    coords <- c("x", "y")
    return(is.data.frame(dat) &&
           identical(colnames(dat), c("x", "y")))
}


clean_data <- function(mask, data){
    XY <- terra::extract(mask, data)
    XY <- cbind(XY, data)
    XY <- na.omit(XY)[, -c(1:2)]
    XY <- unique(XY)

    return(XY)
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


