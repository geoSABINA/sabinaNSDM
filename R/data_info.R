#'@export
summary.hsdm.input <- function(hsdm_input){
    summary_df <- as.data.frame(hsdm_input$Summary)
   
    return(summary_df)

}

#'@export
summary.hsdm.finput <- function(hsdm_finput){
    summary_df <- as.data.frame(hsdm_finput$Summary)
   
    return(summary_df)

}

#'@export
summary.hsdm.predict.g <- function(hsdm_predict){
    summary_df <- as.data.frame(hsdm_predict$Summary)
    
    return(summary_df)

}

#'@export
summary.hsdm.predict.r <- function(hsdm_predict){
    summary_df <- as.data.frame(hsdm_predict$Summary)
    
    return(summary_df)

}

#'@export
summary.hsdm.predict <- function(hsdm_predict){
    summary_df <- as.data.frame(hsdm_predict$Summary)
    
    return(summary_df)

}
