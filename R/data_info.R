#' @export
summary.nsdm.input <- function(nsdm_input){
    summary_df <- as.data.frame(nsdm_input$Summary)
   
    return(summary_df)

}

#' @export
summary.nsdm.finput <- function(nsdm_finput){
    summary_df <- as.data.frame(nsdm_finput$Summary)
   
    return(summary_df)

}

#' @export
summary.nsdm.vinput <- function(nsdm_vinput){
    summary_df <- as.data.frame(nsdm_vinput$Summary)
   
    return(summary_df)

}

#' @export
summary.nsdm.predict <- function(nsdm_predict){
    summary_df <- as.data.frame(nsdm_predict$Summary)
    
    return(summary_df)

}
