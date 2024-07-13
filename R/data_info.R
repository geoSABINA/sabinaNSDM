#' @export
summary.nsdm.input <- function(nsdm_input){
    summary_df <- as.data.frame(nsdm_input$Summary)
   
    return(summary_df)

}

#' @export
summary.nsdm.predict <- function(nsdm_predict){
    summary_df <- as.data.frame(nsdm_predict$Summary)
    
    return(summary_df)

}
