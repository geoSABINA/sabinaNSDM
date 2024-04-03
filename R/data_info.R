#'@export
summary.nshsdm.input <- function(nshsdm_input){
    summary_df <- as.data.frame(nshsdm_input$Summary)
   
    return(summary_df)

}


#'@export
summary.nshsdm.predict <- function(nshsdm_global){
    summary_df <- as.data.frame(nshsdm_global$Summary)
    
    return(summary_df)

}
