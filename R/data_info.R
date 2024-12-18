#' @title Summary for nsdm.input Objects
#' @description Provides a summary for objects of class \code{nsdm.input}.
#' @param object An object of class \code{nsdm.input}.
#' @param ... Additional arguments (not used).
#' @return A summary in a data frame format.
#' @seealso \code{\link{NSDM.InputData}} \code{\link{NSDM.FormattingData}} \code{\link{NSDM.SelectCovariates}}
#' @export
#' @method summary nsdm.input
summary.nsdm.input <- function(object, ...){
    summary_df <- as.data.frame(object$Summary)
    return(summary_df)
}

#' @title Summary for nsdm.predict Objects
#' @description Provides a summary for objects of class \code{nsdm.predict}.
#' @param object An object of class \code{nsdm.predict}.
#' @param ... Additional arguments (not used).
#' @return A summary in a data frame format.
#' @seealso \code{\link{NSDM.Global}} \code{\link{NSDM.Regional}} \code{\link{NSDM.Covariate}} \code{\link{NSDM.Multiply}}
#' @export
#' @method summary nsdm.predict
summary.nsdm.predict <- function(object, ...){
    summary_df <- as.data.frame(object$Summary)
    return(summary_df)
}


