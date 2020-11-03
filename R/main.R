#### data_handling.R ####


#-------------------------------------------------------------------------------
# Published Functions with Documentation
#-------------------------------------------------------------------------------

#' Test run
#'
#' Test run description
#'
#' @usage example test
#'
#' @param test parameter
#'
#' @details Test details
#'
#' @return test output
#'
#' @export

test <- function(data_list, p = 20, niter = 10000, nburn = 1000){
  lasts <- cumsum(sapply(data_list, ncol))
  time <- seq(0,1,len = nrow(data_list[[1]]))
  .g2g_fda(as.matrix(as.data.frame(data_list)), lasts-1, time, p, niter, nburn)
}


