#### data_handling.R ####


#-------------------------------------------------------------------------------
# Published Functions with Documentation
#-------------------------------------------------------------------------------

#' Multi-Group Bayesian FDA
#'
#' Fits a multi-group Bayesian FDA using penalized P-splines
#'
#' @usage multi_group_fda(data_list, p = 20, niter = 10000, nburn = 1000)
#'
#' @param data_list list of group matrices where columns correspond to observed curves
#'
#' @details Fits Bayesian P-splines model. The model can be given by:
#' \deqn{y_{ij} \sim N(H\beta_{ij}, \sigma^2)\\
#' \beta_{ij} \sim N(\beta_j, \tau^2P^{-1})\\
#' p(\beta_j) \propto 1; p(\sigma^2) \propto \frac{1}{\sigma^2}; p(\tau^2) \propto \frac{1}{\tau^2}}
#'
#' @return List of MCMC chains. Primary interest may lie in analyzing beta_grp chains.
#'
#' @export

multi_group_fda <- function(data_list, p = 20, niter = 10000, nburn = 1000){
  lasts <- cumsum(sapply(data_list, ncol))
  time <- seq(0,1,len = nrow(data_list[[1]]))
  out <- .g2g_fda(as.matrix(as.data.frame(data_list)), lasts-1, time, p, niter, nburn)
  names(out) <- c("tau2", "sig2","beta_grp", "beta_curve")
}



