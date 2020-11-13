#### main.R ####


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
  names(out) <- c("tau2", "sig2","beta_grp", "beta_crv")
  return(out)
}

#' Extract areas of significance
#'
#' Given a 3-column matrix corresponding to the lower, average, and upper curves,
#' extract the regions of significance, as well as the locations and value of
#' the largest extrema.
#'
#' @usage extract_sig_area(trip)
#'
#' @param trip Matrix with three columns corresponding to the lower bound (
#' first column), average curve (second column), and upper bound (third column).
#'
#' @return 4xk Matrix where each row is a significance region and the columns
#' correspond to the beginning index, ending index, maximum location, and max value.
#'
#' @export

extract_sig_area <- function(trip){
  regions <- as.numeric(trip[,1]*trip[,3] > 0)
  change_points <- diff(regions)
  starts <- which(change_points == 1) + 1
  ends <- which(change_points == -1)
  if(head(regions,1) != 0){ starts <- c(1, starts) }
  if(tail(regions,1) != 0){ ends <- c(ends, nrow(trip)) }

  sig_mat <- cbind(starts, ends)
  if(length(sig_mat) == 0){
    area_mat <- matrix(c(NA,NA,NA,NA), nrow=1)
  } else {
    maxes <- apply(sig_mat, 1, function(r){
      max_loc <- which.max(abs(trip[r[1]:r[2],2])) + r[1] - 1
      max_val <- trip[,2][max_loc]
      return(c(max_loc, max_val))
    })
    area_mat <- cbind(sig_mat, t(maxes))
  }
  colnames(area_mat) <- c("start","end","max_location","max_value")
  return(area_mat)
}

#' Compute Cohen's d
#'
#' Calculates Cohen's d, a measure of effect size as well as confidence interval
#' bounds.
#'
#' @usage cohend(diff, ns, sds, conf = 0.95)
#'
#' @param diff Numeric difference value
#' @param ns Numeric vector containing the sample sizes of the two groups
#' @param sds Numeris vector containing the sample standard deviations of the
#' two groups (at the point of interest if taken on functional data).
#' @param conf Numeric value giving the desired confidence level.
#'
#' @details Arguments are given in vector form so as to be more amenable to the
#' structure of the package output. Confidence interval calculation comes from
#' Hedge and Olkin (2014).
#'
#' @return A vector of length 3 containing the lower bound, the effect size,
#' and the upper bound in that order.
#'
#' @export

cohend <- function(diff, ns, sds, conf = 0.95){
  poolsd <- sqrt(sum(sds^2*(ns-1))/sum(ns-1))
  effect_size <- diff/poolsd
  effect_se <- sqrt(sum(1/ns) + 0.5*effect_size^2/sum(ns))
  es_vec <- effect_size + c(-1,0,1)*qnorm(conf)*effect_se
  names(es_vec) <- c("es_lower","es","es_upper")
  return(es_vec)
}
