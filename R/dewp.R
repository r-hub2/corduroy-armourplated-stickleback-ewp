#likelihood functions
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


#' Probability mass function of the three-parameter EWP
#'
#' @param x vector of (positive integer) quantiles.
#' @param lambda centrality parameter
#' @param beta1 lower-tail dispersion parameter
#' @param beta2 upper tail dispersion parameter
#' @param sum_limit summation limit for the normalizing factor
#'
#' @return a vector of probabilities
#' @export
#'
dewp3 <- function(x, lambda, beta1, beta2, sum_limit=max(x)*3){
  stopifnot(is.wholenumber(x))
  if (sum_limit<15){
    warning("sum_limit < 15 detected. A sum_limit of 3 times the maximum count value is recommended, or of
            at least 15, in the case of a small maximum count.")
  }
  return(dewp3_cpp(x, lambda, beta1, beta2, sum_limit))
}

#' Random samples from the three-parameter EWP
#'
#' @param n number of observations
#' @param lambda centrality parameter
#' @param beta1 lower-tail dispersion parameter
#' @param beta2 upper tail dispersion parameter
#' @param sum_limit summation limit for the normalizing factor
#'
#' @return random deviates from the EWP_3 distribution
#' @importFrom stats runif
#' @export
#'
rewp3 <- function(n, lambda, beta1, beta2, sum_limit=30){
    if(any(lambda >= sum_limit)) stop('sum_limit must be larger than lambda')
    probs <- vapply(1:sum_limit, dewp3_cpp, numeric(1), lambda, beta1, beta2, sum_limit)
    sample(x = 1:sum_limit, n, replace = T, prob = probs)
}
