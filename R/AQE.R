#' AQE Estimation
#'
#' This function estimates the average quantile difference for two different values
#' of a covariate while controlling for other covariates.
#'
#' @param y A numeric vector of the response variable.
#' @param Xl A matrix of covariates for logistic regression.
#' @param Xq A matrix of covariates for quantile regression.
#' @param indexl An integer index indicating which covariate is being modified for logistic regression.
#' @param indexq An integer index indicating which covariate is being modified for quantile regression.
#' @param value1 The first value of the covariate to be analyzed.
#' @param value2 The second value of the covariate to be analyzed.
#' @param taus A numeric vector of quantiles to estimate.
#' @param delta A numeric value for the delta parameter (default is 0.499).
#' @param m An integer parameter, the order of B-spline in single-index quantile regression.
#' @param u q*1 vector, usually using a vector where all values are equal and the norm is 1, the initial parameter for minimizing the profile likelihood function for quantile regression.
#'
#' @return A numeric vector of quantile estimates.
#' Simulation Example
#'
#' This script demonstrates the simulation of data and the usage of the AQE function.
#'
#' @examples
#' # Set the random seed for reproducibility
#' set.seed(100002)
#'
#' # Sample size
#' n = 500
#'
#' # Generate the true model for simulation
#' ## Probability of Y > 0 given covariate x
#' p = function(x1, x2, x3, x4, x5, gam0 = -0.4, gam1 = -0.480,
#'               gam2 = -0.022, gam3 = 0.021, gam4 = 0.015, gam5 = -0.009) {
#'   lc = gam0 + gam1 * x1 + gam2 * x2 + gam3 * x3 + gam4 * x4 + gam5 * x5
#'   exp(lc) / (1 + exp(lc))
#' }
#'
#' ## Define beta functions
#' bet1 = function(x) {(0.3 * sqrt(x) - x) * 2}
#' bet2 = function(x) {x^2 * 2.2}
#' bet3 = function(x) {(x^2 - 0.5 * x + 0.6) * 2 / 3}
#' bet4 = function(x) {-sin(x * 2 * pi) * 0.1}
#' bet5 = function(x) {-(0.3 * x^2 - x) * 2}
#' bet0 = function(x) {-147.7 * x - 50 * x^2 - 20}
#' bet = function(x,tau){x*tau*20.4}
#' ## G_tau function
#' func = function(x, tau) {
#'   return(bet(x %*% rbind(bet1(tau), bet2(tau), bet3(tau), bet4(tau), bet5(tau)) + bet0(tau), tau))
#' }
#'
#' # Given samples
#' x1 = rbinom(n, 1, 0.5)  # Medicament use
#' x2 = rnorm(n, 28, 2)  # bmi
#' x3 = rnorm(n, 92.5, 13)  # waist
#' x4 = rnorm(n, 80, 12)  # diastolic_bp
#' x5 = rnorm(n, 124, 18.5)  # systolic_bp
#' x0 = rep(1, n)
#'
#' # Covariates of each sample
#' X = cbind(x0, x1, x2, x3, x4, x5)
#' u = runif(n)
#' b = rbinom(n, 1, p(x1, x2, x3, x4, x5))
#' w = bet(bet1(u) * x1 + bet2(u) * x2 + bet3(u) * x3 + bet4(u) * x4 + bet5(u) * x5 + bet0(u), u)
#'
#' # Simulated quantiles of each sample
#' y = b * w
#'
#' P_ZIQSI = AQE(y,X[,-1],X,indexl = 2, indexq = 3, value1 = 23, value2 = 28,taus = 0.5, m = 3, u = rep(sqrt(1/6),6))

#' @export
AQE = function(y, Xl, Xq, indexl, indexq, value1, value2, taus, delta = 0.499, m, u) {

  # Size of data, used in the averaging step
  n = length(y)

  # Create a matrix of new covariates for logistic regression
  Xl.temp1 = Xl.temp2 = Xl
  Xl.temp1[, indexl] = value1
  Xl.temp2[, indexl] = value2
  xl = rbind(Xl.temp1, Xl.temp2)

  # Create a matrix of new covariates for quantile regression
  Xq.temp1 = Xq.temp2 = Xq
  Xq.temp1[, indexq] = value1
  Xq.temp2[, indexq] = value2
  xq = rbind(Xq.temp1, Xq.temp2)

  # Estimate quantiles for the new covariates
  value = proposed.nonsmooth.spline(y = y, Xl = Xl, Xq = Xq, xl = xl, xq = xq, taus = taus, m = m, u = u)

  # AQE estimate = average of Q_y(tau | Xj=value2, X(-j)) - Q_y(tau | Xj=value1, X(-j))
  AQE = mean(value[1:ncol(Xl)]-value[(ncol(Xl)+1):(2*ncol(Xl))])
  return(AQE)
}
