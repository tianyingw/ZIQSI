#' proposed.nonsmooth.spline
#'
#' This function estimates the quantiles for given sample(s) on a grid of quantile levels.
#'
#' @param y n*1 vector, the observed outcome
#' @param Xl n*p matrix, the observed covariates in logistic regression
#' @param Xq n*q matrix, the observed covariates in quantile regression
#' @param xl m*p matrix, the new covariates in logistic regression, with which conditional quantile function are estimated
#' @param xq m*q matrix, the new covariates in quantile regression, with which conditional quantile function are estimated
#' @param taus k*1 vector, the grid of target tau's of y
#' @param delta constant, better to keep the default 0.499
#' @param m numeric variable, the order of B-spline function
#' @param u q*1 vector, usually using a vector where all values are equal and the norm is 1, the initial parameter for minimizing the profile likelihood function for quantile regression.
#' @return quantiles of a m*k matrix, each row is the estimated quantiles for each new case
#' @examples
#' # Set the random seed for reproducibility
#' set.seed(10001)
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
#' bet = function(x,u){x^4*u*10^(-5)/6+x^2*u*0.2/3}
#' ## G_tau function
#' func = function(x, tau) {
#'   return(bet(x %*% rbind(bet1(tau), bet2(tau), bet3(tau), bet4(tau), bet5(tau)) + bet0(tau), tau))
#' }
#'
#' # Given covariates for quantile estimation
#' X1 = c(0, 1)  # sex
#' X2 = qnorm(0.5, 28, 2)  # bmi
#' X3 = qnorm(0.5, 92.5, 13)  # waist
#' X4 = qnorm(0.5, 80, 12)  # diastolic_bp
#' X5 = qnorm(0.5, 124, 18.5)  # systolic_bp
#' X0 = cbind(c(rep(X1[1], 1), rep(X1[2], 1)), rep(X2, 2), rep(X3, 2), rep(X4, 2), rep(X5, 2))
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
#' # Estimated quantile curves for the given covariate x
#' A = proposed.nonsmooth.spline(y = y, Xl = cbind(x1, x2, x3, x4, x5),
#'                                  Xq = cbind(x0, x1, x2, x3, x4, x5),
#'                                  xl = X0,
#'                                  xq = cbind(1, X0),
#'                                  taus = seq(0, 0.99, by = 0.01),
#'                                  delta = 0.499,
#'                                  m = 4,
#'                                  u = rep(1/sqrt(6), 6))
#' # plot the estimated quantile curve for the first sample
#' # plot(seq(0,0.99,0.01),A[1,])
#' # plot the estimated quantile curve for the second sample
#' # plot(seq(0,0.99,0.01),A[2,])
#' @export
proposed.nonsmooth.spline = function(y, Xl, Xq, xl, xq, taus = seq(0,0.99,0.01), delta=0.499, m = 3, u){

  # score function
  rho = function(x,tau){tau*x - x*(x<0)}

  profile_likelihood_beta = function(x,X1,Y1,tau1,Nn,m)
    # to get the score function regarding to single-index quantile regression given beta(x)
  {
    Z1 = vector(length = length(Y1));
    y = vector(length = ncol(X1));
    # normalization of beta
    y = x/rep(sqrt(crossprod(x)),length(x))
    # the single index given normalized beta
    Z1 = X1 %*% y ;

    # the structure of single-index B-spline regression
    S = quantreg::rq(Y1~splines::bs(Z1,df = Nn+m+1, degree = m, intercept = T)-1,tau = tau1);
    # the difference between real value and predicted value
    Res = Y1 - predict(S);
    # calculate the score
    Resq = rho(Y1-predict(S),tau1)

    sumit = sum(Resq) # sum the score
    return(sumit)
  }

  # to calculate BIC to get the coefficient theta
  BIC = function(Z1,Y1,tau1,Nn,m){
    n = length(Y1)
    S = quantreg::rq(Y1~splines::bs(Z1,df = Nn+m+1, degree = m, intercept = T)-1,tau = tau1);
    Resq = (Y1-predict(S))*tau1-(Y1-predict(S)<0)*(Y1-predict(S))
    sumit = log(sum(Resq)/n) +  log(n)/2/n*(Nn+m)# sum the score
    return(sumit)
  }
  # logistic regression
  b = 1*(y > 0)
  logistic = glm(b~Xl,family=binomial)
  gamma = logistic$coef
  q = ncol(Xq)

  # estimate the probability of observing a positive y for each new case in xl
  xl.temp = cbind(1, xl)

  p_hat = exp(xl.temp %*% gamma)/( 1+exp(xl.temp %*% gamma) )

  # size of data, affect the data-driven interpolation window
  n = length(y)

  # estimate conditional quantiles for each new case specifically
  QF = lapply(1:length(p_hat), function(ii){

    # one-to-one mapping of the target tau's of y to the nominal tau's of y|y>0
    p_h = p_hat[ii]

    part1 = length( which(taus < (1 - as.vector(p_h))) ) # length of An
    part2 = ( taus[ which(taus >= (1 - as.vector(p_h)) & taus <= (1 - as.vector(p_h) + n^(-delta))) ] - 1 + as.vector(p_h) ) * n^delta # transform in Bn
    part3 = ( taus[ which(taus > (1 - as.vector(p_h) + n^(-delta))) ] - 1 + as.vector(p_h) ) / as.vector(p_h) # transform in Cn


    taus.s = c(n^(-delta) / p_h, part3) # estimated nominal tau's of y|y>0

    if (any(taus.s > 1)){ # problematic case, the probability of observing a positive y is too small

      quant = c( rep(0, part1) ) # simply 0 for all quantiles

    } else { # ordinary case

      # quantile regression on the estimated nominal tau's of y|y>0 with positive y only
      Y = y[which(y>0)]; X = Xq[which(y>0), ]; hat_Nn = floor(length(Y)^{1/(2*m+1)})+1

      # estimate the coefficient beta for the nominal tau's of y|y>0
      beta_taus = lapply(1:length(taus.s), function(jj){
        tau1 = taus.s[jj]

        profile_score_beta = function(x)
        {
          return(profile_likelihood_beta(x,X,Y,tau1,hat_Nn,m))
        }

        beta_tau1 = optim(u,profile_score_beta)$par
        norm_beta = sqrt(beta_tau1 %*% beta_tau1)
        beta_tau1 = beta_tau1/rep(norm_beta,length(beta_tau1))

        return(beta_tau1)

      }
      )
      # estimate the conditional quantiles on Cn
      fitted = lapply(1:length(taus.s),function(kk){

        tau1 = taus.s[kk]
        beta_tau1 = beta_taus[[kk]][1:q]
        Z = X %*% beta_tau1
        h1 = ceiling(hat_Nn/2)
        h2 = ceiling(hat_Nn*3/2)
        h = h1:h2
        Res = lapply(1:length(h), function(j){BIC(Z, Y, tau1, h[j], m)})
        Res = unlist(Res)
        Nn = h[Res==min(Res)]

        S_tau = quantreg::rq(Y ~ splines::bs(Z ,df = Nn+m, degree = m), tau = tau1)
        fit_tau = predict(S_tau, data.frame(Z = xq[ii,  ]%*% beta_tau1))
        return(fit_tau)
      }
      )

      fittedy = as.vector(unlist(fitted))

      # estimate the conditional quantiles on An, Bn and Cn
      quant = c(rep(0, part1), fittedy[1]*part2, fittedy[-1] )

    }

  })

  Quantiles=do.call(rbind, QF);
  colnames(Quantiles) = taus;
  ll = nrow(Quantiles)
  sample_vector = paste("sample", seq(1, ll))
  rownames(Quantiles) = sample_vector
  return(Quantiles)

}

