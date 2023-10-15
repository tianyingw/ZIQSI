#### demo
# sample size 
set.seed(10001)
n = 500 

# the probability of Y>0 given covariate x
p = function(x1,x2,x3,x4,x5,gam0=-0.4,gam1=-0.480,gam2=-0.022,gam3=0.021,gam4 = 0.015,gam5 = -0.009){
  lc = gam0 + gam1*x1 + gam2*x2 + gam3*x3 + gam4*x4 + gam5*x5
  exp(lc)/(1+exp(lc))
}

# beta_tau
bet1 = function(x){(0.3*sqrt(x)-x)*2}
bet2 = function(x){x^2*2.2}
bet3 = function(x){
  (x^2-0.5*x+0.6)*2/3
}
bet5 = function(x){-(0.3*x^2-x)*2}
bet4 = function(x){-sin(x*2*pi)*0.1}
bet0 = function(x){-147.7*x-50*x^2-20}
bet = function(x,u){x^4*u*10^(-5)/6+x^2*u*0.2/3}

# G_tau function
func <- function(x, tau)
{
  return(bet(x %*% rbind(bet1(tau),bet2(tau),bet3(tau),bet4(tau),bet5(tau))+bet0(tau),tau))
}

# the given covariate x
X1 = c(0,1) # sex
X2 = qnorm(0.5,28,2) # bmi
X3 = qnorm(0.5,92.5,13) # waist
X4 = qnorm(0.5,80,12) # diastolic_bp
X5 = qnorm(0.5,124,18.5) # systolic_bp
X0 =cbind(c(rep(X1[1], 1), rep(X1[2], 1)),  rep(X2, 2), rep(X3, 2), rep(X4, 2), rep(X5,2))

# given samples
x1 = rbinom(n,1,0.5) # Medicament use
x2 = rnorm(n,28,2) # bmi
x3 = rnorm(n,92.5,13) # waist
x4 = rnorm(n,80,12) # diastolic_bp
x5 = rnorm(n,124,18.5) # systolic_bp
x0 = rep(1,n)
X = cbind(x0,x1,x2,x3,x4,x5)
u = runif(n)
b = rbinom(n,1,p(x1,x2,x3,x4,x5))
w = bet(bet1(u)*x1+bet2(u)*x2+bet3(u)*x3+bet4(u)*x4+bet5(u)*x5+bet0(u),u)
y = b*w
# estimated quantile curve for the given covariate x
A = proposed.nonsmooth.spline(y=y, Xl=cbind(x1,x2,x3,x4,x5), Xq=cbind(x0,x1,x2,x3,x4,x5), xl=X0, xq=cbind(1,X0), taus=seq(0, 0.99, by=0.01), delta=0.499, m = 4, u = rep(1/sqrt(6),6))

# hypothesis testing for the covariate waist given the G_tau function, coeffecients beta_tau and gamma
Combination(y,X,m = 4,test_num = 4)

