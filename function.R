library(quantreg)
library(cobs)
library(stats)
library(splines)
library(foreach)
library(doParallel)
library(MASS)
library(aod)
library(openxlsx)
library(lme4)
library(dplyr)
library(caret)
library(PearsonDS)

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

# constructing the score for hypothesis testing and approximated variance-covariance for the score
test_stat_separate = function(Y,X,taus,m,test_num){
  y1 = Y[which(Y>0)]
  y = y1
  #y = quantreg::dither(log(y1))
  x = X[which(Y>0),-c(test_num)]

  n = length(y)
  l = ncol(x)
  u = rep(1/sqrt(l),l)
  beta = matrix(0,nrow = length(taus),ncol = ncol(x))
  hat_Nn = floor(n^{1/(2*m+1)})+1
  
  # estimation for the coefficient beta_tau
  for (i in 1:length(taus))
  {
    profile_score_beta = function(v)
    {
      return(profile_likelihood_beta(v,x,y,taus[i],hat_Nn,m))
    }
    beta1 = optim(u,profile_score_beta)$par
    beta[i, ] = beta1/sqrt(sum(beta1*beta1))
  }
  # estimation for the G_tau function
  theta = lapply(1:length(taus),function(ii){
    Z = x %*% beta[ii, ]
    
    # choosing the knots num Nn using BIC defined above
    h1 = ceiling(hat_Nn/2)
    h2 = ceiling(hat_Nn*3/2)
    h = h1:h2
    Res = lapply(1:length(h), function(j){BIC(Z, y, taus[ii], h[j], m)})
    Res = unlist(Res)
    Nn = h[Res==min(Res)]
    basis = splines::bs(Z,df = Nn+m+1, degree = m, intercept = T)
    # using the Nn with the locally minimum BIC to estimate the G function
    gamma = quantreg::rq(y ~ basis-1, tau = taus[ii])
    u = list(gamma = gamma, Nn = Nn, basis = basis)
    return(u)
  })
  
  result = lapply(1:length(taus),function(ii){
    Z = x %*% beta[ii, ]
    L = theta[[ii]]$gamma
    Nn = as.numeric(theta[[ii]]$Nn)
    
    D = theta[[ii]]$basis
    expect = D %*% solve(t(D) %*% D) %*% t(D) %*% X[which(Y>0),]
    # calculate the estimated value of quantile curve G_tau
    Gtau =  predict(L,data.frame(Z))
    
    # calculate the estimated derivative of quantile curve G_tau
    dev = splines2::dbs(Z, df = Nn+m+1, degree = m, intercept = T) %*% L$coefficients
    dev2 = taus[ii] - (y-Gtau<0)
    s1 = 0
    
    X.star = X[which(Y>0),] - expect
    for (j in 1:n){
      s1 = s1 + (X.star[j,]*dev[j]) %*% t(X.star[j,]*dev[j])
    }
    Omega = s1[test_num,test_num] - s1[test_num,-test_num] %*% MASS::ginv(s1[-test_num,-test_num]) %*% s1[-test_num,test_num]
    dev = matrix(dev,nrow = nrow(X.star),ncol = length(test_num))
    hat_X = dev*as.matrix(X.star[,test_num])
    hat_Omega = hat_X %*% solve(Omega) %*% t(hat_X)
    hat_Omega = (1-taus[ii])^{-1}*taus[ii]^{-1}*hat_Omega
    LLLL = list(score = dev2, W = hat_Omega)
    return(LLLL)
  })
  return(result)
}

# using Pearson Type III distribution to get the p-value for hypothesis testing for single-index quantile regression 
fKMQR.test = function(Y, X, tau, m, test_num, score = NULL, K = NULL) {
  ## fast Kernel machine quantile regression
  ### X nonparamteric var (nxp)
  ## Z parametric var (nxq)
  ## Y response var (nx1)
  ## tau quantile
  
  ## Define some auxiliary functions: tr and KRV function
  KRV=function(K,L){
    n=nrow(K)
    A=scale(K,scale=F) ## that is A=HK
    W=scale(L,scale=F) ## that is W=HL
    Fstar=tr(A%*%W)
    mean.krv=tr(A)*tr(W)/(n-1)	## mean of KRV 
    T=tr(A);T2=tr(A%*%A);S2=sum(diag(A)^2)
    Ts=tr(W);T2s=tr(W%*%W);S2s=sum(diag(W)^2)
    temp1=2*((n-1)*T2-T^2)*((n-1)*T2s-Ts^2)/(n-1)^2/(n+1)/(n-2)
    temp21=n*(n+1)*S2- (n-1)*(T^2+2*T2)
    temp22=n*(n+1)*S2s- (n-1)*(Ts^2+2*T2s)
    temp23=(n+1)*n*(n-1)*(n-2)*(n-3)
    temp2=temp21*temp22/temp23
    variance.krv=temp1+temp2		## variance of KRV
    T3=tr(A%*%A%*%A);S3=sum(diag(A)^3);U=sum(A^3);R=t(diag(A))%*%diag(A%*%A);B=t(diag(A))%*%A%*%diag(A)
    T3s=tr(W%*%W%*%W);S3s=sum(diag(W)^3);Us=sum(W^3);Rs=t(diag(W))%*%diag(W%*%W);Bs=t(diag(W))%*%W%*%diag(W)
    t1=n^2*(n+1)*(n^2+15*n-4)*S3*S3s
    t2=4*(n^4-8*n^3+19*n^2-4*n-16)*U*Us
    t3=24*(n^2-n-4)*(U*Bs+B*Us)
    t4=6*(n^4-8*n^3+21*n^2-6*n-24)*B*Bs
    t5=12*(n^4-n^3-8*n^2+36*n-48)*R*Rs
    t6=12*(n^3-2*n^2+9*n-12)*(T*S2*Rs+R*Ts*S2s)
    t7=3*(n^4-4*n^3-2*n^2+9*n-12)*T*Ts*S2*S2s
    t81=(n^3-3*n^2-2*n+8)*(R*Us+U*Rs);t82=(n^3-2*n^2-3*n+12)*(R*Bs+B*Rs)
    t8=24*(t81+t82)
    t9=12*(n^2-n+4)*(T*S2*Us+U*Ts*S2s)
    t10=6*(2*n^3-7*n^2-3*n+12)*(T*S2*Bs+B*Ts*S2s)
    t11=-2*n*(n-1)*(n^2-n+4)*((2*U+3*B)*S3s+(2*Us+3*Bs)*S3)
    t12=-3*n*(n-1)^2*(n+4)*((T*S2+4*R)*S3s+(Ts*S2s+4*Rs)*S3)
    t13=2*n*(n-1)*(n-2)*((T^3+6*T*T2+8*T3)*S3s+(Ts^3+6*Ts*T2s+8*T3s)*S3)
    t14=T^3*((n^3-9*n^2+23*n-14)*Ts^3+6*(n-4)*Ts*T2s+8*T3s)
    t15=6*T*T2*((n-4)*Ts^3+(n^3-9*n^2+24*n-14)*Ts*T2s+4*(n-3)*T3s)
    t16=8*T3*(Ts^3+3*(n-3)*Ts*T2s+(n^3-9*n^2+26*n-22)*T3s)
    t17=-16*(T^3*Us+U*Ts^3)-6*(T*T2*Us+U*Ts*T2s)*(2*n^2-10*n+16)
    t18=-8*(T3*Us+U*T3s)*(3*n^2-15*n+16)-(T^3*Bs+B*Ts^3)*(6*n^2-30*n+24)
    t19=-6*(T*T2*Bs+B*Ts*T2s)*(4*n^2-20*n+24)-8*(T3*Bs+B*T3s)*(3*n^2-15*n+24)
    t201=24*(T^3*Rs+R*Ts^3)+6*(T*T2*Rs+R*Ts*T2s)*(2*n^2-10*n+24)
    t202=8*(T3*Rs+R*T3s)*(3*n^2-15*n+24)+(3*n^2-15*n+6)*(T^3*Ts*S2s+T*S2*Ts^3)
    t203=6*(T*T2*Ts*S2s+Ts*T2s*T*S2)*(n^2-5*n+6)+48*(T3*Ts*S2s+T3s*T*S2)
    t20=-(n-2)*(t201+t202+t203)
    temp31=t1+t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12+t13+t14+t15+t16+t17+t18+t19+t20
    temp32=n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5)
    mom3=temp31/temp32
    skewness.krv= (mom3-3*mean.krv*variance.krv-mean.krv^3)/variance.krv^1.5 ## skewness of KRV
    m1=mean.krv
    m2=variance.krv
    m3=skewness.krv
    shape=4/m3^2
    scale=sqrt(m2)*m3/2
    location=m1-2*sqrt(m2)/m3
    PIIIpars=list(shape,location,scale)
    pv=1-PearsonDS::ppearsonIII(Fstar, params=PIIIpars) 
    return(pv)
  }
  
  tr=function(x){return(sum(diag(x))) }
  #if(!is.null(K)){
  #  K = IBS.kernel(X) ##kernel matrix
  #}
  RRR = test_stat_separate(Y,X,taus = tau,m = m,test_num = test_num)
  w = RRR[[1]]$score
  K1 = RRR[[1]]$W
  
  Kw=w %*% t(w)
  pv=KRV(Kw,K1)
  return(c(pv))
}

# Combining the p-values for logistic regression and single-index quantile regression
Combination = function(y, X, taus = c(0.1,0.25,0.5,0.75,0.9), m, test_num){
  # getting the p-value for logistic regression
  b = 1*(y>0)
  mod.logistic = glm(b ~ X[,-c(1)], family=binomial(link = 'logit'))
  mod.logistic.null = glm(b ~ X[,-c(1,test_num)], family=binomial(link = 'logit'))
  pvalue.logistic = anova(mod.logistic.null, mod.logistic, test="Rao")$`Pr(>Chi)`[2]
  
  # getting the p-value for single-index quantile regression
  MM = rep(0,length(taus))
  for (j in 1:length(taus)){
    MM[j] = fKMQR.test(y,X,tau = taus[j],m,test_num)
  }
  p_value = MM
  
  # getting the combining proportions for each p-values 
  omega = lapply(taus,function(x){
    a = taus
    s = sum(a*(a-0.5<0) + (1-a)*(1 - (a-0.5<0)))
    (x*(x-0.5<0) + (1-x)*(1 - (x-0.5<0)))/s
  })
  rn = length(y[y==0])/length(y)
  omega = (1-rn)*unlist(omega)
  
  # getting the p-value for hypothesis testing
  Test_stat = rn*tan((0.5-pvalue.logistic)*pi) + sum(omega*tan((0.5-p_value)*pi))
  oo1 = 1-pcauchy(Test_stat)
  return(oo1)
}

# estimating the quantile function given samples
proposed.nonsmooth.spline = function(y, Xl, Xq, xl, xq, taus, delta=0.499, m, u){
  
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
  
  return(Quantiles=do.call(rbind, QF))
  
}

AQE = function(Y,Xl,Xq,Xl_1,Xq_1,indexl, indexq, value1, value2,taus,delta = 0.499,m){
  ZIQSI = function(Y,Xl,Xq,xl,xq,m,taus,delta = 0.499){
    n = length(Y)
    l = ncol(X)
    u = rep(1/sqrt(l),l)
    y = Y[which(Y>0)]
    x = Xq[which(Y>0),]

    
    beta = matrix(0,nrow = 99,ncol = ncol(x))
    n0 = length(y)
    hat_Nn = floor(n0^{1/(2*m+1)})+1
    Tau = seq(0.01,1,0.01)
    
    for (i in 1:99){
      profile_score_beta = function(v)
      {
        return(profile_likelihood_beta(v,x,y,Tau[i],hat_Nn,m))
      }
      beta1 = optim(u,profile_score_beta)$par
      beta[i, ] = beta1/sqrt(sum(beta1*beta1))
    }
    
    Nn = rep(0,length(Tau))
    theta = lapply(1:99,function(ii){
      Z = x %*% beta[ii, ] 
      tau = ii*0.01
      # choosing the knots num Nn using BIC defined above
      h1 = ceiling(hat_Nn/2)
      h2 = ceiling(hat_Nn*3/2)
      h = h1:h2
      Res = lapply(1:length(h), function(j){BIC(Z, y, tau, h[j], m)})
      Res = unlist(Res)
      Nn = h[Res==min(Res)]
      
      # using the Nn with the locally minimum BIC to estimate the G function
      gamma = quantreg::rq(y ~ splines::bs(Z,df = Nn+m, degree = m), tau = Tau[ii])$coefficients
      u = data.frame(gamma = gamma, Nn = Nn)
      return(u)
    })
    
    b = 1*(Y > 0)
    logistic = glm(b~Xl,family=binomial)
    gamma = logistic$coef
    
    q = ncol(x)
    
    # estimate the probability of observing a positive y for each new case in xl
    xl.temp = cbind(1, xl)
    
    p_hat = exp(xl.temp %*% gamma)/( 1+exp(xl.temp %*% gamma) )
    QF = lapply(1:length(p_hat),function(ii){
      p_h = p_hat[ii]
      part1 = length( which(taus < (1 - as.vector(p_h))) ) # length of An
      part2 = ( taus[ which(taus >= (1 - as.vector(p_h)) & taus <= (1 - as.vector(p_h) + n^(-delta))) ] - 1 + as.vector(p_h) ) * n^delta # transform in Bn
      part3 = ( taus[ which(taus > (1 - as.vector(p_h) + n^(-delta))) ] - 1 + as.vector(p_h) ) / as.vector(p_h) # transform in Cn
      
      
      taus.s = c(n^(-delta) / p_h, part3) # estimated nominal tau's of y|y>0
      
      if (any(taus.s > 1)){ # problematic case, the probability of observing a positive y is too small 
        
        quant = c( rep(0, part1) ) # simply 0 for all quantiles
        
      } else {
        tauss = round(taus.s,2)
        
        fitted = lapply(1:length(tauss),function(jj){
          jjj = tauss[jj]*100
          beta1 = beta[jjj, ]
          Nn = theta[[jjj]]$Nn[1]
          theta1 = theta[[jjj]]$gamma
          Z = xq[ii,] %*% beta1
          
          Gtau = cbind(1,splines::bs(Z, df = Nn+m, degree = m)) %*% theta1
          return(Gtau)
        })
        fittedy = as.vector(unlist(fitted))
        quant = c(rep(0, part1), fittedy[1]*part2, fittedy[-1] )
        
      }
      
    })
    return(Quantiles=do.call(rbind, QF))
  }
  
  
  n = length(Y)
  
  Xl.temp1 = Xl.temp2 = Xl_1
  Xl.temp1[, indexl] = value1 
  Xl.temp2[, indexl] = value2
  xl = rbind(Xl.temp1, Xl.temp2)
  
  Xq.temp1 = Xq.temp2 = Xq_1
  Xq.temp1[, indexq] = value1 
  Xq.temp2[, indexq] = value2
  xq = rbind(Xq.temp1, Xq.temp2)
  
  quant = ZIQSI(Y,Xl,Xq,xl,xq,m = 3, tau = seq(0,0.99,0.01))
  
  return(quant)
}


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

