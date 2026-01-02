###############################################################
#                 Test for Maxwell Distribution (nu = 3)
#                       All GOF Tests Included
#                       n = 25,  a = 0.5,1,1.5
###############################################################

library(stats4)
library(VGAM)
library(stargazer)
library(extraDistr) 
library(EnvStats)     # for levy, invgamma, lomax, invgauss etc
library(statmod)
###############################################################

rmaxwell <- function(n, parm) {
  y <- rchisq(n, df = 3)
  X <- sqrt(y) * parm
  return(X)
} 

pmaxwell <- function(x, parm) {
  pchisq((x / parm)^2, df = 3)
}
# Unified CDF for nu = 3 (Maxwell)
pchi_family <- function(x, parm, nu = 3) pmaxwell(x, parm)

###############################################################
#   Parameter Estimation (MLE)
###############################################################
nu=3
estimation <- function(x) sqrt(mean(x^2) / nu)

###############################################################
#   Tn Test Statistic
###############################################################
Tn <- function(X, parm, nu, a) {
  data = X / parm
  n = length(data)
  mdata = matrix(data, n, n)
  pdata = mdata + t(mdata) + a
  multdata = (nu - 1 - mdata^2) * t(nu - 1 - mdata^2)
  
  SUM = 2 / pdata^3 -
    2 * (nu - 1 - mdata^2) / (mdata * pdata^2) +
    multdata / (mdata * t(mdata) * pdata)
  
  sum(SUM) / n
}
###############################################################
#   Comparison Test  GV Avhad, B Ebner (2025)
###############################################################
AEn<-function(X, parm){ # here X is an array
  X = X / parm  
  FH = pmaxwell(X, parm)
  return(mean((X^2-2)*(1-FH)+((X^2-2)/X)*sqrt(2/pi)*(2-(X^2+2)*exp(-X^2/2))+FH-1))
}
###############################################################
#   Comparison Tests (KS, CVM, AD, MA, Watson, ZA, ZB, ZC)
###############################################################
KSn <- function(X, parm) {
  n = length(X)
  FH = pmaxwell(X, parm)
  max(abs(c((1:n)/n - FH, FH - (0:(n-1))/n)))
}

CVn <- function(X, parm ) {
  n = length(X)
  FH = pmaxwell(X, parm)
  sum((FH - (2 * (1:n) - 1)/(2 * n))^2) + 1 / (12 * n)
}

ADn <- function(X, parm) {
  n = length(X)
  z = pmaxwell(X, parm)
  -n - mean((2*(1:n)-1) * (log(z) + log(1 - z[n:1])))
}

MAn <- function(X, parm) {
  n = length(X)
  z = pmaxwell(X, parm)
  n/2 - 2 * sum(z) - sum((2 - (2*(1:n) - 1)/n) * log(1 - z))
}

watson_U2_stat <- function(u) {
  n <- length(u); u_sorted <- sort(u)
  sum((u_sorted - (2*(1:n) - 1)/(2*n))^2) + 1 / (12*n)
}

Wn <- function(X, parm) watson_U2_stat(pmaxwell(X, parm))

ZAn <- function(X, parm, epsilon = 1e-4) {
  n = length(X)
  X[X == min(X)] <- X[X == min(X)] + epsilon
  j = 1:n
  FH =pmaxwell(X, parm)
  out = -sum(log(FH) / (n - j + 0.5) + log(1 - FH) / (j - 0.5))
  return(out)
}

ZBn <- function(X, parm,  epsilon = 1e-4) {
  n = length(X)
  X[X == min(X)] <- X[X == min(X)] + epsilon
  j = 1:n
  FH = pmaxwell(X, parm)
  arg = (1 / FH - 1) / ((n - 0.5) / (j - 0.75) - 1)
  out = sum((log(arg))^2)
  return(out)
}

ZCn <- function(X, parm) {
  n = length(X)
  j = 1:n
  FH = pmaxwell(X, parm)
  T1 = n * (j - 0.5) / (n - j + 0.5)^2 * log((j - 0.5) / (n * FH))
  T2 = n / (n - j + 0.5) * log((n - j + 0.5) / (n * (1 - FH)))
  out = 2 * sum(T1 + T2)
  return(out)
}
###############################################################################
# ========== Alternative Distributions ========= #
rEV <- function(n,parm){
  X = log(1-parm*log(1-runif(n)))
  return(X)
}
rLFR <- function(n,parm){
  U = runif(n)
  X = (-1+sqrt(1-2*parm*log(1-U)))/parm
  return(X)
}
rBE <- function(n,parm){
  B = rbeta(n,parm,1)
  X = -log(1-B)
  return(X)
}
rInvBeta <- function(n,parm){
  U = runif(n)
  X = 1/(1-U^(1/(parm+1)))
  return(X)
}

rTiltedPareto <- function(n,parm){
  out = (1+parm)/(1-runif(n))-parm
  return(out)
}
rBenini <- function(n,parm){
  U = runif(n)
  X = exp((-1+sqrt(1-4*parm*log(1-U)))/(2*parm))
  return(X)
}
dBenini <- function(x,parm){
  T1  = exp(-parm*(log(x))^2)/x^2
  T2  = 1+2*parm*log(x)
  out = T1*T2
  return(out)
}
# Chen's:
rChen <- function(n,parm){
  X = log(1-0.5*log(1-runif(n)))^(1/parm)
  return(X)
}

# Extreme Value:
rEV <- function(n,parm){
  X = log(1-parm*log(1-runif(n)))
  return(X)
}
###############################################################
#   Simulation Settings
###############################################################

n = 50 
MC = 100000
a_values = c(0.5, 1, 1.5)
nu=3
###############################################################
#   Alternatives  
###############################################################

alternatives <- list(  
  max1 = function(n) rmaxwell(n,1),
  max2 = function(n) rmaxwell(n,3),
  max3 = function(n) rmaxwell(n,5),
  
  Halfnormal_1 = function(n) rhalfnormal(n, 1),
  Halfnormal_3 = function(n) rhalfnormal(n, 3),
  
  Ray1 = function(n) rrayleigh(n, 0.5),
  Ray2 = function(n) rrayleigh(n, 1),
  Ray3 = function(n) rrayleigh(n, 3),
  
  Pareto_1  = function(n) rpareto(n, 1, 0.5),
  Pareto_2  = function(n) rpareto(n, 2, 1),
  
  Levy_1 = function(n) rlevy(n, 0, 4),
  Levy_2 = function(n) rlevy(n, 1, 2),
  
  Ev1=function(n) rEV(n,1),
  Ev2=function(n) rEV(n,2),
  
  Ch1  = function(n) rChen(n, 1),
  Ch2  = function(n) rChen(n, 1.5),
  
  LF1=function(n) rLFR(n, 1),
  LF2=function(n) rLFR(n, 2),
  
  ChiSq_1 = function(n) rchisq(n, 1),
  ChiSq_2 = function(n) rchisq(n, 3),
  ChiSq_3 = function(n) rchisq(n, 5),
  
  Gamma_1 = function(n) rgamma(n, 0.5, 1),
  Gamma_2 = function(n) rgamma(n, 1, 1),
  Gamma_3 = function(n) rgamma(n, 1, 3),
  Gamma_4 = function(n) rgamma(n, 2, 5),
  
  InvGauss_1 = function(n) rinvgauss(n, 2, 1),
  InvGauss_2 = function(n) rinvgauss(n, 2, 3),
  InvGauss_3 = function(n) rinvgauss(n, 3, 2),
  
  InvGamma_1 = function(n) rinvgamma(n, 1.5, 1),
  InvGamma_2 = function(n) rinvgamma(n, 3, 1),
  InvGamma_3 = function(n) rinvgamma(n, 3.0, 5),
  
  Weibull_1 = function(n) rweibull(n, 0.7, 1),
  Weibull_2 = function(n) rweibull(n, 1, 1),
  Weibull_3 = function(n) rweibull(n, 2.0, 1) 
)

###############################################################
#   Compute Critical Values (Null)
###############################################################
critical <- matrix(0, nrow = 1, ncol = length(a_values))
colnames(critical) <- paste0("a=", a_values)

for (k in 1:length(a_values)) {
  a = a_values[k]
  stats = numeric(MC)
  
  for (i in 1:MC) {
    x = rmaxwell(n, 1)  #null
    parm = estimation(x)
    stats[i] = Tn(x, parm, nu, a)
  }
  critical[k] = quantile(stats, 0.95)
}

#critical
###############################################################
#   POWER STUDY
###############################################################
results = matrix(0, nrow = length(alternatives), ncol = length(a_values))
rownames(results) = names(alternatives)
colnames(results) = paste0("a=", a_values)

for (alt in names(alternatives)) {
  for (k in 1:length(a_values)) {
    a = a_values[k]
    cv = critical[k]
    power = numeric(MC)
 for (i in 1:MC) {
      x = alternatives[[alt]](n)  #alternatives
      parm = estimation(x)
      power[i] = (Tn(x, parm, nu, a) > cv)
    }
    results[alt, k] = mean(power)
  }
}
###############################################################
#   Final Power Table
###############################################################
#results
###############################################################
#   Classical Tests: Critical Values (Null)
##############################################################
classic_tests <- list(
  AE  = AEn,
  KS  = KSn,
  CVM = CVn,
  AD  = ADn,
  MA  = MAn,
  W   = Wn,
  ZA  = ZAn,
  ZB  = ZBn,
  ZC  = ZCn
)

classic_critical <- matrix(0, nrow = length(classic_tests), ncol = 1)
rownames(classic_critical) <- names(classic_tests)
colnames(classic_critical) <- "CV_0.95"

for (tst in names(classic_tests)) {
  stat = numeric(MC)
  
  for (i in 1:MC) {
    x = rmaxwell(n, 1)
    parm = estimation(x)
    stat[i] = classic_tests[[tst]](x, parm)
  }
  
  classic_critical[tst, 1] = quantile(stat, 0.95)
}
#classic_critical
###############################################################
#   Classical Tests: Power Study
###############################################################
classic_power <- array(0, dim = c(length(alternatives), length(classic_tests)),
  dimnames = list(names(alternatives), names(classic_tests))
)

for (alt in names(alternatives)) {
  for (tst in names(classic_tests)) {
    
    cv = classic_critical[tst, 1]
    reject = numeric(MC)
    
    for (i in 1:MC) {
      x = alternatives[[alt]](n)
      parm = estimation(x)
      reject[i] = (classic_tests[[tst]](x, parm) > cv)
    }
    classic_power[alt, tst] = round(mean(reject),3)
  }
}
#classic_power
###############################################################
#   Combined Power Table
###############################################################
combined_results <- cbind(results, classic_power)
combined_results









