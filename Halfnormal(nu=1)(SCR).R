###############################################################
#                 Test for Half-Normal Distribution (nu = 1)
#                       All GOF Tests Included
#                       n = 25,  a = 0.5,1,1.5
###############################################################
library(stats4)
library(VGAM)
library(stargazer)
library(extraDistr) 
library(EnvStats)   
library(statmod)
###############################################################
#   Half-normal  
###############################################################
rhalfnormal <- function(n, sigma) abs(rnorm(n, 0, sigma))

phalfnormal <- function(x, sigma) {
  erf <- function(z) 2 * pnorm(z * sqrt(2)) - 1
  erf(x / (sigma * sqrt(2)))
}
# Unified CDF for nu = 1 (Half-normal)
 pchi_family <- function(x, parm, nu = 1) phalfnormal(x, parm)

# CDF
#pchi_family <- function(x, parm, nu=1) {
#  pchisq((x / parm)^2, df = nu)
#}
###############################################################
#   Parameter Estimation (MLE)
###############################################################
estimation <- function(x) sqrt(sum(x^2) / (length(x)))
###############################################################
#   Tn Test Statistic
###############################################################
Tn <- function(X, parm, nu = 1, a) {
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
#   Comparison Tests (KS, CVM, AD, MA, Watson, ZA, ZB, ZC)
###############################################################

KSn <- function(X, parm, nu = 1) {
  n = length(X)
  FH = pchi_family(X, parm, nu)
  max(max((1:n)/n - FH), max(FH - (0:(n - 1))/n))
}

CVn <- function(X, parm, nu = 1) {
  n = length(X)
  FH = pchi_family(X, parm, nu)
  sum((FH - (2*(1:n)-1)/(2*n))^2) + 1/(12*n)
}

ADn <- function(X, parm, nu = 1) {
  n = length(X)
  z = pchi_family(X, parm, nu)
  T1 = 2*(1:n) - 1
  T2 = log(z) + log(1 - z[n:1])
  -n - mean(T1 * T2)
}

MAn <- function(X, parm, nu = 1) {
  n = length(X)
  z = pchi_family(X, parm, nu)
  S1 = sum(z)
  T1 = 2 - (2*(1:n)-1)/n
  T2 = log(1 - z)
  S2 = sum(T1 * T2)
  n/2 - 2*S1 - S2
}

Wn <- function(X, parm, nu = 1) {
  u = pchi_family(X, parm, nu)
  n = length(u)
  u_sorted = sort(u)
  i = 1:n
  Ui = (2*i - 1)/(2*n)
  sum((u_sorted - Ui)^2) + 1/(12*n)
}

ZAn <- function(X, parm, nu = 1, eps = 1e-4) {
  n = length(X)
  X[X == min(X)] <- X[X == min(X)] + eps
  j = 1:n
  FH = pchi_family(X, parm, nu)
  -sum(log(FH)/(n - j + 0.5) + log(1 - FH)/(j - 0.5))
}

ZBn <- function(X, parm, nu = 1, eps = 1e-4) {
  n = length(X)
  X[X == min(X)] <- X[X == min(X)] + eps
  j = 1:n
  FH = pchi_family(X, parm, nu)
  arg = (1/FH - 1) / ((n - 0.5)/(j - 0.75) - 1)
  sum(log(arg)^2)
}

ZCn <- function(X, parm, nu = 1) {
  n = length(X)
  j = 1:n
  FH = pchi_family(X, parm, nu)
  T1 = n*(j - 0.5)/(n - j + 0.5)^2 * log((j - 0.5)/(n*FH))
  T2 = n/(n - j + 0.5) * log((n - j + 0.5)/(n*(1 - FH)))
  2 * sum(T1 + T2)
}

###############################################################
#   Simulation Settings
###############################################################

n = 50                       #sample size
MC = 100000                  #replications
a_values = c(0.5, 1, 1.5)    #tuning parameters

###############################################################
#   Alternatives  
###############################################################

alternatives <- list(  
  
  Halfnormal_0.5 = function(n) rhalfnormal(n, 0.5),
  Halfnormal_1 = function(n) rhalfnormal(n, 1),
  Halfnormal_3 = function(n) rhalfnormal(n, 3),
  Halfnormal_5 = function(n) rhalfnormal(n, 5),
  
  Pareto_1  = function(n) rpareto(n, 1, 0.5),
  Pareto_2  = function(n) rpareto(n, 1, 1),
  Pareto_3  = function(n) rpareto(n, 1, 2),
  
  Levy_1 = function(n) rlevy(n, 0, 1),
  Levy_2 = function(n) rlevy(n, 1, 2),
  Levy_3 = function(n) rlevy(n, 2, 5),
  
  ChiSq_1 = function(n) rchisq(n, 1),
  ChiSq_2 = function(n) rchisq(n, 3),
  ChiSq_3 = function(n) rchisq(n, 5),
  
  Gamma_1 = function(n) rgamma(n, 0.5, 1),
  Gamma_2 = function(n) rgamma(n, 1, 1),
  Gamma_3 = function(n) rgamma(n, 1, 3),
  
  InvGauss_1 = function(n) rinvgauss(n, 2, 1),
  InvGauss_2 = function(n) rinvgauss(n, 2, 3),
  InvGauss_3 = function(n) rinvgauss(n, 3, 2),
  
  InvGamma_1 = function(n) rinvgamma(n, 1.5, 1),
  InvGamma_2 = function(n) rinvgamma(n, 3, 1),
  InvGamma_3 = function(n) rinvgamma(n, 3.0, 5),
  
  Weibull_1 = function(n) rweibull(n, 0.7, 1),
  Weibull_2 = function(n) rweibull(n, 1, 1),
  Weibull_3 = function(n) rweibull(n, 2.0, 0.5) 
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
    x = rhalfnormal(n, 1)
    parm = estimation(x)
    stats[i] = Tn(x, parm, nu = 1, a)
  }
  critical[k] = quantile(stats, 0.95)
}
critical
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
      x = alternatives[[alt]](n)
      parm = estimation(x)
      power[i] = (Tn(x, parm, nu = 1, a) > cv)
    }
    
    results[alt, k] = mean(power)
  }
}
#results
###############################################################
#   Classical Tests: Critical Values (Null)
###############################################################
classic_tests <- list(
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
    x = rhalfnormal(n, 1)
    parm = estimation(x)
    stat[i] = classic_tests[[tst]](x, parm, nu = 1)
  }
  
  classic_critical[tst, 1] = quantile(stat, 0.95)
}

###############################################################
#   Classical Tests: Power Study
###############################################################
classic_power <- array(
  0,
  dim = c(length(alternatives), length(classic_tests)),
  dimnames = list(names(alternatives), names(classic_tests))
)

for (alt in names(alternatives)) {
  
  for (tst in names(classic_tests)) {
    
    cv = classic_critical[tst, 1]
    reject = numeric(MC)
    
    for (i in 1:MC) {
      x = alternatives[[alt]](n)
      parm = estimation(x)
      reject[i] = (classic_tests[[tst]](x, parm, nu = 1) > cv)
    }
    
    classic_power[alt, tst] = mean(reject)
  }
}
###############################################################
#   Combined Power Table
###############################################################
combined_results <- cbind(results, classic_power)
combined_results









