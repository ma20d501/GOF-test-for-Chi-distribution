###############################################################
#                 Test for Rayleigh Distribution (nu = 2)
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
# ========== Target Distribution ========== #
rchi_family <- function(n, nu, sigma = 1) {
  return(sigma * sqrt(rchisq(n, df = nu)))
}

pchi_family <- function(x, nu, sigma = 1) {
  return(pchisq((x / sigma)^2, df = nu))
}
# Random number generation for Rayleigh distribution
rrayleigh <- function(n, parm) {
  return(sqrt(-2 * parm^2 * log(runif(n))))
}

# Cumulative distribution function for Rayleigh distribution
prayleigh <- function(x, parm) {
  return(1 - exp(- (x / parm)^2 / 2))
}

###############################################################
#   Parameter Estimation (MLE)
###############################################################
nu=2
estimation <- function(x) sqrt(mean(x^2) / nu)

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
#   Comparison Tests  ("AHn", "MIn", "BKn", "JHn", "PHn", "LAn")
###############################################################

tqt = function(n, x, parm) {
  s = 0; a = 0.5
  g = function(u) (-2*log(1-u))^(0.5*(1-a))
  # muhat = sqrt(sum(x^2)/(2*n))
  g = Vectorize(g)
  for(i in 1:n) {
    s = s + (x[i]^a)*(integrate(g, (i-1)/n, i/n)$value)
  }
  T02 = ((parm^(1-a))*s - (a*mean(x)) - (1-a)*parm*(sqrt(pi/2))) / (mean(x)*(a-1))
  return(T02)
}

AHn = function(x, parm){
  n <- length(x)
  #xx <- rrayleigh(n, parm)
  xx<-x/parm
  out = tqt(n, xx, parm)
  return(out)
}

## Meintanis and Iliopoulos test statistic
MIn <- function(X, parm) {
  n <- length(X)
  data = X / parm  
  a= 2*sqrt(2)
  # Compute the double summation term
  term1 <- 0
  for (j in 1:n) {
    for (k in 1:n) {
      denom <- data[j] + data[k] + a
      term1 <- term1 + (1 / denom) +
        (data[j] + data[k]) / denom^2 +
        2 * (data[j] * data[k] + 2) / denom^3 +
        6 / denom^4 +
        24 / denom^5
    }
  }
  term1 <- sqrt(2) / n * term1
  # Compute the second summation term
  term2 <- 0
  for (j in 1:n) {
    denom <- data[j] + a
    term2 <- term2 + (1 / denom) + (data[j] / denom^2) + (2 / denom^3)
  }
  term2 <- -2 * sqrt(2) * term2
  # Compute the L value
  L <- n / a + term1 + term2
  return(L)
}

## Baratpour and Khodadadi test statistic
BKn <- function(X, parm) {
  # Sort the data
  X1<- X / parm  
  X <- sort(X1)
  n <- length(X)
  if (n < 2) stop("The data must contain at least two elements.")
  # Calculate the mean of X
  X_bar <- mean(X)
  # Compute the first term (A1)
  A1 <- sum((n - 1:(n - 1)) / n * log((n - 1:(n - 1)) / n) * (X[2:n] - X[1:(n - 1)]))
  # Compute the second term (A2)
  A2 <- sqrt(pi / 2) * sqrt(sum(X^3) / (3 * sum(X)))
  # Compute CK_n
  CKn <- (A1 + A2) / X_bar
  return(CKn)
}

##########################################################
jrftest = function(x, parm){
  n <- length(x)
  xneg = rep(NA, n)
  xpos = rep(NA, n)
  #l_hat = sqrt(sum(x^2)/(2*n))
  Xdata = sort(x)
  m = n - 5
  for (k in 1:n){
    neg = k - m
    xneg[k] = ifelse(neg < 1, Xdata[1], Xdata[neg])
    pos = k + m
    xpos[k] = ifelse(pos > n, Xdata[n], Xdata[pos])
  }
  vas = 1 / ((n / (2 * m)) * (xpos - xneg))
  f02 = (Xdata / (parm^2)) * exp(-Xdata^2 / (2 * (parm^2)))
  DHn = (1 / (2 * n)) * sum((sqrt(vas) - sqrt(f02))^2 / vas)
  return(DHn)
}

JHn = function(x, parm){
  n <- length(x)
  #sigma_hat = sqrt(sum(x^2)/(2*n))
  #xx <- rrayleigh(n, parm)
  xx<-x/parm
  out = jrftest(xx, parm)
  return(out)
}

###############################################################
#   PH_n : Hellinger-distance test (Rayleigh)
###############################################################
PHn <- function(Y, sigma_hat, h = NULL) {
  n <- length(Y)
  # Rayleigh density
  g_param <- function(x, sigma) {
    (x / sigma^2) * exp(-x^2 / (2 * sigma^2))
  }
  gY <- g_param(Y, sigma_hat)
  # Nonparametric density (KDE)
  if (is.null(h)) {
    kde <- density(Y, from = min(Y), to = max(Y), n = 4096)
  } else {
    kde <- density(Y, bw = h, from = min(Y), to = max(Y), n = 4096)
  }
  ghat <- approx(kde$x, kde$y, xout = Y, rule = 2)$y
  # Numerical safeguard
  eps <- 1e-8
  ghat <- pmax(ghat, eps)
  PH <- (1 / (2 * n)) * sum((1 - sqrt(gY / ghat))^2)
  return(PH)
}

###############################################################
#   LA_n : Liebenberg & Allison (2019) (Rayleigh)
###############################################################
LAn <- function(Y, phi = 5) {
  n <- length(Y)
  # Precompute sums
  sumY <- sum(Y)
  sumYlogY <- sum(Y * log(Y))
  LA <- 0
  for (k in 1:n) {
    term1 <- Y[k] * sumYlogY
    term2 <- (log(2) + digamma(1 + Y[k]/2)) * Y[k] * sumY
    LA <- LA + (term1 - term2)^2 * exp(-phi * Y[k]^2)
  }
  return(LA)
}

###############################################################
#   Comparison Tests (KS, CVM, AD, MA, Watson, ZA, ZB, ZC)
###############################################################
KSn <- function(X, parm) {
  n = length(X)
  FH = prayleigh(X, parm)
  max(abs(c((1:n)/n - FH, FH - (0:(n-1))/n)))
}

CVn <- function(X, parm ) {
  n = length(X)
  FH = prayleigh(X, parm)
  sum((FH - (2 * (1:n) - 1)/(2 * n))^2) + 1 / (12 * n)
}

ADn <- function(X, parm) {
  n = length(X)
  z = prayleigh(X, parm)
  -n - mean((2*(1:n)-1) * (log(z) + log(1 - z[n:1])))
}

MAn <- function(X, parm) {
  n = length(X)
  z = prayleigh(X, parm)
  n/2 - 2 * sum(z) - sum((2 - (2*(1:n) - 1)/n) * log(1 - z))
}

watson_U2_stat <- function(u) {
  n <- length(u); u_sorted <- sort(u)
  sum((u_sorted - (2*(1:n) - 1)/(2*n))^2) + 1 / (12*n)
}

Wn <- function(X, parm) watson_U2_stat(prayleigh(X, parm))

ZAn <- function(X, parm, epsilon = 1e-4) {
  n = length(X)
  X[X == min(X)] <- X[X == min(X)] + epsilon
  j = 1:n
  FH =prayleigh(X, parm)
  out = -sum(log(FH) / (n - j + 0.5) + log(1 - FH) / (j - 0.5))
  return(out)
}

ZBn <- function(X, parm,  epsilon = 1e-4) {
  n = length(X)
  X[X == min(X)] <- X[X == min(X)] + epsilon
  j = 1:n
  FH = prayleigh(X, parm)
  arg = (1 / FH - 1) / ((n - 0.5) / (j - 0.75) - 1)
  out = sum((log(arg))^2)
  return(out)
}

ZCn <- function(X, parm) {
  n = length(X)
  j = 1:n
  FH = prayleigh(X, parm)
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

n = 25                       #sample size
MC = 100000                  #replications
a_values = c(0.5, 1, 1.5)    #tuning parameters

###############################################################
#   Alternatives (you can add/remove)
###############################################################

alternatives <- list(  
  
  #Ray1 = function(n) rrayleigh(n, 0.5),
  Ray2 = function(n) rrayleigh(n, 1),
  #Ray3 = function(n) rrayleigh(n, 3),
  
  Halfnormal_1 = function(n) rhalfnormal(n, 1),
  Halfnormal_3 = function(n) rhalfnormal(n, 3),
  
  #max1 = function(n) rmaxwell(n,1),
  #max2 = function(n) rmaxwell(n,3),
  #max3 = function(n) rmaxwell(n,5),
  
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
    x = rrayleigh(n, 1)
    parm = estimation(x)
    stats[i] = Tn(x, parm, nu = 2, a)
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
      power[i] = (Tn(x, parm, nu = 2, a) > cv)
    }
    
    results[alt, k] = round(mean(power),3)
  }
}
#results
###############################################################
#   Classical Tests: Critical Values (Null)
###############################################################
classic_tests <- list(
  AH  = AHn,
  MI   = MIn,
  BK  = BKn,
  JH  = JHn,
  PH  = PHn,
  LA  = LAn,
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
    x = rrayleigh(n, 1)
    parm = estimation(x)
    stat[i] = classic_tests[[tst]](x, parm)
  }
  
  classic_critical[tst, 1] = quantile(stat, 0.95)
}
#classic_critical
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










