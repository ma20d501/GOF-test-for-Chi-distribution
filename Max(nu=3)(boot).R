# ======================================================================================================
# ================ Unified Simulation for Chi-family Tests (Half-normal, Rayleigh, Maxwell)  ===========
# ======================================================================================================
###############################################################
#                 Test for Maxwell Distribution (nu = 3)(Warp-speed bootstrap)
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
# ========== chi Distribution ==========
rchi_family <- function(n, nu, sigma = 1) {
  sigma * sqrt(rchisq(n, df = nu))
}

# Unified CDF
pchi_family <- function(x, parm, nu) {
  pchisq((x / parm)^2, df = nu)
}
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

# ========== Alternative Distributions ========= #
r_alternatives <- list(  
  max1 = function(n) rmaxwell(n,1),
  #max2 = function(n) rmaxwell(n,3),
  #max3 = function(n) rmaxwell(n,5),
  
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

# ========== Simulation Parameters ========== #
n_values <- c(25,50)
MC <- 100000      
nu_values <- 3  #  Maxwell   
a_values <- c(0.5, 1, 1.5)
  
# ========== Simulation Loop ========== #
results <- list()

for (nu in nu_values) {
  for (n in n_values) {
    #cat("Running for nu =", nu, "n =", n, "\n")
    mac <- matrix(0, nrow = length(r_alternatives), ncol = 12)
    rownames(mac) <- names(r_alternatives)
    colnames(mac) <- c("Tn(0.5)", "Tn(1)","Tn(1.5)", "AEn", 
                       "KSn", "CVn", "ADn", "MAn", "Wn", "ZAn","ZBn","ZCn" )
    
    for (alt in names(r_alternatives)) {
      s.mat <- matrix(0, nrow = MC, ncol = 12)
      s.boot <- matrix(0, nrow = MC, ncol = 12)
      
      for (i in 1:MC) {
        x <- r_alternatives[[alt]](n)
        x.boot <- rmaxwell(n, 1)   #rchi_family(n, nu)
        
        mu <- estimation(x)
        mu.boot <- estimation(x.boot)
        
        s.mat[i,1] <- Tn(x, mu, nu, a_values[1])
        s.boot[i,1] <- Tn(x.boot, mu.boot, nu, a_values[1])
        
        s.mat[i,2] <- Tn(x, mu, nu, a_values[2])
        s.boot[i,2] <- Tn(x.boot, mu.boot, nu, a_values[2])
        
        s.mat[i,3] <- Tn(x, mu, nu, a_values[3])
        s.boot[i,3] <- Tn(x.boot, mu.boot, nu, a_values[3])
        
        
        s.mat[i,4] <- AEn(x, mu)
        s.boot[i,4] <- AEn(x.boot, mu.boot)
        
        
        s.mat[i,5]<- KSn(x, mu)
        s.boot[i,5]<- KSn(x.boot,mu.boot)
        
        s.mat[i,6]<- CVn(x, mu)
        s.boot[i,6]<-CVn(x.boot,mu.boot)
        
        s.mat[i,7]<- ADn(x, mu)
        s.boot[i,7]<-ADn(x.boot,mu.boot)
        
        s.mat[i,8]<- MAn(x, mu)
        s.boot[i,8]<- MAn(x.boot,mu.boot)
        
        s.mat[i,9]<- Wn(x, mu)
        s.boot[i,9]<- Wn(x.boot,mu.boot)   
        
        s.mat[i,10]<- ZAn(x,mu)
        s.boot[i,10]<- ZAn(x.boot,mu.boot)    
        
        s.mat[i,11]<- ZBn(x,mu)
        s.boot[i,11]<-  ZBn(x.boot,mu.boot)
        
        s.mat[i,12]<- ZCn(x,mu)
        s.boot[i,12]<- ZCn(x.boot,mu.boot) 
        
      }
      
      # Compute empirical power
      for (j in 1:12)
        mac[alt, j] <- mean(s.mat[, j] > quantile(s.boot[, j], 0.95), na.rm = TRUE)
    }
    
    results[[paste0("nu", nu, "_n", n)]] <- round(mac, 3)
  }
}
results[["nu3_n25"]]   
results[["nu3_n50"]]   


