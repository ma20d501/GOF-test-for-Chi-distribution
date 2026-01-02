# ======================================================================================================
# ================ Unified Simulation for Chi-family Tests (Half-normal, Rayleigh, Maxwell)  ===========
# ======================================================================================================
###############################################################
#                 Test for Half-normal Distribution (nu = 1)(Warp-speed bootstrap)
#                       All GOF Tests Included
#                       n = 25,  a = 0.5,1,1.5
#============================================================

library(stats4); library(VGAM); library(extraDistr); library(actuar); library(evd)

# ========== Target Distribution ========== #
rchi_family <- function(n, nu, sigma = 1) {
  return(sigma * sqrt(rchisq(n, df = nu)))
}

pchi_family <- function(x, nu, sigma = 1) {
  return(pchisq((x / sigma)^2, df = nu))
}
rhalfnormal <- function(n, sigma) abs(rnorm(n, 0, sigma))

phalfnormal <- function(x, sigma) {
  erf <- function(z) 2 * pnorm(z * sqrt(2)) - 1
  erf(x / (sigma * sqrt(2)))
}
# ========== Parameter Estimation ========== #
estimation <- function(x) sqrt(mean(x^2))

# ========== Test Statistics ========== #
Tn <- function(X, parm, a, nu) {
  data = X / parm
  n = length(data)
  mdata = matrix(data, n, n)
  pdata = mdata + t(mdata) + a
  multdata = (nu - 1 - mdata^2) * t(nu - 1 - mdata^2)
  SUM = 2 / pdata^3 - 2 * (nu - 1 - mdata^2) / (mdata * pdata^2) + multdata / (mdata * t(mdata) * pdata)
  sum(SUM) / n
}

KSn <- function(X, parm, nu) {
  n = length(X)
  FH = pchi_family(X, nu, parm)
  max(abs(c((1:n)/n - FH, FH - (0:(n-1))/n)))
}

CVn <- function(X, parm, nu) {
  n = length(X)
  FH = pchi_family(X, nu, parm)
  sum((FH - (2 * (1:n) - 1)/(2 * n))^2) + 1 / (12 * n)
}

ADn <- function(X, parm, nu) {
  n = length(X)
  z = pchi_family(X, nu, parm)
  -n - mean((2*(1:n)-1) * (log(z) + log(1 - z[n:1])))
}

MAn <- function(X, parm, nu) {
  n = length(X)
  z = pchi_family(X, nu, parm)
  n/2 - 2 * sum(z) - sum((2 - (2*(1:n) - 1)/n) * log(1 - z))
}

watson_U2_stat <- function(u) {
  n <- length(u); u_sorted <- sort(u)
  sum((u_sorted - (2*(1:n) - 1)/(2*n))^2) + 1 / (12*n)
}

Wn <- function(X, parm, nu) watson_U2_stat(pchi_family(X, nu, parm))

ZAn <- function(X, parm,nu, epsilon = 1e-4) {
  n = length(X)
  X[X == min(X)] <- X[X == min(X)] + epsilon
  j = 1:n
  FH =pchi_family(X, nu, parm)
  out = -sum(log(FH) / (n - j + 0.5) + log(1 - FH) / (j - 0.5))
  return(out)
}

ZBn <- function(X, parm, nu, epsilon = 1e-4) {
  n = length(X)
  X[X == min(X)] <- X[X == min(X)] + epsilon
  j = 1:n
  FH = pchi_family(X, nu, parm)
  arg = (1 / FH - 1) / ((n - 0.5) / (j - 0.75) - 1)
  out = sum((log(arg))^2)
  return(out)
}

ZCn <- function(X, parm, nu) {
  n = length(X)
  j = 1:n
  FH = pchi_family(X, nu, parm)
  T1 = n * (j - 0.5) / (n - j + 0.5)^2 * log((j - 0.5) / (n * FH))
  T2 = n / (n - j + 0.5) * log((n - j + 0.5) / (n * (1 - FH)))
  out = 2 * sum(T1 + T2)
  return(out)
}

###############################################################################
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
# ========== Alternative Distributions ========== #
r_alternatives <- list(
  Halfnormal_0.5 = function(n) rhalfnormal(n, 0.5),
  Halfnormal_1 = function(n) rhalfnormal(n, 1),
  HalfNormal_3 = function(n) rhalfnormal(n, 3),
  
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
  

# ========== Simulation Parameters ========== #
#set.seed(1234)
n_values <- c(25, 50)
MC <- 100000 
nu_values <- 1 #c(1, 2, 3)  # Half-normal, Rayleigh, Maxwell
a_values <- c(0.5, 1, 1.5)

# ========== Simulation Loop ========== #
results <- list()

for (nu in nu_values) {
  for (n in n_values) {
    #cat("Running for nu =", nu, "n =", n, "\n")
    mac <- matrix(0, nrow = length(r_alternatives), ncol = 11)
    rownames(mac) <- names(r_alternatives)
    colnames(mac) <- c("Tn(0.5)", "Tn(1)","Tn(1.5)",  "KSn", "CVn", "ADn", "MAn", 
                       "Wn", "ZAn","ZBn","ZCn" )
    
    for (alt in names(r_alternatives)) {
      s.mat <- matrix(0, nrow = MC, ncol = 11)
      s.boot <- matrix(0, nrow = MC, ncol = 11)
      
      for (i in 1:MC) {
        x <- r_alternatives[[alt]](n)
        x.boot <- rchi_family(n, nu)
        
        mu <- estimation(x)
        mu.boot <- estimation(x.boot)
        
        s.mat[i,1] <- Tn(x, mu, nu, a_values[1])
        s.boot[i,1] <- Tn(x.boot, mu.boot, nu, a_values[1])
        
        s.mat[i,2] <- Tn(x, mu, nu, a_values[2])
        s.boot[i,2] <- Tn(x.boot, mu.boot, nu, a_values[2])
        
        s.mat[i,3] <- Tn(x, mu, nu, a_values[3])
        s.boot[i,3] <- Tn(x.boot, mu.boot, nu, a_values[3])
        
      
        s.mat[i,4] <- KSn(x, mu, nu)
        s.boot[i,4] <- KSn(x.boot, mu.boot, nu)
        
        s.mat[i,5] <- CVn(x, mu, nu)
        s.boot[i,5] <- CVn(x.boot, mu.boot, nu)
        
        s.mat[i,6] <- ADn(x, mu, nu)
        s.boot[i,6] <- ADn(x.boot, mu.boot, nu)
        
        s.mat[i,7] <- MAn(x, mu, nu)
        s.boot[i,7] <- MAn(x.boot, mu.boot, nu)
        
        s.mat[i,8] <- Wn(x, mu, nu)
        s.boot[i,8] <- Wn(x.boot, mu.boot, nu)
        
        s.mat[i,9] <- ZAn(x, mu, nu)
        s.boot[i,9] <- ZAn(x.boot, mu.boot, nu)
        
        s.mat[i,10] <- ZBn(x, mu, nu)
        s.boot[i,10] <- ZBn(x.boot, mu.boot, nu)
        
        s.mat[i,11] <- ZCn(x, mu, nu)
        s.boot[i,11] <- ZCn(x.boot, mu.boot, nu)
      }
      
      # Compute empirical power
      for (j in 1:11)
        mac[alt, j] <- mean(s.mat[, j] > quantile(s.boot[, j], 0.95), na.rm = TRUE)
    }
    
    results[[paste0("nu", nu, "_n", n)]] <- round(mac, 3)
  }
}
results[["nu1_n25"]]  # Half-normal case, n=25
results[["nu1_n50"]]  # Half-normal case, n=50
 
