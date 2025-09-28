# install package ---------------------------------------------------------
library(mvtnorm)
library(piqp)
library(SHT)
library(dbEmpLikeGOF)
library(eummd)
library(reticulate)

# main func ---------------------------------------------------------------
phi <- function(z) {
  sqrt(z)/2
}

T_stat = function(X, Y, KK=2) {
  result = rep(0, KK)
  D = as.matrix(dist(rbind(X,Y), diag=TRUE, upper=TRUE))
  for (k in 1:KK) {
    Dk = D^k
    temp = phi(Dk)
    dXX = temp[1:m, 1:m]
    dXY = temp[1:m, (m+1):(m+n)]
    dYY = temp[(m+1):(m+n), (m+1):(m+n)]
    result[k] = 2*sum(dXY)/(m*n) - sum(dXX)/(m^2) - sum(dYY)/(n^2)
  }
  return(result)
}

A_st <- function(X, Y, s, t) {  
  D <- as.matrix(dist(rbind(X, Y), diag = TRUE, upper = TRUE))
  Ds = D^s
  DXXs <- Ds[1:m, 1:m]
  DXYs <- Ds[1:m, (m+1):(m+n)]
  DYYs <- Ds[(m+1):(m+n), (m+1):(m+n)]
  Dt = D^t
  DXXt <- Dt[1:m, 1:m]
  DXYt <- Dt[1:m, (m+1):(m+n)]
  DYYt <- Dt[(m+1):(m+n), (m+1):(m+n)]
  com_s <- sum(phi(DYYs)) / (n^2)
  com_t <- sum(phi(DYYt)) / (n^2)
  sum_XY_s_row <- rowSums(phi(DXYs)) / n
  sum_XY_t_row <- rowSums(phi(DXYt)) / n
  a1 <- outer(sum_XY_s_row, rep(1, m)) + outer(rep(1, m), sum_XY_s_row) - 
    phi(DXXs/hs) - com_s
  a2 <- outer(sum_XY_t_row, rep(1, m)) + outer(rep(1, m), sum_XY_t_row) - 
    phi(DXXt/ht) - com_t
  sum(a1 * a2) / (m^2)
}

A_hat = function(X, Y, KK=2) {
  Re = matrix(0, KK, KK)
  for (s in 1:KK) {
    for (t in 1:KK) {
      Re[s,t] = A_st(X, Y, s, t)
    }
  }
  return(Re)
}

Tw = function(X, Y, KK=2) {
  res = solve_piqp(P=2*A_hat(X, Y, KK), c=rep(0,KK), 
                   A=matrix(rep(1,KK),ncol=KK), b=1,
                   G=matrix(rep(0,KK*KK),ncol=KK), h=rep(0,KK), x_lb=rep(0,KK))                  
  w_hat = res$x
  re = c(w_hat) %*% T_stat(X, Y, KK)
  return(c(re, w_hat))
}

T_Zhou2017smooth <- function(X, Y, Es, Up) {
  m = nrow(X)  
  n = nrow(Y)
  d = ncol(X)
  T0 = rep(0, len)
  for (i in 1:len) {
    u <- matrix(Up[i,], ncol = 1)  
    X_proj <- c(X %*% u)  
    Y_proj <- c(Y %*% u)  
    Z <- colMeans(outer(X_proj, Y_proj, "<=")) 
    Wavelet_Z <- sqrt(2) * cos(outer(pi * Z, 1:dd))
    Phi1 <- colSums(Wavelet_Z)
    T0[i] <- max(abs(Phi1)) / sqrt(n)
  }
  stats <- max(T0) * sqrt(m / (n + m))
  T1 <- numeric(B)
  for (k in 1:B) {
    e <- Es[,k]  
    Tk = rep(0, len)
    for (i in 1:len) {
      u <- matrix(Up[i,], ncol = 1)
      X_proj <- c(X %*% u)
      Z <- colMeans(outer(X_proj, X_proj, "<="))
      Wavelet_Z <- sqrt(2) * cos(outer(pi * Z, 1:dd))
      Phi1 <- colSums(e * Wavelet_Z) 
      Tk[i] <- max(abs(Phi1)) / sqrt(m)
    }
    T1[k] <- max(Tk)
  }
  critival = quantile(T1, probs=1-alpha)
  return( c( pow = (stats >= unname(critival)),
             p = (1 + sum(T1 >= stats)) / (B + 1),
             stats = stats, T1 = T1, critival = critival ) )
}


# DGP ---------------------------------------------------------------------
Y_ZhouEx1 <- function(n, mu) {
  if (mu == 0) return(runif(n, -1, 1))
  samples <- numeric(n)
  accepted <- logical(n)
  while (sum(accepted) < n) {
    needed <- n - sum(accepted)
    u <- runif(needed, -1, 1)
    g <- ifelse(abs(u) < mu, 0.5 + 2*u*(mu - abs(u))/mu^2, 0.5)
    max_g <- 1  
    accept_new <- runif(needed) < g / max_g 
    samples[!accepted] <- u
    accepted[!accepted] <- accept_new
  }
  return(samples)
}


# implementation ----------------------------------------------------------
m = 20
n = 20
d = 1
Sim = 1000
R = 299 
alpha = 0.05
KK = 2
theta = 0.5
X = matrix(runif(m*d, -1, 1), nrow=m, ncol=d)
Y = matrix(Y_ZhouEx1(n*d, mu = theta), nrow=n, ncol=d)
Z = rbind(X, Y)
ind = 1:(m + n)

# Composite test 
Tw0 = Tw(X, Y, KK)
T0.test = Tw0[1]
w0 = Tw0[-1]
TR.test = rep(0,R)
wR = matrix(0, KK, R)
for (b in 1:R) {
  k <- sample(ind, size = m, replace = FALSE) 
  X1 <- matrix(Z[k, ], ncol = d)
  Y1 <- matrix(Z[-k, ], ncol = d)
  TwR = Tw(X1, Y1, KK)
  TR.test[b] = TwR[1]
  wR[, b] = TwR[-1]
}
p.w = mean(c(T0.test, TR.test) >= T0.test)
pow.w = 1 * (p.w <= alpha)

# Zhou etal.(2017) 
dd = 4 
len = 10 
Up = matrix(rnorm(len * d), ncol = d) 
Up = Up / sqrt(rowSums(Up^2)) 
B = R+1 
Es = matrix(rnorm(m * B), ncol = B)
rezhou = T_Zhou2017smooth(X, Y, Es, Up)
p.Zhou = unname(rezhou["p"])
pow.Zhou = unname(rezhou["pow"])

# Biswas and Ghosh(2014)
if (ncol(X) == 1) {
  p.BG = eqdist.2014BG(c(X), c(Y))$p.value
  pow.BG = 1 * (p.BG <= alpha)
} else {
  p.BG = eqdist.2014BG(X, Y)$p.value
  pow.BG = 1 * (p.BG <= alpha)
}

# Gurevich and Vexler(2011)
p.GV = dbEmpLikeGOF(X, Y, delta=0.1)$pvalue
pow.GV = 1 * (p.GV <= alpha)

# Schrab etal.(2023)
conda_create(envname = "mmdagg-env", python_version = "3.9", 
             packages = c("numpy==1.26.4", "scipy")) 
py_install(c("jax==0.4.23", "jaxlib==0.4.23"), envname = "mmdagg-env", pip = TRUE)
py_install("C:/Users/YOURUSERNAME/Desktop/mmdagg-master", envname = "mmdagg-env", pip = TRUE) 
use_condaenv("mmdagg-env") 
X_np <- np$array(X, dtype = "float32")
Y_np <- np$array(Y, dtype = "float32")
re.Schrab = mmdagg$mmdagg(X_np, Y_np, return_dictionary=TRUE) 
pval_all0 = c(sapply(paste0("Single test 1.", 1:10), function(n) re.Schrab[[2]][[n]][["p-value"]]),
              sapply(paste0("Single test 2.", 1:10), function(n) re.Schrab[[2]][[n]][["p-value"]]))
pval_all = sapply(pval_all0, function(x) x$item())
p.Schrab = min(unname(pval_all)) 
threshold = re.Schrab[[2]][["Single test 1.1"]][["p-value threshold"]]$item()
pow.Schrab = re.Schrab[[1]]$item()

# Ramdas(2015)
p.MMD = eummd::mmd(X, Y, beta = -0.1, kernel = "Gaussian")$pval
pow.MMD = 1 * (p.MMD <= alpha)
band = as.matrix(dist(rbind(X, Y), method = "euclidean", diag = TRUE, upper = TRUE))
band05 = 1/(median(band[band>0])*d^(0.5))^2
p.Ramdas05 = eummd::mmd(X, Y, beta = band05, kernel = "Gaussian")$pval
pow.Ramdas05 = 1 * (p.Ramdas05 <= alpha)
band075 = 1/(median(band[band>0])*d^(0.75))^2
p.Ramdas075 = eummd::mmd(X, Y, beta = band075, kernel = "Gaussian")$pval
pow.Ramdas075 = 1 * (p.Ramdas075 <= alpha)                                        
                                        


