#packages ------------------------------------------------------------

library(MASS)
library(zoo)
library(stats)
library(vars)
library(glmnet)
library(ggplot2)
library(tidyr)
library(dplyr)
library(latex2exp) #write expression with latex in ggplot
library(RSpectra)
library(gridExtra)
library(Matrix)

#defining function------------------------------------------------------
stab_test <- function(kp, A, tol = 1e-8)
{
  if (!is.matrix(A) || nrow(A) != ncol(A)) {
    stop("The matrix is not square")
  }
  eig <- eigen(A, only.values = TRUE)$values  # computing the eigenvalues
  
  for (i in 1:length(eig)) {     
    if (Mod(eig[i]) >= 1 - tol) { 
      return(FALSE)               
    }
  }
  return(TRUE)
}


#function to compute companion matrix
comp_mtrx <- function(AA,p ){ #AA is the cbind of covariate matrices
  K<-nrow(AA)
  Kp<-K*p
  C<-matrix(0, nrow=Kp, ncol=Kp)
  C[1:K,1:Kp]<-AA
  C[(K+1):Kp,1:(K*(p-1))]<- diag(1,nrow=K, ncol=K )
  return(C)
}

#
est_autocov <- function(y_t, Y, Z, T, p=1){
  K <- ncol(y_t)
  Kp <- K * p
  I_t <- diag(T)
  
  # QR decomposition of Z to avoid singularity issues
  qr_decomp <- qr(Z)
  Q <- qr.Q(qr_decomp)
  P_Z <- Q %*% t(Q)  # Projection matrix
  
  # Compute bias-corrected covariance
  bias_sigma <- 1/T * t(Y) %*% (I_t - P_Z) %*% Y
  
  # Degrees of freedom correction
  d.f. <- T / (T - Kp - 1)
  unbiased <- d.f. * bias_sigma  # Corrected covariance estimate
  
  return(unbiased)
}

## VAR(p) process Simulator
var_sim <- function(AA, nu, Sigma_u, nSteps, y0) {
  K <- nrow(Sigma_u)
  Kp <- ncol(AA)
  p <- Kp/K
  
  if (p > 1) {
    C <- comp_mtrx(AA) # form the  companion matrix of the var(p) process
  } else {
    C <- AA  
  }
  y_t <- matrix(0, nrow =  nSteps, ncol=Kp) #trajectories matrix nSteps x Kp
  y_t[1, 1:Kp] <- y0 #add initial value to initiate the simulation
  noise <- mvrnorm(n = nSteps, mu = rep(0, K), Sigma = Sigma_u) #assuming that 
  #residuals follow a multivariate normal distribution    
  
  for (t in 2:nSteps) {
    y_t[t, ] <- C %*% y_t[t-1, ]
    y_t[t, 1:K] <- y_t[t, 1:K] + nu + noise[t,]
  }
  
  y_t <- zoo(y_t[,1:K], 1:nSteps)  
  return(y_t)
}

#compute the largest eigenvalues of A^TA, used to compute l2 norm 

lmda_max<-function(A,k){#k is the number of largest eigenvalues to compute
  s1<-svds(A,k=k, nu=0, nv=0)$d[1]
  return(s1)
}

#we then define a function to scale off-diagonal elements to ensure that 
# ||A||=target
scaled_A_matrix<- function(A,target, tol=1e-6, max_k=1000){
  D<-diag(diag(A))
  B<-A-D
  #check the l2nrom of A, that is the max eigenvalue of A^TA
  #we use hte previousy defined function lmbda_max
  f<-function(k){
    M<-D + k*B
    lmda_max(M, k=1)-target
  }
  #we evaluate at k=0
  f0<-f(0)
  if(abs(f0)<tol){
    return(D)
  }
  k_low <- 0
  k_high <- 1
  fk_high <- f(k_high)
  iter <- 0
  while ((fk_high < 0) && (k_high < max_k)) {
    k_high <- k_high * 2
    fk_high <- f(k_high)
    iter <- iter + 1
    if (iter > 200) stop("Could not bracket the root: try increasing max_k or a different initial guess.")
  }
  if (fk_high < 0) stop("Failed to find k_high that yields operator norm >= target. Increase max_k.")
  root <- uniroot(f, lower = k_low, upper = k_high, tol = tol)$root
  A_scaled <- D + root * B
  return(A_scaled)
}

#single events

K<-200
rho_A<-0.2
nSteps<-600
A_coef <- diag(rho_A, K, K) # initial diagonal matrix
y_0<-rep(0,K)
Sigma_u<-diag(0.5,K,K)
two_upper_indx<- which(col(A_coef)== row(A_coef)+1 | col(A_coef)==row(A_coef)+2,
                       arr.ind = TRUE)
A_coef[two_upper_indx]<-0.000

#compute VAR(1) process
start_time <- Sys.time()
X_t <- matrix(0L, nrow=nSteps, ncol=K)
X_t[1,]<-y_0
#generate the noise
noise <- mvrnorm(n = nSteps, mu = rep(0, K), Sigma = Sigma_u)


for(t in 2:nSteps){
  X_t[t,]<-A_coef%*%X_t[t-1,]
  X_t[t,1:K]<-X_t[t,1:K]+noise[t,]  
  
}

x_t<-X_t[-nrow(X_t),]
y_t<-X_t[-1,]
fit_cv <- cv.glmnet(
  x = as.matrix(x_t),
  y = as.matrix(y_t),
  family = "mgaussian",   
  alpha = 1,
  intercept = FALSE,
  nfolds = 5,
  type.measure = "mse",
  grouped = FALSE
)
beta_star<-coef(fit_cv, s="lambda.min")
beta_star_mtrx<-matrix(0L, nrow=K, ncol=K)
#extract the value and fit in a matrix 
for(i in seq_len(K)){
  beta_vec<-beta_star[paste0("y",i)]@x
  beta_star_mtrx[, i] <- beta_vec[-1] 
}

end_time <- Sys.time()
elapsed  <-difftime(end_time, start_time, units = "mins")
print(paste0("Total elapsed time (mins): ", round(as.numeric(elapsed), 2)))
#mc--------------------------------------------------------------
K<-200
rho_A<-0.2
nSteps<-600
norm_A<-c(0.2,0.92,0.96,1,1.01,1.02,1.03)
gammas<-sapply(norm_A, function(x)(x-0.2)/2) #||A||=alpha+2 gamma
A_coef <- diag(rho_A, K, K)
#n_values<-seq(0,600, 50)#seq of numbes of observations we would like to analyse

## define paramters to generate the VAR(1)
##we assume data is centered, therefore intercepts are a seq of 0
nu<-rep(0,K)
## The var model is gaussian
Sigma_u<-diag(0.5,K,K)
## initial value for the iterations is 0
y_0<-rep(0, K)

#generate a matrix to store lasso errors
#Lasso_error<-matrix(NA, nrow=length(gammas), ncol=length(n_values))
#rownames(Lasso_error)<-paste0("(0.2,", norm_A, ")")
#colnames(Lasso_error)<-as.character(n_values)


##try on a single case to then generalise

A_coef <- diag(rho_A, K, K)

two_upper_indx<- which(col(A_coef)== row(A_coef)+1 | col(A_coef)==row(A_coef)+2,
                       arr.ind = TRUE)
A_coef[two_upper_indx]<-0.415
lasso_error_mc<-rep(NA_real_, 10)
start_time<-Sys.time()
for(j in 1:1){

set.seed(200+j)
X_t <- matrix(0L, nrow=nSteps, ncol=K)
X_t[1,]<-y_0
#generate the noise
noise <- mvrnorm(n = nSteps, mu = rep(0, K), Sigma = Sigma_u)


for(t in 2:nSteps){
  X_t[t,]<-A_coef%*%X_t[t-1,]
  X_t[t,]<-X_t[t,]+noise[t,]  
  
}

x_t<-X_t[-nrow(X_t),]
y_t<-X_t[-1,]
fit_cv <- cv.glmnet(
  x = as.matrix(x_t),
  y = as.matrix(y_t),
  family = "mgaussian",   
  alpha = 1,
  intercept = FALSE,
  nfolds = 5,
  type.measure = "mse",
  grouped = FALSE
)
#compute the beta
beta_star<-coef(fit_cv, s="lambda.min")
beta_star_mtrx<-matrix(0L, nrow=K, ncol=K)
#extract the value and fit in a matrix 
for(i in seq_len(K)){
  beta_vec<-beta_star[paste0("y",i)]@x
  beta_star_mtrx[, i] <- beta_vec[-1] 
  }
lasso_error_mc[j]<-sqrt(sum((beta_star_mtrx-A_coef)^2))

}
end_time <- Sys.time()
elapsed  <-difftime(end_time, start_time, units = "mins")
lasso_error_final<-mean(lasso_error_mc)
print(paste0("Total elapsed time (mins): ", round(as.numeric(elapsed), 2)))

#try with multiple iterations---------------------------------------------
