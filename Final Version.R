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
#matrices--------------------------------------------------------------
K<-200
rho_A<-0.2
norm_A<-c(0.2,0.92,0.96,1,1.01,1.02,1.03)
gammas<-sapply(norm_A, function(x)(x-0.2)/2) #||A||=alpha+2 gamma
A_coef <- diag(rho_A, K, K)
n_values<-seq(0,600, 50)#seq of numbes of observations we would like to analyse

## define paramters to generate the VAR(1)
##we assume data is centered, therefore intercepts are a seq of 0
nu<-rep(0,K)
## The var model is gaussian
sigma_u<-diag(0.5,K,K)
## initial value for the iterations is 0
y_0<-rep(0, K)

#generate a matrix to store lasso errors
Lasso_error<-matrix(NA, nrow=length(gammas), ncol=length(n_values))
rownames(Lasso_error)<-paste0("(0.2,", norm_A, ")")
colnames(Lasso_error)<-as.character(n_values)


##try on a single case to then generalise

A_coef <- diag(rho_A, K, K)
two_upper_indx<- which(col(A_coef)== row(A_coef)+2 | col(A_coef)==row(A_coef)+1,
                       arr.ind = TRUE)
A_coef[two_upper_indx]<-0.360
Y_t<-var_sim(A_coef,nu, Sigma_u=sigma_u, nSteps=600,y0=y_0 )
X_t<-Y_t[-nrow(Y_t),]
y_t<-Y_t[-1,]
fit_lasso<-cv.glmnet(
    x=as.matrix(X_t),
    y=as.matrix(y_t),
    type.measure="mse",
    nfolds=10,
    alpha=1,
    intercept=FALSE
)
  
beta_hat<-coef(fit_lasso, s="lambda.mim")
for(i in seq(0,1000,1)){
  set.seed(2000+i)
  #store values of error for each iterations
  list_error<-c(NA, 1000)
  Y_t<-var_sim(A_coef,1, nu , sigma_u, nSteps=n,y_0 )
  #create lagged value to derive the predictors
  X_t<-as.matrix(Y_t[-nrow(X_t),])
  y_t<-Y_t[-1,]
  #we use built-in function from the package glmnet to compute a lasso 
  #regression
  fit_lasso<-cv.glmnet(
    x=X_t,
    y=y_t,
    type.measure="mse",
    intercept=FALSE,
    alpha=1,
    family="gaussian", 
    
    nfolds=10,
    
    family="gaussian"
  )
  #extract coefficient
  beta_hat<-coef(fit_lasso, s="lambda.mim")
  
  lasso_error[i]<-sqrt(sum(beta_hat))
  
}

for (gamma in seq_along(gammas)) {
  A_coef <- diag(rho_A, K, K)
  two_upper_indx<- which(col(A_coef)== row(A_coef)+2 | col(A_coef)==row(A_coef)+1,
                         arr.ind = TRUE)
  A_coef[two_upper_indx]<-gamma
  
  #we set seed to eállow the process to be replicable
  
  for(i in seq(0,1000,1)){
    set.seed(2000+i)
    #store values of error for each iterations
    list_error<-c(NA, 1000)
    for( n in seq_along(n_values) ){
      
    #we generate the predictors
    Y_t<-var_sim(A_coef,1, nu , sigma_u, nSteps=n,y_0 )
    #create lagged value to derive the predictors
    X_t<-Y_t[-nrow(X_t),]
    y_t<-Y_t[-1,]
    #we use built-in function from the package glmnet to compute a lasso 
    #regression
    fit_lasso<-cv.glmnet(
      x=as.matrix(X_t)
      y=as.matrix(y_t),
      intercept=FALSE,
      alpha=1,
      family="gaussian", 
      grouped=FALSE,
      nfold=10,
      type.measure="mse",
      family="gaussian"
    )
    #extract coefficient
    beta_hat<-coef(fit_lasso, s="lambda.mim")
    
    lasso_error[i]<-sqrt(sum(beta_hat<-))
    
    }
    
  }  
}
#fill the non diagonal values with gamma values
A_coef[two_upper_indx]<-gamma

#check whether the operator norm 



prod_coef<-t(A_coef)%*%A_coef
rho_val<-max(abs(eigen(prod_coef)$values))
rho_val
a<-0.2
scalar<-norm_A/sqrt(rho_val)
A_scalar<-A_coef *scalar
A_scaled_sim<-scaled_A_matrix(A_coef, 1.01)
print(max(abs(eigen(t(A_scalar)%*%A_scalar)$value)))
print(max(abs(eigen(t(A_scaled_sim)%*%A_scaled_sim)$value)))
