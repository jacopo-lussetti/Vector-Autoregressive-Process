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
#single simulation-----------------------------------------------------------
rho_A<-0.2
K<-200
n_Steps<-20
gam<-0.415

## define paramters to generate the VAR(1)
##we assume data is centered, therefore intercepts are a seq of 0
nu<-rep(0,K)
## The var model is gaussian
Sigma_u<-diag(0.5,K,K)
## initial value for the iterations is 0
y_0<-rep(0, K)

#generate the coef matrix

A_coef<-diag(rho_A,K,K)
two_upper_matrix<-which(col(A_coef)==row(A_coef)+1 | col(A_coef)==row(A_coef)+2,
                        arr.ind = TRUE)
A_coef[two_upper_matrix]<-gam

#we generate the VAR(1) process
X_t<-var_sim(A_coef,nu = nu, Sigma_u = Sigma_u,nSteps = n_Steps,y0=y_0)

beta_star<-rep(0,K)

nonzero<-K*0.3
beta_star[sample(1:K,nonzero)]<-rnorm(nonzero, mean=0, sd=0.5)

#we generate the response data from the stochastic regression
noise <-rnorm(n_Steps, mean=0, sd=1)    #noise

y_t<-X_t %*% beta_star  +noise

#lasso regression

lasso_reg<-cv.glmnet(
  x=as.matrix(X_t),
  y=y_t,
  nfolds = 10, 
  family="gaussian",
  intercept=FALSE,
  alpha=1,
  grouped = FALSE
)


beta_lasso_sim <- coef(lasso_reg, s = "lambda.min") 
beta_hat <- as.numeric(beta_lasso_sim[-1])
lasso_iter <- sqrt(sum((beta_hat - beta_star)^2)) 

lasso_iter



#All combinatons----------------------------------
K<-200
rho_A<-0.2
n_values<-seq(30,600,30)
norm_A<-c(0.2,0.92,0.96,1,1.01,1.02,1.03)
iter<- 50#number of iteration for monte carlo
gammas<-sapply(norm_A, function(x)(x-0.2)/2) #||A||=alpha+2 gamma
##we assume data is centered, therefore intercepts are a seq of 0
nu<-rep(0,K)
## The var model is gaussian
Sigma_u<-diag(0.5,K,K)
## initial value for the iterations is 0
y_0<-rep(0, K)
# Initialize matrix to store lasso errors
Lasso_error <- matrix(NA, nrow = length(norm_A), ncol = length(n_values))
colnames(Lasso_error)<-as.character(n_values)
rownames(Lasso_error) <- paste0("(", rho_A, ",", norm_A, ")")
start_time<-Sys.time()
for(gamma_idx in seq_along(gammas)){
  gamma <- gammas[gamma_idx]
  A_coef<-diag(rho_A,K,K)
  two_upper_matrix<-which(col(A_coef)==row(A_coef)+1 | col(A_coef)==row(A_coef)+2,
                              arr.ind = TRUE)
  A_coef[two_upper_matrix]<-gamma
  #compute beta
  
    for(n_idx in seq_along(n_values)){
    n <- n_values[n_idx]
    lasso_error_mc <- rep(NA_real_, iter)
    
    for(i in 1:iter){
        set.seed(400+i)
     
      
      X_t<-var_sim(A_coef,nu = nu, Sigma_u = Sigma_u,nSteps = n,y0=y_0)
      noise_y<- rnorm(n, mean = 0, sd = 1)   #noise
      
      beta_star<-rep(0,K)
      
      nonzero<-round(K*0.3)
      beta_star[sample(1:K,nonzero)]<-rnorm(nonzero, mean=0, sd=0.5)
      y_t<-X_t %*% beta_star  +noise_y
      #lasso regression
      
      lasso_reg<-cv.glmnet(
        x=as.matrix(X_t),
        y=y_t,
        nfolds = 5, 
        family="gaussian",
        intercept=FALSE,
        alpha=1,
        grouped = FALSE
      )
      #compute error
      beta_lasso_sim <- coef(lasso_reg, s = "lambda.min") 
      beta_hat <- as.numeric(beta_lasso_sim[-1,])
      lasso_iter <- sqrt(sum((beta_hat - beta_star)^2)) 
      lasso_error_mc[i]<-lasso_iter
      }#finish mc loop
    #fill a matrix to track average estimator error of lasso
    Lasso_error[gamma_idx,n_idx ] <- mean(lasso_error_mc)    
  }# finish sample size loop
}#finish sample gamma level

end_time <- Sys.time()
elapsed  <-difftime(end_time, start_time, units = "mins")
print(paste("Elapsed time:", elapsed))

#Plotting

# after you built Lasso_error_long and ensured n_obs is numeric

p <- ggplot(
  data = Lasso_error_long,
  mapping = aes(x = n_obs, y = Lasso_Error)
) +
  # map color/linetype to the 'comb' column (the row labels)
  geom_line(aes(color = factor(comb), linetype = factor(comb)), size = 0.5) + 
  scale_color_manual(values = color_palette) +
  scale_linetype_manual(values = line_types) +
  # use actual sample sizes as breaks (not 0..6)
  scale_x_continuous(breaks = seq(0,600, by=100)) +
  scale_y_continuous(breaks = seq(0, 6, by = 1)) +
  labs(
    x = "n",
    y = TeX("$\\|\\hat{\\beta}- \\beta^*\\|_2$"),
    color = TeX("$\\rho(\\alpha), \\|A\\|$")
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    legend.position = c(1, 1),
    legend.justification = c("right", "top"),
    legend.key.size = unit(0.4, 'cm'),
    legend.key.width = unit(0.4, 'cm'),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6)
  ) +
  guides(linetype = "none")  # keep as you had it (shows only color in legend)

p

