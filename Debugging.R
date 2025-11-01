#Debugging
#libraries-------------------------------

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
#Functions--------------------------------
#we define first of all the various functions
##Simulation upper triangular matrix
### To ensure that process is stable, we need to check the abs of eigenvalue
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
#debugging exercise 1---------------------------
#Debugging
n_mc <- 100
K <- 200
n_values <- 600
target_norm <- 1.01

# Fixed parameters for VAR(1) simulation
nu <- rep(0, K)
y_0 <- rep(0, K)
Cov_epsilon <- diag(0.2, K)
beta_star <- rep(0, K)
nonzero <- round(K * 0.2)
beta_star[sample(1:K, nonzero)] <- rnorm(nonzero, mean = 0, sd = 2)

mc_errors <- rep(NA, n_mc)  # Initialize mc_errors outside the loop

for(i in 1:n_mc){
  # Step 1: Generate A_coef with target rho (two upper bands)
  A_coef <- matrix(0, nrow = K, ncol = K)
  two_upper_indx <- which(col(A_coef) == row(A_coef) + 2, arr.ind = TRUE)
  diag(A_coef) <- 0.2
  A_coef[two_upper_indx] <- rnorm(length(two_upper_indx[, 1]), mean = 0, sd = 0.1)
  A_scaled <- scaled_A_matrix(A_coef, target = target_norm)
  
  # Step 2: Simulate VAR(1) data with n = n_values time steps
  pred_sim <- var_sim(A_scaled, nu = nu, Sigma_u = Cov_epsilon, nSteps = n_values, y0 = y_0)
  
  # Step 3: Generate response
  y_t_sim <- as.numeric(pred_sim %*% beta_star) + rnorm(n_values, mean = 0, sd = 3)
  
  # Step 4: Fit LASSO
  lasso_reg_sim <- cv.glmnet(
    x = as.matrix(pred_sim),
    y = y_t_sim,
    type.measure = "deviance",
    alpha = 1,
    nfolds = 10,
    family = "gaussian",
    grouped = FALSE
  )
  
  beta_lasso_sim <- coef(lasso_reg_sim, s = "lambda.1se")
  beta_hat <- beta_lasso_sim[-1]
  mc_errors[i] <- sqrt(sum((beta_hat - beta_star)^2))  # Store error
}

print(mean(mc_errors, na.rm = TRUE))  # Calculate and print the mean error


#debugging for exercise 2---------------------------------------------------

#debug
  alpha<-0.90
  n_value<-1500
  K<-500
  nu <- rep(0, K)
  x_0 <- rep(0, K*2)
  n_mc<-10
  #we compute the covariance matrix for the residuals of VAR(2) process 
  
  
  
  for(i in seq_along(alpha)){
    alpha_i<-alpha[i]
    #we compute the relative covariance_V
    S<-(1+alpha_i)^2/(1-alpha_i^2)^3
    sigma_2<- 1/S
    #now we compute Sigma_V
    
    Sigma_V<-diag(sigma_2, nrow=K, ncol=K)
    
    for(l in seq_along(n_values)){  
      n<-n_values[l]
      errs <- numeric(n_mc)
      nnz<- integer(n_mc)
      bad_runs<-integer(0)
      beta_star <- matrix(0, 2 * K, 1)
      nonzero_first_lag <- round(K * 0.4)
      nonzero_second_lag <- round(K * 0.4)
      beta_star[sample(1:K, nonzero_first_lag)] <- runif(nonzero_first_lag, 2, 3)
      beta_star[sample((K + 1):(2 * K), nonzero_second_lag)] <- runif(nonzero_second_lag, 2, 3)
      
      beta_star <- as.numeric(beta_star)
      
      A1 <- diag(2*alpha_i,K,K)
      A2<-diag(-(alpha_i^2),K,K)
      AA <- cbind(A1, A2)
      #we compute the companion matrix
      AA_comp<-comp_mtrx(AA,2)
      if(max(abs(eigen(AA_comp)$values))>=1){
        warning("The process will be not stable")
      }
      
      # # we compute the Sigma_V based on Yule-Walker equations
      #Sigma_V<-diag(1,nrow(AA_comp),ncol(AA_comp)) - AA_comp %*% t(AA_comp)
      X_t <- matrix(0, nrow = n_value, ncol = 2 * K)  # trajectories matrix nSteps x Kp
      X_t[1, ] <- x_0  # add initial value to initiate the simulation
      
      # Generate noise
      noise <- mvrnorm(n = n, mu = rep(0, K), Sigma = Sigma_V)
      
      # Simulation loop
      for (t in 2:n_value) {
        X_t[t, ] <- AA_comp %*% X_t[t-1, ]
        X_t[t, 1:K] <- X_t[t, 1:K] + nu + noise[t, ]
      }
      
    
      y_t_sim <- as.numeric(X_t %*% beta_star) +rnorm(n_value, mean = 0, sd = 3)
      #we asssume normality of the error term
      # Step 4: Fit LASSO
      
      lasso_reg_sim <- cv.glmnet(
        x = as.matrix(X_t),
        y = as.matrix(y_t_sim),
        intercept=FALSE,
        alpha = 1,
        nfolds = 10, 
        grouped = FALSE
      )
        
        coef_hat<coef(lasso_reg_sim, s="lambda.min")
        coef_hat<-as.numeric(coef_hat[-1])
        errs[i]<-sqrt(sum((coef_hat-beta_star)^2))
        tol_level<- max(abs(coef_hat[!coef_hat==0]))
        nnz[i]<-sum(abs(coef_hat)>tol_level)
    }#loop mc
    valid_errs <- errs[!is.na(errs)]
    Lasso_error_ex_2[k, j] <- mean(valid_errs, na.rm = TRUE)
    cat("alpha=",alpha_i,"observation=",n, "True non-zero", nnz)
    
  }
  }
  print(mean(lasso_iter))
  
  # epsilon<- mvrnorm(n = n_value, mu = rep(0, K*2), Sigma = Sigma_V) #assuming that
  # #now we compute a VAR(2) process in a VAR(1) form
  # 
  #     X_t <- matrix(0, nrow =  n_value, ncol=K*2) #trajectories matrix nSteps x Kp
  #     X_t[1, ] <- y_0 #add initial value to initiate the simulation
  #     
  #     for (t in 2:n_value) {
  #         X_t[t, ] <- AA_comp %*% X_t[t-1, ]
  #         X_t[t, 1:K] <- X_t[t, 1:K] + nu[1:K] +epsilon[t,1:K]
  #     }
  
  
  
    alpha <- c(0.1, 0.3, 0.5,0.6, 0.7, 0.8, 0.9, 0.95)
    n_values <- c(200, 400, 700, 1000, 1300, 1500)
    alpha<-0.1
    n_values<-10
  K <- 500
    n_mc <- 50
    
    nu <- rep(0, K)
    x_0 <- rep(0, K * 2)
    
    # number of true nonzeros in each lag
    nonzero_first_lag  <- round(K * 0.4)
    nonzero_second_lag <- round(K * 0.4)
    #store lasso error
    Lasso_error_ex_2<-matrix(NA,nrow=length(alpha), ncol=length(n_values))
    
    for(a_idx in seq_along(alpha)) {
      alpha_i <- alpha[a_idx]
      
      S <- (1 + alpha_i)^2 / (1 - alpha_i^2)^3
      sigma_2 <- 1 / S
      Sigma_V <- diag(sigma_2, nrow = K, ncol = K)
      
      A1 <- diag(2 * alpha_i, K, K)
      A2 <- diag(-(alpha_i^2), K, K)
      AA <- cbind(A1, A2)
      AA_comp <- comp_mtrx(AA, 2)
      
      if(max(abs(eigen(AA_comp, only.values = TRUE)$values)) >= 1) {
        warning("Alpha ", alpha_i, ": companion matrix spectral radius >= 1 (process may be unstable)")
      }
      
      for(n_idx in seq_along(n_values)) {
        #check how long it takes for each loop 
        start_time <- Sys.time()
        n <- n_values[n_idx]
        
        errs <- numeric(n_mc)
        nnz  <- integer(n_mc)
        tol_level_vec <- numeric(n_mc)
        
        # generate one fixed beta_star for this (alpha,n) combination (you can move outside if you want global fixed)
        set.seed(2025)                          # optional - reproducible beta per (alpha,n)
        beta_star <- numeric(2 * K)
        beta_star[sample(1:K, nonzero_first_lag)] <- runif(nonzero_first_lag, 2, 3)
        beta_star[sample((K + 1):(2 * K), nonzero_second_lag)] <- runif(nonzero_second_lag, 2, 3)
        
        for(mc in 1:n_mc) {
          # simulate X_t for this MC iteration
          X_t <- matrix(0, nrow = n, ncol = 2 * K)
          X_t[1, ] <- x_0
          noise <- mvrnorm(n = n, mu = rep(0, K), Sigma = Sigma_V)
          
          for(t in 2:n) {
            X_t[t, ] <- AA_comp %*% X_t[t - 1, ]
            X_t[t, 1:K] <- X_t[t, 1:K] + nu + noise[t, ]
          }
          
          y_t_sim <- as.numeric(X_t %*% beta_star) + rnorm(n, mean = 0, sd = 3)
          
          lasso_reg_sim <- cv.glmnet(
            x = as.matrix(X_t),
            y = as.numeric(y_t_sim),
            intercept = FALSE,
            alpha = 1,
            nfolds = 10,
            grouped = FALSE
          )
          
          coef_hat <- as.numeric(coef(lasso_reg_sim, s = "lambda.min")[-1])
          errs[mc] <- sqrt(sum((coef_hat - beta_star)^2))
          
  
        } # end MC loop
        end_time <- Sys.time()
        elapsed <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 2)
        valid_errs <- errs[!is.na(errs)]
        Lasso_error_ex_2[ a_idx,n_idx]<-mean(valid_errs, na.rm = TRUE)
        cat("alpha=", alpha_i, " n=", n, " MC valid=", length(valid_errs),
            " median_err=", round(median(valid_errs, na.rm = TRUE), 4),
            " mean_err=", round(mean(valid_errs, na.rm = TRUE), 4),
            " sd=", round(sd(valid_errs, na.rm = TRUE), 4), "\n")
        
        cat("  quantiles (0.1,0.25,0.5,0.75,0.9): ",
            paste(round(quantile(valid_errs, probs = c(.1, .25, .5, .75, .9),
                                 na.rm = TRUE), 3), collapse = ", "), "\n")
        
        cat("Execution time for each (alpha,n) combination", elapsed)
        
        
      } # end n_values loop

    } # end alpha loop
    
#plots
    
    # Convert matrix to data frame and tidy it
    Lasso_error_df2 <- as.data.frame(Lasso_error_ex_2, stringsAsFactors = FALSE)
    Lasso_error_df2$alpha_rho <- rownames( Lasso_error_df2)
    
    Lasso_error_long <- Lasso_error_df %>%
      pivot_longer(
        cols = -c(alpha_rho),  # Use norm_A (not rho) for ||A||
        names_to = "n_obs",
        values_to = "Lasso_Error"
      ) %>%
      mutate(n_obs = as.integer(n_obs))
    
    # Define color and linetype palettes
    color_palette <- c(
      "(0.2,0.2)"   = "black",
      "(0.2,0.92)"  = "red",
      "(0.2,0.96)"  = "blue",
      "(0.2,1)"     = "green",
      "(0.2,1.01)"  = "orange",
      "(0.2,1.02)"  = "purple",
      "(0.2,1.03)"  = "cyan"
    )
    
    line_types <- c(
      "(0.2,0.2)"   = "solid",
      "(0.2,0.92)"  = "dashed",
      "(0.2,0.96)"  = "dashed",
      "(0.2,1)"     = "dotdash",
      "(0.2,1.01)"  = "dotdash",
      "(0.2,1.02)"  = "dotdash",
      "(0.2,1.03)"  = "dotdash"
    )
    
    # Plot
    g <- ggplot(
      data = Lasso_error_long,
      mapping = aes(x = n_obs, y = Lasso_Error)
    ) +
      geom_line(aes(color=factor(alpha_rho), linetype=factor(alpha_rho)), size=0.5) + 
      scale_color_manual(
        values = color_palette
        
      ) +
      scale_linetype_manual(
        values = line_types,  
      ) +
      labs(
        x="n",
        y=TeX("$\\|\\hat{\\beta}- \\beta^*\\|_2$"),
        color = TeX("$\\rho(A),\\,\\|A\\|_2$")  
      )+
      theme_minimal()+
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # make sure axes lines are visible
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.ticks = element_line(colour = "black"),
        
        # add full black border around the plot panel
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        
        # move legend to top-right
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        
        legend.key.size = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        legend.title = element_text(size=8),
        # make legend text a bit smaller
        legend.text = element_text(size = 6)
      )+
      guides(linetype = "none")
    g    