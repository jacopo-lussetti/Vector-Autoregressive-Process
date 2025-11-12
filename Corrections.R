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
library(forecast)
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
set.seed(1245)
rho_A<-0.2
K<-200
n_Steps<-600
gam<-0.400

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
##compute first the lambda as function of dimension 
##we definet he uncodintional variance sigma and tune parameter C
sigma<-1
## recall the definition: C*sigma_epsilon *||Sigma_x ||^(1/2)* sqrt(log(K)/n)
#derive forst sogma
#X is a random design which is zero mean Gaussian
sum_sqrt<-colSums(abs(X_t)^2)/n_Steps # column normalised
sum_sqrt<-sum_sqrt/(sqrt(n_Steps)) # column normalised
Sigma_x<-max(sum_sqrt)
C<-seq(0,10,0.2)
lambda_n<-sapply(C, function(c)(2*(c*sigma*Sigma_x*sqrt(log(K)/n_Steps))))

#we compute the best C by out-of-sample
train_n<-round(n_Steps*0.7)
X_train<-X_t[1:train_n,]
y_train<-y_t[1:train_n]
X_test<-X_t[(train_n+1):n_Steps,]
y_test<-y_t[(train_n+1):n_Steps]

## now we compute 
lasso_reg<-glmnet(
  x=as.matrix(X_train),
  y=as.numeric(y_train),
  family="gaussian",
  intercept=FALSE,
  alpha=1,
  lambda=lambda_n
)
#predict on test set
df<-lasso_reg$df #degrees of freedom
Y_hat<-predict(lasso_reg, newx=as.matrix(X_test))
Y_hat<-as.matrix(Y_hat)
RSS <- colSums((Y_hat - y_test)^2)
sigma_hat<-RSS /(n_Steps-df)
#compute BIC values
BIC<-n_Steps*log(2*pi*sigma_hat) + RSS/sigma_hat + (log(n_Steps))*df
best_C_idx<-which.min(BIC)
best_lambda<-lambda_n[best_C_idx]
#now we compute again the lasso with the best lambda
lasso_reg<-glmnet(
  x=as.matrix(X_t),
  y=as.numeric(y_t),
  family="gaussian",
  intercept=FALSE,
  alpha=1,
  lambda=best_lambda
)

beta_lasso_sim <- coef(lasso_reg) 
beta_hat <- as.numeric(beta_lasso_sim[-1])
lasso_iter <- sqrt(sum((beta_hat - beta_star)^2)) 

lasso_iter




#All combinatons----------------------------------
  
  K <- 200
  rho_A <- 0.2
  n_values <- seq(30, 600, 30)
  norm_A <- c(0.2, 0.92, 0.96, 1, 1.01, 1.02, 1.03)
  iter <- 100 # Monte Carlo iterations
  gammas <- sapply(norm_A, function(x) (x - 0.2) / 2) # ||A|| = alpha + 2*gamma
  nu <- rep(0, K)
  Sigma_u <- diag(0.5, K, K)
  y_0 <- rep(0, K)
  # Storage for coefficient estimates
  Coef_estimate <- vector("list", length(gammas))
  names(Coef_estimate) <- paste0("gamma_", seq_along(gammas))
  for (g in seq_along(Coef_estimate)) {
    Coef_estimate[[g]] <- vector("list", length(n_values))
    names(Coef_estimate[[g]]) <- paste0("n_", n_values)
  }
  # Storage error
  Lasso_error <- list(
    
    Lasso_Error = matrix(NA, nrow = length(norm_A), ncol = length(n_values)),
    Lambda = matrix(NA, nrow = length(norm_A), ncol = length(n_values)),
    Estimates = Coef_estimate
    
  )
  
  colnames(Lasso_error$Lasso_Error) <- as.character(n_values)
  rownames(Lasso_error$Lasso_Error) <- paste0("(", rho_A, ",", norm_A, ")")
  colnames(Lasso_error$Lambda) <- as.character(n_values)
  rownames(Lasso_error$Lambda) <- paste0("(", rho_A, ",", norm_A, ")")
  #storage lambdas
  
  
  start_time <- Sys.time()
  
  for (gamma_idx in seq_along(gammas)) {
    gamma <- gammas[gamma_idx]
    A_coef <- diag(rho_A, K, K)
    two_upper_matrix <- which(col(A_coef) == row(A_coef) + 1 | col(A_coef) == row(A_coef) + 2,
                              arr.ind = TRUE)
    A_coef[two_upper_matrix] <- gamma
    
    for (n_idx in seq_along(n_values)) {
      n <- n_values[n_idx]
      lasso_error_mc <- rep(NA_real_, iter)
      #store the lambas as function of dimension
      lambda_i<-rep(NA_real_, iter)
      beta_hat_mat <- matrix(NA_real_, nrow = K, ncol = iter)
      beta_star_mat <- matrix(NA_real_, nrow = K, ncol = iter)
      
      for (i in 1:iter) {
        set.seed(400 + i)
        
        # simulate VAR and response
        X_t <- var_sim(A_coef, nu = nu, Sigma_u = Sigma_u, nSteps = n, y0 = y_0)
        noise_y <- rnorm(n, mean = 0, sd = 1)   # regression noise; sigma = 1
        sigma <- 1
        
        beta_star <- rep(0, K)
        nonzero <- round(K * 0.3)
        beta_star[sample(1:K, nonzero)] <- rnorm(nonzero, mean = 0, sd = 2)
        beta_star_mat[, i] <- beta_star
        
        y_t <- X_t %*% beta_star + noise_y
        
        # split train/test (time-aware)
        train_n <- round(n * 0.7)
        X_train <- X_t[1:train_n, ]
        y_train <- y_t[1:train_n]
        X_test <- X_t[(train_n + 1):n, ]
        y_test <- y_t[(train_n + 1):n]
        
        ## recall the definition: C*sigma_epsilon *||Sigma_x ||^(1/2)* sqrt(log(K)/n)
        #derive forst sigma
        #X is a random design which is zero mean Gaussian
        sum_sqrt<-colSums(abs(X_train)^2)/train_n 
        sum_sqrt<-sum_sqrt/(sqrt(train_n))# column normalised
        
        #C has to be greater or equal to 
        # grid for C (avoid 0)
        Sigma_x<-max(sum_sqrt)
        C<-seq(Sigma_x,Sigma_x*50,0.2)
        #compute the function of dimensions for lambads
  
        
        
        lambda_n <- sapply(C, function(c) (2*(c * sigma * Sigma_x *
                                                sqrt(log(K) / train_n))))
        
        # fit on training for all lambdas
        lasso_reg_C <- glmnet(
          x = as.matrix(X_train),
          y = as.numeric(y_train),
          family = "gaussian",
          intercept = FALSE,
          alpha = 1,
          standardize = FALSE,
          lambda = lambda_n
        )
        
        # training rss and df for each lambda (aligned)
        coef_mat <- as.matrix(coef(lasso_reg_C))  # (K+1) x L
        L <- ncol(coef_mat)
        rss_train <- numeric(L)
        df_vec <- numeric(L)
        for (j in seq_len(L)) {
          beta_j <- coef_mat[-1, j]  # drop intercept
          preds_train <- as.numeric(X_train %*% beta_j)
          resid_train <- y_train - preds_train
          rss_train[j] <- sum(resid_train^2)
          df_vec[j] <- sum(beta_j != 0)
        }
        
        # compute BIC on training data (classical BIC* simplified)
        # BIC* = n_train * log(RSS / n_train) + log(n_train) * df
        n_tr <- train_n
        BIC_vals <- n_tr * log(rss_train / n_tr) + log(n_tr) * df_vec
        
        best_C_idx <- which.min(BIC_vals)
        best_lambda <- lambda_n[best_C_idx]
        
        # refit on full data with the selected lambda
        lasso_reg_full <- glmnet(
          x = as.matrix(X_t),
          y = as.numeric(y_t),
          family = "gaussian",
          intercept = FALSE,
          alpha = 1,
          standardize = FALSE,
          lambda = best_lambda
        )
        
        # extract coefficients and compute error
        beta_lasso_sim <- as.numeric(coef(lasso_reg_full)[-1])  # vector length K
        beta_hat <- beta_lasso_sim
        beta_hat_mat[, i] <- beta_hat
        
        lasso_iter <- sqrt(sum((beta_hat - beta_star)^2))
        lasso_error_mc[i] <- lasso_iter
      } # end MC loop
      
      Coef_estimate[[gamma_idx]][[n_idx]] <- list(
        beta_hat_mat = beta_hat_mat,
        beta_star_mat = beta_star_mat,
        lasso_error_iter = lasso_error_mc
      )
      
      Lasso_error$Lasso_Error[gamma_idx, n_idx] <- mean(lasso_error_mc)
      Lasso_error$Lambda[gamma_idx, n_idx] <- best_lambda
    } # finish n loop
  } # finish gamma loop
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")
  print(paste("Elapsed time:", elapsed))
  
  #Plotting
  
  
  # --- assume Lasso_error exists and rownames are the "(rho,norm)" strings
  Lasso_error_df <- as.data.frame(Lasso_error$Lasso_Error, stringsAsFactors = FALSE)
  Lasso_error_df$alpha <- rownames(Lasso_error_df)
  
  Lasso_error_long <- Lasso_error_df %>%
    pivot_longer(
      cols = -alpha,
      names_to = "n_obs",
      values_to = "Lasso_Error"
    ) %>%
    mutate(n_obs = as.integer(n_obs),
           alpha = as.factor(alpha))   # ensure factor
  
  # make sure your palette keys exactly match the levels of alpha
  color_palette <- c(
    "(0.2,0.2)"   = "black",
    "(0.2,0.92)"  = "red",
    "(0.2,0.96)"  = "grey40",
    "(0.2,1)"     = "purple",
    "(0.2,1.01)"  = "cyan3",
    "(0.2,1.02)"  = "gold2",
    "(0.2,1.03)"  = "blue"
  )
  
  line_types <- c(
    "(0.2,0.2)"   = "solid",
    "(0.2,0.92)"  = "dashed",
    "(0.2,0.96)"  = "dotdash",
    "(0.2,1)"     = "dotted",
    "(0.2,1.01)"  = "longdash",
    "(0.2,1.02)"  = "dotdash",
    "(0.2,1.03)"  = "solid"
  )
  
  # Defensive: ensure palette keys match the actual alpha levels
  pal_names <- names(color_palette)
  if (!all(levels(Lasso_error_long$alpha) %in% pal_names)) {
    stop("Some alpha levels are missing in color_palette / line_types. Check names.")
  }
  
  # compute a sensible y range / breaks
  ymin <- min(Lasso_error_long$Lasso_Error, na.rm = TRUE)
  ymax <- max(Lasso_error_long$Lasso_Error, na.rm = TRUE)
  y_breaks <- pretty(c(ymin, ymax), n = 6)
  
  p <- ggplot(
    data = Lasso_error_long,
    mapping = aes(x = n_obs, y = Lasso_Error, group = alpha)
  ) +
    geom_line(aes(color = alpha, linetype = alpha), size = 0.8) +
    scale_color_manual(values = color_palette, breaks = names(color_palette)) +
    scale_linetype_manual(values = line_types, breaks = names(line_types)) +
    scale_x_continuous(breaks = seq(0, max(Lasso_error_long$n_obs, na.rm = TRUE),
                                    by = 100)) +
    expand_limits(x=700)+
    scale_y_continuous(breaks = y_breaks) +
    labs( x = "n", y = TeX("$\\|\\hat{\\beta}- \\beta^*\\|_2$"), 
          color = TeX("$\\rho(\\alpha), \\|A\\|$") )+
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black", size = 0.5),
      axis.ticks = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.6),
      legend.position = c(1, 1),
      legend.justification = c("right", "top"),
      legend.key.size = unit(0.45, "cm"),
      legend.key.width = unit(0.6, "cm"),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 5)
    ) +
    guides(linetype = "none")   # hide linetype from legend, keep color
  
  p

#Simulation Chapter 6 --------------------------------------------
#For this simulation, we keep the dependence structure unchanged, by fixing 
#dependence parameter
#single value simulation
p<-1024
n<-100
X_t<-matrix(NA_real_, n+3, p)
X_t[1:2,]<-0

#generate white noise
Sigma_nu<-diag(1,p)    #Sigma_nu
wn<-mvrnorm(n=n+3, mu=rep(0, p),Sigma = Sigma_nu)
eps_t<-matrix(NA_real_, nrow=n+3, ncol=p) #matrix to store the MA component
#generate MA(2) process
for(t in 3:nrow(wn)){
  
  eps_t[t,]<-wn[t,] - 0.8 * wn[(t-1),] + 0.16 * wn[(t-2),]
  X_t[t,]<-1.2*X_t[(t-1),]- 0.36 * X_t[(t-2),] + eps_t[t,]
}
#we remove the first two observations to avoid the effect of initial values
eps_t<-eps_t[-c(1,2,3),]
dim(eps_t)
X_t<-X_t[-c(1,2,3),]
dim(X_t)

#lambda as function of dimension
lambda<-sqrt(log(p)/n)
#generate true beta

beta_star<-rep(0, p)
#determine the the sparsity
non_0<-round(sqrt(p))
non_zero_indices <- sample(1:p, non_0)



# Set the non-zero coefficients
beta_star[non_zero_indices] <- rnorm(non_0, mean = 0, sd = 0.5)

#ensure that SNR is 1.2
SNR_target <- 1.2

# --- compute "noise variance" to be used for the response y ---
# Here we use the average per-scalar variance of eps_t_trim as the noise variance sigma_e^2
sigma_noise2 <- mean(eps_t^2)   # scalar

# --- scale beta so ||beta||_2^2 = SNR_target * sigma_noise2 (because Gamma_X(0)=I) 
beta_norm2_current <- sum(beta_star^2)

beta_norm2_target <- SNR_target * sigma_noise2
c_scale <- sqrt(beta_norm2_target / beta_norm2_current)
beta_star <- beta_star * c_scale
#generate  response
y_t<-X_t %*% beta_star + rnorm(n, mean=0, sd=sqrt(sigma_noise2))
#now we compute the lasso regression
lasso_reg<-glmnet(
  x=as.matrix(X_t),
  y=as.numeric(y_t),
  family="gaussian",
  intercept=FALSE,
  alpha=1,
  lambda=lambda
)
beta_lasso_sim <- coef(lasso_reg)
beta_hat <- as.numeric(beta_lasso_sim[-1])
lasso_error <- sqrt(sum((beta_hat - beta_star)^2))
lasso_error

# All combinations-----------------------------------------------------
p_values <- c(128, 264, 512, 1024)  # Number of variables
n_values <- c(250, 350, 450, 550, 650, 850, 1150, 1350, 1750, 2000)  # Number of observations

iter <- 50  # Number of iterations for the MC simulation

# Storage for lasso errors
lasso_error <- matrix(NA_real_, nrow = length(p_values), ncol = length(n_values))
colnames(lasso_error) <- as.character(n_values)
rownames(lasso_error) <- as.character(p_values)

start_time <- Sys.time()

for (p_index in seq_along(p_values)) {
  p <- p_values[p_index]
  cat("Starting simulations for p =", p, "\n")
  
  Sigma_nu <- diag(1, p)
  
  for (j in seq_along(n_values)) {
    n <- n_values[j]
    cat("  n =", n, "\n");
    
    # Define the matrices
    X_t <- matrix(NA_real_, nrow = n + 3, ncol = p)
    X_t[1:2, ] <- 0  # Initial values
    
    # Create matrix for residuals
    eps_t <- matrix(NA_real_, nrow = n + 3, ncol = p)
    
    results_n <- numeric(iter)  # Initialize results_n as a numeric vector
    
    for (i in 1:iter) {
      set.seed(100 + i)
      
      # White noise component
      wn <-   mvrnorm(n = n + 3, mu = rep(0, p), Sigma = Sigma_nu)
      
      
      # Compute the ARMA model
      for (tt in 3:(n + 3)) {
        
        eps_t[tt, ] <- wn[tt, ] - 0.8 * wn[(tt - 1), ] + 0.16 * wn[(tt - 2), ]
        X_t[tt, ] <- 1.2 * X_t[(tt - 1), ] - 0.36 * X_t[(tt - 2), ] + eps_t[tt, ]
      }
      
      # Remove the first two observations to avoid the effect of initial values
      eps_obs <- eps_t[-c(1, 2, 3), , drop = FALSE]
      X_obs <- X_t[-c(1, 2, 3), , drop = FALSE]
      
      # Compute the sparse true parameter matrix
      beta_star <- rep(0, p)
      non_0 <- round(sqrt(p))
      non_zero_indices <- sample(1:p, non_0)
      
      # Set the non-zero coefficients
      beta_star[non_zero_indices] <- rnorm(non_0, mean = 0, sd = 0.5)
      
      # Rescale beta to ensure that SNR is 1.2
      SNR_target <- 1.2
      sigma_noise2 <- mean(eps_obs^2)  # scalar
      
      beta_norm2_current <- sum(beta_star^2)

      
      beta_norm2_target <- SNR_target * sigma_noise2
      c_scale <- sqrt(beta_norm2_target / beta_norm2_current)
      beta_star <- beta_star * c_scale
      
      # Generate response
      y_t <- as.numeric(X_obs %*% beta_star + rnorm(n, mean = 0, sd = sqrt(sigma_noise2)))
      
      lambda_np <- sqrt(log(p) / n)
      lasso_reg <- 
        glmnet(
          x = as.matrix(X_obs),
          y = as.numeric(y_t),
          family = "gaussian",
          intercept = FALSE,
          alpha = 1,
          lambda = lambda_np
        )
  
      
      # Extract estimates
      cf <- as.vector(coef(lasso_reg, s = lambda_np))
      beta_hat <- as.numeric(cf[-1])
      results_n[i] <- sqrt(sum((beta_hat - beta_star)^2))
    }
    
    avg_coef <- mean(results_n, na.rm = TRUE)
    lasso_error[p_index, j] <- avg_coef
  }
}
save(lasso_error, file="6.1.RData")
Lasso_error_df <- as.data.frame(lasso_error, stringsAsFactors = FALSE)

# if rownames are numeric p-values, make them explicit as alpha label
Lasso_error_df$alpha <- rownames(Lasso_error_df)

# pivot to long format
Lasso_error_long <- Lasso_error_df %>%
  pivot_longer(
    cols = -alpha,
    names_to = "n_obs",
    values_to = "Lasso_Error"
  ) %>%
  mutate(
    n_obs = as.integer(n_obs),
    alpha = factor(alpha, levels = unique(alpha))   # ensure preserved order
  )

# ---- build color & linetype palettes automatically to match alpha levels
alpha_levels <- levels(Lasso_error_long$alpha)
k <- length(alpha_levels)

color_palette <- c(
  "128"   = "black",
  "256"     = "purple",
  "512"  = "gold2",
  "1024"  = "blue"
)

line_types <- c(
  "128"   = "solid",
  "256"  = "dotdash",
  "512"  = "dotdash",
  "1024"  = "dotdash",
)

ymin <- min(Lasso_error_long$Lasso_Error, na.rm = TRUE)
ymax <- max(Lasso_error_long$Lasso_Error, na.rm = TRUE)
y_breaks <- pretty(c(ymin, ymax), n = 6)

# ---- plot
q <- ggplot(
  data = Lasso_error_long,
  mapping = aes(x = n_obs, y = Lasso_Error, group = alpha)
) +
  geom_line(aes(color = alpha, linetype = alpha), size = 0.9) +
  scale_color_manual(values = color_palette, breaks = names(color_palette)) +
  scale_linetype_manual(values = line_types, breaks = names(line_types)) +
  geom_point(aes(color = alpha), size = 1.2) +
  scale_color_manual(values = palette_colors, name = TeX("$p$ (variables)")) +
  scale_linetype_manual(values = linetype_vec, name = NULL) +
  scale_x_continuous(
    breaks = seq(min(Lasso_error_long$n_obs, na.rm = TRUE),
                 max(Lasso_error_long$n_obs, na.rm = TRUE),
                 by = 100),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(breaks = y_breaks) +
  labs(
    x = "n (number of observations)",
    y = TeX("$\\|\\hat{\\beta} - \\beta^*\\|_2$"),
    color = TeX("$p$ (variables)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = 0.5),
    axis.ticks = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.6),
    legend.position = c(1, 1),
    legend.justification = c("right", "top"),
    legend.key.size = unit(0.5, "cm"),
    legend.key.width = unit(0.7, "cm"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  ) +
  guides(color = guide_legend(override.aes = list(linetype = linetype_vec, size = 1.1)),
         linetype = "none")  # hide duplicate linetype legend entry

# Show plot
q

#Exercise 6.2--------------------------------------------
##define functions
### autocovariance matrix from equation 2.1.32 Lukpol
#(I_K^2- A_1 \otimes A_1)^-1 vec(Sigma_u)  
Gamma_y_0<-function(A1, K, Sigma_u){
  I_k2<-diag(1, K^2)
  A_o_A<-kronecker(A1, A1)
  vec_Gamma_0<-solve(I_k2 - A_o_A)%*%as.vector(Sigma_u)
  #now convert back to matrix form
  var_mtrx<-matrix(vec_Gamma_0, nrow=K, ncol=K)
}
## Lasso regression

## Small VAR  

p_a<-10
d_a<-1
T_a<-c(30,50)
rho<-c(0.5,0.7,0.9)

#### Define the transaction Matrix 
set.seed(345)
non_0<-round(0.05 * p_a^2)
element_ii<-runif(p_a, -1,1)
A1_ii<-diag(element_ii)

#### Define the non-diagonal elements
A1_od<-matrix(0, nrow=p_a, ncol=p_a)
off_diag<-which(row(A1_od)!= col(A1_od))
random<-sample(off_diag, size=non_0, replace = FALSE)

A1_od[random] <- runif(non_0, min = -0.1, max = 0.1)
A1_od
A1<-A1_ii + A1_od
View(A1)
#ensure that the process this stable by normalising A1
prod<-t(A1)%*%A1
norm_A1<-max(eigen(prod)$values)
A1<-(A1 / (norm_A1))*0.90
max(eigen(A1)$values) ## 0.9
### Sigmas - we try for a fixed dependence matrix

Sigma_I<-diag(1, p_a)
Sigma_I[which(col(Sigma_I)!= row(Sigma_I) & col(Sigma_I)<= round(p_a/2))]<-0.5
#
Sigma_II<-diag(1, p_a)
Sigma_II[which((Sigma_II!= row(Sigma_I) & col(Sigma_I)<= round(p_a/2)) 
         | round(p_a/2)<= row(Sigma_II))]<-0.5
Sigma_II
Sigma_III
####ensure now that SNR is fixed
varX_0<-Gamma_y_0(A1, p_a, diag(1, p_a))



#### Off-diagonal 

#### We use sample() function to select randomly some off diagonal values



#### select randomly off diago

### Medium VAR
p_b<-50
d_b<-1
T_b<-c(80,120,160)

#Define the function for




## OLS Estimation

## Ridge Regression

## Comparison