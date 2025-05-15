# ------------------------------------------------------------------------------
# ARE 242 Problem Set 2
# Spring 2025
# Author: Yazen Kashlan
# Date: 5/8/2025
# 
# Outline: This script 
# - Sets ups parameters: beta, alpha, omega, price, sigma 
# - creates dataframes: X, M, V, e 
# - creates functions which take in params, data, (and preference shocks):
# 
# - Question 1: estimate logit model using aggregate shares and assuming no heterogeneity
# - Question 2: estimate Mixed Logit RC Model using aggregate data

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Setup
#install.packages("foreach")
#install.packages("evd")
#install.packages("doParallel")
#install.packages("latex2exp")   # if not already installed
#library(EmpiricalIO)
#library(magrittr)
#library(foreach)
#library(ggplot2)
library(latex2exp)
#library(evd)
library(doParallel)
library(dplyr)

registerDoParallel()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Set seed and parameters
rm(list=ls())

folder <- "/Users/yazenkashlan/Documents/GitHub/personal-code/ARE242/"
folder_output <- "/Users/yazenkashlan/Documents/GitHub/personal-code/ARE242/output/"

# seed
set.seed(1234)
# number of products
J <- 10
# dimension of product characteristics including the intercept
K <- 3
# number of markets
T <- 100
# number of consumers per market
N <- 500
# number of Monte Carlo draws (more on this below)
L <- 500

# set parameters of interests: 
#price ones
a <- 0.5
omega <- 1

#Note : the constant beta01 is set to 4 and other two beta0s are drawn; 
beta <- rnorm(K); # 3 betas
beta[1] <- 4      
beta
sigma <- abs(rnorm(K)); # 3 sigmas
sigma
#displayed as output after this cell is executed are betas first and then sigma of betas


#The value of the auxiliary parameters are set as follows:
price_xi <- 1
prop_jt <- 0.6 # share of market to drop in random draws
sd_x <- 0.5
sd_c <- 0.05
sd_p <- 0.05

# make product characteristics data
X <- matrix(sd_x * rnorm(J * (K - 1)), nrow = J) # first two characteristics
X <- cbind(rep(1, J), X)  # third characteristic
colnames(X) <- paste("x", 1:K, sep = "_") # label characteristics
X <- data.frame(j = 1:J, X) %>% # convert array to dataframe with index j
  tibble::as_tibble()
X <- rbind(rep(0, dim(X)[2]), X) # add outside option as first row
head(X)

# make market-product data
M <- expand.grid(j = 1:J, t = 1:T) %>% # jt grid 10 x 100 = 1000 rows
  tibble::as_tibble() %>%
  dplyr::mutate(
    xi = 0, # endogenous component
    #note we will change this later allowing for price to be related to xi
    c = exp(sd_c * rnorm(J*T)),
    p = exp(price_xi * xi + sd_p * rnorm(J*T)) + c) 
# Take 60% sample of M
M <- M %>%
  dplyr::group_by(t) %>%
  dplyr::sample_frac(prop_jt) %>%
  dplyr::ungroup()
# add outside option in every market (100 outside options)
outside <- data.frame(j = 0, t = 1:T, xi = 0, c = 0, p = 0)
M <- rbind( M, outside) %>%
  dplyr::arrange(t, j)
head(M)

# Generate the consumer-level heterogeneity (v and nu)
# make consumer-market data
V      <- matrix(rnorm(N * T * (K + 1)), nrow = N * T) # NT grid 500 x 100 = 50,000 rows
colnames(V) <- c(paste("v_x", 1:K, sep = "_"), "v_p")
V <- data.frame(expand.grid(i = 1:N, t = 1:T), V) %>%
  tibble::as_tibble()
head(V)

# Combine all datasets
# To make choice data, let's create df
df <- expand.grid(t = 1:T, i = 1:N, j = 0:J) %>% # tij = 50,000 * 7 
  tibble::as_tibble() %>%
  dplyr::left_join(V, by = c("i", "t")) %>%
  dplyr::left_join(X, by = c("j")) %>%
  dplyr::left_join(M, by = c("j", "t")) %>%
  dplyr::filter(!is.na(p)) %>%
  dplyr::arrange(t, i, j)
head(df)

# Draw preferene shocks (epsilon_ijt), extreme value type 1
e <- evd::rgev(dim(df)[1])



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# DEFINE FUNCTIONS

# (1) Define indirect utility function
compute_indirect_utility <-  function( df, beta,  sigma, a, omega ) {
  #Recall that even though here xi=0 and no random coeff, 
  #This function is coded to allow for that with xi and v_x and v_p 
  #sigma and omega
  
  # extract matrices
  X <- as.matrix(dplyr::select(df, dplyr::starts_with("x_")))
  p <- as.matrix(dplyr::select(df, p)) 
  v_x <- as.matrix(dplyr::select(df, dplyr::starts_with("v_x")))
  v_p <- as.matrix(dplyr::select(df, v_p))
  xi <- as.matrix(dplyr::select(df, xi))
  
  # random coefficients
  beta_i <- as.matrix(rep(1, dim(v_x)[1])) %*% t(as.matrix(beta)) + v_x %*% diag(sigma) 
  alpha_i <- - exp(a + omega * v_p)
  
  # conditional mean indirect utility
  value <- as.matrix(rowSums(beta_i * X) + p * alpha_i + xi) 
  colnames(value) <- "u"
  
  return(value)  
}

# (2) Define the function compute choice to compute u and q
compute_choice <-function( X, M, V,  e, beta,  sigma,a, omega) {
  #This program takes df and e and returns df with u and q
  
  # constants
  T <- max(M$t)
  N <- max(V$i)
  J <- max(X$j)
  
  # make choice data
  df <- expand.grid(t = 1:T,   i = 1:N,   j = 0:J) %>%
    tibble::as_tibble() %>%
    dplyr::left_join( V, by = c("i", "t")) %>%
    dplyr::left_join(X, by = c("j") ) %>%
    dplyr::left_join( M,  by = c("j", "t")) %>%
    dplyr::filter(!is.na(p)) %>%
    dplyr::arrange( t,  i,   j   )
  
  # compute indirect utility
  u <-  compute_indirect_utility(df, beta,  sigma, a, omega)
  
  # add u and e
  df_choice <- data.frame( df,  u, e) %>%
    tibble::as_tibble()
  
  # make choice
  df_choice <- df_choice %>%
    dplyr::group_by( t,    i  ) %>%
    dplyr::mutate(q = ifelse(u + e == max(u + e), 1, 0)) %>%
    dplyr::ungroup()
  
  # return
  return(df_choice)
}

#compute choice 
df_choice <- compute_choice(X, M, V, e, beta, sigma, a, omega)
head(df_choice)


# (3) Define the function compute_share to compute s and share inversion y
compute_share <-function(X,  M,  V, e, beta, sigma,  a, omega) {
  #This program takes df and e (ijt-level), computes u and q, 
  # then returns market share s and BLP inversion y (jt-level)
  
  # constants
  T <- max(M$t)
  N <- max(V$i)
  J <- max(X$j)
  
  # compute choice
  df_choice <- 
    compute_choice(X,  M, V, e, beta, sigma, a, omega)
  
  # make share data
  df_share <- 
    df_choice %>%
    dplyr::select(-dplyr::starts_with("v_"), -u, -e, -i) %>%
    dplyr::group_by(t, j ) %>%
    dplyr::mutate(q = sum(q)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(t, j, .keep_all = TRUE) %>%
    dplyr::group_by(t) %>%
    dplyr::mutate(s = q/sum(q)) %>%
    dplyr::ungroup()
  
  # log share difference
  df_share <- 
    df_share %>%
    dplyr::group_by(t) %>%
    dplyr::mutate(y = log(s/sum(s * (j == 0)))) %>%
    dplyr::ungroup()
  return(df_share)
}

# we will work with aggregate data df_share (jt-level)
df_share <-compute_share(X, M, V, e, beta, sigma, a, omega)
head(df_share)
summary(df_share)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Assignment questions

# ------------------------------------------------------------------------------
# Question 1
# Estimate Logit model using agg. share data - no heterogeneity (4 params)
result_logit <- lm(data = df_share, formula = y ~ - 1 + x_1 + x_2 + x_3 + p)
summary(result_logit)


# ------------------------------------------------------------------------------
# Question 2 - Mixed Logit RC Model - aggregate data
# Estimate logit model using agg. share data - with heterogeneity (8 params)
# Use V_mcmc and e_mcmc instead of V and e above

## draw mcmc V (This was NT before, now it's LT)
V_mcmc <- matrix(rnorm(L*T*(K + 1)), nrow = L*T) # LT grid 100 x 500 = 50,000 rows
colnames(V_mcmc) <- c(paste("v_x", 1:K, sep = "_"), "v_p")
V_mcmc <- data.frame( expand.grid(i = 1:L, t = 1:T),V_mcmc) %>%
  tibble::as_tibble()
V_mcmc

## draw mcmc e
df_mcmc <- expand.grid(t = 1:T, i = 1:L, j = 0:J) %>%
  tibble::as_tibble() %>%
  dplyr::left_join(V_mcmc, by = c("i", "t")) %>%
  dplyr::left_join(X, by = c("j")) %>%
  dplyr::left_join(M, by = c("j", "t")) %>%
  dplyr::filter(!is.na(p)) %>%
  dplyr::arrange(t, i, j)

# draw idiosyncratic shocks
e_mcmc <- evd::rgev(dim(df_mcmc)[1])
head(e_mcmc)

# compute predicted share
df_share_mcmc <-compute_share(X, M, V_mcmc, e_mcmc, beta, sigma, a, omega)

# Vectorize the parameters to a vector `theta` because `optim` requires the maximiand to be a vector.
theta <- c(beta, sigma, a, omega)


#define the objective function to optimize
compute_nlls_objective_pset2 <- function(theta, df_share, X, M, V_mcmc, e_mcmc) {
  # constants
  K <- length(grep("x_", colnames(X)))
  # extract parameters
  beta <- theta[1:K]
  sigma <- theta[(K + 1):(2 * K)]
  a <- theta[2 * K + 1]
  omega <- theta[2 * K + 2]
  # compute predicted share
  df_share_mcmc <-compute_share(X,  M,  V_mcmc, e_mcmc, beta, sigma, a, omega)
  # compute distance
  distance <- mean((df_share_mcmc$s - df_share$s)^2)
  # return
  return(distance)   
}


# ------------------------------------------------------------------------------
# Hold all but one of the parameters fixed at true value and optimize the one

# 1. Draw a graph of the objective function that varies each parameter 
# from 0.5, 0.6, ..., 1.5 of the true value, 
# while keeping all other parameters at the true values.

#The graphs with the true shocks:
label <- c(paste("\\beta_", 1:K, sep = ""), paste("\\sigma_", 1:K, sep = ""), "a", "\\omega")
label <- paste("$", label, "$", sep = "")

# get the graphs into plots display
# Loop over each parameter in theta, theta[i]
graph_true <- foreach (i = 1:length(theta)) %do% {
  # assign theta in loop
  theta_i <- theta[i]
  # create seq from 0.5-1.5 *theta
  theta_i_list <- theta_i * seq(0.5, 1.5, by = 0.1)
  
  objective_i <- foreach (theta_ij = theta_i_list, .combine = "rbind") %dopar% {
    # collect original theta again
    theta_j <- theta
    # replace theta of interest with the item from the sequence in this loop
    theta_j[i] <- theta_ij
    # estimate the squared distance from NLLS with the new parameters
    objective_ij <- compute_nlls_objective_pset2( theta_j, df_share, X, M, V, e)
    return(objective_ij)
  }
  # collect the squared distance for each parameter value
  df_graph <- data.frame(x = theta_i_list, y = objective_i)
  # plot that optimized value against the sequence around the true value
  g <- ggplot(data = df_graph, aes(x = x, y = y)) +
    geom_point() +
    geom_vline(xintercept = theta_i, linetype = "dotted") +
    ylab("objective function") + xlab(TeX(label[i]))
  
  # save images
  graph_label = i
  ggsave(glue("{folder_output}PSet2_q2_1_theta{graph_label}.png"), g)
  return(g)
}

# plot the graphs
saveRDS(graph_true, file = "/Users/yazenkashlan/Documents/GitHub/personal-code/ARE242/Pset2_graph_true.rds")
graph_true <- readRDS(file = "/Users/yazenkashlan/Documents/GitHub/personal-code/ARE242/Pset2_graph_true.rds")
graph_true # plot

# ------------------------------------------------------------------------------
# Estimate the 8 parameters by NLLS minimization

#define the objective function to optimize. Same as before

# find NLLS estimator
result_NLLS <-
  optim(par = theta, fn = compute_nlls_objective_pset2,
        method = "Nelder-Mead",
        df_share = df_share,
        X = X,
        M = M,
        V_mcmc = V_mcmc,
        e_mcmc = e_mcmc)
saveRDS(result_NLLS, file = glue("{folder}Pset2_result_NLLS.rds"))

# Compare the estimates with the true parameters.
result_NLLS <- readRDS(file = glue("{folder}Pset2_result_NLLS.rds"))
result <- data.frame(true = theta, estimates = result_NLLS$par)


# 3. The optimization takes about 4 minutes (in a mac laptop). 
# What is the mean x2 estimated marginal utility? 
# ANS: 0.2242372 
result$estimates[2]

# What is the estimated standard deviation of the price coeficient, omega?
# ANS: 0.9457209
result$estimates[8]

# 4. Graph the estimated alpha_i result
# First construct alpha_it from v_p_it. 
# Collect the first 500 observations in the simulated V_mcmc
alpha_df <- V_mcmc %>%
  slice(1:500) %>%
  mutate(alpha_i = -exp(result$estimates[7] + result$estimates[8] * v_p))
# Histogram of alpha_i
alpha_i_hist <- ggplot(alpha_df, aes(x = alpha_i)) +
  geom_histogram(binwidth =2.5, fill = "steelblue", color = "white") +
  labs(title = "Histogram of Estimated alpha_i", x = expression(alpha[i]), y = "Count")
ggsave(glue("{folder_output}PSet2_q3_4a.png"), alpha_i_hist)

# Scatter plot of alpha_i
alpha_i_scat <- ggplot(alpha_df, aes(x = i, y = alpha_i)) +
  geom_point(size = 1.5) +
  labs(
    x = "individuals",
    y = "Marginal Utility of Price alpha_i",
    title = "Scatter Plot of αᵢ Across Individuals"
  ) +
  theme_minimal()
ggsave(glue("{folder_output}PSet2_q3_4b.png"), alpha_i_scat)



