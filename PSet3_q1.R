# ------------------------------------------------------------------------------
# ARE 242 Problem Set 3
# Spring 2025
# Author: Yazen Kashlan
# Date: 5/12/2025
# 
# Outline: This script 
# - Sets ups parameters: beta, alpha, omega, price, sigma 
# - creates dataframes: X, M, V, e 
# - creates functions which takes in params, data, (and preference shocks):
#     (1) - compute_indirect_utility which defines utility and alpha_i as
#             rowSums(beta_i * X) + p * alpha_i + xi
#             alpha_i <- - exp(a + omega * v_p)
#     (2) - compute_choice which and calls (1)  
#             where q = ifelse(u + e == max(u + e), 1, 0)
#     (2') - compute_choice_smooth which instead computes q \in (0,1)
#             q = exp(u)/sum(exp(u))
#     (3) - compute_share which calls (2)
#     (4) - compute_nlls_objective which calls (3)
#
# Account for endogenous prices (BLP)
# Repeat the above with 
#     (1) - compute_indirect_utility_delta which defines utility and alpha_i as
#             rowSums(beta_i * X) + p * tildealpha_i + delta_ijt
#             tildealpha_i <- - exp(a + omega * v_p) - (- exp(a + omega^2/2))
#     (2') - compute_choice_smooth_delta which instead computes q \in (0,1)
#             q = exp(u)/sum(exp(u))
#
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
#get rid of scientific display of numbers
options(scipen = 100, digits = 4)


#install.packages("foreach")
#install.packages("evd")
#library(EmpiricalIO)
#library(magrittr)
#library(foreach)
#library(ggplot2)
#library(latex2exp)
#library(evd)

rm(list=ls())
set.seed(1234)

folder <- "/Users/yazenkashlan/Documents/GitHub/personal-code/ARE242/"
folder_output <- "/Users/yazenkashlan/Documents/GitHub/personal-code/ARE242/output/"
folder_data <- "/Users/yazenkashlan/Documents/GitHub/personal-code/ARE242/data/"
folder_tables <- "/Users/yazenkashlan/Documents/GitHub/personal-code/ARE242/tables/"


# ------------------------------------------------------------------------------
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

#Note : the constant beta01 is set to 4 
#and other two beta0s and drawn; 
#Note, there are k=3, three sigmas, one for each k
beta <- rnorm(K); 
beta[1] <- 4      
beta
sigma <- abs(rnorm(K));    
sigma
#displayed as output after this cell is executed are betas first and then sigma of betas

# set auxiliary parameters
price_xi <- 1
sd_x <- 2
sd_xi <- 0.5
sd_c <- 0.05
sd_p <- 0.05



# ------------------------------------------------------------------------------
# make product characteristics data
X <- matrix(sd_x * rnorm(J * (K - 1)), nrow = J)
X <- cbind(rep(1, J), X)
colnames(X) <- paste("x", 1:K, sep = "_")
X <- data.frame(j = 1:J, X) %>%
  tibble::as_tibble()
# add outside option
X <- rbind(rep(0, dim(X)[2]), X) 
head(X)

# make market-product data
M <- expand.grid(j = 1:J, t = 1:T) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    xi = sd_xi * rnorm(J*T),
    c = exp(sd_c * rnorm(J*T)),
    p = exp(price_xi * xi + sd_p * rnorm(J*T)) + c
  )
M <- M %>%
  dplyr::group_by(t) %>%
  dplyr::sample_frac(size = purrr::rdunif(1, J)/J) %>%
  dplyr::ungroup()
# add outside option
outside <- data.frame(j = 0, t = 1:T, xi = 0, c = 0, p = 0)
M <- rbind(M, outside) %>%
  dplyr::arrange(t, j)
M

# make consumer-market data
V <- matrix(rnorm(N * T * (K + 1)), nrow = N * T) 
colnames(V) <- c(paste("v_x", 1:K, sep = "_"), "v_p")
V <- data.frame(expand.grid(i = 1:N, t = 1:T), V) %>%
  tibble::as_tibble()
head(V)


# To make choice data, lets create df
df <- expand.grid(t = 1:T, i = 1:N, j = 0:J) %>%
  tibble::as_tibble() %>%
  dplyr::left_join(V, by = c("i", "t")) %>%
  dplyr::left_join(X, by = c("j")) %>%
  dplyr::left_join(M, by = c("j", "t")) %>%
  dplyr::filter(!is.na(p)) %>%
  dplyr::arrange(t, i, j)
head(df)


# ------------------------------------------------------------------------------
# Define functions
compute_indirect_utility <-  function( df, beta,  sigma, a, omega ) {
  #...  see R file with functions for you to use in Bcourses Pset 2 folder
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
  
  return(value)  }

# compute indirect utility
u <-  compute_indirect_utility(df, beta, sigma, a, omega)



compute_choice_smooth <-
  function(X,M, V, beta, sigma, a, omega) {
    # constants
    T <- max(M$t)
    N <- max(V$i)
    J <- max(X$j)
    # make choice data
    df <- expand.grid(t = 1:T,  i = 1:N, j = 0:J) %>%
      tibble::as_tibble() %>%
      dplyr::left_join(V, by = c("i", "t")) %>%
      dplyr::left_join(X, by = c("j")) %>%
      dplyr::left_join(M, by = c("j", "t")) %>%
      dplyr::filter(!is.na(p)) %>%
      dplyr::arrange(t, i, j)
    # compute indirect utility
    u <- compute_indirect_utility(df, beta, sigma, a, omega)
    # add u 
    df_choice <- data.frame(df, u) %>%
      tibble::as_tibble()
    # make choice
    df_choice <- 
      df_choice %>%
      dplyr::group_by(t, i) %>%
      dplyr::mutate(q = exp(u)/sum(exp(u))) %>%
      dplyr::ungroup()
    # return
    return(df_choice)
  }

df_choice_smooth <-compute_choice_smooth(X, M, V, beta, sigma, a, omega)
summary(df_choice_smooth)



compute_share_smooth <-function(X, M, V, beta, sigma, a, omega) {
  # constants
  T <- max(M$t)
  N <- max(V$i)
  J <- max(X$j)
  # compute choice
  df_choice <- compute_choice_smooth(X, M, V, beta, sigma,a, omega)
  # make share data
  df_share_smooth <- 
    df_choice %>%
    dplyr::select(-dplyr::starts_with("v_"), -u, -i ) %>%
    dplyr::group_by(t, j) %>%
    dplyr::mutate(q = sum(q)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(t, j, .keep_all = TRUE) %>%
    dplyr::group_by(t) %>%
    dplyr::mutate(s = q/sum(q)) %>%
    dplyr::ungroup()
  # log share difference
  df_share_smooth <- 
    df_share_smooth %>%
    dplyr::group_by(t) %>%
    dplyr::mutate(y = log(s/sum(s * (j == 0)))) %>%
    dplyr::ungroup()
  return(df_share_smooth)
}


df_share_smooth <- compute_share_smooth(X, M, V, beta, sigma, a, omega)
head(df_share_smooth)


## draw V and V_mcmc
V_mcmc <- matrix(rnorm(L*T*(K + 1)), nrow = L*T)
colnames(V_mcmc) <- c(paste("v_x", 1:K, sep = "_"), "v_p")
V_mcmc <- data.frame(
  expand.grid(i = 1:L, t = 1:T),
  V_mcmc
) %>%
  tibble::as_tibble()
V_mcmc

df_mcmc <- expand.grid(t = 1:T, i = 1:L, j = 0:J) %>%
  tibble::as_tibble() %>%
  dplyr::left_join(V_mcmc, by = c("i", "t")) %>%
  dplyr::left_join(X, by = c("j")) %>%
  dplyr::left_join(M, by = c("j", "t")) %>%
  dplyr::filter(!is.na(p)) %>%
  dplyr::arrange(t, i, j)

## draw mcmc e
# draw idiosyncratic shocks
e_mcmc <- evd::rgev(dim(df_mcmc)[1])


#set parameters
theta <- c(beta, sigma, a, omega)

# get starting values for NLLS step below
M_no <- M %>%
  dplyr::mutate(xi = 0)


# compute share
compute_share <-function(X,  M,  V, e, beta, sigma,  a, omega) {
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

# compute choice
compute_choice <- function( X, M, V,  e, beta,  sigma, a, omega) {
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
  u <-  compute_indirect_utility( df, beta,  sigma, a, omega)
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

# nlls objective function
compute_nlls_objective <- function(theta, df_share, X, M, V_mcmc, e_mcmc) {
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
# Get starting values for case with xi != 0 using market data where xi = 0
# find NLLS estimator
result_NLLS <- 
  optim(par = theta, fn = compute_nlls_objective,
        method = "Nelder-Mead",
        df_share = df_share_smooth, 
        X = X, 
        M = M_no, 
        V_mcmc = V_mcmc, 
        e_mcmc = e_mcmc)
save(result_NLLS, file = glue("{folder_data}result_NLLS.RData_lecture5"))

load(file = glue("{folder_data}result_NLLS.RData_lecture5"))
result_NLLS
result <- data.frame(true = theta, estimates = result_NLLS$par)
result

save(df_share_smooth, file = glue("{folder_data}df_share_smooth.RData"))


# ------------------------------------------------------------------------------
# BLP
load(file = glue("{folder_data}df_share_smooth.RData"))
head(df_share_smooth)

# compute delta at the true parameters: delta = beta*x + alpha*p + xi
XX <- as.matrix(dplyr::select(df_share_smooth, dplyr::starts_with("x_")))
pp <- as.matrix(dplyr::select(df_share_smooth, p)) 
xi <- as.matrix(dplyr::select(df_share_smooth, xi))
alpha0 <- - exp(a + omega^2/2)
delta <- XX %*% as.matrix(beta) + pp * alpha0 + xi
delta <- dplyr::select(df_share_smooth, t, j) %>%
  dplyr::mutate(delta = as.numeric(delta))
head(delta)

# compute indirect utility from delta
# first define the function
compute_indirect_utility_delta <- function(df,delta, sigma, a,omega) {
  # extract matrices
  X <- as.matrix(dplyr::select(df, dplyr::starts_with("x_")))
  p <- as.matrix(dplyr::select(df, p)) 
  v_x <- as.matrix(dplyr::select(df, dplyr::starts_with("v_x")))
  v_p <- as.matrix(dplyr::select(df, v_p))
  # expand delta
  delta_ijt <- df %>%
    dplyr::left_join(delta, by = c("t", "j")) %>%
    dplyr::select(delta) %>%
    as.matrix()
  # random coefficients
  beta_i <- v_x %*% diag(sigma) 
  tildealpha_i <- - exp(a + omega * v_p) - (- exp(a + omega^2/2))
  # conditional mean indirect utility
  value <- as.matrix(delta_ijt + rowSums(beta_i * X) + p * tildealpha_i) 
  colnames(value) <- "u"
  return(value)
}

u_delta <-compute_indirect_utility_delta(df, delta, sigma,a, omega)
head(u_delta)
summary(u - u_delta)

# compute choice from delta
compute_choice_smooth_delta <-function(X, M, V, delta, sigma, a, omega) {
  # constants
  T <- max(M$t)
  N <- max(V$i)
  J <- max(X$j)
  # make choice data
  df <- expand.grid(t = 1:T,  i = 1:N,  j = 0:J) %>%
    tibble::as_tibble() %>%
    dplyr::left_join(V, by = c("i", "t")) %>%
    dplyr::left_join(X,  by = c("j")) %>%
    dplyr::left_join( M, by = c("j", "t")) %>%
    dplyr::filter(!is.na(p)) %>%
    dplyr::arrange(t, i, j)
  # compute indirect utility
  u <- compute_indirect_utility_delta(df, delta, sigma, a, omega)
  # add u 
  df_choice <- data.frame(df, u) %>%
    tibble::as_tibble()
  # make choice
  df_choice <- df_choice %>%
    dplyr::group_by(t, i) %>%
    dplyr::mutate(q = exp(u)/sum(exp(u))) %>%
    dplyr::ungroup()
  # return
  return(df_choice)
}


df_choice_smooth_delta <- 
  compute_choice_smooth_delta(X, M, V, delta, sigma, a, omega)
head(df_choice_smooth_delta)

#compare data with the estimated one (given true parameters)
summary(df_choice_smooth$q - df_choice_smooth_delta$q)

# compute share from delta
compute_share_smooth_delta <-function(X, M, V, delta, sigma, a, omega) {
  # constants
  T <- max(M$t)
  N <- max(V$i)
  J <- max(X$j)
  # compute choice
  df_choice <- compute_choice_smooth_delta(X, M, V, delta, sigma,a, omega)
  # make share data
  df_share_smooth <-
    df_choice %>%
    dplyr::select(-dplyr::starts_with("v_"),-u,-i) %>%
    dplyr::group_by(t, j) %>%
    dplyr::mutate(q = sum(q)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(t, j, .keep_all = TRUE) %>%
    dplyr::group_by(t) %>%
    dplyr::mutate(s = q/sum(q)) %>%
    dplyr::ungroup()
  # log share difference
  df_share_smooth <- 
    df_share_smooth %>%
    dplyr::group_by(t) %>%
    dplyr::mutate(y = log(s/sum(s * (j == 0)))) %>%
    dplyr::ungroup()
  return(df_share_smooth)
}

#Now compute share
df_share_smooth_delta <-
  compute_share_smooth_delta(X, M, V, delta, sigma, a, omega) 
head(df_share_smooth_delta)
summary(df_share_smooth$s - df_share_smooth_delta$s)

solve_delta <-function(df_share_smooth, X, M, V, delta, sigma, a, 
                       omega, kappa, lambda) {
  # initial distance
  distance <- 10000
  # fixed-point algorithm
  delta_old <- delta
  while (distance > lambda) {
    # save the old delta
    delta_old$delta <- delta$delta
    # compute the share with the old delta
    df_share_smooth_predicted <- compute_share_smooth_delta(X, M, V, delta_old, sigma, a, omega )  
    # update the delta
    delta$delta <-  delta_old$delta + 
      (log(df_share_smooth$s) - log(df_share_smooth_predicted$s)) * kappa
    delta <-  delta %>%
      dplyr::mutate(delta = ifelse(j == 0, 0, delta))
    # update the distance
    distance <- max(abs(delta$delta - delta_old$delta))
    # print(distance)
  }
  return(delta)
}


# 2. Start the algorithm with the true delta_jt and check if the algorithm returns (almost) the same delta_jt when
# the actual and predicted smooth share are equated.
kappa <- 1
lambda <- 1e-3
delta_new <- solve_delta(df_share_smooth, X, M, V, delta, sigma, a, omega, kappa, lambda)
head(delta_new)
summary(delta_new$delta - delta$delta) # pretty good


# 3. Check how long it takes to compute the limit delta under the Monte Carlo shocks starting from the true
# to match with `df share smooth`. This is approximately the time to evaluate the objective function.
delta_new <- solve_delta(df_share_smooth, X, M, V_mcmc, delta, sigma, a,
              omega, kappa, lambda)

save(delta_new, file = glue("{folder_data}Lecture6_delta_new.RData"))
load(file = glue("{folder_data}Lecture6_delta_new.RData"))
summary(delta_new$delta - delta$delta)


# 4. Compute optimal linear parameters associated with data and delta
#Psi is the identity matrix (as weighting matrix 
Psi <- diag(length(beta) + 1)


#Get the optimal linear parameters
compute_theta_linear <-function(df_share_smooth, delta, a, omega, Psi) {
  # extract matrices
  X <-  df_share_smooth %>%
    dplyr::filter(j != 0) %>%
    dplyr::select(dplyr::starts_with("x_")) %>%
    as.matrix()
  p <- df_share_smooth %>%
    dplyr::filter(j != 0) %>%
    dplyr::select(p) %>%
    as.matrix()
  W <- df_share_smooth %>%
    dplyr::filter(j != 0) %>%
    dplyr::select(dplyr::starts_with("x_"), c) %>%
    as.matrix()
  delta_m <- delta %>%
    dplyr::filter(j != 0) %>%
    dplyr::select(delta) %>%
    as.matrix()
  alpha <- - exp(a + omega^2/2)
  # compute the optimal linear parameters
  theta_linear_1 <-crossprod(X, W) %*% solve(Psi, crossprod(W, X))
  theta_linear_2 <- crossprod(X, W) %*% solve(Psi, crossprod(W, delta_m - alpha * p))
  theta_linear <-  solve(theta_linear_1, theta_linear_2)
  return(theta_linear)
}

theta_linear <- compute_theta_linear(df_share_smooth, delta, a, omega, Psi) 
cbind(theta_linear, beta)


# solve xi associated with delta and linear parameters
solve_xi <-function(df_share_smooth, delta, beta, a, omega) {
  # extract matrices
  X1 <- df_share_smooth %>%
    dplyr::filter(j != 0) %>%
    dplyr::select(dplyr::starts_with("x_"),  p) %>%
    as.matrix()
  delta_m <- delta %>%
    dplyr::filter(j != 0) %>%
    dplyr::select(delta) %>%
    as.matrix()
  alpha <- - exp(a + omega^2/2)
  theta_linear <- c(beta, alpha)
  # compute xi
  xi <- delta_m - X1 %*% theta_linear
  colnames(xi) <- "xi"
  # return
  return(xi)
}


xi_new <- solve_xi(df_share_smooth, delta, beta, a, omega)
head(xi_new)
xi_true <-
  df_share_smooth %>%
  dplyr::filter(j != 0) %>%
  dplyr::select(xi)
summary(xi_true - xi_new)

# compute GMM objective function
compute_gmm_objective_Lecture6 <-
  function(theta_nonlinear, delta, df_share_smooth, Psi,X, M, V_mcmc, kappa, 
           lambda) {
    # exctract parameters
    a <- theta_nonlinear[1]
    omega <- theta_nonlinear[2]
    sigma <- theta_nonlinear[3:length(theta_nonlinear)]
    # extract matrix
    W <- df_share_smooth %>%
      dplyr::filter(j != 0) %>%
      dplyr::select(dplyr::starts_with("x_"), c) %>%
      as.matrix()
    # compute the delta that equates the actual and predicted shares
    delta <- solve_delta(
      df_share_smooth, X, M, V_mcmc, delta, sigma, a, 
      omega, kappa, lambda)
    # compute the optimal linear parameters
    beta <-compute_theta_linear(
      df_share_smooth,  delta,  a, omega, Psi) 
    # compute associated xi
    xi <-  solve_xi(df_share_smooth, delta, beta, a, omega)
    # compute objective
    objective <- crossprod(xi, W) %*% solve(Psi, crossprod(W, xi))
    # return
    return(objective)
  }



# non-linear parmaeters
theta_nonlinear <- c(a, omega, sigma)
theta_nonlinear

# Compute GMM objective function
objective <-compute_gmm_objective_Lecture6(theta_nonlinear, delta, df_share_smooth, Psi, 
                           X, M, V_mcmc, kappa, lambda)

save(objective, file = glue("{folder_data}Lecture6_objective.RData"))
load(file = glue("{folder_data}Lecture6_objective.RData"))
objective


# BLP Step
#set timer
ptm <- proc.time()
#result <- optim(par = theta_nonlinear,
#                fn = compute_gmm_objective_Lecture6,
#                method = "BFGS",
#                delta = delta, 
#                df_share_smooth = df_share_smooth, 
#                Psi = Psi, 
#                X = X, M = M, 
#                V_mcmc = V_mcmc, 
#                kappa = kappa,
#                lambda = lambda)
#save(result, file = glue("{folder_data}Lecture6_result.RData"))
load(file = glue("{folder_data}Lecture6_result.RData"))
result
comparison <- cbind(theta_nonlinear, abs(result$par))
colnames(comparison) <- c("true", "estimate")
comparison
proc.time() - ptm
# 6+ hours for me

# ------------------------------------------------------------------------------
# Question 1 
# 1. What is the standard deviation of the attribute x2?
# ANS: 0.4402
result$par[4]

# 2. Given the non linear parameters, what is the value of the GMM objective function?
# ANS: about zero
result$value

# 3. Please graph the estimated alpha_i coefficients. And the histogram of alpha_i
# Collect the first 500 observations in the simulated V_mcmc
alpha_df <- V_mcmc %>%
  slice(1:500) %>%
  mutate(alpha_i = -exp(result$par[1] + result$par[2] * v_p))
# Histogram of alpha_i
alpha_i_hist <- ggplot(alpha_df, aes(x = alpha_i)) +
  geom_histogram(binwidth =2.5, fill = "steelblue", color = "white") +
  labs(title = "Histogram of Estimated alpha_i", x = expression(alpha[i]), y = "Count")
ggsave(glue("{folder_output}PSet3_q1_3a.png"), alpha_i_hist)
# Scatter plot of alpha_i
alpha_i_scat <- ggplot(alpha_df, aes(x = i, y = alpha_i)) +
  geom_point(size = 1.5) +
  labs(
    x = "Individuals",
    y = "Marginal Utility of Price alpha_i",
    title = "Scatter Plot of αᵢ Across Individuals"
  ) +
  theme_minimal()
ggsave(glue("{folder_output}PSet2_q3_4b.png"), alpha_i_scat)


# 4. Please estimate the linear parameters. (beta_0 and alpha_0)
alpha_0 = -exp(result$par[1] + result$par[2]^2/2)
alpha_0
beta_0 = compute_theta_linear(df_share_smooth, delta, result$par[1], result$par[2], Psi)  
beta_0
cbind(beta, beta_0) # compare with original betas

#     beta   beta_0
# x_1 4.0000 4.1648
# x_2 0.2774 0.2746
# x_3 1.0844 1.0826


# What is the average marginal utility for the attribute x2
# 0.2746
beta_0[2]
