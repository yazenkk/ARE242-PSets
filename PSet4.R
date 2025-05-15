# ------------------------------------------------------------------------------
# ARE 242 Problem Set 4
# Spring 2025
# Author: Yazen Kashlan
# Date: 5/14/2025
# Note: this code is adapted from Villas-Boas lecture notes for ARE 242

# ------------------------------------------------------------------------------
options(scipen = 999)  # turn off scientific notation globally

# Setup
#install.packages("foreach")
#install.packages("evd")
#install.packages("ivreg")
#install.packages("pacman") 

# begin preamble of the R script  ------------------------------------------------------
# install and call needed packages
library(pacman) 
library(magrittr)
library(foreach)
library(latex2exp)
library(evd)

p_load(dplyr, haven, readr,cran) #using pacman manager
p_load(ggplot2)
p_load(AER,stargazer) #AER has canned ivreg, fyi, and stargazer for latex
p_load(sandwich,lmtest) #for robust standard errors
p_load(nleqslv) #for solving non linear systems of equations
p_load(np) #for non parametrics in Olley-Pakes method


#get rid of scientific display of numbers
options(scipen = 100, digits = 4)

folder <- "/Users/yazenkashlan/Documents/GitHub/personal-code/ARE242/"
folder_output <- "/Users/yazenkashlan/Documents/GitHub/personal-code/ARE242/output/"


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Question 1 - O-P Production Function Estimation
# ------------------------------------------------------------------------------
# Part 1. Simulate the data
set.seed(1234)
# The values of the auxiliary parameters are set as follows:
# size
J <- 1000
T <- 10
# production function
beta_0 <- 1
beta_l <- 0.2
beta_k <- 0.7
sigma_eta <- 0.2
# anticipated shocks
sigma_nu <- 0.5
alpha <- 0.7
# wage
sigma_w <- 0.1
# capital
delta <- 0.05

# ------------------------------------------------------------------------------
# 2. Write a function that returns the log output given 
# $l_{jt}$, $k_{jt}$, $\omega_{jt}$, and $\eta_{jt}$ 
# log production function
log_production <- function(l,  k, omega, eta, beta_0, beta_l, beta_k) {
  y <- beta_0 + beta_l * l + beta_k * k + omega + eta
  return(y)
}

# ------------------------------------------------------------------------------
# 3. Derive the optimal log labor as a function of omega_jt, eta_jt, k_jt, and wage.
# Define function for static labor choice without optimization error
# labor choice without iota (the error for manager when choosing labor)

log_labor_choice <-function(k, wage, omega, beta_0, beta_l,  beta_k, sigma_eta) {
  L <- ((1 / wage) *exp(beta_0 + omega + sigma_eta^2 / 2) *beta_l *exp(k)^beta_k)^(1 / (1 - beta_l))
  l <- log(L)
  return(l)
}


# ------------------------------------------------------------------------------
# 4. Modify the previous function by including an iota_jt in manager's perceived shock
# remember price p_jt is set to 1
# static labor choice with optimization error iota
log_labor_choice_error <- 
  function(k, wage, omega, beta_0, beta_l, beta_k, iota, sigma_eta) {
    l <- (
      (1/wage) * exp(beta_0 + omega + iota + sigma_eta^2/2) *  beta_l * exp(k)^beta_k
    )^(1/(1 - beta_l))
    l <- log(l)
    return(l)
  }
# optimization error
sigma_iota <- 0.05

# ------------------------------------------------------------------------------
# 5-6. Define the investment process
gamma <- 0.1
investment_choice <-function(k, omega, gamma, delta) {
  I <- (delta + gamma * omega) * exp(k)
  return(I)
}


# ------------------------------------------------------------------------------
# 7. Simulate the data first using the labor choice without optimization error 
# and second using the labor choice with optimization error.

# parameters for the initial state distribution
mean_k_initial <- 1
sigma_k_initial <- 0.5

# Draw omega_j1 from its stationary distribution (AR(1))
sigma_v_initial <- sigma_nu / sqrt(1 - alpha^2)

# Draw a wage
# draw initial state values
wage <- 0.5
set.seed(1234) # resetting seed for consistency with problem set notes

#price<-1 is normalized to 1
df <- 
  expand.grid(
    j = 1:J, 
    t = 1
  ) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    k = rnorm(J, mean_k_initial, sigma_k_initial),
    omega = rnorm(J, 0, sigma_v_initial),
    wage = wage
  )
head(df)


# ------------------------------------------------------------------------------
# 8. Draw optimization error iota $\iota_{jt}$ and compute the labor and investment choice. 

# compute labor and investment
df <- 
  df %>%
  dplyr::mutate(
    iota = rnorm(J, 0, sigma_iota)
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    l = log_labor_choice(k, wage, omega, beta_0, beta_l, beta_k, sigma_eta),
    l_error = log_labor_choice_error(k, wage, omega, beta_0, 
                                     beta_l, beta_k, iota, sigma_eta),
    I = investment_choice(k, omega, gamma, delta)
  ) %>%
  dplyr::ungroup()
head(df)


# ------------------------------------------------------------------------------
# 9. Draw ex post shock and compute the output according to the production function for both labor without optimization error and with optimization error. 

# compute output
df <- 
  df %>%
  dplyr::mutate(eta = rnorm(J, 0, sigma_eta)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    y = log_production(l, k, omega, eta, beta_0, beta_l, beta_k),
    y_error = log_production(l_error, k, omega, eta, beta_0, beta_l, beta_k) 
  ) %>%
  dplyr::ungroup()
head(df)


# ------------------------------------------------------------------------------
# 10. Repeat this procedure for all periods, $t = 1, \cdots 10$ by updating the 
# capital and anticipated shocks, and name the resulting data frame `df\_T`.

df_T <- df
for (t in 2:T) {
  # change time index
  df$t <- t
  
  # draw wage
  wage <- 0.5
  df$wage <- wage
  
  # update capital
  df <- 
    df %>%
    dplyr::mutate(
      k = log((1 - delta) * exp(k) + I) )
  
  # update omega
  df <- 
    df %>%
    dplyr::mutate(
      nu = rnorm(J, 0, sigma_nu),
      omega = alpha * omega + nu )
  
  # compute labor and investment
  df <- df %>%
    dplyr::mutate(iota = rnorm(J, 0, sigma_iota)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      l = log_labor_choice(k, wage, omega, beta_0, beta_l, beta_k, sigma_eta),
      l_error = log_labor_choice_error(k, wage, omega, beta_0, beta_l, beta_k, iota, sigma_eta),
      I = investment_choice(k, omega, gamma, delta)  ) %>%
    dplyr::ungroup()
  # compute output
  df <- df %>%
    dplyr::mutate(eta = rnorm(J, 0, sigma_eta)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      y = log_production(l, k, omega, eta, beta_0, beta_l, beta_k),
      y_error = log_production(l_error, k, omega, eta, beta_0, beta_l, beta_k)  ) %>%
    dplyr::ungroup()
  # append
  df_T <- dplyr::bind_rows(df_T, df)
}

save(df_T, file = glue("{folder_data}pset4df_T.RData"))

#you should get (top rows)
head(df_T)


# ------------------------------------------------------------------------------
# 11. Make a summary stats table of all the variables
summary(df_T)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Part 2- Estimation

# ------------------------------------------------------------------------------
# 1. First, simply regress $y_{jt}$ on $l_{jt}$ and $k_{jt}$ using  the least square method. This is likely to give an upwardly biased estimates on $\beta_l$. Why is it?
result_ols <-lm( data = df_T,formula = y_error ~ l_error + k) 
summary(result_ols)


# ------------------------------------------------------------------------------
# 2.  Second, Panel
df_T_within <- 
  df_T %>%
  dplyr::group_by(j) %>%
  dplyr::mutate(
    dy_error = y_error - mean(y_error),
    dl_error = l_error - mean(l_error),
    dk = k - mean(k)
  ) %>%
  dplyr::ungroup()
result_within <-lm(data = df_T_within,formula = dy_error ~ - 1 + dl_error + dk)
summary(result_within)


# ------------------------------------------------------------------------------
# 3. Estimate the first-step model of Olley-Pakes method:

# obtain optimal bandwidth
result_1st_bw <-npplregbw(data = df_T,
                          formula = y_error ~ l_error + k + I| k + I)
save(result_1st_bw, file = glue("{folder_data}pset4_result_1st_bw.RData"))
load(file = glue("{folder_data}pset4_result_1st_bw.RData"))

#estimate model under optimal bandwidth
result_1st <-npplreg(data = df_T,
                     formula = y_error ~ l_error + k + I| k + I,
                     bws = result_1st_bw)
save(result_1st, file = glue("{folder_data}pset4_result_1st.RData"))

# Return the summary of the first stage estimation and plot the 
# fitted values against the data points.
load(file = glue("{folder_data}pset4df_T.RData"))
load(file = glue("{folder_data}pset4_result_1st.RData"))
summary(result_1st)
result_1st_plot <-data.frame( actual = df_T$y_error,fitted = fitted(result_1st))
PSet4_q2_3 <- result_1st_plot %>%
  ggplot(aes(x = fitted, y = actual) ) +geom_point()
ggsave(glue("{folder_output}PSet4_q2_3.png"), PSet4_q2_3)


# ------------------------------------------------------------------------------
# 4. Estimate the second-step of Olley-Pakes method
df_T_1st <- df_T %>%
  dplyr::mutate(
    y_error_tilde = y_error - result_1st$xcoef["l_error"] * l_error ) %>%
  dplyr::select( j,  t, y_error_tilde)
phi <- fitted(result_1st) - 
  result_1st$xcoef["l_error"] * df_T$l_error
df_T_1st$phi_t <- phi
df_T_1st <- 
  df_T_1st %>%
  dplyr::arrange(j, t ) %>%
  dplyr::group_by(j) %>%
  dplyr::mutate(phi_t_1 = dplyr::lag(phi_t, 1)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-phi_t)
head(df_T_1st)


# ------------------------------------------------------------------------------
# 5. Compute a function that returns the value of $g_{JT}(\alpha, \beta_0, \beta_k)$ 
# given parameter values, data, and `df\_T\_1st`, and name it `moment\_OP\_2nd`. 

moment_OP_2nd <-  function(alpha, beta_0, beta_k, df_T, df_T_1st ) {
  moment <- 
    df_T %>%
    dplyr::group_by(j) %>%
    dplyr::mutate(
      k_t_1 = dplyr::lag(k, 1),
      I_t_1 = dplyr::lag(I, 1)      ) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(df_T_1st, by = c("j", "t") ) %>%
    dplyr::filter(!is.na(phi_t_1)) %>%
    dplyr::mutate(
      g_jt = y_error_tilde - beta_0 - beta_k * k - 
        alpha * (phi_t_1 - beta_0 - beta_k * k_t_1)       )
  moment <-     cbind(  moment$g_jt * moment$k,
                        moment$g_jt * moment$k_t_1,
                        moment$g_jt * moment$I_t_1 )
  moment <- apply(moment, 2, mean)
  return(moment)
}

moment_OP_2nd(alpha, beta_0, beta_k, df_T, df_T_1st)


# ------------------------------------------------------------------------------
# 6. Write a function that returns the value of objective Function
objective_OP_2nd <-function(theta, df_T, df_T_1st, W) {
  alpha <- theta[1]
  beta_0 <- theta[2]
  beta_k <- theta[3]
  moment <-moment_OP_2nd(alpha, beta_0, beta_k, df_T, df_T_1st)
  objective <- t(moment) %*% W %*% moment
  return(objective)
}

#Setting $W$ at the identity matrix, show the value of the 
#objective function evaluated at the true parameters.
W <- diag(3)

theta <- c(alpha, beta_0, beta_k)

objective_OP_2nd(theta, df_T, df_T_1st, W)


# ------------------------------------------------------------------------------
# 7. Find the parameters that minimize the objective function using `optim`.

theta <- c(alpha, beta_0, beta_k)
W <- diag(3)
result_2nd <-
  optim(
    par = theta,
    f = objective_OP_2nd,
    method = "L-BFGS-B",
    df_T = df_T,
    df_T_1st = df_T_1st,
    W = W
  )
result_2nd


# ------------------------------------------------------------------------------
# 8. Draw a graph for the objective function when one of the parameters varies and all others are at their true values
objective_alpha <-
  foreach ( i = seq(0, 1, 0.1),
            .combine = "rbind") %do% {
              objective_i <-
                objective_OP_2nd(c(i, beta_0, beta_k), df_T, df_T_1st, W)
              return(objective_i)
            }
objective_beta_0 <-
  foreach (i = seq(0, 1, 0.1),
           .combine = "rbind") %do% {
             objective_i <-
               objective_OP_2nd(c(alpha, i, beta_k), df_T, df_T_1st, W)
             return(objective_i)
           }
objective_beta_k <-
  foreach (i = seq(0, 1, 0.1),
           .combine = "rbind") %do% {
             objective_i <-
               objective_OP_2nd(c(alpha, beta_0, i), df_T, df_T_1st, W)
             return(objective_i)
           }

objective_plot <-
  data.frame(
    i = seq(0, 1, 0.1),
    alpha = objective_alpha,
    beta_0 = objective_beta_0,
    beta_k = objective_beta_k )

alpha_opt <- objective_plot %>%
  ggplot(aes( x = i, y = alpha )) +geom_point() +
  xlab(TeX("$\\alpha$")) +
  ylab("Objective")
ggsave(glue("{folder_output}PSet4_q2_8a.png"), alpha_opt)

beta0_opt <- objective_plot %>%
  ggplot(aes(x = i, y = beta_0)) +geom_point() +
  xlab(TeX("$\\beta_0$")) +
  ylab("Objective")
ggsave(glue("{folder_output}PSet4_q2_8b0.png"), beta0_opt)

betak_opt <- objective_plot %>%
  ggplot( aes(x = i, y = beta_k)) +geom_point() +
  xlab(TeX("$\\beta_k$")) +
  ylab("Objective")
ggsave(glue("{folder_output}PSet4_q2_8bk.png"), betak_opt)

