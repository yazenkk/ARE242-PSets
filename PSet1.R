# ------------------------------------------------------------------------------
# ARE 242 Problem Set 1
# Spring 2025
# Author: Yazen Kashlan
# Date: 3/8/2025
# 
# Outline: This script 
# - Adapts the script fish.R - Lecture 1
# - Generates figures: PSet1_q2_6.png
# 
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
#install and call needed packages
#install.packages("evd") 
#install.packages("pacman") 
#install.packages("foreach") 
library(evd)
library(pacman) 
library(foreach) # for %do%
p_load(dplyr, haven, readr,cran) 
p_load(ggplot2)
p_load(AER,stargazer)
p_load(glue)

#INSTRUMENTAL VARIABLE IV AND 2SLS – CASE OF ESTIMATING DEMAND IN THE NY FISH MARKET
#First stage
#White Robust standard errors
#Newey West robust std errors
#in fish.R now
rm(list=ls())


#-------------------------------------------
# Read in data
folder <- "/Users/yazenkashlan/Documents/GitHub/personal-code/ARE242/"
folder_output <- "/Users/yazenkashlan/Documents/GitHub/personal-code/ARE242/output/"
mydata <- read_dta(paste(folder, "fishdata.dta", sep=""))


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
## Question 1 part 1
# 1. Report the estimates of the first stage of specification in column (4) and 
# discuss the the R squared and the F test of significance of the weather instrument.
# Do we have a weak instrument problem?

# Replication column 4's IV estimates
reg4<-ivreg(qty~price + day1 + day2 + day3 + day4 + cold + rainy | 
                stormy + day1 + day2 + day3 + day4 + cold + rainy, 
            data=mydata)
summary(reg4, diagnostics = TRUE)

#test weak instrument
regfirststage4<-lm(price~stormy + day1 + day2 + day3 + day4 + cold + rainy, mydata)
summary(regfirststage4)
waldtest(regfirststage4, . ~ . - stormy, test = "F")  # F-test for instrument relevance
#F-statistic = 14.612
# ANS: > 10 (benchmark for weak instruments) hence no weak instrument problem.

# Partial R^2 of instrument Z
SSR_full <- sum(resid(lm(price~day1 + day2 + day3 + day4 + cold + rainy, 
                         data = mydata))^2)  # Sum of squared residuals without Z
SSR_reduced <- sum(resid(regfirststage4)^2)  # With Z
partial_R2 <- 1 - (SSR_reduced / SSR_full)
print(partial_R2)
# ANS: (adjusted) R^2 is lower than regular R^2 as expected given lower F-stat.
# ANS: No weak instrument problem.


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
## Question 1 part 2 - add 2 instruments 
# 2. Instead of generating a dummy for stormy, use as instruments the raw wind 
# speed data and the square of wind speed as instruments in addition to the stormy 
# instrument. Repeat the instrumental variable approach to estimate the coe cient 
# on log price (the elasticity). Do you have a weak instrument problem?


#column 5: IV estimates with windspeed and its square in addition to stormy indicator
reg5<-ivreg(qty~price + day1 + day2 + day3 + day4 + cold + rainy | 
                stormy + windspd + windspd2 + day1 + day2 + day3 + day4 + cold + rainy, 
            data=mydata)
summary(reg5, diagnostics = TRUE)

#test weak instrument
regfirststage5<-lm(price~stormy + windspd + windspd2 + day1 + day2 + day3 + day4 + cold + rainy, 
                   mydata)
summary(regfirststage5)
waldtest(regfirststage5, . ~ . - stormy-windspd - windspd2, test = "F")  # F-test for instrument relevance
#F-statistic: 6.677
# ANS: 
# F < 10 hence we have a weak instrument. 

#column 5.1: try using just windspeed as an IV
reg5_1<-ivreg(qty~price + day1 + day2 + day3 + day4 + cold + rainy |
                windspd + day1 + day2 + day3 + day4 + cold + rainy, 
              data=mydata)
summary(reg5_1, diagnostics = TRUE)

#test weak instrument
regfirststage5_1<-lm(price~windspd + day1 + day2 + day3 + day4 + cold + rainy, 
                     mydata)
summary(regfirststage5_1)
waldtest(regfirststage5_1, . ~ . - windspd, test = "F")  # F-test for instrument relevance
#F-statistic: 17.31
# ANS: Here there is evidence of a many-instruments problem.
# Even though each instrument works on its own, including them all is not recommended.



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
## Question 2: Simulate data - ML Logit estimation exogenous price

# 1. Make the dataframe with 1000 observations for each k; k = 1 or 2; x1 = 0; x2 = 1
df <-
  expand.grid( i = 1:1000, k = 1:2) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(x = ifelse(k == 1, 0, 1)) %>%
  dplyr::arrange(i, k)
df

## 2. Add price and the error term to the dataframe
set.seed(1234)
df <-
  df %>%
  #First, add a random (exogenous) price
  dplyr::mutate(price = runif(dim(df)[1])) %>%
  #Second, draw type-I extreme value random variables.
  dplyr::mutate(e = evd::rgev(dim(df)[1]))
df

# 3. Third, compute the latent value of each option for alpha = -0:6; beta = 0:2
beta <- 0.2
alpha<- -0.6
theta<-c(alpha,beta)
df <-
  df %>%
  dplyr::mutate(latent = alpha*price+ beta * x + e)
df

# 4. Finally, compute y (the choice 0 or 1) by comparing the latent values of k = 1; 2 for each i to obtain the following
df <-
  df %>%
  dplyr::group_by(i) %>%
  dplyr::mutate(y = ifelse(latent == max(latent), 1, 0)) %>%
  dplyr::ungroup()
df


# 5. Generate log-likelihood
loglikelihood_quest2 <-function( temp, df )
{
  l <- df %>%
    dplyr::group_by(i) %>%
    dplyr::mutate(p = exp(temp[1]*price+temp[2]*x)/sum(exp(temp[1]*price+temp[2]*x))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(y == 1)
  l <- mean(log(l$p))
  return(l)
}
ltry <- loglikelihood_quest2( c(-0.3,0.2), df )
ltry

# 6. Compute logL for different betas at true alpha
alpha <- -0.6
b_seq <- seq(0, 1, 0.1)
output <-
  foreach ( b = b_seq, .combine = "rbind" ) %do% {
    l <-
      loglikelihood_quest2(c(alpha, b), df)
    return(l)
  }
#make graph
output <- data.frame(x = b_seq, y = output )
output %>%
  ggplot( aes( x = x, y = y) ) + geom_point() + xlab(("beta")) + ylab("Loglikelihood")

# repeat for alpha at true beta = 0.2
beta = 0.2
a_seq <- seq(-1, 0, 0.1)
output <-
  foreach ( a = a_seq, .combine = "rbind" ) %do% {
    l <-
      loglikelihood_quest2(c(a, beta), df)
    return(l)
  }
#make graph
output <- data.frame(x = a_seq, y = output )
plot_6_1 <- output %>%
          ggplot( aes( x = x, y = y) ) + geom_point() + xlab(("alpha")) + ylab("Loglikelihood")
print(plot_6_1)
ggsave(glue("{folder_output}PSet1_q2_6.png"), plot_6_1)

# 7. Find and report alpha and beta that maximizes the log likelihood for the simulated data. What is the value of the log likelihood function.
result <-
  optim( par = theta, fn = loglikelihood_quest2, df = df, control = list(fnscale = - 1))
result


# 8. Test whether alpha=-1. Find the optimal beta with restriction: alpha = -1.
# Do you reject the null at the 10 percent signi cance level?

alpha_r<- -1
loglikelihood_onlyb <-function( tempb, df ) {
  l <- df %>%
    dplyr::group_by(i) %>%
    dplyr::mutate(p = exp(alpha_r*price+tempb*x)/sum(exp(alpha_r*price+tempb*x))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(y == 1)
  l <- mean(log(l$p))
  return(l) 
}
result_onlyb <-
  optim( par = b, fn = loglikelihood_onlyb,df = df, 
         method = "Brent", lower = -1, upper = 1, 
         control = list(fnscale = - 1))
result_onlyb

# ANS: Compare restricted and unrestricted model for likelihood ratio test.
# Under the null hypothesis (that α = -1), the LR test statistic is asymptotically
# chi-squared distributed with degrees of freedom equal to the number of 
# restrictions (here, 1).
# Compute the LR test stat and get its p-value.
LR_stat <- 2 * (result$value - result_onlyb$value)
p_value <- 1 - pchisq(LR_stat, df = 1)
p_value
# Since the p-value is greater than 0.1, we cannot reject the null that alpha = -1.



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
## Question 3
# Simulate data - ML Logit estimation endogenous price
# We have a new choice model. Let:

# util=alpha*p_ik + beta*x_ik + eps_ik (previously)
# eps_ik ~EVD
# p_ik ~RUNIF[0,1]

# util=alpha*priceE_ik + beta*x_ik + ksi_ik + eps_ik (now, with:)
# priceE = p_ik + x_ik + ksi_ik + v_ik (endogenous)
# ksi_ik ~RUNIF[0,1] (source of endogeneity)
# p_ik ~RUNIF[0,1]
# v_ik ~RUNIF[0,1]

# 1. Make a dataframe with 1000 observations for each k; k = 1 or 2; x1 = 0; 
# x2 = 1 and add an endogenous price
# 1.1 define an omitted variable Xe (the Xes in logit demand)

#define the model of priceE with Xe
#set.seed(1234)
df <-
  df %>%
  dplyr::mutate(Xe = runif(dim(df)[1])) %>%
  dplyr::mutate(priceE = price +Xe+runif(dim(df)[1]))
df

# 2. Generate the new latent2 and y2 choices given Xe and priceE
beta <- 0.2
alpha<- -0.6
theta<-c(alpha,beta)
df <-
  df %>%
  dplyr::mutate(latent2 = alpha*priceE + beta*x + Xe + e)
df

#Compute $y2$ by comparing the latent2 values of $k = 1, 2$ for each $i$ to obtain the following result:
df <-
  df %>%
  dplyr::group_by(i) %>%
  dplyr::mutate(y2 = ifelse(latent2 == max(latent2), 1, 0)) %>%
  dplyr::ungroup()
df

# 3. You observe priceE, xk and yk. Write up the log likelihood function of the new choice model and estimate it by ML
loglikelihood_E <-function( temp, df ) {
  lE <- df %>%
    dplyr::group_by(i) %>%
    dplyr::mutate(p2 = exp(temp[1]*priceE+temp[2]*x)/sum(exp(temp[1]*priceE+temp[2]*x))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(y2 == 1)
  lE <- mean(log(lE$p2))
  return(lE) 
}

result <-
  optim( par = theta, fn = loglikelihood_E,
         df = df, control = list(fnscale = - 1))
result

# Q: Can you explain the Omitted Variable Bias you get for alpha_hat compared to the true alpha?
# ANS: we do not include Xe in the optimization above.

# 4. Suppose you also get data for price. Use as instrument.
reg_firstStage<-lm(priceE ~ x + price - 1, df)
summary(reg_firstStage)
#add first stage residuals to dataframe
df <-
  df %>%
  dplyr::mutate(eFS=reg_firstStage$residuals)

# 5. Estimate alpha and beta by ML using the exogenous portion of priceE
loglikelihood_cf <-function( temp, df ) {
  lcf <- df %>%
    dplyr::group_by(i) %>%
    dplyr::mutate(pcf = exp(temp[1]*priceE+temp[2]*x+temp[3]*eFS)/sum(exp(temp[1]*priceE+temp[2]*x)+temp[3]*eFS))  %>%
    dplyr::ungroup() %>%
    dplyr::filter(y2 == 1)
  lcf <- mean(log(lcf$pcf))
  return(lcf) 
}

result <-
  optim( par = c(theta, 0), fn = loglikelihood_cf,
         df = df, control = list(fnscale = - 1))
result

# ANS: my answers are not exactly the same but close. This is likely due to the
# seed being set at different points along the script.
