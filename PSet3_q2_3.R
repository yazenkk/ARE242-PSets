# ------------------------------------------------------------------------------
# ARE 242 Problem Set 3 Questions 2+3
# Spring 2025
# Author: Yazen Kashlan
# Date: 5/14/2025
# Note: this code is adapted from Villas-Boas lecture notes for ARE 242

# ------------------------------------------------------------------------------
options(scipen = 999)  # turn off scientific notation globally

# Setup
#install.packages("ivreg")
library(dplyr)
library(haven) 
library(ggplot2)
library(AER)
library(stargazer)
library(ivreg)


# ------------------------------------------------------------------------------
# Import data
folder <- "/Users/yazenkashlan/Documents/GitHub/personal-code/ARE242/"
folder_data <- "/Users/yazenkashlan/Documents/GitHub/personal-code/ARE242/data/"
folder_output <- "/Users/yazenkashlan/Documents/GitHub/personal-code/ARE242/output/"
mydata <- read_dta(glue("{folder_data}pset3data.dta"))


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 2.1. Estimate a Logit Discrete Choice Demand differentiated for old and new
# type of products

# ------------------------------------------------------------------------------
# a) Create the lnsj 􀀀 lns0 left hand side variable and label it \log share di erence." Replace the values of
# the variable you generated with NA if it is equal to Inf. Drop the rows for which the variable is NA
mydata <- mydata %>%
  mutate(lnsjs0 = log(sj) - log(s0))
#deal with infinite drop those rows them from the dataset
mydata$lnsjs0[is.infinite(mydata$lnsjs0)] <- NA
#drop those from the data
mydata <- mydata %>% filter(!is.na(mydata$lnsjs0))
#label it
attr(mydata$lnsjs0, "label") <- "log share difference"

# ------------------------------------------------------------------------------
# b) Estimate the above logit model by OLS 
# lnsj - lns0 = beta_0 + alpha*pricejt + beta_1*lncapacityjt + beta_3*newjt +
# beta_4*NonUSjt + beta_5*use1jt + yearFE + ksi_jt. 

# we need year FEs
mydata$feyear <- model.matrix(~ factor(mydata$year) - 1)  # Removes the intercept term
# Price: use p98th
q2_1_b <- lm(lnsjs0 ~ p98th + lncapacity + New + nonUS + use1 + feyear, 
             data = mydata)
summary(q2_1_b)

# What is the estimated marginal utility of price?
# ANS:
q2_1_b$coefficients["p98th"]


# ------------------------------------------------------------------------------
# c) Estimate the same logit model by IV using the Hausman instrument for price 
# (the average price in other markets given by priceIV1). 
# What is the price marginal utility estimate using these instruments ? 
# Did the change make sense given the endogeneity bias concern of the OLS regression?

# Run Hausmann regression
q2_1_iv <- ivreg(
  lnsjs0 ~ p98th    + lncapacity + New + nonUS + use1 + feyear |
           priceIV1 + lncapacity + New + nonUS + use1 + feyear,
  data = mydata
)
summary(q2_1_iv) # View results


# ------------------------------------------------------------------------------
# d) Estimate the same logit model by IV using the BLP instruments for price 
# (the characteristics of other products int he market)

q2_1_blp <- ivreg(
  lnsjs0 ~ p98th    
            + lncapacity + New + nonUS + use1 + feyear |
           frontierTypeCapacity + NofYearsInMktOfTypeCapacity +
           NofFirmsInTypeCapacity + NofProductsInType +
           NofFirmsInType + NofYearsInMktOfCapacity +
            + lncapacity + New + nonUS + use1 + feyear,
  data = mydata
)
summary(q2_1_blp) # View results
summary(q2_1_blp, diagnostics = TRUE)

# Check first stage
first_stage <- lm(
  p98th ~ frontierTypeCapacity + NofYearsInMktOfTypeCapacity +
    NofFirmsInTypeCapacity + NofProductsInType +
    NofFirmsInType + NofYearsInMktOfCapacity +
    lncapacity + New + nonUS + use1 + feyear,
  data = mydata
)
summary(first_stage)$fstatistic # 34>10 no weak IV


# ------------------------------------------------------------------------------
#e) Using Hausman and BLP instruments, estimate the same logit model and 
# recover the unobserved product quality ksi_jt.

reg_iv_combined <- ivreg(lnsjs0 ~ p98th + lncapacity + New + nonUS + use1 + feyear | 
                           priceIV1 + frontierTypeCapacity + NofYearsInMktOfTypeCapacity + 
                           NofFirmsInTypeCapacity + NofProductsInType + 
                           NofFirmsInType + NofYearsInMktOfCapacity + 
                           lncapacity + New + nonUS + use1 + feyear, data = mydata)

mydata$Xis_jt<-residuals(reg_iv_combined)
#label it
attr(mydata$Xis_jt, "label") <- "IV residuals"

#average the Xi by year type first, then plot
#Calculate mean of Xis_jt by year and type
#Collapse data by calculating the mean of Xis_jt for each combination of year and type
collapsed_data <- mydata %>%
  group_by(year, type) %>%
  summarise(mean_Xis_jt = mean(Xis_jt, na.rm = TRUE))
attr(collapsed_data$mean_Xis_jt, "label") <- "Average Quality"
#View the collapsed data
head(collapsed_data)


# ------------------------------------------------------------------------------
# f) In the same graph, please plot the estimated unobserved average 
# (averaged to year type) product quality $\xi_{jt}$ separately for the 
# new and old product types over time. 

delta_plot <- ggplot(collapsed_data, aes(x = year, y = mean_Xis_jt, color = type)) +
  geom_line(linewidth = 1) +     # updated here
  geom_point() +
  labs(
    title = "Estimated Unobserved Quality (Ksi)",
    x = "Year",
    y = "Estimated Unobserved Quality (Ksi)",
    color = "Type"
  ) +
  theme_minimal()
ggsave(glue("{folder_output}PSet3_q2_2.png"), delta_plot)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 2.2. Assuming Cournot Model estimate margins and marginal costs


# a) Load the data
mydata22 <- read_dta(glue("{folder_data}datapset3cournot.dta"))
head(mydata22)


# Extract the needed demand parameters from the demand estimates of IV combined Hausman and BLP instruments
reg_iv_combined <- ivreg(lnsjs0 ~ p98th + lncapacity + New + nonUS + use1 + feyear | 
                           priceIV1 + frontierTypeCapacity + NofYearsInMktOfTypeCapacity + 
                           NofFirmsInTypeCapacity + NofProductsInType + 
                           NofFirmsInType + NofYearsInMktOfCapacity + 
                           lncapacity + New + nonUS + use1 + feyear, data = mydata)

mydata$Xis_jt<-residuals(reg_iv_combined)

# a) cont. Extract the needed demand parameters from the demand estimates of IV combined Hausman and BLP instruments
coefftemp<-coefficients(reg_iv_combined) 
alpha1<-coefftemp[2]
alpha1
#constant
beta0<-coefftemp[1]
beta0
#coeff on lncapacity
beta1<-coefftemp[3]
beta1
#on new
beta3<-coefftemp[4]
beta3
#the estimated year FE
yearFE<-coefftemp[7:21]
# we do not have non-us nor use1 level transactions data so we ignore those here

# b) initialize the matrices to store the MC and prices given the demand and supply model to zero
MCo <- matrix(0, nrow = 16, ncol = 1)
MCn <- matrix(0, nrow = 16, ncol = 1)
Po <- matrix(0, nrow = 16, ncol = 1)
Pn <- matrix(0, nrow = 16, ncol = 1)

#the last year t=16 does not have yearFE estimate, it was the one dropped
# Loop over the time periods from 1 to 15
#the last year 16 has no estimated yearFE estimate it was the one dropped

# Loop over the time periods from 1 to 15 
for (t in 1:15) {
  # Basic data
  m <- mydata22[t, 1]  # Market size
  Qo <- mydata22[t, 2] # Quantity of old products
  Qn <- mydata22[t, 3] # Quantity of new products
  Q0 <- m - Qo - Qn  # Quantity of other products
  sn<-Qn/m
  so<-Qo/m
  s0<-Q0/m
  
  #$\frac{\partial p_o}{\partial Q_o}=\frac{1}{\alpha}[\frac{Q_0+Q_o}{Q_o Q0}]$.
  #$\frac{\partial p_n}{\partial Q_n}=\frac{1}{\alpha}[\frac{Q_0+Q_n}{Q_n Q0}]$.
  # Price derivatives (Cournot model)
  dPoQo <- (Qo + Q0) / (alpha1 * Qo * Q0)
  dPnQn <- (Qn + Q0) / (alpha1 * Qn * Q0)
  dPnQo <- 1 / (alpha1 * Q0)
  dPoQn <- 1 / (alpha1 * Q0)
  
  
  #prices given demand estimated parameters, Q and covariates of logit
  #$p_o=\frac{ln \frac{Q_o}{M}-ln[1-\frac{Q_n}{M}-\frac{Q_o}{M}]-\beta_0- \beta_1 lncapacity_{ot}- yearFE - \xi_{ot}}{\alpha}$
  #  and 
  #$p_n=\frac{ln \frac{Q_n}{M}-ln[1-\frac{Q_n}{M}-\frac{Q_o}{M}]-\beta_0- \beta_1 lncapacity_{nt}-\beta_3 - \beta_4 NonUS_{nt}- \beta_5 fresh_{nt}- yearFE - \xi_{nt}}{\alpha}$
  
  Pn[t]<-(log(sn)-log(s0)-beta0-beta1*mydata22$lncap_n[t]-beta3-mydata22$feyear[t]*yearFE[t])/alpha1
  Po[t]<-(log(so)-log(s0)-beta0-beta1*mydata22$lncap_o[t]-mydata22$feyear[t]*yearFE[t])/alpha1
  
  # Marginal costs
  MCo[t] <- Po[t] + dPoQo * Qo
  MCn[t] <- Pn[t] + dPnQn * Qn
}

#For year t=16 no yearFE

#=================  
#for last t=16
# Basic data
m <- mydata22[16, 1]  # Market size
Qo <- mydata22[16, 2] # Quantity of old products
Qn <- mydata22[16, 3] # Quantity of new products
Q0 <- m - Qo - Qn  # Quantity of other products
sn<-Qn/m
so<-Qo/m
s0<-Q0/m

#$\frac{\partial p_n}{\partial Q_n}=\frac{1}{\alpha}[\frac{Q_0+Q_n}{Q_n Q0}]$.
# Price derivatives (Cournot model)
dPoQo <- (Qo + Q0) / (alpha1 * Qo * Q0)
dPnQn <- (Qn + Q0) / (alpha1 * Qn * Q0)
dPnQo <- 1 / (alpha1 * Q0)
dPoQn <- 1 / (alpha1 * Q0)

#prices given demand estimated parameters, Q and covariates of logit
Pn[16]<-(log(sn)-log(s0)-beta0-beta1*mydata22$lncap_n[16]-beta3)/alpha1
Po[16]<-(log(so)-log(s0)-beta0-beta1*mydata22$lncap_o[16])/alpha1 

# Marginal costs
MCo[16] <- Po[16] + dPoQo * Qo
MCn[16] <- Pn[16] + dPnQn * Qn


# d) In the same graph, plot the marginal costs of the old and new product 
# producing firms over the periods.
mydata22$MCo<-MCo
mydata22$MCn<-MCn

#combined graphs
gCombined <- ggplot(data = mydata22, aes(x = year)) +
  geom_line(aes(y = as.numeric(MCo), color = "MCo")) +
  geom_line(aes(y = as.numeric(MCn), color = "MCn")) +
  labs(y = "Marginal Costs", x = "Year", color = "Legend") +
  theme_minimal()
gCombined
ggsave(glue("{folder_output}PSet3_q2_2d.png"), gCombined)

# Please comment on the evolution of marginal costs over the years for the new and old product.
# ANS: marginal costs for the new products fall more sharply


# e) Plot the evolution of profits for old and new product firms over the years 
# in the same graph. Which profits have been improving over the years?
# construct profits = (p-mc)*Q
mydata22$Pn<-Pn
mydata22$Po<-Po
# convert to numeric
mydata22 <- mydata22 %>%
  mutate(
    Po = as.numeric(as.character(Po)),
    MCo = as.numeric(as.character(MCo)),
    Pn = as.numeric(as.character(Pn)),
    MCn = as.numeric(as.character(MCn))
  )

mydata22 <- mydata22 %>%
  mutate(pi_o = (Po - MCo) * qo) %>%
  mutate(pi_n = (Pn - MCn) * qn)

# plot profits over time
#combined graphs
gCombined_pi <- ggplot(data = mydata22, aes(x = year)) +
  geom_line(aes(y = as.numeric(pi_o), color = "Old product")) +
  geom_line(aes(y = as.numeric(pi_n), color = "New product")) +
  labs(y = "Profits", x = "Year", color = "Legend") +
  theme_minimal()
gCombined_pi
ggsave(glue("{folder_output}PSet3_q2_2e.png"), gCombined_pi)

# f) save for later use
save(mydata22, file=glue("{folder_data}mydata22MC.RData"))


# Question 3
# a) Your task is to find what new post merger marginal costs would need to be 
# such that the prices do not change.

# Compute the marginal costs under the new market structure
#Initialize the matrices
MCos <- matrix(0, nrow = 16, ncol = 1)
MCns <- matrix(0, nrow = 16, ncol = 1)

# Loop over the time periods from 1 to 16 
for (t in 1:16) {
  # Basic data
  m <- mydata22[t, 1]  # Market size
  Qo <- mydata22[t, 2] # Quantity of old products
  Qn <- mydata22[t, 3] # Quantity of new products
  Q0 <- m - Qo - Qn  # Quantity of other products
  sn<-Qn/m
  so<-Qo/m
  s0<-Q0/m
  
  # Price derivatives (Cournot model)
  dPoQo <- (Qo + Q0) / (alpha1 * Qo * Q0)
  dPnQn <- (Qn + Q0) / (alpha1 * Qn * Q0)
  dPnQo <- 1 / (alpha1 * Q0)
  dPoQn <- 1 / (alpha1 * Q0)
  
  #prices are unchanged and so are quantities, therefore
  # $p_o(Q_o^M,Q^M_n) +Q^M_o \frac{\partial p_o}{\partial Q_o}+ Q^M_n \frac{\partial p_n}{\partial Q_o}-mc_o=0$
  # $p_n(Q_o^M,Q^M_n)+Q^M_n \frac{\partial p_n}{\partial Q_n}+Q^M_n \frac{\partial p_o}{\partial Q_n}-mc_n=0$
  # where \frac{\partial p_n}{\partial Q_o} = \frac{\partial p_o}{\partial Q_n} = \frac{1}{alpha Q_o}
  
  # Marginal costs
  MCos[t] <- Po[t] - dPoQo * Qo - dPoQn * Qn
  MCns[t] <- Pn[t] - dPnQn * Qn - dPnQo * Qo
}


# b) Plot the estimated Needed Cost Savings over the years for both types of 
# products on the same graph.
# Hint, use functions and equations coded already from above and add the 
# new parts due to the merger.
CostSavings_o<-as.numeric(MCos)-as.numeric(MCo)
CostSavings_n<-as.numeric(MCns)-as.numeric(MCn)

#combined graphs
gCombinedSavings <- ggplot(data = mydata22, aes(x = year)) +
  geom_line(aes(y = as.numeric(CostSavings_o), color = "MC Savings Old")) +
  geom_line(aes(y = as.numeric(CostSavings_n), color = "MC Savings New")) +
  labs(y = "Marginal Costs Savings such that prices remain unchanged", x = "Year", color = "Legend") +
  theme_minimal()
gCombinedSavings
ggsave(glue("{folder_output}PSet3_q3_b.png"), gCombinedSavings)

