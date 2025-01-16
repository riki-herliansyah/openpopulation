# load the packages needed 
library(jagsUI)
library(mcmcOutput)
library(tidyverse)
library(readxl)

# import the data
group_count <- read_excel("dataset.xlsx", sheet = "group_count")
group_size <- read_excel("dataset.xlsx", sheet = "group_size")

# Build the model using rjags
model_open <- "model{
  
  # log-Poisson model for group size estimation
  for(i in 1:Gn){
    size[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta0 + beta1*group[i] + beta2*year[i]
  }
  # obtain the expected group size for each group and each year
  for(i in 1:2){
    for(t in 1:3){
      mu[i,t] <- exp(beta0 + beta1*(i-1) + beta2*year_[t])
    }
  }
  # Priors
  beta0 ~ dnorm(0, 0.001)
  beta1 ~ dnorm(0, 0.001)
  beta2 ~ dnorm(0, 0.001)
  
 # Open Population Model
  for(i in 1:R) {
     # Define the pdf of N_{i1} ~ Pois(lambda)
    Ni1[i] ~ dpois(lambda_N)
    # combine Ni1 to matrix Nit
    Nit[i,1] <- Ni1[i]
    
    for(t in 1:(T-1)) {
      # recruitment process
      Git[i,t] ~ dpois(gamma)
      # survival process
      Sit[i,t] ~ dbin(omega, Nit[i,t])
      Nit[i,t+1] <- Sit[i,t] + Git[i,t] 
    }
  }
  for(i in 1:R){
    for(t in 1:T){
      #observation process
      Ys[i,t] <- Nit[i,t]*phi[t]*mu[2,t] + Nit[i,t]*(1-phi[t])*mu[1,t]
      Y[i,t] ~ dbin(p, Nit[i,t])
      G[i,t] ~ dbin(phi[t], Y[i,t])
    }
  }
  # Priors for group detection
  prec ~ dgamma(4,3)
  sigma <- 1/prec
  for(t in 1:1){
    alpha ~ dnorm(0, sigma)
    logit(p) <- alpha
  }
  # prior for omega
  logit(omega) <- tau
  tau ~ dnorm(0, s_tau)
  prec_tau ~ dgamma(4,3)
  s_tau <- 1/prec_tau
  
  # prior for gamma
  gamma ~ dgamma(2, 0.2)
  
  # prior for lambda_N
  lambda_N ~ dgamma(0.1, 0.1)
 
  # prior for OMU detection
  phi[1] ~ dbeta(1/2,1/2)
  phi[2] ~ dbeta(1/2,1/2)
  phi[3] ~ dbeta(1/2,1/2)
  
  # obtain the abundance by summing over Ys
  for (t in 1:T){
    N[t] <- sum(Ys[1:R,t])
  }
}"

# pre-processing data
G_count <- group_count %>%
  group_by(Site) %>%
  summarise(T7 = sum(`2007`),
            T12 = sum(`2012`),
            T17 = sum(`2017`),
            .groups = 'drop') %>%
  as.data.frame()
# observed group counts
y <- as.matrix(G_count[,2:4])
# #site (R) and #time (T)
R = 12; T = 3
# observed group size count
group_size <- as.data.frame(group_size)
size = group_size[,4]
# covariate model Poisson
group = group_size[,3] # group cat
year = group_size$Year - 2006 # year
# number of OMU group detected
G <- group_count[group_count$Group_cat == "OMU",3:5]
# import all data required
data <- list(Y = y, R = R, T = T, G = G, size = size, group = group, year = year, 
             Gn = dim(group_size)[1], year_ = c(1, 6, 11))
# list of parameter of interest
parameters <- c("N[1]","N[2]","N[3]","omega", "gamma",
                "p", "beta1", "beta2", "beta0",
                "mu[1,1]","mu[2,1]",
                "mu[1,2]","mu[2,2]",
                "mu[1,3]","mu[2,3]")
# inital state for population at each location (12 river)
inits <- function(){
  list(Ni1 = rpois(12, 100))
}

jagsOut <- jags(data, inits, parameters, textConnection(model_open), DIC=FALSE,
                n.chains=3, n.iter=105000, n.adapt=5000, n.burnin = 5000, parallel=FALSE)
pop.mod <- mcmcOutput(jagsOut)
diagPlot(pop.mod, main ="Open population model parameter")
summary(pop.mod)
