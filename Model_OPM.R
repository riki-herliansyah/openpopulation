# load the packages needed 
library(jagsUI)
library(mcmcOutput)
library(tidyverse)
library(readxl)

# import the data
group_count <- read_excel("dataset.xlsx", sheet = "group_count")
group_size <- read_excel("dataset.xlsx", sheet = "group_size")

model_size <- "model{
# Priors
  beta0 ~ dnorm(0, 0.001)
  beta1 ~ dnorm(0, 0.001)
  beta2 ~ dnorm(0, 0.001)
  for(i in 1:G){
    size[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta0 + beta1*group[i] + beta2*year[i]
  }
  for(i in 1:2){
    for(t in 1:3){
      size_[i,t] ~ dpois(lambda_[i,t])
      log(lambda_[i,t]) <- beta0 + beta1*(i-1) + beta2*year_[t]
    }
  }
}"
group_size <- as.data.frame(group_size)
data <- list(size = group_size[,4], group = group_size[,3], G = dim(group_size)[1],
             year_ = c(1, 6, 11), year = group_size$Year - 2006)
parameters <- c("beta0", "beta1","beta2",
                "size_[1,1]", "size_[1,2]", "size_[1,3]",
                "size_[2,1]", "size_[2,2]", "size_[2,3]")
# fit the model using jags function
jagsOut <- jags(data, NULL, parameters, textConnection(model_size), DIC=FALSE,
                n.chains=3, n.iter=55000, n.adapt=5000, n.burnin = 5000, parallel=FALSE)
pop.mod <- mcmcOutput(jagsOut)
# trace plots
diagPlot(pop.mod, main ="N-Mixture model parameter")
# posterior summaries
summary.mod <- summary(pop.mod)
summary.mod

# Bayesian Analysis for Open Population Model
model_open <- "model{
  # Priors for group detection
  psi[1] ~ dbeta(2, 2)
  psi[2] ~ dbeta(2, 2)
  psi[3] ~ dbeta(2, 2)
  omega ~ dbeta(2, 2)
  gamma ~ dgamma(2, 0.2)
  lambda_N ~ dgamma(2, 0.2)
  
  # Define the pdf of N_{i1} ~ Pois(lambda)
  for(i in 1:R) {
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
      Ys[i,t] <- Nit[i,t]*phi[t]*size_[t,2] + Nit[i,t]*(1-phi[t])*size_[t,1]
      Y[i,t] ~ dbin(psi[t], Nit[i,t])
       G[i,t] ~ dbin(phi[t], Y[i,t])
    }
  }
  # abundance 
  # # prior for group 
  phi[1] ~ dbeta(1/2,1/2)
  phi[2] ~ dbeta(1/2,1/2)
  phi[3] ~ dbeta(1/2,1/2)
  for (t in 1:T){
    Nt[t] <- sum(Nit[1:R,t])
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
y <- as.matrix(G_count[,2:4])
group_size <- as.data.frame(group_size)
# median of population size per group per period (2007, 2012, 2017)
# AMU = (4,5,6), OMU = (9, 11, 13)
size_ <- matrix(c(4,5,6,9,11,13), nrow = 3)
G <- group_count[group_count$Group_cat == "OMU",3:5]
data <- list(Y = y, R = 12, T = 3, G = G, size_ = size_)
# list of parameter of interest
parameters <- c("N[1]","N[2]","N[3]","omega", "gamma",
                "psi[1]","psi[2]","psi[3]",
                "phi[1]","phi[2]","phi[3]")
# inial state for population at each location (12 river)
inits <- function(){
  list(Ni1 = rpois(12, 100))
}

jagsOut <- jags(data, inits, parameters, textConnection(model_open), DIC=FALSE,
                n.chains=3, n.iter=105000, n.adapt=5000, n.burnin = 5000, parallel=FALSE)
pop.mod <- mcmcOutput(jagsOut)
diagPlot(pop.mod, main ="Open population model parameter")
summary(pop.mod)
