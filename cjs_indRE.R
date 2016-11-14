## 7. Estimation of survival probabilities using capture-recapture data
## 7.5. Models with individual variation
## 7.5.3. Individual random effects

library(rstan)
library(shinystan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in getAndPrepareDataWB.R

setwd("/home/ben/stan/cjs/")
#stan_data <- read_rdump("dat.RData")
load(file="dat.RData")

# if subsetting enc data
nI <- 100

dat$y <- dat$y[1:nI,]
dat$nOcc = ncol(dat$y)
dat$nInd = nrow(dat$y)


## Initial values 
inits <- function() list(mean_phi = runif(1, 0, 1),
                         mean_p = runif(1, 0, 1),
                         sigma = runif(1, 0, 1))

## Parameters monitored
params <- c("mean_phi", "mean_p", "sigma", 'epsilon')

## MCMC settings
ni <- 5000
nt <- 4
nb <- 1000
nc <- 4

## Call Stan from R
cjs_indRE  <- stan("cjs_indRE.stan",
                        data = dat, init = inits, pars = params,
                        chains = nc, iter = ni, warmup = nb, thin = nt,
                        seed = 1,
                        open_progress = FALSE)

## Summarize posteriors
print(cjs_indRE, digits = 3)

launch_shinystan(cjs_indRE)
