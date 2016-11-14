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

## Initial values 
inits <- function() list(mean_phi = runif(1, 0, 1),
                         mean_p = runif(1, 0, 1))#,
#                         sigma = runif(1, 0, 2))

## Parameters monitored
params <- c("mean_phi", "mean_p")#, "sigma2","chi")

## MCMC settings
ni <- 5000
nt <- 4
nb <- 1000
nc <- 4

## Call Stan from R
cjs_first  <- stan("cjs_first.stan",
                        data = dat, init = inits, pars = params,
                        chains = nc, iter = ni, warmup = nb, thin = nt,
                        seed = 1,
                        open_progress = FALSE)

## Summarize posteriors
print(cjs_first, digits = 3)

launch_shinystan(cjs_first)
