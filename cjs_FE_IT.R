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
nI <- nrow(dat$y)#100

dat$y <- dat$y[1:nI,]
dat$nOcc = ncol(dat$y)
dat$nInd = nrow(dat$y)


## Initial values 
inits <- function() list(phiBeta = runif(dat$nOcc-1, 0, 1),
                         pBeta = runif(dat$nOcc-1, 0, 1),
                         sigmaI = runif(1, 0, 2)
                         )

## Parameters monitored
params <- c("phiBeta", "pBeta", "sigmaI")

## MCMC settings
ni <- 5000
nt <- 4
nb <- 1000
nc <- 4

## Call Stan from R
cjs_FE_IT  <- stan("cjs_FE_IT.stan",
                        data = dat, init = inits, pars = params,
                        chains = nc, iter = ni, warmup = nb, thin = nt,
                        seed = 1,
                        open_progress = FALSE)

## Summarize posteriors
print(cjs_FE_IT, digits = 3)

launch_shinystan(cjs_FE_IT)
