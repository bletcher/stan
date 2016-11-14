## 7. Estimation of survival probabilities using capture-recapture data
## 7.5. Models with individual variation
## 7.5.3. Individual random effects


# need to add
# cholesky-based mvn
# ind random effects

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

dat$g <- 2
dat$riverOBear <- dat$riverOBear[1:nI]
 
dat$R <- structure(c(5, 0, 0, 1), .Dim = c(2L, 2L))
dat$df <- 3

## Initial values 
inits <- lapply(1:nc, function(i) {
  list(Omega = matrix(c(1, 0, 0, 1), ncol = 2))})


## Parameters monitored
params <- c("eta_phi", "eta_p","eta_i", "sigma", "sigmaP", "sigmaI", "mean_phi", "mean_p", "mean_i")

## MCMC settings
ni <- 2000
nt <- 4
nb <- 1000
nc <- 4

## Call Stan from R
cjs_FE_g_RE_T_cholesky  <- stan("cjs_FE_g_RE_T_cholesky.stan",
                        data = dat, init = inits, pars = params,
                        chains = nc, iter = ni, warmup = nb, thin = nt,
                        seed = 1,
                        open_progress = FALSE)

## Summarize posteriors
print(cjs_FE_g_RE_T_cholesky, digits = 3)

launch_shinystan(cjs_FE_g_RE_T_cholesky)
