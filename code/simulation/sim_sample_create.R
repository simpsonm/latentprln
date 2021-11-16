## create 1000 datasets by sampling the synthetic population
##
## first create a stratified random sample from the population
##  then create tract-level estimates using the sample & sample weights
##  then subsample the original sample to create the synthetic PUMS

library(rstan)
source("sim_fun.R")
load("population.RData")

RNGkind(sample.kind = "Rounding") ## needed for reproducibility in R >= 3.6

nsam <- 1000

## quantiles we want to do inference on within tracts
target.quantiles <- seq(0.05, 0.95, 0.05)

bounds.1 <- c(5, 10, 15, 20, 25, 35, 50, 75, 100, 150)
bounds.2 <- c(10, 15, 25, 35, 50, 75, 100, 150, 200)
bounds <- c(0, sort(unique(c(bounds.1, bounds.2)))*1000)
nbin <- length(bounds)
nknot <- nbin

pums.frac <- 0.385 ## fraction of the full sample to include in the PUMS subsample
ntract <- length(unique(population$tract))

knots <- replicate(ntract, bounds) %>% t()

us.percentages <- c(3.4, 3.8, 5.3, 5.3, 5.3, 10.1, 13.4, 17.8, 12.1, 13.1, 5.1, 5.3)
pknot.tract.prior.loc <- matrix(us.percentages/100, ntract, nknot, byrow = TRUE)

## basic list of data to be fed into rstan. other data will be added to this.
## this is needed to do construct the approximate SE of the inverted quantile estimates
standat.base <- list(
    ntract = ntract,
    nbin = nbin,
    bounds = bounds,
    nknot = nknot,
    knots = knots,
    alpha_prior_sd = 1,
    pknot_tract_prior_loc = pknot.tract.prior.loc,
    pknot_tract_prior_scale = 1/10
)


samples <- list()
standats <- list()

set.seed(71249234)
for(iter in 1:nsam){
    ## create the synthetic PUMS and tract-level estimates
    ## (see sim_fun.R for description of data.create()
    samples[[iter]] <- data.create(population, bounds[-1], target.quantiles, pums.frac)
    sam <- samples[[iter]]
    standat <- standat.base
    standat$nobs <- length(sam$pums.sam$target)
    standat$z_obs <- sam$pums.sam$target
    standat$bin_est <- sam$bin.est
    standat$bin_se <- sam$bin.se
    standat$mean_est <- sam$mean.est
    standat$mean_se <- sam$mean.se
    standat$median_est <- sam$median.est
    standat$median_se <- sam$median.se
    standat$weight <- sam$pums.sam$weight
    standat$tract_pops <- table(population$tract)
    nunif <- sapply(1:ntract,
                    function(x){max(which(standat$median_est[x] > knots[x,]))})
    nalpha <- nknot - nunif
    standat$nalpha <- nalpha
    standat$alpha_prior_mean = rep(2, sum(nalpha))
    standats[[iter]] <- standat
}

save(samples, file = "samples.RData")
save(standats, file = "standats.RData")

