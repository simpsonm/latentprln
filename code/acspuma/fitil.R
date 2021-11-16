library(rstan)
library(tidyverse)
options(warn=1)


load("../data/puma_standats.RData")

standat <- puma_standats$IL$y2015

nknot <- length(standat$knots)
knots <- replicate(standat$ntract, standat$knots) %>% t()
nunif <- sapply(1:standat$ntract,
                function(x){max(which(standat$median_est[x] > knots[x,]))})
nalpha <- nknot - nunif

state <- 'il'

standat.il <- list(
    ntract = standat$ntract,
    nbin = standat$nbin,
    nalpha = nalpha,
    nknot = nknot,
    knots = knots,
    bounds = standat$knots,
    bin_est = standat$bin_est,
    bin_se  = standat$bin_se,
    mean_est = c(standat$mean_est),
    mean_se  = c(standat$mean_se),
    median_est = c(standat$median_est),
    median_se  = c(standat$median_se),
    alpha_prior_mean = rep(2, sum(nalpha)),
    alpha_prior_sd = 1,
    pknot_tract_prior_loc = t(replicate(standat$ntract, standat$pknot_prior_loc)),
    pknot_tract_prior_scale = 1/10,
    tract_pops = c(standat$total_est)
)

save(standat.il, file = paste('standat.', state, '.RData', sep = ""))

set.seed(8713482)

fit0 <- stan("../models/prln_manytracts.stan",
             data = standat.il, chains = 1, iter = 1)

fit <- stan(fit = fit0, data = standat.il, chains = 4, cores = 4,
            iter = 8000, warmup = 4000,
            open_progress = FALSE)

file.name <- paste("fit.", state, ".RData", sep = "")
save(fit, file = file.name)
