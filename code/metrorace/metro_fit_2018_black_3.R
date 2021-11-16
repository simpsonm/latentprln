library(tidyverse)
library(rstan)
library(foreach)
library(doParallel)
source("postprocess_functions.R")
##library(shinystan)

nwarm = 2000
niter = nwarm + 2000
nthin = 1
nchain = 1
ncore = 1
delta = 0.995
maxtree = 18
initr = 0.5
nclustercore = 28
nklsim = 500
nhrsim = 1000

chain = 3
race = "black"
year = "y2018"
logfile_name = "metro_2018_black_3.txt"
outfile_name = "metro_kls_black_3.RData"

pknot_prior_scale = 1/10
alpha_prior_mean = 2
alpha_prior_sd = 1

load("../data/metro_race_standats.RData")

niter_per_chain = (niter - nwarm) / nthin

set.seed(374239109)
seeds = list()
for(metro in names(metro_race_standats)){
    seeds[[metro]] = list()
    metro_standat = metro_race_standats[[metro]][[year]][[race]]
    for(tract in 1:metro_standat$ntract){
        seeds[[metro]][[tract]] =
            sum(sample(0:9, 8, TRUE)*10^c(1:8))
    }
}

model_init = stan_model(file = "../models/prln_tract.stan")
model_init_nomedian = stan_model(file = "../models/prln_tract_nomedian.stan")

cl = makeCluster(nclustercore, outfile = logfile_name)
registerDoParallel(cl)

system.time({
metro_kls = foreach(metro = names(metro_race_standats),
                    .packages = c("tidyverse", "rstan")) %dopar% {
    metro_warnings = list()
    metro_standat = metro_race_standats[[metro]][[year]][[race]]
    metro_samples = list()
    for(tract in 1:metro_standat$ntract){
        set.seed(seeds[[metro]][[tract]])
        standat = metro_standat
        standat$bin_est = metro_standat$bin_est[tract,]
        standat$bin_se = metro_standat$bin_se[tract,]
        standat$mean_est = metro_standat$mean_est[tract]
        standat$mean_se = metro_standat$mean_se[tract]
        standat$median_est = metro_standat$median_est[tract]
        standat$median_se = metro_standat$median_se[tract]
        standat$total_est = metro_standat$total_est[tract]
        standat$total_se = metro_standat$total_se[tract]
        standat$pknot_prior_scale = pknot_prior_scale
        standat$alpha_prior_mean = alpha_prior_mean
        standat$alpha_prior_sd = alpha_prior_sd
        cat("\n")
        cat("\n")
        cat("Metro: ")
        cat(metro)
        cat("\n")
        cat("Year: ")
        cat(year)
        cat("\n")
        cat("Race: ")
        cat(race)
        cat("\n")
        cat("Tract: ")
        cat(tract)
        cat("\n")
        cat("Starting")
        cat("\n")
        cat("\n")
        if(is.na(standat$median_se)){
            mymodel = model_init_nomedian
        } else {
            mymodel = model_init
        }
        if(!is.na(standat$mean_se) && standat$total_est >= 100){
            quiet_fit = quietly(sampling)(
                mymodel,
                data = standat,
                chains = nchain, cores = ncore, thin = nthin,
                warmup = nwarm, iter = niter,
                init_r = initr,
                open_progress = FALSE,
                control = list(adapt_delta = delta,
                               max_treedepth = maxtree))
            metro_warnings[[tract]] = quiet_fit$warnings
            cat(quiet_fit$warnings)
            metro_samples[[tract]] =
                rstan::extract(quiet_fit$result, pars = c("pknot", "alpha"),
                               permuted = FALSE)
            rm(quiet_fit)
        } else {
            metro_samples[[tract]] = NA
            metro_warnings[[tract]] = NA
        }
        cat("\n")
        cat("\n")
        cat("Metro: ")
        cat(metro)
        cat("\n")
        cat("Year: ")
        cat(year)
        cat("\n")
        cat("Race: ")
        cat(race)
        cat("\n")
        cat("Tract: ")
        cat(tract)
        cat("\n")
        cat("Finished")
        cat("\n")
        cat("\n")
    }
    weights = metro_standat$total_est %>% as.vector
    weights[c(which(is.na(metro_samples)),
              which(sapply(metro_samples, is.null)))] = 0
    weights = weights / sum(weights)
    knots = c(metro_standat$knots, Inf)
    if(sum(!is.na(metro_samples)) > 1){
        kls = compute_kl_index(nklsim, metro_samples, weights, knots)
        hrs = compute_hr_index(nhrsim, metro_samples, weights, knots)
        out = tibble(metro = metro,
                     year = year,
                     race = race,
                     ntract = standat$ntract,
                     nas = sum(is.na(metro_samples)),
                     chain = chain,
                     iter = 1:niter_per_chain,
                     klest = as.vector(kls$kl_index_ests),
                     klse = as.vector(kls$kl_index_ses),
                     hrest = as.vector(hrs$hr_index_ests),
                     hrse = as.vector(hrs$hr_index_ses))
    } else {
        out = tibble(metro = metro,
                     year = year,
                     race = race,
                     ntract = standat$ntract,
                     nas = sum(is.na(metro_samples)),                     
                     chain = chain,
                     iter = 1:niter_per_chain,
                     klest = NA,
                     klse = NA,
                     hrest = NA,
                     hrse = NA)
    }
    rm(metro_samples)
    list(warnings = metro_warnings, kls = out)
}
})

save(metro_kls, file = outfile_name)

stopCluster(cl)
