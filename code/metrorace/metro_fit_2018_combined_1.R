library(tidyverse)
library(rstan)
library(foreach)
library(doParallel)
source("postprocess_functions.R")
##library(shinystan)

load("../data/metro_standats.RData")

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

chain = 1
race = "combined"
year = "y2018"
logfile_name = "metro_2018_combined_1.txt"
outfile_name = "metro_kls_combined_1.RData"

pknot_prior_scale = 1/10
alpha_prior_mean = 2
alpha_prior_sd = 1

niter_per_chain = (niter - nwarm) / nthin

set.seed(892394010)
seeds = list()
init_rs = list()
for(metro in names(metro_standats)){
    seeds[[metro]] = list()
    init_rs[[metro]] = list()
    metro_standat = metro_standats[[metro]][[year]]
    for(tract in 1:metro_standat$ntract){
        seeds[[metro]][[tract]] =
            sum(sample(0:9, 8, TRUE)*10^c(1:8))
        init_rs[[metro]][[tract]] = initr
    }
}

model_init =
    stan_model(file = "../models/prln_tract_full.stan")
model_init_nomedian =
    stan_model(file = "../models/prln_tract_full_nomedian.stan")
model_init_noshare =
    stan_model(file = "../models/prln_tract_full_noshare.stan")
model_init_noquant =
    stan_model(file = "../models/prln_tract_gini.stan")
model_init_noquant_nomedian =
    stan_model(file = "../models/prln_tract_gini_nomedian.stan")


cl = makeCluster(nclustercore, outfile = logfile_name)
registerDoParallel(cl)

system.time({
metro_kls = foreach(metro = names(metro_standats),
                    .packages = c("tidyverse", "rstan")) %dopar% {
    metro_warnings = list()
    metro_standat = metro_standats[[metro]][[year]]
    metro_samples = list()
    for(tract in 1:metro_standat$ntract){
        set.seed(seeds[[metro]][[tract]])
        standat = metro_standat
        standat$bin_est = metro_standat$bin_est[tract,]
        standat$bin_se = metro_standat$bin_se[tract,]
        standat$quant_est = metro_standat$quant_est[tract,]
        standat$quant_se = metro_standat$quant_se[tract,]
        standat$share_est = metro_standat$share_est[tract,]
        standat$share_se = metro_standat$share_se[tract,]
        standat$share_lb = metro_standat$share_lb[tract,]
        standat$share_ub = metro_standat$share_ub[tract,]
        standat$mean_est = metro_standat$mean_est[tract]
        standat$mean_se = metro_standat$mean_se[tract]
        standat$median_est = metro_standat$median_est[tract]
        standat$median_se = metro_standat$median_se[tract]
        standat$gini_est = metro_standat$gini_est[tract]
        standat$gini_se = metro_standat$gini_se[tract]
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
        cat("Tract: ")
        cat(tract)
        cat("\n")
        cat("Starting")
        cat("\n")
        cat("\n")
        ## catch NAs in quantiles and fit the model without those estimates
        ##  as well as any corresponding share estimates
        quantNAidx = unique(c(which(is.na(standat$quant_est)),
                              which(is.na(standat$quant_se))))
        if(length(quantNAidx) > 0){
            oldnquant = standat$nquant
            standat$nquant = standat$nquant - length(quantNAidx)
            if(standat$nquant > 0){
                standat$quant_est = as.array(standat$quant_est[-quantNAidx])
                standat$quant_se = as.array(standat$quant_se[-quantNAidx])
                standat$quantiles = as.array(standat$quantiles[-quantNAidx])
                ## also remove share estimates corresponding to any
                ## missing quantile estimates
                ## first 4 quantile estimates are quintile estimates,
                ##  each missing quintile estimate removes the
                ##  income share estimate just below and just above it
                quintNAidx = quantNAidx[quantNAidx < oldnquant]
                ## the last quantile estimate is the 95%ile estimate
                ##  if it is missing, it only removes the last share est
                ninefiveNAidx = quantNAidx[quantNAidx == oldnquant]
                shareNAidx = unique(c(quintNAidx,
                                      quintNAidx + 1,
                                      ninefiveNAidx + 1))
                standat$nshare = standat$nshare - length(shareNAidx)
                standat$share_est = as.array(standat$share_est[-shareNAidx])
                standat$share_se = as.array(standat$share_se[-shareNAidx])
                standat$share_lb = as.array(standat$share_lb[-shareNAidx])
                standat$share_ub = as.array(standat$share_ub[-shareNAidx])
            } else {
                standat$nshare = 0
            }
        }
        ## catch any additional NAs in share_lb/ub/se
        ##  and fit the model without those estimates
        if(standat$nshare > 0){
            shareNAidx = unique(c(which(is.na(standat$share_est)),
                                  which(is.na(standat$share_lb)),
                                  which(is.na(standat$share_ub)),
                                  which(is.na(standat$share_se))))
            if(length(shareNAidx) > 0){
                standat$nshare = standat$nshare - length(shareNAidx)
                if(standat$nshare > 0){
                    standat$share_est = as.array(standat$share_est[-shareNAidx])
                    standat$share_se = as.array(standat$share_se[-shareNAidx])
                    standat$share_lb = as.array(standat$share_lb[-shareNAidx])
                    standat$share_ub = as.array(standat$share_ub[-shareNAidx])
                } else {
                    standat$nquant = 0
                }
            }
        }
        if(standat$nquant > 0){
            if(is.na(standat$median_se)){
                mymodel = model_init_nomedian
            } else {
                mymodel = model_init
            }
        } else {
            if(!is.na(standat$median_se)){
                mymodel = model_init_noquant
            } else {
                mymodel = model_init_noquant_nomedian
            }
        }
        if(!is.na(standat$mean_se) && standat$total_est >= 100){
            quiet_fit = quietly(sampling)(
                mymodel,
                data = standat,
                chains = nchain, cores = ncore, thin = nthin,
                warmup = nwarm, iter = niter,
                init_r = init_rs[[metro]][[tract]],
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
