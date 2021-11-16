source("../pareto_linear/pareto_linear.R")
source("../postprocess_functions.R")
library(tidyverse)
library(sf)
library(rstan)
library(sp)
library(maptools)
library(Hmisc) ## for wtd quantiles / means / etc
library(survey) ## for Hadamard matrix
library(reldist)
library(reshape2)
library(xtable)

load("standat.co.RData")
load("standat.il.RData")
load("standat.mo.RData")
load("standat.mt.RData")
load("standat.ny.RData")
load("../data/full_state_datasets.RData")

co.data <- full_state_datasets %>%
    filter(state == "CO", as.numeric(puma) == 821)
co.ests <- estimate_summary(co.data)

il.data <- full_state_datasets %>%
    filter(state == "IL", as.numeric(puma) == 3502)
il.ests <- estimate_summary(il.data)

mo.data <- full_state_datasets %>%
    filter(state == "MO", as.numeric(puma) == 600)
mo.ests <- estimate_summary(mo.data)

mt.data <- full_state_datasets %>%
    filter(state == "MT", as.numeric(puma) == 600)
mt.ests <- estimate_summary(mt.data)

ny.data <- full_state_datasets %>%
    filter(state == "NY", as.numeric(puma) == 3706)
ny.ests <- estimate_summary(ny.data)

## count the number of missing tract-level estimates of each type in each PUMA
## mostly missing tract-level estimates of the 95th percentile of income
## these numbers are reported in the captions of the tables.
## NOTE: estimate_summary treats estimates with missing SEs as missing
##    since these estimates are often deliberately poor due to disclosure limitations
co.ests %>% select(contains("est")) %>% apply(2, function(x){sum(is.na(x))})
il.ests %>% select(contains("est")) %>% apply(2, function(x){sum(is.na(x))})
mo.ests %>% select(contains("est")) %>% apply(2, function(x){sum(is.na(x))})
mt.ests %>% select(contains("est")) %>% apply(2, function(x){sum(is.na(x))})
ny.ests %>% select(contains("est")) %>% apply(2, function(x){sum(is.na(x))})

## for each model/PUMA combination, compute the RMSE of the posterior median estimates
## of each target percentile, treating the direct estimate as the true value,
## as a percentage of the direct estimate
target_quantiles <- c(0.2, 0.4, 0.6, 0.8, 0.95)
  ## quantiles we wish to do inference at
states <- c("CO", "IL", "MO", "MT", "NY")

set.seed(923429)
ncore <- 1  ## can set > 1 when not on windows
results <- list()
for(state in states){
    print(state)
    file.name <- paste("fit.", tolower(state), ".RData", sep = "")
    load(file.name)
    standat <- switch(state,
                      CO = standat.co,
                      IL = standat.il,
                      MO = standat.mo,
                      MT = standat.mt,
                      NY = standat.ny)
    ests <- switch(state,
                   CO = co.ests,
                   IL = il.ests,
                   MO = mo.ests,
                   MT = mt.ests,
                   NY = ny.ests)
    ## if the SE is missing, treat the estimate as missing

    pred <- latentprln_predict(target_quantiles, standat, fit, ncore)
    save(pred, file = paste("pred.", tolower(state), ".RData", sep = ""))

    prln_est <- matrix(0, nrow = standat$ntract, ncol = 1 + length(target_quantiles))
    colnames(prln_est) <- c("Gini", paste("Q", c(20, 40, 60, 80, 95), sep=""))
    alphalist <- list()
    for(i in 1:standat$ntract){
        temp <- pareto.linear(standat$bounds, standat$bin_est[i,], target_quantiles)
        alphalist[[i]] <- temp$alpha
        prln_est[i,] <- c(temp$giniq, temp$quants)
    }
    ## note: alpha = 1 implies the bin is uniform
    ## or if it's the top bin it implies that the bin is a point mass on the bin minimum
    ## the 12th bin is the top bin (which is what we are printing)
    print(sort(simplify2array(alphalist)[12,]))

    estimates <- ests %>% select(GEOID, model, contains("est"))

    ## evaluate model based latent PRLN predictions using held out estimates
    pred_lprln_metric <- eval_lprln_preds_acspuma(pred, estimates) %>%
        mutate(State = toupper(state)) %>%
        select(State, Estimator, Metric, everything())

    ## evaluate PRLN based predictions
    pred_prln_metric <- eval_prln_preds_acspuma(prln_est, estimates) %>%
        mutate(State = state) %>%
        select(State, Estimator, Metric, everything())

    pred_metric <- bind_rows(pred_lprln_metric, pred_prln_metric) %>%
        select(State, Metric, Estimator, everything()) %>%
        mutate(Estimator = factor(Estimator, levels = c("PRLN", "Mean", "Median"))) %>%
        arrange(State, Metric, Estimator)

    results[[state]] <- pred_metric
}

save(results, file = "results.RData")


## Create tables D.1, D.2, D.3, D.4, and D.5
load("results.RData")
for(state in states){
    results[[state]] %>% xtable %>% print
}
