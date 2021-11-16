library(sp)
library(maptools)
library(Hmisc) ## for wtd quantiles / means / etc
library(survey) ## for Hadamard matrix
library(reldist)
library(dplyr)

## Construct a stratified (ACS) sample from a given population, then given that sample:
##  1) create tract level estimates of various quantities
##  2) create SEs of those estimates using successive difference replication
##  3) create an artificial PUMS as a subset of the original sample
##
## population = data frame of the population
##    each row is a unit (person/household) of the population
##    several columns:
##      tract:         tract indicator for that unit
##      stratum:       stratum indicator for that unit
##      target:        value of the target variable for that unit
##      weight:        weight associated with their stratum
##      stratum.nobs:  number of observations for their stratum in the full sample
## bounds = bounds of desired bin estimates for each tract
## target.quantiles = target quantiles we will estimate with the model
## sam.frac = proportion of each stratum to sample in the ACS sample
## pums.frac = proprotion of each stratum to subsample for the PUMS
##
## Returns a named list of:
##   The full ACS sample (data frame)
##   The PUMS subsample  (data frame)
##   Various tract level estimates (matrices, or when appropriate vectors)
##     (separate elements for bin, mean, median, gini coefficient, and target quantile estimates)
##   SEs of those tract level estimates (matrices, or when appropriate vectors)
data.create <- function(population, bounds, target.quantiles, pums.frac = 0.385){
  n.tract <- length(unique(population$tract))
  n.strata <- length(unique(population$stratum))
  n.bounds <- length(bounds)
  n.bins <- n.bounds + 1
  n.target.quantiles <- length(target.quantiles)
  q.target <- matrix(0, nrow = n.tract, ncol = n.target.quantiles)
  q.target.se <- matrix(0, nrow = n.tract, ncol = n.target.quantiles)
  median.est <- rep(0, n.tract)
  median.se <- median.est
  gini.est <- rep(0, n.tract)
  gini.se <- median.est
  mean.est <- rep(0, n.tract)
  mean.se <- mean.est
  bin.est <- matrix(0, n.tract, n.bins)
  bin.se <- bin.est
  Had <- 2*hadamard(79) - 1  ## hadamard matrix for computing replicate weights
  ## create the initial ACS sample
  strata.samples <- lapply(1:n.strata, function(i){
    pop.stratum <- filter(population, stratum == i)
    pop.stratum.size <- nrow(pop.stratum)
    n.sam <- pop.stratum$stratum.nobs[1]
    sam.stratum <- pop.stratum[sample(1:pop.stratum.size, n.sam),]
    return(sam.stratum)
  })
  combined.acs <- Reduce(rbind, strata.samples, NULL)
  ## subsample to create the pums
  pums.samples <- lapply(1:n.strata, function(i){
    sam.stratum <- filter(combined.acs, stratum == i)
    sam.stratum.size <- nrow(sam.stratum)
    n.sam <- ceiling(pums.frac * sam.stratum.size)
    pums.stratum <- sam.stratum[sample(1:sam.stratum.size, n.sam),]
    return(pums.stratum)
  })
  combined.pums <- Reduce(rbind, pums.samples, NULL)
  ## create tract-level estimates and SEs
  for(i.tract in 1:n.tract){
    ## pull out the ACS sample for current tract
    sam.tract <- filter(combined.acs, tract == i.tract)
    n.sam <- nrow(sam.tract)
    ## create replicate weights for current tract
    ## (using successive difference replication for variance estimation)
    rep.wts <- matrix(0, nrow = n.sam, ncol = 80)
    RI <- 1; RII <- 2; inc <- 1; cyc <- 1
    for(i in 1:n.sam){
      RI <- RI + inc
      RII <- RI + inc
      if(RII > 80){
        if(cyc < inc){
          cyc <- cyc + 1
        } else {
          cyc <- 1
          inc <- inc + 1
        }
        RI <- inc + cyc - 1
        RII <- RI + inc
      }
      rep.wts[i,] <- (1 + (2^(-3/2)*Had[RI,]) - (2^(-3/2)*Had[RII,]))*sam.tract$weight[i]
    }
    colnames(rep.wts) <- paste("RepWt", 1:80, sep ="")
    sam.tract <- cbind(sam.tract, rep.wts)
    ## construct ACS estimates for the tract
    median.est[i.tract] <- Hmisc::wtd.quantile(sam.tract$target, sam.tract$weight, 0.5)
    gini.est[i.tract] <- gini(sam.tract$target, sam.tract$weight)
    q.target[i.tract,] <- Hmisc::wtd.quantile(sam.tract$target, sam.tract$weight, target.quantiles)
    mean.est[i.tract] <- Hmisc::wtd.mean(sam.tract$target, sam.tract$weight)
    cdf.ests <- sapply(1:n.bounds, function(x){
      sum((sam.tract$target < bounds[x])*sam.tract$weight)/ sum(sam.tract$weight)})
    bin.est[i.tract,] <- diff(c(0, cdf.ests, 1))
    ## construct SEs of ACS estimates for the tract
    #### first construct replicate estimates
    q.target.reps <- sapply(1:80, function(x){
      Hmisc::wtd.quantile(sam.tract$target, rep.wts[,x], target.quantiles)})
    median.reps <- sapply(1:80, function(x){
      Hmisc::wtd.quantile(sam.tract$target, rep.wts[,x], 0.5)})
    gini.reps <- sapply(1:80, function(x){
      gini(sam.tract$target, rep.wts[,x])})
    mean.reps <- sapply(1:80, function(x){
      Hmisc::wtd.mean(sam.tract$target, rep.wts[,x])})
    bin.reps <- sapply(1:80, function(x){
      cdf.ests <- sapply(1:n.bounds, function(y){
        sum((sam.tract$target < bounds[y])*rep.wts[,x])/ sum(rep.wts[,x])})
      diff(c(0, cdf.ests, 1))})
    #### use replicate estimates to compute standard errors
    median.se[i.tract]  <- sqrt(4/80*sum((median.reps  - median.est[i.tract] )^2))
    gini.se[i.tract] <- sqrt(4/80*sum((gini.reps  - gini.est[i.tract] )^2))
    mean.se[i.tract] <- sqrt(4/80*sum((mean.reps - mean.est[i.tract])^2))
    q.target.se[i.tract,] <- sapply(1:n.target.quantiles, function(x){
      sqrt(4/80*sum((q.target[i.tract,x] - q.target.reps[x,])^2))})
    bin.se[i.tract,] <- sapply(1:n.bins, function(x){
      sqrt(4/80*sum((bin.est[i.tract,x] - bin.reps[x,])^2))})
    ## correction for bin estimates that are exactly zero
    idx0 <- which(bin.se[i.tract,] <= 0)
    if(length(idx0) > 0){
      ## using the true population average weight and denominator
      pop.tract <- filter(population, tract == i.tract)
      avg.weight <- mean(pop.tract$weight)
      denom <- nrow(pop.tract)
      pstar <- min(2.3*avg.weight/denom, 0.5)
      bin.se[i.tract, idx0] <- sqrt(pstar*(1 - pstar)*avg.weight/denom)
    }
  }
  out <- list(acs.sam = combined.acs, pums.sam = combined.pums,
              median.est = median.est, q.target = q.target, mean.est = mean.est, bin.est = bin.est,
              median.se = median.se, q.target.se = q.target.se, mean.se = mean.se, bin.se = bin.se,
              gini.est = gini.est, gini.se = gini.se)
  return(out)
}

## Given a sample from the joint posterior of the tract-level SDs of the
## inverted median estimates, fit an inverse gamma model to the draws from
## each tract using maximum likelihood and return the parameter estimates
##
## SDs = n.sim x n.tract matrix MCMC samples of the SD of the inverted median estimate by tract
##
## Returns a list of two elements:
##   dfs:    each tract's df estimate (df = alpha/2 in terms of the usual IG parameterization)
##   scales: each tract's scale estimate (scale = sqrt(beta/alpha) in terms of the usual IG parameterization)
ig.fit <- function(SDs){
  n.tract <- ncol(SDs)
  vars <- SDs^2
  dfs <- rep(0, n.tract)
  scales <- rep(0, n.tract)
  for(tract in 1:n.tract){
    draws <- vars[,tract]
    scale.squared <- 1/mean(1/draws) ## can derive optimal scale parameter analytically
    log.scale.squared <- log(scale.squared)
    scales[tract] <- sqrt(scale.squared)
    logmean <- mean(log(draws))
    n <- length(draws)
    ## find optimal df parameter using line search to set d loglikelihood / d df = 0
    ## plugging in estimate of scale parameter
    uniout <- uniroot(function(x){n*log(x/2)/2 + n*log.scale.squared/2 - n*digamma(x/2)/2 - n*logmean/2},
                      c(10^-10, 10^10))
    dfs[tract] <- uniout$root
    ## very occasionally in simulation studies, we get an estimated df < 0, which
    ## implies the standard deviation approximation does not exist.
    ## This is a hack to prevent it from occurring and breaking the simulation
    ## Note: in the both the simulation study reported in the paper and
    ##   in the models for the ACS data it never happened.
    if(dfs[tract] <= 2)
      dfs[tract] <- 2.1
  }
  return(list(dfs = dfs, scales = scales))
}

## Given draws from the posterior distribution of the tract-level parameters,
## construct the posterior predictive distribution of a variety of tract-level estimates.
##
## pars =              list of MCMC samples of parameters
##                      minimally contains named elements mu_tract and sigma_tract
##                      (both n.sim x n.tract x n.mix arrays)
##                      mixture models (when n.mix > 1) also contain pmix_tract
##                      (n.xim x n.tract x n.mix array)
## tract.size =        vector of tract populations
## bounds =            vector of bin estimate bounds
## target.quantiles =  vector of target quantiles we are interested in
## pop.summary =       data frame of the true values of each estimate in each tract in the population
##
## Returns: a data frame containing psoterior draws for:
##          1) the difference between the posterior median of each tract-level estimand and
##             the true value of the estimand in the population
##          2) whether the 95% credible intervals of each tract-level estimand covered
##             the true value of the estimand in the population
##          quantile estimate (including median), mean estimate, and gini coefficient estimate.
stat.mix.sim <- function(pars, tract.size, bounds, target.quantiles, pop.summary){
  mu <- pars$mu_tract
  sigma <- pars$sigma_tract
  n.tract <- length(tract.size)
  n.mix <- dim(mu)[3]
  ## check to see if parameters came from a mixture model
  if(is.na(n.mix)){
    n.mix <- 1
  } else {
    pmix <- pars$pmix_tract
  }
  n.sim <- dim(mu)[1]
  ## For each tract:
  ##   1) generate the joint posterior predictive distribution of the target variable
  ##      for each unit (person/household) in that tract
  ##   2) Compute the various estimates (bin, quantile, mean, gini) for each
  ##      draw from that joint posterior predictive
  ##   3) Summarize the distribution of those estimates with 95% credival intervals and posterior medians
  ##   4) Return the difference between the true values and the posterior medians,
  ##      and also return whether the intervals covered the true values
  stat.dfs <- lapply(1:n.tract, function(i.tract){
    ## generate the joint posterior predictive
    if(n.mix > 1){
      pred.pop <- sapply(1:n.sim, function(i.sim){
        ids <- sample(1:n.mix, tract.size[i.tract], TRUE, pmix[i.sim, i.tract, ])
        mus <- mu[cbind(i.sim, i.tract, ids)]
        sigmas <- sigma[cbind(i.sim, i.tract, ids)]
        exp(rnorm(tract.size[i.tract], mus, sigmas))
      }) %>% t()
    } else {
      mus <- mu[,i.tract]
      sigmas <- sigma[,i.tract]
      pred.pop <- replicate(tract.size[i.tract], exp(mus + sigmas*rnorm(n.sim)))
    }
    ## compute estimates for each draw from the posterior predictive
    stat.out <- apply(pred.pop, 1, function(x){
      cdfs <- c(0, sapply(bounds, function(y){mean(x<y)}), 1)
      out <- c(quantile(x, probs = target.quantiles),
               mean(x), diff(cdfs), gini(x))
      names(out) <- c(paste("quant", target.quantiles, sep=".."), "mean",
                      paste("bin", c("below", bounds), c(bounds, "above"), sep = ".."),
                      "gini")
      out
    })
    ## summarise posterior predictive of estimates
    stat.median <- apply(stat.out, 1, median)
    stat.lower <- apply(stat.out, 1, quantile, probs = 0.025)
    stat.upper <- apply(stat.out, 1, quantile, probs = 0.975)
    stat.pop <- filter(pop.summary, tract == i.tract) %>% select(-tract) %>% unlist()
    stat.diff <- stat.median - stat.pop
    stat.cover <- (stat.lower < stat.pop & stat.pop < stat.upper)*1
    stat.combined <- c(stat.diff, stat.cover)
    names(stat.combined) <- c(paste(names(stat.diff), "diff", sep = "..."),
                              paste(names(stat.cover), "cover", sep = "..."))
    data.frame(tract = i.tract, name = names(stat.combined), value = stat.combined)
  })
  out <- Reduce(rbind, stat.dfs)
  return(out)
}

## Given the values of the target variable in a tract in the population,
## summarize it computing bin, quantile, mean, and gini estimates
##
## x =                vector target variable for a given tract in the population
## target quantiles = vector of target quantiles
## bounds =           vector of bounds
##
## Returns: a vector values of bin, quantile, mean, and gini estimates
pop.summarise <- function(x, target.quantiles, bounds){
  q.out <- quantile(x, probs = target.quantiles)
  cdfs <- c(0, sapply(bounds, function(y){mean(x<y)}), 1)
  bin.out <- diff(cdfs)
  mean.out <- mean(x)
  gini.out <- gini(x)
  out <- data.frame(t(c(q.out, mean.out, bin.out, gini.out)))
  names(out) <- c(paste("quant", target.quantiles, sep=".."), "mean",
                  paste("bin", c("below", bounds), c(bounds, "above"), sep = ".."),
                  "gini")
  return(out)
}
