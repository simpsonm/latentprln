## R functions for running the Pareto-linear procedure.
## Translated from Pascal code avaialble in prln04.pas available at http://www.unc.edu/~nielsen/data/data.htm

## test on lower.bounds and bin.ests below.
## should yield:
##  ginir = 33.41, ginirq = 33.35
## lower.bounds <- c(0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
##                   10000, 12000, 15000, 25000, 50000)
## bin.ests <- c(496, 730, 824, 1221, 1499, 1704, 1712, 1888, 2149,
##               2122, 3808, 3713, 3084, 642, 133)

## lb = lower bound of category
## ub = upper bound of category
## plower = P(X > lb) (or counts)
## pupper = P(X > ub) (or counts)
estimate.alpha <- function(lb, ub, n.low, n.up){
  log(n.low/n.up)/log(ub/lb)
}

## lower.bounds = lower bounds of each income category
##  (assume first category's lower bound is 0)
## bin.ests = estimates of number (or percentage) of people
##  in each income category
pareto.linear <- function(lower.bounds, bin.ests, target.quantiles){
  n.bin <- length(bin.ests)
  rev.cum.counts <- rev(cumsum(rev(bin.ests)))
  count.total <- rev.cum.counts[1]
  cdf.ests <- c(0, cumsum(bin.ests))/count.total
  income.shares <- rep(0, n.bin)
  cum.income.shares <- rep(0, n.bin + 1)
  alphas <- rep(1, n.bin)
  median.idx <- min(which(cdf.ests > 0.5)) - 1 ## because we added zero
  n.quant <- 100
  quant.idx <- 0
  quant.interval <- 1/n.quant
  quant.cum.temp <- rep(0, n.quant)
  for(i in 1:n.bin){
    if( bin.ests[i] == 0 ) {
      ## if the category has no observations, easy
      income.shares[i] <- 0
      cum.income.shares[i+1] <- cum.income.shares[i] + income.shares[i]
    } else if( (i <= median.idx | bin.ests[i + 1] == 0) & i < n.bin ) {
      ## if it has observations and...
      ## if the category is completely < median, or
      ## if the bound above me has no observations,
      ## set income share to midpoint * percentage in bin
      income.shares[i] <- bin.ests[i]*(lower.bounds[i+1] + lower.bounds[i])/2
      cum.income.shares[i+1] <- cum.income.shares[i] + income.shares[i]
      alphas[i] <- 1
      ## compute quantiles in this bound
      while(quant.interval*quant.idx < cdf.ests[i+1] &
            (quant.idx+1)*quant.interval <= cdf.ests[i+1]){
              quant.idx <- quant.idx + 1
              quant.frac <- (quant.idx*quant.interval - cdf.ests[i])/
                (cdf.ests[i+1] - cdf.ests[i])
              quant.x <- lower.bounds[i] + quant.frac*(lower.bounds[i+1] - lower.bounds[i])
              quant.cum.temp[quant.idx] <- quant.frac*bin.ests[i]*(lower.bounds[i] + quant.x)/2 +
                cum.income.shares[i]
            }
    } else if(i < n.bin & i > median.idx){
      ## otherwise, estimate alpha
      alpha <- estimate.alpha(lower.bounds[i], lower.bounds[i+1], rev.cum.counts[i], rev.cum.counts[i+1])
      if(alpha <= 1){
        ## if alpha is estimated <= 1, revert to midpoint method
        income.shares[i] <- (lower.bounds[i+1] + lower.bounds[i])/2 * bin.ests[i]
        cum.income.shares[i+1] <- cum.income.shares[i] + income.shares[i]
        alphas[i] <- 1
        ## compute quantiles in this category
        while(quant.interval*quant.idx < cdf.ests[i+1] &
              (quant.idx+1)*quant.interval <= cdf.ests[i+1]){
      quant.idx <- quant.idx + 1
      quant.frac <- (quant.idx*quant.interval - cdf.ests[i])/
        (cdf.ests[i+1] - cdf.ests[i])
      quant.x <- lower.bounds[i] + quant.frac*(lower.bounds[i+1] - lower.bounds[i])
      quant.cum.temp[quant.idx] <- quant.frac*bin.ests[i]*(lower.bounds[i] + quant.x)/2 +
        cum.income.shares[i]
    }
      } else {
        ## otherwise, use pareto distribution to get the income share
        alphas[i] <- alpha
        income.shares[i] <- (alpha/(alpha - 1))*
          (lower.bounds[i]*rev.cum.counts[i] - lower.bounds[i+1]*rev.cum.counts[i+1])
        cum.income.shares[i+1] <- cum.income.shares[i] + income.shares[i]
        ## compute quantiles in this category
        while(quant.interval*quant.idx < cdf.ests[i+1] &
              (quant.idx+1)*quant.interval <= cdf.ests[i+1]){
                quant.idx <- quant.idx + 1
                quant.n <- count.total*(1 - quant.idx * quant.interval)
                quant.x <- lower.bounds[i]*exp(log(rev.cum.counts[i]/quant.n)/alpha)
                quant.cum.temp[quant.idx] <- cum.income.shares[i] +
                  (alpha/(alpha-1))*(lower.bounds[i]*rev.cum.counts[i] - quant.x*quant.n)
              }
      }
    } else if(i == n.bin){
      ## now deal with the upper category
      ## first use alpha from category just below
      back.idx <- n.bin - 1
      alpha <- estimate.alpha(lower.bounds[back.idx], lower.bounds[n.bin],
                              rev.cum.counts[back.idx], rev.cum.counts[n.bin])
      while(alpha <= 1 & median.idx < back.idx - 1){
      back.idx <-  back.idx - 1
      alpha <- estimate.alpha(lower.bounds[back.idx], lower.bounds[n.bin],
                              rev.cum.counts[back.idx], rev.cum.counts[n.bin])
    }
      if(alpha <= 1){
        ## if we can't get an alpha > 1, then put everything at the category min
        income.shares[n.bin] <- lower.bounds[n.bin]*bin.ests[n.bin]
        cum.income.shares[n.bin+1] <- cum.income.shares[n.bin] + income.shares[n.bin]
        alphas[n.bin] <- 1
        ## compute quantiles in this category
        while(quant.interval*quant.idx < cdf.ests[i+1] &
              (quant.idx+1)*quant.interval <= cdf.ests[i+1]){
      quant.idx <- quant.idx + 1
      quant.frac <- (quant.interval*quant.idx - cdf.ests[i])/
        (cdf.ests[i+1] - cdf.ests[i])
      quant.x <- lower.bounds[n.bin]
      quant.cum.temp[quant.idx] <- cum.income.shares[i] + quant.frac*bin.ests[i]*quant.x
    }
      } else {
        ## if we can get an alpha > 1, then use the pareto with a huge new upper bound
        alphas[n.bin] <- alpha
        end.pt <- 2*lower.bounds[n.bin]
        count.above.end.pt <- bin.ests[n.bin]*exp(alpha * log(lower.bounds[n.bin]/end.pt))
        perc.below.end.pt <- 1 - count.above.end.pt / count.total
        income.shares[n.bin] <-
          (alpha/(alpha-1))*(lower.bounds[n.bin]*bin.ests[n.bin] - end.pt*count.above.end.pt) +
          end.pt*count.above.end.pt
        cum.income.shares[n.bin+1] <- cum.income.shares[n.bin] + income.shares[n.bin]
        ## compute quantiles in this bound
        while(quant.interval*quant.idx < perc.below.end.pt &
              (quant.idx+1)*quant.interval <= perc.below.end.pt){
                quant.idx <- quant.idx + 1
                quant.n <- count.total*(1 - quant.idx * quant.interval)
                quant.x <- lower.bounds[n.bin]*
                  exp(log(rev.cum.counts[n.bin]/quant.n)/alpha)
                quant.cum.temp[quant.idx] <- cum.income.shares[n.bin] +
                  (alpha/(alpha-1))*(lower.bounds[n.bin]*rev.cum.counts[n.bin] - quant.x*quant.n)
              }
        while(quant.interval*quant.idx < 1.0 &
              quant.interval*(quant.idx+1) <= 1.0){
                quant.idx <- quant.idx + 1
                quant.frac <- (quant.idx*quant.interval - perc.below.end.pt)/
                  (1.0 - perc.below.end.pt)
                quant.x <- end.pt
                quant.cum.temp[quant.idx] <- cum.income.shares[n.bin] +
                  alpha/(alpha-1)*(lower.bounds[n.bin]*rev.cum.counts[n.bin] -
                                   end.pt*count.above.end.pt) +
                  quant.frac*quant.x*count.above.end.pt
              }
      }
    }
  }
  ## estimate the gini coefficient using the lorenz curve
  ## at the points of the income categories
  bin.probs <- rep(0, n.bin+1)
  lorenz <- rep(0, n.bin+1)
  bin.probs[n.bin+1] <- 1
  lorenz[n.bin+1] <- 1
  ginir <- bin.probs[n.bin+1]
  for(i in 2:n.bin){
    bin.probs[i] <- (count.total - rev.cum.counts[i])/count.total
    lorenz[i] <- cum.income.shares[i]/cum.income.shares[n.bin+1]
  }
  for(i in 2:n.bin){
    ginir <- ginir - lorenz[i]*(bin.probs[i+1] - bin.probs[i-1])
  }
  ## estimate the gini coefficient using 100 points of the implied lorenz curve
  quant.cum <- c(0, quant.cum.temp/cum.income.shares[n.bin + 1])
  giniq <- quant.cum[1]*(5/12) + quant.cum[2]*(13/12)
  for(i in 3:(n.quant-1)){
    giniq <- giniq + quant.cum[i]
  }
  giniq <- giniq + quant.cum[n.quant]*13/12 + quant.cum[n.quant+1]*5/12
  giniq <- giniq * quant.interval
  giniq <- (0.5 - giniq)/0.5
  n.quant <- length(target.quantiles)
  out.quantiles <- rep(0, n.quant)
  probs <- bin.ests/sum(bin.ests)
  cum.probs <- c(0, cumsum(probs))
  for(i in 1:n.quant){
    p.idx <- max(which(target.quantiles[i] > cum.probs))
    p.me <- probs[p.idx]
    p.other <- cum.probs[p.idx]
    tau.me <- (target.quantiles[i] - p.other)/p.me
    alpha <- alphas[p.idx]
    if(p.idx < n.bin){
      lb <- lower.bounds[p.idx]
      ub <- lower.bounds[p.idx + 1]
      if(alpha == 1){
        out.quantiles[i] <- lb + (ub - lb)*tau.me
      } else {
        out.quantiles[i] <- lb*exp(-log(1 - (1 - exp(alpha*log(lb/ub)))*tau.me)/alpha)
      }
    } else {
      lb <- lower.bounds[p.idx]
      if(alpha == 1){
        out.quantiles[i] <- lb
      } else {
        out.quantiles[i] <- lb*exp(-log(1 - tau.me)/alpha)
      }
    }
  }
  return(list(ginir = ginir, giniq = giniq, alphas = alphas, quants = out.quantiles))
}
