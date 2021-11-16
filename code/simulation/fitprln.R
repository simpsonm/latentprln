load("standats.RData")
source("../pareto_linear/pareto_linear.R")


iters <- 1:1000

target.quantiles = seq(0.05, 0.95, by = 0.05)

prln.out <- NULL

set.seed(7812384)

for(iter in iters){
    bounds <- standats[[iter]]$bounds
    for(i in 1:standats[[iter]]$ntract){
        temp <- pareto.linear(bounds, standats[[iter]]$bin_est[i,], target.quantiles)
        prln.est <- matrix(c(temp$quants, temp$giniq), nrow = 1)
        colnames(prln.est) <- c(paste("P", target.quantiles*100, sep=""), 'gini')
        prln.out <- rbind(prln.out, data.frame(iter = iter, tract = i, prln.est))
    }
}

save(prln.out, file = "prln.out.RData")
