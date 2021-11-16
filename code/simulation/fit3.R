library(rstan)
source("../postprocess_functions.R")
load("standats.RData")

iters <- 334 + 333 + 1:333

model.init <- stan_model(file = "../models/prln_manytracts.stan")

target.quantiles = seq(0.05, 0.95, by = 0.05)
ncore <- 4
ncore2 <- 1 ## if not on windows, set this to 4

stat.out <- list()
set.seed(34325425)
for(iter in iters){
    print(iter)
    fit <- sampling(model.init, data = standats[[iter]],
                    cores = 4, chains = 4,
                    warmup = 4000, iter = 8000,
                    control = list(adapt_delta = 0.9, max_treedepth = 13))
    stat <- prln.tract.predict(target.quantiles, standats[[iter]], fit, ncore2)
    stat$iter <- iter
    stat.out[[iter]] <- stat
    rm(fit)
    save(stat.out, file = "stat.out3.RData")
}
