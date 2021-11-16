library(tidyverse)

logsumexp = function(x){
    idx = which.max(x)
    log1p(sum(exp(x[-idx] - x[idx]))) + x[idx]
}

## assumes that knots[length(knots)] = Inf
rprln = function(n, alpha, pknot, knots){
    nknot = length(pknot)
    ids = sample(1:nknot, n, TRUE, pknot)
    lb = knots[ids]
    ub = knots[ids + 1]
    u = runif(n, 0, 1)
    nalpha = length(alpha)
    nunif = nknot - nalpha
    out = sapply(1:n, function(i){
        if(ids[i] > nunif){
            lb[i] /
                (1 - u[i]*
                 (1 - (lb[i]/ub[i])^alpha[ids[i] - nunif])
                )^(1/alpha[ids[i] - nunif])
        } else {
            lb[i] + (ub[i] - lb[i])*u[i]
        }
    })
    out
}

## assumes that knots[length(knots)] = Inf
ldprln = function(x, alpha, pknot, knots){
    nknot = length(pknot)
    nalpha = length(alpha)
    nunif = nknot - nalpha
    out = sapply(x, function(y){
        idx = max(which(knots <= y))
        lb = knots[idx]
        ub = knots[idx + 1]
        mypknot = pknot[idx]
        if(idx > nunif){
            myalpha = alpha[idx - nunif]
            log(mypknot) + log(myalpha) + myalpha*log(lb) -
                (myalpha + 1)*log(y) -
                log(1 - (lb/ub)^myalpha)
        } else {
            log(mypknot) - log(ub - lb)
        }
    })
    out
}

## assumes that knots[length(knots)] = Inf
cdfprln = function(x, alpha, pknot, knots){
    nknot = length(pknot)
    nalpha = length(alpha)
    nunif = nknot - nalpha
    out = sapply(x, function(y){
        idx = max(which(knots <= y))
        lb = knots[idx]
        ub = knots[idx + 1]
        pknot = pknot[idx]
        if(idx > nunif){
            myalpha = alpha[idx - nunif]
            pknot * (1 -(lb / y)^(myalpha)) /
                (1 - (lb / ub)^(myalpha))
        } else {
            pknot * (y - lb) / (ub - lb)
        }
    })
    out
}

## simulates n samples per bin, does not need bin probs
rprln_bin = function(n, alpha, knots, iknot){
    nknot = length(knots)
    if(knots[nknot] == Inf) nknot = nknot - 1
    nalpha = length(alpha)
    nunif = nknot - nalpha
    u = runif(n, 0, 1)
    lb = knots[iknot]
    ub = knots[iknot + 1]
    if(iknot > nunif){
        lb / ( 1 - u*(1 - (lb/ub)^alpha[iknot - nunif]) )^(1/alpha[iknot - nunif])
    } else {
        lb + (ub - lb)*u
    }
}

## computes the log bin density for each observation,
## does not need bin probs
ldprln_bin = function(x, alpha, knots){
    nknot = length(knots)
    if(knots[nknot] == Inf) nknot = nknot - 1
    nalpha = length(alpha)
    nunif = nknot - nalpha
    out = sapply(x, function(y){
        idx = max(which(knots <= y))
        lb = knots[idx]
        ub = knots[idx + 1]
        if(idx > nunif){
            myalpha = alpha[idx - nunif]
            log(myalpha) + myalpha*log(lb) -
                (myalpha + 1)*log(y) -
                log(1 - (lb/ub)^myalpha)
        } else {
            - log(ub - lb)
        }
    })
    out
}

compute_kl_index = function(nsim, fit, weights, knots){
    ntract = length(fit)
    niters = sapply(fit, function(x){dim(x)[[1]]})
    nchains = sapply(fit, function(x){dim(x)[[2]]})
    if(length(unique(unlist(niters))) == 1){
        niter = unique(unlist(niters))
    } else {
        stop("Unbalanced sample")
    }
    if(length(unique(unlist(nchains))) == 1){
        nchain = unique(unlist(nchains))
    } else {
        stop("Unbalanced sample")
    }
    parnames = sapply(fit, function(x){attr(x, "dimnames")$parameters})
    nalphas = sapply(parnames, function(x){
        sum(grepl("alpha", x))
    })
    nknot = length(knots) - 1
    nunifs = nknot - nalphas
    nonzero_weights = weights[weights > 0]
    nonzero_tracts = which(weights > 0)
    ntract0 = length(nonzero_tracts)
    logweights = log(nonzero_weights)
    neasy = min(nunifs)   ## number of easy bins to integrate
    nhard = nknot - neasy ## number of hard bins to integrate
    divsl = sapply(1:nchain, function(ichain){
        divs = sapply(1:niter, function(iter){
            alphas = sapply(1:ntract0, function(itract0){
                itract = nonzero_tracts[itract0]
                alpha = fit[[itract]][iter, ichain, nknot + nalphas[itract]]
            })
            idx_tract = which.min(alphas)
            itract_idx = nonzero_tracts[idx_tract]
            alpha_idx = fit[[itract_idx]][iter, ichain, nknot + 1:nalphas[itract_idx]]
            ests = lapply(1:nhard, function(iknot0){
                iknot = neasy + iknot0
                sims = rprln_bin(nsim, alpha_idx, knots, iknot)
                logdens = sapply(1:ntract0, function(itract0){
                    itract = nonzero_tracts[itract0]
                    alpha = fit[[itract]][iter, ichain, nknot + 1:nalphas[itract]]
                    pknot = fit[[itract]][iter, ichain, 1:nknot]
                    logdens = ldprln_bin(sims, alpha, knots)
                })
                logprobs = sapply(1:ntract0, function(itract0){
                    itract = nonzero_tracts[itract0]
                    log(fit[[itract]][iter, ichain, iknot])
                })
                logdiff = apply(apply(logdens, 1,
                                      function(x){x + logprobs}),
                                2,
                                function(x){x - logsumexp(x + logweights) })
                logimpweight = apply(logdens, 2,
                                     function(x){x - logdens[,idx_tract]})
                integrand = logdiff*exp(apply(logimpweight, 1,
                                              function(x){x + logprobs}))
                mcests = apply(integrand, 1, mean)
                covest = cov(t(integrand))
                list(mcests = mcests, covest = covest)
            })
            easy_pknot = sapply(1:ntract0, function(itract0){
                itract = nonzero_tracts[itract0]
                fit[[itract]][iter, ichain, 1:neasy]
            })
            log_weighted_pknot = log(drop(easy_pknot %*% nonzero_weights))
            if(length(dim(easy_pknot)) > 0){
                easyests = apply(easy_pknot, 2, function(x){
                    x * (log(x) - log_weighted_pknot)}) %>% t()
            } else {
                easyests = easy_pknot * (log(easy_pknot) - log_weighted_pknot)
            }
            mcests = sapply(ests, function(x){x$mcests})
            fullests = apply(cbind(mcests, easyests), 1, sum)
            est = nonzero_weights %*% fullests
            covest = ests[[1]]$covest
            if(nhard > 1){
                for(i in 2:nhard){
                    covest = covest + ests[[i]]$covest
                }
            }
            se = sqrt((t(nonzero_weights) %*% covest %*% nonzero_weights) / nsim)
            c(est, se)
        })
    }, simplify="array")
    out = list(kl_index_ests = divsl[1,,],
               kl_index_ses = divsl[2,,])
    out
}

compute_hr_index = function(nsim, fit, weights, knots){
    ntract = length(fit)
    niters = sapply(fit, function(x){dim(x)[[1]]})
    nchains = sapply(fit, function(x){dim(x)[[2]]})
    if(length(unique(unlist(niters))) == 1){
        niter = unique(unlist(niters))
    } else {
        stop("Unbalanced sample")
    }
    if(length(unique(unlist(nchains))) == 1){
        nchain = unique(unlist(nchains))
    } else {
        stop("Unbalanced sample")
    }
    parnames = sapply(fit, function(x){attr(x, "dimnames")$parameters})
    nalphas = sapply(parnames, function(x){
        sum(grepl("alpha", x))
    })
    nknot = length(knots) - 1
    nunifs = nknot - nalphas
    nonzero_weights = weights[weights > 0]
    nonzero_tracts = which(weights > 0)
    ntract0 = length(nonzero_tracts)
    logweights = log(nonzero_weights)
    bdivsl = sapply(1:nchain, function(ichain){
        bdivs = sapply(1:niter, function(iter){
            alphas = sapply(1:ntract0, function(itract0){
                itract = nonzero_tracts[itract0]
                alpha = fit[[itract]][iter, ichain, nknot + nalphas[itract]]
            })
            idx_tract = which.min(alphas)
            itract_idx = nonzero_tracts[idx_tract]
            pknot_idx = rep(1/nknot, nknot) ## need to ensure pknots are nonzero
            alpha_idx = fit[[itract_idx]][iter, ichain, nknot + 1:nalphas[itract_idx]]
            sims = rprln(nsim, alpha_idx, pknot_idx, knots)
            logdens_idx = ldprln(sims, alpha_idx, pknot_idx, knots)
            logdens = sapply(1:ntract0, function(itract0){
                itract = nonzero_tracts[itract0]
                alpha = fit[[itract]][iter, ichain, nknot + 1:nalphas[itract]]
                pknot = fit[[itract]][iter, ichain, 1:nknot]
                ldprln(sims, alpha, pknot, knots)
            })
            logmetrodens = apply(logdens, 1,
                                 function(x){logsumexp(x + logweights)})
            impweight = exp(logmetrodens - logdens_idx)
            integrand = sapply(1:ntract0, function(itract0){
                itract = nonzero_tracts[itract0]
                alpha = fit[[itract]][iter, ichain, nknot + 1:nalphas[itract]]
                pknot = fit[[itract]][iter, ichain, 1:nknot]
                cdfs = cdfprln(sims, alpha, pknot, knots)
                (cdfs * log(cdfs) + (1 - cdfs) * log(1 - cdfs)) * impweight
            })
            mcest = apply(integrand, 2, mean)
            covest = cov(integrand)
            est = 1 + 2 * (mcest %*% nonzero_weights)
            se = 2 * sqrt((t(nonzero_weights) %*% covest %*% nonzero_weights) / nsim)
            c(est, se)
        })
    }, simplify = "array")
    out = list(hr_index_ests = bdivsl[1,,],
               hr_index_ses = bdivsl[2,,])    
    out
}



reg_model_fits = function(standat){
    reg_model = stan_model("../models/regression_noprior.stan")
    nreg = 1 + standat$nx
    out = list()
    for(i in 1:standat$nsample){
        mystandat = standat
        mystandat$index = standat$index[i,]
        mystandat$index_se = standat$index_se[i,]
        X = cbind(1, mystandat$x)
        Y = mystandat$index
        olsest = chol2inv(chol(crossprod(X, X)))%*%crossprod(X,Y)
        initlist = list(alpha = olsest[1], beta = olsest[-1],
                        sigma2 = mean((Y - X%*%olsest)^2))
        o = optimizing(reg_model, data = mystandat, hessian=TRUE,
                       init = initlist,
                       tol_obj = 1e-20,
                       tol_rel_obj = 1e0,
                       tol_grad = 1e-16,
                       tol_rel_grad = 1e-1,
                       tol_param = 1e-5,
                       iter = 10000)
        print(c(i, o$return_code, o$value))
        cnt = 1
        while(o$return_code){
            cnt = cnt + 1
            initlist = list(
                alpha = o$par[1],
                beta = o$par[1 + 1:mystandat$nx],
                sigma2 = o$par[1 + mystandat$nx + 1]
            )
            o = optimizing(reg_model, data = mystandat, hessian=TRUE,
                           init = initlist,
                           tol_obj = 1e-20,
                           tol_rel_obj = 1e0,
                           tol_grad = 1e-16,
                           tol_rel_grad = 1e-1,
                           tol_param = 1e-5,
                           iter = 10000)
            print(c(i, cnt, o$return_code, o$value))
        }
        out[[i]] = list(est = o$par[1:nreg],
                        hess = o$hess[1:nreg, 1:nreg],
                        return = o$return_code)
        rm(o)
    }
    out
}    


eiv_reg_model_fits = function(standat){
    reg_model = stan_model("../models/eiv_regression_noprior.stan")
    nreg = 1 + standat$nx
    out = list()
    for(i in 1:standat$nsample){
        mystandat = standat
        mystandat$index = standat$index[i,]
        mystandat$index_se = standat$index_se[i,]
        X = cbind(1, mystandat$x)
        Y = mystandat$index
        olsest = chol2inv(chol(crossprod(X, X)))%*%crossprod(X,Y)
        initlist = list(alpha = olsest[1], beta = olsest[-1],
                        mu = apply(mystandat$x, 2, mean),
                        sigma2 = apply(mystandat$x, 2, var),
                        tau2 = mean((Y - X%*%olsest)^2))
        o = optimizing(reg_model, data = mystandat, hessian=TRUE,
                       init = initlist,
                       tol_obj = 1e-25,
                       tol_rel_obj = 1e-3,
                       tol_grad = 1e-20,
                       tol_rel_grad = 1e-3,
                       tol_param = 1e-8,
                       iter = 10000)
        print(c(i, o$return_code, o$value))
        cnt = 1
        while(o$return_code){
            cnt = cnt + 1
            initlist = list(
                alpha = o$par[1],
                beta = o$par[1 + 1:mystandat$nx],
                mu = o$par[1 + mystandat$nx + 1:mystandat$nx],
                sigma2 = o$par[1 + 2*mystandat$nx + 1:mystandat$nx],
                tau2 = o$par[1 + 3*mystandat$nx + 1]
            )
            o = optimizing(reg_model, data = mystandat, hessian=TRUE,
                           init = initlist,
                           tol_obj = 1e-25,
                           tol_rel_obj = 1e-3,
                           tol_grad = 1e-20,
                           tol_rel_grad = 1e-3,
                           tol_param = 1e-8,
                           iter = 10000)
            print(c(i, cnt, o$return_code, o$value))
        }
        out[[i]] = list(est = o$par[1:nreg],
                        hess = o$hess[1:nreg, 1:nreg],
                        return = o$return_code)
        rm(o)
    }
    out
}    

mysummary = function(x){
    c(mean = mean(x), sd = sd(x),
      quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))}

regsummary = function(regfits, regnames){
    out = sapply(regfits, function(x){
        sigma = chol2inv(chol(-x$hess))
        chol_sigma = chol(sigma)
        mu = x$est
        mu + crossprod(chol_sigma, rnorm(length(x$est)))
    }) %>% apply(., 1, mysummary) %>% t()
    rownames(out) = regnames
    out
}
