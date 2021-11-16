library(reldist)
library(parallel)
library(MASS)
library(data.table)
library(abind)
library(rstan)
library(tidyverse)

logsumexp <- function(x){
    idx <- which.max(x)
    log1p(sum(exp(x[-idx] - x[idx]))) + x[idx]
}

rprln <- function(n, alpha, pknot, knots){
    nknot <- length(pknot)
    ids <- sample(1:nknot, n, TRUE, pknot)
    lb <- knots[ids]
    ub <- knots[ids + 1]
    u <- runif(n, 0, 1)
    nalpha <- length(alpha)
    nunif <- nknot - nalpha
    out <- sapply(1:n, function(i){
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

rprlnpuma <- function(n, alpha_tract, pknot_tract, weights, nalphas, knots){
    nknot <- ncol(pknot_tract)
    ntract <- nrow(pknot_tract)
    nunifs <- nknot - nalphas
    tracts <- sample(1:ntract, n, TRUE, weights)
    sapply(tracts, function(itract){
        nunif <- nunifs[itract]
        nalpha <- nalphas[itract]
        alphaIDX <- sum(nalphas[1:itract]) - nalpha
        pknot <- pknot_tract[itract,]
        alpha <- alpha_tract[alphaIDX + 1:nalpha]
        rprln(1, alpha, pknot, knots)
    })
}

cdfprln <- function(x, alpha, pknot, knots){
    nknot <- length(pknot)
    nalpha <- length(alpha)
    nunif <- nknot - nalpha
    out <- sapply(x, function(y){
        idx <- max(which(knots <= y))
        lb <- knots[idx]
        ub <- knots[idx + 1]
        pknot <- pknot[idx]
        if(idx > nunif){
            myalpha <- alpha[idx - nunif]
            pknot * (1 -(lb / y)^(myalpha)) /
                (1 - (lb / ub)^(myalpha))
        } else {
            pknot * (y - lb) / (ub - lb)
        }
    })
    out
}

ldprln <- function(x, alpha, pknot, knots){
    nknot <- length(pknot)
    nalpha <- length(alpha)
    nunif <- nknot - nalpha
    out <- sapply(x, function(y){
        idx <- max(which(knots <= y))
        lb <- knots[idx]
        ub <- knots[idx + 1]
        mypknot <- pknot[idx]
        if(idx > nunif){
            myalpha <- alpha[idx - nunif]
            log(mypknot) + log(myalpha) + myalpha*log(lb) -
                (myalpha + 1)*log(y) -
                log(1 - (lb/ub)^myalpha)
        } else {
            log(mypknot) - log(ub - lb)
        }
    })
    out
}

ldprlnpuma <- function(x, alpha_tract, pknot_tract, weights, nalphas, knots){
    ntract <- nrow(pknot_tract)
    nknot <- ncol(pknot_tract)
    nunifs <- nknot - nalphas
    logdens <- sapply(1:ntract, function(itract){
        nunif <- nunifs[itract]
        nalpha <- nalphas[itract]
        alphaIDX <- sum(nalphas[1:itract]) - nalpha
        pknot <- pknot_tract[itract,]
        alpha <- alpha_tract[alphaIDX + 1:nalpha]
        ldprln(x, alpha, pknot, knots)
    })
    apply(logdens, 1, function(y){
        logsumexp(y + log(weights))
    })
}

compute_quantiles <- function(quants, fit, standat, ncore){
    pknot_tract <- sapply(fit, function(x){
        rstan::extract(x$result, pars = "pknot")$pknot
    }, simplify="array") %>% aperm(c(1,3,2))
    alpha_tract <- sapply(fit, function(x){
        rstan::extract(x$result, pars = "alpha")$alpha
    })
    nalphas <- sapply(alpha_tract, ncol)
    alpha_tract <- do.call(cbind, alpha_tract)

    knots <- c(standat$knots, Inf)
    nknot <- dim(pknot_tract)[3]
    ntract <- dim(pknot_tract)[2]
    niter <- dim(pknot_tract)[1]
    nquant <- length(quants)
    nunifs <- nknot - nalphas

    quantsl <- mclapply(1:niter, function(iter){
        sapply(1:ntract, function(itract){
            nunif <- nunifs[itract]
            nalpha <- nalphas[itract]
            alphaIDX <- sum(nalphas[1:itract]) - nalpha
            pknots <- pknot_tract[iter, itract,]
            cumpknots <- cumsum(pknots)
            cumpknots0 <- c(0, cumpknots)
            outs <- NULL
            for(iquant in 1:nquant){
                if(quants[iquant] < 1){
                    idx <- min(which(quants[iquant] < cumpknots))
                    padj <- (quants[iquant] - cumpknots0[idx]) / pknots[idx]
                    lb <- knots[idx]
                    ub <- knots[idx + 1]
                    if(idx > nunif){
                        alpha <- alpha_tract[iter, alphaIDX + idx - nunif]
                        out <- lb*(1 - padj*(1 - (lb/ub)^alpha))^(-1/alpha)
                    } else {
                        out <- lb + (ub - lb)*padj
                    }
                } else {
                    out <- Inf
                }
                outs <- c(outs, out)
            }
            outs
        })
    }, mc.cores = ncore)

    aperm(simplify2array(quantsl))
}

compute_gini <- function(fit, standat, ncore){
    pknot_tract <- sapply(fit, function(x){
        rstan::extract(x$result, pars = "pknot")$pknot
    }, simplify="array") %>% aperm(c(1,3,2))
    alpha_tract <- sapply(fit, function(x){
        rstan::extract(x$result, pars = "alpha")$alpha
    })
    nalphas <- sapply(alpha_tract, ncol)
    alpha_tract <- do.call(cbind, alpha_tract)

    knots <- c(standat$knots, Inf)
    nknot <- dim(pknot_tract)[3]
    ntract <- dim(pknot_tract)[2]
    niter <- dim(pknot_tract)[1]
    nunifs <- nknot - nalphas

    lower_knot_means <- (knots[1:nknot] + knots[-1])/2
    lower_knot_ilorenz <- (1 + knots[1:nknot] / (knots[1:nknot] + knots[-1]))/3

    ginisl <- mclapply(1:niter, function(iter){
        sapply(1:ntract, function(itract){
            nunif <- nunifs[itract]
            nalpha <- nalphas[itract]
            alphaIDX <- sum(nalphas[1:itract]) - nalpha + 1:nalpha
            pknots <- pknot_tract[iter, itract,]
            alphas <- alpha_tract[iter, alphaIDX]
            revcumpknots0 <- c(rev(cumsum(rev(pknots)))[-1], 0)
            lbs <- knots[nunif + 1:nalpha]
            ubs <- knots[nunif + 1 + 1:nalpha]
            upper_knot_means <- alphas / (alphas - 1) * lbs *
                (1 - (lbs/ubs)^(alphas - 1)) /
                (1 - (lbs/ubs)^(alphas))
            upper_knot_ilorenz <- 1 - alphas / (2*alphas - 1) *
                (1 - (lbs/ubs)^(2*alphas - 1)) /
                (1 - (lbs/ubs)^(alphas))
            knot_means <- c(lower_knot_means[1:nunif], upper_knot_means)
            knot_ilorenz <- c(lower_knot_ilorenz[1:nunif], upper_knot_ilorenz)
            ilorenz =
                sum(knot_means*pknots*(pknots*knot_ilorenz + revcumpknots0)) /
                sum(knot_means*pknots)
            gini <- 1 - 2*ilorenz
            gini
        })
    }, mc.cores = ncore)
    aperm(simplify2array(ginisl))
}

compute_shares <- function(quant_lb, quant_ub, fit, standat, ncore){
    pknot_tract <- sapply(fit, function(x){
        rstan::extract(x$result, pars = "pknot")$pknot
    }, simplify="array") %>% aperm(c(1,3,2))
    alpha_tract <- sapply(fit, function(x){
        rstan::extract(x$result, pars = "alpha")$alpha
    })
    nalphas <- sapply(alpha_tract, ncol)
    alpha_tract <- do.call(cbind, alpha_tract)

    knots <- c(standat$knots, Inf)
    nknot <- dim(pknot_tract)[3]
    ntract <- dim(pknot_tract)[2]
    niter <- dim(pknot_tract)[1]
    nunifs <- nknot - nalphas
    nquant <- length(quant_lb)

    lower_knot_means <- (knots[1:nknot] + knots[-1])/2
    lower_knot_ilorenz <- (1 + knots[1:nknot] / (knots[1:nknot] + knots[-1]))/3

    qlbs <- compute_quantiles(quant_lb, fit, standat, ncore)
    qubs <- compute_quantiles(quant_ub, fit, standat, ncore)

    sharesl <- mclapply(1:niter, function(iter){
        sapply(1:ntract, function(itract){
            nunif <- nunifs[itract]
            nalpha <- nalphas[itract]
            alphaIDX <- sum(nalphas[1:itract]) - nalpha + 1:nalpha
            pknots <- pknot_tract[iter, itract,]
            alphas <- alpha_tract[iter, alphaIDX]
            revcumpknots0 <- c(rev(cumsum(rev(pknots)))[-1], 0)
            lbs <- knots[nunif + 1:nalpha]
            ubs <- knots[nunif + 1 + 1:nalpha]
            upper_knot_means <- alphas / (alphas - 1) * lbs *
                (1 - (lbs/ubs)^(alphas - 1)) /
                (1 - (lbs/ubs)^(alphas))
            knot_means <- c(lower_knot_means[1:nunif], upper_knot_means)
            outs <- NULL
            for(iquant in 1:nquant){
                qlb <- qlbs[iter, itract, iquant]
                qub <- qubs[iter, itract, iquant]
                lbIDX <- max(which(standat$knots <= qlb))
                ubIDX <- max(which(standat$knots <= qub))
                IDXs <- lbIDX:ubIDX
                ## lorenz curve at lower quantile
                lb <- knots[lbIDX]
                ub <- knots[lbIDX + 1]
                muk <- knot_means[lbIDX]
                if(lbIDX <= nunif){
                    lorenz_lb <- ((qlb^2 - lb^2)/(ub - lb)) / (2*muk)
                } else {
                    alpha <- alphas[lbIDX - nunif]
                    lorenz_lb <- (lb/muk)*(alpha / (alpha - 1)) *
                        (1 - (lb/qlb)^(alpha - 1)) /
                        (1 - (lb/ub)^(alpha))
                }
                ## lorenz curve at upper quantile
                lb <- knots[ubIDX]
                ub <- knots[ubIDX + 1]
                muk <- knot_means[ubIDX]
                if(ubIDX <= nunif){
                    lorenz_ub <- ((qub^2 - lb^2)/(ub - lb)) / (2*muk)
                } else if(qub < Inf){
                    alpha <- alphas[ubIDX - nunif]
                    lorenz_ub <- (lb/muk)*(alpha / (alpha - 1)) *
                        (1 - (lb/qub)^(alpha - 1)) /
                        (1 - (lb/ub)^(alpha))
                } else {
                    lorenz_ub <- 1
                }

                out <- sum(pknots[IDXs] * knot_means[IDXs]) -
                    pknots[ubIDX] * knot_means[ubIDX] * (1 - lorenz_ub) -
                    pknots[lbIDX] * knot_means[lbIDX] * lorenz_lb

                outs <- c(outs, out / sum(pknots * knot_means))
            }
            outs
        })
    }, mc.cores = ncore)
    aperm(simplify2array(sharesl))
}

compute_lpdfs <- function(true_values, est_mat, se_mat, ncore){
    niter <- dim(true_values)[1]
    out <- mclapply(1:niter, function(iter){
        dnorm(true_values[iter,,], est_mat, se_mat, log = TRUE)
    }, mc.cores = ncore)
    aperm(simplify2array(out), c(3,1,2))
}

compute_waic <- function(standat, fit, ncore){
    target_quantiles <- c(0.2, 0.4, 0.6, 0.8, 0.95)
    true_quants <- compute_quantiles(c(0.5, target_quantiles), fit, standat, ncore)
    true_bins <- sapply(fit, function(x){
        rstan::extract(x$result,
                       pars = "tract_true_values"
                       )$tract_true_values[,1:standat$nbin]
    }, simplify="array") %>% aperm(c(1,3,2))
    true_means <- sapply(fit, function(x){
        rstan::extract(x$result,
                       pars = "tract_true_values"
                       )$tract_true_values[,1 + standat$nbin]
    }, simplify="array")
    true_gini <- compute_gini(fit, standat, ncore)
    quant_lb <- c(0, 0.2, 0.4, 0.6, 0.8, 0.95)
    quant_ub <- c(0.2, 0.4, 0.6, 0.8, 1, 1)
    true_shares <- compute_shares(quant_lb, quant_ub, fit, standat, ncore)

    true_values <- abind(true_bins, 1000*true_means,
                        1000*true_quants, true_shares, true_gini)

    est_mat <- cbind(standat$bin_est, 1000*standat$mean_est,
                    1000*standat$median_est, 1000*standat$quant_est,
                    standat$share_est, standat$gini_est)

    se_mat <- cbind(standat$bin_se, 1000*standat$mean_se,
                   1000*standat$median_se, 1000*standat$quant_se,
                   standat$share_se, standat$gini_se)

    lpdfs <- compute_lpdfs(true_values, est_mat, se_mat, ncore)
    mean_pdfs <- apply(exp(lpdfs), c(2, 3), mean, na.rm = TRUE)
    var_lpdfs <- apply(lpdfs, c(2, 3), var, na.rm = TRUE)

    full_lpdfs <- apply(lpdfs, c(1, 2), sum, na.rm = TRUE)
    mean_full_pdfs <- apply(exp(full_lpdfs), c(2), mean, na.rm = TRUE)
    var_full_lpdfs <- apply(full_lpdfs, c(2), var, na.rm = TRUE)

    model_lpdfs <- apply(lpdfs[,,1:(standat$nbin + 2)],
                        c(1, 2), sum, na.rm = TRUE)
    mean_model_pdfs <- apply(exp(model_lpdfs), c(2), mean, na.rm = TRUE)
    var_model_lpdfs <- apply(model_lpdfs, c(2), var, na.rm = TRUE)

    pred_lpdfs <- apply(lpdfs[,,-c(1:(standat$nbin + 2))],
                       c(1, 2), sum, na.rm = TRUE)
    mean_pred_pdfs <- apply(exp(pred_lpdfs), c(2), mean, na.rm = TRUE)
    var_pred_lpdfs <- apply(pred_lpdfs, c(2), var, na.rm = TRUE)

    quant_lpdfs <- apply(lpdfs[,,standat$nbin + 2 + 1:standat$nquant],
                        c(1, 2), sum, na.rm = TRUE)
    mean_quant_pdfs <- apply(exp(quant_lpdfs), c(2), mean, na.rm = TRUE)
    var_quant_lpdfs <- apply(quant_lpdfs, c(2), var, na.rm = TRUE)

    share_lpdfs <- apply(lpdfs[,,standat$nbin + 2 + standat$nshare +
                                1:standat$nshare],
                        c(1, 2), sum, na.rm = TRUE)
    mean_share_pdfs <- apply(exp(share_lpdfs), c(2), mean, na.rm = TRUE)
    var_share_lpdfs <- apply(share_lpdfs, c(2), var, na.rm = TRUE)

    waics <- log(cbind(mean_model_pdfs, mean_pred_pdfs, mean_full_pdfs,
                      mean_quant_pdfs, mean_share_pdfs, mean_pdfs)) -
        cbind(var_model_lpdfs, var_pred_lpdfs, var_full_lpdfs,
              var_quant_lpdfs, var_share_lpdfs, var_lpdfs)

    waic_est <- apply(waics, 2, sum, na.rm = TRUE)
    waic_est_se <- apply(waics, 2, function(x){
        sd(x, na.rm = TRUE)*sqrt(sum(1-is.na(x)))})

    name_bounds <- c( standat$knots, "Up")
    bin_names <- paste("Bin", name_bounds[-length(name_bounds)],
                      "--",  name_bounds[-1], sep = "")
    quant_names <- paste("P", c(20, 40, 60, 80, 95), sep="")
    share_names <- paste("Share", c(0, 20, 40, 60, 80, 95),
                        c(20, 40, 60, 80, 100, 100), sep = "--")
    waic_df <- tibble(Est = c("Base", "Extra", "Full", "Percentiles", "Shares",
                             bin_names, "Mean", "Median",
                             quant_names, share_names, "Gini"),
                     WAIC = waic_est,
                     SE = waic_est_se)
    waic_df
}

compute_kl_index <- function(nsim, fit, weights, knots, ncore){
    pknot_tract <- sapply(fit, function(x){
        rstan::extract(x$result, pars = "pknot")$pknot
    }, simplify="array") %>% aperm(c(1,3,2))
    alpha_tract <- sapply(fit, function(x){
        rstan::extract(x$result, pars = "alpha")$alpha
    })
    nalphas <- sapply(alpha_tract, ncol)
    alpha_tract <- do.call(cbind, alpha_tract)
    nknot <- dim(pknot_tract)[3]
    ntract <- dim(pknot_tract)[2]
    niter <- dim(pknot_tract)[1]
    nunifs <- nknot - nalphas
    divsl <- mclapply(1:niter, function(iter){
        sims <- rprlnpuma(nsim, alpha_tract[iter,],
                         pknot_tract[iter,,],
                         weights, nalphas, knots)
        logdens <- sapply(1:ntract, function(itract){
            alphaIDX <- sum(nalphas[1:itract]) - nalphas[itract]
            logdens <- ldprln(sims,
                             alpha_tract[iter, alphaIDX + 1:nalphas[itract]],
                             pknot_tract[iter, itract,], knots)
        })
        logdiff <- apply(logdens, 1, function(y){
            y - logsumexp(y + log(weights))
        })
        apply(logdiff, 1, function(y){
            mean(y*exp(y))
        })
    }, mc.cores = ncore)
    tract_kls <- do.call(rbind, divsl)
    weighted_kl <- tract_kls %*% weights
    out <- list(tract_kls = tract_kls, weighted_kl = weighted_kl)
    out
}

compute_hr_index <- function(nsim, fit, weights, knots, ncore){
    pknot_tract <- sapply(fit, function(x){
        rstan::extract(x$result, pars = "pknot")$pknot
    }, simplify="array") %>% aperm(c(1,3,2))

    alpha_tract <- sapply(fit, function(x){
        rstan::extract(x$result, pars = "alpha")$alpha
    })
    nalphas <- sapply(alpha_tract, ncol)
    alpha_tract <- do.call(cbind, alpha_tract)
    nknot <- dim(pknot_tract)[3]
    ntract <- dim(pknot_tract)[2]
    niter <- dim(pknot_tract)[1]
    nunifs <- nknot - nalphas
    bdivsl <- mclapply(1:niter, function(iter){
        sims <- rprlnpuma(nsim, alpha_tract[iter,],
                         pknot_tract[iter,,],
                         weights, nalphas, knots)
        tract_integrals <- sapply(1:ntract, function(itract){
            nunif <- nunifs[itract]
            nalpha <- nalphas[itract]
            alphaIDX <- sum(nalphas[1:itract]) - nalpha
            pknots <- pknot_tract[iter, itract,]
            alphas <- alpha_tract[iter, alphaIDX + 1:nalpha]
            cdfs <- cdfprln(sims, alphas, pknots, knots)
            mean(cdfs*log(cdfs) + (1-cdfs)*log(1-cdfs))
        })
    }, mc.cores = ncore)
    tract_hrs <- 1 + do.call(rbind, bdivsl)
    weighted_hr <- tract_hrs %*% weights
    out <- list(tract_hrs = tract_hrs, weighted_hr = weighted_hr)
    out
}

eval.model <- function(pred, trues){
    pred.cast <- dcast(pred, iter + tract ~ name)

    pred.diffs <- pred.cast %>%
        group_by(iter) %>%
        arrange(tract) %>%
        transmute(tract = tract,
                  q0.05.med.diff = quant..0.05...median - trues$quant..0.05[tract],
                  q0.1.med.diff  = quant..0.1...median  - trues$quant..0.1[tract],
                  q0.15.med.diff = quant..0.15...median - trues$quant..0.15[tract],
                  q0.2.med.diff  = quant..0.2...median  - trues$quant..0.2[tract],
                  q0.25.med.diff = quant..0.25...median - trues$quant..0.25[tract],
                  q0.3.med.diff  = quant..0.3...median  - trues$quant..0.3[tract],
                  q0.35.med.diff = quant..0.35...median - trues$quant..0.35[tract],
                  q0.4.med.diff  = quant..0.4...median  - trues$quant..0.4[tract],
                  q0.45.med.diff = quant..0.45...median - trues$quant..0.45[tract],
                  q0.5.med.diff  = quant..0.5...median  - trues$quant..0.5[tract],
                  q0.55.med.diff = quant..0.55...median - trues$quant..0.55[tract],
                  q0.6.med.diff  = quant..0.6...median  - trues$quant..0.6[tract],
                  q0.65.med.diff = quant..0.65...median - trues$quant..0.65[tract],
                  q0.7.med.diff  = quant..0.7...median  - trues$quant..0.7[tract],
                  q0.75.med.diff = quant..0.75...median - trues$quant..0.75[tract],
                  q0.8.med.diff  = quant..0.8...median  - trues$quant..0.8[tract],
                  q0.85.med.diff = quant..0.85...median - trues$quant..0.85[tract],
                  q0.9.med.diff  = quant..0.9...median  - trues$quant..0.9[tract],
                  q0.95.med.diff = quant..0.95...median - trues$quant..0.95[tract],
                  gini.med.diff  = gini...median        - trues$gini[tract],
                  q0.05.mean.diff = quant..0.05...mean - trues$quant..0.05[tract],
                  q0.1.mean.diff  = quant..0.1...mean  - trues$quant..0.1[tract],
                  q0.15.mean.diff = quant..0.15...mean - trues$quant..0.15[tract],
                  q0.2.mean.diff  = quant..0.2...mean  - trues$quant..0.2[tract],
                  q0.25.mean.diff = quant..0.25...mean - trues$quant..0.25[tract],
                  q0.3.mean.diff  = quant..0.3...mean  - trues$quant..0.3[tract],
                  q0.35.mean.diff = quant..0.35...mean - trues$quant..0.35[tract],
                  q0.4.mean.diff  = quant..0.4...mean  - trues$quant..0.4[tract],
                  q0.45.mean.diff = quant..0.45...mean - trues$quant..0.45[tract],
                  q0.5.mean.diff  = quant..0.5...mean  - trues$quant..0.5[tract],
                  q0.55.mean.diff = quant..0.55...mean - trues$quant..0.55[tract],
                  q0.6.mean.diff  = quant..0.6...mean  - trues$quant..0.6[tract],
                  q0.65.mean.diff = quant..0.65...mean - trues$quant..0.65[tract],
                  q0.7.mean.diff  = quant..0.7...mean  - trues$quant..0.7[tract],
                  q0.75.mean.diff = quant..0.75...mean - trues$quant..0.75[tract],
                  q0.8.mean.diff  = quant..0.8...mean  - trues$quant..0.8[tract],
                  q0.85.mean.diff = quant..0.85...mean - trues$quant..0.85[tract],
                  q0.9.mean.diff  = quant..0.9...mean  - trues$quant..0.9[tract],
                  q0.95.mean.diff = quant..0.95...mean - trues$quant..0.95[tract],
                  gini.mean.diff  = gini...mean        - trues$gini[tract],
                  q0.05.med.pctdiff = q0.05.med.diff / trues$quant..0.05[tract],
                  q0.1.med.pctdiff  = q0.1.med.diff  / trues$quant..0.1[tract],
                  q0.15.med.pctdiff = q0.15.med.diff / trues$quant..0.15[tract],
                  q0.2.med.pctdiff  = q0.2.med.diff  / trues$quant..0.2[tract],
                  q0.25.med.pctdiff = q0.25.med.diff / trues$quant..0.25[tract],
                  q0.3.med.pctdiff  = q0.3.med.diff  / trues$quant..0.3[tract],
                  q0.35.med.pctdiff = q0.35.med.diff / trues$quant..0.35[tract],
                  q0.4.med.pctdiff  = q0.4.med.diff  / trues$quant..0.4[tract],
                  q0.45.med.pctdiff = q0.45.med.diff / trues$quant..0.45[tract],
                  q0.5.med.pctdiff  = q0.5.med.diff  / trues$quant..0.5[tract],
                  q0.55.med.pctdiff = q0.55.med.diff / trues$quant..0.55[tract],
                  q0.6.med.pctdiff  = q0.6.med.diff  / trues$quant..0.6[tract],
                  q0.65.med.pctdiff = q0.65.med.diff / trues$quant..0.65[tract],
                  q0.7.med.pctdiff  = q0.7.med.diff  / trues$quant..0.7[tract],
                  q0.75.med.pctdiff = q0.75.med.diff / trues$quant..0.75[tract],
                  q0.8.med.pctdiff  = q0.8.med.diff  / trues$quant..0.8[tract],
                  q0.85.med.pctdiff = q0.85.med.diff / trues$quant..0.85[tract],
                  q0.9.med.pctdiff  = q0.9.med.diff  / trues$quant..0.9[tract],
                  q0.95.med.pctdiff = q0.95.med.diff / trues$quant..0.95[tract],
                  gini.med.pctdiff  = gini.med.diff  / trues$gini[tract],
                  q0.05.mean.pctdiff = q0.05.mean.diff / trues$quant..0.05[tract],
                  q0.1.mean.pctdiff  = q0.1.mean.diff  / trues$quant..0.1[tract],
                  q0.15.mean.pctdiff = q0.15.mean.diff / trues$quant..0.15[tract],
                  q0.2.mean.pctdiff  = q0.2.mean.diff  / trues$quant..0.2[tract],
                  q0.25.mean.pctdiff = q0.25.mean.diff / trues$quant..0.25[tract],
                  q0.3.mean.pctdiff  = q0.3.mean.diff  / trues$quant..0.3[tract],
                  q0.35.mean.pctdiff = q0.35.mean.diff / trues$quant..0.35[tract],
                  q0.4.mean.pctdiff  = q0.4.mean.diff  / trues$quant..0.4[tract],
                  q0.45.mean.pctdiff = q0.45.mean.diff / trues$quant..0.45[tract],
                  q0.5.mean.pctdiff  = q0.5.mean.diff  / trues$quant..0.5[tract],
                  q0.55.mean.pctdiff = q0.55.mean.diff / trues$quant..0.55[tract],
                  q0.6.mean.pctdiff  = q0.6.mean.diff  / trues$quant..0.6[tract],
                  q0.65.mean.pctdiff = q0.65.mean.diff / trues$quant..0.65[tract],
                  q0.7.mean.pctdiff  = q0.7.mean.diff  / trues$quant..0.7[tract],
                  q0.75.mean.pctdiff = q0.75.mean.diff / trues$quant..0.75[tract],
                  q0.8.mean.pctdiff  = q0.8.mean.diff  / trues$quant..0.8[tract],
                  q0.85.mean.pctdiff = q0.85.mean.diff / trues$quant..0.85[tract],
                  q0.9.mean.pctdiff  = q0.9.mean.diff  / trues$quant..0.9[tract],
                  q0.95.mean.pctdiff = q0.95.mean.diff / trues$quant..0.95[tract],
                  gini.mean.pctdiff  = gini.mean.diff  / trues$gini[tract]) %>%
        ungroup()

    med.rmse <- pred.diffs %>%
        summarise_at(vars(contains("med.diff")), function(x){
            sqrt(mean(x^2, na.rm = TRUE))
        }) %>%
        setnames(old = names(.),
                 new = c(paste("P", seq(5, 95, by = 5), sep=""), 'Gini')) %>%
        mutate(Estimator = "Median", Metric = "RMSE")

    mean.rmse <- pred.diffs %>%
        summarise_at(vars(contains("mean.diff")), function(x){
            sqrt(mean(x^2, na.rm = TRUE))
        }) %>%
        setnames(old = names(.),
                 new = c(paste("P", seq(5, 95, by = 5), sep=""), 'Gini')) %>%
        mutate(Estimator = "Mean", Metric = "RMSE")

    med.mad <- pred.diffs %>%
        summarise_at(vars(contains("med.diff")), function(x){
            mean(abs(x), na.rm = TRUE)
        }) %>%
        setnames(old = names(.),
                 new = c(paste("P", seq(5, 95, by = 5), sep=""), 'Gini')) %>%
        mutate(Estimator = "Median", Metric = "MAD")

    mean.mad <- pred.diffs %>%
        summarise_at(vars(contains("mean.diff")), function(x){
            mean(abs(x), na.rm = TRUE)
        }) %>%
        setnames(old = names(.),
                 new = c(paste("P", seq(5, 95, by = 5), sep=""), 'Gini')) %>%
        mutate(Estimator = "Mean", Metric = "MAD")

    med.rmspe <- pred.diffs %>%
        summarise_at(vars(contains("med.pctdiff")), function(x){
            sqrt(mean(x^2, na.rm = TRUE))
        }) %>%
        setnames(old = names(.),
                 new = c(paste("P", seq(5, 95, by = 5), sep=""), 'Gini')) %>%
        mutate(Estimator = "Median", Metric = "RMSPE")

    mean.rmspe <- pred.diffs %>%
        summarise_at(vars(contains("mean.pctdiff")), function(x){
            sqrt(mean(x^2, na.rm = TRUE))
        }) %>%
        setnames(old = names(.),
                 new = c(paste("P", seq(5, 95, by = 5), sep=""), 'Gini')) %>%
        mutate(Estimator = "Mean", Metric = "RMSPE")

    med.mape <- pred.diffs %>%
        summarise_at(vars(contains("med.pctdiff")), function(x){
            mean(abs(x), na.rm = TRUE)
        }) %>%
        setnames(old = names(.),
                 new = c(paste("P", seq(5, 95, by = 5), sep=""), 'Gini')) %>%
        mutate(Estimator = "Median", Metric = "MAPE")

    mean.mape <- pred.diffs %>%
        summarise_at(vars(contains("mean.pctdiff")), function(x){
            mean(abs(x), na.rm = TRUE)
        }) %>%
        setnames(old = names(.),
                 new = c(paste("P", seq(5, 95, by = 5), sep=""), 'Gini')) %>%
        mutate(Estimator = "Mean", Metric = "MAPE")

    out <- bind_rows(med.rmse, mean.rmse, med.rmspe, mean.rmspe,
                     med.mad, mean.mad, med.mape, mean.mape) %>%
        select(Estimator, Metric, everything())
    out
}

eval.prln <- function(prln, trues){
    pred.diffs <- prln %>%
        group_by(iter) %>%
        arrange(tract) %>%
        mutate(P5.diff  = P5  - trues$quant..0.05,
               P10.diff = P10 - trues$quant..0.1,
               P15.diff = P15 - trues$quant..0.15,
               P20.diff = P20 - trues$quant..0.2,
               P25.diff = P25 - trues$quant..0.25,
               P30.diff = P30 - trues$quant..0.3,
               P35.diff = P35 - trues$quant..0.35,
               P40.diff = P40 - trues$quant..0.4,
               P45.diff = P45 - trues$quant..0.45,
               P50.diff = P50 - trues$quant..0.5,
               P55.diff = P55 - trues$quant..0.55,
               P60.diff = P60 - trues$quant..0.6,
               P65.diff = P65 - trues$quant..0.65,
               P70.diff = P70 - trues$quant..0.7,
               P75.diff = P75 - trues$quant..0.75,
               P80.diff = P80 - trues$quant..0.8,
               P85.diff = P85 - trues$quant..0.85,
               P90.diff = P90 - trues$quant..0.9,
               P95.diff = P95 - trues$quant..0.95,
               gini.diff = gini - trues$gini,
               P5.pctdiff  = P5.diff  / trues$quant..0.05,
               P10.pctdiff = P10.diff / trues$quant..0.1,
               P15.pctdiff = P15.diff / trues$quant..0.15,
               P20.pctdiff = P20.diff / trues$quant..0.2,
               P25.pctdiff = P25.diff / trues$quant..0.25,
               P30.pctdiff = P30.diff / trues$quant..0.3,
               P35.pctdiff = P35.diff / trues$quant..0.35,
               P40.pctdiff = P40.diff / trues$quant..0.4,
               P45.pctdiff = P45.diff / trues$quant..0.45,
               P50.pctdiff = P50.diff / trues$quant..0.5,
               P55.pctdiff = P55.diff / trues$quant..0.55,
               P60.pctdiff = P60.diff / trues$quant..0.6,
               P65.pctdiff = P65.diff / trues$quant..0.65,
               P70.pctdiff = P70.diff / trues$quant..0.7,
               P75.pctdiff = P75.diff / trues$quant..0.75,
               P80.pctdiff = P80.diff / trues$quant..0.8,
               P85.pctdiff = P85.diff / trues$quant..0.85,
               P90.pctdiff = P90.diff / trues$quant..0.9,
               P95.pctdiff = P95.diff / trues$quant..0.95,
               gini.pctdiff = gini.diff / trues$gini) %>%
        ungroup()

    rmse <- pred.diffs %>%
        summarise_at(vars(contains(".diff")), function(x){
            sqrt(mean(x^2, na.rm = TRUE))
        }) %>%
        setnames(old = names(.),
                 new = c(paste("P", seq(5, 95, by = 5), sep=""), 'Gini')) %>%
        mutate(Estimator = "PRLN", Metric = "RMSE")

    mad <- pred.diffs %>%
        summarise_at(vars(contains(".diff")), function(x){
            mean(abs(x), na.rm = TRUE)
        }) %>%
        setnames(old = names(.),
                 new = c(paste("P", seq(5, 95, by = 5), sep=""), 'Gini')) %>%
        mutate(Estimator = "PRLN", Metric = "MAD")

    rmspe <- pred.diffs %>%
        summarise_at(vars(contains(".pctdiff")), function(x){
            sqrt(mean(x^2, na.rm = TRUE))
        }) %>%
        setnames(old = names(.),
                 new = c(paste("P", seq(5, 95, by = 5), sep=""), 'Gini')) %>%
        mutate(Estimator = "PRLN", Metric = "RMSPE")

    mape <- pred.diffs %>%
        summarise_at(vars(contains(".pctdiff")), function(x){
            mean(abs(x), na.rm = TRUE)
        }) %>%
        setnames(old = names(.),
                 new = c(paste("P", seq(5, 95, by = 5), sep=""), 'Gini')) %>%
        mutate(Estimator = "PRLN", Metric = "MAPE")

    out <- bind_rows(rmse, rmspe, mad, mape) %>%
        select(Estimator, Metric, everything())
    out
}

eval.direct <- function(directs, trues){
    target.quantiles <- seq(0.05, 0.95, by = 0.05)
    trues.mat <- trues %>%
        select(contains('quant'), gini) %>%
        as.matrix()

    rmse <- lapply(directs,
                   function(x){
                       (cbind(x$q.target, x$gini.est) - trues.mat)^2
                   }) %>%
        simplify2array(.) %>%
        apply(., 2, function(x){sqrt(mean(x))}) %>% t %>%
        data.frame(.) %>%
        setNames(c(paste("P", 100*target.quantiles, sep = ""), 'Gini')) %>%
        mutate(Estimator = 'Direct', Metric = 'RMSE')

    rmspe <- lapply(directs,
                   function(x){
                       ((cbind(x$q.target, x$gini.est) - trues.mat)/trues.mat)^2
                   }) %>%
        simplify2array(.) %>%
        apply(., 2, function(x){sqrt(mean(x))}) %>% t %>%
        data.frame(.) %>%
        setNames(c(paste("P", 100*target.quantiles, sep = ""), 'Gini')) %>%
        mutate(Estimator = 'Direct', Metric = 'RMSPE')

    mad <- lapply(directs,
                   function(x){
                       abs(cbind(x$q.target, x$gini.est) - trues.mat)
                   }) %>%
        simplify2array(.) %>%
        apply(., 2, function(x){mean(x)}) %>% t %>%
        data.frame(.) %>%
        setNames(c(paste("P", 100*target.quantiles, sep = ""), 'Gini')) %>%
        mutate(Estimator = 'Direct', Metric = 'MAD')

    mape <- lapply(directs,
                   function(x){
                       abs((cbind(x$q.target, x$gini.est) - trues.mat)/trues.mat)
                   }) %>%
        simplify2array(.) %>%
        apply(., 2, function(x){mean(x)}) %>% t %>%
        data.frame(.) %>%
        setNames(c(paste("P", 100*target.quantiles, sep = ""), 'Gini')) %>%
        mutate(Estimator = 'Direct', Metric = 'MAPE')

    out <- bind_rows(rmse, rmspe, mad, mape) %>%
        select(Estimator, Metric, everything())
    out
}


cover.pop <- function(pred, trues){
    pred.cast <- dcast(pred, iter + tract ~ name)

    pred.diffs <- pred.cast %>%
        group_by(iter) %>%
        arrange(tract) %>%
        transmute(tract = tract,
                  q0.05.cover =
                      1 * ((quant..0.05...min < trues$quant..0.05[tract]) &
                           (quant..0.05...max > trues$quant..0.05[tract])),
                  q0.1.cover  =
                      1 * ((quant..0.1...min < trues$quant..0.1[tract]) &
                           (quant..0.1...max > trues$quant..0.1[tract])),
                  q0.15.cover =
                      1 * ((quant..0.15...min < trues$quant..0.15[tract]) &
                           (quant..0.15...max > trues$quant..0.15[tract])),
                  q0.2.cover  =
                      1 * ((quant..0.2...min < trues$quant..0.2[tract]) &
                           (quant..0.2...max > trues$quant..0.2[tract])),
                  q0.25.cover =
                      1 * ((quant..0.25...min < trues$quant..0.25[tract]) &
                           (quant..0.25...max > trues$quant..0.25[tract])),
                  q0.3.cover  =
                      1 * ((quant..0.3...min < trues$quant..0.3[tract]) &
                           (quant..0.3...max > trues$quant..0.3[tract])),
                  q0.35.cover =
                      1 * ((quant..0.35...min < trues$quant..0.35[tract]) &
                           (quant..0.35...max > trues$quant..0.35[tract])),
                  q0.4.cover  =
                      1 * ((quant..0.4...min < trues$quant..0.4[tract]) &
                           (quant..0.4...max > trues$quant..0.4[tract])),
                  q0.45.cover =
                      1 * ((quant..0.45...min < trues$quant..0.45[tract]) &
                           (quant..0.45...max > trues$quant..0.45[tract])),
                  q0.5.cover  =
                      1 * ((quant..0.5...min < trues$quant..0.5[tract]) &
                           (quant..0.5...max > trues$quant..0.5[tract])),
                  q0.55.cover =
                      1 * ((quant..0.55...min < trues$quant..0.55[tract]) &
                           (quant..0.55...max > trues$quant..0.55[tract])),
                  q0.6.cover  =
                      1 * ((quant..0.6...min < trues$quant..0.6[tract]) &
                           (quant..0.6...max > trues$quant..0.6[tract])),
                  q0.65.cover =
                      1 * ((quant..0.65...min < trues$quant..0.65[tract]) &
                           (quant..0.65...max > trues$quant..0.65[tract])),
                  q0.7.cover  =
                      1 * ((quant..0.7...min < trues$quant..0.7[tract]) &
                           (quant..0.7...max > trues$quant..0.7[tract])),
                  q0.75.cover =
                      1 * ((quant..0.75...min < trues$quant..0.75[tract]) &
                           (quant..0.75...max > trues$quant..0.75[tract])),
                  q0.8.cover  =
                      1 * ((quant..0.8...min < trues$quant..0.8[tract]) &
                           (quant..0.8...max > trues$quant..0.8[tract])),
                  q0.85.cover =
                      1 * ((quant..0.85...min < trues$quant..0.85[tract]) &
                           (quant..0.85...max > trues$quant..0.85[tract])),
                  q0.9.cover  =
                      1 * ((quant..0.9...min < trues$quant..0.9[tract]) &
                           (quant..0.9...max > trues$quant..0.9[tract])),
                  q0.95.cover =
                      1 * ((quant..0.95...min < trues$quant..0.95[tract]) &
                           (quant..0.95...max > trues$quant..0.95[tract])),
                  gini.cover =
                      1 * ((gini...min < trues$gini[tract]) &
                           (gini...max > trues$gini[tract]))) %>%
        ungroup()
    coverage.rate <- pred.diffs %>%
        summarise_at(vars(contains("cover")), function(x){
            mean(x, na.rm = TRUE)
        }) %>%
        setnames(old = names(.),
                 new = c(paste("P", seq(5, 95, by = 5), sep=""), 'Gini'))
    coverage.rate
}

cover.prln <- function(pred, prln){
    pred.cast <- dcast(pred, iter + tract ~ name)
    pred.diffs <- pred.cast %>%
        group_by(iter) %>%
        arrange(tract) %>%
        transmute(tract = tract,
                  q0.05.cover =
                      1 * ((quant..0.05...min < prln$P5[tract]) &
                           (quant..0.05...max > prln$P5[tract])),
                  q0.1.cover  =
                      1 * ((quant..0.1...min < prln$P10[tract]) &
                           (quant..0.1...max > prln$P10[tract])),
                  q0.15.cover =
                      1 * ((quant..0.15...min < prln$P15[tract]) &
                           (quant..0.15...max > prln$P15[tract])),
                  q0.2.cover  =
                      1 * ((quant..0.2...min < prln$P20[tract]) &
                           (quant..0.2...max > prln$P20[tract])),
                  q0.25.cover =
                      1 * ((quant..0.25...min < prln$P25[tract]) &
                           (quant..0.25...max > prln$P25[tract])),
                  q0.3.cover  =
                      1 * ((quant..0.3...min < prln$P30[tract]) &
                           (quant..0.3...max > prln$P30[tract])),
                  q0.35.cover =
                      1 * ((quant..0.35...min < prln$P35[tract]) &
                           (quant..0.35...max > prln$P35[tract])),
                  q0.4.cover  =
                      1 * ((quant..0.4...min < prln$P40[tract]) &
                           (quant..0.4...max > prln$P40[tract])),
                  q0.45.cover =
                      1 * ((quant..0.45...min < prln$P45[tract]) &
                           (quant..0.45...max > prln$P45[tract])),
                  q0.5.cover  =
                      1 * ((quant..0.5...min < prln$P50[tract]) &
                           (quant..0.5...max > prln$P50[tract])),
                  q0.55.cover =
                      1 * ((quant..0.55...min < prln$P55[tract]) &
                           (quant..0.55...max > prln$P55[tract])),
                  q0.6.cover  =
                      1 * ((quant..0.6...min < prln$P60[tract]) &
                           (quant..0.6...max > prln$P60[tract])),
                  q0.65.cover =
                      1 * ((quant..0.65...min < prln$P65[tract]) &
                           (quant..0.65...max > prln$P65[tract])),
                  q0.7.cover  =
                      1 * ((quant..0.7...min < prln$P70[tract]) &
                           (quant..0.7...max > prln$P70[tract])),
                  q0.75.cover =
                      1 * ((quant..0.75...min < prln$P75[tract]) &
                           (quant..0.75...max > prln$P75[tract])),
                  q0.8.cover  =
                      1 * ((quant..0.8...min < prln$P80[tract]) &
                           (quant..0.8...max > prln$P80[tract])),
                  q0.85.cover =
                      1 * ((quant..0.85...min < prln$P85[tract]) &
                           (quant..0.85...max > prln$P85[tract])),
                  q0.9.cover  =
                      1 * ((quant..0.9...min < prln$P90[tract]) &
                           (quant..0.9...max > prln$P90[tract])),
                  q0.95.cover =
                      1 * ((quant..0.95...min < prln$P95[tract]) &
                           (quant..0.95...max > prln$P95[tract])),
                  gini.cover =
                      1 * ((gini...min < prln$gini[tract]) &
                           (gini...max > prln$gini[tract]))) %>%
        ungroup()
    coverage.rate <- pred.diffs %>%
        summarise_at(vars(contains("cover")), function(x){
            mean(x, na.rm = TRUE)
        }) %>%
        setnames(old = names(.),
                 new = c(paste("P", seq(5, 95, by = 5), sep=""), 'Gini'))
    coverage.rate
}


estimate_summary <- function(puma.df){
    est.df  <-  puma.df %>%
        transmute(GEOID = GEOID,
                  model = "Design",
                  q0.2.est = quant.20.est,
                  q0.4.est = quant.40.est,
                  q0.6.est = quant.60.est,
                  q0.8.est = quant.80.est,
                  q0.95.est = quant.95.est,
                  gini.est = gini.est)
    se.df <- puma.df %>%
        transmute(GEOID = GEOID,
                  model = "Design",
                  q0.2.se = quant.20.se,
                  q0.4.se = quant.40.se,
                  q0.6.se = quant.60.se,
                  q0.8.se = quant.80.se,
                  q0.95.se = quant.95.se,
                  gini.se = gini.se)
    ## if the SE is missing, treat the estimate as missing too
    ##  since it is low quality due to disclosure limitation
    naidx <- which(is.na(se.df), TRUE)
    if(nrow(naidx) > 0){
        for(i in 1:nrow(naidx)){
            est.df[naidx[i,1], naidx[i,2]] <- NA
        }
    }
    se.df$geometry <- NULL
    full_join(est.df, se.df) %>%
        mutate(
            q0.2.lower = q0.2.est - q0.2.se*1.96,
            q0.2.upper = q0.2.est + q0.2.se*1.96,
            q0.4.lower = q0.4.est - q0.4.se*1.96,
            q0.4.upper = q0.4.est + q0.4.se*1.96,
            q0.6.lower = q0.6.est - q0.6.se*1.96,
            q0.6.upper = q0.6.est + q0.6.se*1.96,
            q0.8.lower = q0.8.est - q0.8.se*1.96,
            q0.8.upper = q0.8.est + q0.8.se*1.96,
            q0.95.lower = q0.95.est - q0.95.se*1.96,
            q0.95.upper = q0.95.est + q0.95.se*1.96,
            gini.lower = gini.est - gini.se*1.96,
            gini.upper = gini.est + gini.se*1.96)
}

latentprln_predict <- function(target.quantiles, standat, fit, ncore){
    pars <- rstan::extract(fit, pars = c("pknot_tract", "alpha_tract"))
    tract.size <- standat$tract_pops
    bounds <- standat$bounds
    nalphas <- standat$nalpha
    knots <- standat$knots
    rplnknots <- cbind(knots, Inf)
    ntract <- length(tract.size)
    pbin <- pars$pknot_tract
    alpha <- pars$alpha_tract
    nbin <- dim(pbin)[3]
    nsim <- dim(pbin)[1]
    stat.dfs <- mclapply(1:ntract, function(tract){
        stat.out <- sapply(1:nsim, function(sim){
            tractids <- sum(nalphas[1:tract]) - nalphas[tract] + 1:nalphas[tract]
            pred.pop <- rprln(tract.size[tract], alpha[sim,tractids],
                              pbin[sim,tract,], rplnknots[tract,])
            cdfs <- c(0, sapply(bounds, function(y){mean(pred.pop<y)}), 1)
            out <- c(quantile(pred.pop, probs = target.quantiles),
                     mean(pred.pop), diff(cdfs), gini(pred.pop))
            names(out) <- c(paste("quant", target.quantiles, sep=".."), "mean",
                            paste("bin", c("below", bounds),
                                  c(bounds, "above"), sep = ".."),
                            "gini")
            out
        })
        stat.median <- apply(stat.out, 1, median)
        stat.mean <- apply(stat.out, 1, mean)
        stat.min <- apply(stat.out, 1, quantile, probs = 0.025)
        stat.max <- apply(stat.out, 1, quantile, probs = 0.975)
        stat.lower <- apply(stat.out, 1, quantile, probs = 0.25)
        stat.upper <- apply(stat.out, 1, quantile, probs = 0.75)

        stat.combined <- c(stat.median, stat.mean, stat.lower, stat.upper,
                           stat.min, stat.max)
        names(stat.combined) <- c(paste(names(stat.median), "median", sep = "..."),
                                  paste(names(stat.mean), "mean", sep = "..."),
                                  paste(names(stat.lower), "lower", sep = "..."),
                                  paste(names(stat.upper), "upper", sep = "..."),
                                  paste(names(stat.min), "min", sep = "..."),
                                  paste(names(stat.max), "max", sep = "..."))
        data.frame(tract = tract, name = names(stat.combined), value = stat.combined)
    }, mc.cores = ncore)
    out <- Reduce(rbind, stat.dfs)
    out
}

eval_lprln_preds_acspuma <- function(pred, ests){
    pred.cast <- dcast(pred, tract ~ name)

    pred.diffs <- pred.cast %>%
        transmute(tract = factor(tract),
                  q0.2.med.diff = quant..0.2...median - ests$q0.2.est[tract],
                  q0.4.med.diff = quant..0.4...median - ests$q0.4.est[tract],
                  q0.6.med.diff = quant..0.6...median - ests$q0.6.est[tract],
                  q0.8.med.diff = quant..0.8...median - ests$q0.8.est[tract],
                  q0.95.med.diff = quant..0.95...median - ests$q0.95.est[tract],
                  gini.med.diff = gini...median - ests$gini.est[tract],
                  q0.2.mean.diff = quant..0.2...mean - ests$q0.2.est[tract],
                  q0.4.mean.diff = quant..0.4...mean - ests$q0.4.est[tract],
                  q0.6.mean.diff = quant..0.6...mean - ests$q0.6.est[tract],
                  q0.8.mean.diff = quant..0.8...mean - ests$q0.8.est[tract],
                  q0.95.mean.diff = quant..0.95...mean - ests$q0.95.est[tract],
                  gini.mean.diff = gini...mean - ests$gini.est[tract],
                  q0.2.med.pctdiff = 100 * q0.2.med.diff / ests$q0.2.est[tract],
                  q0.4.med.pctdiff = 100 * q0.4.med.diff / ests$q0.4.est[tract],
                  q0.6.med.pctdiff = 100 * q0.6.med.diff / ests$q0.6.est[tract],
                  q0.8.med.pctdiff = 100 * q0.8.med.diff / ests$q0.8.est[tract],
                  q0.95.med.pctdiff = 100 * q0.95.med.diff / ests$q0.95.est[tract],
                  gini.med.pctdiff = 100 * gini.med.diff / ests$gini.est[tract],
                  q0.2.mean.pctdiff = 100 * q0.2.mean.diff / ests$q0.2.est[tract],
                  q0.4.mean.pctdiff = 100 * q0.4.mean.diff / ests$q0.4.est[tract],
                  q0.6.mean.pctdiff = 100 * q0.6.mean.diff / ests$q0.6.est[tract],
                  q0.8.mean.pctdiff = 100 * q0.8.mean.diff / ests$q0.8.est[tract],
                  q0.95.mean.pctdiff = 100 * q0.95.mean.diff / ests$q0.95.est[tract],
                  gini.mean.pctdiff = 100 * gini.mean.diff / ests$gini.est[tract])

    med.rmse <- pred.diffs %>%
        summarise_at(vars(contains("med.diff")), function(x){
            1000*sqrt(mean(x^2, na.rm = TRUE))
        }) %>%
        setnames(old = names(.),
                 new = c(paste("Q", c(0.2, 0.4, 0.6, 0.8, 0.95), sep=""), "Gini")) %>%
        mutate(Estimator = "Median", Metric = "RMSE")

    mean.rmse <- pred.diffs %>%
        summarise_at(vars(contains("mean.diff")), function(x){
            1000*sqrt(mean(x^2, na.rm = TRUE))
        }) %>%
        setnames(old = names(.),
                 new = c(paste("Q", c(0.2, 0.4, 0.6, 0.8, 0.95), sep=""), "Gini")) %>%
        mutate(Estimator = "Mean", Metric = "RMSE")

    med.mad <- pred.diffs %>%
        summarise_at(vars(contains("med.diff")), function(x){
            1000*mean(abs(x), na.rm = TRUE)
        }) %>%
        setnames(old = names(.),
                 new = c(paste("Q", c(0.2, 0.4, 0.6, 0.8, 0.95), sep=""), "Gini")) %>%
        mutate(Estimator = "Median", Metric = "MAD")

    mean.mad <- pred.diffs %>%
        summarise_at(vars(contains("mean.diff")), function(x){
            1000*mean(abs(x), na.rm = TRUE)
        }) %>%
        setnames(old = names(.),
                 new = c(paste("Q", c(0.2, 0.4, 0.6, 0.8, 0.95), sep=""), "Gini")) %>%
        mutate(Estimator = "Mean", Metric = "MAD")

    med.rmspe <- pred.diffs %>%
        summarise_at(vars(contains("med.pctdiff")), function(x){
            sqrt(mean(x^2, na.rm = TRUE))
        }) %>%
        setnames(old = names(.),
                 new = c(paste("Q", c(0.2, 0.4, 0.6, 0.8, 0.95), sep=""), "Gini")) %>%
        mutate(Estimator = "Median", Metric = "RMSPE")

    mean.rmspe <- pred.diffs %>%
        summarise_at(vars(contains("mean.pctdiff")), function(x){
            sqrt(mean(x^2, na.rm = TRUE))
        }) %>%
        setnames(old = names(.),
                 new = c(paste("Q", c(0.2, 0.4, 0.6, 0.8, 0.95), sep=""), "Gini")) %>%
        mutate(Estimator = "Mean", Metric = "RMSPE")

    med.mape <- pred.diffs %>%
        summarise_at(vars(contains("med.pctdiff")), function(x){
            mean(abs(x), na.rm = TRUE)
        }) %>%
        setnames(old = names(.),
                 new = c(paste("Q", c(0.2, 0.4, 0.6, 0.8, 0.95), sep=""), "Gini")) %>%
        mutate(Estimator = "Median", Metric = "MAPE")

    mean.mape <- pred.diffs %>%
        summarise_at(vars(contains("mean.pctdiff")), function(x){
            mean(abs(x), na.rm = TRUE)
        }) %>%
        setnames(old = names(.),
                 new = c(paste("Q", c(0.2, 0.4, 0.6, 0.8, 0.95), sep=""), "Gini")) %>%
        mutate(Estimator = "Mean", Metric = "MAPE")

    out <- bind_rows(med.rmse, mean.rmse, med.rmspe, mean.rmspe,
                     med.mad, mean.mad, med.mape, mean.mape) %>%
        select(Estimator, Metric, everything())
    out
}

eval_prln_preds_acspuma <- function(prln.est, ests){
    pred.diffs <- data.frame(prln.est, tract = 1:nrow(prln.est)) %>%
        transmute(tract = factor(tract),
                  q0.2.diff  = Q20  - ests$q0.2.est[tract],
                  q0.4.diff  = Q40  - ests$q0.4.est[tract],
                  q0.6.diff  = Q60  - ests$q0.6.est[tract],
                  q0.8.diff  = Q80  - ests$q0.8.est[tract],
                  q0.95.diff = Q95  - ests$q0.95.est[tract],
                  gini.diff  = Gini - ests$gini.est[tract],
                  q0.2.pctdiff  = 100 * q0.2.diff / ests$q0.2.est[tract],
                  q0.4.pctdiff  = 100 * q0.4.diff / ests$q0.4.est[tract],
                  q0.6.pctdiff  = 100 * q0.6.diff / ests$q0.6.est[tract],
                  q0.8.pctdiff  = 100 * q0.8.diff / ests$q0.8.est[tract],
                  q0.95.pctdiff = 100 * q0.95.diff / ests$q0.95.est[tract],
                  gini.pctdiff  = 100 * gini.diff / ests$gini.est[tract])

    rmse <- pred.diffs %>%
        summarise_at(vars(contains(".diff")), function(x){
            1000*sqrt(mean(x^2, na.rm = TRUE))
        }) %>%
        setnames(old = names(.),
                 new = c(paste("Q", c(0.2, 0.4, 0.6, 0.8, 0.95), sep=""), "Gini")) %>%
        mutate(Estimator = "PRLN", Metric = "RMSE")

    mad <- pred.diffs %>%
        summarise_at(vars(contains(".diff")), function(x){
            1000*mean(abs(x), na.rm = TRUE)
        }) %>%
        setnames(old = names(.),
                 new = c(paste("Q", c(0.2, 0.4, 0.6, 0.8, 0.95), sep=""), "Gini")) %>%
        mutate(Estimator = "PRLN", Metric = "MAD")

    rmspe <- pred.diffs %>%
        summarise_at(vars(contains(".pctdiff")), function(x){
            sqrt(mean(x^2, na.rm = TRUE))
        }) %>%
        setnames(old = names(.),
                 new = c(paste("Q", c(0.2, 0.4, 0.6, 0.8, 0.95), sep=""), "Gini")) %>%
        mutate(Estimator = "PRLN", Metric = "RMSPE")

    mape <- pred.diffs %>%
        summarise_at(vars(contains(".pctdiff")), function(x){
            mean(abs(x), na.rm = TRUE)
        }) %>%
        setnames(old = names(.),
                 new = c(paste("Q", c(0.2, 0.4, 0.6, 0.8, 0.95), sep=""), "Gini")) %>%
        mutate(Estimator = "PRLN", Metric = "MAPE")

    out <- bind_rows(rmse, rmspe, mad, mape) %>%
        select(Estimator, Metric, everything())
    out
}
