library(tidyverse)
library(rstan)
library(xtable)
library(sf)
library(rstan)
library(reldist)  ## for computing gini
#source("postprocess_functions.R")
load("../data/metro_dataset.RData")
load("../data/metro_race_dataset.RData")
load("../data/metro_race_income_dataset.RData")
load("../data/metro_race_income_standats.RData")

metro_kls_black = list()
metro_kls_white = list()
metro_kls_combined = list()
for(i in 1:4){
    load(paste("metro_kls_black_", i, ".RData", sep = ""))
    metro_kls_black[[i]] = metro_kls
    load(paste("metro_kls_white_", i, ".RData", sep = ""))
    metro_kls_white[[i]] = metro_kls
    load(paste("metro_kls_combined_", i, ".RData", sep = ""))
    metro_kls_combined[[i]] = metro_kls
}

warnings = list()
indicies = list()
for(race in c("black", "white", "combined")){
    warnings[[race]] = list()
    indicies[[race]] = list()
}

for(i in 1:4){
    warnings[["black"]][[i]] = lapply(metro_kls_black[[i]],
                                      function(x){x$warnings})
    indicies[["black"]][[i]] = lapply(metro_kls_black[[i]],
                                      function(x){x$kls}) %>%
        reduce(., rbind)
    warnings[["white"]][[i]] = lapply(metro_kls_white[[i]],
                                      function(x){x$warnings})
    indicies[["white"]][[i]] = lapply(metro_kls_white[[i]],
                                      function(x){x$kls}) %>%
        reduce(., rbind)
    warnings[["combined"]][[i]] = lapply(metro_kls_combined[[i]],
                                         function(x){x$warnings})
    indicies[["combined"]][[i]] = lapply(metro_kls_combined[[i]],
                                         function(x){x$kls}) %>%
        reduce(., rbind)
}

## only a few fits with a small number of divergent transitions
for(race in c("black", "white", "combined")){
    for(i in 1:4){
        for(j in 1:length(warnings[[race]][[i]])){
            if(length(warnings[[race]][[i]][[1]][[j]]) &&
               !is.na(warnings[[race]][[i]][[1]][[j]])){
                print(warnings[[race]][[i]][[1]][[j]])
            }
        }
    }
}




full_indicies = lapply(indicies, function(x){reduce(x, rbind)}) %>%
    reduce(., rbind)
rm(indicies)

eiv_indicies = full_indicies %>%
    select(-chain, -iter) %>%
    group_by(metro, race) %>%
    summarise(kl.est = mean(klest),
              hr.est = mean(hrest),
              kl.se = sqrt(mean(klse^2) + var(klest)),
              hr.se = sqrt(mean(klse^2) + var(klest)),
              ntract = mean(ntract),
              ntract_excluded = mean(nas),
              ntract_included = ntract - ntract_excluded) %>%
    ungroup()

kl_summary = full_indicies %>%
    group_by(metro, year, race) %>%
    summarize(p025 = quantile(klest, probs = 0.025, na.rm = TRUE),
              p250 = quantile(klest, probs = 0.25, na.rm = TRUE),
              p500 = quantile(klest, probs = 0.5, na.rm = TRUE),
              p750 = quantile(klest, probs = 0.75, na.rm = TRUE),
              p975 = quantile(klest, probs = 0.975, na.rm = TRUE)) %>%
    mutate(index = "Divergence Index") %>%
    select(metro, year, race, index, everything())

hr_summary = full_indicies %>%
    group_by(metro, year, race) %>%
    summarize(p025 = quantile(hrest, probs = 0.025, na.rm = TRUE),
              p250 = quantile(hrest, probs = 0.25, na.rm = TRUE),
              p500 = quantile(hrest, probs = 0.5, na.rm = TRUE),
              p750 = quantile(hrest, probs = 0.75, na.rm = TRUE),
              p975 = quantile(hrest, probs = 0.975, na.rm = TRUE)) %>%
    mutate(index = "Rank-Order Information Theory Index") %>%
    select(metro, year, race, index, everything())

combined_summary = rbind(kl_summary, hr_summary)

level_order = kl_summary %>%
    filter(race == "combined") %>%
    (function(x){x$metro[order(x$p500)]})

## Create Figures H.4, H.5, and H.6
index_plot_combined = ggplot(filter(combined_summary, race == "combined"),
                    aes(x = factor(metro, levels = level_order))) +
    geom_boxplot(aes(ymin = p025, lower = p250, middle = p500,
                     upper = p750, ymax = p975), stat = "identity") +
    facet_wrap(~index, scales="free", nrow = 2) +
    xlab("Metro Areas Ordered by Divergence Index") +
    theme(axis.ticks.x=element_blank(),
          axis.text.x=element_blank())

index_plot_combined

ggsave(filename = "../../doc/index_plot_combined.png",
       plot = index_plot_combined, width = 8, height = 4)

level_order = kl_summary %>%
    filter(race == "black") %>%
    (function(x){x$metro[order(x$p500)]})

index_plot_black = ggplot(filter(combined_summary, race == "black"),
                    aes(x = factor(metro, levels = level_order))) +
    geom_boxplot(aes(ymin = p025, lower = p250, middle = p500,
                     upper = p750, ymax = p975), stat = "identity") +
    facet_wrap(~index, scales="free", nrow = 2) +
    xlab("Metro Areas Ordered by Divergence Index") +
    theme(axis.ticks.x=element_blank(),
          axis.text.x=element_blank())

index_plot_black

ggsave(filename = "../../doc/index_plot_black.png",
       plot = index_plot_black, width = 8, height = 4)

level_order = kl_summary %>%
    filter(race == "white") %>%
    (function(x){x$metro[order(x$p500)]})

index_plot_white = ggplot(filter(combined_summary, race == "white"),
                    aes(x = factor(metro, levels = level_order))) +
    geom_boxplot(aes(ymin = p025, lower = p250, middle = p500,
                     upper = p750, ymax = p975), stat = "identity") +
    facet_wrap(~index, scales="free", nrow = 2) +
    xlab("Metro Areas Ordered by Divergence Index") +
    theme(axis.ticks.x=element_blank(),
          axis.text.x=element_blank())

index_plot_white

ggsave(filename = "../../doc/index_plot_white.png",
       plot = index_plot_white, width = 8, height = 4)



#### now set up for EIV regression

### first: estimate gini coefficients for each race in each metro area

## to do that, first fit metro-level LPRLN income distributionsa
metros = unique(metro_race_income_dataset$metro)
years = c(2018)
races = c("black", "white")

model_init = stan_model(file = "../models/prln_tract.stan")
gini_fits = list()
set.seed(238993490)
for(metroName in metros){
    gini_fits[[metroName]] = list()
    for(curr_year in years){
        yearString = paste("y", curr_year, sep="")
        gini_fits[[metroName]][[yearString]] = list()
        for(race in races){
            standat = metro_race_income_standats[[metroName]][[yearString]][[race]]
            standat$pknot_prior_scale = 1/10
            standat$alpha_prior_mean = 2
            standat$alpha_prior_sd = 1
            quiet_fit = quietly(sampling)(
                model_init,
                data = standat,
                chains = 4, cores = 4,
                init_r = 0.5,
                warmup = 5000, iter = 10000,
                open_progress = FALSE,
                control = list(adapt_delta = 0.9999,
                               max_treedepth = 15))

            cat("\n")
            cat("\n")
            cat("Metro: ")
            cat(metroName)
            cat("\n")
            cat("Year: ")
            cat(curr_year)
            cat("\n")
            cat("Race: ")
            cat(race)
            cat("\n")
            cat("\n")
            cat("Warnings: ")
            cat("\n")
            cat("\n")
            cat(quiet_fit$warnings)
            cat("\n")
            cat("\n")
            cat("\n")
            cat("\n")

            gini_fits[[metroName]][[yearString]][[race]] =
                extract(quiet_fit$result,
                        pars = c("pknot", "alpha"))

            rm(quiet_fit)
        }
    }
}
## very few metro areas have a very small number of divergent transitions
save(gini_fits, file = "gini_fits.RData")

load("gini_fits.RData")

## now compute gini from samples
ginis = tibble(
    metro = character(),
    year = numeric(),
    race = character(),
    gini.est = numeric(),
    gini.se = numeric()
)

for(metroName in metros){
    for(curr_year in years){
        yearString = paste("y", curr_year, sep="")
        for(race in races){
            standat = metro_race_income_standats[[metroName]][[yearString]][[race]]
            draws = gini_fits[[metroName]][[yearString]][[race]]
            lb = standat$knots
            ub = c(standat$knots[-1], Inf)
            alpha = draws$alpha
            pknot = draws$pknot
            niter = nrow(pknot)
            nalpha = ncol(alpha)
            nunif = ncol(pknot) - nalpha
            lorenz_unif = (1 + lb / (lb + ub))/3
            mean_unif = (lb + ub) / 2
            lorenz_unif = lorenz_unif[1:nunif]
            mean_unif = mean_unif[1:nunif]
            lb = lb[nunif + 1:nalpha]
            ub = ub[nunif + 1:nalpha]
            gini_sample = sapply(1:niter, function(iter){
                myalpha = alpha[iter,]
                mypknot = pknot[iter,]
                lorenz_pareto = 1 - myalpha / (2 * myalpha - 1) *
                    (1 - (lb/ub) ^ (2*myalpha - 1)) /
                    (1 - (lb/ub) ^ myalpha)
                mean_pareto = myalpha / (myalpha - 1) * lb *
                    (1 - (lb/ub)^(myalpha - 1)) /
                    (1 - (lb/ub)^myalpha)
                lorenz_bin = c(lorenz_unif, lorenz_pareto)
                mean_bin = c(mean_unif, mean_pareto)
                pmu = mypknot * mean_bin
                plorenz = mypknot * lorenz_bin
                psums = c(rev(cumsum(rev(mypknot)))[-1], 0)
                lorenz = sum(pmu * (plorenz + psums)) / sum(pmu)
                1 - 2*lorenz
            })
            ginis = ginis %>%
                add_row(tibble_row(metro = metroName,
                                   year = curr_year,
                                   race = race,
                                   gini.est = mean(gini_sample),
                                   gini.se = sd(gini_sample)))
        }
    }
}

save(ginis, file = "ginis.RData")

rm(gini_fits)
load("ginis.RData")

### next, construct covariate + se matrices for EIV regressions
combined_xest = metro_dataset %>%
    st_drop_geometry %>%
    filter(year == 2018) %>%
    arrange(NAME) %>%
    select(contains(".est"), NAME) %>%
    select(-contains("employ.pop"), -contains("labor.part")) %>%
    rename(Gini = 1, Unemp = 2, AgeOver65 = 3, AgeUnder18 = 4, Edu = 5,
           Income = 6, Foreign = 7, IndustyConstruct = 8, IndustryManuf = 9,
           IndustryFIRE = 10, IndustryProf = 11, FemaleHHer = 12,
           Population = 13, SameHouse = 14, SameCounty = 15,
           NewHouse = 16, Name = 17) %>%
    select(Name, Gini, Population, Unemp, Edu, Income, everything())

combined_xse = metro_dataset %>%
    st_drop_geometry %>%
    filter(year == 2018) %>%
    arrange(NAME) %>%
    select(contains(".se"), NAME) %>%
    select(-contains("employ.pop"), -contains("labor.part")) %>%
    rename(Gini = 1, Unemp = 2, AgeOver65 = 3, AgeUnder18 = 4, Edu = 5,
           Income = 6, Foreign = 7, IndustyConstruct = 8, IndustryManuf = 9,
           IndustryFIRE = 10, IndustryProf = 11, FemaleHHer = 12,
           Population = 13, SameHouse = 14, SameCounty = 15,
           NewHouse = 16, Name = 17) %>%
    select(Name, Gini, Population, Unemp, Edu, Income, everything())

combined_xest_for_race = combined_xest %>%
    select(-Gini, -Population, -Unemp, - Edu, -Income, - FemaleHHer)

combined_xse_for_race = combined_xse %>%
    select(-Gini, -Population, -Unemp, - Edu, -Income, - FemaleHHer)

gini_black = ginis %>%
    filter(year == 2018) %>%
    filter(race == "black") %>%
    mutate(NAME = metro) %>%
    select(NAME, gini.est, gini.se)

black_xest = metro_race_dataset %>%
    st_drop_geometry %>%
    filter(year == 2018) %>%
    select(NAME, contains(".black")) %>%
    full_join(gini_black) %>%
    select(NAME, contains(".est")) %>%
    select(-contains("employ.pop"), -contains("labor.part")) %>%
    rename(Name = 1, Population = 2, Unemp = 3, Edu = 4, Income = 5, Gini = 6) %>%
    select(Name, Gini, everything()) %>%
    full_join(combined_xest_for_race)

black_xse = metro_race_dataset %>%
    st_drop_geometry %>%
    filter(year == 2018) %>%
    select(NAME, contains(".black")) %>%
    full_join(gini_black) %>%
    select(NAME, contains(".se")) %>%
    select(-contains("employ.pop"), -contains("labor.part")) %>%
    rename(Name = 1, Population = 2, Unemp = 3, Edu = 4, Income = 5, Gini = 6) %>%
    select(Name, Gini, everything()) %>%
    full_join(combined_xse_for_race)

gini_white = ginis %>%
    filter(year == 2018) %>%
    filter(race == "white") %>%
    mutate(NAME = metro) %>%
    select(NAME, gini.est, gini.se)

white_xest = metro_race_dataset %>%
    st_drop_geometry %>%
    filter(year == 2018) %>%
    select(NAME, contains(".white")) %>%
    full_join(gini_white) %>%
    select(NAME, contains(".est")) %>%
    select(-contains("employ.pop"), -contains("labor.part")) %>%
    rename(Name = 1, Population = 2, Unemp = 3, Edu = 4, Income = 5, Gini = 6) %>%
    select(Name, Gini, everything()) %>%
    full_join(combined_xest_for_race)

white_xse = metro_race_dataset %>%
    st_drop_geometry %>%
    filter(year == 2018) %>%
    select(NAME, contains(".white")) %>%
    full_join(gini_white) %>%
    select(NAME, contains(".se")) %>%
    select(-contains("employ.pop"), -contains("labor.part")) %>%
    rename(Name = 1, Population = 2, Unemp = 3, Edu = 4, Income = 5, Gini = 6) %>%
    select(Name, Gini, everything()) %>%
    full_join(combined_xse_for_race)

combined_indicies = eiv_indicies %>%
    filter(race == "combined") %>%
    mutate(Name = paste(metro, "Metro Area", sep = " ")) %>%
    select(-metro, -race) %>%
    select(Name, everything())

combined_data = combined_xest %>%
    full_join(combined_xse, by = "Name", suffix = c(".est", ".se")) %>%
    full_join(combined_indicies) %>%
    arrange(Name) %>%
    drop_na() %>%
    filter(ntract_included > 4)


black_indicies = eiv_indicies %>%
    filter(race == "black") %>%
    mutate(Name = paste(metro, "Metro Area", sep = " ")) %>%
    select(-metro, -race) %>%
    select(Name, everything())

black_data = black_xest %>%
    full_join(black_xse, by = "Name", suffix = c(".est", ".se")) %>%
    full_join(black_indicies) %>%
    arrange(Name) %>%
    drop_na() %>%
    filter(ntract_included > 4)


white_indicies = eiv_indicies %>%
    filter(race == "white") %>%
    mutate(Name = paste(metro, "Metro Area", sep = " ")) %>%
    select(-metro, -race) %>%
    select(Name, everything())

white_data = white_xest %>%
    full_join(white_xse, by = "Name", suffix = c(".est", ".se")) %>%
    full_join(white_indicies) %>%
    arrange(Name) %>%
    drop_na() %>%
    filter(ntract_included > 4)



## now beginning running regressions
olsfun = function(standat){
    yy = standat$index
    XX = cbind(1, standat$x)
    XpX = t(XX)%*%XX
    out = chol2inv(chol(XpX))%*%t(XX)%*%yy
    regnames = colnames(standat$x) %>%
        strsplit("\\.est") %>%
        unlist %>%
        c("Intercept", .)
    rownames(out) = regnames
    out
}

## all households regressions
nx = ncol(combined_xest) - 1
combined_standat = list(
    nmetro = nrow(combined_data),
    nx = nx,
    x = combined_data %>%
        select(contains(".est")) %>%
        select(-kl.est, -hr.est) %>%
        as.matrix,
    x_se = combined_data %>%
        select(contains(".se")) %>%
        select(-kl.se, -hr.se) %>%
        as.matrix,
    alpha_cs_prior_mean = 0,
    alpha_cs_prior_sd = 100,
    beta_cs_prior_mean = rep(0, nx),
    beta_cs_prior_sd = rep(3, nx),
    mu_cs_prior_mean = rep(0, nx),
    mu_cs_prior_sd = rep(3, nx),
    tau_cs_prior_mean = 0,
    tau_cs_prior_sd = 0.8,
    sigma_cs_prior_mean = rep(0, nx),
    sigma_cs_prior_sd = rep(2, nx)
)

combined_kl_standat = combined_standat
combined_kl_standat$index = combined_data$kl.est
combined_kl_standat$index_se = combined_data$kl.se

combined_hr_standat = combined_standat
combined_hr_standat$index = combined_data$hr.est
combined_hr_standat$index_se = combined_data$hr.se

eiv_model_init = stan_model(file = "../models/eiv_latent_regression.stan")

## each fit takes a couple minutes
set.seed(9089432)
combined_hr_fit = sampling(eiv_model_init, data = combined_hr_standat,
                           chains = 4, cores = 4, warmup = 2000,
                           iter = 4000, open_progress = FALSE)
combined_hr_summary = summary(combined_hr_fit,
                              pars = c("alpha", "beta"))$summary
combined_hr_cs_summary = summary(combined_hr_fit,
                                 pars = c("alpha_cs", "beta_cs"))$summary
rm(combined_hr_fit)

set.seed(72349345)
combined_kl_fit = sampling(eiv_model_init, data = combined_kl_standat,
                           chains = 4, cores = 4, warmup = 2000,
                           iter = 4000, open_progress = FALSE,
                           control = list(max_treedepth = 12))
combined_kl_summary = summary(combined_kl_fit,
                              pars = c("alpha", "beta"))$summary
combined_kl_cs_summary = summary(combined_kl_fit,
                                 pars = c("alpha_cs", "beta_cs"))$summary
rm(combined_kl_fit)


regnames = colnames(combined_hr_standat$x) %>%
    strsplit("\\.est") %>%
    unlist %>%
    c("Intercept", .)

rownames(combined_hr_summary) = regnames
rownames(combined_hr_cs_summary) = regnames
rownames(combined_kl_summary) = regnames
rownames(combined_kl_cs_summary) = regnames

save(combined_hr_summary, file = "combined_hr_summary.RData")
save(combined_kl_summary, file = "combined_kl_summary.RData")
save(combined_hr_cs_summary, file = "combined_hr_cs_summary.RData")
save(combined_kl_cs_summary, file = "combined_kl_cs_summary.RData")


## black households regressions
nx = ncol(black_xest) - 1
black_standat = list(
    nmetro = nrow(black_data),
    nx = nx,
    x = black_data %>%
        select(contains(".est")) %>%
        select(-kl.est, -hr.est) %>%
        as.matrix,
    x_se = black_data %>%
        select(contains(".se")) %>%
        select(-kl.se, -hr.se) %>%
        as.matrix,
    alpha_cs_prior_mean = 0,
    alpha_cs_prior_sd = 100,
    beta_cs_prior_mean = rep(0, nx),
    beta_cs_prior_sd = rep(3, nx),
    mu_cs_prior_mean = rep(0, nx),
    mu_cs_prior_sd = rep(3, nx),
    tau_cs_prior_mean = 0,
    tau_cs_prior_sd = 0.8,
    sigma_cs_prior_mean = rep(0, nx),
    sigma_cs_prior_sd = rep(2, nx)
)

black_kl_standat = black_standat
black_kl_standat$index = black_data$kl.est
black_kl_standat$index_se = black_data$kl.se

black_hr_standat = black_standat
black_hr_standat$index = black_data$hr.est
black_hr_standat$index_se = black_data$hr.se

eiv_model_init = stan_model(file = "../models/eiv_latent_regression_repar.stan")

## each fit takes about an hour
set.seed(7823489)
black_hr_fit = sampling(eiv_model_init, data = black_hr_standat,
                        chains = 4, cores = 4, warmup = 1000,
                        iter = 12000, open_progress = FALSE,
                        control = list(adapt_delta = 0.9, max_treedepth = 15))
black_hr_summary = summary(black_hr_fit,
                              pars = c("alpha", "beta"))$summary
black_hr_cs_summary = summary(black_hr_fit,
                              pars = c("alpha_cs", "beta_cs"))$summary
rm(black_hr_fit)

set.seed(8792349)
black_kl_fit = sampling(eiv_model_init, data = black_kl_standat,
                        chains = 4, cores = 4, warmup = 1000,
                        iter = 12000, open_progress = FALSE,
                        control = list(adapt_delta = 0.9, max_treedepth = 15))
black_kl_summary = summary(black_kl_fit,
                              pars = c("alpha", "beta"))$summary
black_kl_cs_summary = summary(black_kl_fit,
                              pars = c("alpha_cs", "beta_cs"))$summary
rm(black_kl_fit)

regnames = colnames(black_hr_standat$x) %>%
    strsplit("\\.est") %>%
    unlist %>%
    c("Intercept", .)

rownames(black_hr_summary) = regnames
rownames(black_hr_cs_summary) = regnames
rownames(black_kl_summary) = regnames
rownames(black_kl_cs_summary) = regnames

save(black_hr_summary, file = "black_hr_summary.RData")
save(black_kl_summary, file = "black_kl_summary.RData")
save(black_hr_cs_summary, file = "black_hr_cs_summary.RData")
save(black_kl_cs_summary, file = "black_kl_cs_summary.RData")



## white households regressions
nx = ncol(white_xest) - 1
white_standat = list(
    nmetro = nrow(white_data),
    nx = nx,
    x = white_data %>%
        select(contains(".est")) %>%
        select(-kl.est, -hr.est) %>%
        as.matrix,
    x_se = white_data %>%
        select(contains(".se")) %>%
        select(-kl.se, -hr.se) %>%
        as.matrix,
    alpha_cs_prior_mean = 0,
    alpha_cs_prior_sd = 100,
    beta_cs_prior_mean = rep(0, nx),
    beta_cs_prior_sd = rep(3, nx),
    mu_cs_prior_mean = rep(0, nx),
    mu_cs_prior_sd = rep(3, nx),
    tau_cs_prior_mean = 0,
    tau_cs_prior_sd = 0.8,
    sigma_cs_prior_mean = rep(0, nx),
    sigma_cs_prior_sd = rep(2, nx)
)

white_kl_standat = white_standat
white_kl_standat$index = white_data$kl.est
white_kl_standat$index_se = white_data$kl.se

white_hr_standat = white_standat
white_hr_standat$index = white_data$hr.est
white_hr_standat$index_se = white_data$hr.se

eiv_model_init = stan_model(file = "../models/eiv_latent_regression_repar.stan")

## each fit takes about 2.5-3 hours
set.seed(8923429)
white_hr_fit = sampling(eiv_model_init, data = white_hr_standat,
                        chains = 4, cores = 4, warmup = 2000,
                        iter = 12000, open_progress = FALSE,
                        init_r = 1,
                        control = list(max_treedepth = 16,
                                       adapt_delta = 0.9))
white_hr_summary = summary(white_hr_fit,
                              pars = c("alpha", "beta"))$summary
white_hr_cs_summary = summary(white_hr_fit,
                                 pars = c("alpha_cs", "beta_cs"))$summary
rm(white_hr_fit)

set.seed(4892342)
white_kl_fit = sampling(eiv_model_init, data = white_kl_standat,
                        chains = 4, cores = 4, warmup = 2000,
                        iter = 12000, open_progress = FALSE,
                        init_r = 1,
                        control = list(max_treedepth = 16,
                                       adapt_delta = 0.9))
white_kl_summary = summary(white_kl_fit,
                              pars = c("alpha", "beta"))$summary
white_kl_cs_summary = summary(white_kl_fit,
                                 pars = c("alpha_cs", "beta_cs"))$summary
rm(white_kl_fit)

regnames = colnames(white_hr_standat$x) %>%
    strsplit("\\.est") %>%
    unlist %>%
    c("Intercept", .)

rownames(white_hr_summary) = regnames
rownames(white_hr_cs_summary) = regnames
rownames(white_kl_summary) = regnames
rownames(white_kl_cs_summary) = regnames

save(white_hr_summary, file = "white_hr_summary.RData")
save(white_kl_summary, file = "white_kl_summary.RData")
save(white_hr_cs_summary, file = "white_hr_cs_summary.RData")
save(white_kl_cs_summary, file = "white_kl_cs_summary.RData")



## Create Tables H.1, H.2, and H.3. Table 4 is constructed from these tables.
load("combined_hr_summary.RData")
load("black_hr_summary.RData")
load("white_hr_summary.RData")
digit_matrix = matrix(3, nrow = 17, ncol = 9)
digit_matrix[3,] = -1
digit_matrix[6,] = -1
digit_matrix2 = digit_matrix[-nrow(digit_matrix),]

tabcolnames = c("OLS", "Mean", "SD", "2.5%", "25%", "50%", "75%", "97.5%")

combined_hr_table = cbind(olsfun(combined_hr_standat),
                          combined_hr_summary[,-c(2,9,10)])
colnames(combined_hr_table) = tabcolnames

black_hr_table = cbind(olsfun(black_hr_standat),
                          black_hr_summary[,-c(2,9,10)])
colnames(black_hr_table) = tabcolnames

white_hr_table = cbind(olsfun(white_hr_standat),
                          white_hr_summary[,-c(2,9,10)])
colnames(white_hr_table) = tabcolnames

xtable(combined_hr_table, digits = digit_matrix)
xtable(black_hr_table, digits = digit_matrix2)
xtable(white_hr_table, digits = digit_matrix2)


## Create Tables H.4, H.5, and H.6. Table 5 is constructed from these tables.
load("combined_kl_summary.RData")
load("black_kl_summary.RData")
load("white_kl_summary.RData")
combined_kl_table = cbind(olsfun(combined_kl_standat),
                          combined_kl_summary[,-c(2,9,10)])
colnames(combined_kl_table) = tabcolnames

black_kl_table = cbind(olsfun(black_kl_standat),
                          black_kl_summary[,-c(2,9,10)])
colnames(black_kl_table) = tabcolnames

white_kl_table = cbind(olsfun(white_kl_standat),
                          white_kl_summary[,-c(2,9,10)])
colnames(white_kl_table) = tabcolnames

xtable(combined_kl_table, digits = digit_matrix)
xtable(black_kl_table, digits = digit_matrix2)
xtable(white_kl_table, digits = digit_matrix2)
