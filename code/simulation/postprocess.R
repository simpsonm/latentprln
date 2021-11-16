library(data.table)
library(reshape2)
library(rstan)
library(xtable)
source("../postprocess_functions.R")
source("sim_fun.R")
library(tidyverse)

load('prln.out.RData')
load("samples.RData")
load("population.RData")

target.quantiles <- seq(0.05, 0.95, 0.05)
bounds <- c(5, 10, 15, 20, 25, 35, 50, 75, 100, 150, 200)*1000

pop.summary <- population %>% group_by(tract) %>%
    do(pop.summarise(.$target,
                     target.quantiles = target.quantiles,
                     bounds = bounds)) %>%
    ungroup() %>%
    arrange(tract)

pop.summary.mat <- pop.summary %>%
    select(contains('quant'))

load('stat.out1.RData')
stat.out1 <- stat.out
load('stat.out2.RData')
stat.out2 <- stat.out
load('stat.out3.RData')
stat.out3 <- stat.out

stat.out <- rbind(do.call(rbind, stat.out1),
                  do.call(rbind, stat.out2),
                  do.call(rbind, stat.out3))
rm(stat.out1)
rm(stat.out2)
rm(stat.out3)

model.metrics <- eval.model(stat.out, pop.summary)
prln.metrics <- eval.prln(prln.out, pop.summary)
direct.metrics <- eval.direct(samples, pop.summary)

mean.mat <- model.metrics %>% filter(Estimator == 'Mean') %>%
    arrange(Metric) %>%
    select(-Metric, -Estimator) %>%
    as.matrix()

median.mat <- model.metrics %>% filter(Estimator == 'Median') %>%
    arrange(Metric) %>%
    select(-Metric, -Estimator) %>%
    as.matrix()

prln.mat <- prln.metrics %>%
    arrange(Metric) %>%
    select(-Metric, -Estimator) %>%
    as.matrix()

direct.mat <- direct.metrics %>%
    arrange(Metric) %>%
    select(-Metric, -Estimator) %>%
    as.matrix()

mean.rel.mat <- 100 * (mean.mat - direct.mat) / direct.mat
median.rel.mat <- 100 * (median.mat - direct.mat) / direct.mat
prln.rel.mat <- 100 * (prln.mat - direct.mat) / direct.mat

mean.rel.df <- model.metrics %>% filter(Estimator == 'Mean') %>%
    arrange(Metric) %>%
    select(Metric, Estimator) %>%
    bind_cols(data.frame(mean.rel.mat))
median.rel.df <- model.metrics %>% filter(Estimator == 'Median') %>%
    arrange(Metric) %>%
    select(Metric, Estimator) %>%
    bind_cols(data.frame(median.rel.mat))
prln.rel.df <- prln.metrics %>%
    arrange(Metric) %>%
    select(Metric, Estimator) %>%
    bind_cols(data.frame(prln.rel.mat))


rel.df <- bind_rows(mean.rel.df, median.rel.df, prln.rel.df) %>%
    arrange(Metric, Estimator)


## create Table 1
rel.df %>%
    select(Metric, Estimator, P5, P10, P15, P20, P25, P30, P35, P40, P45, P50) %>%
    xtable %>% print(., include.rownames = FALSE)

## create Table 2
rel.df %>%
    select(Metric, Estimator, P55, P60, P65, P70, P75, P80, P85, P90, P95, Gini) %>%
    xtable %>% print(., include.rownames = FALSE)




covertrue <- cover.pop(stat.out, pop.summary) %>%
    mutate(Reference = "Population")

coverprln <- cover.prln(stat.out, prln.out) %>%
    mutate(Reference = "PRLN")

## create Table 3
bind_rows(covertrue, coverprln) %>%
    select(Reference, everything()) %>%
    melt %>%
    dcast(variable ~ Reference) %>%
    mutate(Percentile = variable) %>%
    select(Percentile, Population, PRLN) %>%
    xtable %>%
    print(., include.rownames = FALSE)
