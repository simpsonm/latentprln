## This script creates the synthetic population for the simulation study in Section 4.
## It additionally creates Figure A.2 in Appendix A.
## NOTE: FOR REASONS WE DO NOT UNDERSTAND, IN NEWER VERSIONS OF R, THIS SCRIPT
##       NO LONGER RECREATES OUR ORIGINAL POPULATION.
##       We include the generated file, population.RData, in the GitHub repository
##       for this reason, and have commented out the line which saves the population
##       in this script.
library(tidyverse)
library(tidycensus)
library(sf)
library(sp)
library(maptools)
library(reldist) ## for gini calculation
library(broom)
library(ggthemes)
library(gridExtra)

RNGkind(sample.kind = "Rounding") ## needed for reproducibility in R >= 3.6
set.seed(225613478)

## load the 2010-2014 MO household PUMS and MO tract shapefiles
## then subset both to Boone county and format the data appropriately
mo.pums <- read.csv("../data/MO/ss14hmo.csv") %>% select(WGTP, PUMA00, PUMA10, FINCP, HINCP)
mo.pums2tract <- read.csv("../data/tract2puma2010.txt") %>% filter(STATEFP == 29)
mo.shape <- readShapePoly("../data/MO/shapefile/14000.shp")
boone.tracts <- filter(mo.pums2tract, PUMA5CE == 600)
boone.tracts$TRACT <- ifelse(nchar(boone.tracts$TRACTCE) == 3,
                            paste("000", boone.tracts$TRACTCE, sep=""),
                            paste("00", boone.tracts$TRACTCE, sep=""))
boone.shape <- mo.shape[as.character(mo.shape@data$TRACT) %in% boone.tracts$TRACT &
                        mo.shape@data$COUNTY == "019",]

boone.pops <- get_acs(geography = "tract", table = "S1901", year = 2019,
                      state = "MO", county = "Boone", geometry = TRUE) %>%
    filter(grepl("C01", variable)) %>%
    mutate(variable = recode(variable, "S1901_C01_001" = "Total")) %>%
    filter(variable == "Total") %>%
    mutate(Total = estimate) %>%
    select(-variable, -estimate, -moe)
boone.pops$tract <- substr(boone.pops$GEOID, 6, 12)

boone.shape$TRACT <- as.character(boone.shape$TRACT)
boone.shape@data <- full_join(boone.shape@data, select(boone.pops, tract, Total),
                              c("TRACT" = "tract"))
boone.shape@data$geometry <- NULL
## we will use the observed incomes to create reasonable looking synthetic income distributions
boone.pums <- mo.pums %>% filter(PUMA00 == 600 | PUMA10 == 600) %>% na.omit() %>%
  transmute(weight = WGTP, income = HINCP) %>% filter(weight > 0)

## extract the strata, their weights, and the sample sizes from the PUMS
## and the estimated total population of boone county from the 2010-2014 ACS estimates
nobs.strata.df <- boone.pums %>% group_by(weight) %>% summarise(nobs = n())
wts.unique <- nobs.strata.df$weight
n.strata <- length(wts.unique)
N.pop <- sum(boone.shape@data$Total)

## use stratum weights and sample sizes to create stratum population sizes
pop.probs <- nobs.strata.df$nobs*nobs.strata.df$weight
pop.probs <- pop.probs/sum(pop.probs)
pop.raw <- pop.probs*N.pop
pop.floor <- floor(pop.raw)
pop.diff <- pop.raw - pop.floor
ids <- order(pop.diff, decreasing = TRUE)[1:(N.pop - sum(pop.floor))]
pop.strata <- rep(0, n.strata)
pop.strata[ids] <- 1
pop.strata <- pop.strata + pop.floor
strata.sample.size <- ceiling((nobs.strata.df$nobs / sum(nobs.strata.df$nobs))*(0.1 * N.pop))
ids <- which(strata.sample.size > pop.strata)
strata.sample.size[ids] <- pop.strata[ids] - 1

## consolidate synthetic population info known so far into a list
pop.dat <- list(tract.pops = boone.shape@data$Total, strata.pops = pop.strata,
                stratum.nobs = strata.sample.size,
                wts = pop.strata / strata.sample.size)


## now to assign strata to tracts
## assumes a strong correlation between strata and tracts
n.tract <- length(boone.shape)
tract.stratum.pop.fixed <- matrix(0, n.tract, n.strata)
strata.ids <- 1:n.strata
strata.pops.left <- pop.dat$strata.pops
for(tract in 1:n.tract){
  tract.pop.left <- pop.dat$tract.pops[tract]
  while(tract.pop.left > 0){
    ## if/else because sample() is poorly designed
    sample.space <- strata.ids[strata.pops.left > 0]
    if(length(sample.space) > 1)
      strata.id <- sample(sample.space, 1)
    else
      strata.id <- sample.space
    n.assign <- min(tract.pop.left, strata.pops.left[strata.id])
    tract.stratum.pop.fixed[tract, strata.id] <- tract.stratum.pop.fixed[tract, strata.id] + n.assign
    strata.pops.left[strata.id] <- strata.pops.left[strata.id] - n.assign
    tract.pop.left <- tract.pop.left - n.assign
  }
}

## check to make sure population margins line up for tracts and strata
## (both should evaluate to TRUE)
all.equal(apply(tract.stratum.pop.fixed, 1, sum), pop.dat$tract.pops)
all.equal(apply(tract.stratum.pop.fixed, 2, sum), pop.dat$strata.pops)


## create an in-out gradient for the spatial random effect to generate synthetic incomes
bbox <- boone.shape@bbox
center <- (bbox[,2] - bbox[,1])/2 + bbox[,1]
nsim <- 10000
mean.distance <- rep(0, n.tract)
for(tract in 1:n.tract){
  sims <- spsample(boone.shape@polygons[[tract]], nsim, "regular")
  distance <- apply(sims@coords, 1, function(x){
     sqrt(sum((x - center)^2))
  })
  mean.distance[tract] <- mean(distance)
}
scaled.distance <- (mean.distance - mean(mean.distance))/sd(mean.distance)
scaled.distance[6] <- 2

## add the gradient info to the spatial data frame for plotting
boone.shape@data$Distance <- scaled.distance
boone.geom <- tidy(boone.shape)
centroids <- as.data.frame(coordinates(boone.shape))
names(centroids) <- paste("center", c("x", "y"), sep=".")
IDs <- sapply(boone.shape@polygons, function(ll){ll@ID})
centroids$id <- IDs
centroids$tract <- 1:n.tract
boone.shape@data$id <- IDs
boone.geom.full <- left_join(boone.geom, boone.shape@data, by="id") %>% left_join(centroids)

## plot the spatial random effect's basis function
## basis function has a nice center-out gradient
ggplot(data = boone.geom.full) +
  geom_map(map = boone.geom.full, aes(map_id = id, fill = Distance), color = "black", size = .25) +
  expand_limits(x = boone.geom.full$long, y = boone.geom.full$lat) + coord_equal() + theme_map() +
  scale_fill_gradient(low = "white", high = "blue")

## simulate the population using a two component mixture of lognormals for each
##  tract/stratum combination, using 2010-2014 PUMS incomes as a guideline
set.seed(1465489134)
population <- data.frame(tract = NULL, stratum = NULL, target = NULL)
scaled.log.wts <- (log(wts.unique) - mean(log(wts.unique)))/sd(log(wts.unique))
omega.st <- matrix(0, n.tract, n.strata)      ## mixture probabilities for each tract/stratum combo
mu.st <- array(0, c(n.tract, n.strata, 2))    ## mean parameters for each tract/stratum combo
sigma.st <- array(0, c(n.tract, n.strata, 2)) ## SD parameters for each tract/stratum combo
m.st <- array(0, c(n.tract, n.strata, 2))     ## mean of each lognormal component for each tract/stratum combo
s.st <- array(0, c(n.tract, n.strata, 2))     ## SD of each lognormal component for each tract/stratum combo
## use the average difference between stratum incomes and the overal mean and
##  the SD of incomes within strata to help come up with reasonable parameter
##  values for the lognormal mixture
mu.boone <- wtd.mean(log(boone.pums$income + 1), boone.pums$weight)
sig.boone <- sqrt(wtd.var(log(boone.pums$income + 1), boone.pums$weight))
boone.strata.stats <- data.frame(income = log(boone.pums$income + 1), weight = boone.pums$weight,
                                 stratum = 1:length(boone.pums$income)) %>%
  group_by(weight) %>% summarise(n = n(), mn.diff = mean(income - mu.boone)/(n + 5),
                                 sd = sqrt((mean((income - mean(income))^2)*n + sig.boone^2*500)/(n + 500)))
for(tract in 1:n.tract){
  for(stratum in 1:n.strata){
    ## pick parameter values for the current tract/stratum combo
    omega.st[tract, stratum] <- 1/(1 + exp(0 + .2*scaled.distance[tract] + .2*scaled.log.wts[stratum]))
    mu.st[tract, stratum, 1] <- mu.boone*0.87 +
      (- .3*scaled.distance[tract] + boone.strata.stats$mn.diff[stratum])
    mu.st[tract, stratum, 2] <- mu.boone*1.05 +
      (- .2*scaled.distance[tract] + 1.5*boone.strata.stats$mn.diff[stratum])
    sigma.st[tract, stratum, 1] <- 1*exp( scaled.distance[tract]/5 -
                                            log(boone.strata.stats$sd[stratum]/1)/5)
    sigma.st[tract, stratum, 2] <- 0.6*exp( scaled.distance[tract]/5 -
                                            log(boone.strata.stats$sd[stratum]/0.6)/5)
    n.pop.st <- tract.stratum.pop.fixed[tract, stratum]
    if(n.pop.st > 0){
      ## if the tract/stratum combo has a positive population,
      ## simulate that population from its 2 component mixture of lognormals
      id.st <- 1 + rbinom(n.pop.st, 1, omega.st[tract, stratum])
      sim.st <- exp(mu.st[cbind(tract, stratum, id.st)] +
                       sigma.st[cbind(tract, stratum, id.st)]*rnorm(n.pop.st))
      if(max(sim.st) > 1000000){
        ## if anyone in the tract/stratum combo has a very large income, print info about it
        ## useful for checking how realistic the synthetic distribution looks
        idx <- which(sim.st > 1000000)
        print(sim.st[idx]/1000)
        print(id.st[idx])
        print(mu.st[cbind(tract, stratum, unique(id.st[idx]))])
        print(sigma.st[cbind(tract, stratum, unique(id.st[idx]))])
      }
      population <- rbind(population, data.frame(tract = tract, stratum = stratum, target = sim.st))
    }
  }
}

## collect all of the population info into the population list and save it
population$weight <- pop.dat$wts[population$stratum]
population$stratum.nobs <- pop.dat$stratum.nobs[population$stratum]

## Commented out because this script will not generate the population
##  we used in current versions of R. See note at the top.
## save(population, file = "population.RData")


## create plots to look at various features of the population
load("population.RData")
pop.summary <- population %>% group_by(tract) %>%
  summarize(mean = mean(target), sd = sd(target), median = median(target))
boone.shape@data$median <- pop.summary$median
boone.shape@data$mean <- pop.summary$mean
boone.shape@data$sd <- pop.summary$sd
boone.geom.full.2 <- left_join(boone.geom.full, boone.shape@data, by="id")
p.mean <- ggplot(data = boone.geom.full.2) +
  geom_map(map = boone.geom.full.2, aes(map_id = id, fill = mean),
           color = "black", size = .25) +
  expand_limits(x = boone.geom.full.2$long, y = boone.geom.full.2$lat) + coord_equal() + theme_map() +
  scale_fill_gradient(low = "white", high = "blue") +
  geom_text(aes(label = tract, x = center.x, y = center.y))
p.sd <- ggplot(data = boone.geom.full.2) +
  geom_map(map = boone.geom.full.2, aes(map_id = id, fill = sd), color = "black", size = .25) +
  expand_limits(x = boone.geom.full.2$long, y = boone.geom.full.2$lat) + coord_equal() + theme_map() +
  scale_fill_gradient(low = "white", high = "blue") +
  geom_text(aes(label = tract, x = center.x, y = center.y))
p.median <- ggplot(data = boone.geom.full.2) +
  geom_map(map = boone.geom.full.2, aes(map_id = id, fill = median),
           color = "black", size = .25) +
  expand_limits(x = boone.geom.full.2$long, y = boone.geom.full.2$lat) + coord_equal() + theme_map() +
  scale_fill_gradient(low = "white", high = "blue") +
  geom_text(aes(label = tract, x = center.x, y = center.y))
grid.arrange(p.mean, p.median, p.sd, ncol = 3)

## look at the right tail of the poulation
mean(population$target > 250000)
mean(population$target > 1000000)
sort(population$target, TRUE)[1:10]/1000000

## look at a histogram of the population by tract
ggplot(population) + geom_histogram(aes(target, ..density..), bins = 50) +
  facet_wrap(~tract) + xlim(c(0,250000))

## create density plots of the bins from the bin estimates for the population
bounds.1 <- c(5, 10, 15, 20, 25, 35, 50, 75, 100, 150)
bounds.2 <- c(10, 15, 25, 35, 50, 75, 100, 150, 200)
bounds <- sort(unique(c(bounds.1, bounds.2)))*1000
ggplot(population) + geom_histogram(aes(target, ..density..), breaks = c(0, bounds, 250000)) +
  facet_wrap(~tract) + xlim(c(0,250000))

## look at the histogram of the population by stratum, only for stratums with enough households
ggplot(filter(population, stratum.nobs > 50)) +
  geom_histogram(aes(target, ..density..), bins = 50) + facet_wrap(~stratum) +
  xlim(c(0,250000))


## create Figure A.2 in Appendix A
boone.map.data <- full_join(boone.pops, boone.shape@data, by = c("tract" = "TRACT"))

mean.map <- ggplot(boone.map.data) +
    geom_sf(aes(geometry = geometry, fill = mean/1000),
            color = "black", size = .25) +
    theme_map() +
    scale_fill_continuous(low = "white", high = "darkgreen", limits = c(20, 100)) +
    guides(fill = guide_colourbar(title = "Mean\nIncome\n($1,000)")) +
    ggtitle("Synthetic population\ntract-level mean income") +
    theme(plot.title = element_text(hjust = 0.5), legend.position = c(-0.05,0.05))

median.map <- ggplot(boone.map.data) +
    geom_sf(aes(geometry = geometry, fill = median/1000),
            color = "black", size = .25) +
    theme_map() +
    scale_fill_continuous(low = "white", high = "darkgreen", limits = c(20, 100)) +
    guides(fill = guide_colourbar(title = "Median\nIncome\n($1,000)")) +
    ggtitle("Synthetic population\ntract-level median income") +
    theme(plot.title = element_text(hjust = 0.5), legend.position = c(-0.05,0.05))

sd.map <- ggplot(boone.map.data) +
    geom_sf(aes(geometry = geometry, fill = sd/1000),
            color = "black", size = .25) +
    theme_map() +
    scale_fill_continuous(low = "white", high = "darkgreen", limits = c(60, 80)) +
    guides(fill = guide_colourbar(title = "SD of\nIncome\n($1,000)")) +
    ggtitle("Synthetic population\ntract-level SD of income") +
    theme(plot.title = element_text(hjust = 0.5), legend.position = c(-0.05,0.05))

w <- 3
h <- 5
ggsave("../../doc/mean.png", mean.map, width = w, height = h)
ggsave("../../doc/median.png", median.map, width = w, height = h)
ggsave("../../doc/sd.png", sd.map, width = w, height = h)
