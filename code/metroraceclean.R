library(tidyverse)
library(tidycensus)
library(sf)
source('clean_functions.R')

years = c(2013, 2018)

metro_tract_race_dataset = create_metro_tract_race_dataset(years)
save(metro_tract_race_dataset, file = "data/metro_tract_race_dataset.RData")

metro_race_standats = create_metro_race_standats(st_drop_geometry(metro_tract_race_dataset))
save(metro_race_standats, file = "data/metro_race_standats.RData")

metros = get_acs(
    geography = "metropolitan statistical area/micropolitan statistical area",
    table = "B01003",
    cache_table = TRUE,
    year = 2018) %>%    ## using top 100 metros of 2018
    top_n(100, estimate) %>%
    arrange(desc(estimate))

metro_race_dataset = create_metro_race_dataset(years) %>%
    filter(GEOID %in% metros$GEOID)
save(metro_race_dataset, file = "data/metro_race_dataset.RData")

metro_race_income_dataset = create_metro_race_income_dataset(years) %>%
    filter(GEOID %in% metros$GEOID)
save(metro_race_income_dataset, file = "data/metro_race_income_dataset.RData")

metro_race_income_standats =
    create_metro_race_income_standats(years, metro_race_income_dataset)
save(metro_race_income_standats, file = "data/metro_race_income_standats.RData")
