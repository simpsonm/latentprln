library(tidyverse)
library(tidycensus)
library(sf)
source('clean_functions.R')

## If you don't already have a census API key, get one and input it as follows
## census_api_key("key", install=TRUE)

years = c(2013, 2018)

metro_tract_dataset = create_metro_tract_dataset(years)
save(metro_tract_dataset, file = "data/metro_tract_dataset.RData")

metro_standats = create_metro_standats(st_drop_geometry(metro_tract_dataset))
save(metro_standats, file = "data/metro_standats.RData")

metros = get_acs(
    geography = "metropolitan statistical area/micropolitan statistical area",
    table = "B01003",
    cache_table = TRUE,
    year = 2018) %>%    ## using top 100 metros of 2018
    top_n(100, estimate) %>%
    arrange(desc(estimate))

metro_dataset = create_metro_dataset(years) %>%
    filter(GEOID %in% metros$GEOID)
save(metro_dataset, file = "data/metro_dataset.RData")


