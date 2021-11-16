source('clean_functions.R')

## If you don't already have a census API key, get one and input it as follows
## census_api_key("key", install=TRUE)

states = c("CO", "IL", "MO", "MT", "NY")
years = 2015
pumas = c(821, 3502, 600, 600, 3706)

full_state_datasets = create_state_tract_dataset(states, years)

save(full_state_datasets, file = "data/full_state_datasets.RData")

puma_standats = create_puma_standats(states, pumas, years, full_state_datasets)

save(puma_standats, file = "data/puma_standats.RData")
