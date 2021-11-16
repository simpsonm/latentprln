library(tidyverse)
library(tidycensus)
library(sf)
options(tigris_use_cache = TRUE)

## var_tab <- load_variables(2015, "acs5", cache = TRUE)
## var_sub <- load_variables(2015, "acs5/subject", cache = TRUE)
## var_tab %>% filter(grepl("B19082", name))
## var_sub %>% filter(grepl("S1901", name))

## ggplot(out_geom, aes(fill = mean.est, geometry = geometry)) +
##     geom_sf(color = NA) +
##     coord_sf(crs = 26911) +
##     scale_fill_viridis_c(option = "magma")

mysign = function(x){
    1*(x > 0) - 1*(x <= 0)
}

countse2propse = function(sx, sn, p, n){
    sqrt(sx^2 - mysign(sx - p*sn) * p^2 * sn^2) / n
}

aggse2meanse = function(sx, sn, m, n){
    sqrt(sx^2 + m^2 * sn^2) / n
}

create_metro_race_income_standats = function(years, metro_race_income_dataset){
    pknot_prior_loc = list()
    for(curr_year in years){
        yearString = paste("y", curr_year, sep="")
        pknot_prior_loc[[yearString]] =
            create_us_race_bin_proportions(curr_year)
    }

    metro_standats = list()

    metros = unique(metro_race_income_dataset$metro)
    races = c("black", "white")

    for(metroName in metros){
        metro_standats[[metroName]] = list()
        for(curr_year in years){
            yearString = paste("y", curr_year, sep="")
            metro_standats[[metroName]][[yearString]] = list()
            for(race in races){
                metro_data = metro_race_income_dataset %>%
                    filter(year == 2018, metro == metroName) %>%
                    select(contains(race))

                ## total population
                total_est = metro_data %>%
                    select(contains("total")) %>%
                    select(contains("est")) %>%
                    as.matrix %>%
                    drop

                total_se = metro_data %>%
                    select(contains("total")) %>%
                    select(contains("se"))
                    as.matrix %>%
                    drop

                metro_standats[[metroName]][[yearString]][[race]]$total_est =
                    total_est

                metro_standats[[metroName]][[yearString]][[race]]$total_se =
                    total_se

                ## bin estimates and knots
                bin_ests = metro_data %>%
                    select(contains("bin")) %>%
                    select(contains("est")) %>%
                    as.matrix %>%
                    drop

                bin_ses = metro_data %>%
                    select(contains("bin")) %>%
                    select(contains("se")) %>%
                    as.matrix %>%
                    drop

                knots = names(bin_ests) %>%
                    str_split("\\.", simplify = TRUE) %>%
                    .[-1,2] %>%
                    as.numeric %>%
                    c(0, .)

                metro_standats[[metroName]][[yearString]][[race]]$nbin =
                    length(knots)
                metro_standats[[metroName]][[yearString]][[race]]$knots =
                    knots
                metro_standats[[metroName]][[yearString]][[race]]$bin_est =
                    bin_ests
                metro_standats[[metroName]][[yearString]][[race]]$bin_se =
                    bin_ses

                ## mean ests
                metro_standats[[metroName]][[yearString]][[race]]$mean_est =
                    metro_data %>%
                    select(contains("mean")) %>%
                    select(contains("est")) %>%
                    as.matrix %>%
                    drop

                metro_standats[[metroName]][[yearString]][[race]]$mean_se =
                    metro_data %>%
                    select(contains("mean")) %>%
                    select(contains("se")) %>%
                    as.matrix %>%
                    drop

                ## median ests
                metro_standats[[metroName]][[yearString]][[race]]$median_est =
                    metro_data %>%
                    select(contains("median")) %>%
                    select(contains("est")) %>%
                    as.matrix %>%
                    drop

                metro_standats[[metroName]][[yearString]][[race]]$median_se =
                    metro_data %>%
                    select(contains("median")) %>%
                    select(contains("se")) %>%
                    as.matrix %>%
                    drop

                ## prior location for pknot_tract
                metro_standats[[metroName]][[yearString]][[race]]$pknot_prior_loc =
                    pknot_prior_loc[[yearString]][[race]]
            }
        }
    }
    metro_standats
}

create_metro_race_income_dataset = function(years){
    metro_geog = "metropolitan statistical area/micropolitan statistical area"
    out_all = NULL

    for(curr_year in years){
        bin_ests_black = get_acs(geography = metro_geog,
                                 table = "B19001B",
                                 cache_table = TRUE,
                                 year = curr_year) %>%
            mutate(variable = recode(variable,
                                     "B19001B_001" = "total",
                                     "B19001B_002" = "bin.below.10",
                                     "B19001B_003" = "bin.10.15",
                                     "B19001B_004" = "bin.15.20",
                                     "B19001B_005" = "bin.20.25",
                                     "B19001B_006" = "bin.25.30",
                                     "B19001B_007" = "bin.30.35",
                                     "B19001B_008" = "bin.35.40",
                                     "B19001B_009" = "bin.40.45",
                                     "B19001B_010" = "bin.45.50",
                                     "B19001B_011" = "bin.50.60",
                                     "B19001B_012" = "bin.60.75",
                                     "B19001B_013" = "bin.75.100",
                                     "B19001B_014" = "bin.100.125",
                                     "B19001B_015" = "bin.125.150",
                                     "B19001B_016" = "bin.150.200",
                                     "B19001B_017" = "bin.200.above")) %>%
            mutate(se = moe / 1.645)

        bin_ests_white = get_acs(geography = metro_geog,
                                 table = "B19001A",
                                 cache_table = TRUE,
                                 year = curr_year) %>%
            mutate(variable = recode(variable,
                                     "B19001A_001" = "total",
                                     "B19001A_002" = "bin.below.10",
                                     "B19001A_003" = "bin.10.15",
                                     "B19001A_004" = "bin.15.20",
                                     "B19001A_005" = "bin.20.25",
                                     "B19001A_006" = "bin.25.30",
                                     "B19001A_007" = "bin.30.35",
                                     "B19001A_008" = "bin.35.40",
                                     "B19001A_009" = "bin.40.45",
                                     "B19001A_010" = "bin.45.50",
                                     "B19001A_011" = "bin.50.60",
                                     "B19001A_012" = "bin.60.75",
                                     "B19001A_013" = "bin.75.100",
                                     "B19001A_014" = "bin.100.125",
                                     "B19001A_015" = "bin.125.150",
                                     "B19001A_016" = "bin.150.200",
                                     "B19001A_017" = "bin.200.above")) %>%
            mutate(se = moe / 1.645)

        median_ests_black = get_acs(geography = metro_geog,
                                    table = "B19013B",
                                    cache_table = TRUE,
                                    year = curr_year) %>%
            mutate(estimate = estimate / 1000, moe = moe / 1000,
                   se = moe / 1.645) %>%
            mutate(variable = recode(variable,
                                     "B19013B_001" = "median"))

        median_ests_white = get_acs(geography = metro_geog,
                                    table = "B19013A",
                                    cache_table = TRUE,
                                    year = curr_year) %>%
            mutate(estimate = estimate / 1000, moe = moe / 1000,
                   se = moe / 1.645) %>%
            mutate(variable = recode(variable,
                                     "B19013A_001" = "median"))

        agg_ests_black = get_acs(geography = metro_geog,
                                 table = "B19013B",
                                 cache_table = TRUE,
                                 year = curr_year) %>%
            mutate(estimate = estimate / 1000, moe = moe / 1000,
                   se = moe / 1.645) %>%
            mutate(variable = recode(variable,
                                     "B19013B_001" = "agg"))

        agg_ests_white = get_acs(geography = metro_geog,
                                 table = "B19013A",
                                 cache_table = TRUE,
                                 year = curr_year) %>%
            mutate(estimate = estimate / 1000, moe = moe / 1000,
                   se = moe / 1.645) %>%
            mutate(variable = recode(variable,
                                     "B19013A_001" = "median"))

        agg_ests_black = get_acs(geography = metro_geog,
                                 table = "B19025B",
                                 cache_table = TRUE,
                                 year = curr_year) %>%
            mutate(variable = recode(variable,
                                     "B19025B_001" = "aggregate")) %>%
            mutate(estimate = estimate / 1000, moe = moe / 1000,
                   se = moe / 1.645)

        agg_ests_white = get_acs(geography = metro_geog,
                                 table = "B19025A",
                                 cache_table = TRUE,
                                 year = curr_year) %>%
            mutate(variable = recode(variable,
                                     "B19025A_001" = "aggregate")) %>%
            mutate(estimate = estimate / 1000, moe = moe / 1000,
                   se = moe / 1.645)

        bin_ests_black_wide = bin_ests_black %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable,
                        values_from = estimate) %>%
            arrange(GEOID) %>%
            mutate_at(vars(matches("bin")), list(~ ./total))

        bin_ses_black_wide = bin_ests_black %>% select(-moe, -estimate) %>%
            pivot_wider(names_from = variable,
                        values_from = se) %>%
            arrange(GEOID) %>%
            mutate(bin.below.10 =
                       countse2propse(bin.below.10, total,
                                      bin_ests_black_wide$bin.below.10,
                                      bin_ests_black_wide$total),
                   bin.10.15 =
                       countse2propse(bin.10.15, total,
                                      bin_ests_black_wide$bin.10.15,
                                      bin_ests_black_wide$total),
                   bin.15.20 =
                       countse2propse(bin.15.20, total,
                                      bin_ests_black_wide$bin.15.20,
                                      bin_ests_black_wide$total),
                   bin.20.25 =
                       countse2propse(bin.20.25, total,
                                      bin_ests_black_wide$bin.20.25,
                                      bin_ests_black_wide$total),
                   bin.25.30 =
                       countse2propse(bin.25.30, total,
                                      bin_ests_black_wide$bin.25.30,
                                      bin_ests_black_wide$total),
                   bin.30.35 =
                       countse2propse(bin.30.35, total,
                                      bin_ests_black_wide$bin.30.35,
                                      bin_ests_black_wide$total),
                   bin.35.40 =
                       countse2propse(bin.35.40, total,
                                      bin_ests_black_wide$bin.35.40,
                                      bin_ests_black_wide$total),
                   bin.40.45 =
                       countse2propse(bin.40.45, total,
                                      bin_ests_black_wide$bin.40.45,
                                      bin_ests_black_wide$total),
                   bin.45.50 =
                       countse2propse(bin.45.50, total,
                                      bin_ests_black_wide$bin.45.50,
                                      bin_ests_black_wide$total),
                   bin.50.60 =
                       countse2propse(bin.50.60, total,
                                      bin_ests_black_wide$bin.50.60,
                                      bin_ests_black_wide$total),
                   bin.60.75 =
                       countse2propse(bin.60.75, total,
                                      bin_ests_black_wide$bin.60.75,
                                      bin_ests_black_wide$total),
                   bin.75.100 =
                       countse2propse(bin.75.100, total,
                                      bin_ests_black_wide$bin.75.100,
                                      bin_ests_black_wide$total),
                   bin.100.125 =
                       countse2propse(bin.100.125, total,
                                      bin_ests_black_wide$bin.100.125,
                                      bin_ests_black_wide$total),
                   bin.125.150 =
                       countse2propse(bin.125.150, total,
                                      bin_ests_black_wide$bin.125.150,
                                      bin_ests_black_wide$total),
                   bin.150.200 =
                       countse2propse(bin.150.200, total,
                                      bin_ests_black_wide$bin.150.200,
                                      bin_ests_black_wide$total),
                   bin.200.above =
                       countse2propse(bin.200.above, total,
                                      bin_ests_black_wide$bin.200.above,
                                      bin_ests_black_wide$total))

        bin_ests_white_wide = bin_ests_white %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable,
                        values_from = estimate) %>%
            arrange(GEOID) %>%
            mutate_at(vars(matches("bin")), list(~ ./total))

        bin_ses_white_wide = bin_ests_white %>% select(-moe, -estimate) %>%
            pivot_wider(names_from = variable,
                        values_from = se) %>%
            arrange(GEOID) %>%
            mutate(bin.below.10 =
                       countse2propse(bin.below.10, total,
                                      bin_ests_white_wide$bin.below.10,
                                      bin_ests_white_wide$total),
                   bin.10.15 =
                       countse2propse(bin.10.15, total,
                                      bin_ests_white_wide$bin.10.15,
                                      bin_ests_white_wide$total),
                   bin.15.20 =
                       countse2propse(bin.15.20, total,
                                      bin_ests_white_wide$bin.15.20,
                                      bin_ests_white_wide$total),
                   bin.20.25 =
                       countse2propse(bin.20.25, total,
                                      bin_ests_white_wide$bin.20.25,
                                      bin_ests_white_wide$total),
                   bin.25.30 =
                       countse2propse(bin.25.30, total,
                                      bin_ests_white_wide$bin.25.30,
                                      bin_ests_white_wide$total),
                   bin.30.35 =
                       countse2propse(bin.30.35, total,
                                      bin_ests_white_wide$bin.30.35,
                                      bin_ests_white_wide$total),
                   bin.35.40 =
                       countse2propse(bin.35.40, total,
                                      bin_ests_white_wide$bin.35.40,
                                      bin_ests_white_wide$total),
                   bin.40.45 =
                       countse2propse(bin.40.45, total,
                                      bin_ests_white_wide$bin.40.45,
                                      bin_ests_white_wide$total),
                   bin.45.50 =
                       countse2propse(bin.45.50, total,
                                      bin_ests_white_wide$bin.45.50,
                                      bin_ests_white_wide$total),
                   bin.50.60 =
                       countse2propse(bin.50.60, total,
                                      bin_ests_white_wide$bin.50.60,
                                      bin_ests_white_wide$total),
                   bin.60.75 =
                       countse2propse(bin.60.75, total,
                                      bin_ests_white_wide$bin.60.75,
                                      bin_ests_white_wide$total),
                   bin.75.100 =
                       countse2propse(bin.75.100, total,
                                      bin_ests_white_wide$bin.75.100,
                                      bin_ests_white_wide$total),
                   bin.100.125 =
                       countse2propse(bin.100.125, total,
                                      bin_ests_white_wide$bin.100.125,
                                      bin_ests_white_wide$total),
                   bin.125.150 =
                       countse2propse(bin.125.150, total,
                                      bin_ests_white_wide$bin.125.150,
                                      bin_ests_white_wide$total),
                   bin.150.200 =
                       countse2propse(bin.150.200, total,
                                      bin_ests_white_wide$bin.150.200,
                                      bin_ests_white_wide$total),
                   bin.200.above =
                       countse2propse(bin.200.above, total,
                                      bin_ests_white_wide$bin.200.above,
                                      bin_ests_white_wide$total))

        median_ests_white_wide = median_ests_white %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable,
                        values_from = estimate)

        median_ests_black_wide = median_ests_black %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable,
                        values_from = estimate)

        median_ses_white_wide = median_ests_white %>% select(-moe, -estimate) %>%
            pivot_wider(names_from = variable,
                        values_from = se)

        median_ses_black_wide = median_ests_black %>% select(-moe, -estimate) %>%
            pivot_wider(names_from = variable,
                        values_from = se)

        mean_ests_white_wide = agg_ests_white %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable,
                        values_from = estimate) %>%
            full_join(bin_ests_white_wide) %>%
            arrange(GEOID) %>%
            transmute(GEOID = GEOID,
                      NAME = NAME,
                      total = total,
                      mean = aggregate / total)

        mean_ests_black_wide = agg_ests_black %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable,
                        values_from = estimate) %>%
            full_join(bin_ests_black_wide) %>%
            arrange(GEOID) %>%
            transmute(GEOID = GEOID,
                      NAME = NAME,
                      total = total,
                      mean = aggregate / total)

        mean_ses_white_wide = agg_ests_white %>% select(-moe, -estimate) %>%
            pivot_wider(names_from = variable,
                        values_from = se) %>%
            full_join(bin_ses_black_wide) %>%
            arrange(GEOID) %>%
            transmute(GEOID = GEOID,
                      NAME = NAME,
                      mean = aggse2meanse(aggregate,
                                          total,
                                          mean_ests_white_wide$mean,
                                          mean_ests_white_wide$total))

        mean_ses_black_wide = agg_ests_black %>% select(-moe, -estimate) %>%
            pivot_wider(names_from = variable,
                        values_from = se) %>%
            full_join(bin_ses_black_wide) %>%
            arrange(GEOID) %>%
            transmute(GEOID = GEOID,
                      NAME = NAME,
                      mean = aggse2meanse(aggregate,
                                          total,
                                          mean_ests_black_wide$mean,
                                          mean_ests_black_wide$total))


        full_ests_white = bin_ests_white_wide %>%
            full_join(mean_ests_white_wide) %>%
            full_join(median_ests_white_wide) %>%
            select(GEOID, NAME, total, mean, median, contains("bin"))

        full_ests_black = bin_ests_black_wide %>%
            full_join(mean_ests_black_wide) %>%
            full_join(median_ests_black_wide) %>%
            select(GEOID, NAME, total, mean, median, contains("bin"))

        full_ses_white = bin_ses_white_wide %>%
            full_join(mean_ses_white_wide) %>%
            full_join(median_ses_white_wide) %>%
            select(GEOID, NAME, total, mean, median, contains("bin"))

        full_ses_black = bin_ses_black_wide %>%
            full_join(mean_ses_black_wide) %>%
            full_join(median_ses_black_wide) %>%
            select(GEOID, NAME, total, mean, median, contains("bin"))


        full_ests = full_join(full_ests_white, full_ests_black,
                              by = c('GEOID', 'NAME'),
                              suffix = c(".white", ".black")) %>%
            select(GEOID, NAME, contains("white"), contains("black"))

        full_ses = full_join(full_ses_white, full_ses_black,
                             by = c('GEOID', 'NAME'),
                             suffix = c(".white", ".black")) %>%
            select(GEOID, NAME, contains("white"), contains("black"))

        out = full_join(full_ests, full_ses,
                        by = c('GEOID', 'NAME'), suffix = c(".est", ".se")) %>%
            mutate(metro = NAME, year = curr_year) %>%
            select(GEOID, NAME, metro, year, everything())

        out_all = rbind(out_all, out)
    }
    out_all
}


create_metro_race_dataset = function(years){
    metro_geog = "metropolitan statistical area/micropolitan statistical area"
    out_all = NULL

    for(curr_year in years){
        agg_ests_white = get_acs(geography = metro_geog,
                                 table = "B19025A",
                                 cache_table = TRUE,
                                 year = curr_year) %>%
            mutate(variable = recode(variable,
                                     "B19025A_001" = "aggregate")) %>%
            mutate(estimate = estimate / 1000, moe = moe / 1000,
                   se = moe / 1.645)

        agg_ests_black = get_acs(geography = metro_geog,
                                 table = "B19025B",
                                 cache_table = TRUE,
                                 year = curr_year) %>%
            mutate(variable = recode(variable,
                                     "B19025B_001" = "aggregate")) %>%
            mutate(estimate = estimate / 1000, moe = moe / 1000,
                   se = moe / 1.645)

        race_pop_ests = get_acs(geography = metro_geog,
                                table = "B02001",
                                cache_table = TRUE,
                                year = curr_year) %>%
            mutate(variable = recode(variable,
                                     "B02001_001" = "pop",
                                     "B02001_002" = "pop.white",
                                     "B02001_003" = "pop.black",
                                     "B02001_004" = "pop.native",
                                     "B02001_005" = "pop.asian",
                                     "B02001_006" = "pop.islander",
                                     "B02001_007" = "pop.other",
                                     "B02001_008" = "pop.twoplus")) %>%
            filter(!grepl("B02001", variable)) %>%
            mutate(se = moe / 1.645)

        unemp_ests_new = get_acs(geography = metro_geog,
                                 table = "S2301",
                                 cache_table = TRUE,
                                 year = curr_year) %>%
            mutate(variable = recode(variable,
                                     "S2301_C02_012" = "labor.part.rate.white",
                                     "S2301_C02_013" = "labor.part.rate.black",
                                     "S2301_C03_012" = "employ.pop.ratio.white",
                                     "S2301_C03_013" = "employ.pop.ratio.black",
                                     "S2301_C04_012" = "unemploy.rate.white",
                                     "S2301_C04_013" = "unemploy.rate.black"
                                     )) %>%
            filter(!grepl("S2301", variable)) %>%
            mutate(estimate = estimate / 100, moe = moe / 100,
                   se = moe / 1.645)

        edu_ests_new = get_acs(geography = metro_geog,
                               table = "S1501",
                               cache_table = TRUE,
                               year = curr_year) %>%
            mutate(variable = recode(variable,
                                     "S1501_C02_029" =
                                         "edu.25up.hs.plus.white",
                                     "S1501_C02_035" =
                                         "edu.25up.hs.plus.black")) %>%
            filter(!grepl("S1501", variable)) %>%
            mutate(estimate = estimate / 100, moe = moe / 100,
                   se = moe / 1.645)

        femalehead_ests_new = get_acs(geography = metro_geog,
                                      table = "B09019",
                                      cache_table = TRUE,
                                      year = curr_year) %>%
            mutate(variable = recode(variable,
                                     "B09019_005" =
                                         "total.family.male.householder",
                                     "B09019_006" =
                                         "total.family.female.householder"
                                     )) %>%
            filter(!grepl("B09019", variable)) %>%
            mutate(se = moe / 1.645)

        race_pop_ests_wide = race_pop_ests %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable, values_from = estimate)

        race_pop_ses_wide = race_pop_ests %>% select(-moe, -estimate) %>%
            pivot_wider(names_from = variable, values_from = se)

        unemp_ests_wide = unemp_ests_new %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable, values_from = estimate)
        unemp_ses_wide = unemp_ests_new %>% select(-estimate, -moe) %>%
            pivot_wider(names_from = variable, values_from = se)

        edu_ests_wide = edu_ests_new %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable, values_from = estimate)
        edu_ses_wide = edu_ests_new %>% select(-estimate, -moe) %>%
            pivot_wider(names_from = variable, values_from = se)

        percapita_ests_white_wide = agg_ests_white %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable,
                        values_from = estimate) %>%
            full_join(race_pop_ests_wide) %>%
            arrange(GEOID) %>%
            transmute(GEOID = GEOID,
                      NAME = NAME,
                      pop = pop.white,
                      percapita.income = aggregate / pop.white)
        percapita_ses_white_wide = agg_ests_white %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable,
                        values_from = estimate) %>%
            full_join(race_pop_ses_wide) %>%
            arrange(GEOID) %>%
            transmute(GEOID = GEOID,
                      NAME = NAME,
                      pop = pop.white,
                      percapita.income =
                          aggse2meanse(aggregate, pop,
                                       percapita_ests_white_wide$percapita.income,
                                       percapita_ests_white_wide$pop))
        percapita_ests_white_wide = percapita_ests_white_wide %>%
            select(-pop)
        percapita_ests_white_wide = percapita_ses_white_wide %>%
            select(-pop)

        percapita_ests_black_wide = agg_ests_black %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable,
                        values_from = estimate) %>%
            full_join(race_pop_ests_wide) %>%
            arrange(GEOID) %>%
            transmute(GEOID = GEOID,
                      NAME = NAME,
                      pop = pop.black,
                      percapita.income = aggregate / pop.black)
        percapita_ses_black_wide = agg_ests_black %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable,
                        values_from = estimate) %>%
            full_join(race_pop_ses_wide) %>%
            arrange(GEOID) %>%
            transmute(GEOID = GEOID,
                      NAME = NAME,
                      pop = pop.black,
                      percapita.income =
                          aggse2meanse(aggregate, pop,
                                       percapita_ests_black_wide$percapita.income,
                                       percapita_ests_black_wide$pop))
        percapita_ests_black_wide = percapita_ests_black_wide %>%
            select(-pop)
        percapita_ests_black_wide = percapita_ses_black_wide %>%
            select(-pop)

        full_ests_white = race_pop_ests_wide %>%
            select(GEOID, NAME, pop.white) %>%
            full_join(unemp_ests_wide) %>%
            full_join(edu_ests_wide) %>%
            full_join(percapita_ests_white_wide) %>%
            select(-contains("black")) %>%
            select(GEOID, NAME, everything())

        full_ests_black = race_pop_ests_wide %>%
            select(GEOID, NAME, pop.black) %>%
            full_join(unemp_ests_wide) %>%
            full_join(edu_ests_wide) %>%
            full_join(percapita_ests_black_wide) %>%
            select(-contains("white")) %>%
            select(GEOID, NAME, everything())

        full_ses_white = race_pop_ses_wide %>%
            select(GEOID, NAME, pop.white) %>%
            full_join(unemp_ses_wide) %>%
            full_join(edu_ses_wide) %>%
            full_join(percapita_ses_white_wide) %>%
            select(-contains("black")) %>%
            select(GEOID, NAME, everything())

        full_ses_black = race_pop_ses_wide %>%
            select(GEOID, NAME, pop.black) %>%
            full_join(unemp_ses_wide) %>%
            full_join(edu_ses_wide) %>%
            full_join(percapita_ses_black_wide) %>%
            select(-contains("white")) %>%
            select(GEOID, NAME, everything())

        full_ests = full_join(full_ests_white, full_ests_black,
                              by = c('GEOID', 'NAME'),
                              suffix = c(".white", ".black")) %>%
            select(GEOID, NAME, contains("white"), contains("black"))

        full_ses = full_join(full_ses_white, full_ses_black,
                             by = c('GEOID', 'NAME'),
                             suffix = c(".white", ".black")) %>%
            select(GEOID, NAME, contains("white"), contains("black"))

        out = full_join(full_ests, full_ses,
                        by = c('GEOID', 'NAME'), suffix = c(".est", ".se"))

        geom = get_acs(geography = metro_geog,
                       variables = "B19013_001",
                       cache_table = TRUE,
                       year = curr_year,
                       geometry = TRUE) %>%
            select(-estimate, -moe, -variable)

        out_geom = full_join(geom, out) %>%
            mutate(year = curr_year) %>%
            select(GEOID, NAME, year, everything())

        out_all = rbind(out_all, out_geom)
    }
    out_all
}

create_metro_dataset = function(years){
    metro_geog = "metropolitan statistical area/micropolitan statistical area"
    out_all = NULL

    for(curr_year in years){
        gini_ests_new = get_acs(geography = metro_geog,
                                table = "B19083",
                                cache_table = TRUE,
                                year = curr_year) %>%
            mutate(variable = recode(variable, "B19083_001" = "gini")) %>%
            mutate(se = moe / 1.645)
        unemp_ests_new = get_acs(geography = metro_geog,
                                 table = "S2301",
                                 cache_table = TRUE,
                                 year = curr_year) %>%
            mutate(variable = recode(variable,
                                     "S2301_C02_001" = "labor.part.rate",
                                     "S2301_C03_001" = "employ.pop.ratio",
                                     "S2301_C04_001" = "unemploy.rate")) %>%
            filter(!grepl("S2301", variable)) %>%
            mutate(estimate = estimate / 100, moe = moe / 100,
                   se = moe / 1.645)

        age_ests_new = get_acs(geography = metro_geog,
                               table = "S0101",
                               cache_table = TRUE,
                               year = curr_year) %>%
            mutate(variable = recode(variable,
                                     "S0101_C02_026" = "age.18.over",
                                     "S0101_C02_030" = "age.65.over")) %>%
            filter(!grepl("S0101", variable)) %>%
            mutate(estimate = estimate / 100, moe = moe / 100,
                   se = moe / 1.645)

        edu_ests_new = get_acs(geography = metro_geog,
                               table = "S1501",
                               cache_table = TRUE,
                               year = curr_year) %>%
            mutate(variable = recode(variable,
                                     "S1501_C02_014" = "edu.25up.hs.plus")) %>%
            filter(!grepl("S1501", variable)) %>%
            mutate(estimate = estimate / 100, moe = moe / 100,
                   se = moe / 1.645)

        percapita_ests_new = get_acs(geography = metro_geog,
                                     table = "B19301",
                                     cache_table = TRUE,
                                     year = curr_year) %>%
            mutate(variable = recode(variable,
                                     "B19301_001" = "per.capita.income")) %>%
            mutate(estimate = estimate / 1000, moe = moe / 1000,
                   se = moe / 1.645)

        foreign_ests_new = get_acs(geography = metro_geog,
                                   table = "DP02",
                                   cache_table = TRUE,
                                   year = curr_year) %>%
            mutate(variable = recode(variable,
                                     "DP02_0092P" = "foreign.born")) %>%
            filter(!grepl("DP02", variable)) %>%
            mutate(estimate = estimate / 100, moe = moe / 100,
                   se = moe / 1.645)

        occu_ests_new = get_acs(geography = metro_geog,
                                table = "S2405",
                                cache_table = TRUE,
                                year = curr_year) %>%
            mutate(variable = recode(variable,
                                     "S2405_C01_001" = "total.civilian",
                                     "S2405_C01_003" = "total.construction",
                                     "S2405_C01_004" = "total.manufacturing",
                                     "S2405_C01_009" = "total.FIRE",
                                     "S2405_C01_008" = "total.information",
                                     "S2405_C01_010" = "total.professional",
                                     "S2405_C01_011" = "total.edu.health",
                                     "S2405_C01_014" = "total.public.admin")) %>%
            filter(!grepl("S2405", variable)) %>%
            mutate(se = moe / 1.645)

        femalehead_ests_new = get_acs(geography = metro_geog,
                                      table = "B09019",
                                      cache_table = TRUE,
                                      year = curr_year) %>%
            mutate(variable = recode(variable,
                                     "B09019_005" =
                                         "total.family.male.householder",
                                     "B09019_006" =
                                         "total.family.female.householder"
                                     )) %>%
            filter(!grepl("B09019", variable)) %>%
            mutate(se = moe / 1.645)

        mobility_ests_new = get_acs(geography = metro_geog,
                                    table = "B07204",
                                    cache_table = TRUE,
                                    year = curr_year) %>%
            mutate(variable = recode(variable,
                                     "B07204_001" = "pop",
                                     "B07204_002" = "pop.same.house",
                                     "B07204_005" =
                                         "pop.same.town.same.county",
                                     "B07204_008" =
                                         "pop.diff.town.same.county"
                                     )) %>%
            filter(!grepl("B07204", variable)) %>%
            mutate(se = moe / 1.645)

        housing_ests_new = get_acs(geography = metro_geog,
                                   table = "DP04",
                                   cache_table = TRUE,
                                   year = curr_year) %>%
            mutate(variable = recode(variable,
                                     "DP04_0017P" = "house.built.less.5",
                                     )) %>%
            filter(!grepl("DP04", variable)) %>%
            mutate(estimate = estimate / 100, moe = moe / 100,
                   se = moe / 1.645)

        gini_ests_wide = gini_ests_new %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable, values_from = estimate)
        gini_ses_wide = gini_ests_new %>% select(-estimate, -moe) %>%
            pivot_wider(names_from = variable, values_from = se)

        unemp_ests_wide = unemp_ests_new %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable, values_from = estimate)
        unemp_ses_wide = unemp_ests_new %>% select(-estimate, -moe) %>%
            pivot_wider(names_from = variable, values_from = se)

        age_ests_wide = age_ests_new %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable, values_from = estimate) %>%
            mutate(age.under.18 = 1 - age.18.over) %>%
            select(-age.18.over)
        age_ses_wide = age_ests_new %>% select(-estimate, -moe) %>%
            pivot_wider(names_from = variable, values_from = se) %>%
            mutate(age.under.18 = age.18.over) %>%
            select(-age.18.over)

        edu_ests_wide = edu_ests_new %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable, values_from = estimate)
        edu_ses_wide = edu_ests_new %>% select(-estimate, -moe) %>%
            pivot_wider(names_from = variable, values_from = se)

        percapita_ests_wide = percapita_ests_new %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable, values_from = estimate)
        percapita_ses_wide = percapita_ests_new %>% select(-estimate, -moe) %>%
            pivot_wider(names_from = variable, values_from = se)

        foreign_ests_wide = foreign_ests_new %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable, values_from = estimate)
        foreign_ses_wide = foreign_ests_new %>% select(-estimate, -moe) %>%
            pivot_wider(names_from = variable, values_from = se)

        occu_ests_wide = occu_ests_new %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable, values_from = estimate) %>%
            arrange(GEOID) %>%
            transmute(
                GEOID = GEOID,
                NAME = NAME,
                total = total.civilian,
                construction = total.construction / total.civilian,
                manufacturing = total.manufacturing / total.civilian,
                fire = total.FIRE / total.civilian,
                professional =
                    (total.FIRE +
                     total.information +
                     total.professional +
                     total.edu.health +
                     total.public.admin) / total.civilian)
        occu_ses_wide = occu_ests_new %>% select(-estimate, -moe) %>%
            pivot_wider(names_from = variable, values_from = se) %>%
            arrange(GEOID) %>%
            transmute(
                GEOID = GEOID,
                NAME = NAME,
                construction =
                    countse2propse(total.construction,
                                   total.civilian,
                                   occu_ests_wide$construction,
                                   occu_ests_wide$total),

                manufacturing =
                    countse2propse(total.manufacturing,
                                   total.civilian,
                                   occu_ests_wide$manufacturing,
                                   occu_ests_wide$total),
                fire =
                    countse2propse(total.FIRE,
                                   total.civilian,
                                   occu_ests_wide$fire,
                                   occu_ests_wide$total),
                professional =
                    countse2propse(
                        sqrt(total.FIRE^2 + total.information^2 +
                             total.professional^2 + total.edu.health^2),
                        total.civilian,
                        occu_ests_wide$professional,
                        occu_ests_wide$total))
        occu_ests_wide = occu_ests_wide %>% select(-total)

        femalehead_ests_wide = femalehead_ests_new %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable, values_from = estimate) %>%
            arrange(GEOID) %>%
            transmute(GEOID = GEOID,
                      NAME = NAME,
                      total = total.family.male.householder +
                          total.family.female.householder,
                      female.hher = total.family.female.householder / total)

        femalehead_ses_wide = femalehead_ests_new %>% select(-estimate, -moe) %>%
            pivot_wider(names_from = variable, values_from = se) %>%
            arrange(GEOID) %>%
            transmute(GEOID = GEOID,
                      NAME = NAME,
                      female.hher =
                          countse2propse(
                              total.family.female.householder,
                              sqrt(total.family.female.householder^2 +
                                   total.family.male.householder^2),
                              femalehead_ests_wide$female.hher,
                              femalehead_ests_wide$total))
        femalehead_ests_wide = femalehead_ests_wide %>% select(-total)

        mobility_ests_wide = mobility_ests_new %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable, values_from = estimate) %>%
            arrange(GEOID) %>%
            transmute(GEOID = GEOID,
                      NAME = NAME,
                      pop = pop,
                      same.house = pop.same.house / pop,
                      same.county =
                          (pop.same.town.same.county +
                           pop.diff.town.same.county) / pop)
        mobility_ses_wide = mobility_ests_new %>% select(-estimate, -moe) %>%
            pivot_wider(names_from = variable, values_from = se) %>%
            arrange(GEOID) %>%
            transmute(GEOID = GEOID,
                      NAME = NAME,
                      pop = pop,
                      same.house =
                          countse2propse(
                              pop.same.house,
                              pop,
                              mobility_ests_wide$same.house,
                              mobility_ests_wide$pop),
                      same.county =
                          countse2propse(
                              sqrt(pop.same.town.same.county^2 +
                                   pop.diff.town.same.county^2),
                              pop,
                              mobility_ests_wide$same.county,
                              mobility_ests_wide$pop))

        housing_ests_wide = housing_ests_new %>% select(-moe, -se) %>%
            pivot_wider(names_from = variable, values_from = estimate)
        housing_ses_wide = housing_ests_new %>% select(-estimate, -moe) %>%
            pivot_wider(names_from = variable, values_from = se)

        full_ests = full_join(gini_ests_wide, unemp_ests_wide) %>%
            full_join(age_ests_wide) %>%
            full_join(edu_ests_wide) %>%
            full_join(percapita_ests_wide) %>%
            full_join(foreign_ests_wide) %>%
            full_join(occu_ests_wide) %>%
            full_join(femalehead_ests_wide) %>%
            full_join(mobility_ests_wide) %>%
            full_join(housing_ests_wide) %>%
            select(GEOID, NAME, gini, everything())

        full_ses = full_join(gini_ses_wide, unemp_ses_wide) %>%
            full_join(age_ses_wide) %>%
            full_join(edu_ses_wide) %>%
            full_join(percapita_ses_wide) %>%
            full_join(foreign_ses_wide) %>%
            full_join(occu_ses_wide) %>%
            full_join(femalehead_ses_wide) %>%
            full_join(mobility_ses_wide) %>%
            full_join(housing_ses_wide) %>%
            select(GEOID, NAME, gini, everything())

        out = full_join(full_ests, full_ses,
                        by = c('GEOID', 'NAME'), suffix = c(".est", ".se"))

        geom = get_acs(geography = metro_geog,
                       variables = "B19013_001",
                       cache_table = TRUE,
                       year = curr_year,
                       geometry = TRUE) %>%
            select(-estimate, -moe, -variable)

        out_geom = full_join(geom, out) %>%
            mutate(year = curr_year) %>%
            select(GEOID, NAME, year, everything())

        out_all = rbind(out_all, out_geom)

    }
    out_all
}


create_metro_tract_race_dataset = function(years){

    pums2tract = read_csv("data/tract2puma2010.txt")
    out_all = NULL

    metros = get_acs(
        geography = "metropolitan statistical area/micropolitan statistical area",
        table = "B01003",
        cache_table = TRUE,
        year = 2018) %>%    ## using top 100 metros of 2018
        top_n(100, estimate) %>%
        arrange(desc(estimate))

    metro_deliniation = read_csv("data/metro_deliniation_4-2018.csv") %>%
        filter(as.character(`CBSA Code`) %in% metros$GEOID) %>%
        mutate(GEOID = as.character(`CBSA Code`),
               MetroName = `CBSA Title`,
               StateName = `State Name`,
               County = `County/County Equivalent`,
               StateFips = `FIPS State Code`,
               CountyFips = `FIPS County Code`) %>%
        select(GEOID, MetroName, StateName, County, StateFips, CountyFips)

    for(curr_metro in metros$GEOID){
        metroName = metro_deliniation %>%
            filter(GEOID == curr_metro) %>%
            .[['MetroName']] %>%
            unique

        stateFips = metro_deliniation %>%
            filter(GEOID == curr_metro) %>%
            .[['StateFips']] %>%
            unique

        states = fips_codes %>%
            filter(state_code %in% stateFips) %>%
            .[['state']] %>%
            unique

        for(curr_state in states){

            fips = fips_codes %>%
                filter(state == curr_state) %>%
                select(state_code) %>% .[1,1]

            countyFips = metro_deliniation %>%
                filter(GEOID == curr_metro, StateFips == fips) %>%
                .[['CountyFips']] %>%
                unique

            tractFips = pums2tract %>%
                filter(STATEFP == fips, COUNTYFP %in% countyFips) %>%
                mutate(TractFips = paste(fips, COUNTYFP, TRACTCE, sep = "")) %>%
                .[['TractFips']] %>%
                unique

            for(curr_year in years){

                bin_ests_white = get_acs(geography = "tract",
                                         table = "B19001A",
                                         state = curr_state,
                                         county = countyFips,
                                         cache_table = TRUE,
                                         year = curr_year) %>%
                    mutate(variable = recode(variable,
                                             "B19001A_001" = "total.white",
                                             "B19001A_002" = "bin.below.10",
                                             "B19001A_003" = "bin.10.15",
                                             "B19001A_004" = "bin.15.20",
                                             "B19001A_005" = "bin.20.25",
                                             "B19001A_006" = "bin.25.30",
                                             "B19001A_007" = "bin.30.35",
                                             "B19001A_008" = "bin.35.40",
                                             "B19001A_009" = "bin.40.45",
                                             "B19001A_010" = "bin.45.50",
                                             "B19001A_011" = "bin.50.60",
                                             "B19001A_012" = "bin.60.75",
                                             "B19001A_013" = "bin.75.100",
                                             "B19001A_014" = "bin.100.125",
                                             "B19001A_015" = "bin.125.150",
                                             "B19001A_016" = "bin.150.200",
                                             "B19001A_017" = "bin.200.above")) %>%
                    mutate(se = moe / 1.645)

                bin_ests_black = get_acs(geography = "tract",
                                         table = "B19001B",
                                         state = curr_state,
                                         county = countyFips,
                                         cache_table = TRUE,
                                         year = curr_year) %>%
                    mutate(variable = recode(variable,
                                             "B19001B_001" = "total.black",
                                             "B19001B_002" = "bin.below.10",
                                             "B19001B_003" = "bin.10.15",
                                             "B19001B_004" = "bin.15.20",
                                             "B19001B_005" = "bin.20.25",
                                             "B19001B_006" = "bin.25.30",
                                             "B19001B_007" = "bin.30.35",
                                             "B19001B_008" = "bin.35.40",
                                             "B19001B_009" = "bin.40.45",
                                             "B19001B_010" = "bin.45.50",
                                             "B19001B_011" = "bin.50.60",
                                             "B19001B_012" = "bin.60.75",
                                             "B19001B_013" = "bin.75.100",
                                             "B19001B_014" = "bin.100.125",
                                             "B19001B_015" = "bin.125.150",
                                             "B19001B_016" = "bin.150.200",
                                             "B19001B_017" = "bin.200.above")) %>%
                    mutate(se = moe / 1.645)

                median_ests_white = get_acs(geography = "tract",
                                            table = "B19013A",
                                            state = curr_state,
                                            county = countyFips,
                                            cache_table = TRUE,
                                            year = curr_year) %>%
                    mutate(variable = recode(variable,
                                             "B19013A_001" = "median")) %>%
                    mutate(estimate = estimate / 1000, moe = moe / 1000,
                           se = moe / 1.645)

                median_ests_black = get_acs(geography = "tract",
                                            table = "B19013B",
                                            state = curr_state,
                                            county = countyFips,
                                            cache_table = TRUE,
                                            year = curr_year) %>%
                    mutate(variable = recode(variable,
                                             "B19013B_001" = "median")) %>%
                    mutate(estimate = estimate / 1000, moe = moe / 1000,
                           se = moe / 1.645)

                agg_ests_white = get_acs(geography = "tract",
                                         table = "B19025A",
                                         state = curr_state,
                                         county = countyFips,
                                         cache_table = TRUE,
                                         year = curr_year) %>%
                    mutate(variable = recode(variable,
                                             "B19025A_001" = "aggregate")) %>%
                    mutate(estimate = estimate / 1000, moe = moe / 1000,
                           se = moe / 1.645)

                agg_ests_black = get_acs(geography = "tract",
                                         table = "B19025B",
                                         state = curr_state,
                                         county = countyFips,
                                         cache_table = TRUE,
                                         year = curr_year) %>%
                    mutate(variable = recode(variable,
                                             "B19025B_001" = "aggregate")) %>%
                    mutate(estimate = estimate / 1000, moe = moe / 1000,
                           se = moe / 1.645)

                race_pop_ests = get_acs(geography = "tract",
                                        table = "B02001",
                                        state = curr_state,
                                        county = countyFips,
                                        cache_table = TRUE,
                                        year = curr_year) %>%
                    mutate(variable = recode(variable,
                                             "B02001_001" = "pop",
                                             "B02001_002" = "pop.white",
                                             "B02001_003" = "pop.black",
                                             "B02001_004" = "pop.native",
                                             "B02001_005" = "pop.asian",
                                             "B02001_006" = "pop.islander",
                                             "B02001_007" = "pop.other",
                                             "B02001_008" = "pop.twoplus")) %>%
                    filter(!grepl("B02001", variable)) %>%
                    mutate(se = moe / 1.645)

                unemp_ests_new = get_acs(geography = "tract",
                                         table = "S2301",
                                         state = curr_state,
                                         county = countyFips,
                                         cache_table = TRUE,
                                         year = curr_year) %>%
                    mutate(variable = recode(variable,
                                             "S2301_C02_012" = "labor.part.rate.white",
                                             "S2301_C02_013" = "labor.part.rate.black",
                                             "S2301_C03_012" = "employ.pop.ratio.white",
                                             "S2301_C03_013" = "employ.pop.ratio.black",
                                             "S2301_C04_012" = "unemploy.rate.white",
                                             "S2301_C04_013" = "unemploy.rate.black"
                                             )) %>%
                    filter(!grepl("S2301", variable)) %>%
                    mutate(estimate = estimate / 100, moe = moe / 100,
                           se = moe / 1.645)

                edu_ests_new = get_acs(geography = "tract",
                                       table = "S1501",
                                       state = curr_state,
                                       county = countyFips,
                                       cache_table = TRUE,
                                       year = curr_year) %>%
                    mutate(variable = recode(variable,
                                             "S1501_C02_029" =
                                                 "edu.25up.hs.plus.white",
                                             "S1501_C02_035" =
                                                 "edu.25up.hs.plus.black")) %>%
                    filter(!grepl("S1501", variable)) %>%
                    mutate(estimate = estimate / 100, moe = moe / 100,
                           se = moe / 1.645)

                femalehead_ests_new = get_acs(geography = "tract",
                                              table = "B09019",
                                              state = curr_state,
                                              county = countyFips,
                                              cache_table = TRUE,
                                              year = curr_year) %>%
                    mutate(variable = recode(variable,
                                             "B09019_005" =
                                                 "total.family.male.householder",
                                             "B09019_006" =
                                                 "total.family.female.householder"
                                             )) %>%
                    filter(!grepl("B09019", variable)) %>%
                    mutate(se = moe / 1.645)


                bin_ests_white_wide = bin_ests_white %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable,
                                values_from = estimate) %>%
                    arrange(GEOID) %>%
                    mutate_at(vars(matches("bin")), list(~ ./total.white))

                bin_ses_white_wide = bin_ests_white %>% select(-moe, -estimate) %>%
                    pivot_wider(names_from = variable,
                                values_from = se) %>%
                    arrange(GEOID) %>%
                    mutate(bin.below.10 =
                               countse2propse(bin.below.10, total.white,
                                              bin_ests_white_wide$bin.below.10,
                                              bin_ests_white_wide$total.white),
                           bin.10.15 =
                               countse2propse(bin.10.15, total.white,
                                              bin_ests_white_wide$bin.10.15,
                                              bin_ests_white_wide$total.white),
                           bin.15.20 =
                               countse2propse(bin.15.20, total.white,
                                              bin_ests_white_wide$bin.15.20,
                                              bin_ests_white_wide$total.white),
                           bin.20.25 =
                               countse2propse(bin.20.25, total.white,
                                              bin_ests_white_wide$bin.20.25,
                                              bin_ests_white_wide$total.white),
                           bin.25.30 =
                               countse2propse(bin.25.30, total.white,
                                              bin_ests_white_wide$bin.25.30,
                                              bin_ests_white_wide$total.white),
                           bin.30.35 =
                               countse2propse(bin.30.35, total.white,
                                              bin_ests_white_wide$bin.30.35,
                                              bin_ests_white_wide$total.white),
                           bin.35.40 =
                               countse2propse(bin.35.40, total.white,
                                              bin_ests_white_wide$bin.35.40,
                                              bin_ests_white_wide$total.white),
                           bin.40.45 =
                               countse2propse(bin.40.45, total.white,
                                              bin_ests_white_wide$bin.40.45,
                                              bin_ests_white_wide$total.white),
                           bin.45.50 =
                               countse2propse(bin.45.50, total.white,
                                              bin_ests_white_wide$bin.45.50,
                                              bin_ests_white_wide$total.white),
                           bin.50.60 =
                               countse2propse(bin.50.60, total.white,
                                              bin_ests_white_wide$bin.50.60,
                                              bin_ests_white_wide$total.white),
                           bin.60.75 =
                               countse2propse(bin.60.75, total.white,
                                              bin_ests_white_wide$bin.60.75,
                                              bin_ests_white_wide$total.white),
                           bin.75.100 =
                               countse2propse(bin.75.100, total.white,
                                              bin_ests_white_wide$bin.75.100,
                                              bin_ests_white_wide$total.white),
                           bin.100.125 =
                               countse2propse(bin.100.125, total.white,
                                              bin_ests_white_wide$bin.100.125,
                                              bin_ests_white_wide$total.white),
                           bin.125.150 =
                               countse2propse(bin.125.150, total.white,
                                              bin_ests_white_wide$bin.125.150,
                                              bin_ests_white_wide$total.white),
                           bin.150.200 =
                               countse2propse(bin.150.200, total.white,
                                              bin_ests_white_wide$bin.150.200,
                                              bin_ests_white_wide$total.white),
                           bin.200.above =
                               countse2propse(bin.200.above, total.white,
                                              bin_ests_white_wide$bin.200.above,
                                              bin_ests_white_wide$total.white))

                bin_ests_black_wide = bin_ests_black %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable,
                                values_from = estimate) %>%
                    arrange(GEOID) %>%
                    mutate_at(vars(matches("bin")), list(~ ./total.black))

                bin_ses_black_wide = bin_ests_black %>% select(-moe, -estimate) %>%
                    pivot_wider(names_from = variable,
                                values_from = se) %>%
                    arrange(GEOID) %>%
                    mutate(bin.below.10 =
                               countse2propse(bin.below.10, total.black,
                                              bin_ests_black_wide$bin.below.10,
                                              bin_ests_black_wide$total.black),
                           bin.10.15 =
                               countse2propse(bin.10.15, total.black,
                                              bin_ests_black_wide$bin.10.15,
                                              bin_ests_black_wide$total.black),
                           bin.15.20 =
                               countse2propse(bin.15.20, total.black,
                                              bin_ests_black_wide$bin.15.20,
                                              bin_ests_black_wide$total.black),
                           bin.20.25 =
                               countse2propse(bin.20.25, total.black,
                                              bin_ests_black_wide$bin.20.25,
                                              bin_ests_black_wide$total.black),
                           bin.25.30 =
                               countse2propse(bin.25.30, total.black,
                                              bin_ests_black_wide$bin.25.30,
                                              bin_ests_black_wide$total.black),
                           bin.30.35 =
                               countse2propse(bin.30.35, total.black,
                                              bin_ests_black_wide$bin.30.35,
                                              bin_ests_black_wide$total.black),
                           bin.35.40 =
                               countse2propse(bin.35.40, total.black,
                                              bin_ests_black_wide$bin.35.40,
                                              bin_ests_black_wide$total.black),
                           bin.40.45 =
                               countse2propse(bin.40.45, total.black,
                                              bin_ests_black_wide$bin.40.45,
                                              bin_ests_black_wide$total.black),
                           bin.45.50 =
                               countse2propse(bin.45.50, total.black,
                                              bin_ests_black_wide$bin.45.50,
                                              bin_ests_black_wide$total.black),
                           bin.50.60 =
                               countse2propse(bin.50.60, total.black,
                                              bin_ests_black_wide$bin.50.60,
                                              bin_ests_black_wide$total.black),
                           bin.60.75 =
                               countse2propse(bin.60.75, total.black,
                                              bin_ests_black_wide$bin.60.75,
                                              bin_ests_black_wide$total.black),
                           bin.75.100 =
                               countse2propse(bin.75.100, total.black,
                                              bin_ests_black_wide$bin.75.100,
                                              bin_ests_black_wide$total.black),
                           bin.100.125 =
                               countse2propse(bin.100.125, total.black,
                                              bin_ests_black_wide$bin.100.125,
                                              bin_ests_black_wide$total.black),
                           bin.125.150 =
                               countse2propse(bin.125.150, total.black,
                                              bin_ests_black_wide$bin.125.150,
                                              bin_ests_black_wide$total.black),
                           bin.150.200 =
                               countse2propse(bin.150.200, total.black,
                                              bin_ests_black_wide$bin.150.200,
                                              bin_ests_black_wide$total.black),
                           bin.200.above =
                               countse2propse(bin.200.above, total.black,
                                              bin_ests_black_wide$bin.200.above,
                                              bin_ests_black_wide$total.black))

                median_ests_white_wide = median_ests_white %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable,
                                values_from = estimate)

                median_ests_black_wide = median_ests_black %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable,
                                values_from = estimate)

                median_ses_white_wide = median_ests_white %>% select(-moe, -estimate) %>%
                    pivot_wider(names_from = variable,
                                values_from = se)

                median_ses_black_wide = median_ests_black %>% select(-moe, -estimate) %>%
                    pivot_wider(names_from = variable,
                                values_from = se)

                mean_ests_white_wide = agg_ests_white %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable,
                                values_from = estimate) %>%
                    full_join(bin_ests_white_wide) %>%
                    arrange(GEOID) %>%
                    transmute(GEOID = GEOID,
                              NAME = NAME,
                              total = total.white,
                              mean = aggregate / total.white)

                mean_ests_black_wide = agg_ests_black %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable,
                                values_from = estimate) %>%
                    full_join(bin_ests_black_wide) %>%
                    arrange(GEOID) %>%
                    transmute(GEOID = GEOID,
                              NAME = NAME,
                              total = total.black,
                              mean = aggregate / total.black)

                mean_ses_white_wide = agg_ests_white %>% select(-moe, -estimate) %>%
                    pivot_wider(names_from = variable,
                                values_from = se) %>%
                    full_join(bin_ses_black_wide) %>%
                    arrange(GEOID) %>%
                    transmute(GEOID = GEOID,
                              NAME = NAME,
                              mean = aggse2meanse(aggregate,
                                                  total.black,
                                                  mean_ests_white_wide$mean,
                                                  mean_ests_white_wide$total))

                mean_ses_black_wide = agg_ests_black %>% select(-moe, -estimate) %>%
                    pivot_wider(names_from = variable,
                                values_from = se) %>%
                    full_join(bin_ses_black_wide) %>%
                    arrange(GEOID) %>%
                    transmute(GEOID = GEOID,
                              NAME = NAME,
                              mean = aggse2meanse(aggregate,
                                                  total.black,
                                                  mean_ests_black_wide$mean,
                                                  mean_ests_black_wide$total))

                mean_ests_white_wide = mean_ests_white_wide %>% select(-total)
                mean_ests_black_wide = mean_ests_black_wide %>% select(-total)

                race_pop_ests_wide = race_pop_ests %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable, values_from = estimate)

                race_pop_ses_wide = race_pop_ests %>% select(-moe, -estimate) %>%
                    pivot_wider(names_from = variable, values_from = se)

                unemp_ests_wide = unemp_ests_new %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable, values_from = estimate)
                unemp_ses_wide = unemp_ests_new %>% select(-estimate, -moe) %>%
                    pivot_wider(names_from = variable, values_from = se)

                edu_ests_wide = edu_ests_new %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable, values_from = estimate)
                edu_ses_wide = edu_ests_new %>% select(-estimate, -moe) %>%
                    pivot_wider(names_from = variable, values_from = se)

                percapita_ests_white_wide = agg_ests_white %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable,
                                values_from = estimate) %>%
                    full_join(race_pop_ests_wide) %>%
                    arrange(GEOID) %>%
                    transmute(GEOID = GEOID,
                              NAME = NAME,
                              pop = pop.white,
                              percapita.income = aggregate / pop.white)
                percapita_ses_white_wide = agg_ests_white %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable,
                                values_from = estimate) %>%
                    full_join(race_pop_ses_wide) %>%
                    arrange(GEOID) %>%
                    transmute(GEOID = GEOID,
                              NAME = NAME,
                              pop = pop.white,
                              percapita.income =
                                  aggse2meanse(aggregate, pop,
                                               percapita_ests_white_wide$percapita.income,
                                               percapita_ests_white_wide$pop))
                percapita_ests_white_wide = percapita_ests_white_wide %>%
                    select(-pop)
                percapita_ests_white_wide = percapita_ses_white_wide %>%
                    select(-pop)

                percapita_ests_black_wide = agg_ests_black %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable,
                                values_from = estimate) %>%
                    full_join(race_pop_ests_wide) %>%
                    arrange(GEOID) %>%
                    transmute(GEOID = GEOID,
                              NAME = NAME,
                              pop = pop.black,
                              percapita.income = aggregate / pop.black)
                percapita_ses_black_wide = agg_ests_black %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable,
                                values_from = estimate) %>%
                    full_join(race_pop_ses_wide) %>%
                    arrange(GEOID) %>%
                    transmute(GEOID = GEOID,
                              NAME = NAME,
                              pop = pop.black,
                              percapita.income =
                                  aggse2meanse(aggregate, pop,
                                               percapita_ests_black_wide$percapita.income,
                                               percapita_ests_black_wide$pop))
                percapita_ests_black_wide = percapita_ests_black_wide %>%
                    select(-pop)
                percapita_ests_black_wide = percapita_ses_black_wide %>%
                    select(-pop)

                full_ests_white = race_pop_ests_wide %>%
                    select(GEOID, NAME, pop.white) %>%
                    full_join(bin_ests_white_wide) %>%
                    full_join(mean_ests_white_wide) %>%
                    full_join(median_ests_white_wide) %>%
                    full_join(unemp_ests_wide) %>%
                    full_join(edu_ests_wide) %>%
                    full_join(percapita_ests_white_wide) %>%
                    select(-contains("black")) %>%
                    select(GEOID, NAME, total.white, mean, median,
                           contains("bin"), everything())

                full_ests_black = race_pop_ests_wide %>%
                    select(GEOID, NAME, pop.black) %>%
                    full_join(bin_ests_black_wide) %>%
                    full_join(mean_ests_black_wide) %>%
                    full_join(median_ests_black_wide) %>%
                    full_join(unemp_ests_wide) %>%
                    full_join(edu_ests_wide) %>%
                    full_join(percapita_ests_black_wide) %>%
                    select(-contains("white")) %>%
                    select(GEOID, NAME, total.black, mean, median,
                           contains("bin"), everything())

                full_ses_white = race_pop_ses_wide %>%
                    select(GEOID, NAME, pop.white) %>%
                    full_join(bin_ses_white_wide) %>%
                    full_join(mean_ses_white_wide) %>%
                    full_join(median_ses_white_wide) %>%
                    full_join(unemp_ses_wide) %>%
                    full_join(edu_ses_wide) %>%
                    full_join(percapita_ses_white_wide) %>%
                    select(-contains("black")) %>%
                    select(GEOID, NAME, total.white, mean, median,
                           contains("bin"), everything())

                full_ses_black = race_pop_ses_wide %>%
                    select(GEOID, NAME, pop.black) %>%
                    full_join(bin_ses_black_wide) %>%
                    full_join(mean_ses_black_wide) %>%
                    full_join(median_ses_black_wide) %>%
                    full_join(unemp_ses_wide) %>%
                    full_join(edu_ses_wide) %>%
                    full_join(percapita_ses_black_wide) %>%
                    select(-contains("white")) %>%
                    select(GEOID, NAME, total.black, mean, median,
                           contains("bin"), everything())

                full_ests = full_join(full_ests_white, full_ests_black,
                                      by = c('GEOID', 'NAME'),
                                      suffix = c(".white", ".black")) %>%
                    select(GEOID, NAME, contains("white"), contains("black"))

                full_ses = full_join(full_ses_white, full_ses_black,
                                     by = c('GEOID', 'NAME'),
                                     suffix = c(".white", ".black")) %>%
                    select(GEOID, NAME, contains("white"), contains("black"))

                out = full_join(full_ests, full_ses,
                                by = c('GEOID', 'NAME'), suffix = c(".est", ".se"))

                geom = get_acs(geography = "tract",
                               variables = "B19013_001",
                               state = curr_state,
                               county = countyFips,
                               cache_table = TRUE,
                               year = curr_year,
                               geometry = TRUE) %>%
                    select(-estimate, -moe, -variable)

                out_geom = full_join(geom, out) %>%
                    mutate(metro = metroName, metroFips = curr_metro,
                           state = curr_state, year = curr_year) %>%
                    select(GEOID, NAME, metro, metroFips, state, year, everything())

                out_all = rbind(out_all, out_geom)
            }
        }
    }
    out_all
}

create_metro_tract_dataset = function(years){

    pums2tract = read_csv("data/tract2puma2010.txt")
    out_all = NULL

    metros = get_acs(
        geography = "metropolitan statistical area/micropolitan statistical area",
        table = "B01003",
        cache_table = TRUE,
        year = 2018) %>%    ## using top 100 metros of 2018
        top_n(100, estimate) %>%
        arrange(desc(estimate))

    metro_deliniation = read_csv("data/metro_deliniation_4-2018.csv") %>%
        filter(as.character(`CBSA Code`) %in% metros$GEOID) %>%
        mutate(GEOID = as.character(`CBSA Code`),
               MetroName = `CBSA Title`,
               StateName = `State Name`,
               County = `County/County Equivalent`,
               StateFips = `FIPS State Code`,
               CountyFips = `FIPS County Code`) %>%
        select(GEOID, MetroName, StateName, County, StateFips, CountyFips)

    for(curr_metro in metros$GEOID){
        metroName = metro_deliniation %>%
            filter(GEOID == curr_metro) %>%
            .[['MetroName']] %>%
            unique

        stateFips = metro_deliniation %>%
            filter(GEOID == curr_metro) %>%
            .[['StateFips']] %>%
            unique

        states = fips_codes %>%
            filter(state_code %in% stateFips) %>%
            .[['state']] %>%
            unique

        for(curr_state in states){

            fips = fips_codes %>%
                filter(state == curr_state) %>%
                select(state_code) %>% .[1,1]

            countyFips = metro_deliniation %>%
                filter(GEOID == curr_metro, StateFips == fips) %>%
                .[['CountyFips']] %>%
                unique

            tractFips = pums2tract %>%
                filter(STATEFP == fips, COUNTYFP %in% countyFips) %>%
                mutate(TractFips = paste(fips, COUNTYFP, TRACTCE, sep = "")) %>%
                .[['TractFips']] %>%
                unique

            for(curr_year in years){
                mean_ests_new = get_acs(geography = "tract",
                                        table = "S1901",
                                        state = curr_state,
                                        county = countyFips,
                                        cache_table = TRUE,
                                        year = curr_year) %>%
                    filter(grepl("C01", variable)) %>%
                    mutate(variable = recode(variable,
                                             "S1901_C01_001" = "total",
                                             "S1901_C01_002" = "bin.below.10",
                                             "S1901_C01_003" = "bin.10.15",
                                             "S1901_C01_004" = "bin.15.25",
                                             "S1901_C01_005" = "bin.25.35",
                                             "S1901_C01_006" = "bin.35.50",
                                             "S1901_C01_007" = "bin.50.75",
                                             "S1901_C01_008" = "bin.75.100",
                                             "S1901_C01_009" = "bin.100.150",
                                             "S1901_C01_010" = "bin.150.200",
                                             "S1901_C01_011" = "bin.200.above",
                                             "S1901_C01_012" = "median",
                                             "S1901_C01_013" = "mean")) %>%
                    filter(!grepl("S1901", variable)) %>%
                    mutate(estimate = ifelse(grepl("bin", variable),
                                             estimate / 100, estimate),
                           moe = ifelse(grepl("bin", variable),
                                        moe / 100, moe)) %>%
                    mutate(estimate = ifelse(grepl("mean", variable),
                                             estimate / 1000, estimate),
                           moe = ifelse(grepl("mean", variable),
                                        moe / 1000, moe)) %>%
                    mutate(estimate = ifelse(grepl("median", variable),
                                             estimate / 1000, estimate),
                           moe = ifelse(grepl("median", variable),
                                        moe / 1000, moe)) %>%
                    mutate(se = moe / 1.645)

                bin_ests_new = get_acs(geography = "tract",
                                       table = "S2503",
                                       state = curr_state,
                                       county = countyFips,
                                       cache_table = TRUE,
                                       year = curr_year) %>%
                    filter(grepl("C01", variable)) %>%
                    mutate(variable = recode(variable,
                                             "S2503_C01_001" = "total",
                                             "S2503_C01_002" = "bin.below.5",
                                             "S2503_C01_003" = "bin.5.10",
                                             "S2503_C01_004" = "bin.10.15",
                                             "S2503_C01_005" = "bin.15.20",
                                             "S2503_C01_006" = "bin.20.25",
                                             "S2503_C01_007" = "bin.25.35",
                                             "S2503_C01_008" = "bin.35.50",
                                             "S2503_C01_009" = "bin.50.75",
                                             "S2503_C01_010" = "bin.75.100",
                                             "S2503_C01_011" = "bin.100.150",
                                             "S2503_C01_012" = "bin.150.above",
                                             "S2503_C01_013" = "median")) %>%
                    filter(!grepl("S2503", variable)) %>%
                    mutate(estimate = ifelse(grepl("bin", variable),
                                             estimate / 100, estimate),
                           moe = ifelse(grepl("bin", variable),
                                        moe / 100, moe)) %>%
                    mutate(estimate = ifelse(grepl("median", variable),
                                             estimate / 1000, estimate),
                           moe = ifelse(grepl("median", variable),
                                        moe / 1000, moe)) %>%
                    mutate(se = moe / 1.645)

                gini_ests_new = get_acs(geography = "tract",
                                        table = "B19083",
                                        state = curr_state,
                                        county = countyFips,
                                        cache_table = TRUE,
                                        year = curr_year) %>%
                    mutate(variable = recode(variable, "B19083_001" = "gini")) %>%
                    mutate(se = moe / 1.645)

                quant_ests_new = get_acs(geography = "tract",
                                         table = "B19080",
                                         state = curr_state,
                                         county = countyFips,
                                         cache_table = TRUE,
                                         year = curr_year) %>%
                    mutate(variable = recode(variable,
                                             "B19080_001" = "quant.20",
                                             "B19080_002" = "quant.40",
                                             "B19080_003" = "quant.60",
                                             "B19080_004" = "quant.80",
                                             "B19080_005" = "quant.95")) %>%
                    mutate(estimate = estimate / 1000, moe = moe / 1000,
                           se = moe / 1.645)

                share_ests_new = get_acs(geography = "tract",
                                         table = "B19082",
                                         state = curr_state,
                                         county = countyFips,
                                         cache_table = TRUE,
                                         year = curr_year) %>%
                    mutate(variable = recode(variable,
                                             "B19082_001" = "share.0.20",
                                             "B19082_002" = "share.20.40",
                                             "B19082_003" = "share.40.60",
                                             "B19082_004" = "share.60.80",
                                             "B19082_005" = "share.80.100",
                                             "B19082_006" = "share.95.100")) %>%
                    mutate(estimate = estimate / 100, moe = moe / 100,
                           se = moe / 1.645)

                unemp_ests_new = get_acs(geography = "tract",
                                         table = "S2301",
                                         state = curr_state,
                                         county = countyFips,
                                         cache_table = TRUE,
                                         year = curr_year) %>%
                    mutate(variable = recode(variable,
                                             "S2301_C02_001" = "labor.part.rate",
                                             "S2301_C03_001" = "employ.pop.ratio",
                                             "S2301_C04_001" = "unemploy.rate")) %>%
                    filter(!grepl("S2301", variable)) %>%
                    mutate(estimate = estimate / 100, moe = moe / 100,
                           se = moe / 1.645)

                age_ests_new = get_acs(geography = "tract",
                                       table = "S0101",
                                       state = curr_state,
                                       county = countyFips,
                                       cache_table = TRUE,
                                       year = curr_year) %>%
                    mutate(variable = recode(variable,
                                             "S0101_C02_026" = "age.18.over",
                                             "S0101_C02_030" = "age.65.over")) %>%
                    filter(!grepl("S0101", variable)) %>%
                    mutate(estimate = estimate / 100, moe = moe / 100,
                           se = moe / 1.645)

                edu_ests_new = get_acs(geography = "tract",
                                       table = "S1501",
                                       state = curr_state,
                                       county = countyFips,
                                       cache_table = TRUE,
                                       year = curr_year) %>%
                    mutate(variable = recode(variable,
                                             "S1501_C02_014" = "edu.25up.hs.plus")) %>%
                    filter(!grepl("S1501", variable)) %>%
                    mutate(estimate = estimate / 100, moe = moe / 100,
                           se = moe / 1.645)

                percapita_ests_new = get_acs(geography = "tract",
                                             table = "B19301",
                                             state = curr_state,
                                             county = countyFips,
                                             cache_table = TRUE,
                                             year = curr_year) %>%
                    mutate(variable = recode(variable,
                                             "B19301_001" = "per.capita.income")) %>%
                    mutate(estimate = estimate / 1000, moe = moe / 1000,
                           se = moe / 1.645)

                foreign_ests_new = get_acs(geography = "tract",
                                           table = "DP02",
                                           state = curr_state,
                                           county = countyFips,
                                           cache_table = TRUE,
                                           year = curr_year) %>%
                    mutate(variable = recode(variable,
                                             "DP02_0092P" = "foreign.born")) %>%
                    filter(!grepl("DP02", variable)) %>%
                    mutate(estimate = estimate / 100, moe = moe / 100,
                           se = moe / 1.645)

                occu_ests_new = get_acs(geography = "tract",
                                        table = "S2405",
                                        state = curr_state,
                                        county = countyFips,
                                        cache_table = TRUE,
                                        year = curr_year) %>%
                    mutate(variable = recode(variable,
                                             "S2405_C01_001" = "total.civilian",
                                             "S2405_C01_003" = "total.construction",
                                             "S2405_C01_004" = "total.manufacturing",
                                             "S2405_C01_009" = "total.FIRE",
                                             "S2405_C01_008" = "total.information",
                                             "S2405_C01_010" = "total.professional",
                                             "S2405_C01_011" = "total.edu.health",
                                             "S2405_C01_014" = "total.public.admin")) %>%
                    filter(!grepl("S2405", variable)) %>%
                    mutate(se = moe / 1.645)

                femalehead_ests_new = get_acs(geography = "tract",
                                              table = "B09019",
                                              state = curr_state,
                                              county = countyFips,
                                              cache_table = TRUE,
                                              year = curr_year) %>%
                    mutate(variable = recode(variable,
                                             "B09019_005" =
                                                 "total.family.male.householder",
                                             "B09019_006" =
                                                 "total.family.female.householder"
                                             )) %>%
                    filter(!grepl("B09019", variable)) %>%
                    mutate(se = moe / 1.645)

                mobility_ests_new = get_acs(geography = "tract",
                                            table = "B07204",
                                            state = curr_state,
                                            county = countyFips,
                                            cache_table = TRUE,
                                            year = curr_year) %>%
                    mutate(variable = recode(variable,
                                             "B07204_001" = "pop",
                                             "B07204_002" = "pop.same.house",
                                             "B07204_005" =
                                                 "pop.same.town.same.county",
                                             "B07204_008" =
                                                 "pop.diff.town.same.county"
                                             )) %>%
                    filter(!grepl("B07204", variable)) %>%
                    mutate(se = moe / 1.645)

                housing_ests_new = get_acs(geography = "tract",
                                           table = "DP04",
                                           state = curr_state,
                                           county = countyFips,
                                           cache_table = TRUE,
                                           year = curr_year) %>%
                    mutate(variable = recode(variable,
                                             "DP04_0017P" = "house.built.less.5",
                                             )) %>%
                    filter(!grepl("DP04", variable)) %>%
                    mutate(estimate = estimate / 100, moe = moe / 100,
                           se = moe / 1.645)

                bin_ests_wide = bin_ests_new %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable,
                                values_from = estimate) %>%
                    select(-total, -median, -bin.150.above)
                bin_ses_wide = bin_ests_new %>% select(-moe, -estimate) %>%
                    pivot_wider(names_from = variable,
                                values_from = se) %>%
                    select(-total, -median, -bin.150.above)

                mean_ests_wide = mean_ests_new %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable,
                                values_from = estimate) %>%
                    select(GEOID, NAME, total, mean, median,
                           bin.150.200, bin.200.above)
                mean_ses_wide = mean_ests_new %>% select(-moe, -estimate) %>%
                    pivot_wider(names_from = variable,
                                values_from = se) %>%
                    select(GEOID, NAME, total, mean, median,
                           bin.150.200, bin.200.above)

                quant_ests_wide = quant_ests_new %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable, values_from = estimate)
                quant_ses_wide = quant_ests_new %>% select(-estimate, -moe) %>%
                    pivot_wider(names_from = variable, values_from = se)

                share_ests_wide = share_ests_new %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable, values_from = estimate)
                share_ses_wide = share_ests_new %>% select(-estimate, -moe) %>%
                    pivot_wider(names_from = variable, values_from = se)

                gini_ests_wide = gini_ests_new %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable, values_from = estimate)
                gini_ses_wide = gini_ests_new %>% select(-estimate, -moe) %>%
                    pivot_wider(names_from = variable, values_from = se)

                unemp_ests_wide = unemp_ests_new %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable, values_from = estimate)
                unemp_ses_wide = unemp_ests_new %>% select(-estimate, -moe) %>%
                    pivot_wider(names_from = variable, values_from = se)

                age_ests_wide = age_ests_new %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable, values_from = estimate) %>%
                    mutate(age.under.18 = 1 - age.18.over) %>%
                    select(-age.18.over)
                age_ses_wide = age_ests_new %>% select(-estimate, -moe) %>%
                    pivot_wider(names_from = variable, values_from = se) %>%
                    mutate(age.under.18 = age.18.over) %>%
                    select(-age.18.over)

                edu_ests_wide = edu_ests_new %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable, values_from = estimate)
                edu_ses_wide = edu_ests_new %>% select(-estimate, -moe) %>%
                    pivot_wider(names_from = variable, values_from = se)

                percapita_ests_wide = percapita_ests_new %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable, values_from = estimate)
                percapita_ses_wide = percapita_ests_new %>% select(-estimate, -moe) %>%
                    pivot_wider(names_from = variable, values_from = se)

                foreign_ests_wide = foreign_ests_new %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable, values_from = estimate)
                foreign_ses_wide = foreign_ests_new %>% select(-estimate, -moe) %>%
                    pivot_wider(names_from = variable, values_from = se)

                occu_ests_wide = occu_ests_new %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable, values_from = estimate) %>%
                    arrange(GEOID) %>%
                    transmute(
                        GEOID = GEOID,
                        NAME = NAME,
                        total = total.civilian,
                        construction = total.construction / total.civilian,
                        manufacturing = total.manufacturing / total.civilian,
                        fire = total.FIRE / total.civilian,
                        professional =
                            (total.FIRE +
                             total.information +
                             total.professional +
                             total.edu.health +
                             total.public.admin) / total.civilian)
                occu_ses_wide = occu_ests_new %>% select(-estimate, -moe) %>%
                    pivot_wider(names_from = variable, values_from = se) %>%
                    arrange(GEOID) %>%
                    transmute(
                        GEOID = GEOID,
                        NAME = NAME,
                        construction =
                            countse2propse(total.construction,
                                           total.civilian,
                                           occu_ests_wide$construction,
                                           occu_ests_wide$total),

                        manufacturing =
                            countse2propse(total.manufacturing,
                                           total.civilian,
                                           occu_ests_wide$manufacturing,
                                           occu_ests_wide$total),
                        fire =
                            countse2propse(total.FIRE,
                                           total.civilian,
                                           occu_ests_wide$fire,
                                           occu_ests_wide$total),
                        professional =
                            countse2propse(
                                sqrt(total.FIRE^2 + total.information^2 +
                                     total.professional^2 + total.edu.health^2),
                                total.civilian,
                                occu_ests_wide$professional,
                                occu_ests_wide$total))
                occu_ests_wide = occu_ests_wide %>% select(-total)

                femalehead_ests_wide = femalehead_ests_new %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable, values_from = estimate) %>%
                    arrange(GEOID) %>%
                    transmute(GEOID = GEOID,
                              NAME = NAME,
                              total = total.family.male.householder +
                                  total.family.female.householder,
                              female.hher = total.family.female.householder / total)

                femalehead_ses_wide = femalehead_ests_new %>% select(-estimate, -moe) %>%
                    pivot_wider(names_from = variable, values_from = se) %>%
                    arrange(GEOID) %>%
                    transmute(GEOID = GEOID,
                              NAME = NAME,
                              female.hher =
                                  countse2propse(
                                      total.family.female.householder,
                                      sqrt(total.family.female.householder^2 +
                                           total.family.male.householder^2),
                                      femalehead_ests_wide$female.hher,
                                      femalehead_ests_wide$total))
                femalehead_ests_wide = femalehead_ests_wide %>% select(-total)

                mobility_ests_wide = mobility_ests_new %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable, values_from = estimate) %>%
                    arrange(GEOID) %>%
                    transmute(GEOID = GEOID,
                              NAME = NAME,
                              pop = pop,
                              same.house = pop.same.house / pop,
                              same.county =
                                  (pop.same.town.same.county +
                                   pop.diff.town.same.county) / pop)
                mobility_ses_wide = mobility_ests_new %>% select(-estimate, -moe) %>%
                    pivot_wider(names_from = variable, values_from = se) %>%
                    arrange(GEOID) %>%
                    transmute(GEOID = GEOID,
                              NAME = NAME,
                              pop = pop,
                              same.house =
                                  countse2propse(
                                      pop.same.house,
                                      pop,
                                      mobility_ests_wide$same.house,
                                      mobility_ests_wide$pop),
                              same.county =
                                  countse2propse(
                                      sqrt(pop.same.town.same.county^2 +
                                           pop.diff.town.same.county^2),
                                      pop,
                                      mobility_ests_wide$same.county,
                                      mobility_ests_wide$pop))

                housing_ests_wide = housing_ests_new %>% select(-moe, -se) %>%
                    pivot_wider(names_from = variable, values_from = estimate)
                housing_ses_wide = housing_ests_new %>% select(-estimate, -moe) %>%
                    pivot_wider(names_from = variable, values_from = se)

                bounds1 = c(5, 10, 15, 20, 25, 35, 50, 75, 100, 150)
                bounds2 = c(10, 15, 25, 35, 50, 75, 100, 150, 200)
                bounds = unique(c(bounds1, bounds2))
                lbs = c("below", bounds)
                ubs = c(bounds, "above")

                full_ests = full_join(mean_ests_wide, bin_ests_wide) %>%
                    full_join(quant_ests_wide) %>%
                    full_join(share_ests_wide) %>%
                    full_join(gini_ests_wide) %>%
                    full_join(unemp_ests_wide) %>%
                    full_join(age_ests_wide) %>%
                    full_join(edu_ests_wide) %>%
                    full_join(percapita_ests_wide) %>%
                    full_join(foreign_ests_wide) %>%
                    full_join(occu_ests_wide) %>%
                    full_join(femalehead_ests_wide) %>%
                    full_join(mobility_ests_wide) %>%
                    full_join(housing_ests_wide) %>%
                    select(GEOID, NAME, total, mean, median,
                           paste("bin", lbs, ubs, sep="."),
                           contains("quant"), contains("share"), gini,
                           everything())

                full_ses = full_join(mean_ses_wide, bin_ses_wide) %>%
                    full_join(quant_ses_wide) %>%
                    full_join(share_ses_wide) %>%
                    full_join(gini_ses_wide) %>%
                    full_join(unemp_ses_wide) %>%
                    full_join(age_ses_wide) %>%
                    full_join(edu_ses_wide) %>%
                    full_join(percapita_ses_wide) %>%
                    full_join(foreign_ses_wide) %>%
                    full_join(occu_ses_wide) %>%
                    full_join(femalehead_ses_wide) %>%
                    full_join(mobility_ses_wide) %>%
                    full_join(housing_ses_wide) %>%
                    select(GEOID, NAME, total, mean, median,
                           paste("bin", lbs, ubs, sep="."),
                           contains("quant"), contains("share"), gini,
                           everything())

                out = full_join(full_ests, full_ses,
                                by = c('GEOID', 'NAME'), suffix = c(".est", ".se"))

                geom = get_acs(geography = "tract",
                               variables = "B19013_001",
                               state = curr_state,
                               county = countyFips,
                               cache_table = TRUE,
                               year = curr_year,
                               geometry = TRUE) %>%
                    select(-estimate, -moe, -variable)

                out_geom = full_join(geom, out) %>%
                    mutate(metro = metroName, metroFips = curr_metro,
                           state = curr_state, year = curr_year) %>%
                    select(GEOID, NAME, metro, metroFips, state, year, everything())

                out_all = rbind(out_all, out_geom)
            }
        }
    }
    out_all
}

create_state_tract_dataset = function(states, years){

    pums2tract = read_csv("data/tract2puma2010.txt")
    out_all = NULL

    for(curr_state in states){
        for(curr_year in years){

            mean_ests_new = get_acs(geography = "tract",
                                    table = "S1901",
                                    state = curr_state,
                                    cache_table = TRUE,
                                    year = curr_year) %>%
                filter(grepl("C01", variable)) %>%
                mutate(variable = recode(variable,
                                         "S1901_C01_001" = "total",
                                         "S1901_C01_002" = "bin.below.10",
                                         "S1901_C01_003" = "bin.10.15",
                                         "S1901_C01_004" = "bin.15.25",
                                         "S1901_C01_005" = "bin.25.35",
                                         "S1901_C01_006" = "bin.35.50",
                                         "S1901_C01_007" = "bin.50.75",
                                         "S1901_C01_008" = "bin.75.100",
                                         "S1901_C01_009" = "bin.100.150",
                                         "S1901_C01_010" = "bin.150.200",
                                         "S1901_C01_011" = "bin.200.above",
                                         "S1901_C01_012" = "median",
                                         "S1901_C01_013" = "mean")) %>%
                filter(!grepl("S1901", variable)) %>%
                mutate(estimate = ifelse(grepl("bin", variable),
                                         estimate / 100, estimate),
                       moe = ifelse(grepl("bin", variable), moe / 100, moe)) %>%
                mutate(estimate = ifelse(grepl("mean", variable),
                                         estimate / 1000, estimate),
                       moe = ifelse(grepl("mean", variable), moe / 1000, moe)) %>%
                mutate(estimate = ifelse(grepl("median", variable),
                                         estimate / 1000, estimate),
                       moe = ifelse(grepl("median", variable), moe / 1000, moe)) %>%
                mutate(se = moe / 1.645)

            bin_ests_new = get_acs(geography = "tract",
                                   table = "S2503",
                                   state = curr_state,
                                   cache_table = TRUE,
                                   year = curr_year) %>%
                filter(grepl("C01", variable)) %>%
                mutate(variable = recode(variable,
                                         "S2503_C01_001" = "total",
                                         "S2503_C01_002" = "bin.below.5",
                                         "S2503_C01_003" = "bin.5.10",
                                         "S2503_C01_004" = "bin.10.15",
                                         "S2503_C01_005" = "bin.15.20",
                                         "S2503_C01_006" = "bin.20.25",
                                         "S2503_C01_007" = "bin.25.35",
                                         "S2503_C01_008" = "bin.35.50",
                                         "S2503_C01_009" = "bin.50.75",
                                         "S2503_C01_010" = "bin.75.100",
                                         "S2503_C01_011" = "bin.100.150",
                                         "S2503_C01_012" = "bin.150.above",
                                         "S2503_C01_013" = "median")) %>%
                filter(!grepl("S2503", variable)) %>%
                mutate(estimate = ifelse(grepl("bin", variable),
                                         estimate / 100, estimate),
                       moe = ifelse(grepl("bin", variable), moe / 100, moe)) %>%
                mutate(estimate = ifelse(grepl("median", variable),
                                         estimate / 1000, estimate),
                       moe = ifelse(grepl("median", variable), moe / 1000, moe)) %>%
                mutate(se = moe / 1.645)

            gini_ests_new = get_acs(geography = "tract",
                                    table = "B19083",
                                    state = curr_state,
                                    cache_table = TRUE,
                                    year = curr_year) %>%
                mutate(variable = recode(variable, "B19083_001" = "gini")) %>%
                mutate(se = moe / 1.645)

            quant_ests_new = get_acs(geography = "tract",
                                     table = "B19080",
                                     state = curr_state,
                                     cache_table = TRUE,
                                     year = curr_year) %>%
                mutate(variable = recode(variable,
                                         "B19080_001" = "quant.20",
                                         "B19080_002" = "quant.40",
                                         "B19080_003" = "quant.60",
                                         "B19080_004" = "quant.80",
                                         "B19080_005" = "quant.95")) %>%
                mutate(estimate = estimate / 1000, moe = moe / 1000, se = moe / 1.645)

            share_ests_new = get_acs(geography = "tract",
                                     table = "B19082",
                                     state = curr_state,
                                     cache_table = TRUE,
                                     year = curr_year) %>%
                mutate(variable = recode(variable,
                                         "B19082_001" = "share.0.20",
                                         "B19082_002" = "share.20.40",
                                         "B19082_003" = "share.40.60",
                                         "B19082_004" = "share.60.80",
                                         "B19082_005" = "share.80.100",
                                         "B19082_006" = "share.95.100")) %>%
                mutate(estimate = estimate / 100, moe = moe / 100, se = moe / 1.645)

            bin_ests_wide = bin_ests_new %>% select(-moe, -se) %>%
                pivot_wider(names_from = variable, values_from = estimate) %>%
                select(-total, -median, -bin.150.above)
            bin_ses_wide = bin_ests_new %>% select(-moe, -estimate) %>%
                pivot_wider(names_from = variable, values_from = se) %>%
                select(-total, -median, -bin.150.above)

            mean_ests_wide = mean_ests_new %>% select(-moe, -se) %>%
                pivot_wider(names_from = variable, values_from = estimate) %>%
                select(GEOID, NAME, total, mean, median, bin.150.200, bin.200.above)
            mean_ses_wide = mean_ests_new %>% select(-moe, -estimate) %>%
                pivot_wider(names_from = variable, values_from = se) %>%
                select(GEOID, NAME, total, mean, median, bin.150.200, bin.200.above)

            quant_ests_wide = quant_ests_new %>% select(-moe, -se) %>%
                pivot_wider(names_from = variable, values_from = estimate)
            quant_ses_wide = quant_ests_new %>% select(-estimate, -moe) %>%
                pivot_wider(names_from = variable, values_from = se)

            share_ests_wide = share_ests_new %>% select(-moe, -se) %>%
                pivot_wider(names_from = variable, values_from = estimate)
            share_ses_wide = share_ests_new %>% select(-estimate, -moe) %>%
                pivot_wider(names_from = variable, values_from = se)

            gini_ests_wide = gini_ests_new %>% select(-moe, -se) %>%
                pivot_wider(names_from = variable, values_from = estimate)
            gini_ses_wide = gini_ests_new %>% select(-estimate, -moe) %>%
                pivot_wider(names_from = variable, values_from = se)

            bounds1 = c(5, 10, 15, 20, 25, 35, 50, 75, 100, 150)
            bounds2 = c(10, 15, 25, 35, 50, 75, 100, 150, 200)
            bounds = unique(c(bounds1, bounds2))
            lbs = c("below", bounds)
            ubs = c(bounds, "above")

            full_ests = full_join(mean_ests_wide, bin_ests_wide) %>%
                full_join(quant_ests_wide) %>%
                full_join(share_ests_wide) %>%
                full_join(gini_ests_wide) %>%
                select(GEOID, NAME, total, mean, median, paste("bin", lbs, ubs, sep="."),
                       contains("quant"), contains("share"), gini)

            full_ses = full_join(mean_ses_wide, bin_ses_wide) %>%
                full_join(quant_ses_wide) %>%
                full_join(share_ses_wide) %>%
                full_join(gini_ses_wide) %>%
                select(GEOID, NAME, total, mean, median, paste("bin", lbs, ubs, sep="."),
                       contains("quant"), contains("share"), gini)

            out = full_join(full_ests, full_ses,
                            by = c('GEOID', 'NAME'), suffix = c(".est", ".se"))

            geom = get_acs(geography = "tract",
                           variables = "B19013_001",
                           state = curr_state,
                           cache_table = TRUE,
                           year = curr_year,
                           geometry = TRUE) %>%
                select(-estimate, -moe, -variable)

            fips = fips_codes %>%
                filter(state == curr_state) %>%
                select(state_code) %>% .[1,1]

            state_pums2tract = pums2tract %>%
                filter(STATEFP == fips) %>%
                mutate(TRACT = str_pad(TRACTCE, 6, "left", "0"),
                       COUNTY = str_pad(COUNTYFP, 3, "left", "0"),
                       GEOID = paste(fips, COUNTY, TRACT, sep = "")) %>%
                transmute(GEOID = GEOID, puma = PUMA5CE)

            out_geom = full_join(geom, out) %>%
                mutate(state = curr_state, year = curr_year) %>%
                full_join(state_pums2tract) %>%
                select(GEOID, NAME, state, puma, year, everything())

            out_all = rbind(out_all, out_geom)
        }
    }
    out_all
}

create_us_bin_proportions = function(year){
    mean_ests_new = get_acs(geography = "US",
                            table = "S1901",
                            cache_table = TRUE,
                            year = year) %>%
        filter(grepl("C01", variable)) %>%
        mutate(variable = recode(variable,
                                 "S1901_C01_001" = "total",
                                 "S1901_C01_002" = "bin.below.10",
                                 "S1901_C01_003" = "bin.10.15",
                                 "S1901_C01_004" = "bin.15.25",
                                 "S1901_C01_005" = "bin.25.35",
                                 "S1901_C01_006" = "bin.35.50",
                                 "S1901_C01_007" = "bin.50.75",
                                 "S1901_C01_008" = "bin.75.100",
                                 "S1901_C01_009" = "bin.100.150",
                                 "S1901_C01_010" = "bin.150.200",
                                 "S1901_C01_011" = "bin.200.above",
                                 "S1901_C01_012" = "median",
                                 "S1901_C01_013" = "mean")) %>%
        filter(grepl("bin", variable)) %>%
        mutate(estimate = estimate / 100,
               moe = moe / 100,
               se = moe / 1.645)

    bin_ests_new = get_acs(geography = "US",
                           table = "S2503",
                           cache_table = TRUE,
                           year = year)
    if(year <= 2016){
        bin_ests_new = bin_ests_new %>%
            mutate(variable = recode(variable,
                                     "S2503_C01_001" = "total",
                                     "S2503_C01_002" = "bin.below.5",
                                     "S2503_C01_003" = "bin.5.10",
                                     "S2503_C01_004" = "bin.10.15",
                                     "S2503_C01_005" = "bin.15.20",
                                     "S2503_C01_006" = "bin.20.25",
                                     "S2503_C01_007" = "bin.25.35",
                                     "S2503_C01_008" = "bin.35.50",
                                     "S2503_C01_009" = "bin.50.75",
                                     "S2503_C01_010" = "bin.75.100",
                                     "S2503_C01_011" = "bin.100.150",
                                     "S2503_C01_012" = "bin.150.above",
                                     "S2503_C01_013" = "median")) %>%
            filter(grepl("bin", variable)) %>%
            mutate(estimate = estimate / 100,
                   moe = moe / 100,
                   se = moe / 1.645)
    } else {
        bin_ests_new = bin_ests_new %>%
            mutate(variable = recode(variable,
                                     "S2503_C01_001" = "total",
                                     "S2503_C02_002" = "bin.below.5",
                                     "S2503_C02_003" = "bin.5.10",
                                     "S2503_C02_004" = "bin.10.15",
                                     "S2503_C02_005" = "bin.15.20",
                                     "S2503_C02_006" = "bin.20.25",
                                     "S2503_C02_007" = "bin.25.35",
                                     "S2503_C02_008" = "bin.35.50",
                                     "S2503_C02_009" = "bin.50.75",
                                     "S2503_C02_010" = "bin.75.100",
                                     "S2503_C02_011" = "bin.100.150",
                                     "S2503_C02_012" = "bin.150.above",
                                     "S2503_C02_013" = "median")) %>%
            filter(grepl("bin", variable)) %>%
            mutate(estimate = estimate / 100,
                   moe = moe / 100,
                   se = moe / 1.645)
    }

    bin_ests_wide = bin_ests_new %>% select(-moe, -se) %>%
        pivot_wider(names_from = variable, values_from = estimate) %>%
        select(-bin.150.above)
    bin_ses_wide = bin_ests_new %>% select(-moe, -estimate) %>%
        pivot_wider(names_from = variable, values_from = se) %>%
        select(-bin.150.above)

    mean_ests_wide = mean_ests_new %>% select(-moe, -se) %>%
        pivot_wider(names_from = variable, values_from = estimate) %>%
        select(GEOID, NAME, bin.150.200, bin.200.above)
    mean_ses_wide = mean_ests_new %>% select(-moe, -estimate) %>%
        pivot_wider(names_from = variable, values_from = se) %>%
        select(GEOID, NAME, bin.150.200, bin.200.above)

    bounds1 = c(5, 10, 15, 20, 25, 35, 50, 75, 100, 150)
    bounds2 = c(10, 15, 25, 35, 50, 75, 100, 150, 200)
    bounds = unique(c(bounds1, bounds2))
    lbs = c("below", bounds)
    ubs = c(bounds, "above")

    full_ests = full_join(mean_ests_wide, bin_ests_wide) %>%
        select(GEOID, NAME, paste("bin", lbs, ubs, sep="."))

    full_ses = full_join(mean_ses_wide, bin_ses_wide) %>%
        select(GEOID, NAME, paste("bin", lbs, ubs, sep="."))

    us_proportions = full_ests %>%
        select(-GEOID, -NAME) %>%
        as.matrix %>%
        as.vector

    us_proportions
}

create_us_race_bin_proportions = function(year){
    bin_ests_white = get_acs(geography = "US",
                             table = "B19001A",
                             cache_table = TRUE,
                             year = year) %>%
        mutate(variable = recode(variable,
                                 "B19001A_001" = "total.white",
                                 "B19001A_002" = "bin.below.10",
                                 "B19001A_003" = "bin.10.15",
                                 "B19001A_004" = "bin.15.20",
                                 "B19001A_005" = "bin.20.25",
                                 "B19001A_006" = "bin.25.30",
                                 "B19001A_007" = "bin.30.35",
                                 "B19001A_008" = "bin.35.40",
                                 "B19001A_009" = "bin.40.45",
                                 "B19001A_010" = "bin.45.50",
                                 "B19001A_011" = "bin.50.60",
                                 "B19001A_012" = "bin.60.75",
                                 "B19001A_013" = "bin.75.100",
                                 "B19001A_014" = "bin.100.125",
                                 "B19001A_015" = "bin.125.150",
                                 "B19001A_016" = "bin.150.200",
                                 "B19001A_017" = "bin.200.above")) %>%
        mutate(se = moe / 1.645)

    bin_ests_black = get_acs(geography = "US",
                             table = "B19001B",
                             cache_table = TRUE,
                             year = year) %>%
        mutate(variable = recode(variable,
                                 "B19001B_001" = "total.black",
                                 "B19001B_002" = "bin.below.10",
                                 "B19001B_003" = "bin.10.15",
                                 "B19001B_004" = "bin.15.20",
                                 "B19001B_005" = "bin.20.25",
                                 "B19001B_006" = "bin.25.30",
                                 "B19001B_007" = "bin.30.35",
                                 "B19001B_008" = "bin.35.40",
                                 "B19001B_009" = "bin.40.45",
                                 "B19001B_010" = "bin.45.50",
                                 "B19001B_011" = "bin.50.60",
                                 "B19001B_012" = "bin.60.75",
                                 "B19001B_013" = "bin.75.100",
                                 "B19001B_014" = "bin.100.125",
                                 "B19001B_015" = "bin.125.150",
                                 "B19001B_016" = "bin.150.200",
                                 "B19001B_017" = "bin.200.above")) %>%
        mutate(se = moe / 1.645)

    white_prop = bin_ests_white %>% select(-moe, -se) %>%
        pivot_wider(names_from = variable,
                    values_from = estimate) %>%
        arrange(GEOID) %>%
        mutate_at(vars(matches("bin")), list(~ ./total.white)) %>%
        select(-total.white, -GEOID, -NAME) %>%
        as.matrix %>%
        as.vector

    black_prop = bin_ests_black %>% select(-moe, -se) %>%
        pivot_wider(names_from = variable,
                    values_from = estimate) %>%
        arrange(GEOID) %>%
        mutate_at(vars(matches("bin")), list(~ ./total.black)) %>%
        select(-total.black, -GEOID, -NAME) %>%
        as.matrix %>%
        as.vector

    out = list()
    out[["white"]] = white_prop
    out[["black"]] = black_prop
    out
}


create_puma_standats = function(states, pumas, years, full_state_datasets){

    pknot_prior_loc = list()
    for(curr_year in years){
        yearString = paste("y", curr_year, sep="")
        pknot_prior_loc[[yearString]] =
            create_us_bin_proportions(curr_year)
    }

    puma_standats = list()

    for(i in 1:length(states)){
        curr_state = states[i]
        curr_puma = pumas[i]

        puma_standats[[curr_state]] = list()

        for(curr_year in years){
            yearString = paste("y", curr_year, sep="")
            puma_standats[[curr_state]][[yearString]] = list()

            state_puma_data = full_state_datasets %>%
                filter(state == curr_state,
                       as.numeric(puma) == curr_puma,
                       year == curr_year)
            state_puma_data$geometry = NULL

            puma_standats[[curr_state]][[yearString]]$ntract =
                nrow(state_puma_data)

            ## total population
            total_ests = state_puma_data %>%
                select(contains("total")) %>%
                select(contains("est"))

            total_ses = state_puma_data %>%
                select(contains("total")) %>%
                select(contains("se"))

            puma_standats[[curr_state]][[yearString]]$total_est =
                as.matrix(total_ests)

            puma_standats[[curr_state]][[yearString]]$total_se =
                as.matrix(total_ses)

            ## bin estimates and knots
            bin_ests = state_puma_data %>%
                select(contains("bin")) %>%
                select(contains("est"))

            bin_ses = state_puma_data %>%
                select(contains("bin")) %>%
                select(contains("se"))

            knots = names(bin_ests) %>%
                str_split("\\.", simplify = TRUE) %>%
                .[-1,2] %>%
                as.numeric %>%
                c(0, .)

            puma_standats[[curr_state]][[yearString]]$nbin =
                length(knots)
            puma_standats[[curr_state]][[yearString]]$knots =
                knots
            puma_standats[[curr_state]][[yearString]]$bin_est =
                as.matrix(bin_ests)
            puma_standats[[curr_state]][[yearString]]$bin_se =
                as.matrix(bin_ses)

            ## quantile estimates
            quant_ests = state_puma_data %>%
                select(contains("quant")) %>%
                select(contains("est"))

            quant_ses = state_puma_data %>%
                select(contains("quant")) %>%
                select(contains("se"))

            quantiles = names(quant_ests) %>%
                str_split("\\.", simplify = TRUE) %>%
                .[,2] %>%
                as.numeric %>%
                (function(x){x/100})

            puma_standats[[curr_state]][[yearString]]$nquant =
                ncol(quant_ests)
            puma_standats[[curr_state]][[yearString]]$quantiles =
                quantiles
            puma_standats[[curr_state]][[yearString]]$quant_est =
                as.matrix(quant_ests)
            puma_standats[[curr_state]][[yearString]]$quant_se =
                as.matrix(quant_ses)

            ## share estimates
            share_ests = state_puma_data %>%
                select(contains("share")) %>%
                select(contains("est"))

            share_ses = state_puma_data %>%
                select(contains("share")) %>%
                select(contains("se"))

            share_lb = cbind(0, quant_ests)
            share_ub = quant_ests %>%
                select(-quant.95.est) %>%
                cbind(., -1, -1)

            puma_standats[[curr_state]][[yearString]]$nshare =
                ncol(share_ests)
            puma_standats[[curr_state]][[yearString]]$share_est =
                as.matrix(share_ests)
            puma_standats[[curr_state]][[yearString]]$share_se =
                as.matrix(share_ses)
            puma_standats[[curr_state]][[yearString]]$share_lb =
                as.matrix(share_lb)
            puma_standats[[curr_state]][[yearString]]$share_ub =
                as.matrix(share_ub)

            ## mean ests
            puma_standats[[curr_state]][[yearString]]$mean_est =
                state_puma_data %>%
                select(contains("mean")) %>%
                select(contains("est")) %>%
                as.matrix

            puma_standats[[curr_state]][[yearString]]$mean_se =
                state_puma_data %>%
                select(contains("mean")) %>%
                select(contains("se")) %>%
                as.matrix

            ## median ests
            puma_standats[[curr_state]][[yearString]]$median_est =
                state_puma_data %>%
                select(contains("median")) %>%
                select(contains("est")) %>%
                as.matrix

            puma_standats[[curr_state]][[yearString]]$median_se =
                state_puma_data %>%
                select(contains("median")) %>%
                select(contains("se")) %>%
                as.matrix

            ## gini ests
            puma_standats[[curr_state]][[yearString]]$gini_est =
                state_puma_data %>%
                select(contains("gini")) %>%
                select(contains("est")) %>%
                as.matrix

            puma_standats[[curr_state]][[yearString]]$gini_se =
                state_puma_data %>%
                select(contains("gini")) %>%
                select(contains("se")) %>%
                as.matrix

            ## prior location for pknot_tract
            puma_standats[[curr_state]][[yearString]]$pknot_prior_loc =
                pknot_prior_loc[[yearString]]
        }
    }
    puma_standats
}


create_metro_standats = function(full_metro_datasets){

    years = unique(full_metro_datasets$year)

    pknot_prior_loc = list()
    for(curr_year in years){
        yearString = paste("y", curr_year, sep="")
        pknot_prior_loc[[yearString]] =
            create_us_bin_proportions(curr_year)
    }

    metros = unique(full_metro_datasets$metro)
    metro_standats = list()

    for(curr_metro in metros){
        metro_standats[[curr_metro]] = list()
        for(curr_year in years){
            yearString = paste("y", curr_year, sep="")
            metro_standats[[curr_metro]][[yearString]] = list()

            metro_year_data = full_metro_datasets %>%
                filter(metro == curr_metro,
                       year == curr_year)

            ## total population
            total_ests = metro_year_data %>%
                select(contains("total")) %>%
                select(contains("est"))

            total_ses = metro_year_data %>%
                select(contains("total")) %>%
                select(contains("se"))

            metro_standats[[curr_metro]][[yearString]]$ntract =
                nrow(metro_year_data)

            metro_standats[[curr_metro]][[yearString]]$total_est =
                as.matrix(total_ests)

            metro_standats[[curr_metro]][[yearString]]$total_se =
                as.matrix(total_ses)

            ## bin estimates and knots
            bin_ests = metro_year_data %>%
                select(contains("bin")) %>%
                select(contains("est"))

            bin_ses = metro_year_data %>%
                select(contains("bin")) %>%
                select(contains("se"))

            knots = names(bin_ests) %>%
                str_split("\\.", simplify = TRUE) %>%
                .[-1,2] %>%
                as.numeric %>%
                c(0, .)

            metro_standats[[curr_metro]][[yearString]]$nbin =
                length(knots)
            metro_standats[[curr_metro]][[yearString]]$knots =
                knots
            metro_standats[[curr_metro]][[yearString]]$bin_est =
                as.matrix(bin_ests)
            metro_standats[[curr_metro]][[yearString]]$bin_se =
                as.matrix(bin_ses)

            ## quantile estimates
            quant_ests = metro_year_data %>%
                select(contains("quant")) %>%
                select(contains("est"))

            quant_ses = metro_year_data %>%
                select(contains("quant")) %>%
                select(contains("se"))

            quantiles = names(quant_ests) %>%
                str_split("\\.", simplify = TRUE) %>%
                .[,2] %>%
                as.numeric %>%
                (function(x){x/100})

            metro_standats[[curr_metro]][[yearString]]$nquant =
                ncol(quant_ests)
            metro_standats[[curr_metro]][[yearString]]$quantiles =
                quantiles
            metro_standats[[curr_metro]][[yearString]]$quant_est =
                as.matrix(quant_ests)
            metro_standats[[curr_metro]][[yearString]]$quant_se =
                as.matrix(quant_ses)

            ## share estimates
            share_ests = metro_year_data %>%
                select(contains("share")) %>%
                select(contains("est"))

            share_ses = metro_year_data %>%
                select(contains("share")) %>%
                select(contains("se"))

            share_lb = cbind(0, quant_ests)
            share_ub = quant_ests %>%
                select(-quant.95.est) %>%
                cbind(., -1, -1)

            metro_standats[[curr_metro]][[yearString]]$nshare =
                ncol(share_ests)
            metro_standats[[curr_metro]][[yearString]]$share_est =
                as.matrix(share_ests)
            metro_standats[[curr_metro]][[yearString]]$share_se =
                as.matrix(share_ses)
            metro_standats[[curr_metro]][[yearString]]$share_lb =
                as.matrix(share_lb)
            metro_standats[[curr_metro]][[yearString]]$share_ub =
                as.matrix(share_ub)

            ## mean ests
            metro_standats[[curr_metro]][[yearString]]$mean_est =
                metro_year_data %>%
                select(contains("mean")) %>%
                select(contains("est")) %>%
                as.matrix

            metro_standats[[curr_metro]][[yearString]]$mean_se =
                metro_year_data %>%
                select(contains("mean")) %>%
                select(contains("se")) %>%
                as.matrix

            ## median ests
            metro_standats[[curr_metro]][[yearString]]$median_est =
                metro_year_data %>%
                select(contains("median")) %>%
                select(contains("est")) %>%
                as.matrix

            metro_standats[[curr_metro]][[yearString]]$median_se =
                metro_year_data %>%
                select(contains("median")) %>%
                select(contains("se")) %>%
                as.matrix

            ## gini ests
            metro_standats[[curr_metro]][[yearString]]$gini_est =
                metro_year_data %>%
                select(contains("gini")) %>%
                select(contains("est")) %>%
                as.matrix

            metro_standats[[curr_metro]][[yearString]]$gini_se =
                metro_year_data %>%
                select(contains("gini")) %>%
                select(contains("se")) %>%
                as.matrix

            ## prior location for pknot_tract
            metro_standats[[curr_metro]][[yearString]]$pknot_prior_loc =
                pknot_prior_loc[[yearString]]
        }
    }
    metro_standats
}


create_metro_race_standats = function(full_metro_race_datasets){

    years = unique(full_metro_race_datasets$year)

    pknot_prior_loc = list()
    for(curr_year in years){
        yearString = paste("y", curr_year, sep="")
        pknot_prior_loc[[yearString]] =
            create_us_race_bin_proportions(curr_year)
    }

    metros = unique(full_metro_race_datasets$metro)
    metro_race_standats = list()

    for(curr_metro in metros){
        metro_race_standats[[curr_metro]] = list()
        for(curr_year in years){
            yearString = paste("y", curr_year, sep="")
            metro_race_standats[[curr_metro]][[yearString]] = list()
            metro_race_standats[[curr_metro]][[yearString]][["white"]] = list()
            metro_race_standats[[curr_metro]][[yearString]][["black"]] = list()
            metro_year_data = full_metro_race_datasets %>%
                filter(metro == curr_metro, year == curr_year)

            ## total population
            total_ests_white = metro_year_data %>%
                select(contains("white")) %>%
                select(contains("total")) %>%
                select(contains("est"))

            total_ses_white = metro_year_data %>%
                select(contains("white")) %>%
                select(contains("total")) %>%
                select(contains("se"))

            total_ests_black = metro_year_data %>%
                select(contains("black")) %>%
                select(contains("total")) %>%
                select(contains("est"))

            total_ses_black = metro_year_data %>%
                select(contains("black")) %>%
                select(contains("total")) %>%
                select(contains("se"))

            metro_race_standats[[curr_metro]][[yearString]][["white"]]$ntract =
                nrow(metro_year_data)
            metro_race_standats[[curr_metro]][[yearString]][["black"]]$ntract =
                nrow(metro_year_data)

            metro_race_standats[[curr_metro]][[yearString]][["white"]]$total_est =
                as.matrix(total_ests_white)
            metro_race_standats[[curr_metro]][[yearString]][["white"]]$total_se =
                as.matrix(total_ses_white)

            metro_race_standats[[curr_metro]][[yearString]][["black"]]$total_est =
                as.matrix(total_ests_black)
            metro_race_standats[[curr_metro]][[yearString]][["black"]]$total_se =
                as.matrix(total_ses_black)


            ## bin estimates and knots
            bin_ests_white = metro_year_data %>%
                select(contains("white")) %>%
                select(contains("bin")) %>%
                select(contains("est"))

            bin_ses_white = metro_year_data %>%
                select(contains("white")) %>%
                select(contains("bin")) %>%
                select(contains("se"))

            bin_ests_black = metro_year_data %>%
                select(contains("black")) %>%
                select(contains("bin")) %>%
                select(contains("est"))

            bin_ses_black = metro_year_data %>%
                select(contains("black")) %>%
                select(contains("bin")) %>%
                select(contains("se"))

            knots = names(bin_ests_white) %>%
                str_split("\\.", simplify = TRUE) %>%
                .[-1,2] %>%
                as.numeric %>%
                c(0, .)

            metro_race_standats[[curr_metro]][[yearString]][["white"]]$nbin =
                length(knots)
            metro_race_standats[[curr_metro]][[yearString]][["white"]]$knots =
                knots
            metro_race_standats[[curr_metro]][[yearString]][["white"]]$bin_est =
                as.matrix(bin_ests_white)
            metro_race_standats[[curr_metro]][[yearString]][["white"]]$bin_se =
                as.matrix(bin_ses_white)

            metro_race_standats[[curr_metro]][[yearString]][["black"]]$nbin =
                length(knots)
            metro_race_standats[[curr_metro]][[yearString]][["black"]]$knots =
                knots
            metro_race_standats[[curr_metro]][[yearString]][["black"]]$bin_est =
                as.matrix(bin_ests_black)
            metro_race_standats[[curr_metro]][[yearString]][["black"]]$bin_se =
                as.matrix(bin_ses_black)

            ## mean ests
            metro_race_standats[[curr_metro]][[yearString]][["white"]]$mean_est =
                metro_year_data %>%
                select(contains("white")) %>%
                select(contains("mean")) %>%
                select(contains("est")) %>%
                as.matrix

            metro_race_standats[[curr_metro]][[yearString]][["white"]]$mean_se =
                metro_year_data %>%
                select(contains("white")) %>%
                select(contains("mean")) %>%
                select(contains("se")) %>%
                as.matrix

            metro_race_standats[[curr_metro]][[yearString]][["black"]]$mean_est =
                metro_year_data %>%
                select(contains("black")) %>%
                select(contains("mean")) %>%
                select(contains("est")) %>%
                as.matrix

            metro_race_standats[[curr_metro]][[yearString]][["black"]]$mean_se =
                metro_year_data %>%
                select(contains("black")) %>%
                select(contains("mean")) %>%
                select(contains("se")) %>%
                as.matrix

            ## median ests
            metro_race_standats[[curr_metro]][[yearString]][["white"]]$median_est =
                metro_year_data %>%
                select(contains("white")) %>%
                select(contains("median")) %>%
                select(contains("est")) %>%
                as.matrix

            metro_race_standats[[curr_metro]][[yearString]][["white"]]$median_se =
                metro_year_data %>%
                select(contains("white")) %>%
                select(contains("median")) %>%
                select(contains("se")) %>%
                as.matrix

            metro_race_standats[[curr_metro]][[yearString]][["black"]]$median_est =
                metro_year_data %>%
                select(contains("black")) %>%
                select(contains("median")) %>%
                select(contains("est")) %>%
                as.matrix

            metro_race_standats[[curr_metro]][[yearString]][["black"]]$median_se =
                metro_year_data %>%
                select(contains("black")) %>%
                select(contains("median")) %>%
                select(contains("se")) %>%
                as.matrix

            ## prior location for pknot_tract
            metro_race_standats[[curr_metro]][[yearString]][["white"]]$pknot_prior_loc =
                pknot_prior_loc[[yearString]][["white"]]
            metro_race_standats[[curr_metro]][[yearString]][["black"]]$pknot_prior_loc =
                pknot_prior_loc[[yearString]][["black"]]

        }
    }
    metro_race_standats
}
