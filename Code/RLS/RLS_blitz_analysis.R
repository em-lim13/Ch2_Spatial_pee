# Code to look at all the spatial data together
# June 2, 2023
# Em Lim

# Load packages and functions ----
library(tidyverse)
library(visreg)
library(ggplot2)
library(PNWColors)
library(ggeffects)
library(lubridate)
library(sf)
library(lmerTest)
library(vegan) # for diversity indexes
library(MuMIn) # for dredge?
library(TMB)
library(glmmTMB) # better for random effects?

# Source pretty functions
source("Code/theme_black.R")
source("Code/Functions.R") # Length to weight function here!

# Set default plotting
theme_set(theme_bw())

# RLS data: Load data + manipulate  ------
# RLS data from the website!

# Just the pelagic and cryptic fish
fish <- read_csv("Data/RLS/RLS_data/reef_fish_abundance_and_biomass.csv",
                 show_col_types = FALSE) %>%
  rbind(read_csv("Data/RLS/RLS_data/cryptobenthic_fish_abundance_sizes_added.csv",
                 show_col_types = FALSE)) %>%
  as.data.frame() %>%
  mutate(survey_date = ymd(survey_date),
         year = as.factor(year(survey_date)),
         site_code = ifelse(site_name == "Swiss Boy", "BMSC24", site_code)) %>%
  filter(month(survey_date) == 4 | month(survey_date) == 5) %>% # Just the RLS blitz data for now
  mutate(size_class = case_when(
    species_name == "Rhinogobiops nicholsii" & size_class == 0 ~ 7.5,
    species_name == "Artedius harringtoni" & size_class == 0 ~ 5,
    species_name == "Jordania zonope" & size_class == 0 ~ 5,
    species_name == "Oxylebius pictus" & size_class == 0 ~ 12.5,
    TRUE ~ as.numeric(size_class)))

# Rhinogobiops nicholsii 7.5 avg
# Artedius harringtoni avg 5
# Jordania zonope avg 5

# Just the mobile inverts
invert <- read_csv("Data/RLS/RLS_data/mobile_macroinvertebrate_abundance.csv",
                   show_col_types = FALSE) %>%
  as.data.frame() %>%
  mutate(survey_date = ymd(survey_date),
         year = as.factor(year(survey_date)),
         site_code = ifelse(site_name == "Swiss Boy", "BMSC24", site_code),
         species_name = case_when(
           species_name == "Montereina nobilis" ~ "Peltodoris nobilis",
           species_name == "Parastichopus californicus" ~ "Apostichopus californicus",
           species_name == "Berthella californica" ~ "Berthella chacei",
           species_name == "Henricia leviuscula" ~ "Henricia spp.",
                        TRUE ~ as.character(species_name))) %>%
  filter(month(survey_date) == 4 | month(survey_date) == 5) # Just the RLS blitz data for now


# RLS Biomass calculations -------

# A note on length to weight relationships
# The overall relationship is W = a*L^b
# log form of the above formula is log(W) = log(a) + b*log(L)
# log(W) = log(a) + b*log(L) is the same as log10(W) = log10(a) + b*log10(L) so don't worry about the base as long as it's the same across the formula

# Use the a and b parameters from FishBase
# W = exp(log(a) + b*log(L))
# Cite Froese R, Thorson JT, Reyes Jr RB. A Bayesian approach for estimating length‐weight relationships in fishes. Journal of Applied Ichthyology. 2014 Feb;30(1):78-85.

# Use WL relationships from fishbase to estimate weight of fish in RLS surveys
fishes <- fish %>%
  length_to_weight() %>% # Use nice length to weight function!
          # careful, the function shrinks the big wolf eel
  filter(species_name != "Bolinopsis infundibulum") %>%
  filter(species_name != "Pleuronichthys coenosus") %>%
  filter(species_name != "Pleurobrachia bachei") %>%
  filter(species_name != "Polyorchis penicillatus") %>% # Filter inverts
  filter(species_name != "Phoca vitulina") %>% # Remove seal
  filter(species_name != "Actinopterygii spp.") %>% # Remove unidentif fish
  filter(species_name != "Myoxocephalus aenaeus") %>% # Remove east coast fish
  home_range() %>% # calculate each fish's home range
  mutate(biomass_per_indiv = biomass/total, # see how the RLS biomass calc estimated each fish size
         range_weight = 1/(range*100), # scale smallest range to a weight of 1, everything else is fraction
         weight_weighted = weight_size_class_sum*range_weight) # weight weights by range

  
# That one huge wolf eel can't be right
# Fishbase: max size = 240 cm, max weight = 18.4 kg
# Ours is apparently 187 cm and 63 kg
# That formula must be off
# I also have size data for inverts "Haliotis kamtschatkana", "Crassadoma gigantea", "Pycnopodia helianthoides", "Polyorchis penicillatus", "Bolinopsis infundibulum", Pleurobrachia bachei, Pleuronichthys coenosus


# Some people counted M2 fishes on M1, but not always so I don't actually want goby and sculpin counts from M1
# Figure out what to do here

# now calc invert weights
inverts <- invert %>%
  invert_length_to_weight() %>%
  mutate(biomass_per_indiv = biomass/total, # see how the RLS biomass calc estimated each fish size
         range = 0,
         range_weight = 1, 
         weight_weighted = weight_size_class_sum*range_weight) 
    # For now the invert ranges are just NAs, because I'm pretty sure they'll all be really small ranges and also the biomass estimates are already so vague


# Join all rls data together
rls <- rbind(fishes, inverts)

# extract just one row per survey to join with the pee data and tide data
rls_survey_info <- rls %>%
  select(site_code, year, depth, survey_id, survey_date, hour) %>%
  rename(survey_depth = depth) %>%
  unique() %>%
  mutate(date_time = ymd_hms(paste(survey_date, hour)))


# Pivot the data wider for diversity metrics
# currently the df I want to tack the metrics onto is sorted by date, oldest at top
# Spreading long data for taxa abundance to columns for each species

rls_wider <- rls %>% # Spreading the data to wide format
  dplyr::select(survey_id, species_name, total) %>% 
  group_by(survey_id, species_name) %>%
  summarise(total = sum(total)) %>%
  ungroup() %>%
  # Removing any columns not needed
  spread(key = species_name, value = total) %>%
  replace(is.na(.), 0)


rls_wide <- rls_wider %>%
  mutate(shannon = (diversity((rls_wider %>% select(-survey_id)), index = "shannon")),
         simpson = (diversity((rls_wider %>% select(-survey_id)), index = "simpson"))) %>%
  select(survey_id, shannon, simpson)


# NH+ data: Load + manipulate  ----

# RLS nh4 from 2021
# This is the spring blitz, the june samples, and the July samples
# The may samples have estimated matrix effects, not great
# The June and July samples had a own matrix spike from each site but were still compared to DI

# RLS from 2022
# These had a own matrix spike from each site but were still compared to DI
# They came out quite negative maybe because of SFU DI, and maybe because of temperature difference between the standards and samples
# I took the lowest nh4 reading and added it to everything else to "set" that sample to 0 and bump everything up
# Maybe an underestimation

# RLS from 2023
# Did the full "proper" Taylor protocol with standard bottles + BF from each site

# combine these three years into one!
rls_nh4 <- rbind(read_csv("Output/Output_data/RLS_nh4_2021.csv"),
                        read_csv("Output/Output_data/RLS_nh4_2022.csv"),
                        read_csv("Output/Output_data/RLS_nh4_2023.csv")) %>%
  mutate(year = as.factor(year)) %>%
  filter(month == "May") %>% 
  rename(site_code = site_ID) %>%
  left_join(rls_survey_info, by = c("site_code", "year")) %>%
  depth_function() # only keep the RLS survey from the transect where the pee is from


# Tide exchange: Load + manipulate data ----

# Downloaded tide height data (m) from http://tbone.biol.sc.edu/tide/tideshow.cgi?site=Bamfield%2C+British+Columbia in 1 min intervals over the two week spanning each RLS spring blitz, in 24 hour time
# 2021 April 26 start
# added May 20, 2021 for the two Dixon sites
# 2022 April 25 start
# 2023 May 8 start

tide <- read_csv("Data/RLS/tides_1_min.csv") %>%
  mutate(date_time = ymd_hms(paste(date, time)))

# build an empty dataframe
# Do this each time you run the for loop!!!
tide_exchange <- data.frame()

# then write the loop
for (x in 1:nrow(rls_survey_info)) {
  
  survey_start <- ymd_hms(rls_survey_info$date_time[x:x])
  survey_end <- survey_start + hours(1)
  
  output = tide %>%
    filter(between(date_time, survey_start, survey_end)) %>%
    mutate(rate = 100 * (tide_m - lag(tide_m))/lag(tide_m),
           max = rate[which.max(abs(rate))]) %>%
    slice(-1) %>%
    summarise(avg_exchange_rate = mean(rate),
              max_exchange_rate = mean(max)
              ) %>%
    mutate(survey_id = rls_survey_info$survey_id[x:x])
  
  tide_exchange = rbind(tide_exchange, output)
}
# Now I have the average and max rate of change of tide height for each survey!!!
# no diff in models between avg and max tide exchange, just use avg


# Site level averaging -----
# One row per transect
  # Per site, per depth, per year
  # Same as grouping by survey_ID
  # NH4+, fish biomass, richness, abundance, and tide exchange for each survey

rls_final <- 
  # average nh4 for each site
  rls_nh4 %>% 
  group_by(survey_id) %>%
  mutate(nh4_avg = mean(nh4_conc),
         temp_avg = mean(temp_est),
         depth_avg = mean(depth),
         year = as.factor(year)) %>%
  ungroup() %>% 
  select(c(site, site_code, survey_id, date_time, year, nh4_avg, temp_avg, depth_avg)) %>%
  unique() %>%
  # rls fish biomass 
  left_join(fishes %>%
              group_by(survey_id) %>%
              summarize(fish_weight_sum = sum(weight_size_class_sum),
                        fish_weight_weighted = sum(weight_weighted))) %>%
  # rls survey richness and abundance and total biomass
  left_join(rls %>%
              group_by(survey_id) %>%
              summarize(weight_sum = sum(weight_size_class_sum),
                        all_weight_weighted = sum(weight_weighted),
                        species_richness = n_distinct(species_name),
                        abundance = sum(total))) %>%
  # diversity indexes
  left_join(rls_wide, by = "survey_id") %>%
  # tide exchange 
  left_join(tide_exchange, by = "survey_id") %>%
  # scale and centre variables here
  mutate(
    # the biomass + abundance variables
    weight_sum_stand = c(scale(weight_sum)),
    fish_weight_sum_stand = c(scale(fish_weight_sum)),
    fish_weight_weighted_stand = c(scale(fish_weight_weighted)),
    all_weight_weighted_stand = c(scale(all_weight_weighted)),
    abundance_stand = c(scale(abundance)),
    # biodiversity variables
    species_richness_stand = c(scale(species_richness)),
    shannon_stand = c(scale(shannon)),
    simpson_stand = c(scale(simpson)),
    # abiotic variables I should control for
    depth_avg_stand = c(scale(depth_avg)),
    avg_exchange_rate_stand = c(scale(avg_exchange_rate)),
    tide_cat = case_when(avg_exchange_rate < -0.1958187 ~ "Ebb",
                         avg_exchange_rate < 0.1958187 ~ "Slack",
                         avg_exchange_rate > 0.1958187 ~ "Flood")
  )

# can i make a flood, ebb, slack var for exchange?
max(rls_final$avg_exchange_rate)
min(rls_final$avg_exchange_rate)
hist(rls_final$avg_exchange_rate)

# spread/3 = 0.6034917
# flood = 0.6195273 - 1.223019
# slack = 0.0160356 - 0.6195273
# ebb = -0.5874561 - 0.0160356

# but we know slack is actually centered around 0 soooo
# ebb = < -0.1958187 
# slack = -0.1958187 - 0.1958187
# ebb = > 0.1958187


# Stats -------

# Step 0: Ask my question
  # Is there a link between biodiversity or biomass and ammonium concentration?

# Response variable: nh4_avg
  # this is the mean of 3 measurements, how do I incorporate that error in the models?

# Step 1: Choose a distribution
ggplot(rls_final, aes(x = nh4_avg)) +
  geom_histogram(bins = 30) 

# NH4 is positive only and continuous (not intergers) 
# Gamma is likely the best bet for this!

# Step 2: Choose your predictors

# Biological predictors:
# 1) Biomass
    # a) weight_sum = total wet weight of all the animals observed on RLS transects
    # b) all_weight_weighted = weighted fish weight + invert weight
        # Weight sum and weight weighted are barely different
    # c) fish_weight_sum = total wet weight of just the fishes (kg)
    # d) fish_weight_weighted = total weight of fishes weighted by home range size
# I'm going with weight_sum

# 2) Biodiversity
    # a) species_richness = total # species on each transect
    # b) shannon = shannon diversity of each transect
    # c) simpson = simpson diversity of each transect
# I'm going with shannon
    
# These are the variables I really care about!!!! The abiotics are the things I might need to control for

# Continuous abiotic predictors
  # 1) depth_avg = average nh4 sample depths
  # 2) avg_exchange_rate = average rate of change of the tide height over the 1 hour survey

# Categorical abiotic predictors
  # 3) year = 2021, 2022, or 2023
  # 4) site_code = each unique site
    # should year and site be fixed or random effects?
      # interested trends across a broad spatial scale (and not site-specific trends) so site should be random
      # so year should also be random???

# Interactions
  # I think whatever biomass and biodiversity variable I choose should have an interaction with exchange rate

# Build full model with gaussian distribution
mod_norm_full <- glmmTMB(nh4_avg ~ weight_sum_stand*avg_exchange_rate_stand + shannon_stand*avg_exchange_rate_stand + depth_avg_stand + (1|year) + (1|site_code), 
                         family = 'gaussian',
                         data = rls_final)
summary(mod_norm_full)

# Build full model with gamma distribution
mod_gam_full <- glmmTMB(nh4_avg ~ weight_sum_stand*avg_exchange_rate_stand + shannon_stand*avg_exchange_rate_stand + depth_avg_stand + (1|year) + (1|site_code), 
                         family = Gamma(link = 'log'),
                         data = rls_final)
summary(mod_gam_full)

# Step 3: Model residuals 
  # Compare the two distrubutions 

# Gaussian
mod_norm_full_resids <- simulateResiduals(mod_norm_full)
plot(mod_norm_full_resids) 

# Gamma
mod_norm_gam_resids <- simulateResiduals(mod_gam_full)
plot(mod_norm_gam_resids) 

AIC(mod_norm_full, mod_gam_full)

# Gamma residuals look better and that mod has a lower AIC!

# Step 4: Check for collinearity of predictors

### THIS IS WHERE THE WHEELS FALL OFF THE BUS

# car can't handle random effects so make a simplified mod
car::vif(lm(nh4_avg ~ weight_sum_stand + avg_exchange_rate_stand + shannon_stand + depth_avg_stand, data = rls_final))
# Allll goooood

# Double checking my variables!!!

# Would a diff biomass variable would be better than weight sum!

# all weight summed
mod_weight_sum <- glmmTMB(nh4_avg ~ weight_sum_stand*avg_exchange_rate_stand + shannon_stand*avg_exchange_rate_stand + depth_avg_stand + (1|year) + (1|site_code), 
                          family = Gamma(link = 'log'),
                          data = rls_final)

# all_weight_weighted
mod_all_weighted_sum <- glmmTMB(nh4_avg ~ all_weight_weighted_stand*avg_exchange_rate_stand + shannon_stand*avg_exchange_rate_stand + depth_avg_stand + (1|year) + (1|site_code), 
                          family = Gamma(link = 'log'),
                          data = rls_final)

# fish weight
mod_fish_weight <- glmmTMB(nh4_avg ~ fish_weight_sum_stand*avg_exchange_rate_stand + shannon_stand*avg_exchange_rate_stand + depth_avg_stand + (1|year) + (1|site_code), 
                           family = Gamma(link = 'log'),
                           data = rls_final)

# fish weight weighted
mod_fish_weighted <- glmmTMB(nh4_avg ~ fish_weight_weighted_stand*avg_exchange_rate_stand + shannon_stand*avg_exchange_rate_stand + depth_avg_stand + (1|year) + (1|site_code), 
                           family = Gamma(link = 'log'),
                           data = rls_final)


# abundance
mod_abund <- glmmTMB(nh4_avg ~ abundance_stand*avg_exchange_rate_stand + shannon_stand*avg_exchange_rate_stand + depth_avg_stand + (1|year) + (1|site_code), 
                     family = Gamma(link = 'log'),
                     data = rls_final)
summary(mod_abund)

# compare
aic_biomass <- AIC(mod_weight_sum, mod_all_weighted, mod_fish_weight, mod_fish_weighted, mod_abund) %>%
  mutate(delta = AIC - min(AIC)) # abundance is best, then weight sum 


# Which richness variable explains data best

# richness
mod_rich <- glmmTMB(nh4_avg ~ abundance_stand*avg_exchange_rate_stand + species_richness_stand*avg_exchange_rate_stand + depth_avg_stand + (1|year) + (1|site_code), 
                    family = Gamma(link = 'log'),
                    data = rls_final)

# shannon
mod_shan <- glmmTMB(nh4_avg ~ abundance_stand*avg_exchange_rate_stand + shannon_stand*avg_exchange_rate_stand + depth_avg_stand + (1|year) + (1|site_code), 
                    family = Gamma(link = 'log'),
                    data = rls_final)

# simpsons
mod_simp <- glmmTMB(nh4_avg ~ abundance_stand*avg_exchange_rate_stand + simpson_stand*avg_exchange_rate_stand + depth_avg_stand + (1|year) + (1|site_code), 
                    family = Gamma(link = 'log'),
                    data = rls_final)

AIC(mod_rich, mod_shan, mod_simp) # richness is best


# So the best mod would be:
mod_best <- glmmTMB(nh4_avg ~ abundance_stand*avg_exchange_rate_stand + species_richness_stand*avg_exchange_rate_stand + depth_avg_stand + (1|year) + (1|site_code), 
                    family = Gamma(link = 'log'),
                    data = rls_final)
summary(mod_best)
# I double checked, there's no triple interaction or abundance:richness interaction! Phewwww

# Double check the resids on this
mod_best_resids <- simulateResiduals(mod_best)
plot(mod_best_resids) # YEP that's fine!


# The model I liked better was
mod_weight_sum <- glmmTMB(nh4_avg ~ weight_sum_stand*avg_exchange_rate_stand + shannon_stand*avg_exchange_rate_stand + depth_avg_stand + (1|year) + (1|site_code), 
                          family = Gamma(link = 'log'),
                          data = rls_final)
summary(mod_weight_sum)

# But let's go with mod best!!! 


# quick plot to figure out this interaction
abund_plot <- ggplot(rls_final, aes(abundance, nh4_avg, colour = tide_cat)) +
  geom_point() +
  geom_smooth(method = lm) +
  theme(legend.position = "NULL")
# So nh4 increases a bit with abundance at slack and ebb tide, but at flood the trend flips. neat!

weight_plot <- ggplot(rls_final, aes(weight_sum, nh4_avg, colour = tide_cat)) +
  geom_point() +
  geom_smooth(method = lm)

abund_plot + weight_plot

# does this mean weight and abundance are the same?
ggplot(rls_final, aes(abundance, weight_sum)) +
  geom_point() +
  geom_smooth(method = lm) # WOOOAH they basically are... 
# that's a cool result in itself....


# Old stats -----
# Does pee vary by site and year?
simple_model <- lm(nh4_conc ~ site_code + year, data = rls_nh4)
summary(simple_model)
visreg(simple_model)


# Put all predictors in a model
mod_full <- lmer(nh4_avg ~ scale(weight_sum) * scale(shannon) * scale(abundance) * scale(avg_exchange_rate) + (1|site_code), rls_final)
summary(mod_full)

# Put all predictors in a model + depth
# Cut the triple + interactions
mod_fuller <- lmer(nh4_avg ~ scale(weight_sum) + scale(simpson) + scale(abundance) + scale(avg_exchange_rate) + scale(depth_avg) +
                   scale(weight_sum):scale(simpson) + scale(weight_sum):scale(abundance) + 
                   scale(weight_sum):scale(avg_exchange_rate) + scale(weight_sum):scale(depth_avg) +
                   scale(simpson):scale(abundance) + scale(simpson):scale(avg_exchange_rate) + scale(simpson):scale(depth_avg) +
                   scale(abundance):scale(avg_exchange_rate) + scale(abundance):scale(depth_avg) +
                   scale(avg_exchange_rate):scale(depth_avg) +
                   (1|site_code), rls_final)
summary(mod_fuller)


# dredge to compare all models
options(na.action = "na.fail")
dred <- dredge(mod_fuller)
# best model is just the intercept, LOL
# 2nd best (<delta AIC 2) is just richness
# 3rd is just weight, but that's > delta AIC 5
# top 3 models and delta AIC is the same when full mod has all interactions vs when it has none

# When I add depth and only the 2-way interactions top mod = intercept, 2nd = richness (delta 0.73), 3rd is just depth (delta 5.6)

# look at this "best mod"
mod_rich <- lmer(nh4_avg ~ species_richness + (1|site_code), rls_final)
summary(mod_rich)
visreg(mod_rich)
# Negative relationship between species richness and NH4+...... cool cool cool cool cool



# What if I put allll the data together?????
mod_all <- lm(nh4_avg ~ scale(weight_sum) + scale(avg_exchange_rate)  + scale(BiomassM) + 
                scale(weight_sum):avg_exchange_rate + scale(weight_sum):scale(BiomassM), big_rls)

summary(kelp_mod_best)
visreg(kelp_mod_best)


#Data exploration ------
  
  # Data checks
  goby <- fish %>%
  filter(species_name == "Rhinogobiops nicholsii") %>%
  filter(site_name == "Goby Town") %>%
  filter(depth == "5.7") %>%
  select(site_name, depth, method, block, species_name, total, size_class)

# New nudis?
nudi <- invert %>%
  filter(class == "Gastropoda") %>%
  filter(family != "Acmaeidae") %>%
  filter(family != "Calliostomatidae") %>%
  filter(family != "Epitoniidae" ) %>%
  filter(family != "Fissurellidae" ) %>%
  filter(family != "Haliotidae" ) %>%
  filter(family != "Lottiidae" ) %>%
  filter(family != "Muricidae" ) %>%
  filter(family != "Naticidae" ) %>%
  filter(family != "Tegulidae" ) %>%
  filter(family != "Turbinidae") %>%
  count(species_name, year)

# New RLS species?
rls_new <- rls %>%
  count(species_name, year) %>%
  select(-n)

prev <- rls_new %>%
  filter(year != "2023") %>%
  mutate(seen_before = "yes") %>%
  select(-year) %>%
  unique()

new <- rls_new %>%
  filter(year == "2023") %>%
  mutate(seen_this_year = "yes") %>%
  select(-year) %>%
  left_join(prev, by = "species_name") %>%
  replace_na(list(seen_before = "no")) %>%
  filter(seen_before == "no") 

#write_csv(new, "Output/Output_data/new_species.csv")

# missing size data
no_sizes <- fishes %>%
  filter(size_class == "0") %>%
  filter(species_name != "Bolinopsis infundibulum") %>%
  filter(species_name != "Pleuronichthys coenosus") %>%
  filter(species_name != "Pleurobrachia bachei") %>%
  filter(species_name != "Polyorchis penicillatus") %>%
  select(-year)

# just look at goby sizes
goby <- fishes %>%
  filter(species_name == "Oxylebius pictus") %>%
  summarize(mean = mean(size_class))

hist(goby$size_class) 

# Gobies 7.5 avg
# Artedius harringtoni avg 5
# Jordania zonope avg 5
# Oxylebius pictus avg 12.5

  
#write_csv(no_sizes, "Output/Output_data/missing_fish_sizes.csv")

# What's the most abundant species?
rls_abundant_species <- rls %>% 
  count(species_name)

# where are the most urchins???
# Can swap for any species
rls_urchins <- rls %>%
  filter(species_name == "Mesocentrotus franciscanus") %>%
  group_by(site_name) %>%
  summarize(urchins = sum(total)) %>%
  arrange(desc(urchins))

# So basically you can manipulate the data to summarize whatever thing you're interested in, and then you can join that summarized data with another df of interest by the common site code/site name :)

# rank each site by the year
# I'm sure I can do this a better way but I'll just keep this for now
rank <- rls_final %>%
  # 2012
  filter(year(date_time) == "2021") %>%
  group_by(site) %>%
  summarise(nh4_avg2021 = mean(nh4_avg)) %>%
  arrange(desc(nh4_avg2021)) %>% 
  mutate(rank2021 = 1:22,
         grade2021 = 100 - (100*rank2021/22)) %>% # 2021
  #2022
  left_join(rls_final %>%
              filter(year(date_time) == "2022") %>%
              group_by(site) %>%
              summarise(nh4_avg2022 = mean(nh4_avg)) %>%
              arrange(desc(nh4_avg2022)) %>% 
              mutate(rank2022 = 1:19,
                     grade2022 = 100 - (100*rank2022/19)), by= "site") %>%
  # 2023
  left_join(rls_final %>%
              filter(year(date_time) == "2023") %>%
              group_by(site) %>%
              summarise(nh4_avg2023 = mean(nh4_avg)) %>%
              arrange(desc(nh4_avg2023)) %>% 
              mutate(rank2023 = 1:20,
                     grade2023 = 100 - (100*rank2023/20)), by= "site") %>%
  select(grade2021, grade2022, grade2023, site) %>%
  rowwise() %>%
  mutate(avg_grade = mean(c(grade2021, grade2022, grade2023), na.rm = TRUE))


# What were the matrix effects like in 2022?
rls_nh4 %>%
  filter(year == "2022") %>%
  summarise(matrix = mean(matrix))
# 14.6

# how about the end of 2021?
rbind(read_csv("Output/Output_data/RLS_nh4_2021.csv"),
      read_csv("Output/Output_data/RLS_nh4_2022.csv"),
      read_csv("Output/Output_data/RLS_nh4_2023.csv")) %>%
  mutate(year = as.factor(year)) %>%
  filter(year == "2021") %>%
  filter(month != "May") %>%
  summarise(matrix = mean(matrix))
# 7.82
# Go back to the 2021 data and set the ME to 7.82


# Correlation analysis -----
# Make a DF where we only have overlap between sites
b_2021 <- rls_nh4 %>%
  filter(year == "2021") %>%
  transmute(site = site,
            nh4_2021 = nh4_avg) %>%
  unique()

b_2022 <- rls_nh4 %>%
  filter(year == "2022") %>%
  transmute(site = site,
            nh4_2022 = nh4_avg) %>%
  unique()

b_2023 <- rls_nh4 %>%
  filter(year == "2023") %>%
  transmute(site = site,
            nh4_2023 = nh4_avg) %>%
  unique()

# throw them together and drop the NAs
cor_data <- b_2021 %>%
  left_join(b_2022, by = "site") %>%
  left_join(b_2023, by = "site") %>%
  drop_na()

# OK so now do the correlations!

# 2021 vs 2022
cor.test(cor_data$nh4_2021, cor_data$nh4_2022, 
         method = "spearman")
# Yes correlated?

# 2022 vs 2023
cor.test(cor_data$nh4_2022, cor_data$nh4_2023, 
         method = "spearman")
# Yes correlated?

# 2023 vs 2021
cor.test(cor_data$nh4_2023, cor_data$nh4_2021, 
         method = "spearman")
# Yes correlated?

# Graphing ----
pal <- pnw_palette("Sailboat", 3)

rls_nh4 %>%
  group_by(site) %>%
  mutate(nh4_avg = mean(nh4_conc)) %>%
  ungroup() %>%
  ggplot() +
  geom_point(aes(x = reorder(site, -nh4_avg), 
                 nh4_conc, colour = year, fill = year, pch = month),
             size = 3, alpha = 0.75) +
  geom_point(aes(x = reorder(site, -nh4_avg), 
                 nh4_avg),
             size = 3, colour = "black", pch = 20) +
  theme(axis.text.x = element_text(angle = 75, vjust = 0.5)) +
  labs(x = "Site", y = "NH4+ Concentration (umol/L)") +
  scale_colour_manual(values = pal) 


# Plot the ranking of each site year to year
# Plot rank
ggplot(data = rank) +
  geom_point(aes(reorder(site, -avg_grade), grade2021, colour = "2021"), size = 3) +
  geom_point(aes(reorder(site, -avg_grade), grade2022, colour = "2022"), size = 3) +
  geom_point(aes(reorder(site, -avg_grade), grade2023, colour = "2023"), size = 3) +
  labs(x= "Site", y = "Rank", colour = "Year") +
  scale_colour_manual(values = pal) + 
  theme(axis.text.x = element_text(angle = 75, vjust = 0.5)) 

# ggsave("Output/Figures/rank_all_years.png", device = "png",
#        height = 5, width = 8, dpi = 400)


# Plot biomass vs nh4
ggplot(rls_final, aes(weight_sum, nh4_avg, label = site)) +
  geom_point(aes(colour = site_code)) +
  geom_smooth(method = lm) +
  geom_text(check_overlap = TRUE, hjust = 1,  size = 3, ) +
  labs(x = "Total fish biomass (g)", y = "Ammonium concentration (umol)") +
  theme_classic() +
  theme(legend.position = "null")

# richness
ggplot(rls_final, aes(species_richness, nh4_avg, label = site)) +
  geom_point(aes(colour = site_code)) +
  geom_smooth(method = lm) +
  geom_text(check_overlap = TRUE, hjust = 1,  size = 3, ) +
  labs(x = "Richness", y = "Ammonium concentration (umol)") +
  theme_classic() +
  theme(legend.position = "null")

# shannon
ggplot(rls_final, aes(shannon, nh4_avg, label = site)) +
  geom_point(aes(colour = site_code)) +
  geom_smooth(method = lm) +
  geom_text(check_overlap = TRUE, hjust = 1,  size = 3, ) +
  labs(x = "Richness", y = "Ammonium concentration (umol)") +
  theme_classic() +
  theme(legend.position = "null")

# simpson
ggplot(rls_final, aes(simpson, nh4_avg, label = site)) +
  geom_point(aes(colour = site_code)) +
  geom_smooth(method = lm) +
  geom_text(check_overlap = TRUE, hjust = 1,  size = 3, ) +
  labs(x = "Richness", y = "Ammonium concentration (umol)") +
  theme_classic() +
  theme(legend.position = "null")

# abundance
ggplot(rls_final, aes(abundance, nh4_avg, label = site)) +
  geom_point(aes(colour = site_code)) +
  geom_smooth(method = lm) +
  geom_text(check_overlap = TRUE, hjust = 1,  size = 3, ) +
  labs(x = "Abundance", y = "Ammonium concentration (umol)") +
  theme_classic() +
  theme(legend.position = "null")

# tide exchange
ggplot(rls_final, aes(avg_exchange_rate, nh4_avg, label = site)) +
  geom_point(aes(colour = site_code)) +
  geom_smooth(method = lm) +
  geom_text(check_overlap = TRUE, hjust = 1,  size = 3, ) +
  labs(x = "Rate of tidal exchange", y = "Ammonium concentration (umol)") +
  theme_classic() +
  theme(legend.position = "null")

# Site vs pee
ggplot(rls_nh4, aes(nh4_conc, site, colour = year)) +
  geom_boxplot() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15)) +
  labs(x= "Ammonium concentration (umol)", y = "Site") +
  theme_classic()

# black background
ggplot(rls_nh4, aes(nh4_conc, site, fill = site)) +
  geom_boxplot(colour = "white") +
  theme_black() +
  labs(x= "Ammonium concentration (umol/L)", y = "Site") +
  theme(legend.position = "none")

#ggsave("Output/Figures/RLS_sites_pee.png", device = "png",
#       height = 9, width = 16, dpi = 400)



# Visualize temperature
ggplot(rls_nh4, aes(x = temp_est, y = site)) +
  geom_point() +
  labs(x = "Salinity Estimate", y = "Site")

ggplot(bottles, aes(x = temp_est)) +
  geom_histogram(colour="black", fill="white")


# Make a map -----
# Load not great shapefile
potato_map <- sf::st_read("Data/Shapefiles/eez.shp") %>%
  st_sf() %>%
  st_set_crs(4326)

# Define colours
blue <- paste("#b9d1df", sep="")

# Make map without pies, just scaling size of point to %
sf_use_s2(FALSE)

# coords
rls_coords <- rls_final %>%
  left_join(
    rls %>% select(site_code, site_name, latitude, longitude, survey_latitude, survey_longitude,) %>% 
      unique()
    ) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(4326) %>%
  mutate(xjit = ifelse(site_code == "BMSC6" | 
                         site_code == "BMSC3" |
                         site_code == "BMSC21" |
                         site_code == "BMSC15" |
                         site_code == "BMSC11" |
                         site_code == "BMSC19", -0.015, 0.015),
         weight_kg = weight_sum/1000)

## Create legend
#site_names <- as.list(paste(coords$site_num, coords$`Site name`, sep = ": "))

# Make map
ggplot() +
  geom_sf(data = potato_map, fill = blue, colour = "white") +
  geom_sf(data = rls_coords, 
          colour = "black",
          pch = 21,
          alpha = 0.9,
          aes(size = weight_kg,
              fill = nh4_avg)) +
  coord_sf(xlim = c(-125.4, -125.0), ylim = c(48.80, 49), expand = FALSE)  +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 13, colour = "black"),
        legend.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 11, colour = "black"),
        panel.grid.major = element_line(color = "white")) +
  viridis::scale_fill_viridis(option="magma", direction = -1) +
  labs(x = "Longitude", y = "Latitude",
       fill = expression(paste("NH"[4]^" +","(umol)")),
       size = "Fish biomass (kg)")  +
  geom_text(data = rls_coords,
            aes(x = survey_longitude, y = survey_latitude, 
                label = site_code),
            size = 3,
            nudge_x = rls_coords$xjit)

            
#ggsave("Pub_figs/Fig.1.png", device = "png",
#       height = 150, width = 250, units = c("mm"), dpi = 600)


# Graveyard -----

# reduce the df so i can join it with kelp rls data
rls_final_reduced <- rls_final %>%
  mutate(BiomassM = 0) %>%
  select(site, site_code, survey_id, nh4_avg, depth_avg, avg_exchange_rate, BiomassM, weight_sum, species_richness, abundance, avg_exchange_rate,  depth_avg) %>%
  rename(site_code = site_code)

big_rls <- rbind(rls_final_reduced, data_s_reduced) %>%
  mutate(weight_stand = scale(weight_sum),
         rich_stand = scale(species_richness),
         abundance_stand = scale(abundance),
         tide_stand = scale(avg_exchange_rate),
         depth_stand = scale(depth_avg),
         biomass_mean_stand = scale(BiomassM))

# Just look at temp 
temp_df <- rls_final %>%
  drop_na(temp_avg)

mod_temp <- lm(nh4_avg ~ scale(weight_sum) + scale(species_richness) + scale(abundance) + scale(avg_exchange_rate) + scale(temp_avg) + scale(depth_avg), temp_df)
summary(mod_temp)
# Interesting. In 2021, there was a slightly negative relationship between temperature and nh4+ concentration


# Dot and whisker?
sum_stats_pee <- ggpredict(simple_model, terms = c("site_code", "year")) %>% 
  #and then we'll just rename one of the columns so it's easier to plot
  rename(site_code = x,
         nh4_conc = predicted,
         year = group)
# plot
ggplot() +
  geom_point(data = sum_stats_pee, 
             aes(y = site_code, x = nh4_conc, colour = year),
             size = 4) +
  geom_errorbar(data = sum_stats_pee, 
                aes(y = site_code,
                    x = nh4_conc,
                    colour = year,
                    # and you can decide which type of error to show here
                    # we're using 95% CI
                    xmin = conf.low,
                    xmax = conf.high),
                width = 0.2,
                size = 1.2)  +
  geom_point(data = rls_nh4, aes (y = site_code, x = nh4_conc, colour = year), alpha = 0.5, height = 0, size = 2) +
  labs(x= "Ammonium concentration (umol/L)", y = "Site") +
  theme_black() +
  theme(legend.position="none") 

#ggsave("Output/Figures/RLS_pee_black.png", device = "png",
#       height = 9, width = 16, dpi = 400)