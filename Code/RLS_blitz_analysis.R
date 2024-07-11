# Code to look at all the spatial data together
# June 2, 2023
# Em Lim

# Load packages and functions ----
library(tidyverse)
library(visreg)
library(ggplot2)
library(ggeffects)
library(lubridate)
library(vegan) # for diversity indexes
library(MuMIn) # for dredge?
library(glmmTMB) # better for random effects?
library(patchwork)
library(emmeans)
# for R2 value
library(insight)
library(performance)
#library(Matrix)

# package version control
library(renv)

# Source pretty functions
source("Code/Functions.R") # Length to weight function here!

# Set default plotting
theme_set(theme_bw())

# RLS data: Load data + manipulate  ------
# RLS data from the website!
fish_invert <- read_csv("Data/RLS/RLS_data/reef_fish_abundance_and_biomass.csv",
                        show_col_types = FALSE) %>%
  rbind(read_csv("Data/RLS/RLS_data/cryptobenthic_fish_abundance_sizes_added.csv",
                 show_col_types = FALSE)) %>%
  rbind(read_csv("Data/RLS/RLS_data/mobile_macroinvertebrate_abundance.csv",
                 show_col_types = FALSE)) %>%
  as.data.frame() %>%
  clean_phylo_names() %>% # function to fix naming errors
  mutate(survey_date = ymd(survey_date),
         year = as.factor(year(survey_date)),
         site_code = ifelse(site_name == "Swiss Boy", "BMSC24", site_code),
         # correct for rectangle area
         survey_den = case_when(method == 1 ~ total/500,
                                method == 2 ~ total/100)) %>%
  filter(month(survey_date) == 4 | month(survey_date) == 5) %>% # Just the RLS blitz data for now
  mutate(size_class = case_when(
    species_name == "Rhinogobiops nicholsii" & size_class == 0 ~ 7.5,
    species_name == "Artedius harringtoni" & size_class == 0 ~ 5,
    species_name == "Jordania zonope" & size_class == 0 ~ 5,
    species_name == "Oxylebius pictus" & size_class == 0 ~ 12.5,
    TRUE ~ as.numeric(size_class))) %>%
  filter(species_name != "Myoxocephalus aenaeus") %>% # Remove east coast fish
  filter(species_name != "Phoca vitulina") %>% # Remove seal, lets just focus on inverts and fish
  # remove fish seen on wrong method
  filter(!(family == "Gobiidae" & method == 1)) %>% 
  filter(!(family == "Cottidae" & method == 1)) %>% 
  filter(!(family == "Hexagrammidae" & method == 2)) 

# Rhinogobiops nicholsii 7.5 avg
# Artedius harringtoni avg 5
# Jordania zonope avg 5


# RLS Biomass calculations -------

# A note on length to weight relationships
# The overall relationship is W = a*L^b
# log form of the above formula is log(W) = log(a) + b*log(L)
# log(W) = log(a) + b*log(L) is the same as log10(W) = log10(a) + b*log10(L) so don't worry about the base as long as it's the same across the formula

# Use the a and b parameters from FishBase
# W = exp(log(a) + b*log(L))
# W = a*L^b
# Cite Froese R, Thorson JT, Reyes Jr RB. A Bayesian approach for estimating length‐weight relationships in fishes. Journal of Applied Ichthyology. 2014 Feb;30(1):78-85.

# Black rockfish in a paper ranged from 35.4 to 92.5 mm standard length (SL) and 1.10 to 17.78 g wet weight

# Use WL relationships from fishbase to estimate weight of fish in RLS surveys
rls <- fish_invert %>%
  length_to_weight() %>% # Use nice length to weight function!
  # careful, the function shrinks the big wolf eel
  home_range() %>% # calculate each fish's home range function
  mutate(biomass_per_indiv = biomass/survey_den) # see how the RLS biomass calc estimated each fish size
#  mutate(size_class_sum_g = weight_size_class_sum*1000)


# That one huge wolf eel can't be right
# Fishbase: max size = 240 cm, max weight = 18.4 kg
# Ours is apparently 187 cm and 63 kg
# That formula must be off
# I also have size data for inverts "Haliotis kamtschatkana", "Crassadoma gigantea", "Pycnopodia helianthoides", "Polyorchis penicillatus", "Bolinopsis infundibulum", Pleurobrachia bachei, Pleuronichthys coenosus

# Some people counted M2 fishes on M1, but not always so I don't actually want goby and sculpin counts from M1
# Figure out what to do here

# Megathura crenulata probably isn't found here????

#write_csv(rls, "Output/Output_data/rls_species.csv")

# extract just one row per survey to join with the pee data and tide data
rls_survey_info <- rls %>%
  select(site_code, year, depth, survey_id, survey_date, hour) %>%
  rename(survey_depth = depth) %>%
  unique() %>%
  mutate(date_time = ymd_hms(paste(survey_date, hour)),
         date = survey_date)

# Extract phylo data! 
# write_csv(rls %>% 
#   select(phylum, class, order, family, species_name) %>% 
#   unique(), "Output/Output_data/rls_phylo.csv")


# Pivot the data wider for diversity metrics
# currently the df I want to tack the metrics onto is sorted by date, oldest at top
# Spreading long data for taxa abundance to columns for each species
rls_wider <- rls %>% 
  dplyr::select(survey_id, species_name, survey_den) %>% 
  group_by(survey_id, species_name) %>%
  summarise(survey_den = sum(survey_den)) %>%
  ungroup() %>%
  spread(key = species_name, value = survey_den) %>%
  replace(is.na(.), 0)

# write_csv(rls_wider, "Output/Output_data/rls_wider.csv")

# then calculate biodiversity metrics
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
  rename(site_code = site_ID) %>%
  filter(month == "May") %>%
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

# what's the overall rate of exchange for this period? How to define slack vs ebb vs flood
tide_overall <- tide %>%
  mutate(rate = 100 * (tide_m - lag(tide_m))/lag(tide_m)
  ) %>%
  group_by(date) %>%
  slice(-1) %>%
  ungroup()
# mean = 0.0003725301
# sd = 0.379465
# max = +/- 1.511335

# tide cat:
# looked at rate of change over whole 6 week period, SD = 0.379465
# 1/2 SD would make slack period one SD wide
# that makes slack 1 hour long at low tide (ebb to flood) but 2 hours at high tide
# but not all high and lows bc some swing more than others

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
  left_join(rls %>%
              filter(phylum == "Chordata") %>%
              group_by(survey_id) %>%
              summarize(fish_weight_sum = sum(weight_size_class_sum),
                        fish_weight_weighted = sum(weight_weighted),
                        fish_abund = sum(survey_den))) %>%
  # invert biomass
  left_join(rls %>%
              filter(phylum != "Chordata") %>%
              group_by(survey_id) %>%
              summarize(invert_weight_sum = sum(weight_size_class_sum),
                        invert_weight_weighted = sum(weight_weighted),
                        invert_abund = sum(survey_den))) %>%
  # rls survey richness and abundance and total biomass
  left_join(rls %>%
              group_by(survey_id) %>%
              summarize(weight_sum = sum(weight_size_class_sum),
                        all_weight_weighted = sum(weight_weighted),
                        species_richness = n_distinct(species_name),
                        abundance = sum(survey_den))
  ) %>%
  # diversity indexes
  left_join(rls_wide, by = "survey_id") %>%
  # tide exchange 
  left_join(tide_exchange, by = "survey_id") %>%
  # scale and center variables here
  mutate(
    # the biomass + abundance variables
    weight_sum_scale = c(scale(weight_sum)),
    fish_weight_sum_scale = c(scale(fish_weight_sum)),
    fish_weight_weighted_scale = c(scale(fish_weight_weighted)),
    fish_abund_scale = c(scale(fish_abund)),
    invert_weight_sum_scale = c(scale(invert_weight_sum)),
    invert_weight_weighted_scale = c(scale(invert_weight_weighted)),
    invert_abund_scale = c(scale(invert_abund)),
    all_weighted_scale = c(scale(all_weight_weighted)),
    abundance_scale = c(scale(abundance)),
    # biodiversity variables
    rich_scale = c(scale(species_richness)),
    shannon_scale = c(scale(shannon)),
    simpson_scale = c(scale(simpson)),
    # abiotic variables I should control for
    depth_avg_scale = c(scale(depth_avg)),
    tide_scale = c(scale(avg_exchange_rate)),
    tide_cat = factor(as.factor(case_when(avg_exchange_rate < -0.1897325 ~ "Ebb",
                                          avg_exchange_rate < 0.1897325 ~ "Slack",
                                          avg_exchange_rate > 0.1897325 ~ "Flood")),
                      levels = c("Ebb", "Slack", "Flood")),
    # center instead of scale
    abundance_center = c(scale(abundance, scale = FALSE)),
    tide_center = c(scale(avg_exchange_rate, scale = FALSE)),
    shannon_center = c(scale(shannon, scale = FALSE)),
    depth_center = c(scale(depth_avg, scale = FALSE))
  )


#write_csv(rls_final, "Output/Output_data/rls_final.csv")

# can i make a flood, ebb, slack var for exchange?
#max(rls_final$avg_exchange_rate)
#min(rls_final$avg_exchange_rate)
#hist(rls_final$avg_exchange_rate)

# spread/3 = 0.6034917
# flood = 0.6195273 - 1.223019
# slack = 0.0160356 - 0.6195273
# ebb = -0.5874561 - 0.0160356

# but we know slack is actually centered around 0 soooo
# what's the min value/3
# that divides the min - 0 into three parts, two parts = ebb, one part on either side of 0 = slack
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

# Build full model with gaussian distribution
# mod_norm_full <- glmmTMB(nh4_avg ~ weight_sum_scale*tide_scale + shannon_scale + depth_avg_scale + (1|year) + (1|site_code), 
#                          family = 'gaussian',
#                          data = rls_final)
# plot(simulateResiduals(mod_norm_full)) 


# Build full model with gamma distribution
# mod_gam_full <- glmmTMB(nh4_avg ~ weight_sum_scale*tide_scale + shannon_scale + depth_avg_scale + (1|year) + (1|site_code), 
#                         family = Gamma(link = 'log'),
#                         data = rls_final)
# plot(simulateResiduals(mod_gam_full)) 
# 
# compare via AIC
# AIC(mod_norm_full, mod_gam_full) 

# Gamma residuals look better and that mod has a lower AIC!

# Step 2: Choose your predictors

# Biological predictors:
# 1) Biomass
# a) weight_sum = total wet weight of all the animals observed on RLS transects
# b) all_weight_weighted = weighted fish weight + invert weight
# Weight sum and weight weighted are barely different
# c) fish_weight_sum = total wet weight of just the fishes (kg)
# d) fish_weight_weighted = total weight of fishes weighted by home range size
# e) abundance
# Abundance is preferred by AIC!

# 2) Biodiversity
# a) species_richness = total # species on each transect
# b) shannon = shannon diversity of each transect
# c) simpson = simpson diversity of each transect
# Richness is slightly preferred by AIC but I'm going with shannon

# These are the variables I really care about!!!! The abiotics are the things I might need to control for

# Continuous abiotic predictors
# 1) depth_avg = average nh4 sample depths
# 2) avg_exchange_rate = average rate of change of the tide height over the 1 hour survey
# Depth doesn't matter but avg_exchange_rate is important!

# Categorical abiotic predictors
# 3) year = 2021, 2022, or 2023
# 4) site_code = each unique site
# should year and site be fixed or random effects?
# interested trends across a broad spatial scale (and not site-specific trends) so site should be random
# so year should also be random???

# What if year is absent or fixed
# Year = fixed no change, years are diff from 2021
# Drop year = interaction still signif but species richness is now signif negative slope
# What if site is absent or fixed
# Site = fixed no change, some sites are diff from each other
# Drop site no change
# What if both are fixed or dropped?
# Both fixed = no change.
# Both dropped = interaction still there but species richness is now signif negative slope

# Interactions
# I think abundance:tide will have an interaction
# I double checked, there's no triple interaction or abundance:richness interaction! Phewwww


# I ran through dredging this model:
mod_all <- glmmTMB(nh4_avg ~ abundance_scale*tide_scale + shannon_scale + depth_avg_scale + (1|year) + (1|site_code), 
                   family = Gamma(link = 'log'),
                   data = rls_final,
                   na.action = na.fail)
summary(mod_all)
plot(DHARMa::simulateResiduals(mod_all))

# dredge <- as.data.frame(dredge(mod_all)) %>% filter(delta < 3)

# Best according to AIC: just intercept, lol
# Second best = abundance_scale*tide_scale
# Third best = just tide

# so the best AIC model
mod_aic <- glmmTMB(nh4_avg ~ abundance_scale*tide_scale + (1|year) + (1|site_code), 
                   family = Gamma(link = 'log'),
                   data = rls_final)
summary(mod_aic)
plot(DHARMa::simulateResiduals(mod_aic))

# I prefer biomass over abundance bc that's theoretically more linked to NH4
# But abundance is a straight up recorded measure, no biomass proxy BS
# And I prefer shannon diversity over richness bc it's a "better" biodiversity measure
# And I want to drop depth bc it doesn't seem to matter...

mod_brain <- glmmTMB(nh4_avg ~ abundance_scale*tide_scale + shannon_scale + depth_avg_scale + (1|year) + (1|site_code), 
                     family = Gamma(link = 'log'),
                     data = rls_final)
summary(mod_brain)
plot(DHARMa::simulateResiduals(mod_brain))

# just look at fish weights
mod_fish <- glmmTMB(nh4_avg ~ fish_weight_sum_scale*tide_scale + shannon_scale + depth_avg_scale + (1|year) + (1|site_code), 
                    family = Gamma(link = 'log'),
                    data = rls_final)
summary(mod_fish)
plot(DHARMa::simulateResiduals(mod_fish)) # ok so neither fish abundance or fish biomass significantly predict nh4

# just look at invert weights
mod_invert <- glmmTMB(nh4_avg ~ invert_abund_scale*tide_scale + shannon_scale + depth_avg_scale + (1|year) + (1|site_code), 
                      family = Gamma(link = 'log'),
                      data = rls_final)
summary(mod_invert)
plot(DHARMa::simulateResiduals(mod_invert)) # invert abundance has a negative abundance:tide interaction and the model output looks the same as the total abundance model


# weight instead of abundance
mod_weight <- glmmTMB(nh4_avg ~ weight_sum_scale*tide_scale + shannon_scale + depth_avg_scale + (1|year) + (1|site_code), 
                      family = Gamma(link = 'log'),
                      data = rls_final)
summary(mod_weight)

# simpson instead of shannon
mod_simp <- glmmTMB(nh4_avg ~ abundance_scale*tide_scale + simpson_scale + depth_avg_scale + (1|year) + (1|site_code), 
                    family = Gamma(link = 'log'),
                    data = rls_final)

mod_richness <- glmmTMB(nh4_avg ~ abundance_scale*tide_scale + rich_scale + depth_avg_scale + (1|year) + (1|site_code), 
                        family = Gamma(link = 'log'),
                        data = rls_final)

summary(mod_richness)

mod_simp_weight <- glmmTMB(nh4_avg ~ weight_sum_scale*tide_scale + simpson_scale + depth_avg_scale + (1|year) + (1|site_code), 
                           family = Gamma(link = 'log'),
                           data = rls_final)


# table for publication
AIC_tab_rls <- AIC(mod_brain, mod_weight, mod_simp, mod_simp_weight) %>%
  rownames_to_column() %>%
  mutate(best = min(AIC),
         delta = round((AIC - best), digits = 2),
         likelihood = exp( -0.5*delta),
         sum = sum(likelihood),
         AICw = round((likelihood/sum), digits = 2),
         AIC = round(AIC, digits = 2)) %>%
  select(rowname, df, AIC, delta, AICw)

# Use AIC to get predictors, then show the coefficients and a model output

# centered variables instead of scaled for estimates
mod_brain_c <- glmmTMB(nh4_avg ~ abundance_center*tide_center + shannon_center + depth_center + (1|year) + (1|site_code), 
                       family = Gamma(link = 'log'),
                       data = rls_final)

summary(mod_brain_c)


# Step 4: Check for collinearity of predictors

# car can't handle random effects so make a simplified mod
car::vif(lm(nh4_avg ~ abundance_scale + tide_scale + shannon_scale + depth_avg_scale, data = rls_final))
# All good, shannon is a little high


# Family data manipulation ------

# Make a dataframe where each family we saw on a survey gets a row
# If we didn't see a certain family on that survey, it doesn't get a row
# No zeros for abundance!
family_df_no0_a <- rls_final %>%
  select(site, site_code, survey_id, year, nh4_avg, depth_avg_scale, weight_sum_scale, abundance_scale, shannon_scale, tide_scale, tide_cat) %>%
  left_join((rls %>%
               mutate(family = as.factor(family))%>%
               group_by(survey_id, method, phylum, family) %>% # if i want methods split up, add it back here
               summarise(fam_den = sum(survey_den),
                         weight_fam_sum_g = 1000*sum(weight_size_class_sum))),
            by = "survey_id") 


# what are the top families?

# df for just the surveys I'm looking at
rls_sm <- rls_final %>%
  select(survey_id) %>%
  unique() %>%
  left_join(rls, by = "survey_id")

# CHOOSE THE FAMILIES TO INCLUDE
# rank of families by total abundance (density)
fam_list_total <- rls_sm %>%
  group_by(family) %>%
  summarise(sum = sum(survey_den)) %>%
  drop_na(family) %>%
  arrange(desc(sum)) %>%
  transmute(family = family, 
            sum_total = sum,
            rank_total = 1:58)

# but which families show up on the most transects?
fam_list_count <- rls_sm %>%
  select(survey_id, family) %>%
  unique() %>%
  count(family) %>%
  drop_na(family) %>%
  arrange(desc(n)) %>%
  transmute(family = family, 
            sum_count = n,
            rank_count = 1:58)

# rank of families by biomass
fam_list_bio <- rls_sm %>%
  group_by(family) %>%
  summarise(sum = sum(weight_size_class_sum)) %>%
  drop_na(family) %>%
  arrange(desc(sum)) %>%
  transmute(family = family, 
            sum_bio = sum,
            rank_bio = 1:58)

# join all that fam stuff together
fam_big_list <- fam_list_total %>%
  left_join(fam_list_count, by = "family") %>%
  left_join(fam_list_bio, by = "family") %>%
  mutate(rank_all = (rank_total + rank_count + rank_bio)/3)

# OK so how do I make sense of this
# plot all rankings together
ggplot(fam_big_list) +
  geom_point(aes(y = reorder(family, -rank_all), 
                 x = rank_total, colour = "density"), alpha = 0.5) +
  geom_point(aes(y = reorder(family, -rank_all), 
                 x = rank_count, colour = "count"), alpha = 0.5) +
  geom_point(aes(y = reorder(family, -rank_all), 
                 x = rank_bio, colour = "biomass"), alpha = 0.5) +
  geom_point(aes(y = reorder(family, -rank_all), 
                 x = rank_all, colour = "all"), alpha = 0.5) +
  geom_vline(xintercept = 15, lty = "dashed") +
  geom_hline(yintercept = "Sebastidae", lty = "dashed") +
  labs(y = "Family", x = "Rank", colour = "Rank Type")

# just the families I want to keep
fam_list_cut <- fam_big_list %>%
  arrange(rank_all) %>%
  head(15) %>%
  select(family) 

# OK and now I want to make a final df where the top 15 families are named, and everything else is "other"

# make this final df for the df without 0's
family_df_no0 <- family_df_no0_a %>%
  filter(family %in% fam_list_cut$family)
# I no longer want to do this!!!!
# now make all the other families "other"
# rbind(family_df_no0_a %>%
#         filter(!family %in% fam_list_cut$family) %>%
#         mutate(family = "other"))

# since I'm working with each family separately I want to change density back to abundance
fam_df_no0 <- family_df_no0 %>%
  mutate(total_fam = case_when(method == 1 ~ fam_den*500,
                               method == 2 ~ fam_den*100),
         weight_abund_fam_sum_g = case_when(method == 1 ~ weight_fam_sum_g*500,
                                            method == 2 ~ weight_fam_sum_g*100),
         abund_fam_scale = c(scale(total_fam)),
         weight_abund_fam_scale = c(scale(weight_abund_fam_sum_g)),
         fam_den_scale = c(scale(fam_den)),
         weight_den_fam_scale = c(scale(weight_fam_sum_g))) %>%
  as.data.frame() %>%
  droplevels()

# mess around with orders next time!!!!!! 

# rank of orders by abundance
order_list_total <- rls_sm %>%
  group_by(order) %>%
  summarise(sum = sum(total)) %>%
  drop_na(order) %>%
  arrange(desc(sum))

# rank orders by biomass
order_list_bio <- rls_sm %>%
  group_by(order) %>%
  summarise(sum = sum(weight_size_class_sum)) %>%
  drop_na(order) %>%
  arrange(desc(sum))

# Family stats ----

# I want to run one model for each family
# list family names
fam_levels <- levels(fam_df_no0$family)

# lapply fam_fun over each family name to create df of predictions for each fam mod
fam_predictions <- lapply(fam_levels, function(family_name) {
  fam_fun_combo(fam_df_no0, family_name)  
})

# join up those predictions
rls_fam_predicts <- bind_rows(fam_predictions)


# what are the top 6 inverts and fish?
top_fam_r2 <- rls_fam_predicts %>%
  select(c(family, r2)) %>%
  unique() %>%
  arrange(desc(r2)) %>%
  head(7) %>% # only want top 3 inverts and 3 fish but urchins are above the 3rd fish fam
  filter(family != "Strongylocentrotidae")

# make df of just the tops
rls_top_fam_predict <- rls_fam_predicts %>%
  filter(family %in% top_fam_r2$family) %>%
  # optional reorder of families
  arrange(desc(r2)) %>%
  mutate(family = factor(family, unique(family)))


# Outputs for these families
# list top 6 family names
top_fam_levels <- levels(rls_top_fam_predict$family)

# get model outputs for each family model
lapply(top_fam_levels, function(family_name) {
  diagnose_fun(fam_df_no0, family_name) 
})

# now filter full family df to just include those top 6 families
rls_top_fam <- top_fam_r2 %>%
  left_join(fam_df_no0, by = "family") %>%
  arrange(desc(r2)) %>%
  mutate(family = factor(family, unique(family)))


# Family community analysis -----
# multivariate space for the big fam model??? are communities dominated by certain fams when there's a lot of nh4 in the water?
# make matrix from communities and then see which direction the nh4 loads???
# envfit
# if there's a site with a lots of greenlings, do those sites have higher nh4

# create the wide form family data

# make wide for families
fam_wide <- family_df_no0 %>% 
  dplyr::select(survey_id, family, fam_den) %>% 
  group_by(survey_id, family) %>%
  summarise(total = sum(fam_den)) %>%
  ungroup() %>%
  spread(key = family, value = total) %>%
  replace(is.na(.), 0) %>%
  select(-17) # cut the last column, it's NA

# first use the rls final to filter for included surveys and make sure the order of rows is the same  
com <- rls_final %>%
  left_join(fam_wide, by = "survey_id") %>%
  select(34:48)

# make env data
env <- rls_final %>%
  select(nh4_avg, depth_avg, avg_exchange_rate)

#convert com to a matrix
m_com = as.matrix(com)

# Perform the NMDS ordination
set.seed(123)
nmds = metaMDS(m_com, distance = "bray")
nmds

# Now we run the envfit function with our environmental data frame, env
en = envfit(nmds, env, permutations = 999, na.rm = TRUE)
# The first parameter is the metaMDS object from the NMDS ordination we just performed. Next is env, our environmental data frame. Then we state we want 999 permutations, and to remove any rows with missing data.
en

plot(nmds)
plot(en)

plot(nmds, type = "t")
plot(en)


# my old code
# make the dissim matrix
dissim_mat <- vegdist(fam_wider, method = "horn")
mod <- adonis(dissim_mat ~ nh4_avg, data = rls_final, permutations = 9999)
# this looks at species distribution as function of environmental data
summary(mod)

# Ordination: nMDS
myNMDS <- metaMDS(fam_wider, k = 2)
myNMDS #most important: is the stress low? Here it is 0.22 which is a bit on the high side
stressplot(myNMDS) #low stress means that the observed dissimilarity between site pairs matches that on the 2-D plot fairly well (points hug the line)

my_envfit <- envfit(myNMDS, site_data, permutations = 999)
spp_fit <- envfit(myNMDS, fam_wider, permutations = 999)

# try to plot?

#save NMDS results into dataframe
site_scrs <- as.data.frame(scores(myNMDS, display = "sites")) 

spp_scrs <- as.data.frame(scores(spp_fit, display = "vectors")) #save species intrinsic values into dataframe
spp_scrs <- cbind(spp_scrs, Species = rownames(spp_scrs)) #add species names to dataframe
spp_scrs <- cbind(spp_scrs, pval = spp_fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant

# order species
signif_spp_scrs1 <- spp_scrs[order(spp_scrs$pval),]

#cut species that weren't significant
sig_spp_scrs <- subset(signif_spp_scrs1, pval<0.137) #subset data to show species significant at 0.05

#ugly plot
ordiplot(myNMDS, type = "n") 
ordihull(myNMDS, groups = site_data$beach,draw = "polygon",col = "grey99",label = T)
orditorp(myNMDS, display = "species", col = "purple4",air = 0.01, cex = 0.9)
orditorp(myNMDS, display = "sites", cex = 0.75, air = 0.01)


# Graphing ----
# sort out some palettes!
pal20 <- viridis::viridis(20)
pie(rep(1, 20), col = pal20)

pal30 <- viridis::viridis(30)


pal1 <- pal20[14] # fam plots
#pal3 <- c(pal20[20], pal20[15], pal20[11]) # pred plot
pal3 <- c(pal30[16], pal30[25], pal30[29]) # tide colours

pal5 <- c(pal20[1], pal20[3], pal20[5], pal20[7], pal20[9]) # coeff plot


# Fig 3a: Coefficient plot ----

# generate coefficients
# try using emtrends to generate standard error/ confidence intervals around slope estimates
rls_coeffs <- confint(mod_brain, level = 0.95, method = c("wald"), component = c("all", "cond", "zi", "other"), estimate = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(variable = rowname,
         lower_CI = `2.5 %`,
         upper_CI = `97.5 %`,
         estimate = Estimate) %>%
  head(- 2)  %>% # remove random effects
  mutate(estimate = ifelse(variable == "(Intercept)", exp(estimate), estimate),
         lower_CI = ifelse(variable == "(Intercept)", exp(lower_CI), lower_CI),
         upper_CI = ifelse(variable == "(Intercept)", exp(upper_CI), upper_CI)) %>%
  tail(-1) %>%
  mutate(
    variable = factor(as.factor(variable), 
                      labels = c("Abundance", "Tide", "Biodiversity", "Depth", "Abundance:tide"))) %>%
  arrange(desc(estimate)) %>%
  mutate(variable = factor(variable, unique(variable)))



# Coefficient plot 
rls_coeff_plot <- coeff_plot(coeff_df = rls_coeffs, 
                             pal = pal5) +
  place_label("(a)")


# ggsave("Output/Figures/rls_mod_coeff.png", device = "png", height = 9, width = 12, dpi = 400)


# Fig 3b: Abundance vs nh4 -----
mod_pred_plot <- glmmTMB(nh4_avg ~ abundance*tide_scale + shannon_scale + depth_avg_scale + (1|year) + (1|site_code), 
                         family = Gamma(link = 'log'),
                         data = rls_final)

# get mean tide for each level
tide_means <- rls_final %>%
  group_by(tide_cat) %>%
  summarise(tide = mean(tide_scale))

v <- c(as.numeric(tide_means[1,2]), as.numeric(tide_means[2,2]), as.numeric(tide_means[3,2]))

# ggpredict
predict <- ggpredict(mod_pred_plot, terms = c("abundance", "tide_scale [v]")) %>% 
  as.data.frame() %>%
  mutate(abundance = x,
         tide_cat = factor(as.factor(case_when(group == as.numeric(tide_means[1,2]) ~ "Ebb",
                                               group == as.numeric(tide_means[2,2]) ~ "Slack",
                                               group == as.numeric(tide_means[3,2]) ~ "Flood")),
                           levels = c("Ebb", "Slack", "Flood")))

# now plot these predictions
rls_pred_plot <- 
  plot_pred(raw_data = rls_final,
            predict_data = predict, 
            plot_type = "rls",
            x_var = abundance, y_var = nh4_avg, 
            lty_var = tide_cat,
            pal = pal3,
            theme = "white") +
  place_label("(b)")


# Fig 3c: Family plots -----

fam_plot <- ggplot() + 
  geom_point(data = rls_top_fam, 
             aes(x = fam_den, y = nh4_avg), colour = pal1,
             alpha = 0.8) +
  labs(y = expression(paste("Ammonium ", (mu*M))), 
       x = expression(paste("Animal abundance/m"^2))) +
  facet_wrap(~family, scales = 'free_x', ncol = 3) +
  geom_line(data = rls_top_fam_predict,
            aes(x = fam_den, y = predicted), colour = pal1,
            linewidth = 1) +
  geom_ribbon(data = rls_top_fam_predict,
              aes(x = fam_den, y = predicted, 
                  ymin = conf.low, ymax = conf.high), fill = pal1,
              alpha = 0.15) +
  theme_white() +
  theme(strip.background = element_rect(fill = "grey", color = "grey"))
#  geom_text(
#    data = rls_fam %>% select(family, slope) %>% unique(),
#    mapping = aes(x = Inf, y = Inf, label = slope),
#    hjust   = 1.1,
#    vjust   = 1.5,
#    size = 9)

fam_plot

#ggsave("Output/Pub_figs/Fig3d.png", device = "png", height = 18, width = 12, dpi = 400)

#ggsave("Output/Pub_figs/Fig3d.png", device = "png", height = 9, width = 16, dpi = 400)

# Fig 3 -----
rls_coeff_plot + rls_pred_plot  &
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(hjust = -1, vjust = 0))

# ggsave("Output/Pub_figs/Fig3.png", device = "png", height = 9, width = 16, dpi = 400)


# Alt Fig 3 tri-panel -----
squish <- theme(axis.title.y = element_text(margin = margin(r = -120, unit = "pt")))

(rls_coeff_plot + rls_pred_plot)/(fam_plot + squish) &
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(hjust = -1, vjust = 0))

# ggsave("Output/Pub_figs/Fig3.png", device = "png", height = 18, width = 16, dpi = 400)


# Figs for presentations -----
# just plot each line one by one for the model prediction vs raw data plot
# just ebb
plot_pred(raw_data = rls_final %>% filter(tide_cat == "Ebb"),
          predict_data = predict %>% filter(tide_cat == "Ebb"),
          plot_type = "rls",
          x_var = abundance, y_var = nh4_avg, 
          pal = pal3,
          theme = "black") +
  theme(legend.position = "right") 

#ggsave("Output/Pres_figs/Fig3ba.png", device = "png", height = 9, width = 12, dpi = 400)

# add slack
plot_pred(raw_data = rls_final %>% filter(tide_cat != "Flood"), 
          predict_data = predict %>% filter(tide_cat != "Flood"),
          plot_type = "rls",
          x_var = abundance, y_var = nh4_avg, 
          pal = pal3,
          theme = "black") +
  theme(legend.position = "right") 

#ggsave("Output/Pres_figs/Fig3bb.png", device = "png", height = 9, width = 12, dpi = 400)

# all of it
plot_pred(raw_data = rls_final, 
          predict_data = predict,
          plot_type = "rls",
          x_var = abundance, y_var = nh4_avg, 
          pal = pal3,
          theme = "black") +
  theme(legend.position = "right") 

#ggsave("Output/Pres_figs/Fig3b.png", device = "png", height = 9, width = 12, dpi = 400)


# Data exploration ------

# abundant across all surveys
all_rls <- (rls %>% select(phylum, class, order, family, species_name, total, weight_per_indiv_g)) %>%
  rbind(kelp_rls %>% select(phylum, class, order, family, species_name, total, weight_per_indiv_g))

d <- all_rls %>% count(phylum, class, order, family, species_name, weight_per_indiv_g)
#write_csv(d %>% filter(phylum != "Chordata"), "Output/Output_data/species_list.csv")

e <- all_rls %>% filter(weight_per_indiv_g == 0.5) %>% count(species_name)

a <- all_rls %>% filter(species_name == "Cryptolithodes sitchensis")

all_abund <- all_rls %>%
  group_by(family) %>%
  summarise(total = sum(total))

# can I save a species list from this?

# pycno size
(rls %>% select(survey_date, species_name, size_class, total) %>% filter(species_name == "Pycnopodia helianthoides")) %>%
  rbind(kelp_rls %>% select(Date, species_name, size_class, total) %>% rename(survey_date = Date) %>% filter(species_name == "Pycnopodia helianthoides"))  %>%
  uncount(total) %>%
  mutate(year = year(survey_date)) %>%
  ggplot(aes(year, size_class, colour = year)) +
  geom_jitter(width = 0.2) +
  stat_summary(fun = "mean", geom = "point", size = 5, mapping = aes(group = year)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1.5, mapping = aes(group = year)) +
  labs(x = "Date", y = "Pycno diameter (cm)") +
  theme(legend.position = "null") +
  scale_x_continuous(breaks = c(2021, 2022, 2023))

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
no_sizes <- rls %>%
  filter(size_class == "0") %>%
  filter(species_name != "Bolinopsis infundibulum") %>%
  filter(species_name != "Pleuronichthys coenosus") %>%
  filter(species_name != "Pleurobrachia bachei") %>%
  filter(species_name != "Polyorchis penicillatus") %>%
  select(-year)

# just look at goby sizes
goby <- rls %>%
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

# What's the most abundant species?
rls_abundant <- fish %>% 
  group_by(species_name) %>%
  summarise(total = sum(total))



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


# Are the sites spatially autocorrelated?
resids <- DHARMa::simulateResiduals(mod_brain)
new_resids <- recalculateResiduals(resids, group = rls_final$site_code)

rls_auto <- rls_final %>%
  left_join(read_csv("Data/RLS/RLS_data/true_coords.csv"), by = "site_code") %>%
  select(longitude, latitude) %>%
  unique()

testSpatialAutocorrelation(new_resids, x = rls_auto$longitude, y = rls_auto$latitude)

# Fuck around ----

abalone <- rls %>%
  filter(species_name == "Haliotis kamtschatkana") %>%
  group_by(survey_id) %>%
  summarize(density = sum(total)*100) 


# Graveyard -----
# wendy code for R2 values 
performance::r2(beta_reg_model, tolerance = 0.0000000000001)
performance::r2_nakagawa(beta_reg_model)
tidy(beta_reg_model)
summary(beta_reg_model)
performance::r2(beta_reg_model, tolerance = 0.0000000000001)

# Calculate the marginal R2
marginal_R2 <- 1 - (mod_brain$deviance / mod_brain$null_deviance)
# Calculate the conditional R2
conditional_R2 <- 1 - (mod_brain$deviance / mod_brain$null_deviance)
marginal_R2
conditional_R2

# The big family model

# On surveys that we saw a family, does that family abundance ~ nh4 ???
mod_fam_no0 <- glmmTMB(nh4_avg ~ abund_fam_scale * family +
                         (1|year) + (1|site_code), 
                       family = Gamma(link = 'log'),
                       data = family_df_no0)
summary(mod_fam_no0)
plot(DHARMa::simulateResiduals(mod_fam_no0)) # less bad with the triple,

# can I emtrends to get the slopes
df <- emtrends(mod_fam_no0, pairwise ~ family, var = "abund_fam_scale")$emtrends %>%
  as.data.frame()


ggplot(df, aes(x = abund_fam_scale.trend, y = reorder(family, abund_fam_scale.trend), xmin = asymp.LCL, xmax = asymp.UCL)) +
  geom_point(size = 2.7) +
  geom_errorbar(width = 0, linewidth = 0.5) +
  geom_vline(xintercept=0, color="black", linetype="dashed")


# random effect of family
mod_fam_no0 <- glmmTMB(nh4_avg ~ abund_fam_scale +
                         (1|year) + (1|site_code) + (1+abund_fam_scale|family), 
                       family = Gamma(link = 'log'),
                       data = family_df_no0)

d <- ranef(mod_fam_no0)

# OK
# No random effect of family because the same nh4 value is replicated for each family, so I can't really have a random effect of family too