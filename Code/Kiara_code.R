# Script to look Kiara's honours
# July 2024
# Em Lim

# Load packages and functions ----
library(tidyverse)
library(visreg)
library(ggplot2)
library(ggeffects)
library(lubridate)
library(vegan) # for diversity indexes
library(glmmTMB) # better for random effects?
library(patchwork)
library(emmeans)
# for R2 value
library(insight)
library(performance)

# Load functions
source("Code/Functions.R")

# Load your nh4 data here!!!! ----

# RLS data ----
kelp_rls <- read_csv("Data/Team_kelp/RLS_KCCA_2022.csv") %>%
  # processing required to get this df into the RLS data format
  filter(Method != 0) %>% # get rid of all method 0's
  slice(2:n()) %>% # cuts the first blank row
  # rename columns
  rename(
    site_code = `Site No.`,
    site_name = `Site Name`, 
    common_name = `Common name`,
    `0` = Inverts,
    species_name = Species,
    method = Method
  )  %>% 
  # Rename columns with spaces
  mutate(species_name = str_to_sentence(species_name),
         common_name = str_to_sentence(common_name),
         Date = dmy(Date),
         date_time_survey = ymd_hms(paste(Date, Time))) %>%
  # use function to fix species naming errors 
  clean_sp_names() %>%
  # Pivot longer for biomass
  pivot_longer( cols = `0`:`400`, names_to = "size_class", values_to = "total") %>%
  drop_na(total) %>%
  filter(total > 0) %>%
  select(-Total) %>%
  # group blocks 1 and 2
  group_by(site_code, site_name, Date, date_time_survey, Depth, method, species_name, common_name, size_class) %>%
  summarise(survey_total = sum(total)) %>% # sum blocks 1 and 2
  ungroup() %>%
  # join phylo names from the rls blitz data
  left_join(read_csv("Output/Output_data/rls_phylo.csv"), by = "species_name") %>%
  clean_phylo_names() %>% # function to fix naming errors
  filter(species_name != "Tonicella spp.") %>% # remove chitons
  filter(species_name != "Cryptochiton stelleri") %>%
  # filter fish seen on wrong survey
  filter(!(family == "Gobiidae" & method == 1)) %>% 
  filter(!(family == "Cottidae" & method == 1)) %>% 
  filter(!(family == "Hexagrammidae" & method == 2)) %>%
  mutate(size_class = as.numeric(size_class),
         # correct for rectangle area
         survey_den = case_when(method == 1 ~ survey_total/500,
                                method == 2 ~ survey_total/100)) %>%
  length_to_weight() %>% # fish length to weight function
  home_range() %>% # calculate each fish's home range function
  as.data.frame() 


# pivot back to wide for biodiversity
kelp_rls_wider <- kelp_rls %>% 
  dplyr::select(site_code, species_name, survey_den) %>% 
  group_by(site_code, species_name) %>%
  summarise(survey_den = sum(survey_den)) %>%
  ungroup() %>%
  spread(key = species_name, value = survey_den) %>%
  replace(is.na(.), 0)

# then calculate biodiversity metrics
kelp_rls_wide <- kelp_rls_wider %>%
  mutate(shannon = (vegan::diversity((kelp_rls_wider %>% select(-survey_id)), index = "shannon")),
         simpson = (vegan::diversity((kelp_rls_wider %>% select(-survey_id)), index = "simpson"))) %>%
  select(site_code, shannon, simpson)


# extract just one row per survey to join with the pee data and tide data
kcca_surveys <- kelp_rls %>%
              select(site_code, site_name, Date, date_time_survey) %>%
              unique() %>%
              mutate(survey_id = 101:117) %>%
  mutate(date_time_kelp = ymd_hms(paste(Date, Time_start)))

# Tide data ------

tide_kcca <- read_csv("Data/Team_kelp/kcca_tide_data.csv") %>%
  mutate(date_time = ymd_hms(paste(date, time)))

# build an empty dataframe
# Do this each time you run the for loop!!!
tide_exchange_kcca <- data.frame()

# then write the loop
for (x in 1:nrow(kcca_survey_info)) {
  
  survey_start <- ymd_hms(kcca_survey_info$date_time_kelp[x:x]) 
  # tie the times to the time the nh4 samples were collected, during kelp surveys
  survey_end <- survey_start + hours(1)
  
  output = tide_kcca %>%
    filter(between(date_time, survey_start, survey_end)) %>%
    mutate(rate = 100 * (tide_m - lag(tide_m))/lag(tide_m),
           max = rate[which.max(abs(rate))]) %>%
    slice(-1) %>%
    summarise(avg_exchange_rate = mean(rate),
              max_exchange_rate = mean(max)
    ) %>%
    mutate(survey_id = kcca_survey_info$survey_id[x:x])
  
  tide_exchange_kcca = rbind(tide_exchange_kcca, output)
}
# Now I have the average and max rate of change of tide height for each survey!!!


# Put them all together! ----
data <- pee %>% # put your nh4 df here
  # join site survey info for the survey ID
  left_join(kcca_survey_info %>%
              select(-c(Date, site_name)), by = "site_code") %>%
  # Join biomass, richness and abundance
  left_join(kelp_rls %>%
              group_by(site_code) %>%
              summarize(weight_sum = sum(weight_size_class_sum),
                        species_richness = n_distinct(species_name),
                        abundance = sum(survey_den)) %>%
              ungroup()) %>%
  # rls fish biomass 
  left_join(kelp_rls %>%
              filter(phylum == "Chordata") %>%
              group_by(survey_id) %>%
              summarize(fish_weight_sum = sum(weight_size_class_sum),
                        fish_abund = sum(survey_den))) %>%
  # rls CUKES!
  left_join(kelp_rls %>%
              filter(species_name == "Apostichopus californica") %>%
              group_by(survey_id) %>%
              summarize(cuke_weight_sum = sum(weight_size_class_sum),
                        cuke_abund = sum(survey_den))) %>%
  # diversity indexes
  left_join(kelp_rls_wide, by = "survey_id") %>%
  # join tide exchange data
  left_join(tide_exchange_kcca, by = "survey_id") %>%
  as.data.frame()%>%
  # scale and center variables here
  mutate(
    # the biomass + abundance variables
    weight_sum_scale = c(scale(weight_sum)),
    abundance_scale = c(scale(abundance)),
    # fish
    fish_weight_sum_scale = c(scale(fish_weight_sum)),
    fish_abund_scale = c(scale(fish_abund)),
    # cukes
    cuke_weight_sum_scale = c(scale(cuke_weight_sum)),
    cuke_abund_scale = c(scale(cuke_abund)),
    # biodiversity variables
    rich_scale = c(scale(species_richness)),
    shannon_scale = c(scale(shannon)),
    simpson_scale = c(scale(simpson)),
    # abiotic variables I should control for
    depth_avg_scale = c(scale(depth_avg)),
    tide_scale = c(scale(avg_exchange_rate)),
    # tide_cat = factor(as.factor(case_when(avg_exchange_rate < -0.1897325 ~ "Ebb",
    #                                       avg_exchange_rate < 0.1897325 ~ "Slack",
    #                                       avg_exchange_rate > 0.1897325 ~ "Flood")),
    #                   levels = c("Ebb", "Slack", "Flood")),
    
    # center instead of scale
    abundance_center = c(scale(abundance, scale = FALSE)),
    tide_center = c(scale(avg_exchange_rate, scale = FALSE)),
    shannon_center = c(scale(shannon, scale = FALSE)),
    depth_center = c(scale(depth_avg, scale = FALSE))
  )


# Stats -------

# Step 0: ask question
# Step 1: choose distribution
# Step 2: choose predictors

# overall abundance
mod_brain <- glmmTMB(nh4_avg ~ abundance_scale*tide_scale + shannon_scale + depth_avg_scale, 
                     family = Gamma(link = 'log'),
                     data = data)

summary(mod_brain)
plot(DHARMa::simulateResiduals(mod_brain))


# just look at fish weights
mod_fish <- glmmTMB(nh4_avg ~ fish_weight_sum_scale*tide_scale + shannon_scale + depth_avg_scale, 
                    family = Gamma(link = 'log'),
                    data = rls_final)
summary(mod_fish)
plot(DHARMa::simulateResiduals(mod_fish))

# just look at cukes abundance
mod_cuke <- glmmTMB(nh4_avg ~ cuke_abund_scale*tide_scale + shannon_scale + depth_avg_scale, 
                    family = Gamma(link = 'log'),
                    data = rls_final)
summary(mod_cuke)
plot(DHARMa::simulateResiduals(mod_cuke))

# AIC table
AIC_tab_rls <- AIC(mod_brain, mod_weight, mod_simp, mod_simp_weight) %>% # inset models here
  rownames_to_column() %>%
  mutate(best = min(AIC),
         delta = round((AIC - best), digits = 2),
         likelihood = exp( -0.5*delta),
         sum = sum(likelihood),
         AICw = round((likelihood/sum), digits = 2),
         AIC = round(AIC, digits = 2)) %>%
  select(rowname, df, AIC, delta, AICw)


# centered variables instead of scaled for estimates
mod_brain_c <- glmmTMB(nh4_avg ~ abundance_center*tide_center + shannon_center + depth_center, 
                       family = Gamma(link = 'log'),
                       data = rls_final)
summary(mod_brain_c)


# Step 4: Check for collinearity of predictors

# car can't handle random effects so make a simplified mod
car::vif(lm(nh4_avg ~ abundance_scale + tide_scale + shannon_scale + depth_avg_scale, data = rls_final))
# All good, shannon is a little high
# above 2 is BAD



