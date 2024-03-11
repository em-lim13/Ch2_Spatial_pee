# Script to look into the variation in pee inside vs outside kelp forests
# March 3, 2023
# Em Lim

# Load packages and functions -----
library(tidyverse)
library(MuMIn) # for dredging
library(visreg)
library(ggplot2)
library(TMB)
library(glmmTMB) # better for random effects
library(lubridate)
library(DHARMa)
library(broom.mixed)
library(dotwhisker)
library(ggeffects)


# set theme and load functions
theme_set(theme_bw())
source("Code/Functions.R")

# Kelp data ----

#Should I drop the no kelp control from the model and just say we found no diff in ammonium on the sand? Basically using the no kelp control as a methods control but not using it in the model
#Also double check that less dangerous bay was the no kelp control sand site???
#Bc now I can’t remember which site we didn’t write down the inside vs outside and probably got it messed up.. 

# load site names
names <- read_csv("Data/Team_kelp/Output_data/site_names.csv")

# transect level kelp biomass + density data from Claire
kelp <- read_csv("Data/Team_kelp/Output_data/transect_biomass.csv") %>%
  as.data.frame() %>% 
  mutate(sample = case_when(Transect == 0 ~ 1, 
                            Transect == 5 ~ 2,
                            Transect == 10 ~ 3),
         kelp_sp = as.factor(case_when(
           Kelp == 0 ~ "none", # no kelp = none
           SiteName == "Wizard Islet North" ~ "none", # no kelp site that had 2 nereo
           Macro_5m2 == 0 ~ "nereo", # no macro = nereo
           TRUE ~ as.character("macro"))) # everything else = macro
         ) %>% 
  left_join(names, by = "SiteName") %>% 
  replace(is.na(.), 0) %>%
  # Add the averaged site level variables from Claire!
  left_join(read_csv("Data/Team_kelp/Output_data/kelpmetrics_2022.csv"), by = "SiteName") %>%
  group_by(SiteName) %>%
  mutate(BiomassM = mean(Biomassm2kg),
         Area_m2 = ifelse(kelp_sp == "none", 0, Area_m2)) %>%
  ungroup() %>%
  rename(kelp_den = Kelp,
         site_name = SiteName) %>%
filter(Transect != "15")
  
    #select(site_code, SiteName, HeightT, BiomassTkg, Biomassm2kg, sample) 

# Less dangerous bay was the "real" no kelp control!!!!
# Sand town was a super weird outlier, I think that was the site we messed up the samples from, that's AOKAY!!! phew
# both dropped from Claire's stuff because she's not looking at no kelp sites


# Transect level variables:
  # BiomassTkg is kelp density x average biomass for each transect
  # Biomassm2kg is just that number / 5 bc the transect was 5 m2

# Site level variables
# BiomassM is the average biomass/m2 at each site across all 4 transects
# DensityM is the average density at each site across all 4 transects


# NH4+ data -----
# kelp pee inside vs outside data from Kelp_pee_nh4_calc.R
pee <- read_csv("Data/Team_kelp/Output_data/kelp_pee.csv")


# RLS data ----
kelp_rls1 <- read_csv("Data/Team_Kelp/RLS_KCCA_2022.csv") %>%
  as.data.frame() %>%
  filter(Method != 0) %>% # get rid of all method 0's
  slice(2:n()) %>% # cuts the first blank row
  rename(
    site_code = `Site No.`,
    site_name = `Site Name`, 
    common_name = `Common name`,
    `0` = Inverts,
    species_name = Species
  )  %>% # Rename columns with spaces
  mutate(species_name = str_to_sentence(species_name),
         common_name = str_to_sentence(common_name),
         Date = dmy(Date),
         date_time_survey = ymd_hms(paste(Date, Time))) %>%
  clean_sp_names() %>%
  left_join(read_csv("Output/Output_data/rls_phylo.csv"), by = "species_name") %>%
  clean_phylo_names() %>% # function to fix naming errors
  filter(species_name != "Tonicella spp.") %>%
  filter(species_name != "Cryptochiton stelleri")

# Pivot longer for biomass
kelp_rls <- kelp_rls1 %>%
  pivot_longer( cols = `0`:`400`, names_to = "size_class", values_to = "total") %>%
  mutate(size_class = as.numeric(size_class)) %>%
  drop_na(total) %>%
  filter(total > 0) %>%
  select(-Total) %>%
  length_to_weight() %>% # fish length to weight function
  home_range() # calculate each fish's home range function
  

# save csv for mapping 
#  kelp_rls_csv <- kelp_rls %>%
#    transmute(site_code = site_code,
#              latitude = Latitude,
#              longitude = Longitude) %>%
#    unique()
#  
#  write_csv(kelp_rls_csv, "Output/Output_data/kelp_rls.csv")


# extract just one row per survey to join with the pee data and tide data
# Only keep the surveys I have pee samples for!
kcca_survey_info <- pee %>%
  transmute(Date = date,
            site_code = site_code) %>%
  unique() %>%
  left_join(kelp_rls %>%
              select(site_code, site_name, Date, date_time_survey) %>%
              unique() %>%
              mutate(survey_id = 101:127), by = c("Date", "site_code")) %>%
  left_join((kelp %>%
              select(c(site_code, Time_start)) %>%
              unique()), by = "site_code") %>%
  mutate(date_time_kelp = ymd_hms(paste(Date, Time_start)))

# pivot back to wide for biodiversity
kelp_rls_wider <- kelp_rls %>% 
  left_join(kcca_survey_info, by = c("site_code", "date_time_survey")) %>%
  dplyr::select(survey_id, species_name, total) %>% 
  group_by(survey_id, species_name) %>%
  summarise(total = sum(total)) %>%
  ungroup() %>%
  spread(key = species_name, value = total) %>%
  replace(is.na(.), 0)

# then calculate biodiversity metrics
kelp_rls_wide <- kelp_rls_wider %>%
  mutate(shannon = (vegan::diversity((kelp_rls_wider %>% select(-survey_id)), index = "shannon")),
         simpson = (vegan::diversity((kelp_rls_wider %>% select(-survey_id)), index = "simpson"))) %>%
  select(survey_id, shannon, simpson)

# Tide data ------
# July 7, 4 days
# July 24, 4 days
# Aug 3, 4 days
# Aug 7, 1 day
# Aug 18, 4 days
# Aug 22, 1 day
# Sept 1 + 5, each 1 day

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
# no diff in models between avg and max tide exchange, just use avg


# Put them all together! ----
data_all <- pee %>%
  # so some site level averaging
  group_by(site_code) %>%
  mutate(nh4_avg = mean(c(nh4_outside, nh4_inside)),
         nh4_in_avg = mean(nh4_inside),
         nh4_out_avg = mean(nh4_outside),
         in_out_avg = mean(in_minus_out),
         percent_diff = 100*(nh4_inside-nh4_outside)/nh4_outside,
         depth_avg = mean(depth_m)) %>%
  ungroup() %>%
  # join kelp data with both transect + site avgs
  left_join(kelp, by = c("site_code", "sample")) %>%
  # join site survey info for the survey ID
  left_join(kcca_survey_info %>%
              select(-c(Date, site_name)), by = "site_code"
            ) %>%
  # Join biomass, richness and abundance
  left_join(kelp_rls %>%
              group_by(site_code) %>%
              summarize(weight_sum = sum(weight_size_class_sum),
                        all_weight_weighted = sum(weight_weighted),
                        species_richness = n_distinct(species_name),
                        abundance = sum(total))) %>%
  # diversity indexes
  left_join(kelp_rls_wide, by = "survey_id") %>%
  # join tide exchange data
  left_join(tide_exchange_kcca, by = "survey_id")


# Remove Second Beach South for in vs out
data <- data_all %>%
  filter(site_name != "Second Beach South") %>%
  # scale factors using function
  scale_vars() 

# reduced data to export a file for mapping
data_map <- data %>%
  select(site, site_code, survey_id, in_out_avg, nh4_in_avg, nh4_out_avg, nh4_avg, depth_avg, avg_exchange_rate, kelp_sp, BiomassM, kelp_bio_scale, forest_bio_scale, weight_sum, weight_sum_scale, all_weighted_scale, abundance_scale, rich_scale, shannon_scale, simpson_scale, tide_scale, depth_scale, tide_cat) %>%
  unique()
#  filter(site != "Wizard_I_North" | kelp_sp != "none") # just getting rid of the duplicate row, bc there wasn't kelp on one transect unique misses this one

# write_csv(data_map, "Output/Output_data/kelp_final.csv")
  


# Stats for in minus out -----

# Step 0: Ask my question
  # Does kelp density affect the retention of ammonium in kelp forests?

# Response variable = in_minus_out
  # how much ammonium is inside the forest vs outside
  # I'm using the transect level data instead of the averaged site data, the residuals look better and it makes more sense

# Step 1: Choose a distribution
ggplot(data, aes(x = in_minus_out)) +
  geom_histogram(bins = 30) 
# Hey that's roughly Gaussian!
# Positive and negative values roughly centered around 0

# Step 2: Choose your predictors
# And Step 3: Model residuals!

# Biological
  # Kelp forest
    # forest_bio_scale (site level = area*biomass*density)
    # kelp_sp (site level)
    # kelp_bio_scale (site level, density*biomass)
    # log_kelp_scale log(bio_mean), also scaled
    # log_den (site level)
    # bio_tran_scale (mini transect level)
    # den_tran_scale (mini transect level)
# I'm using kelp_bio_scale, best for AIC and better than log(kelp) too
# Tried including one transect level variable and it wasn't as good

  # RLS community (site level)
    # Biomass:
      # weight_sum_scale 
      # all_weighted_scale 
      # abundance 
    # Biodiversity
      # rich_scale
      # shannon_scale
      # simpson_scale
# Using weight_sum_scale and shannon, best AIC

# Abiotic to control for
  # depth_avg (depth_scale)
  # avg_exchange_rate (tide_scale)
# Depth doesn't matter but tide does!

# Random variables
  # side_code: if I use mini transect as the level of study than I need a random effect of site

# Interactions:
  # interaction between animal biomass:tide_exchange
  # also interaction between kelp:tide_exchange
  # also interaction between animal biomass:kelp
  # maybe a triple interaction between kelp:animal_biomass:tide??
# Only kelp:animal interaction! 

# Build full model with gaussian distribution
# This is the transect level model
# One row per transect, 3 per site
mod_tran <- glmmTMB(in_minus_out ~ -1 + kelp_sp + 
                      kelp_bio_scale*tide_scale*weight_sum_scale + 
                      bio_tran_scale + shannon_scale + depth_scale + 
                      (1|site_code),
                    family = 'gaussian',
                    data = data,
                    na.action = na.fail)

plot(simulateResiduals(mod_tran)) # looks fine
summary(mod_tran)

# Dredge to find best combo of variables
# dredge <- as.data.frame(dredge(mod_tran)) %>% filter(delta < 3)

# Without second beach:
# best = kelp_bio_scale*tide_scale + weight_sum_scale + shannon_scale + kelp_sp 
# second best = same + weight_sum_scale:tide_scale

mod_best <- glmmTMB(in_minus_out ~ kelp_sp + kelp_bio_scale*tide_scale + weight_sum_scale + shannon_scale + (1|site_code), 
                    family = 'gaussian',
                    data = data)

plot(simulateResiduals(mod_best)) # looks fine
summary(mod_best)

# use my stupid brain to think of a good model
mod_in_out <- glmmTMB(in_minus_out ~ kelp_sp +
                        kelp_bio_scale*tide_scale*weight_sum_scale + 
                        shannon_scale + depth_scale - 
                        kelp_bio_scale:tide_scale:weight_sum_scale +
                        (1|site_code),
                      family = 'gaussian',
                      data = data)

plot(simulateResiduals(mod_in_out)) # looks fine
summary(mod_in_out)

# biomass and abundance, shannon vs simpson checks
mod_abund <- glmmTMB(in_minus_out ~ kelp_sp +
                        kelp_bio_scale*tide_scale*abundance_scale + 
                        shannon_scale + depth_scale - 
                        kelp_bio_scale:tide_scale:abundance_scale +
                        (1|site_code),
                      family = 'gaussian',
                      data = data)

mod_simp <- glmmTMB(in_minus_out ~ kelp_sp +
                        kelp_bio_scale*tide_scale*weight_sum_scale + 
                        simpson_scale + depth_scale - 
                        kelp_bio_scale:tide_scale:weight_sum_scale +
                        (1|site_code),
                      family = 'gaussian',
                      data = data)

mod_abund_simp <- glmmTMB(in_minus_out ~ kelp_sp +
                       kelp_bio_scale*tide_scale*abundance_scale + 
                         simpson_scale + depth_scale - 
                       kelp_bio_scale:tide_scale:abundance_scale +
                       (1|site_code),
                     family = 'gaussian',
                     data = data)

mod_tran <- glmmTMB(in_minus_out ~ kelp_sp +
                      bio_tran_scale*tide_scale*weight_sum_scale + 
                      shannon_scale + depth_scale - 
                      bio_tran_scale:tide_scale:weight_sum_scale +
                      (1|site_code),
                      family = 'gaussian',
                      data = data)

mod_forest <- glmmTMB(in_minus_out ~ kelp_sp +
                        forest_bio_scale*tide_scale*weight_sum_scale + 
                        shannon_scale + depth_scale - 
                        forest_bio_scale:tide_scale:weight_sum_scale +
                        (1|site_code),
                    family = 'gaussian',
                    data = data)


AIC_tab_kelp <- AIC(mod_in_out, mod_abund, mod_simp, mod_abund_simp) %>%
  rownames_to_column() %>%
  mutate(best = min(AIC),
         delta = AIC - best,
         likelihood = exp( -0.5*delta),
         sum = sum(likelihood),
         AICw = likelihood/sum) %>%
  select(rowname, df, AIC, delta, AICw)

# run the model without an intercept
mod_in_out2 <- glmmTMB(in_minus_out ~ -1 + kelp_sp +
                        kelp_bio_scale*tide_scale*weight_sum_scale + 
                        shannon_scale + depth_scale - 
                        kelp_bio_scale:tide_scale:weight_sum_scale +
                        (1|site_code),
                      family = 'gaussian',
                      data = data)
summary(mod_in_out2)

# center instead of scale for estimates in normal units
mod_in_out_c <- glmmTMB(in_minus_out ~ - 1 + kelp_sp +
                          kelp_bio_center*tide_center*weight_sum_center + 
                          shannon_center + depth_center - 
                          kelp_bio_center:tide_center:weight_sum_center +
                          (1|site_code),
                        family = 'gaussian',
                        data = data)

summary(mod_in_out_c)

# Step 4: Check for collinearity of predictors

# car can't handle random effects so make a simplified mod
car::vif(lm(in_minus_out ~ kelp_sp + kelp_bio_scale + tide_scale + weight_sum_scale + shannon_scale +  depth_scale, data = data))
# Yep this looks fine!


# Notes from stats beers!
#  Double check that the negative in – out sites aren’t the places where the transect was on the other side of the forest.
# Maybe try to PCA????? see if there's clustering?
  

# Graphing ------

# Palettes
pal12 <- viridis::viridis(11)
pal_k <- viridis::viridis(10)
pal2 <- c(pal_k[8], pal_k[5])

# Save model coefficients 
# use the no intercept model for plotting
df <- confint(mod_in_out2, level = 0.95, method = c("wald"), component = c("all", "cond", "zi", "other"), estimate = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(variable = rowname,
         lower_CI = `2.5 %`,
         upper_CI = `97.5 %`,
         estimate = Estimate) %>%
  head(- 1)  %>%
  mutate(variable = factor(as.factor(variable), 
                           levels = c("kelp_spmacro", "kelp_spnereo", "kelp_spnone", "shannon_scale", "depth_scale", "kelp_bio_scale", "weight_sum_scale", "tide_scale", "kelp_bio_scale:tide_scale", "kelp_bio_scale:weight_sum_scale", "tide_scale:weight_sum_scale"),
                           labels = c("Macro", "Nereo", "No kelp", "Biodiversity", "Depth",  "Kelp biomass", "Animal biomass", "Tide", "Kelp:tide", "Kelp:animals", "Tide:animals")),
         se = (upper_CI - estimate)/1.96
  )

# Use function to plot coefficients
kelp_coeff_plot <- coeff_plot(coeff_df = df,
                              pal = pal12) +
  place_label("(a)")

# ggsave("Output/Figures/kelp_in_out_mod_coeff.png", device = "png", height = 8, width = 12, dpi = 400)


# Model predictions against raw data plot
# make new mod with untransformed vars

# create range vector
tide_means_kelp <- data %>%
  group_by(tide_cat) %>%
  summarise(tide = mean(tide_scale))

v2 <- c(as.numeric(tide_means_kelp[2,2]), as.numeric(tide_means_kelp[1,2]))

# now make predictions
predict_kelp <- ggpredict(mod_in_out, terms = c("kelp_bio_scale", "tide_scale [v2]")) %>% 
  mutate(kelp_bio_scale = x,
         tide_cat = factor(as.factor(ifelse(group == as.character(tide_means_kelp[1,2]), "Slack", "Flood")),
         levels = c("Ebb", "Slack", "Flood"))
         ) %>%
  filter(tide_cat != "Flood" | kelp_bio_scale < 0.21) %>%
  filter(tide_cat != "Flood" | kelp_bio_scale > -0.6)


# now plot these predictions
kelp_pred_plot <- plot_kelp_pred(raw_data = data, 
                                 predict_data = predict_kelp, 
                                 pal = pal2,
                                 theme_white = TRUE) +
  place_label("(b)")

# try again with new function
plot_pred(raw_data = (data %>%
                        mutate(tide = ifelse(avg_exchange_rate < 0, "Slack", "Flood"))),
          predict_data = predict_kelp, 
          plot_type = "kelp",
          x_var = kelp_bio_scale, y_var = in_minus_out, 
          lty_var = tide_cat,
          size_var = weight_sum, 
          pal = pal2) +
  place_label("(b)")

# ggsave("Output/Figures/kelp_in_out_mod_predict.png", device = "png", height = 9, width = 12, dpi = 400)


# Figure 3 for pub with white -----
kelp_coeff_plot + kelp_pred_plot

#ggsave("Output/Pub_figs/Fig4.png", device = "png", height = 9, width = 16, dpi = 400)



# Stats: site to site variation ----

# This is not complete!!!!

# Step 0: Ask my question
  # Is there a link between biodiversity or biomass and ammonium concentration? 

# Step 1: Choose a distribution

# response variables:
  # nh4_avg (average of in and out)
  # nh4_in_avg (inside kelp avg)
  # nh4_out_avg (outside kelp avg)
    # most comparable to the RLS blitz work?

ggplot(data_s, aes(x = nh4_out_avg)) +
  geom_histogram(bins = 30) 
# Either way we're using Gamma!

# Yes I confirmed that gamma is better than gaussian

# Step 2: Choose your predictors

# Biological predictors:
  # 1) Animal biomass
    # a) weight_sum = total wet weight of all the animals observed on RLS transects
    # b) all_weight_weighted = weighted fish weight + invert weight
    # c) abundance
      # I'm going with weight_sum
  # 2) Kelp forest
    # a) forest_bio_scale (site level = area*biomass*density)
    # b) kelp_bio_scale (site level, density*biomass)
    # c) kelp_sp (site level)

  # 3) Biodiversity
   # a) species_richness = total # species on each transect
   # b) shannon = shannon diversity of each transect
   # c) simpson = simpson diversity of each transect
     # I'm going with shannon

# Continuous abiotic predictors
  # 1) depth_avg = average nh4 sample depths
  # 2) avg_exchange_rate = average rate of change of the tide height over the 1 hour survey

# Interactions
  # I think whatever biomass and biodiversity variable I choose should have an interaction with exchange rate
  # maybe triple kelp:biomass:tide interaction

# And Step 3: Model residuals 

# Build full model with gamma distribution
mod_kelp_site <- glmmTMB(nh4_out_avg ~ weight_sum_scale*tide_scale*kelp_bio_scale + shannon_scale + kelp_sp + depth_scale - weight_sum_scale:tide_scale:kelp_bio_scale, 
                        family = Gamma(link = 'log'),
                        data = data_s,
                        na.action = na.fail)
summary(mod_kelp_site)
plot(simulateResiduals(mod_kelp_site))

# try making a model with MORE interactions!
mod_kelp_site_ints <- glmmTMB(nh4_out_avg ~ weight_sum_scale*tide_scale*kelp_bio_scale*shannon_scale + kelp_sp + depth_scale - 
                           weight_sum_scale:tide_scale:kelp_bio_scale:shannon_scale -
                           weight_sum_scale:tide_scale:kelp_bio_scale - 
                           weight_sum_scale:tide_scale:shannon_scale -
                           weight_sum_scale:kelp_bio_scale:shannon_scale -
                           tide_scale:kelp_bio_scale:shannon_scale, 
                         family = Gamma(link = 'log'),
                         data = data_s)
summary(mod_kelp_site_ints) 
plot(simulateResiduals(mod_kelp_site_ints)) # when I put Second beach south back in, the residuals look fine. Why is that. why. 

AIC(mod_kelp_site, mod_kelp_site_ints) # mod_kelp_site is better

# Finalize which predictors are best
# dredge <- as.data.frame(dredge(mod_gam_kelp)) %>% filter(delta < 3)
# best model: just kelp_sp + tide_scale

mod_top <- glmmTMB(nh4_out_avg ~ kelp_sp + tide_scale, 
                   family = Gamma(link = 'log'),
                   data = data_s)
summary(mod_top)
plot(simulateResiduals(mod_top)) # looks good.....

# But I'm sticking with my first model!


# Step 4: Check for collinearity of predictors

# car can't handle random effects so make a simplified mod
car::vif(lm(nh4_out_avg ~ weight_sum_scale + tide_scale + kelp_bio_scale + shannon_scale + kelp_sp + depth_scale, data = data_s))
# all good


# Plot kelp site mod ------
# generate df for plotting
 
df2 <- confint(mod_kelp_site, level = 0.95, method = c("wald"), component = c("all", "cond", "zi", "other"), estimate = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(variable = rowname,
         lower_CI = `2.5 %`,
         upper_CI = `97.5 %`,
         estimate = Estimate) %>%
  mutate(estimate = ifelse(variable == "(Intercept)" |
                             variable == "kelp_spmixed" |
                             variable == "kelp_spnereo" |
                             variable == "kelp_spnone",
                           exp(estimate), estimate),
         lower_CI = ifelse(variable == "(Intercept)" |
                             variable == "kelp_spmixed" |
                             variable == "kelp_spnereo" |
                             variable == "kelp_spnone",
                           exp(lower_CI), lower_CI),
         upper_CI = ifelse(variable == "(Intercept)" |
                             variable == "kelp_spmixed" |
                             variable == "kelp_spnereo" |
                             variable == "kelp_spnone",
                           exp(upper_CI), upper_CI),
         variable = factor(as.factor(variable), 
                           levels = c("(Intercept)", "kelp_spmixed", "kelp_spnereo", "kelp_spnone", "shannon_scale", "depth_scale", "kelp_bio_scale", "weight_sum_scale", "tide_scale", "tide_scale:kelp_bio_scale", "weight_sum_scale:kelp_bio_scale", "weight_sum_scale:tide_scale"),
                           labels = c("Intercept", "Mixed kelp", "Nereo", "No kelp", "Biodiversity", "Depth",  "Kelp biomass", "Animal biomass", "Tide", "Kelp:tide", "Kelp:animals", "Tide:animals")))

# Coefficient plot
pal12 <- viridis::viridis(12)

# plot un-logged cat variables
ggplot(df2, aes(x = estimate, y = (variable), xmin = lower_CI, xmax = upper_CI, colour = variable)) +
  geom_point(size = 10) +
  geom_errorbar(width = 0, linewidth = 3) +
  geom_vline(xintercept=0, color="white", linetype="dashed") +
  labs(x = "Coefficient", y = " ") +
  scale_y_discrete(limits = rev(levels(df2$variable))) +
  theme_black() +
  theme(legend.position = "none") + 
  scale_colour_manual(values = pal12)

# ggsave("Output/Figures/kelp_site_mod_coeff.png", device = "png", height = 9, width = 12, dpi = 400)

# plot model predictions
ggplot(data = data_s, 
       aes(x = tide_scale, y = nh4_out_avg,
           colour = kelp_sp,
           fill = kelp_sp)) + 
  geom_point(alpha = 0.8) +
  geom_smooth(method = lm) +
  labs(y = expression(paste(Delta, " Ammonium ", (mu*M))), 
       x = "Tide exchange (m/s)",
       colour = "Kelp sp", fill = "Kelp sp") +
  theme_black() + 
  scale_colour_manual(values = pal5) +
  scale_fill_manual(values = pal5)

# fucking around -----

kcca_table <- pee %>%
  select(site_code, date) %>%
  unique() %>%
  filter(site_code != "KCCA13")



# Graveyard ------
data2 <- data %>%
  select(site, kelp_sp, Composition)
