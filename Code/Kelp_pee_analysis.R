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

# set theme and load functions
theme_set(theme_bw())
source("Code/Functions.R")

# Kelp data ----

# load site names
names <- read_csv("Data/Team_kelp/Output_data/site_names.csv")

# transect level kelp biomass + density data from Claire
kelp <- read_csv("Data/Team_kelp/Output_data/transect_biomass.csv") %>%
  as.data.frame() %>%
  mutate(sample = ifelse(Transect == 0, 1, 
                  ifelse(Transect == 5, 2, 3)),
         kelp_sp = as.factor(ifelse(Kelp == 0, "none",
                          ifelse(Macro_5m2 == 0, "nereo", 
                                 ifelse(Nereo_5m2 == "0", "macro", "mixed")))),
         # but the "mixed" kelps are actually mostly macro
         kelp_sp = case_when(kelp_sp == "mixed" ~ "macro", 
                             TRUE ~ as.character(kelp_sp))) %>%
  left_join(names, by = "SiteName") %>% 
  replace(is.na(.), 0) %>%
  # Add the averaged site   level variables from Claire!
  left_join(read_csv("Data/Team_kelp/Output_data/kelpmetrics_2022.csv"), by = "SiteName") %>%
  group_by(SiteName) %>%
  mutate(BiomassM = mean(Biomassm2kg),
         Area_m2 = ifelse(Composition == "None", 0, 
                          ifelse(site_code == "KCCA19", 0.1, Area_m2))) %>%
  ungroup() %>%
  rename(kelp_den = Kelp,
         site_name = SiteName) %>% 
  filter(Transect != "15") 
    #select(site_code, SiteName, HeightT, BiomassTkg, Biomassm2kg, sample) 

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
    common_name = `Common name`
  )  %>% # Rename columns with spaces
  mutate(Species = str_to_sentence(Species),
         common_name = str_to_sentence(common_name),
         Date = dmy(Date),
         date_time_survey = ymd_hms(paste(Date, Time))) %>% 
  filter(Species != "Debris - Metal") %>%
  filter(Species != "Debris - Other") %>%
  filter(Species != "Debris - Wood") %>%
  filter(Species != "Debris - Glass") %>%
  filter(Species != "Debris - Fishing Gear") %>%
  filter(Species != "Debris - Zero")  # Cut the non-animal species

# Pivot longer for biomass
kelp_rls <- kelp_rls1 %>%
  rename(`0` = Inverts,
         species_name = Species) %>%
  pivot_longer( cols = `0`:`400`, names_to = "size_class", values_to = "total") %>%
  mutate(size_class = as.numeric(size_class)) %>%
  drop_na(total) %>%
  filter(total > 0) %>%
  select(-Total) %>%
  invert_length_to_weight() %>% # invert length to weight
  length_to_weight() %>% # fish length to weight function
  home_range() # calculate each fish's home range

# save csv for mapping 
  #kelp_rls_csv <- kelp_rls %>%
  #  transmute(site_code = site_code,
  #            latitude = Latitude,
  #            longitude = Longitude) %>%
  #  unique()
  #
  #write_csv(kelp_rls_csv, "Output/Output_data/kelp_rls.csv")


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

# reduced data for the kelp site model that includes second beach
data_s <- data_all %>%
  # scale factors using function
  scale_vars() %>%
  select(site, site_code, survey_id, in_out_avg, nh4_in_avg, nh4_out_avg, nh4_avg, depth_avg, avg_exchange_rate, kelp_sp, BiomassM, kelp_bio_scale, forest_bio_scale, weight_sum, weight_sum_scale, all_weighted_scale, abundance_scale, rich_scale, shannon_scale, simpson_scale, tide_scale, depth_scale, tide_cat) %>%
  unique() %>%
  filter(site != "Wizard_I_North" | kelp_sp != "none") # just getting rid of the duplicate row, bc there wasn't kelp on one transect unique misses this one

#write_csv(data_s, "Output/Output_data/kelp_final.csv")
  


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
mod_tran <- glmmTMB(in_minus_out ~ kelp_sp + kelp_bio_scale*tide_scale*weight_sum_scale + bio_tran_scale + shannon_scale + depth_scale + (1|site_code), 
                         family = 'gaussian',
                         data = data,
                    na.action = na.fail)

# Dredge to find best combo of variables
# dredge <- as.data.frame(dredge(mod_tran)) %>% filter(delta < 3)

# Without second beach:
# best = kelp_bio_scale*tide_scale + weight_sum_scale + shannon_scale + kelp_sp 
# second best = same + weight_sum_scale:tide_scale

# with Second beach:
# best = kelp_bio_scale*tide_scale*weight_sum_scale + kelp_sp + shannon_scale 
# second best is the same + bio_tran_scale

mod_best <- glmmTMB(in_minus_out ~ +kelp_sp + kelp_bio_scale*tide_scale + weight_sum_scale + shannon_scale + (1|site_code), 
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

# run the model without an intercept
mod_in_out2 <- glmmTMB(in_minus_out ~ -1 + kelp_sp +
                        kelp_bio_scale*tide_scale*weight_sum_scale + 
                        shannon_scale + depth_scale - 
                        kelp_bio_scale:tide_scale:weight_sum_scale +
                        (1|site_code),
                      family = 'gaussian',
                      data = data)
summary(mod_in_out2)



# Step 4: Check for collinearity of predictors

# car can't handle random effects so make a simplified mod
car::vif(lm(in_minus_out ~ kelp_bio_scale + tide_scale + weight_sum_scale + shannon_scale + kelp_sp + depth_scale, data = data))
# Yep this looks fine!


# Notes from stats beers!
#  Double check that the negative in – out sites aren’t the places where the transect was on the other side of the forest.
# Maybe try to PCA????? see if there's clustering?
  

# Try to plot in minus out ------

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

# Coefficient plot
pal12 <- viridis::viridis(12)
pal <- viridis::viridis(10)
pal2 <- c(pal[8], pal[5])

# marginal means for cats
ggplot(df, aes(x = estimate, y = variable, 
               xmin = lower_CI, xmax = upper_CI, 
               colour = variable)) +
  geom_point(size = 10) +
  geom_errorbar(width = 0, linewidth = 3) +
  geom_vline(xintercept=0, color="white", linetype="dashed") +
  labs(x = "Coefficient", y = " ") +
  scale_y_discrete(limits = rev(levels(df$variable))) +
  theme_black() +
  theme(legend.position = "none") + 
  scale_colour_manual(values = pal12)

# ggsave("Output/Figures/kelp_in_out_mod_coeff.png", device = "png", height = 8, width = 12, dpi = 400)


# remake the shitty asymptote fig with the log(kelp) model and actual predictions, with geom_point(aes(colour = tide, size = animal weight))! See what that looks like???

# make new mod with untransformed vars

# create range vector
tide_means_kelp <- data %>%
  group_by(tide_cat) %>%
  summarise(tide = mean(tide_scale))

v2 <- c(1.5020449, -0.5006816)

# now make predictions
predict_kelp <- ggpredict(mod_in_out, terms = c("kelp_bio_scale", "tide_scale [v2]")) %>% 
  mutate(kelp_bio_scale = x,
         tide_cat = factor(as.factor(ifelse(group == "-0.5006816", "Slack", "Flood")),
         levels = c("Ebb", "Slack", "Flood"))
         ) %>%
  filter(tide_cat != "Flood" | kelp_bio_scale < 0.21) %>%
  filter(tide_cat != "Flood" | kelp_bio_scale > -0.6)
  
# now plot these predictions
plot_kelp_pred(data, predict_kelp) +
  geom_hline(yintercept= 0, linetype = "dashed", color = "white", linewidth = 0.5) +
  guides(size = guide_legend(override.aes = 
                               list(colour = "white")),
         lty = guide_legend(override.aes = list(linewidth = 0.5))) 



# get original axis
#ggplot(data, aes(BiomassM, kelp_bio_scale)) +
#  geom_smooth() +
#  geom_vline(xintercept= 1.8,  color = "red", linewidth = 0.5) +
#  scale_y_continuous(n.breaks = 50)


# ggsave("Output/Figures/kelp_in_out_mod_predict.png", device = "png", height = 9, width = 12, dpi = 400)


# Figure 3 for pub with white -----
# plot coeffs on white
kelp_coeff_plot <- ggplot(df, aes(x = estimate, y = variable, 
               xmin = lower_CI, xmax = upper_CI, 
               colour = variable)) +
  geom_point(size = 10) +
  geom_errorbar(width = 0, linewidth = 3) +
  geom_vline(xintercept=0, color="black", linetype="dashed") +
  labs(x = "Coefficient", y = " ") +
  scale_y_discrete(limits = rev(levels(df$variable))) +
  theme_white() +
  theme(legend.position = "none") + 
  scale_colour_manual(values = pal12)

# plot predictions on white
kelp_pred_plot <- plot_kelp_pred(data, predict_kelp) +
  geom_hline(yintercept= 0, linetype = "dashed", color = "white", linewidth = 0.5) +
  guides(size = guide_legend(override.aes = 
                               list(colour = "black")),
         lty = guide_legend(override.aes = list(linewidth = 0.5))) +
  theme_white()

# Put them together for pub
kelp_coeff_plot + kelp_pred_plot

# ggsave("Output/Pub_figs/Fig3.png", device = "png", height = 9, width = 16, dpi = 400)


# # Ammonium at the site level?
data_s %>%
  mutate(site = fct_reorder(site, in_out_avg, .fun='median')) %>%
  ggplot() +
  geom_point(aes(x = reorder(site, in_out_avg), nh4_out_avg, 
                 colour = "blue")) +
  geom_point(aes(x = reorder(site, in_out_avg), nh4_in_avg, 
                 colour = "green")) +
  labs(y = "NH4 concentration", x = "Site") + 
  theme(axis.text.x=element_text(angle = 65, hjust = 1))+
  scale_color_identity(name = "Sample",
                       breaks = c("blue", "green"),
                       labels = c("Outside kelp", "Inside kelp"),
                       guide = "legend")
# Ok this is sort of neat, you can see the base ammonium levels and then how different they are, by site. might be neat to arrange these with kelp density instead of site on the x
ggplot(data) +
  geom_point(aes(x = kelp_bio_scale, nh4_out_avg, 
                 colour = "blue")) +
  geom_point(aes(x = kelp_bio_scale, nh4_in_avg, 
                 colour = "green")) +
  labs(y = "NH4 concentration", x = "Kelp biomass") + 
  theme(axis.text.x=element_text(angle = 65, hjust = 1))+
  scale_color_identity(name = "Sample",
                       breaks = c("blue", "green"),
                       labels = c("Outside kelp", "Inside kelp"),
                       guide = "legend")


# Stats: site to site variation ----

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

 ggsave("Output/Figures/kelp_site_mod_coeff.png", device = "png", height = 9, width = 12, dpi = 400)

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




# Graveyard ------

# are the continuous predictors the same for all kelp species :|
# what if I make nereo first instead of macro?
n_data <- data %>%
  mutate(kelp_sp = factor(kelp_sp, levels = c("nereo", "macro", "mixed", "none")))

mod_in_out_nereo <- glmmTMB(in_minus_out ~ -1 + kelp_sp + kelp_bio_scale*tide_scale*weight_sum_scale + shannon_scale + depth_scale - kelp_bio_scale:tide_scale:weight_sum_scale +
                              (1|site_code), 
                            family = 'gaussian',
                            data = n_data)

summary(mod_in_out_nereo)

# plot model
tt <- (list(no_mix = mod_in_out2, mix = mod_in_out_nereo)
       %>% purrr::map_dfr(tidy, effects = "fixed", conf.int = TRUE,
                          .id = "model")
       %>% select(model, component, term, estimate, conf.low, conf.high)
       ## create new 'term' that combines component and term
       %>% mutate(term_orig = term,
                  term = forcats::fct_inorder(paste(term, component, sep = "_")))
)

dwplot(tt)

# it doesn't seem to matter how I order the kelp_sp factor, all of the estimates are EXACTLY the same



# old kelp pee coeff plot using mod with intercept
int <- confint(mod_in_out, estimate = TRUE)[1,3]

df2 <- confint(mod_in_out, level = 0.95, method = c("wald"), component = c("all", "cond", "zi", "other"), estimate = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(variable = rowname,
         lower_CI = `2.5 %`,
         upper_CI = `97.5 %`,
         estimate = Estimate) %>%
  head(- 1)  %>%
  mutate(variable = factor(as.factor(variable), 
                           levels = c("(Intercept)", "kelp_spnereo", "kelp_spnone", "shannon_scale", "depth_scale", "kelp_bio_scale", "weight_sum_scale", "tide_scale", "kelp_bio_scale:tide_scale", "kelp_bio_scale:weight_sum_scale", "tide_scale:weight_sum_scale"),
                           labels = c("Macro", "Nereo", "No kelp", "Biodiversity", "Depth",  "Kelp biomass", "Animal biomass", "Tide", "Kelp:tide", "Kelp:animals", "Tide:animals")),
         adj_estimate = case_when(variable == "Nereo" ~ estimate + int,
                                  variable == "No kelp" ~ estimate + int,
                                  TRUE ~ estimate),
         se = (upper_CI - estimate)/1.96,
         ci.lb_adjust = adj_estimate - (1.96*se),
         ci.up_adjust = adj_estimate + (1.96*se)
  )


# just urchins
urchins_kelp <- kelp_rls %>%
  filter(species_name == "Mesocentrotus franciscanus") %>%
  group_by(site_code) %>%
  summarize(urchins = sum(total)) %>%
  left_join(kelp %>% select(site_code, Composition) %>% unique())


# try making a more simple model???
mod_gam_kelp3 <- glmmTMB(nh4_out_avg ~ weight_sum_scale, 
                         family = Gamma(link = 'log'),
                         data = data_s)
plot(simulateResiduals(mod_gam_kelp3)) # two zig zags but no red

mod_gam_kelp4 <- glmmTMB(nh4_out_avg ~ tide_scale, 
                         family = Gamma(link = 'log'),
                         data = data_s)
plot(simulateResiduals(mod_gam_kelp4)) # super fucked

mod_gam_kelp5 <- glmmTMB(nh4_out_avg ~ kelp_bio_scale, 
                         family = Gamma(link = 'log'),
                         data = data_s)
plot(simulateResiduals(mod_gam_kelp5)) # three zigzags, 2 red

mod_gam_kelp6 <- glmmTMB(nh4_out_avg ~ shannon_scale, 
                         family = Gamma(link = 'log'),
                         data = data_s)
plot(simulateResiduals(mod_gam_kelp6)) # super fucked

mod_gam_kelp7 <- glmmTMB(nh4_out_avg ~ kelp_sp, 
                         family = Gamma(link = 'log'),
                         data = data_s)
plot(simulateResiduals(mod_gam_kelp7)) # no signif problems, useable

mod_gam_kelp8 <- glmmTMB(nh4_out_avg ~ depth_scale, 
                         family = Gamma(link = 'log'),
                         data = data_s)
plot(simulateResiduals(mod_gam_kelp8))  # a little wavy but ok

# mod with just the variables that didn't throw errors for nh4_outside?
mod_gam_kelp9 <- glmmTMB(nh4_out_avg ~ weight_sum_scale + kelp_sp + depth_scale, 
                         family = Gamma(link = 'log'),
                         data = data_s)
plot(simulateResiduals(mod_gam_kelp9)) # looks ok

aic <- AIC(mod_gam_kelp, mod_gam_kelp2, mod_gam_kelp3, mod_gam_kelp4, mod_gam_kelp5, mod_gam_kelp6, mod_gam_kelp7, mod_gam_kelp8, mod_gam_kelp9)
# ok so the two model with all the terms and lots of interactions are preferred


# so i accidentally ran the models with nh4_in_avg, and kelp_bio_scale*tide_scale + kelp_sp was the "best model" with low residual issues and AIC..... 

# mod with just the variables that didn't throw errors for nh4_inside?
mod_gam_kelp_in <- glmmTMB(nh4_in_avg ~ kelp_bio_scale*tide_scale + kelp_sp, 
                           family = Gamma(link = 'log'),
                           data = data_s)
summary(mod_gam_kelp_in) 
plot(simulateResiduals(mod_gam_kelp_in)) # ehhh? U-shape but less fucked
# dropping the interaction makes it worse
# forest_bio_scale makes it worse
# dropping bio mean scale makes it worse
# dropping tide makes it worse
# dropping kelp_sp gets rid of red line but all results non signif
# can't do tide*kelp_sp interaction or kelp_bio_scale*kelp_sp


# try to plot visreg output nicer
p = visreg(mod_best, "kelp_bio_scale", by="weight_sum_scale",
           overlay = TRUE, partial = FALSE, rug = FALSE,
           plot=FALSE)

ggplot(p$fit, aes(kelp_bio_scale, visregFit, 
                  fill = factor(weight_sum_scale))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha = 0.2) +
  geom_line(aes(colour = factor(weight_sum_scale)), linewidth = 1) +
  labs(linetype="weight_sum_scale", fill="weight_sum_scale", colour = "weight_sum_scale")


# Old stats -----
# I don't think this is linear I'm going to have to fit some kind of curve to it.....

# linear model though for fun
kelp_pee_mod <- lmer(in_minus_out ~ kelp_den + (1|site), data = data)
summary(kelp_pee_mod)
visreg(kelp_pee_mod)

# full model for dredging
model_all <- lmer(in_minus_out ~ den_scale + bio_tran_scale + kelp_bio_scale + area_scale + forest_bio_scale + kelp_sp + (1|site), data = data, na.action = na.fail)
summary(model_all)
visreg(model_all)

dredge <- as.data.frame(dredge(model_all)) %>%
  filter(delta < 3)

# the best model has the mean biomass of the whole forest
# best has mean bio + total forest biomass
# second best is just mean bio
# third is mean bio + total forest biomass + total area

bio_mod <- lmer(in_minus_out ~ kelp_bio_scale + (1|site), data = data)
summary(bio_mod)
visreg(bio_mod)
# Forests with more mean biomass/m2 retain more pee!

# redo dredge with the rest of the vars
model_all <- lmer(in_minus_out ~ den_scale + bio_tran_scale + kelp_bio_scale + area_scale + forest_bio_scale + kelp_sp +
                    weight_stand + rich_scale + abundance_scale + tide_scale + depth_stand + (1|site), data = data, na.action = na.fail)
summary(model_all)

dredge <- as.data.frame(dredge(model_all)) %>%
  filter(delta < 3)

# plot?
ggplot(data, aes(kelp_bio_scale, in_minus_out)) +
  geom_point() +
  geom_smooth(method = lm) 

# Let's be intelligent and build the model I think would be best
# I'd guess the biomass/m2 (density x biomass) = bio_tran_scale would matter bc that's the info about the kelp on the transect the samples were taken on
# And I'd guess the mean biomass/m2 (den x bio) of the kelp forest (how much kelp is around) will matter
mod1 <- lmer(in_minus_out ~ bio_tran_scale* kelp_bio_scale  + (1|site), data = data, na.action = na.fail)
summary(mod1)
visreg(mod1, "bio_tran_scale", by = "kelp_bio_scale")
# It looks like at high mean biomass, the relationship levels out
# So maybe there's an asymptote!

# let's try taking the log of the  in_minus_out and see if that improves fit

mod2 <- lmer(log_pee_diff ~ bio_tran_scale* kelp_bio_scale  + (1|site), data = data, na.action = na.fail)
summary(mod1)
visreg(mod1, "bio_tran_scale", by = "kelp_bio_scale")


# go back to the preferred model with just mean biomass
bio_mod2 <- lmer(log_pee_diff ~ kelp_bio_scale  + (1|site), data = data)
summary(bio_mod2)
visreg(bio_mod2)

# check resid
mod_cont_resids <- simulateResiduals(bio_mod)
plot(mod_cont_resids)
# weight fucking shit



# Use AIC to see if this improves things
AIC(bio_mod, bio_mod2, mod1)
# Yes bio_mod2 is preferred, taking the log of the response variable helps a lot



# Site level model!

kelp_mod_full <- lm(nh4_out_avg ~ weight_stand + rich_scale + abundance_scale + tide_scale + depth_stand + kelp_bio_scale + 
                      weight_stand:tide_scale + weight_stand:kelp_bio_scale +
                      rich_scale:tide_scale + rich_scale:kelp_bio_scale +
                      abundance_scale:tide_scale +
                      abundance_scale:kelp_bio_scale +
                       tide_scale:kelp_bio_scale, data_s)
summary(kelp_mod_full)

options(na.action = "na.fail")
dredge <- as.data.frame(dredge(kelp_mod_full)) %>%
  filter(delta < 3)
# nh4 avg second best mod was just tide
# nh4 inside avg second best mod was just bio mean, then just tide
# nh4 outside avg second best mod just tide

kelp_mod_mean <- lm(nh4_avg ~ weight_stand + tide_scale  + kelp_bio_scale + 
                      weight_stand:tide_scale + weight_stand:kelp_bio_scale, data_s)


summary(kelp_mod_best)

visreg(kelp_mod_best)


# Plots ----

ggplot(data, aes(kelp_bio_scale, log_pee_diff)) +
  geom_point(aes(pch = kelp_sp, colour = site_code), size =3 )+ 
  geom_smooth(method = lm, colour = "blue")+
  geom_smooth(method = loess, colour = "red")+
  
labs(y = "Inside - outside kelp forest ammonium (log(uM))", x = "Forest biomass") 
  

# each point is a single transect with a pee difference and a kelp density 
# linear model
ggplot(data, aes(BiomassM, in_minus_out)) +
  geom_point(aes(pch = kelp_sp, colour = site_code), size =3 )+ 
  geom_hline(yintercept= 0, linetype = "dashed", color = "red", size = 1.5) +
  labs(y = "Inside - outside kelp forest ammonium (uM)", x = "Kelp Density") +
  geom_smooth(method = lm,
              alpha = 0.25)

# try plotting the same data but with a better asymptote?
ggplot(data, aes(kelp_den, in_minus_out)) +
  geom_point(aes(pch = kelp_sp, colour = site_code), size =3 )+ 
  geom_hline(yintercept= 0, linetype = "dashed", color = "red", size = 1.5) +
  labs(y = "Inside - outside kelp forest ammonium (uM)", x = "Kelp Density") +
  geom_smooth(method = "loess",
              span = 1,
              alpha = 0.25)

# Percent difference 
# With a truly heinious curve, thanks loess
ggplot(data, aes(kelp_den, percent_diff)) +
  geom_point(aes(pch = kelp_sp, colour = site_code), size =3 )+ 
  geom_hline(yintercept= 0, linetype = "dashed", color = "red", size = 1.5) +
  labs(y = "Inside - outside kelp forest ammonium (uM)", x = "Kelp Density") +
  geom_smooth(method = "loess",
              span = 1,
              alpha = 0.25)

#ggsave("Output/Figures/in_out_kelp_pee.png", device = "png",
#       height = 9, width = 16, dpi = 400)



# Claire has area, density, biomass
# Got average individual biomass for each kelp, which you can scale up to the average biomass of each transect


# Old Kelp density data -----
# Also Claire's biomass estimates
# Basically it would be nice to have a single metric of how much kelp is in each forest
kcca_summary <- kcca_final %>%
  group_by(site_code) %>%
  summarise(kelp_den = mean(kelp_den),
            in_minus_out = mean(in_minus_out),
            kelp_sp = kelp_sp)