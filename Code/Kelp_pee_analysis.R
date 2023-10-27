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
         kelp_sp = ifelse(Kelp == 0, "none",
                          ifelse(Macro_5m2 == 0, "nereo", 
                                 ifelse(Nereo_5m2 == "0", "macro", "mixed")))) %>%
  left_join(names, by = "SiteName") %>% 
  replace(is.na(.), 0) %>%
  # Add the averaged site level variables from Claire!
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
data <- pee %>%
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
  left_join(tide_exchange_kcca, by = "survey_id") %>%
  # scale factors
  mutate(# kelp forest level variables
         forest_biomass = BiomassM*Area_m2,
         den_scale = c(scale(DensityM)), # make sure density is right
         bio_mean_scale = c(scale(BiomassM)),
         forest_bio_scale = c(scale(forest_biomass)),
         area_scale = c(scale(Area_m2)),
         # transect level variables
         bio_tran_scale = c(scale(Biomassm2kg)),
         den_tran_scale = c(scale(kelp_den)), # make sure density is right
         # log the pee diff?
         log_pee_diff = log(in_minus_out + 1),
         # the biomass + abundance variables
         weight_sum_stand = c(scale(weight_sum)),
         all_weighted_stand = c(scale(all_weight_weighted)),
         abundance_stand = c(scale(abundance)),
         # biodiversity variables
         rich_stand = c(scale(species_richness)),
         shannon_stand = c(scale(shannon)),
         simpson_stand = c(scale(simpson)),
         # abiotic variables I should control for
         depth_avg_stand = c(scale(depth_avg)),
         tide_stand = c(scale(avg_exchange_rate)),
         tide_cat = ifelse(avg_exchange_rate < -0.1897325, "Ebb",
                           ifelse(avg_exchange_rate < 0.1897325, "Slack", "Flood")),
                    # only slack and flood
         log_kelp = log(BiomassM + 0.001),
         log_kelp_scale = c(scale(log_kelp)),
         log_kelp_tran = log(Biomassm2kg + 0.001)) %>%
  filter(site_name != "Second Beach South")

# I don't think I might actually want to remove Second beach south.....
  # Tons of nereo but sort of an outlier site??? Claire drops this site for her stuff


# Each row is a mini transect into the kelp forest
  # That has an individual kelp density + biomass estimate + inside vs outside ammonium value

# new df for one row per site
data_s <- data %>%
  select(site, site_code, survey_id, in_out_avg, nh4_in_avg, nh4_out_avg, nh4_avg, depth_avg, avg_exchange_rate, kelp_sp, BiomassM, bio_mean_scale, weight_sum, weight_sum_stand, all_weighted_stand, abundance_stand, rich_stand, shannon_stand, simpson_stand, tide_stand, depth_avg_stand, tide_cat) %>%
  unique() %>%
  filter(site != "Wizard_I_North" | kelp_sp != "none") # just getting rid of the duplicate row, bc there wasn't kelp on one transect unique misses this one

#write_csv(data_s, "Output/Output_data/kelp_final.csv")


# not sure what this is
data_s_reduced <- data %>%
  select(site, site_code, survey_id, nh4_in_avg, depth_avg, avg_exchange_rate, BiomassM, weight_sum, species_richness, abundance, avg_exchange_rate,  depth_avg) %>%
  unique() %>%
  rename(nh4_avg = nh4_in_avg)
  

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

# Biological
  # Kelp forest
    # forest_bio_scale (site level = area*biomass*density)
    # kelp_sp (site level)
    # bio_mean_scale (site level, density*biomass)
    # log_kelp_scale log(bio_mean), also scaled
    # log_den (site level)
    # bio_tran_scale (mini transect level)
    # den_tran_scale (mini transect level)
# I'm using bio_mean_scale, best for AIC and better than log(kelp) too
# Tried including one transect level variable and it wasn't as good


  # RLS community (site level)
    # Biomass:
      # weight_sum_stand 
      # all_weighted_stand 
      # abundance 
    # Biodiversity
      # rich_stand
      # shannon_stand
      # simpson_stand
# Using weight_sum_stand and shannon, best AIC

# Abiotic to control for
  # depth_avg (depth_avg_stand)
  # avg_exchange_rate (tide_stand)
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
mod_tran <- glmmTMB(in_minus_out ~ bio_mean_scale*tide_stand*weight_sum_stand + bio_tran_scale + shannon_stand + depth_avg_stand + kelp_sp + (1|site_code), 
                         family = 'gaussian',
                         data = data,
                    na.action = na.fail)

# Dredge to find best combo of variables
# dredge <- as.data.frame(dredge(mod_tran)) %>% filter(delta < 3)

# Without second beach:
# best = bio_mean_scale*tide_stand + weight_sum_stand + shannon_stand + kelp_sp 
# second best = same + weight_sum_stand:tide_stand

# with Second beach:
# best = bio_mean_scale*tide_stand*weight_sum_stand + kelp_sp + shannon_stand 
# second best is the same + bio_tran_scale

mod_best <- glmmTMB(in_minus_out ~ bio_mean_scale*tide_stand + weight_sum_stand + shannon_stand + kelp_sp + (1|site_code), 
                    family = 'gaussian',
                    data = data)

summary(mod_best)

# Step 3: Model residuals!
plot(simulateResiduals(mod_best)) # looks fine

# Step 4: Check for collinearity of predictors

# car can't handle random effects so make a simplified mod
car::vif(lm(in_minus_out ~ bio_mean_scale + tide_stand + weight_sum_stand + shannon_stand + kelp_sp, data = data))
# Yep this looks fine!


# Notes from stats beers!
#  Double check that the negative in – out sites aren’t the places where the transect was on the other side of the forest.
# Maybe try to PCA????? see if there's clustering?
  

# Try to plot in minus out ------
visreg(mod_best, "bio_mean_scale", by = "tide_stand", overlay =TRUE)
# so when the tide comes in the difference in ammonium in thicker kelp forests is even greater than you'd expect at slack tide


# generate df for plotting
df <- confint(mod_best, level = 0.95, method = c("wald"), component = c("all", "cond", "zi", "other"), estimate = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(variable = rowname,
         lower_CI = `2.5 %`,
         upper_CI = `97.5 %`,
         estimate = Estimate) %>%
  head(- 1)  %>%
  mutate(variable = factor(as.factor(variable), 
                           levels = c("(Intercept)", "bio_mean_scale", "tide_stand", "weight_sum_stand",  "shannon_stand", "kelp_spmixed", "kelp_spnereo", "kelp_spnone", "bio_mean_scale:tide_stand"),
                           labels = c("Intercept", "Kelp biomass", "Tide", "Animal biomass",  "Shannon", "Mixed kelp", "Nereo", "None", "Kelp:tide"))
  )

# Coefficient plot
pal9 <- viridis::viridis(9)

ggplot(df, aes(x = estimate, y = (variable), xmin = lower_CI, xmax = upper_CI, colour = variable)) +
  geom_point(size = 10) +
  geom_errorbar(width = 0, linewidth = 3) +
  geom_vline(xintercept=0, color="white", linetype="dashed") +
  labs(x = "Coefficient", y = " ") +
  scale_y_discrete(limits = rev(levels(df$variable))) +
  theme_black() +
  theme(legend.position = "none") + 
  scale_colour_manual(values = pal9)

# ggsave("Output/Figures/kelp_in_out_mod_coeff.png", device = "png", height = 9, width = 16, dpi = 400)


# remake the shitty asymptote fig with the log(kelp) model and actual predictions, with geom_point(aes(colour = tide, size = animal weight))! See what that looks like???

# make new mod with untransformed vars
mod_plot <- glmmTMB(in_minus_out ~ BiomassM*avg_exchange_rate + weight_sum + shannon + kelp_sp + (1|site_code), 
                   family = 'gaussian',
                   data = data)
summary(mod_plot)
plot(simulateResiduals(mod_log))

# what should the tide be set to?
ggplot(data, aes(avg_exchange_rate, 1)) + geom_point()
slack <- mean(data$tide_stand[data$avg_exchange_rate<0])
flood <- mean(data$tide_stand[data$avg_exchange_rate>0])

# make predictions
# create range vector
v <- seq(-6.97, 0.65, length.out = 100)

v <- c(slack, flood)

# now make predictions
predict <- ggpredict(mod_best, terms = c("bio_mean_scale", "tide_stand [v]")) %>% 
  mutate(bio_mean_scale = x,
         tide = ifelse(group == "-0.750048942265706", "slack", "flood")) %>%
  filter(tide != "flood" | bio_mean_scale < 0) %>%
  filter(tide != "flood" | bio_mean_scale > -0.78)
  
# now plot these predictions

ggplot() + 
  geom_point(data = data %>%
               mutate(tide = ifelse(avg_exchange_rate < 0, "slack", "flood")) , 
             aes(x = bio_mean_scale, y = in_minus_out,
             colour = tide,
             fill = tide,
             size = weight_sum_stand), alpha = 0.8) +
  geom_line(data = predict,
            aes(x = bio_mean_scale, y = predicted, lty = tide, colour = tide),
            linewidth = 1.5) +
  geom_ribbon(data = predict,
              aes(x = bio_mean_scale, y = predicted, 
                  fill = tide,
                  ymin = conf.low, ymax = conf.high), 
              alpha = 0.2) +
  geom_hline(yintercept= 0, linetype = "dashed", color = "white", linewidth = 0.5) +
labs(y = "Delta ammonium (uM)", x = "Kelp Biomass (kg/m2)",
     colour = "Tide", fill = "Tide",
     lty = "Tide",
     size = "Animals (kg)") +
  guides(size = guide_legend(override.aes = 
                              list(colour = "white"))) +
  theme_black()

#ggsave("Output/Figures/kelp_in_out_mod_predict.png", device = "png", height = 9, width = 16, dpi = 400)

  


# try again on transect level ---
# mod for plotting indiv transects?
mod_tran <- glmmTMB(in_minus_out ~ log_kelp_tran*avg_exchange_rate + weight_sum + shannon_stand + kelp_sp + (1|site_code), 
                    family = 'gaussian',
                    data = data)

# make predictions
# create range vector
v2 <- seq(-6.97, 1.2, length.out = 100)
# now make predictions
predict2 <- ggpredict(mod_tran, terms = "log_kelp_tran [v2]") %>% 
  mutate(Biomassm2kg = exp(x) - 0.001) %>%
  filter(predicted > -1)

# now plot these predictions
ggplot() + 
  geom_point(data = data, aes(x = Biomassm2kg, y = in_minus_out,
                              colour = avg_exchange_rate,
                              fill = avg_exchange_rate,
                              size = weight_sum,
                              pch = kelp_sp), alpha = 0.8) +
  #now plot the  model output
  geom_ribbon(data = predict2,
              aes(x = Biomassm2kg, y = predicted, ymin = conf.low, ymax = conf.high), 
              alpha = 0.3) +
  geom_line(data = predict2,
            aes(x = Biomassm2kg, y = predicted),
            linewidth = 1.5) +
  theme_classic() +
  geom_hline(yintercept= 0, linetype = "dashed", color = "red", linewidth = 0.5) +
  labs(y = "Inside - outside kelp forest ammonium (uM)", x = "Kelp Biomass (kg/m2)",
       colour = "Tide exchange (m/s)", fill = "Tide exchange (m/s)",
       size = "Animal biomass (kg)",
       pch = "Kelp species") +
  ylim(c(-1, 1)) +
  scale_shape_manual(values = c(21, 24, 22, 25))




# Stats: site to site variation ----

# Step 0: Ask my question
  # Is there a link between biodiversity or biomass and ammonium concentration? 

# Step 1: Choose a distribution

# response variables:
  # nh4_avg (average of in and out)
  # nh4_in_avg (inside kelp avg)
  # nh4_out_avg (outside kelp avg)
    # most comperable to the RLS blitz work?

ggplot(data_s, aes(x = nh4_out_avg)) +
  geom_histogram(bins = 30) 
# Either way we're using Gamma!

# Step 2: Choose your predictors

# Biological predictors:
  # 1) Animal biomass
    # a) weight_sum = total wet weight of all the animals observed on RLS transects
    # b) all_weight_weighted = weighted fish weight + invert weight
    # c) abundance
      # I'm going with weight_sum
  # 2) Kelp forest
    # a) forest_bio_scale (site level = area*biomass*density)
    # b) bio_mean_scale (site level, density*biomass)
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

# Build full model with gaussian distribution
mod_norm_kelp <- glmmTMB(nh4_out_avg ~ weight_sum_stand*bio_mean_scale*tide_stand + shannon_stand + kelp_sp + depth_avg_stand,
                         family = 'gaussian',
                         data = data_s)
summary(mod_norm_kelp)

# Build full model with gamma distribution
mod_gam_kelp <- glmmTMB(nh4_in_avg ~ weight_sum_stand*bio_mean_scale*tide_stand + shannon_stand + kelp_sp + depth_avg_stand, 
                        family = Gamma(link = 'log'),
                        data = data_s)
summary(mod_gam_kelp)

# Step 3: Model residuals 
  # Compare the two distributions 

# Gaussian
mod_norm_kelp_resids <- simulateResiduals(mod_norm_kelp)
plot(mod_norm_kelp_resids) 
# no signif problems found......

# Gamma
mod_gam_kelp_resids <- simulateResiduals(mod_gam_kelp)
plot(mod_gam_kelp_resids) # not technically wrong but a little funky

AIC(mod_norm_kelp, mod_gam_kelp)
# Gam is better but there's some S-shape stuff going on 


# Step 4: Check for collinearity of predictors

# car can't handle random effects so make a simplified mod
car::vif(lm(nh4_out_avg ~ weight_sum_stand + bio_mean_scale + tide_stand + shannon_stand + kelp_sp + depth_avg_stand, data = data_s))
# all good


# try to plot this to see what the heck is happening
data_s %>%
  mutate(kelp_cat = case_when(BiomassM < 0.565827 ~ "Low", # this does include the 1 no kelp site
                              BiomassM < 1.131654 ~ "Mid",
                              BiomassM >= 1.131654 ~ "Tons")) %>%
  
ggplot(aes(weight_sum, BiomassM, colour = nh4_out_avg)) +
  geom_point() +
  geom_smooth(method = lm)
# In low kelp density when the tide is flooding, instead of getting lower the nh4 outside the kelp forest is increasing? only 2 data points through
# but for mid kelp density i mean I don't think I can say anything there's only 4 points

# Graveyard ------

# try to plot visreg output nicer
p = visreg(mod_best, "bio_mean_scale", by="weight_sum_stand",
           overlay = TRUE, partial = FALSE, rug = FALSE,
           plot=FALSE)

ggplot(p$fit, aes(bio_mean_scale, visregFit, 
                  fill = factor(weight_sum_stand))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha = 0.2) +
  geom_line(aes(colour = factor(weight_sum_stand)), linewidth = 1) +
  labs(linetype="weight_sum_stand", fill="weight_sum_stand", colour = "weight_sum_stand")


# Old stats -----
# I don't think this is linear I'm going to have to fit some kind of curve to it.....

# linear model though for fun
kelp_pee_mod <- lmer(in_minus_out ~ kelp_den + (1|site), data = data)
summary(kelp_pee_mod)
visreg(kelp_pee_mod)

# full model for dredging
model_all <- lmer(in_minus_out ~ den_scale + bio_tran_scale + bio_mean_scale + area_scale + forest_bio_scale + kelp_sp + (1|site), data = data, na.action = na.fail)
summary(model_all)
visreg(model_all)

dredge <- as.data.frame(dredge(model_all)) %>%
  filter(delta < 3)

# the best model has the mean biomass of the whole forest
# best has mean bio + total forest biomass
# second best is just mean bio
# third is mean bio + total forest biomass + total area

bio_mod <- lmer(in_minus_out ~ bio_mean_scale + (1|site), data = data)
summary(bio_mod)
visreg(bio_mod)
# Forests with more mean biomass/m2 retain more pee!

# redo dredge with the rest of the vars
model_all <- lmer(in_minus_out ~ den_scale + bio_tran_scale + bio_mean_scale + area_scale + forest_bio_scale + kelp_sp +
                    weight_stand + rich_stand + abundance_stand + tide_stand + depth_stand + (1|site), data = data, na.action = na.fail)
summary(model_all)

dredge <- as.data.frame(dredge(model_all)) %>%
  filter(delta < 3)

# plot?
ggplot(data, aes(bio_mean_scale, in_minus_out)) +
  geom_point() +
  geom_smooth(method = lm) 

# Let's be intelligent and build the model I think would be best
# I'd guess the biomass/m2 (density x biomass) = bio_tran_scale would matter bc that's the info about the kelp on the transect the samples were taken on
# And I'd guess the mean biomass/m2 (den x bio) of the kelp forest (how much kelp is around) will matter
mod1 <- lmer(in_minus_out ~ bio_tran_scale* bio_mean_scale  + (1|site), data = data, na.action = na.fail)
summary(mod1)
visreg(mod1, "bio_tran_scale", by = "bio_mean_scale")
# It looks like at high mean biomass, the relationship levels out
# So maybe there's an asymptote!

# let's try taking the log of the  in_minus_out and see if that improves fit

mod2 <- lmer(log_pee_diff ~ bio_tran_scale* bio_mean_scale  + (1|site), data = data, na.action = na.fail)
summary(mod1)
visreg(mod1, "bio_tran_scale", by = "bio_mean_scale")


# go back to the preferred model with just mean biomass
bio_mod2 <- lmer(log_pee_diff ~ bio_mean_scale  + (1|site), data = data)
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

kelp_mod_full <- lm(nh4_out_avg ~ weight_stand + rich_stand + abundance_stand + tide_stand + depth_stand + bio_mean_scale + 
                      weight_stand:tide_stand + weight_stand:bio_mean_scale +
                      rich_stand:tide_stand + rich_stand:bio_mean_scale +
                      abundance_stand:tide_stand +
                      abundance_stand:bio_mean_scale +
                       tide_stand:bio_mean_scale, data_s)
summary(kelp_mod_full)

options(na.action = "na.fail")
dredge <- as.data.frame(dredge(kelp_mod_full)) %>%
  filter(delta < 3)
# nh4 avg second best mod was just tide
# nh4 inside avg second best mod was just bio mean, then just tide
# nh4 outside avg second best mod just tide

kelp_mod_mean <- lm(nh4_avg ~ weight_stand + tide_stand  + bio_mean_scale + 
                      weight_stand:tide_stand + weight_stand:bio_mean_scale, data_s)


summary(kelp_mod_best)

visreg(kelp_mod_best)


# Plots ----

ggplot(data, aes(bio_mean_scale, log_pee_diff)) +
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

# Ammonium at the site level?
data %>%
  mutate(site = fct_reorder(site, nh4_avg, .fun='median')) %>%
  ggplot() +
  geom_point(aes(x = reorder(site, nh4_avg), nh4_outside, 
                 colour = "blue")) +
  geom_point(aes(x = reorder(site, nh4_avg), nh4_inside, 
                 colour = "green")) +
  labs(y = "NH4 concentration", x = "Site") + 
  theme(axis.text.x=element_text(angle = 65, hjust = 1))+
  scale_color_identity(name = "Sample",
                       breaks = c("blue", "green"),
                       labels = c("Outside kelp", "Inside kelp"),
                       guide = "legend")
# Ok this is sort of neat, you can see the base ammonium levels and then how different they are, by site. might be neat to arrange these with kelp density instead of site on the x



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