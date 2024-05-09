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
library(patchwork)

# set theme and load functions
theme_set(theme_bw())
source("Code/Functions.R")

# Kelp data ----

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
           Nereo_5m2 == 0 ~ "macro",
           TRUE ~ as.character("macro"))) # everything else = macro
         ) %>% 
  left_join(names, by = "SiteName") %>% 
  replace(is.na(.), 0) %>%
  # Add the averaged site level variables from Claire!
  left_join(read_csv("Data/Team_kelp/Output_data/kelpmetrics_2022.csv"), by = "SiteName") %>%
  group_by(SiteName) %>%
  mutate(BiomassM = ifelse(kelp_sp == "none", 0, mean(Biomassm2kg)),
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
kelp_rls <- read_csv("Data/Team_Kelp/RLS_KCCA_2022.csv") %>%
  # processing required to get this df into the RLS data format
  as.data.frame() %>%
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
  # join phylo names from the rls blitz data
  left_join(read_csv("Output/Output_data/rls_phylo.csv"), by = "species_name") %>%
  clean_phylo_names() %>% # function to fix naming errors
  filter(species_name != "Tonicella spp.") %>% # remove chitons
  filter(species_name != "Cryptochiton stelleri") %>%
  # filter fish seen on wrong survey
  filter(!(family == "Gobiidae" & method == 1)) %>% 
  filter(!(family == "Cottidae" & method == 1)) %>% 
  filter(!(family == "Hexagrammidae" & method == 2)) %>%
  # Pivot longer for biomass
  pivot_longer( cols = `0`:`400`, names_to = "size_class", values_to = "total") %>%
  drop_na(total) %>%
  filter(total > 0) %>%
  select(-Total) %>%
  mutate(size_class = as.numeric(size_class),
         # correct for rectangle area
         total = case_when(method == 1 ~ total/500,
                           method == 2 ~ total/100)) %>%
  length_to_weight() %>% # fish length to weight function
  home_range()  # calculate each fish's home range function


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
  scale_vars() %>%
  mutate(exp_kelp = exp(BiomassM),
         exp_abund = exp(abundance),
         kelp_cat = factor(as.factor(
           case_when(
             kelp_bio_scale < -0.5 ~ "Low",
             kelp_bio_scale < 1 ~ "Mid",
             kelp_bio_scale > 1 ~ "High")),
           levels = c("Low", "Mid", "High")))

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
      # abundance_scale
    # Biodiversity
      # rich_scale
      # shannon_scale
      # simpson_scale
# Using abundance_scale and simpson_scale, best AIC

# Abiotic to control for
  # depth_avg (depth_scale)
  # avg_exchange_rate (tide_scale)
# Tide doesn't matter but depth does????

# Random variables
  # side_code: if I use mini transect as the level of study than I need a random effect of site

# Interactions:
  # interaction between animal biomass:tide_exchange
  # interaction between kelp:tide_exchange
  # interaction between animal biomass:kelp
  # maybe a triple interaction between kelp:animal_biomass:tide??


# Build full model with gaussian distribution
# This is the transect level model
# One row per transect, 3 per site
mod_tran <- glmmTMB(in_minus_out ~ kelp_sp + 
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

# best = kelp_bio_scale*tide_scale + weight_sum_scale + shannon_scale + kelp_sp 
# new best = abundance_scale + depth_scale + kelp_bio_scale + kelp_sp + tide_scale + abundance_scale:kelp_bio_scale + abundance_scale:tide_scale + kelp_bio_scale:tide_scale

mod_best <- glmmTMB(in_minus_out ~ kelp_sp + kelp_bio_scale + abundance_scale + depth_scale + tide_scale + abundance_scale:kelp_bio_scale + abundance_scale:tide_scale + kelp_bio_scale:tide_scale + (1|site_code), 
                    family = 'gaussian',
                    data = data )

plot(simulateResiduals(mod_best)) # looks fine
summary(mod_best)


# use my stupid brain to think of a good model
mod_in_out <- glmmTMB(in_minus_out ~ kelp_sp + 
                        kelp_bio_scale*tide_scale +
                        kelp_bio_scale*weight_sum_scale +
                        weight_sum_scale*tide_scale +
                        shannon_scale + depth_scale + 
                        (1|site_code),
                      family = 'gaussian',
                      data = data) 
# when I remove no kelp sites the coeffs hardly change

plot(simulateResiduals(mod_in_out)) # looks ok!
summary(mod_in_out)

# Let's look at those interactions some more:

# Kelp:Tide
visreg(mod_in_out, "kelp_bio_scale", by = "tide_scale", overlay=TRUE)
# the range for low and high tide exchange rates is pretty small
# Hannah said it's ok

# Abundance:Kelp
visreg(mod_in_out, "weight_sum_scale", by = "kelp_bio_scale", overlay=TRUE)
# looks like an ok spread, 2 groups only have 3 sites, tho
# this basically looks like: at high kelp or high animals, you get a large delta nh4. But if one of the variables is med/low, increasing the other will increase delta
# I think this is my "asympote"

# Abundance:tide
visreg(mod_in_out, "weight_sum_scale", by = "tide_scale", overlay=TRUE)
# this looks fine

# run the model without an intercept
mod_in_out2 <- glmmTMB(in_minus_out ~ -1 + kelp_sp + 
                         kelp_bio_scale*tide_scale +
                         kelp_bio_scale*weight_sum_scale +
                         weight_sum_scale*tide_scale +
                         shannon_scale + depth_scale + 
                         (1|site_code),
                       family = 'gaussian',
                       data = data)
summary(mod_in_out2)

# what if I add kelp sp interaction
# I can't do a triple kelp_sp*kelp_bio_scale*tide_scale interaction bc no nereo flood
# I want to know if the tide:kelp interaction is the same for nereo and macro
    # there is no tide:nereo interaction bc nereo only measured at slack tide
    # and there's no kelp:tide interaction for no kelp bc there's no kelp biomass
# basically the model is guessing at what's happening with nereo based on what's happening with macro and that's find


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
pal8 <- viridis::viridis(8)
pal_k <- viridis::viridis(10)
pal2 <- c(pal_k[8], pal_k[5])
pal <- viridis::viridis(10)
pal3 <- c(pal[10], pal[8], pal[5])

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
                           levels = c("kelp_spmacro", "kelp_spnereo", "kelp_spnone", "kelp_bio_scale", "weight_sum_scale", "depth_scale", "tide_scale", "shannon_scale", "kelp_bio_scale:tide_scale", "kelp_bio_scale:weight_sum_scale", "tide_scale:weight_sum_scale"),
                           labels = c("Macro", "Nereo", "None", "Kelp biomass", "Animal biomass", "Depth", "Tide", "Biodiversity", "Kelp:tide", "Kelp:animals", "Tide:animals")),
         se = (upper_CI - estimate)/1.96
  )

# Use function to plot coefficients
kelp_coeff_plot <- coeff_plot(coeff_df = df,
                              pal = pal12) +
  place_label("(a)")

#ggsave("Output/Pres_figs/Fig4a.png", device = "png", height = 8, width = 12, dpi = 400)

# plot just cont vars
kelp_coeff_plot <- coeff_plot(coeff_df = (df %>% tail(-3) %>% droplevels()),
                              pal = pal8) +
  place_label("(a)")


# Plot kelp species predictions ----
# make predictions for each kelp sp at the group mean kelp bio for that species
d <- data %>%
  mutate(kelp = ifelse(kelp_sp == "none", "none", "kelp")) %>%
  group_by(kelp) %>%
  summarise(mean = mean(kelp_bio_scale))

v <- c(d$mean[1], d$mean[2])
  
predict_sp <- ggpredict(mod_in_out2, terms = c("kelp_bio_scale[v]", "kelp_sp")) %>% 
  dplyr::rename(kelp_bio_scale = x,
                in_minus_out = predicted,
                kelp_sp = group
         ) %>%
  mutate(kelp = ifelse(kelp_sp == "none", "none", "kelp")) %>%
  left_join(d, by = "kelp") %>%
  filter(mean ==  kelp_bio_scale)

# old way, model estimates at mean levels of all other variables
#predict_sp2 <- ggpredict(mod_in_out2, terms = "kelp_sp") %>% 
#  dplyr::rename(kelp_sp = x,
#                in_minus_out = predicted) 

sp_labs <- c(macro = "Macro", nereo = "Nereo", none = "None")

set.seed(234444)
kelp_sp_plot <- 
  dot_whisker(sum_data = predict_sp, 
              all_data = data,
              x_var = kelp_sp,
              y_var = in_minus_out,
              pch_var = kelp_sp,
              labels = sp_labs,
              pal = c(pal12[3], pal12[2],pal12[1])) +
  labs(y = expression(paste(Delta, " Ammonium ", (mu*M))),
       x = "Kelp species") +
    geom_hline(yintercept = 0, lty = "dashed") +
    place_label("(b)") +
  ylim(c(-0.79, 1.05))

set.seed(234444)
kelp_sp_plot <- 
alt_dot_whisker(data, x_var = kelp_sp,
                y_var = in_minus_out,
                group = kelp_sp,
                labels = sp_labs,
                pal = c(pal12[3], pal12[2],pal12[1]))+
  labs(y = expression(paste(Delta, " Ammonium ", (mu*M))),
       x = "Kelp species") +
  geom_hline(yintercept = 0, lty = "dashed") +
  place_label("(b)")

  
# Plot kelp_biomass predictions ----
predict_kelp_bio <- ggpredict(mod_in_out2, terms = "kelp_bio_scale") %>% 
  dplyr::rename(kelp_bio_scale = x)

kelp_pred_plot <- plot_pred(raw_data = data,
                            predict_data = predict_kelp_bio, 
                            plot_type = "new_kelp",
                            x_var = kelp_bio_scale, y_var = in_minus_out, 
                            pch_var = kelp_sp,
                            x_axis_lab = expression(paste("Kelp biomass (kg/m"^2,")")),
                            pal = pal12[4],
                            lty_var = pal12[4]) +
  place_label("(c)")


# Plot animal abundance predictions ----
predict_abund <- ggpredict(mod_in_out2, terms = "weight_sum_scale") %>% 
  dplyr::rename(weight_sum_scale = x)

abund_pred_plot <- plot_pred(raw_data = data,
                             predict_data = predict_abund, 
                             plot_type = "new_kelp",
                             x_var = weight_sum_scale, y_var = in_minus_out,
                             pch_var = kelp_sp,
                             x_axis_lab = "Animal Biomass (kg)",
                             pal = pal12[5],
                             lty_var = pal12[5]) +  
  place_label("(d)")

# Plot depth predictions ---- 
predict_depth <- ggpredict(mod_in_out2, terms = "depth_scale") %>% 
  dplyr::rename(depth_scale = x)

depth_pred_plot <- plot_pred(raw_data = data,
                             predict_data = predict_depth, 
                             plot_type = "new_kelp",
                             x_var = depth_scale, y_var = in_minus_out, 
                             pch_var = kelp_sp,
                             x_axis_lab = "Depth",
                             pal = pal12[6],
                             lty_var = pal12[6]) +  
  place_label("(e)")

# Plot kelp biomass vs tide interaction ------
visreg(mod_in_out, "kelp_bio_scale", by = "tide_scale", overlay=TRUE,
       xlab = "Kelp biomass (kg/m2)")

# create range vector
tide_means_kelp <- data %>%
  group_by(tide_cat) %>%
  summarise(tide = mean(tide_scale))

v2 <- c(as.numeric(tide_means_kelp[2,2]), as.numeric(tide_means_kelp[1,2]))

# now make predictions
predict_kelp <- ggpredict(mod_in_out2, terms = c("kelp_bio_scale", "tide_scale [v2]")) %>% 
  mutate(kelp_bio_scale = x,
         tide_cat = factor(as.factor(ifelse(
           group == as.character(tide_means_kelp[1,2]), "Slack", "Flood")),
           levels = c("Ebb", "Slack", "Flood"))) %>%
  filter(tide_cat != "Flood" | kelp_bio_scale < 0.21) %>%
  filter(tide_cat != "Flood" | kelp_bio_scale > -0.6)


# now plot these predictions
kelp_tide_int_plot <- 
  plot_pred(raw_data = (data %>% mutate(
    tide = ifelse(avg_exchange_rate < 0, "Slack", "Flood"))),
    predict_data = predict_kelp, 
    plot_type = "kelp",
    x_var = kelp_bio_scale, y_var = in_minus_out, 
    lty_var = tide_cat,
     pal = pal2) +
  theme(legend.position = "none") +
#  theme(legend.position = c(0.8, 0.17)) +
  guides(size = "none") +
#  ylim(c(-0.5, 0.8))
  place_label("(c)")


# Plot kelp biomass vs abundance interaction -----
visreg(mod_in_out2, "weight_sum_scale", by = "kelp_bio_scale", 
       overlay=TRUE,
       xlab = "Animal biomass (kg/m2)") # looks better
# create range vector based on visreg numbers

hist(data$kelp_bio_scale)

d <- data %>%
  group_by(kelp_cat) %>%
  summarise(mean = mean(kelp_bio_scale),
            min = min(weight_sum_scale) - 0.2,
            max = max(weight_sum_scale) + 0.5)

v6 <- c(d$mean[1], d$mean[2], d$mean[3])

v3 <- c(-1.176, -0.47, 1.897)

v4 <- c(-0.822849, 0, 1.813896)
v4 <- c(-0.822849, 0, 1.5)
v5 <- c(-1, 0, 1.75) # based on eyeballed means of the three bins looking at the hist


# now make predictions
predict_abund_kelp <- ggpredict(mod_in_out2, terms = c("weight_sum_scale", "kelp_bio_scale [v6]")) %>% 
  mutate(weight_sum_scale = x,
         kelp_cat = factor(as.factor(case_when(
           group == d$mean[1] ~ "Low",
           group == d$mean[2] ~ "Mid",
           group == d$mean[3] ~ "High")),
                           levels = c("Low", "Mid", "High"))) %>%
  left_join(d, by = "kelp_cat") %>%
  rowwise() %>%
  filter(between(weight_sum_scale, min, max))
  

# now plot these predictions
abund_kelp_int_plot <- 
  plot_pred(raw_data = data,
            predict_data = predict_abund_kelp, 
            plot_type = "new_kelp",
            x_var = weight_sum_scale, y_var = in_minus_out, 
            lty_var = kelp_cat,
            x_axis_lab = expression(paste("Animal biomass (kg/m"^2,")")),
            pal = pal3) +
  guides(size = "none") +
  theme(#legend.position = c(0.58, 0.95), legend.direction="horizontal",
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  labs(colour = "Kelp", fill = "Kelp", lty = "Kelp") +
  place_label("(d)") +
  ylim(c(-0.79, 1.05))


# Plot abundance vs tide interaction ------
visreg(mod_in_out, "weight_sum_scale", by = "tide_scale", overlay=TRUE,
       xlab = "Animal biomass (kg/m2)")
visreg(mod_in_out, "tide_scale", by = "weight_sum_scale", overlay=TRUE)

# create range vector
tide_means_kelp <- data %>%
  group_by(tide_cat) %>%
  summarise(tide = mean(tide_scale))

v2 <- c(as.numeric(tide_means_kelp[2,2]), as.numeric(tide_means_kelp[1,2]))

# range of animal vs tide
d2 <- data %>%
  group_by(tide_cat) %>%
  summarise(min = min(weight_sum_scale) - 0.3,
            max = max(weight_sum_scale) + 0.5)

# now make predictions
predict_abund_tide <- ggpredict(mod_in_out2, terms = c("weight_sum_scale", "tide_scale [v2]")) %>% 
  mutate(weight_sum_scale = x,
         tide_cat = factor(as.factor(ifelse(
           group == as.character(tide_means_kelp[1,2]), "Slack", "Flood")),
           levels = c("Ebb", "Slack", "Flood"))) %>%
  left_join(d2, by = "tide_cat") %>%
  rowwise() %>%
  filter(between(weight_sum_scale, min, max))

# now plot these predictions
abund_tide_int_plot <- 
  plot_pred(raw_data = (data %>% mutate(
                            tide = ifelse(avg_exchange_rate < 0, "Slack", "Flood"))),
            predict_data = predict_abund_tide, 
            plot_type = "new_kelp",
            x_var = weight_sum_scale, y_var = in_minus_out, 
            lty_var = tide_cat,
            x_axis_lab = expression(paste("Animal biomass (kg/m"^2,")")),
            pal = pal2) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())+
  place_label("(e)") +
  guides(size = "none")

# Put them allll together???? ----

# plot indiv relationships
plots <- kelp_sp_plot + kelp_pred_plot + abund_pred_plot + depth_pred_plot

kelp_coeff_plot + plots


# plot coeffs + interactions
squish <- theme(axis.title.y = element_text(margin = margin(r = -200, unit = "pt")))

kelp_coeff_plot/ ((kelp_sp_plot + squish) +
  abund_kelp_int_plot + 
  (kelp_tide_int_plot + squish) +
  abund_tide_int_plot )



# ggsave("Output/Pub_figs/Fig4.png", device = "png", height = 16, width = 14, dpi = 400)
  
# comparing the coeff plots, making depth a random effect doesn't really change anything REAL. It just makes the no kelp coeff LOOK like its not signif diff from 0, but when you estimate the delta at the mean no kelp biomass (0) it's basically the same estimate as the model with depth
  
# makes me wonder if we just forgo the coeff plot and just show the four prediction plots: kelp sp, kelp bio, animals, depth??


# Figure 4 for pub with white -----
kelp_coeff_plot + kelp_pred_plot

#ggsave("Output/Pub_figs/Fig4b.png", device = "png", height = 9, width = 16, dpi = 400)




# IDK if this is useful anymore Plot kelp species vs kelp biomass vs nh4 ------
data_kelp <- data %>%
  filter(kelp_sp != "none") %>%
  mutate(var = case_when(kelp_sp == "macro" & tide_cat == "Slack" ~ "macro_slack",
                         kelp_sp == "macro" & tide_cat == "Flood" ~ "macro_flood",
                         kelp_sp == "nereo" & tide_cat == "Slack" ~ "nereo_slack",
                         kelp_sp == "nereo" & tide_cat == "Flood" ~ "nereo_flood"))

mod_sp <- glmmTMB(in_minus_out ~ kelp_sp + kelp_bio_scale + abundance_scale + 
                    shannon_scale + depth_scale +
                    (1|site_code),
                  family = 'gaussian',
                  data = data_kelp)

plot(simulateResiduals(mod_sp)) # looks fine
summary(mod_sp)


visreg(mod_sp, "kelp_bio_scale", by = "abundance_scale", overlay = TRUE)

# max min
max_min_sp <- data %>%
  group_by(kelp_sp) %>%
  summarise(min = min(kelp_bio_scale) - 0.1,
            max = max(kelp_bio_scale) + 0.1)

# now make predictions
predict_kelp_sp <- ggpredict(mod_in_out2, terms = c("kelp_bio_scale", "kelp_sp")) %>% 
  mutate(kelp_bio_scale = x,
         kelp_sp = as.factor(group)) %>%
  left_join(max_min_sp, by = "kelp_sp") %>%
  rowwise() %>%
  filter(between(kelp_bio_scale, min, max)) 
 #filter(tide_cat != "Flood" | kelp_bio_scale < 0.21) %>%
 #filter(tide_cat != "Flood" | kelp_bio_scale > -0.6)

# plot predictions
ggplot() + 
  geom_point(data = data, 
             aes(x = kelp_bio_scale, y = in_minus_out, 
                 colour = kelp_sp, fill = kelp_sp,
                 pch = kelp_sp), 
             alpha = 0.8) +
  geom_line(data = predict_kelp_sp,
            aes(x = kelp_bio_scale, y = predicted, 
                colour = kelp_sp, lty = kelp_sp),
            linewidth = 2) +
  geom_ribbon(data = predict_kelp_sp,
              aes(x = kelp_bio_scale, y = predicted, fill = kelp_sp,
                  ymin = conf.low, ymax = conf.high), 
              alpha = 0.15) +
  stat_summary(data = (data %>% filter(kelp_sp == "none")), 
               aes(x = kelp_bio_scale, y = in_minus_out, group = kelp_sp), 
               fun = "mean",
               fun.data = "mean_cl_boot") +
  labs(colour = "Kelp species", fill = "Kelp species", lty = "Kelp species", pch = "Kelp species") +
  theme_white() +
  scale_colour_manual(values = (pal3)) +
  scale_fill_manual(values = (pal3)) +
  geom_hline(yintercept= 0, linetype = "dashed", color = "black", linewidth = 0.5)+
  guides(lty = guide_legend(override.aes = list(linewidth = 0.5)),
         size = guide_legend(override.aes = list(colour = "black")),
         colour = guide_legend(override.aes = list(size = 2)))  +
    labs(y = expression(paste(Delta, " Ammonium ", (mu*M))), 
         x = expression(paste("Kelp biomass (kg/m"^2,")"))) +
  theme(legend.position = c(0.75,0.25))

# now plot these predictions
plot_pred(raw_data = data_kelp,
          predict_data = predict_kelp_sp, 
          plot_type = "kelp",
          x_var = kelp_bio_scale, y_var = in_minus_out, 
          lty_var = kelp_sp,
          size_var = weight_sum, 
          pal = pal2)


# so basically for Macrocystis, there's a positive interaction between kelp biomass and tide exchange. But that's only driven by two high in - out sites and one low in - out site at flood tide.... so is that interaction real


# Family level ------

# No zeros for abundance!

data_fam_no0_a <- data %>%
  select(site, site_code, in_out_avg, kelp_bio_scale, kelp_sp, depth_scale, weight_sum_scale, abundance_scale, shannon_scale, tide_scale, tide_cat) %>%
  unique() %>%
  left_join((kelp_rls %>%
               mutate(family = as.factor(family))%>%
               group_by(site_code, method, phylum, family) %>% # if i want methods split up, add it back here
               summarise(total_fam = sum(total),
                         weight_fam_sum_g = 1000*sum(weight_size_class_sum))),
            by = "site_code") %>%
  mutate(abund_fam_scale = c(scale(total_fam)),
         weight_fam_scale = c(scale(weight_fam_sum_g)))


# what are the top families?

# df for just the surveys I'm looking at
kelp_sm <- data %>%
  select(site_code) %>%
  unique() %>%
  left_join(kelp_rls, by = "site_code")

# CHOOSE THE FAMILIES TO INCLUDE
# rank of families by total abundance (density)
kelp_fam_list_total <- kelp_sm %>%
  group_by(family) %>%
  summarise(sum = sum(total)) %>%
  drop_na(family) %>%
  arrange(desc(sum)) %>%
  transmute(family = family, 
            sum_total = sum,
            rank_total = 1:43)

# but which families show up on the most transects?
kelp_fam_list_count <- kelp_sm %>%
  select(site_code, family) %>%
  unique() %>%
  count(family) %>%
  drop_na(family) %>%
  arrange(desc(n)) %>%
  transmute(family = family, 
            sum_count = n,
            rank_count = 1:43)

# rank of families by biomass
kelp_fam_list_bio <- kelp_sm %>%
  group_by(family) %>%
  summarise(sum = sum(weight_size_class_sum)) %>%
  drop_na(family) %>%
  arrange(desc(sum)) %>%
  transmute(family = family, 
            sum_bio = sum,
            rank_bio = 1:43)

# join all that fam stuff together
kelp_fam_big_list <- kelp_fam_list_total %>%
  left_join(kelp_fam_list_count, by = "family") %>%
  left_join(kelp_fam_list_bio, by = "family") %>%
  mutate(rank_all = (rank_total + rank_count + rank_bio)/3)

# OK so how do I make sense of this
# plot all rankings together
ggplot(kelp_fam_big_list) +
  geom_point(aes(y = reorder(family, -rank_all), 
                 x = rank_total, colour = "density"), alpha = 0.5) +
  geom_point(aes(y = reorder(family, -rank_all), 
                 x = rank_count, colour = "count"), alpha = 0.5) +
  geom_point(aes(y = reorder(family, -rank_all), 
                 x = rank_bio, colour = "biomass"), alpha = 0.5) +
  geom_point(aes(y = reorder(family, -rank_all), 
                 x = rank_all, colour = "all"), alpha = 0.5) +
  geom_vline(xintercept = 15, lty = "dashed") +
  geom_hline(yintercept = "Muricidae", lty = "dashed") +
  labs(y = "Family", x = "Rank", colour = "Rank Type")

# just the families I want to keep
kelp_fam_list_cut <- kelp_fam_big_list %>%
  arrange(rank_all) %>%
  head(15) %>%
  select(family) 

# OK and now I want to make a final df where the top 15 families are named, and everything else is "other"

# make this final df for the df without 0's
data_fam_no0 <- data_fam_no0_a %>%
  filter(family %in% kelp_fam_list_cut$family) %>%
  # now make all the other families "other"
  rbind(data_fam_no0_a %>%
          filter(!family %in% kelp_fam_list_cut$family) %>%
          mutate(family = "other"))

# Family stats ----
mod_fam <- glmmTMB(in_out_avg ~ kelp_sp + kelp_bio_scale + tide_scale + 
                        family*abund_fam_scale + depth_scale ,
                      family = 'gaussian',
                      data = data_fam_no0)
summary(mod_fam)
plot(DHARMa::simulateResiduals(mod_fam)) 

df <- emtrends(mod_fam, pairwise ~ family, var = "abund_fam_scale")$emtrends %>%
  as.data.frame()

# look at those slopes
ggplot(df, aes(x = abund_fam_scale.trend, y = reorder(family, abund_fam_scale.trend), xmin = lower.CL, xmax = upper.CL)) +
  geom_point(size = 2.7) +
  geom_errorbar(width = 0, linewidth = 0.5) +
  geom_vline(xintercept=0, color="black", linetype="dashed")

# For fish: Hexagrammidae, Sebastidae, Embiotocidae, Cottidae
# For inverts: Echinasteridae, Asteropseidae, Muricidae, Asteriidae,  Stichopodidae by size of slope

# plot each fish family model -----
# since I'm working with each family seperately I want to change density back to abundance
data_fam_no0s <- data_fam_no0 %>%
  mutate(total_fam = case_when(method == 1 ~ total_fam*500,
                               method == 2 ~ total_fam*100),
         abund_fam_scale = c(scale(total_fam)))

# should I do a model for each species?
v_fam <- v_fun(data_fam_no0s, "Hexagrammidae")
predict_green <- fam_fun2(data_fam_no0s, "Hexagrammidae", diagnose = FALSE)
# total_fam  9.4316940  3.9924145   2.362  0.01816 * 

v_fam <- v_fun(data_fam_no0s, "Sebastidae")
predict_rock <- fam_fun2(data_fam_no0s, "Sebastidae", diagnose = FALSE)

v_fam <- v_fun(data_fam_no0s, "Embiotocidae")
predict_per <- fam_fun2(data_fam_no0s, "Embiotocidae", diagnose = FALSE)

v_fam <- v_fun(data_fam_no0s, "Cottidae")
predict_scul <- fam_fun2(data_fam_no0s, "Cottidae", diagnose = FALSE)

predict_fish <- rbind(predict_green, predict_rock, predict_per, predict_scul)%>%
  mutate(family = factor(family, levels = 
                           c("Hexagrammidae", "Sebastidae", "Embiotocidae", "Cottidae")))

# make a df of just those species
fish_fam_data <- data_fam_no0s %>%
  filter(family == "Hexagrammidae" |
           family == "Sebastidae" |
           family == "Cottidae" |
           family == "Embiotocidae") %>%
  mutate(family = factor(family, levels = 
                           c("Hexagrammidae", "Sebastidae", "Embiotocidae", "Cottidae")),
         slope = case_when(family == "Hexagrammidae" ~ "slope = 0.009, p = 0.0002",
                           family == "Sebastidae" ~ "slope = 0.0009, p = 0.10",
                           family == "Cottidae" ~ "slope = 0.003, p = 0.01",
                           family == "Embiotocidae" ~ "slope = 0.0007, p = 0.0001"))

# plot these curves for the fish 
pal <- viridis::viridis(10)
pal1 <- pal[5]

ggplot() + 
  geom_point(data = fish_fam_data, 
             aes(x = total_fam, y = in_out_avg), colour = pal1,
             alpha = 0.8) +
  labs(y = expression(paste("Delta ammonium ", (mu*M))), 
       x = expression(paste("Abundance"))) +
  facet_wrap(~family, scales = 'free_x') +
  geom_line(data = predict_fish,
            aes(x = total_fam, y = predicted), colour = pal1,
            linewidth = 1) +
  geom_ribbon(data = predict_fish,
              aes(x = total_fam, y = predicted, 
                  ymin = conf.low, ymax = conf.high), fill = pal1,
              alpha = 0.15) +
  theme_white() +
  theme(strip.background = element_rect(fill = "grey", color = "grey")) +
  geom_text(
    data = fish_fam_data %>% select(family, slope) %>% unique(),
    mapping = aes(x = Inf, y = Inf, label = slope),
    hjust   = 1.1,
    vjust   = 1.5,
    size = 9)

 ggsave("Output/Figures/fish_families_kelp.png", device = "png", height = 9, width = 12, dpi = 400)


# Next do inverts! -----

# Model selection checks------
# biomass and abundance, shannon vs simpson checks
mod_abund <- glmmTMB(in_minus_out ~ kelp_sp +
                       kelp_bio_scale*tide_scale +
                       kelp_bio_scale*abundance_scale +
                       abundance_scale*tide_scale +
                       shannon_scale + depth_scale + 
                       (1|site_code),
                     family = 'gaussian',
                     data = data) 

mod_simp <- glmmTMB(in_minus_out ~ kelp_sp + 
                      kelp_bio_scale*tide_scale +
                      kelp_bio_scale*weight_sum_scale +
                      weight_sum_scale*tide_scale +
                      simpson_scale + depth_scale + 
                      (1|site_code),
                      family = 'gaussian',
                      data = data) 

# Abundance simpson's is the best AIC mod
mod_abund_simp <- glmmTMB(in_minus_out ~ kelp_sp + 
                            kelp_bio_scale*tide_scale +
                            kelp_bio_scale*abundance_scale +
                            abundance_scale*tide_scale +
                            simpson_scale + depth_scale + 
                            (1|site_code),
                          family = 'gaussian',
                          data = data) 

AIC_tab_kelp <- AIC(mod_in_out, mod_abund, mod_simp, mod_abund_simp) %>%
  rownames_to_column() %>%
  mutate(best = min(AIC),
         delta = round((AIC - best), digits = 2),
         likelihood = exp( -0.5*delta),
         sum = sum(likelihood),
         AICw = round((likelihood/sum), digits = 2),
         AIC = round(AIC, digits = 2)) %>%
  select(rowname, df, AIC, delta, AICw) 


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

ggplot(data) +
  geom_point(aes(BiomassM, nh4_inside, shape = kelp_sp, colour = "inside")) +
  geom_point(aes(BiomassM, nh4_outside, shape = kelp_sp, colour = "outside")) +
  geom_point(aes(BiomassM, in_minus_out, shape = kelp_sp, colour = "delta"), size = 3) +
  geom_hline(yintercept= 0, linetype = "dashed", linewidth = 0.5) 
  
  
  


kcca_table <- pee %>%
  select(site_code, date) %>%
  unique() %>%
  filter(site_code != "KCCA13")



# Graveyard ------



data2 <- data %>%
  select(site, kelp_sp, Composition)
