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

#Should I drop the no kelp control from the model and just say we found no diff in ammonium on the sand? Basically using the no kelp control as a methods control but not using it in the model

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
         exp_abund = exp(abundance))

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
  # interaction between kelp:tide_exchange
  # interaction between animal biomass:kelp
  # maybe a triple interaction between kelp:animal_biomass:tide??
# Only kelp:tide interaction! 

# Let's plot out some of these interactions to make sure everything makes sense
ggplot(data, aes(kelp_bio_scale, in_out_avg, colour = kelp_sp)) +
  geom_point(aes(pch = kelp_sp))+
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, lty = "dashed")
# looks like a positive delta vs kelp biomass

ggplot(data, aes(BiomassM, in_minus_out)) +
  geom_point()+
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, lty = "dashed")

ggplot(data, aes(bio_tran_scale, in_minus_out, colour = kelp_sp)) +
  geom_point(aes(pch = kelp_sp))+
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, lty = "dashed")
# on the transect scale it also looks like there's a positive delta vs kelp trend

ggplot(data, aes(weight_sum_scale, in_out_avg, colour = kelp_sp)) +
  geom_point(aes(pch = kelp_sp))+
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, lty = "dashed")
# looks like a positive delta vs animal biomass/abundance

ggplot(data, aes(depth_scale, in_out_avg, colour = kelp_sp)) +
  geom_point(aes(pch = kelp_sp))+
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, lty = "dashed")
# maybe no depth relationship, but it looks like it's positive for nereo

ggplot(data, aes(tide_scale, in_out_avg, colour = kelp_sp)) +
  geom_point(aes(pch = kelp_sp))+
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, lty = "dashed")
# maybe negative delta vs tide relationship 

# looks like there might be a kelp:tide interaction but I'm worried about the kelp species
  # there's only one nereo site at flood tide, and only three macro sites at flood
  # there's two low delta nh4 low kelp sites and two high delta nh4 high kelp sites
# there's only three nereo sites
  # so I don't think I want to mess around with kelp_sp interactions! Not enough data for nereo or no kelp sites
# the no kelp sites happened to be the deepest sites too

# Build full model with gaussian distribution
# This is the transect level model
# One row per transect, 3 per site
mod_tran <- glmmTMB(in_minus_out ~ kelp_sp + 
                      kelp_bio_scale*tide_scale*abundance_scale + 
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
# I now know abundance is actually better than weight, the resids look better with abundance
# no interactions bc I don't think I have the data to do that
mod_in_out <- glmmTMB(in_minus_out ~ kelp_sp + kelp_bio_scale + tide_scale + 
                        abundance_scale + shannon_scale + depth_scale + 
                        (1|site_code),
                      family = 'gaussian',
                      data = data) # when I remove no kelp sites the coeffs hardly change

plot(simulateResiduals(mod_in_out)) # looks ok but not perfect
summary(mod_in_out)

# what if i just add the abundance:kelp interaction
mod_int <- glmmTMB(in_minus_out ~ kelp_sp + kelp_bio_scale*abundance_scale +
                     tide_scale + shannon_scale + depth_scale + 
                     (1|site_code),
                      family = 'gaussian',
                      data = data) # when I remove no kelp sites the coeffs hardly change

plot(simulateResiduals(mod_int)) # looks better! no red
summary(mod_int)
# seems like the other 2 2-way interactions only come out when the kelp:tide interaction is also involved. if the kelp:tide interaction isn't valid to include, maybe none of them are???


# what if I add kelp sp interaction
# I can't do a triple kelp_sp*kelp_bio_scale*tide_scale interaction bc no nereo flood
# I want to know if the tide:kelp interaction is the same for nereo and macro
    # there is no tide:nereo interaction bc nereo only measured at slack tide
    # and there's no kelp:tide interaction for no kelp bc there's no kelp biomass
# and i want to know if the kelp biomass slope is same for macro and nereo
# looks like it is???? but there's not enough data points to run the model with just nereo

# biomass and abundance, shannon vs simpson checks
mod_weight <- glmmTMB(in_minus_out ~ kelp_sp + kelp_bio_scale + tide_scale + 
                       weight_sum_scale + shannon_scale + depth_scale + 
                       (1|site_code),
                     family = 'gaussian',
                     data = data)

mod_simp <- glmmTMB(in_minus_out ~ kelp_sp + kelp_bio_scale + tide_scale + 
                      abundance_scale + simpson_scale + depth_scale + 
                      (1|site_code),
                    family = 'gaussian',
                    data = data)

# Abundance simpson's is the best AIC mod
mod_weight_simp <- glmmTMB(in_minus_out ~ kelp_sp + kelp_bio_scale + tide_scale + 
                            weight_sum_scale + simpson_scale + depth_scale + 
                            (1|site_code),
                          family = 'gaussian',
                          data = data)

mod_tran <- glmmTMB(in_minus_out ~ kelp_sp + bio_tran_scale + tide_scale +
                      abundance_scale + shannon_scale + depth_scale +
                      (1|site_code),
                      family = 'gaussian',
                      data = data)

mod_forest <- glmmTMB(in_minus_out ~ kelp_sp + forest_bio_scale + tide_scale +
                        abundance_scale + shannon_scale + depth_scale +
                        (1|site_code),
                      family = 'gaussian',
                      data = data)


AIC_tab_kelp <- AIC(mod_in_out, mod_weight, mod_simp, mod_weight_simp, mod_forest, mod_tran, mod_best, mod_int) %>%
  rownames_to_column() %>%
  mutate(best = min(AIC),
         delta = AIC - best,
         likelihood = exp( -0.5*delta),
         sum = sum(likelihood),
         AICw = likelihood/sum) %>%
  select(rowname, df, AIC, delta, AICw)
# ooof the model with the three t-way interactions are super preferred 

# let's look at those interactions some more:

# Kelp:Tide
visreg(mod_best, "kelp_bio_scale", by = "tide_scale", overlay=TRUE)
# the range for low and high tide exchange rates is pretty small
visreg(mod_best, "tide_scale", by = "kelp_bio_scale", overlay=TRUE)
# looking at it the other way, there's a good range of tides across the med kelp density, but not at low or high kelp
# I don't think this tide:kelp interaction is valid!

# Abundance:Kelp
visreg(mod_best, "abundance_scale", by = "kelp_bio_scale", overlay=TRUE)
visreg(mod_best, "kelp_bio_scale", by = "abundance_scale", overlay=TRUE)
# looks like an ok spread, 2 groups only have 3 sites, tho
# this basically looks like: at high kelp or high animals, you get a large delta nh4. But if one of the variables is med/low, increasing the other will increase delta
# I think this is my "asympote"

# Abundance:tide
visreg(mod_int, "abundance_scale", by = "tide_scale", overlay=TRUE)

# OK so tide isn't actually a reasonable thing to include!!!!! the spead of kelp biomass across tides is pittiful, it's trying to estimate slopes for two clumps of points and I don't trust that those slopes are actually different

# what does the visreg look like with no interactions
visreg(mod_in_out)


# run the model without an intercept
mod_in_out2 <- glmmTMB(in_minus_out ~ -1 + kelp_sp + kelp_bio_scale + tide_scale + 
                         abundance_scale + simpson_scale + depth_scale +
                         (1|site_code),
                      family = 'gaussian',
                      data = data)
summary(mod_in_out2)

# don't scale variables for estimates in normal units
mod_in_out_c <- glmmTMB(in_minus_out ~ - 1 + kelp_sp + kelp_bio_center + tide_center + 
                          abundance_center + simpson_center + depth_center +
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
pal12 <- viridis::viridis(8)
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
                           levels = c("kelp_spmacro", "kelp_spnereo", "kelp_spnone", "kelp_bio_scale", "abundance_scale", "depth_scale", "simpson_scale", "tide_scale"),
                           labels = c("Macro", "Nereo", "None", "Kelp biomass", "Animal abund", "Depth", "Biodiversity", "Tide")),
         se = (upper_CI - estimate)/1.96
  )

# Use function to plot coefficients
kelp_coeff_plot <- coeff_plot(coeff_df = df,
                              pal = pal12) +
  place_label("(a)")

#ggsave("Output/Pres_figs/Fig4a.png", device = "png", height = 8, width = 12, dpi = 400)

# Plot kelp species predictions ----
predict_sp <- ggpredict(mod_in_out, terms = "kelp_sp") %>% 
  dplyr::rename(kelp_sp = x,
                in_minus_out = predicted
         ) 

sp_labs <- c(macro = "Macro", nereo = "Nereo", none = "None")

kelp_sp_plot <- 
  dot_whisker(sum_data = predict_sp, 
              all_data = data,
              x_var = kelp_sp,
              y_var = in_minus_out,
              pch_var = kelp_sp,
              labels = sp_labs,
              pal = c(pal12[3], pal12[2],pal12[1])) +
    geom_hline(yintercept = 0, lty = "dashed") +
    place_label("(b)")

# Plot kelp_biomass predictions ----
predict_kelp_bio <- ggpredict(mod_in_out, terms = "kelp_bio_scale") %>% 
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
predict_abund <- ggpredict(mod_in_out, terms = "abundance_scale") %>% 
  dplyr::rename(abundance_scale = x)

abund_pred_plot <- plot_pred(raw_data = data,
                             predict_data = predict_abund, 
                             plot_type = "new_kelp",
                             x_var = abundance_scale, y_var = in_minus_out,
                             pch_var = kelp_sp,
                             x_axis_lab = "Abundance",
                             pal = pal12[5],
                             lty_var = pal12[5]) +  
  place_label("(d)")

# Plot depth? ---- 
predict_depth <- ggpredict(mod_in_out, terms = "depth_scale") %>% 
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

# Put them allll together????

plots <- kelp_sp_plot + kelp_pred_plot + abund_pred_plot + depth_pred_plot

kelp_coeff_plot + plots
  
#ggsave("Output/Figures/Fig4.png", device = "png", height = 13.5, width = 24, dpi = 400)
  
  




# OLD Model predictions against raw data plot ----
# make new mod with untransformed vars

# create range vector
tide_means_kelp <- data %>%
  group_by(tide_cat) %>%
  summarise(tide = mean(tide_scale))

v2 <- c(as.numeric(tide_means_kelp[2,2]), as.numeric(tide_means_kelp[1,2]))

# now make predictions
predict_kelp <- ggpredict(mod_abund, terms = c("kelp_bio_scale", "tide_scale [v2]")) %>% 
  mutate(kelp_bio_scale = x,
         tide_cat = factor(as.factor(ifelse(group == as.character(tide_means_kelp[1,2]), "Slack", "Flood")),
         levels = c("Ebb", "Slack", "Flood"))
         ) %>%
  filter(tide_cat != "Flood" | kelp_bio_scale < 0.21) %>%
  filter(tide_cat != "Flood" | kelp_bio_scale > -0.6)


# now plot these predictions
kelp_pred_plot <- plot_pred(raw_data = (data %>%
                        mutate(tide = ifelse(avg_exchange_rate < 0, "Slack", "Flood"))),
          predict_data = predict_kelp, 
          plot_type = "kelp",
          x_var = kelp_bio_scale, y_var = in_minus_out, 
          lty_var = tide_cat,
          size_var = weight_sum, 
          pal = pal2) +
#  ylim(c(-0.5, 0.8))
  place_label("(b)")

#ggsave("Output/Pres_figs/Fig4b.png", device = "png", height = 9, width = 12, dpi = 400)


# Figure 4 for pub with white -----
kelp_coeff_plot + kelp_pred_plot

#ggsave("Output/Pub_figs/Fig4b.png", device = "png", height = 9, width = 16, dpi = 400)


# Plot kelp species vs kelp biomass vs nh4 ------
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
max_min_sp <- data_kelp %>%
  group_by(kelp_sp) %>%
  summarise(min = min(kelp_bio_scale) - 0.1,
            max = max(kelp_bio_scale) + 0.1)

# tide_stuff
tide_means_kelp <- data_kelp %>%
  group_by(tide_cat) %>%
  summarise(tide = mean(tide_scale))

v2 <- c(as.numeric(tide_means_kelp[2,2]), as.numeric(tide_means_kelp[1,2]))

# now make predictions
predict_kelp_sp <- ggpredict(mod_sp, terms = c("kelp_bio_scale", "kelp_sp")) %>% 
  mutate(kelp_bio_scale = x,
         kelp_sp = as.factor(group)) %>%
  left_join(max_min_sp, by = "kelp_sp") %>%
  rowwise %>%
  filter(between(kelp_bio_scale, min, max)) 
 #filter(tide_cat != "Flood" | kelp_bio_scale < 0.21) %>%
 #filter(tide_cat != "Flood" | kelp_bio_scale > -0.6)

# plot predictions
ggplot() + 
  geom_point(data = data_kelp, 
             aes(x = kelp_bio_scale, y = in_minus_out, 
                 colour = kelp_sp, fill = kelp_sp,
                 pch = kelp_sp), 
             alpha = 0.8) +
  geom_smooth(data = data_kelp, 
              aes(x = kelp_bio_scale, y = in_minus_out, 
                  colour = kelp_sp, fill = kelp_sp),
              method = lm)
  
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
               fun.data = "mean_cl_boot") 

  labs(colour = "Tide", fill = "Tide", lty = "Tide") +
  theme +
  scale_colour_manual(values = (pal)) +
  scale_fill_manual(values = (pal)) +

    geom_hline(yintercept= 0, linetype = "dashed", color = features, linewidth = 0.5)+
    guides(lty = guide_legend(override.aes = list(linewidth = 0.5)),
           size = guide_legend(override.aes = list(colour = features)),
           colour = guide_legend(override.aes = list(size = 2)))  +
    scale_x_continuous(breaks = c(-1.17905227, -0.1, 1, 2.05),
                       labels = c("0", "0.6", "1.2", "1.8")) +
    labs(y = expression(paste(Delta, " Ammonium ", (mu*M))), 
         x = expression(paste("Kelp biomass (kg/m"^2,")")),
         size = "Animals (kg)")

# now plot these predictions
plot_pred(raw_data = data_kelp,
          predict_data = predict_kelp_sp, 
          plot_type = "kelp",
          x_var = kelp_bio_scale, y_var = in_minus_out, 
          lty_var = kelp_sp,
          size_var = weight_sum, 
          pal = pal2)


# so basically for Macrocystis, there's a positive interaction between kelp biomass and tide exchange. But that's only driven by two high in - out sites and one low in - out site at flood tide.... so is that interaction real

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
