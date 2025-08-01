# Script to look into the variation in pee inside vs outside kelp forests
# March 3, 2023, last updated Dec 2024
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
library(vegan)
library(viridis)
library(emmeans)
library(ggtext)

# set theme and load functions
theme_set(theme_bw())
source("Code/Functions.R")
source("Code/Fig_theme_functions.R")


# NH4+ data -----
# kelp pee inside vs outside data from Kelp_pee_nh4_calc.R
pee <- read_csv("Output/Output_data/kelp_pee.csv") %>%
  as.data.frame() %>%
  filter(site_code != "KCCA13") %>% # remove second beach south outlier
  left_join(
    read_csv("Data/Team_kelp/Output_data/site_names.csv") %>%
               rename(site_name = SiteName), by = "site_code") %>%
  select(-site) %>%
  mutate(x_change = nh4_inside/nh4_outside)
      

# just save the sites we surveyed
kcca_surveys <- pee %>%
  transmute(Date = date,
            site_code = site_code) %>%
  unique() 

# Kelp data ----

# transect level kelp biomass + density data from Claire
kelp <- read_csv("Data/Team_kelp/Output_data/kelp_metrics2024.csv") %>%
  as.data.frame() %>%
  # rename vars now that all the data is linked up
  rename(kelp_den = Kelp,
         site_name = SiteName) %>%
  # only keep sites I measured pee 
  filter(site_code %in% kcca_surveys$site_code) %>%
  mutate(sample = case_when(Transect == 0 ~ 1, 
                            Transect == 5 ~ 2,
                            Transect == 10 ~ 3),
         kelp_sp = as.factor(case_when(
           site_code == "KCCA19" ~ "none", # no kelp site that had 2 nereo
           kelp_den == 0 ~ "none", # no kelp = none
           Macro_5m2 == 0 ~ "nereo", # no macro = nereo
           Nereo_5m2 == 0 ~ "macro", # no nereo = macro
           TRUE ~ as.character("macro")))) %>% # everything else = macro
  replace(is.na(.), 0)  %>%
  group_by(site_code) %>%
  # site level summaries
  mutate(BiomassM = ifelse(kelp_sp == "none", 0, mean(biomass_trans_mean)), # used to be BiomassTkg
         Area_m2 = ifelse(kelp_sp == "none", 0, Area_m2),
         DensityM = mean(kelp_den)) %>%
  ungroup()  %>%
  filter(Transect != "15") # this means the fourth transect contributes to the kelp density metrics but is cut bc there's no pee sample on that tran

# Less dangerous bay was the "real" no kelp control


# Transect level variables:
# biomass_trans_mean is kelp density x average biomass for each transect

# Site level variables
# BiomassM is the average biomass/m2 at each site across all 4 transects
# DensityM is the average density at each site across all 4 transects


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
  # only keep sites I measured pee 
  filter(site_code %in% kcca_surveys$site_code) %>% 
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


# extract just one row per survey to join with the pee data and tide data
# Only keep the surveys I have pee samples for!
kcca_survey_info <- kcca_surveys %>%
  left_join(kelp_rls %>%
              select(site_code, site_name, Date, date_time_survey) %>%
              unique() %>%
              mutate(survey_id = 1:n()), by = c("Date", "site_code")) %>%
  left_join((kelp %>%
               select(c(site_code, Time_start)) %>%
               unique()), by = "site_code") %>%
  mutate(date_time_kelp = ymd_hms(paste(Date, Time_start)))

# pivot back to wide for biodiversity
kelp_rls_wider <- kelp_rls %>% 
  left_join(kcca_survey_info, by = c("site_code", "date_time_survey")) %>%
  dplyr::select(survey_id, species_name, survey_den) %>% 
  group_by(survey_id, species_name) %>%
  summarise(survey_den = sum(survey_den)) %>%
  ungroup() %>%
  spread(key = species_name, value = survey_den) %>%
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
         percent_diff = 100*((nh4_in_avg-nh4_out_avg)/nh4_out_avg),
         depth_avg = mean(depth_m)) %>%
  ungroup() %>%
  # join kelp data with both transect + site avgs
  left_join(kelp, by = c("site_code", "site_name", "sample")) %>%
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
                        abundance = sum(survey_den)) %>%
              ungroup()) %>%
  # diversity indexes
  left_join(kelp_rls_wide, by = "survey_id") %>%
  # join tide exchange data
  left_join(tide_exchange_kcca, by = "survey_id") %>%
  # scale factors using function
  scale_vars() %>%
  mutate(exp_kelp = exp(BiomassM),
         exp_abund = exp(abundance),
         kelp_cat = factor(as.factor(
           case_when(
             kelp_bio_scale < -0.5 ~ "Low",
             kelp_bio_scale < 1 ~ "Mid",
             kelp_bio_scale > 1 ~ "High")),
           levels = c("Low", "Mid", "High"))) %>%
  as.data.frame()


# reduced data to export a file for mapping

# save csv for mapping 
kelp_map_data <- data %>%
  select(site_name, site_code, in_out_avg, nh4_in_avg, nh4_out_avg, nh4_avg) %>%
  left_join(read_csv("Data/Team_kelp/RLS_KCCA_2022.csv") %>%
              transmute(site_code = `Site No.`,
                        latitude = Latitude,
                        longitude = Longitude) %>%
              unique(), by = "site_code") %>%
  unique()
# write_csv(kelp_map_data, "Output/Output_data/kelp_map_data.csv")



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

# Abiotic to control for
# depth_avg (depth_scale)
# avg_exchange_rate (tide_scale)

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

mod_best <- glmmTMB(in_minus_out ~ kelp_sp + kelp_bio_scale + abundance_scale + depth_scale + tide_scale + 
                      abundance_scale:kelp_bio_scale + abundance_scale:tide_scale + kelp_bio_scale:tide_scale + (1|site_code), 
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

# what about looking at % change
# This tells me what x increase of nh4+ we had in each forest
mod_xchange <- glmmTMB(x_change ~ -1 + kelp_sp + 
                        kelp_bio_scale*tide_scale +
                        kelp_bio_scale*weight_sum_scale +
                        weight_sum_scale*tide_scale +
                        shannon_scale + depth_scale + 
                        (1|site_code),
                      family = 'gaussian',
                      data = data) 
# when I remove no kelp sites the coeffs hardly change

plot(simulateResiduals(mod_xchange)) # looks ok!
summary(mod_xchange)

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

# centered variables for interpretation in text
mod_c <- glmmTMB(in_minus_out ~ -1 + kelp_sp + 
                   kelp_bio_center*tide_center +
                   kelp_bio_center*weight_sum_center +
                   weight_sum_center*tide_center +
                   shannon_center + depth_center +
                   (1|site_code),
                 family = 'gaussian',
                 data = data) 
summary(mod_c)

# untransformed variables for plotting
mod_in_out3 <- glmmTMB(in_minus_out ~ -1 + kelp_sp + 
                         BiomassM*avg_exchange_rate +
                         BiomassM*weight_sum +
                         weight_sum*avg_exchange_rate +
                         shannon_scale + depth_scale + 
                         (1|site_code),
                       family = 'gaussian',
                       data = data)

# what if I add kelp sp interaction
# I can't do a triple kelp_sp*kelp_bio_scale*tide_scale interaction bc no nereo flood
# I want to know if the tide:kelp interaction is the same for nereo and macro
# there is no tide:nereo interaction bc nereo only measured at slack tide
# and there's no kelp:tide interaction for no kelp bc there's no kelp biomass
# basically the model is guessing at what's happening with nereo based on what's happening with macro and that's fine


# Step 4: Check for collinearity of predictors

# car can't handle random effects so make a simplified mod
car::vif(lm(in_minus_out ~ kelp_sp + kelp_bio_scale + tide_scale + weight_sum_scale + shannon_scale +  depth_scale, data = data))
# Yep this looks fine!


# Notes from stats beers!
# Maybe try to PCA????? see if there's clustering?

# is there a kelp ~ animal relationship?
mod_kelp <- glmmTMB(abundance ~ BiomassM,
                      family = 'gaussian',
                      data = data %>% select(BiomassM, abundance, depth_avg) %>% unique) 

plot(DHARMa::simulateResiduals(mod_kelp))
summary(mod_kelp)

# Graphing ------

# remake with 30 colour palette
pal30 <- viridis::viridis(30)
pie(rep(1, 30), col = pal30)

pal8c <- viridis::viridis(30)[22:30] # coeff plot
pal_spc <- c(pal30[19], pal30[20], pal30[21]) # kelp sp colours
pal3c <- c(pal30[7], pal30[12], pal30[17]) # kelp biomass colours 
pal2c <- c(pal30[1], pal30[9]) # tide colours


# Fig 3a: model coefficients ----
# Save model coefficients 
# use the no intercept model for plotting
df <- confint(mod_in_out, level = 0.95, method = c("wald"), component = c("all", "cond", "zi", "other"), estimate = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(variable = rowname,
         lower_CI = `2.5 %`,
         upper_CI = `97.5 %`,
         estimate = Estimate) %>%
  head(- 1) %>%
  tail(-3) %>%
  mutate(variable = factor(as.factor(variable), 
                           labels = c("Depth", "Kelp biomass", "Kelp:tide", "Kelp:animals", "Biodiversity", "Tide", "Tide:animals", "Animal biomass")),
         se = (upper_CI - estimate)/1.96
  )%>%
  arrange(desc(estimate)) %>%
  mutate(variable = factor(variable, unique(variable)))



# plot just cont vars
kelp_coeff_plot <- coeff_plot(coeff_df = df,
                              pal = rev(pal8c)) +
  place_label("(a)") # this is a function I create in the Functions.R file



# Fig 3b: Kelp species predictions ----
# make predictions for macro and nereo at mean kelp bio, make prediction for the no kelp at the mean of the no kelp sites
kelp_means <- data %>%
  mutate(kelp = ifelse(kelp_sp == "none", "none", "kelp")) %>%
  group_by(kelp) %>%
  summarise(mean = mean(BiomassM))

v <- c(kelp_means$mean[1], kelp_means$mean[2])

predict_sp <- ggpredict(mod_in_out3, terms = c("BiomassM[v]", "kelp_sp")) %>% 
  as.data.frame() %>%
  dplyr::rename(BiomassM = x,
                in_minus_out = predicted,
                kelp_sp = group) %>%
  mutate(kelp = ifelse(kelp_sp == "none", "none", "kelp"),
         kelp_sp = factor(kelp_sp, levels = c("nereo", "macro", "none"))) %>%
  left_join(kelp_means, by = "kelp") %>%
  filter(mean ==  BiomassM)

data2 <- data %>%
  mutate(kelp_sp = factor(kelp_sp, levels = c("nereo", "macro", "none")))

sp_labs <- c(macro = "Macro", nereo = "Nereo", none = "No kelp")

set.seed(234444)
kelp_sp_plot <- 
  dot_whisker(sum_data = predict_sp, 
              all_data = data2,
              x_var = kelp_sp,
              y_var = in_minus_out,
              #pch_var = kelp_sp,
              labels = sp_labs,
              pal = pal_spc) +
  labs(y = expression(paste(Delta, " Ammonium ", (mu*M))),
       x = "Kelp species") +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.25) +
  place_label("(b)") +
  ylim(c(-0.79, 1.05))


# Fig 3c: Kelp biomass x tide interaction ------
visreg(mod_in_out, "kelp_bio_scale", by = "tide_scale", overlay=TRUE,
       xlab = "Kelp biomass (kg/m2)")

# create range vector
tide_means_kelp <- data %>%
  group_by(tide_cat) %>%
  summarise(tide = mean(avg_exchange_rate))

v2 <- c(as.numeric(tide_means_kelp[2,2]), as.numeric(tide_means_kelp[1,2]))

v_kelp <- seq(from = min(data$BiomassM), 
              to = max(data$BiomassM),
              length.out = 100)

# now make predictions
predict_kelp <- ggpredict(mod_in_out3, terms = c("BiomassM [v_kelp]", "avg_exchange_rate [v2]")) %>% 
  as.data.frame() %>%
  mutate(BiomassM = x,
         tide_cat = factor(as.factor(ifelse(
           group == as.character(tide_means_kelp[1,2]), "Slack", "Flood")),
           levels = c("Ebb", "Slack", "Flood"))) %>%
  filter(tide_cat != "Flood" | BiomassM < 3.8) %>%
  filter(tide_cat != "Flood" | BiomassM > 1.8)


# now plot these predictions
kelp_tide_int_plot <- 
  plot_pred(raw_data = (data %>% mutate(
    tide = ifelse(avg_exchange_rate < 0, "Slack", "Flood"))),
    predict_data = predict_kelp, 
    plot_type = "new_kelp",
    x_var = BiomassM, y_var = in_minus_out, 
    lty_var = tide_cat,
    pch_var = tide_cat,
    pal = rev(pal2c),
    x_axis_lab = expression(paste("Kelp biomass (kg/m"^2,")"))) +
  theme(legend.position = "none") +
  place_label("(c)") +
  ylim(c(-0.79, 1.05)) 



# Fig 3d: Kelp biomass x abundance interaction -----
visreg(mod_in_out2, "weight_sum_scale", by = "kelp_bio_scale", 
       overlay=TRUE,
       xlab = "Animal biomass (kg/m2)") # looks better
# create range vector based on visreg numbers

hist(data$kelp_bio_scale)

# get mean kelp bio for each kelp category, and find min and max animal biomass to limit the lines to the available data
kelp_sp_means <- data %>%
  group_by(kelp_cat) %>%
  summarise(mean = mean(BiomassM),
            min = min(weight_sum),
            max = max(weight_sum))

# mean kelp bio for low, mid, and high help biomass
v6 <- c(kelp_sp_means$mean[1], kelp_sp_means$mean[2], kelp_sp_means$mean[3])


# now make predictions
predict_abund_kelp <- ggpredict(mod_in_out3, terms = c("weight_sum", "BiomassM [v6]")) %>% 
  as.data.frame() %>%
  mutate(weight_sum = x,
         kelp_cat = factor(as.factor(case_when(
           group == kelp_sp_means$mean[1] ~ "Low",
           group == kelp_sp_means$mean[2] ~ "Mid",
           group == kelp_sp_means$mean[3] ~ "High")),
           levels = c("High", "Mid", "Low"))) %>%
  left_join(kelp_sp_means, by = "kelp_cat") %>%
  rowwise() %>%
  filter(between(weight_sum, min, max))


# now plot these predictions
abund_kelp_int_plot <- 
  plot_pred(raw_data = data %>% 
              mutate(kelp_cat = factor(kelp_cat, 
                                       levels = c("High", "Mid", "Low"))),
            predict_data = predict_abund_kelp, 
            plot_type = "new_kelp",
            x_var = weight_sum, y_var = in_minus_out, 
            lty_var = kelp_cat,
            pch_var = kelp_cat,
            x_axis_lab = expression(paste("Animal biomass (kg/m"^2,")")),
            pal = rev(pal3c)) +
  guides(size = "none") +
  ylim(c(-0.79, 1.05)) +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank()) +
  labs(colour = "Kelp biomass", fill = "Kelp biomass", pch = "Kelp biomass") +
  place_label("(d)") 


# Fig 3e: Abundance x tide interaction ------
visreg(mod_in_out, "weight_sum_scale", by = "tide_scale", overlay=TRUE,
       xlab = "Animal biomass (kg/m2)")
visreg(mod_in_out, "tide_scale", by = "weight_sum_scale", overlay=TRUE)

# create range vector
tide_means_kelp <- data %>%
  group_by(tide_cat) %>%
  summarise(tide = mean(avg_exchange_rate))

v2 <- c(as.numeric(tide_means_kelp[2,2]), as.numeric(tide_means_kelp[1,2]))

# range of animal vs tide
d2 <- data %>%
  group_by(tide_cat) %>%
  summarise(min = min(weight_sum) - 0.025,
            max = max(weight_sum) + 0.01)

# now make predictions
predict_abund_tide <- ggpredict(mod_in_out3, terms = c("weight_sum", "avg_exchange_rate [v2]")) %>% 
  as.data.frame() %>%
  mutate(weight_sum = x,
         tide_cat = factor(as.factor(ifelse(
           group == as.character(tide_means_kelp[1,2]), "Slack", "Flood")),
           levels = c("Ebb", "Slack", "Flood"))) %>%
  left_join(d2, by = "tide_cat") %>%
  rowwise() %>%
  filter(between(weight_sum, min, max))

# now plot these predictions
abund_tide_int_plot <- 
  plot_pred(raw_data = (data %>% 
                          mutate(tide = ifelse(avg_exchange_rate < 0, "Slack", "Flood"))),
            predict_data = predict_abund_tide, 
            plot_type = "new_kelp",
            x_var = weight_sum, y_var = in_minus_out, 
            lty_var = tide_cat,
            pch_var = tide_cat,
            x_axis_lab = expression(paste("Animal biomass (kg/m"^2,")")),
            pal = rev(pal2c)) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())+
  place_label("(e)") +
  guides(size = "none") + 
  ylim(c(-0.79, 1.05))

# Fig 3 panels ----
# plot coeffs + interactions
squish <- theme(axis.title.y = element_text(margin = margin(r = -90, unit = "pt")))
# -200 for theme black

kelp_coeff_plot/ ((kelp_sp_plot + squish) +
                    abund_kelp_int_plot + 
                    (kelp_tide_int_plot + squish) +
                    abund_tide_int_plot ) & theme(legend.justification = "left")

# correct size
# ggsave("Output/Pub_figs/Fig3.png", device = "png", height = 7, width = 6.5, units = "in", dpi = 400)

# old size
# ggsave("Output/Pub_figs/Fig3.png", device = "png", height = 16, width = 16, dpi = 400)



# Family manipulations ------
data_fam_no0_a <- data %>%
  select(site_name, site_code, in_out_avg, kelp_bio_scale, kelp_sp, depth_scale, weight_sum_scale, abundance_scale, shannon_scale, tide_scale, tide_cat) %>%
  unique() %>%
  left_join((kelp_rls %>%
               mutate(family = as.factor(family))%>%
               group_by(site_code, method, phylum, family) %>% # if i want methods split up, add it back here
               summarise(fam_den = sum(survey_den),
                         weight_fam_sum_g = 1000*sum(weight_size_class_sum))),
            by = "site_code") 


# what are the top families?

# CHOOSE THE FAMILIES TO INCLUDE
# rank of families by total abundance (density)
kelp_fam_list_total <- kelp_rls %>%
  group_by(family) %>%
  summarise(sum = sum(survey_den)) %>%
  drop_na(family) %>%
  arrange(desc(sum)) %>%
  transmute(family = family, 
            sum_total = sum,
            rank_total = 1:43)

# but which families show up on the most transects?
kelp_fam_list_count <- kelp_rls %>%
  select(site_code, family) %>%
  unique() %>%
  count(family) %>%
  drop_na(family) %>%
  arrange(desc(n)) %>%
  transmute(family = family, 
            sum_count = n,
            rank_count = 1:43)

# rank of families by biomass
kelp_fam_list_bio <- kelp_rls %>%
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

# since I'm working with each family separately I want to change density back to abundance
data_fam_no0s <- data_fam_no0 %>%
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


# Family stats ----
mod_fam <- glmmTMB(in_out_avg ~ kelp_sp + kelp_bio_scale + tide_scale + 
                     family*weight_den_fam_scale + depth_scale ,
                   family = 'gaussian',
                   data = data_fam_no0s)
summary(mod_fam)
plot(DHARMa::simulateResiduals(mod_fam)) 

df <- emmeans::emtrends(mod_fam, pairwise ~ family, var = "weight_den_fam_scale")$emtrends %>%
  as.data.frame()

# look at those slopes
ggplot(df, aes(x = weight_den_fam_scale.trend, y = reorder(family, weight_den_fam_scale.trend), xmin = lower.CL, xmax = upper.CL)) +
  geom_point(size = 2.7) +
  geom_errorbar(width = 0, linewidth = 0.5) +
  geom_vline(xintercept=0, color="black", linetype="dashed")

# For fish: Hexagrammidae, Sebastidae, Embiotocidae, Cottidae
# For inverts: Echinasteridae, Asteropseidae, Muricidae, Asteriidae by size of slope

# FAMILY MODELS
# I want to run one model for each family
# list family names
fam_kelp_levels <- levels(data_fam_no0s$family)

# lapply fam_fun over each family name to create df of predictions for each fam mod
fam_kelp_predictions <- lapply(fam_kelp_levels, function(family_name) {
  # This model has untransformed weight var for ease of plotting
  fam_fun_kelp_combo(data_fam_no0s, family_name)  
})

# join up those predictions
kelp_fam_predicts <- bind_rows(fam_kelp_predictions)

# what are the top 6 inverts and fish?
top_fam_kelp_r2 <- kelp_fam_predicts %>%
  select(c(family, r2)) %>%
  unique() %>%
  arrange(desc(r2)) %>%
  head(6) 
# Echinasteridae, Gobiidae, Cottidae, Strongylocentrotidae, Hexagrammidae, Turbinidae top 6

# make df of just the tops
top_fam_kelp_predict <- kelp_fam_predicts %>%
  filter(family %in% top_fam_kelp_r2$family) %>%
  # optional reorder of families
  arrange(desc(r2)) %>%
  mutate(family = factor(family, unique(family)))

# Outputs for these families
# list top 6 family names
top_fam_kelp_levels <- levels(top_fam_kelp_predict$family)

# get model outputs for each family model
# this function uses the scaled weight variable, which is what I want for the 'real' model
lapply(top_fam_kelp_levels, function(family_name) {
  diagnose_kelp_fun(data_fam_no0s, family_name) 
})
# Turbinidae doesn't have the nicest DHARMA plots

# now filter full family df to just include those top 6 families
kelp_fam <- top_fam_kelp_r2 %>%
  left_join(data_fam_no0s, by = "family") %>%
  arrange(desc(r2)) %>%
  mutate(family = factor(family, unique(family)))

# what % do these top 6 families make up
all_fam <- data_fam_no0s %>%
  summarise(type = "all",
            total_den = sum(fam_den),
            total_weight = sum(weight_fam_sum_g)) %>%
  rbind(kelp_fam %>%
          summarise(type = "top",
                    total_den = sum(fam_den),
                    total_weight = sum(weight_fam_sum_g)))

# Fig S2.2: Family plots ----
pal_k <- viridis::viridis(10)
pal1 <- pal_k[4]

fam_plot <- ggplot() + 
  geom_point(data = kelp_fam, 
             aes(x = weight_fam_sum_g, y = in_out_avg), colour = pal1,
             alpha = 0.8, size = 3) +
  labs(y = expression(paste(Delta, " Ammonium ", (mu*M))), 
       x = expression(paste("Biomass (g/m"^2,")"))) +
  facet_wrap(~family, scales = 'free_x', ncol = 3) +
  geom_line(data = top_fam_kelp_predict,
            aes(x = weight_fam_sum_g, y = predicted), colour = pal1,
            linewidth = 1.5) +
  geom_ribbon(data = top_fam_kelp_predict,
              aes(x = weight_fam_sum_g, y = predicted, 
                  ymin = conf.low, ymax = conf.high), fill = pal1,
              alpha = 0.2) +
  theme_white() +
  theme(strip.background = element_rect(fill = "grey", color = "grey"),
        axis.text = element_text(size = 18, color = "black", lineheight = 0.9)) + 
  theme(panel.spacing = unit(1, "lines")) +
  geom_hline(yintercept= 0, linetype = "dashed", color = "black", linewidth = 0.5)


fam_plot

# ggsave("Output/Pub_figs/Supp2Fig2.png", device = "png", height = 9, width = 16, dpi = 400)

# Fig S1.2 shannon vs biomass plots -----
pal10 <- viridis::viridis(10)
pal_sp2 <- c(pal10[8], pal10[6], pal10[4])

# and there's a positive relationship between kelp biomass and shannon?
kelp_shannon <- ggplot(data, aes(BiomassM, shannon, colour = kelp_sp, pch = kelp_sp, fill = kelp_sp)) +
  geom_point(alpha = 0.8, size = 4) +
  geom_smooth(method = lm, linewidth = 2, alpha = 0.2) +
  theme_white() +
  ylim(c(1.2, 2.8)) +
  theme(legend.position = "none") +
  scale_colour_manual(labels = sp_labs, values = (pal_sp2)) +
  scale_fill_manual(labels = sp_labs, values = (pal_sp2)) +
  labs(x = expression(paste("Kelp biomass (kg/m"^2,")")), y = "Shannon diversity") +
  place_label("(a)")


# what about weight? positive weight relationship
animal_shannon <- ggplot(data, aes(weight_sum, shannon, colour = kelp_sp, pch = kelp_sp, fill = kelp_sp)) +
  geom_point(alpha = 0.8, size = 4) +
  geom_smooth(method = lm, linewidth = 2, alpha = 0.2) +
  theme_white() +
  ylim(c(1.2, 2.8)) +
  labs(x = expression(paste("Animal biomass (kg/m"^2,")")), y = "Shannon diversity",
       colour = "Kelp species", pch = "Kelp species", fill = "Kelp species") +
  scale_colour_manual(labels = sp_labs, values = (pal_sp2)) +
  scale_shape_discrete(labels = sp_labs) +
  scale_fill_manual(labels = sp_labs, values = (pal_sp2)) +
  guides(lty = guide_legend(override.aes = list(linewidth = 0.5)),
         size = guide_legend(override.aes = list(colour = "black")),
         colour = guide_legend(override.aes = list(size = 3, linewidth = 1)))+
  place_label("(b)")


kelp_shannon + animal_shannon

# ggsave("Output/Pub_figs/Supp1Fig2.png", device = "png", height = 9, width = 16, dpi = 400)


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

mod_abund_simp <- glmmTMB(in_minus_out ~ kelp_sp + 
                            kelp_bio_scale*tide_scale +
                            kelp_bio_scale*abundance_scale +
                            abundance_scale*tide_scale +
                            simpson_scale + depth_scale + 
                            (1|site_code),
                          family = 'gaussian',
                          data = data) 

mod_rich <- glmmTMB(in_minus_out ~ kelp_sp + 
                        kelp_bio_scale*tide_scale +
                        kelp_bio_scale*weight_sum_scale +
                        weight_sum_scale*tide_scale +
                        rich_scale + depth_scale + 
                        (1|site_code),
                      family = 'gaussian',
                      data = data) 

mod_rich_abund <- glmmTMB(in_minus_out ~ kelp_sp + 
                      kelp_bio_scale*tide_scale +
                      kelp_bio_scale*abundance_scale +
                      abundance_scale*tide_scale +
                      rich_scale + depth_scale + 
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


# Summary stats -----
sum_kelp_pee <- data %>%
  group_by(site_code) %>%
  reframe(kelp_sp = kelp_sp,
          nh4_in = nh4_in_avg,
          nh4_out_avg = nh4_out_avg,
          x_change = nh4_in_avg/nh4_out_avg) %>%
  unique() %>%
  arrange(desc(x_change))

x_ch <- sum_kelp_pee %>%
  group_by(kelp_sp) %>% 
  summarise(mean_x_change = mean(x_change))


# Graveyard ------
# old analyses and scraps of code

# Community stuff -----
# make wide for families
fam_kelp_wide <- data_fam_no0s %>% 
  dplyr::select(site_code, family, weight_fam_sum_g) %>% 
  group_by(site_code, family) %>%
  summarise(total = sum(weight_fam_sum_g)) %>%
  ungroup() %>%
  spread(key = family, value = total) %>%
  replace(is.na(.), 0) %>%
  select(-17) # cut the last column, it's NA

# first use the rls final to filter for included surveys and make sure the order of rows is the same  
kelp_com <- data %>%
  select(site_code) %>%
  unique() %>%
  left_join(fam_kelp_wide, by = "site_code") %>%
  select(-1)

# make env data
kelp_env <- data %>%
  select(in_out_avg, kelp_bio_scale, kelp_sp, depth_scale, tide_scale) %>%
  unique() 

#convert com to a matrix
kelp_m_com = as.matrix(kelp_com)

# Perform the NMDS ordination
set.seed(123)
kelp_nmds = metaMDS(kelp_m_com, distance = "bray")
kelp_nmds # stress is fine

# Now we run the envfit function with our environmental data frame, env
kelp_en = envfit(kelp_nmds, kelp_env, permutations = 999, na.rm = TRUE)
# The first parameter is the metaMDS object from the NMDS ordination we just performed. Next is env, our environmental data frame. Then we state we want 999 permutations, and to remove any rows with missing data.
kelp_en

plot(kelp_nmds)
plot(kelp_en)

plot(kelp_nmds, type = "t")
plot(kelp_en)


# try plotting in ggplot? this code doesn't work
#save NMDS results into dataframe
site_scrs <- as.data.frame(scores(kelp_nmds, display = "sites")) 

spp_scrs <- as.data.frame(scores(spp_fit, display = "vectors")) #save species intrinsic values into dataframe
spp_scrs <- cbind(spp_scrs, Species = rownames(spp_scrs))


# Stats: site to site variation ----

# This was mostly exploratory, we didn't design the kelp study to ask questions about among-site differences

data_site <- data %>%
  select(site_code, site_name, nh4_avg, nh4_in_avg, nh4_out_avg, in_out_avg, Depth_m, mean_bio_og, mean_bio_recalc, BiomassM, DensityM, Area_m2, kelp_sp, weight_sum, abundance, shannon, simpson, avg_exchange_rate, forest_biomass, den_scale, kelp_bio_scale, forest_bio_scale, area_scale, weight_sum_scale, abundance_scale, shannon_scale, simpson_scale, depth_scale, tide_scale, tide_cat) %>%
  mutate(depth_scale = c(scale(Depth_m))) %>%
  unique()

# Step 0: Ask my question
# Is there a link between biodiversity or biomass and ammonium concentration? 

# Step 1: Choose a distribution

# response variables:
# nh4_avg (average of in and out)
# nh4_in_avg (inside kelp avg)
# nh4_out_avg (outside kelp avg)
# most comparable to the RLS blitz work?

ggplot(data, aes(x = nh4_out_avg)) +
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

mod_site <- glmmTMB(nh4_out_avg ~ kelp_sp + 
                        kelp_bio_scale*tide_scale +
                        kelp_bio_scale*weight_sum_scale +
                        weight_sum_scale*tide_scale +
                        shannon_scale + depth_scale,
                      family = 'gaussian',
                      data = data_site) 
summary(mod_site)
plot(simulateResiduals(mod_site))
visreg(mod_site, "kelp_bio_scale", by = "tide_scale")

# try fewer interactions, idk i have the data for that
mod_site <- glmmTMB(nh4_out_avg ~ kelp_sp + 
                      kelp_bio_scale +
                      weight_sum_scale*tide_scale +
                      shannon_scale + depth_scale,
                    family = 'gaussian',
                    data = data_site) 
summary(mod_site)
plot(simulateResiduals(mod_site))

# Build full model with gamma distribution
mod_kelp_site <- glmmTMB(nh4_out_avg ~ weight_sum_scale*tide_scale*kelp_bio_scale + shannon_scale + kelp_sp + depth_scale - weight_sum_scale:tide_scale:kelp_bio_scale, 
                         family = Gamma(link = 'log'),
                         data = data_site,
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
plot(simulateResiduals(mod_kelp_site_ints))

AIC(mod_kelp_site, mod_kelp_site_ints) # mod_kelp_site is better

# Finalize which predictors are best
# dredge <- as.data.frame(dredge(mod_gam_kelp)) %>% filter(delta < 3)
# best model: just kelp_sp + tide_scale

mod_top <- glmmTMB(nh4_out_avg ~ kelp_sp + tide_scale, 
                   family = Gamma(link = 'log'),
                   data = data_site)
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


# messing around -----

# pivot longer to plot kelp biomass of each species on same plot
long <- data %>%
  pivot_longer(cols = Macro:Nereo, values_to = "kelp_sp_density", names_to = "kelp_sp_sp")

# look at kelp biomass at each site
ggplot(long, aes(kelp_sp_density, site_name, colour = kelp_sp_sp)) +
  geom_point() 

# look at this another way 
short <- data %>% select(site_name, kelp_sp, Macro, Nereo)
# all Nereo sites are true Nereo, the Macro sites have very little Nereo
# Ed King East Inside is the only site I'd call truely mixed

ggplot(data %>% mutate(mixed = ifelse(site_name == "Ed King East Inside", "mixed", kelp_sp)), aes(kelp_bio_scale, in_minus_out, colour = mixed)) +
  geom_point() # ok so the one "mixed site" is also the highest kelp biomass.... idk i think I'm gonna have to just accept this


# looks like there's a positive relationship between shannon and delta nh4+
ggplot(data, aes(y = in_minus_out, x=shannon, colour = kelp_sp)) +
  geom_point() +
  geom_smooth(method = lm)


# try dropping not scaling the variables
mod_div <- glmmTMB(in_minus_out ~ kelp_sp + 
                     BiomassM*avg_exchange_rate +
                     BiomassM*weight_sum +
                     weight_sum*avg_exchange_rate +
                     BiomassM*shannon + depth_avg + 
                     (1|site_code),
                      family = 'gaussian',
                      data = data) 

summary(mod_div) # but when I remove the no kelp sites the trend is still negative
# but when you remove the kelp biomass shannon becomes positive









# when i remove weight shannon becomes + and NOT signif
# when I remove kelp biomass it becomes + and not signif

d <- data %>% filter(kelp_cat != "Mid") %>% as.data.frame()
mod <- glmmTMB(in_minus_out ~ kelp_sp + 
                 kelp_bio_scale*tide_scale +
                 kelp_bio_scale*weight_sum_scale +
                 weight_sum_scale*tide_scale +
                 shannon_scale + depth_scale + 
                 (1|site_code),
               family = 'gaussian',
               data = data) 

summary(mod)

# what's happening with diversity
ggplot(data = data,
       aes(x = species_richness, y = in_minus_out), colour = "red") +
  geom_point() +
  geom_smooth(method = 'lm')


visreg(mod)


# old
ggplot(data) +
  geom_point(aes(BiomassM, nh4_inside, shape = kelp_sp, colour = "inside")) +
  geom_point(aes(BiomassM, nh4_outside, shape = kelp_sp, colour = "outside")) +
  geom_point(aes(BiomassM, in_minus_out, shape = kelp_sp, colour = "delta"), size = 3) +
  geom_hline(yintercept= 0, linetype = "dashed", linewidth = 0.5) 


kcca_table <- pee %>%
  select(site_code, date) %>%
  unique() %>%
  filter(site_code != "KCCA13")


# Can i fit a asymptote??

data1 <- data %>%
  mutate(invert_in_out = 1/in_minus_out,
         invert_max_pee = 1/max(in_minus_out),
         invert_biomass = 1/BiomassM) %>%
  filter(BiomassM > 0)

lm(invert_mass ~ (k/maxC)*invert_Hct + invert_maxC)


mod_invert <- glmmTMB(invert_in_out ~ constant*invert_biomass + invert_max_pee, data = data1)

mod_invert2 <- glmmTMB(invert_in_out ~ invert_biomass, data = data1)
summary(mod_invert2)

data2 <- data1 %>%
  mutate(k = 8.4565/2)

mod_invert3 <- glmmTMB(invert_in_out ~ k*invert_max_pee*invert_biomass + invert_max_pee, data = data2)
summary(mod_invert3)
visreg(mod_invert3)

ggplot(data, aes(BiomassM, in_minus_out)) +
  geom_point() +
  geom_smooth()


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
          pal = pal3)


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



# try to figure out code diffs -----

data_check <- data %>%
  select(in_minus_out, kelp_sp, kelp_bio_scale, tide_scale, weight_sum_scale, shannon_scale, depth_scale,site_code, BiomassM)

data_old_check <- data_old %>%
  select(in_minus_out, kelp_sp, kelp_bio_scale, tide_scale, weight_sum_scale, shannon_scale, depth_scale,site_code, BiomassM)

d <- data_check == data_old_check
