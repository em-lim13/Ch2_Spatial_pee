# Script to look into the variation in pee inside vs outside kelp forests
# March 3, 2023
# Em Lim

# Load packages and functions -----
library(tidyverse)
library(MuMIn)
library(visreg)
library(ggplot2)
library(lmerTest)
library(lme4)
library(lubridate)
library(DHARMa)

# set theme and load functions
theme_set(theme_bw())
source("Code/theme_black.R")
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
  mutate(BiomassM = mean(Biomassm2kg)) %>%
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


max(data$nh4_avg)
min(data$nh4_avg) #[rls_nh4$nh4_conc > 0])

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
  length_to_weight() # length to weight function

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
         depth_avg = mean(depth_m)) %>%
  ungroup() %>%
  # join kelp data with both transect + site avgs
  left_join(kelp, by = c("site_code", "sample")) %>%
  # join site survey info for the survey ID
  left_join(kcca_survey_info %>%
              select(-c(Date, site_name)), by = "site_code"
            ) %>%
  # join tide exchange data
  left_join(tide_exchange_kcca, by = "survey_id") %>%
  # Join fish biomass
  left_join(kelp_rls %>%
              drop_na(weight_size_class_sum) %>%
              group_by(site_code) %>%
              summarize(weight_sum = sum(weight_size_class_sum))) %>%
  # Join richness and abundance
  left_join(kelp_rls %>%
              group_by(site_code) %>%
              summarize(species_richness = n_distinct(species_name),
                        abundance = sum(total))) %>%
  # 
  mutate(forest_biomass = BiomassM*Area_m2,
         percent_diff = 100*(nh4_inside-nh4_outside)/nh4_outside,
         den_scale = scale(kelp_den), # make sure density is right
         bio_tran_scale = scale(Biomassm2kg),
         bio_mean_scale = scale(BiomassM),
         forest_bio_scale = scale(forest_biomass),
         area_scale = scale(Area_m2),
         log_pee_diff = log(in_minus_out + 1),
         weight_stand = scale(weight_sum),
         rich_stand = scale(species_richness),
         abundance_stand = scale(abundance),
         tide_stand = scale(avg_exchange_rate),
         depth_stand = scale(depth_avg)) %>%
  filter(site_name != "Second Beach South")

# Each row is a mini transect into the kelp forest
  # That has an individual kelp density + biomass estimate + inside vs outside ammonium value

# new df for one row per site

data_s <- data %>%
  select(site, site_code, survey_id, nh4_in_avg, nh4_out_avg, nh4_avg, depth_avg, avg_exchange_rate, BiomassM, bio_mean_scale, weight_stand, rich_stand, abundance_stand, tide_stand,  depth_stand) %>%
  unique()

data_s_reduced <- data %>%
  select(site, site_code, survey_id, nh4_in_avg, depth_avg, avg_exchange_rate, BiomassM, weight_sum, species_richness, abundance, avg_exchange_rate,  depth_avg) %>%
  unique() %>%
  rename(nh4_avg = nh4_in_avg)
  

# Stats -----

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
  

# boxplot site vs pee diff
ggplot(data, aes(site, in_minus_out)) +
  geom_boxplot()

# each point is a single transect with a pee difference and a kelp density 
# linear model
ggplot(data, aes(kelp_den, in_minus_out)) +
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