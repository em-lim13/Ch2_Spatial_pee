# Script to look into the variation in pee inside vs outside kelp forests
# March 3, 2023
# Em Lim

# Load packages
library(tidyverse)
library(MuMIn)
library(visreg)
library(ggplot2)
library(lmerTest)
library(lme4)
theme_set(theme_bw())
source("Code/theme_black.R")

# Load data

# load site names
names <- read_csv("Data/Team_kelp/Output_data/site_names.csv")

# kelp pee inside vs outside data from Kelp_pee_nh4_calc.R
pee <- read_csv("Data/Team_kelp/Output_data/kelp_pee.csv")

# kelp density data from Claire
den <- read_csv("Data/Team_kelp/Output_data/kelp_density_2022_KDC_CMA.csv") %>%
  as.data.frame() %>%
  filter(Transect_dist != "15") %>% # Remove the transect that's not paired with a pee sample
  mutate(sample = ifelse(Transect_dist == 0, 1, 
                  ifelse(Transect_dist == 5, 2, 3)),
         kelp_den = Macro_5m2 + Nereo_5m2,
         kelp_sp = ifelse(kelp_den == 0, "none",
                   ifelse(Macro_5m2 == 0, "nereo", 
                   ifelse(Nereo_5m2 == "0", "macro", "mixed"))))%>%
  select(site_code, sample, kelp_sp, kelp_den)

# kelp biomass data from Claire
bio <- read_csv("Data/Team_kelp/Output_data/transect_biomass.csv") %>%
  as.data.frame() %>%
  filter(Transect != "15") %>% 
  mutate(sample = ifelse(Transect == 0, 1, 
                  ifelse(Transect == 5, 2, 3))) %>%
  left_join(names, by = "SiteName") %>%
  select(site_code, SiteName, HeightT, BiomassTkg, Biomassm2kg, sample)
# BiomassTkg is kelp density x average biomass for each transect
# Biomassm2kg is just that number / 5 bc the transect was 5 m2

# Site level stuff from Claire
site <- read_csv("Data/Team_kelp/Output_data/kelpmetrics_2022.csv") %>%
  select(SiteName, BiomassM, DensityM, Composition, Area_m2)
# BiomassM is the average biomass/m2 at each site across all 4 transects
# DensityM is the average density at each site across all 4 transects


# Load RLS data
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
         date_time = ymd_hms(paste(Date, Time))) %>% 
  filter(Species != "Debris - Metal") %>%
  filter(Species != "Debris - Other") %>%
  filter(Species != "Debris - Wood") %>%
  filter(Species != "Debris - Glass") %>%
  filter(Species != "Debris - Fishing Gear") %>%
  filter(Species != "Debris - Zero")  # Cut the non-animal species

# Pivot longer
kelp_rls <- kelp_rls1 %>%
  rename(`0` = Inverts) %>%
  pivot_longer( cols = `0`:`400`, names_to = "size_class", values_to = "total") %>%
  drop_na(total) %>%
  filter(total > 0) %>%
  select(-Total)

# extract just one row per survey to join with the pee data and tide data
# Only keep the surveys I have pee samples for!
kcca_survey_info <- pee %>%
  transmute(Date = date,
            site_code = site_code) %>%
  unique() %>%
  left_join(kelp_rls %>%
              select(site_code, site_name, Date, date_time) %>%
              unique() %>%
              mutate(survey_id = 101:127), by = c("Date", "site_code"))
  
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
  
  survey_start <- ymd_hms(kcca_survey_info$date_time[x:x])
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
  left_join(den, by = c("site_code", "sample")) %>%
  left_join(bio, by = c("site_code", "sample")) %>%
  left_join(site, by = "SiteName") %>%
  left_join(kcca_survey_info %>%
              select(-c(Date, site_name)), by = "site_code"
            ) %>%
  left_join(tide_exchange_kcca, by = "survey_id") %>%
  replace(is.na(.), 0) %>% # replace NA with 0 bc those sites had no kelp
  mutate(forest_biomass = BiomassM*Area_m2,
         nh4_avg = mean(c(nh4_outside, nh4_inside)),
         percent_diff = 100*(nh4_inside-nh4_outside)/nh4_outside,
         den_scale = base::scale(kelp_den),
         bio_tran_scale = base::scale(Biomassm2kg),
         bio_mean_scale = base::scale(BiomassM),
         forest_bio_scale = base::scale(forest_biomass),
         area_scale = base::scale(Area_m2),
         log_pee_diff = log(in_minus_out + 1)) %>%
  filter(SiteName != "Second Beach South")


# OK let's do some friggin model dredging -----

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
bio_mod2 <- lmer(log_pee_diff ~ bio_mean_scale + (1|site), data = data)
summary(bio_mod2)
visreg(bio_mod2)

# Use AIC to see if this improves things
AIC(bio_mod, bio_mod2, mod1)
# Yes bio_mod2 is preferred, taking the log of the response variable helps a lot

# Ok let's plot this




# Old Kelp density data -----

# Also Claire's biomass estimates
# Basically it would be nice to have a single metric of how much kelp is in each forest
kcca_summary <- kcca_final %>%
  group_by(site_code) %>%
  summarise(kelp_den = mean(kelp_den),
            in_minus_out = mean(in_minus_out),
            kelp_sp = kelp_sp)


# Plots ----

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

# stats! ----
# I don't think this is linear I'm going to have to fit some kind of curve to it.....

# linear model though for fun
kelp_pee_mod <- lmer(in_minus_out ~ kelp_den + (1|site), data = data)
summary(kelp_pee_mod)
visreg(kelp_pee_mod)



# Claire has area, density, biomass
# Got average individual biomass for each kelp, which you can scale up to the average biomass of each transect


