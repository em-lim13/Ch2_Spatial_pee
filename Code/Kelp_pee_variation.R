# Script to determine whether ammonium is higher or lower in a kelp forest
# July 4, 2022
# Em Lim

# Load packages -----
library(tidyverse)
library(visreg)
library(ggplot2)
library(lmerTest)
library(lme4)
theme_set(theme_bw())
source("Code/theme_black.R")

get.ammonium <- function(data_path) {
  data_kc1_bf <- read_csv(data_path) %>%
    filter(sample == "bf") %>% # just the BF sample
    transmute(bf = (FLU1 + FLU2 + FLU3)/3,
              site_code = site_code) 

  # Load the standards and the sample data
  data_kc1 <- read_csv(data_path) %>%
    as.data.frame() %>%
    filter(sample != "bf") %>% # cut the bf sample
    left_join(data_kc1_bf, by = "site_code") %>% # Add the background flu as a column
    mutate(mean_FLU_unadj = (FLU1 + FLU2 + FLU3)/3,
           mean_FLU = mean_FLU_unadj - bf) # correct the samples with this bf value

  # Build the standard curve for protocol 1!
  # We need to estimate how much ammonium is in the standard water
  standard_kc1 <- data_kc1 %>%
    filter(sample == "standard") %>% 
    mutate(nh4_added_umol = nh4_vol_uL/1e6 * 200, #amount of NH4
           total_vol_L = nh4_vol_uL/1e6 + 0.04, #new volume of sample + NH4
           nh4_conc_final_umol_L = nh4_added_umol / total_vol_L) #concentration
  #of NH4 in seawater sample

  #graph standard curve
  print(ggplot(standard_kc1, aes(nh4_conc_final_umol_L, mean_FLU)) +
    geom_point() +
    geom_smooth(method = lm, se = FALSE))
  # Inspect this curve to make sure nothing is wonky!

  #model for protocol 1
  mod_kc1 <- lm(mean_FLU ~ nh4_conc_final_umol_L, data = standard_kc1)
  summary(mod_kc1)
  # Check Adjusted R-squared value, we're looking for at least 0.98

  #save coefficients
  int_kc1 <- coef(mod_kc1)[1]
  slope_kc1 <- coef(mod_kc1)[2]

  # calculate ammonium in standards using these coefficients
  standard_b_kc1 <- standard_kc1 %>%
    mutate(int = int_kc1, 
           slope = slope_kc1, 
           nh4_conc = abs(int/slope),
           true_nh4_conc = nh4_conc + nh4_conc_final_umol_L)

  # Now move into protocol 2! Use the standards to calculate ammonium in samples
  # BF-corrected FLU of standard curve against calculated NH4 concentration
  mod2_kc1 <- lm(mean_FLU ~ true_nh4_conc, data = standard_b_kc1)
  summary(mod2_kc1)

  # pull coefficients from model
  int2_kc1 <- coef(mod2_kc1)[1]
  slope2_kc1 <- coef(mod2_kc1)[2]

  # calculate how much ammonium is in the samples
  samples_kc1 <- data_kc1 %>%
    filter(sample != "standard") %>%
    mutate(nh4_conc = (mean_FLU - int2_kc1)/slope2_kc1)
}

samples_kc1 <- get.ammonium("Data/Team_kelp/2022_07_04_KCCA1.csv")
samples_kc2 <- get.ammonium("Data/Team_kelp/2022_07_05_KCCA2.csv")
samples_kc3 <- get.ammonium("Data/Team_kelp/2022_07_06_KCCA3.csv")
samples_kc4 <- get.ammonium("Data/Team_kelp/2022_07_07_KCCA4.csv")
samples_kc6 <- get.ammonium("Data/Team_kelp/2022_07_24_KCCA6.csv")
samples_kc7 <- get.ammonium("Data/Team_kelp/2022_07_25_KCCA7.csv")
samples_kc9 <- get.ammonium("Data/Team_kelp/2022_07_27_KCCA9.csv")
samples_kc12 <- get.ammonium("Data/Team_kelp/2022_08_03_KCCA12.csv")
samples_kc13 <- get.ammonium("Data/Team_kelp/2022_08_05_KCCA13.csv")
samples_kc14 <- get.ammonium("Data/Team_kelp/2022_08_06_KCCA14.csv")
samples_kc15 <- get.ammonium("Data/Team_kelp/2022_08_07_KCCA15.csv")
samples_kc16 <- get.ammonium("Data/Team_kelp/2022_08_18_KCCA16.csv")
samples_kc17 <- get.ammonium("Data/Team_kelp/2022_08_20_KCCA17.csv")
samples_kc18 <- get.ammonium("Data/Team_kelp/2022_08_21_KCCA18.csv")
samples_kc19 <- get.ammonium("Data/Team_kelp/2022_08_22_KCCA19.csv")
samples_kc21 <- get.ammonium("Data/Team_kelp/2022_09_01_KCCA21.csv")
samples_kc22 <- get.ammonium("Data/Team_kelp/2022_09_05_KCCA22.csv")

# Pull all the data together ----
kcca <- rbind(samples_kc1, samples_kc2, samples_kc3, samples_kc4, 
              samples_kc6, samples_kc7, samples_kc9, samples_kc12, 
              samples_kc13, samples_kc14, samples_kc15, samples_kc16,
              samples_kc17, samples_kc18, samples_kc19, samples_kc21,
              samples_kc22) %>%
  select(date, site, site_code, sample, kelp, depth_m, nh4_conc)

# separate the inside vs outside so I can put them back together side by side as separate columns
kcca_outside <- kcca %>%
  filter(kelp == "outside") %>%
  select(site_code, sample, nh4_conc) %>%
  rename(nh4_outside = nh4_conc)

kcca_not_final <- kcca %>%
  filter(kelp == "inside") %>%
  rename(nh4_inside = nh4_conc) %>%
  left_join(kcca_outside, by= c("site_code", "sample")) %>%
  mutate(in_minus_out = nh4_inside - nh4_outside,
         sample = as.numeric(sample))

# Kelp density data -----
# Load the kelp density Data from Claire
kelp <- read_csv("Data/Team_kelp/kelp_density_2022_KDC_CMA.csv") %>%
  as.data.frame() %>%
  filter(Transect_dist != "15") %>% # Remove the transect that's not paired with a pee sample
  mutate(sample = ifelse(Transect_dist == 0, 1, 
                         ifelse(Transect_dist == 5, 2, 3)),
         kelp_den = Macro_5m2 + Nereo_5m2,
         kelp_sp = ifelse(kelp_den == 0, "none",
                          ifelse(Macro_5m2 == 0, "nereo", 
                          ifelse(Nereo_5m2 == "0", "macro", "mixed"))))%>%
  select(site_code, sample, kelp_sp, kelp_den)

kelp_summary <- kelp %>%
  group_by(site_code) %>%
  summarise(avg_kelp_den = mean(kelp_den)) 

# merge kelp data with the pee data
kcca_final <- kcca_not_final %>%
  left_join(kelp, by = c("site_code", "sample")) %>%
  group_by(site) %>%
  mutate(nh4_avg = mean(c(nh4_outside, nh4_inside)),
         site = as.factor(site),
         percent_diff = 100*(nh4_inside-nh4_outside)/nh4_outside) %>%
  left_join(kelp_summary, by = "site_code")

# Average kelp density for each forest... I should probably go back and add the fourth transect to this avg
# Also Claire's biomass estimates
# Basically it would be nice to have a single metric of how much kelp is in each forest
kcca_summary <- kcca_final %>%
  group_by(site_code) %>%
  summarise(kelp_den = mean(kelp_den),
            in_minus_out = mean(in_minus_out),
            kelp_sp = kelp_sp)


# Plots ----

# boxplot site vs pee diff
ggplot(kcca_final, aes(site, in_minus_out)) +
  geom_boxplot()

# each point is a single transect with a pee difference and a kelp density 
# linear model
ggplot(kcca_final, aes(kelp_den, in_minus_out)) +
  geom_point(aes(pch = kelp_sp, colour = site_code), size =3 )+ 
  geom_hline(yintercept= 0, linetype = "dashed", color = "red", size = 1.5) +
  labs(y = "Inside - outside kelp forest ammonium (uM)", x = "Kelp Density") +
  geom_smooth(method = lm,
              alpha = 0.25)

# try plotting the same data but with a better asymptote?
ggplot(kcca_final, aes(kelp_den, in_minus_out)) +
  geom_point(aes(pch = kelp_sp, colour = site_code), size =3 )+ 
  geom_hline(yintercept= 0, linetype = "dashed", color = "red", size = 1.5) +
  labs(y = "Inside - outside kelp forest ammonium (uM)", x = "Kelp Density") +
  geom_smooth(method = "loess",
              span = 1,
              alpha = 0.25)

# Percent difference 
# With a truly heinious curve, thanks loess
ggplot(kcca_final, aes(kelp_den, percent_diff)) +
  geom_point(aes(pch = kelp_sp, colour = site_code), size =3 )+ 
  geom_hline(yintercept= 0, linetype = "dashed", color = "red", size = 1.5) +
  labs(y = "Inside - outside kelp forest ammonium (uM)", x = "Kelp Density") +
  geom_smooth(method = "loess",
              span = 1,
              alpha = 0.25)

#ggsave("Output/Figures/in_out_kelp_pee.png", device = "png",
#       height = 9, width = 16, dpi = 400)

# Ammonium at the site level?
kcca_final %>%
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
kelp_pee_mod <- lmer(in_minus_out ~ kelp_den + (1|site), data = kcca_final)
summary(kelp_pee_mod)
visreg(kelp_pee_mod)
