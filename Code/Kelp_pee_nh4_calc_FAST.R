# NEW Script to calculate ammonium in vs outside a kelp forest
# July 4, 2022
# Em Lim

# This code only calculates the concentration of ammonium in each sample and generates a csv
# The actual analysis is completed in "Kelp_pee_analysis.R"

# Load packages -----
library(tidyverse)
library(visreg)

# Source super cool function
source("Code/Functions.R")

# Now read allllllll the site data and let the magic function do all the work
samples_kc1 <- pee_calc("Data/Team_kelp/2022_07_04_KCCA1.csv")
samples_kc2 <- pee_calc("Data/Team_kelp/2022_07_05_KCCA2.csv")
samples_kc3 <- pee_calc("Data/Team_kelp/2022_07_06_KCCA3.csv")
samples_kc4 <- pee_calc("Data/Team_kelp/2022_07_07_KCCA4.csv")
samples_kc6 <- pee_calc("Data/Team_kelp/2022_07_24_KCCA6.csv")
samples_kc7 <- pee_calc("Data/Team_kelp/2022_07_25_KCCA7.csv")
samples_kc9 <- pee_calc("Data/Team_kelp/2022_07_27_KCCA9.csv")
samples_kc12 <- pee_calc("Data/Team_kelp/2022_08_03_KCCA12.csv")
samples_kc13 <- pee_calc("Data/Team_kelp/2022_08_05_KCCA13.csv")
samples_kc14 <- pee_calc("Data/Team_kelp/2022_08_06_KCCA14.csv")
samples_kc15 <- pee_calc("Data/Team_kelp/2022_08_07_KCCA15.csv")
samples_kc16 <- pee_calc("Data/Team_kelp/2022_08_18_KCCA16.csv")
samples_kc17 <- pee_calc("Data/Team_kelp/2022_08_20_KCCA17.csv")
samples_kc18 <- pee_calc("Data/Team_kelp/2022_08_21_KCCA18.csv") # The R2 is only 0.9683819, look into KC18!!!
samples_kc19 <- pee_calc("Data/Team_kelp/2022_08_22_KCCA19.csv")
samples_kc21 <- pee_calc("Data/Team_kelp/2022_09_01_KCCA21.csv")
samples_kc22 <- pee_calc("Data/Team_kelp/2022_09_05_KCCA22.csv")


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
         sample = as.numeric(sample)) %>%
  select(-kelp)

# write_csv(kcca_not_final, "Data/Team_kelp/Output_data/kelp_pee.csv")
