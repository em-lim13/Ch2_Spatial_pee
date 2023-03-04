# Script to calculate ammonium in vs outside a kelp forest
# July 4, 2022
# Em Lim

# This code only calculates the concentration of ammonium in each sample and generates a csv
# The actual analysis is completed in "Kelp_pee_analysis.R"


# Load packages -----
library(tidyverse)

# KCCA1 Slug Island -----
# this data was read on the fluorometer twice
# once 5 hours after the OPA spike and again the next morning
# this code only deals with the immediate reading
# for the next morning reading see ME_BF_analysis 

# Load the background fluorometry reading
# This is how much a sample glows when the OPA is added immediately before being read on the fluorometer
data_kc1_bf <- read_csv("Data/Team_kelp/2022_07_04_KCCA1.csv") %>%
  filter(sample == "bf") %>% # just the BF sample
  transmute(bf = (FLU1 + FLU2 + FLU3)/3,
            site_code = site_code) 

# Load the standards and the sample data
data_kc1 <- read_csv("Data/Team_kelp/2022_07_04_KCCA1.csv") %>%
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
ggplot(standard_kc1, aes(nh4_conc_final_umol_L, mean_FLU)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)
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


# KCCA2 Between Brady's and Scotts -----
# this data was read on the fluorometer twice
# once 5 hours after the OPA spike and again the next morning
# this code only deals with the immediate reading
# for the next morning reading see ME_BF_analysis 

# Load the background fluorometry reading
data_kc2_bf <- read_csv("Data/Team_kelp/2022_07_05_KCCA2.csv") %>%
  filter(sample == "bf") %>% # just the BF sample
  transmute(bf = (FLU1 + FLU2 + FLU3)/3,
            site_code = site_code) 

# Load the standards and the sample data
data_kc2 <- read_csv("Data/Team_kelp/2022_07_05_KCCA2.csv") %>%
  as.data.frame() %>%
  left_join(data_kc2_bf, by = "site_code") %>%
  mutate(mean_FLU_unadj = (FLU1 + FLU2 + FLU3)/3,
         mean_FLU = mean_FLU_unadj - bf) %>%
  filter(sample != "bf")

# Build the standard curve for protocol 1!
standard_kc2 <- data_kc2 %>%
  filter(sample == "standard") %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * 200, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + 0.04, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L) #concentration
#of NH4 in seawater sample

#graph standard curve
ggplot(standard_kc2, aes(nh4_conc_final_umol_L, mean_FLU)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#model for protocol 1
mod_kc2 <- lm(mean_FLU ~ nh4_conc_final_umol_L, data = standard_kc2)
summary(mod_kc2)

#save coefficients
int_kc2 <- coef(mod_kc2)[1]
slope_kc2 <- coef(mod_kc2)[2]

# calculate ammonium in standards
standard_b_kc2 <- standard_kc2 %>%
  mutate(int = int_kc2, 
         slope = slope_kc2, 
         nh4_conc = abs(int/slope),
         true_nh4_conc = nh4_conc + nh4_conc_final_umol_L)

# Now move into protocol 2! Use the standards to calculate ammonium in samples
# Now BF-corrected FLU of standard curve against calculated NH4 concentration
mod2_kc2 <- lm(mean_FLU ~ true_nh4_conc, data = standard_b_kc2)
summary(mod2_kc2)

# pull coefficients from model
int2_kc2 <- coef(mod2_kc2)[1]
slope2_kc2 <- coef(mod2_kc2)[2]

# calculate how much ammonium is in the samples
samples_kc2 <- data_kc2 %>%
  filter(sample != "standard") %>%
  mutate(nh4_conc = (mean_FLU - int2_kc2)/slope2_kc2)

# KCCA3 Flemming 112 -----
# this data was read on the fluorometer twice
# once 5 hours after the OPA spike and again the next morning
# this code only deals with the immediate reading
# for the next morning reading see ME_BF_analysis 

# Load the background fluorometry reading
data_kc3_bf <- read_csv("Data/Team_kelp/2022_07_06_KCCA3.csv") %>%
  filter(sample == "bf") %>% # just the BF sample
  transmute(bf = (FLU1 + FLU2 + FLU3)/3,
            site_code = site_code) 

# Load the standards and the sample data
data_kc3 <- read_csv("Data/Team_kelp/2022_07_06_KCCA3.csv") %>%
  as.data.frame() %>%
  left_join(data_kc3_bf, by = "site_code") %>%
  mutate(mean_FLU_unadj = (FLU1 + FLU2 + FLU3)/3,
         mean_FLU = mean_FLU_unadj - bf) %>%
  filter(sample != "bf")

# Build the standard curve for protocol 1!
standard_kc3 <- data_kc3 %>%
  filter(sample == "standard") %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * 200, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + 0.04, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L) #concentration
#of NH4 in seawater sample

#graph standard curve
ggplot(standard_kc3, aes(nh4_conc_final_umol_L, mean_FLU)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#model for protocol 1
mod_kc3 <- lm(mean_FLU ~ nh4_conc_final_umol_L, data = standard_kc3)
summary(mod_kc3)

#save coefficients
int_kc3 <- coef(mod_kc3)[1]
slope_kc3 <- coef(mod_kc3)[2]

# calculate ammonium in standards
standard_b_kc3 <- standard_kc3 %>%
  mutate(int = int_kc3, 
         slope = slope_kc3, 
         nh4_conc = abs(int/slope),
         true_nh4_conc = nh4_conc + nh4_conc_final_umol_L)

# Now move into protocol 2! Use the standards to calculate ammonium in samples
# Now BF-corrected FLU of standard curve against calculated NH4 concentration
mod2_kc3 <- lm(mean_FLU ~ true_nh4_conc, data = standard_b_kc3)
summary(mod2_kc3)

# pull coefficients from model
int2_kc3 <- coef(mod2_kc3)[1]
slope2_kc3 <- coef(mod2_kc3)[2]

# calculate how much ammonium is in the samples
samples_kc3 <- data_kc3 %>%
  filter(sample != "standard") %>%
  mutate(nh4_conc = (mean_FLU - int2_kc3)/slope2_kc3)

# KCCA4 Fleming112 -----
# Load the background fluorometry reading
data_kc4_bf <- read_csv("Data/Team_kelp/2022_07_07_KCCA4.csv") %>%
  filter(sample == "bf") %>% # just the BF sample
  transmute(bf = (FLU1 + FLU2 + FLU3)/3,
            site_code = site_code) 

# Load the standards and the sample data
data_kc4 <- read_csv("Data/Team_kelp/2022_07_07_KCCA4.csv") %>%
  as.data.frame() %>%
  left_join(data_kc4_bf, by = "site_code") %>%
  mutate(mean_FLU_unadj = (FLU1 + FLU2 + FLU3)/3,
         mean_FLU = mean_FLU_unadj - bf) %>%
  filter(sample != "bf")

# Build the standard curve for protocol 1!
standard_kc4 <- data_kc4 %>%
  filter(sample == "standard") %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * 200, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + 0.04, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L) #concentration
#of NH4 in seawater sample

#graph standard curve
ggplot(standard_kc4, aes(nh4_conc_final_umol_L, mean_FLU)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#model for protocol 1
mod_kc4 <- lm(mean_FLU ~ nh4_conc_final_umol_L, data = standard_kc4)
summary(mod_kc4)

#save coefficients
int_kc4 <- coef(mod_kc4)[1]
slope_kc4 <- coef(mod_kc4)[2]

# calculate ammonium in standards
standard_b_kc4 <- standard_kc4 %>%
  mutate(int = int_kc4, 
         slope = slope_kc4, 
         nh4_conc = abs(int/slope),
         true_nh4_conc = nh4_conc + nh4_conc_final_umol_L)

# Now move into protocol 2! Use the standards to calculate ammonium in samples
# Now BF-corrected FLU of standard curve against calculated NH4 concentration
mod2_kc4 <- lm(mean_FLU ~ true_nh4_conc, data = standard_b_kc4)
summary(mod2_kc4)

# pull coefficients from model
int2_kc4 <- coef(mod2_kc4)[1]
slope2_kc4 <- coef(mod2_kc4)[2]

# calculate how much ammonium is in the samples
samples_kc4 <- data_kc4 %>%
  filter(sample != "standard") %>%
  mutate(nh4_conc = (mean_FLU - int2_kc4)/slope2_kc4)

# KCCA6 Less_dangerous_bay -----
# Load the background fluorometry reading
data_kc6_bf <- read_csv("Data/Team_kelp/2022_07_24_KCCA6.csv") %>%
  filter(sample == "bf") %>% # just the BF sample
  transmute(bf = (FLU1 + FLU2 + FLU3)/3,
            site_code = site_code) 

# Load the standards and the sample data
data_kc6 <- read_csv("Data/Team_kelp/2022_07_24_KCCA6.csv") %>%
  as.data.frame() %>%
  left_join(data_kc6_bf, by = "site_code") %>%
  mutate(mean_FLU_unadj = (FLU1 + FLU2 + FLU3)/3,
         mean_FLU = mean_FLU_unadj - bf) %>%
  filter(sample != "bf")

# Build the standard curve for protocol 1!
standard_kc6 <- data_kc6 %>%
  filter(sample == "standard") %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * 200, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + 0.04, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L) #concentration
#of NH4 in seawater sample

#graph standard curve
ggplot(standard_kc6, aes(nh4_conc_final_umol_L, mean_FLU)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#model for protocol 1
mod_kc6 <- lm(mean_FLU ~ nh4_conc_final_umol_L, data = standard_kc6)
summary(mod_kc6)

#save coefficients
int_kc6 <- coef(mod_kc6)[1]
slope_kc6 <- coef(mod_kc6)[2]

# calculate ammonium in standards
standard_b_kc6 <- standard_kc6 %>%
  mutate(int = int_kc6, 
         slope = slope_kc6, 
         nh4_conc = abs(int/slope),
         true_nh4_conc = nh4_conc + nh4_conc_final_umol_L)

# Now move into protocol 2! Use the standards to calculate ammonium in samples
# Now BF-corrected FLU of standard curve against calculated NH4 concentration
mod2_kc6 <- lm(mean_FLU ~ true_nh4_conc, data = standard_b_kc6)
summary(mod2_kc6)

# pull coefficients from model
int2_kc6 <- coef(mod2_kc6)[1]
slope2_kc6 <- coef(mod2_kc6)[2]

# calculate how much ammonium is in the samples
samples_kc6 <- data_kc6 %>%
  filter(sample != "standard") %>%
  mutate(nh4_conc = (mean_FLU - int2_kc6)/slope2_kc6)

# KCCA7 Ed_king_east_inside ----
# Load the background fluorometry reading
data_kc7_bf <- read_csv("Data/Team_kelp/2022_07_25_KCCA7.csv") %>%
  filter(sample == "bf") %>% # just the BF sample
  transmute(bf = (FLU1 + FLU2 + FLU3)/3,
            site_code = site_code) 

# Load the standards and the sample data
data_kc7 <- read_csv("Data/Team_kelp/2022_07_25_KCCA7.csv") %>%
  as.data.frame() %>%
  left_join(data_kc7_bf, by = "site_code") %>%
  mutate(mean_FLU_unadj = (FLU1 + FLU2 + FLU3)/3,
         mean_FLU = mean_FLU_unadj - bf) %>%
  filter(sample != "bf")

# Build the standard curve for protocol 1!
standard_kc7 <- data_kc7 %>%
  filter(sample == "standard") %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * 200, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + 0.04, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L) #concentration
#of NH4 in seawater sample

#graph standard curve
ggplot(standard_kc7, aes(nh4_conc_final_umol_L, mean_FLU)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#model for protocol 1
mod_kc7 <- lm(mean_FLU ~ nh4_conc_final_umol_L, data = standard_kc7)
summary(mod_kc7)

#save coefficients
int_kc7 <- coef(mod_kc7)[1]
slope_kc7 <- coef(mod_kc7)[2]

# calculate ammonium in standards
standard_b_kc7 <- standard_kc7 %>%
  mutate(int = int_kc7, 
         slope = slope_kc7, 
         nh4_conc = abs(int/slope),
         true_nh4_conc = nh4_conc + nh4_conc_final_umol_L)

# Now move into protocol 2! Use the standards to calculate ammonium in samples
# Now BF-corrected FLU of standard curve against calculated NH4 concentration
mod2_kc7 <- lm(mean_FLU ~ true_nh4_conc, data = standard_b_kc7)
summary(mod2_kc7)

# pull coefficients from model
int2_kc7 <- coef(mod2_kc7)[1]
slope2_kc7 <- coef(mod2_kc7)[2]

# calculate how much ammonium is in the samples
samples_kc7 <- data_kc7 %>%
  filter(sample != "standard") %>%
  mutate(nh4_conc = (mean_FLU - int2_kc7)/slope2_kc7)

# KCCA9 Wizard -----
# Load the background fluorometry reading
data_kc9_bf <- read_csv("Data/Team_kelp/2022_07_27_KCCA9.csv") %>%
  filter(sample == "bf") %>% # just the BF sample
  transmute(bf = (FLU1 + FLU2 + FLU3)/3,
            site_code = site_code) 

# Load the standards and the sample data
data_kc9 <- read_csv("Data/Team_kelp/2022_07_27_KCCA9.csv") %>%
  as.data.frame() %>%
  left_join(data_kc9_bf, by = "site_code") %>%
  mutate(mean_FLU_unadj = (FLU1 + FLU2 + FLU3)/3,
         mean_FLU = mean_FLU_unadj - bf) %>%
  filter(sample != "bf")

# Build the standard curve for protocol 1!
standard_kc9 <- data_kc9 %>%
  filter(sample == "standard") %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * 200, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + 0.04, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L) #concentration
#of NH4 in seawater sample

#graph standard curve
ggplot(standard_kc9, aes(nh4_conc_final_umol_L, mean_FLU)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#model for protocol 1
mod_kc9 <- lm(mean_FLU ~ nh4_conc_final_umol_L, data = standard_kc9)
summary(mod_kc9)

#save coefficients
int_kc9 <- coef(mod_kc9)[1]
slope_kc9 <- coef(mod_kc9)[2]

# calculate ammonium in standards
standard_b_kc9 <- standard_kc9 %>%
  mutate(int = int_kc9, 
         slope = slope_kc9, 
         nh4_conc = abs(int/slope),
         true_nh4_conc = nh4_conc + nh4_conc_final_umol_L)

# Now move into protocol 2! Use the standards to calculate ammonium in samples
# Now BF-corrected FLU of standard curve against calculated NH4 concentration
mod2_kc9 <- lm(mean_FLU ~ true_nh4_conc, data = standard_b_kc9)
summary(mod2_kc9)

# pull coefficients from model
int2_kc9 <- coef(mod2_kc9)[1]
slope2_kc9 <- coef(mod2_kc9)[2]

# calculate how much ammonium is in the samples
samples_kc9 <- data_kc9 %>%
  filter(sample != "standard") %>%
  mutate(nh4_conc = (mean_FLU - int2_kc9)/slope2_kc9)

# KCCA12 North_helby_rock -----
# Load the background fluorometry reading
data_kc12_bf <- read_csv("Data/Team_kelp/2022_08_03_KCCA12.csv") %>%
  filter(sample == "bf") %>% # just the BF sample
  transmute(bf = (FLU1 + FLU2 + FLU3)/3,
            site_code = site_code) 

# Load the standards and the sample data
data_kc12 <- read_csv("Data/Team_kelp/2022_08_03_KCCA12.csv") %>%
  as.data.frame() %>%
  left_join(data_kc12_bf, by = "site_code") %>%
  mutate(mean_FLU_unadj = (FLU1 + FLU2 + FLU3)/3,
         mean_FLU = mean_FLU_unadj - bf) %>%
  filter(sample != "bf")

# Build the standard curve for protocol 1!
standard_kc12 <- data_kc12 %>%
  filter(sample == "standard") %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * 200, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + 0.04, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L) #concentration
#of NH4 in seawater sample

#graph standard curve
ggplot(standard_kc12, aes(nh4_conc_final_umol_L, mean_FLU)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#model for protocol 1
mod_kc12 <- lm(mean_FLU ~ nh4_conc_final_umol_L, data = standard_kc12)
summary(mod_kc12)

#save coefficients
int_kc12 <- coef(mod_kc12)[1]
slope_kc12 <- coef(mod_kc12)[2]

# calculate ammonium in standards
standard_b_kc12 <- standard_kc12 %>%
  mutate(int = int_kc12, 
         slope = slope_kc12, 
         nh4_conc = abs(int/slope),
         true_nh4_conc = nh4_conc + nh4_conc_final_umol_L)

# Now move into protocol 2! Use the standards to calculate ammonium in samples
# Now BF-corrected FLU of standard curve against calculated NH4 concentration
mod2_kc12 <- lm(mean_FLU ~ true_nh4_conc, data = standard_b_kc12)
summary(mod2_kc12)

# pull coefficients from model
int2_kc12 <- coef(mod2_kc12)[1]
slope2_kc12 <- coef(mod2_kc12)[2]

# calculate how much ammonium is in the samples
samples_kc12 <- data_kc12 %>%
  filter(sample != "standard") %>%
  mutate(nh4_conc = (mean_FLU - int2_kc12)/slope2_kc12)

# KCCA13 Second_beach_south ----
# Load the background fluorometry reading
data_kc13_bf <- read_csv("Data/Team_kelp/2022_08_05_KCCA13.csv") %>%
  filter(sample == "bf") %>% # just the BF sample
  transmute(bf = (FLU1 + FLU2 + FLU3)/3,
            site_code = site_code) 

# Load the standards and the sample data
data_kc13 <- read_csv("Data/Team_kelp/2022_08_05_KCCA13.csv") %>%
  as.data.frame() %>%
  left_join(data_kc13_bf, by = "site_code") %>%
  mutate(mean_FLU_unadj = (FLU1 + FLU2 + FLU3)/3,
         mean_FLU = mean_FLU_unadj - bf) %>%
  filter(sample != "bf")

# Build the standard curve for protocol 1!
standard_kc13 <- data_kc13 %>%
  filter(sample == "standard") %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * 200, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + 0.04, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L) #concentration
#of NH4 in seawater sample

#graph standard curve
ggplot(standard_kc13, aes(nh4_conc_final_umol_L, mean_FLU)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#model for protocol 1
mod_kc13 <- lm(mean_FLU ~ nh4_conc_final_umol_L, data = standard_kc13)
summary(mod_kc13)

#save coefficients
int_kc13 <- coef(mod_kc13)[1]
slope_kc13 <- coef(mod_kc13)[2]

# calculate ammonium in standards
standard_b_kc13 <- standard_kc13 %>%
  mutate(int = int_kc13, 
         slope = slope_kc13, 
         nh4_conc = abs(int/slope),
         true_nh4_conc = nh4_conc + nh4_conc_final_umol_L)

# Now move into protocol 2! Use the standards to calculate ammonium in samples
# Now BF-corrected FLU of standard curve against calculated NH4 concentration
mod2_kc13 <- lm(mean_FLU ~ true_nh4_conc, data = standard_b_kc13)
summary(mod2_kc13)

# pull coefficients from model
int2_kc13 <- coef(mod2_kc13)[1]
slope2_kc13 <- coef(mod2_kc13)[2]

# calculate how much ammonium is in the samples
samples_kc13 <- data_kc13 %>%
  filter(sample != "standard") %>%
  mutate(nh4_conc = (mean_FLU - int2_kc13)/slope2_kc13)

# KCCA14 Danvers_danger_rock -----
# Load the background fluorometry reading
data_kc14_bf <- read_csv("Data/Team_kelp/2022_08_06_KCCA14.csv") %>%
  filter(sample == "bf") %>% # just the BF sample
  transmute(bf = (FLU1 + FLU2 + FLU3)/3,
            site_code = site_code) 

# Load the standards and the sample data
data_kc14 <- read_csv("Data/Team_kelp/2022_08_06_KCCA14.csv") %>%
  as.data.frame() %>%
  left_join(data_kc14_bf, by = "site_code") %>%
  mutate(mean_FLU_unadj = (FLU1 + FLU2 + FLU3)/3,
         mean_FLU = mean_FLU_unadj - bf) %>%
  filter(sample != "bf")

# Build the standard curve for protocol 1!
standard_kc14 <- data_kc14 %>%
  filter(sample == "standard") %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * 200, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + 0.04, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L) #concentration
#of NH4 in seawater sample

#graph standard curve
ggplot(standard_kc14, aes(nh4_conc_final_umol_L, mean_FLU)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#model for protocol 1
mod_kc14 <- lm(mean_FLU ~ nh4_conc_final_umol_L, data = standard_kc14)
summary(mod_kc14)

#save coefficients
int_kc14 <- coef(mod_kc14)[1]
slope_kc14 <- coef(mod_kc14)[2]

# calculate ammonium in standards
standard_b_kc14 <- standard_kc14 %>%
  mutate(int = int_kc14, 
         slope = slope_kc14, 
         nh4_conc = abs(int/slope),
         true_nh4_conc = nh4_conc + nh4_conc_final_umol_L)

# Now move into protocol 2! Use the standards to calculate ammonium in samples
# Now BF-corrected FLU of standard curve against calculated NH4 concentration
mod2_kc14 <- lm(mean_FLU ~ true_nh4_conc, data = standard_b_kc14)
summary(mod2_kc14)

# pull coefficients from model
int2_kc14 <- coef(mod2_kc14)[1]
slope2_kc14 <- coef(mod2_kc14)[2]

# calculate how much ammonium is in the samples
samples_kc14 <- data_kc14 %>%
  filter(sample != "standard") %>%
  mutate(nh4_conc = (mean_FLU - int2_kc14)/slope2_kc14)

# KCCA15 Cable_beach ----
# Load the background fluorometry reading
data_kc15_bf <- read_csv("Data/Team_kelp/2022_08_07_KCCA15.csv") %>%
  filter(sample == "bf") %>% # just the BF sample
  transmute(bf = (FLU1 + FLU2 + FLU3)/3,
            site_code = site_code) 

# Load the standards and the sample data
data_kc15 <- read_csv("Data/Team_kelp/2022_08_07_KCCA15.csv") %>%
  as.data.frame() %>%
  left_join(data_kc15_bf, by = "site_code") %>%
  mutate(mean_FLU_unadj = (FLU1 + FLU2 + FLU3)/3,
         mean_FLU = mean_FLU_unadj - bf) %>%
  filter(sample != "bf")

# Build the standard curve for protocol 1!
standard_kc15 <- data_kc15 %>%
  filter(sample == "standard") %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * 200, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + 0.04, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L) #concentration
#of NH4 in seawater sample

#graph standard curve
ggplot(standard_kc15, aes(nh4_conc_final_umol_L, mean_FLU)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#model for protocol 1
mod_kc15 <- lm(mean_FLU ~ nh4_conc_final_umol_L, data = standard_kc15)
summary(mod_kc15)

#save coefficients
int_kc15 <- coef(mod_kc15)[1]
slope_kc15 <- coef(mod_kc15)[2]

# calculate ammonium in standards
standard_b_kc15 <- standard_kc15 %>%
  mutate(int = int_kc15, 
         slope = slope_kc15, 
         nh4_conc = abs(int/slope),
         true_nh4_conc = nh4_conc + nh4_conc_final_umol_L)

# Now move into protocol 2! Use the standards to calculate ammonium in samples
# Now BF-corrected FLU of standard curve against calculated NH4 concentration
mod2_kc15 <- lm(mean_FLU ~ true_nh4_conc, data = standard_b_kc15)
summary(mod2_kc15)

# pull coefficients from model
int2_kc15 <- coef(mod2_kc15)[1]
slope2_kc15 <- coef(mod2_kc15)[2]

# calculate how much ammonium is in the samples
samples_kc15 <- data_kc15 %>%
  filter(sample != "standard") %>%
  mutate(nh4_conc = (mean_FLU - int2_kc15)/slope2_kc15)

# KCCA16 Tzartus 116 ----
# Load the background fluorometry reading
data_kc16_bf <- read_csv("Data/Team_kelp/2022_08_18_KCCA16.csv") %>%
  filter(sample == "bf") %>% # just the BF sample
  transmute(bf = (FLU1 + FLU2 + FLU3)/3,
            site_code = site_code) 

# Load the standards and the sample data
data_kc16 <- read_csv("Data/Team_kelp/2022_08_18_KCCA16.csv") %>%
  as.data.frame() %>%
  left_join(data_kc16_bf, by = "site_code") %>%
  mutate(mean_FLU_unadj = (FLU1 + FLU2 + FLU3)/3,
         mean_FLU = mean_FLU_unadj - bf) %>%
  filter(sample != "bf")

# Build the standard curve for protocol 1!
standard_kc16 <- data_kc16 %>%
  filter(sample == "standard") %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * 200, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + 0.04, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L) #concentration
#of NH4 in seawater sample

#graph standard curve
ggplot(standard_kc16, aes(nh4_conc_final_umol_L, mean_FLU)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#model for protocol 1
mod_kc16 <- lm(mean_FLU ~ nh4_conc_final_umol_L, data = standard_kc16)
summary(mod_kc16)

#save coefficients
int_kc16 <- coef(mod_kc16)[1]
slope_kc16 <- coef(mod_kc16)[2]

# calculate ammonium in standards
standard_b_kc16 <- standard_kc16 %>%
  mutate(int = int_kc16, 
         slope = slope_kc16, 
         nh4_conc = abs(int/slope),
         true_nh4_conc = nh4_conc + nh4_conc_final_umol_L)

# Now move into protocol 2! Use the standards to calculate ammonium in samples
# Now BF-corrected FLU of standard curve against calculated NH4 concentration
mod2_kc16 <- lm(mean_FLU ~ true_nh4_conc, data = standard_b_kc16)
summary(mod2_kc16)

# pull coefficients from model
int2_kc16 <- coef(mod2_kc16)[1]
slope2_kc16 <- coef(mod2_kc16)[2]

# calculate how much ammonium is in the samples
samples_kc16 <- data_kc16 %>%
  filter(sample != "standard") %>%
  mutate(nh4_conc = (mean_FLU - int2_kc16)/slope2_kc16)

# KCCA17 Turf Island ----
# Load the background fluorometry reading
data_kc17_bf <- read_csv("Data/Team_kelp/2022_08_20_KCCA17.csv") %>%
  filter(sample == "bf") %>% # just the BF sample
  transmute(bf = (FLU1 + FLU2 + FLU3)/3,
            site_code = site_code) 

# Load the standards and the sample data
data_kc17 <- read_csv("Data/Team_kelp/2022_08_20_KCCA17.csv") %>%
  as.data.frame() %>%
  left_join(data_kc17_bf, by = "site_code") %>%
  mutate(mean_FLU_unadj = (FLU1 + FLU2 + FLU3)/3,
         mean_FLU = mean_FLU_unadj - bf) %>%
  filter(sample != "bf")

# Build the standard curve for protocol 1!
standard_kc17 <- data_kc17 %>%
  filter(sample == "standard") %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * 200, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + 0.04, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L) #concentration
#of NH4 in seawater sample

#graph standard curve
ggplot(standard_kc17, aes(nh4_conc_final_umol_L, mean_FLU)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#model for protocol 1
mod_kc17 <- lm(mean_FLU ~ nh4_conc_final_umol_L, data = standard_kc17)
summary(mod_kc17)

#save coefficients
int_kc17 <- coef(mod_kc17)[1]
slope_kc17 <- coef(mod_kc17)[2]

# calculate ammonium in standards
standard_b_kc17 <- standard_kc17 %>%
  mutate(int = int_kc17, 
         slope = slope_kc17, 
         nh4_conc = abs(int/slope),
         true_nh4_conc = nh4_conc + nh4_conc_final_umol_L)

# Now move into protocol 2! Use the standards to calculate ammonium in samples
# Now BF-corrected FLU of standard curve against calculated NH4 concentration
mod2_kc17 <- lm(mean_FLU ~ true_nh4_conc, data = standard_b_kc17)
summary(mod2_kc17)

# pull coefficients from model
int2_kc17 <- coef(mod2_kc17)[1]
slope2_kc17 <- coef(mod2_kc17)[2]

# calculate how much ammonium is in the samples
samples_kc17 <- data_kc17 %>%
  filter(sample != "standard") %>%
  mutate(nh4_conc = (mean_FLU - int2_kc17)/slope2_kc17)

# KCCA18 Second Beach ----
# Load the background fluorometry reading
data_kc18_bf <- read_csv("Data/Team_kelp/2022_08_21_KCCA18.csv") %>%
  filter(sample == "bf") %>% # just the BF sample
  transmute(bf = (FLU1 + FLU2 + FLU3)/3,
            site_code = site_code) 

# Load the standards and the sample data
data_kc18 <- read_csv("Data/Team_kelp/2022_08_21_KCCA18.csv") %>%
  as.data.frame() %>%
  left_join(data_kc18_bf, by = "site_code") %>%
  mutate(mean_FLU_unadj = (FLU1 + FLU2 + FLU3)/3,
         mean_FLU = mean_FLU_unadj - bf) %>%
  filter(sample != "bf")

# Build the standard curve for protocol 1!
standard_kc18 <- data_kc18 %>%
  filter(sample == "standard") %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * 200, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + 0.04, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L) #concentration
#of NH4 in seawater sample

#graph standard curve
ggplot(standard_kc18, aes(nh4_conc_final_umol_L, mean_FLU)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#model for protocol 1
mod_kc18 <- lm(mean_FLU ~ nh4_conc_final_umol_L, data = standard_kc18)
summary(mod_kc18)
# DAMN THE CURVE LOOKS LIKE SHIT
# something is weird with the 800 spike, the readings aren't within 3 units
# Double check the data sheet

#save coefficients
int_kc18 <- coef(mod_kc18)[1]
slope_kc18 <- coef(mod_kc18)[2]

# calculate ammonium in standards
standard_b_kc18 <- standard_kc18 %>%
  mutate(int = int_kc18, 
         slope = slope_kc18, 
         nh4_conc = abs(int/slope),
         true_nh4_conc = nh4_conc + nh4_conc_final_umol_L)

# Now move into protocol 2! Use the standards to calculate ammonium in samples
# Now BF-corrected FLU of standard curve against calculated NH4 concentration
mod2_kc18 <- lm(mean_FLU ~ true_nh4_conc, data = standard_b_kc18)
summary(mod2_kc18)

# pull coefficients from model
int2_kc18 <- coef(mod2_kc18)[1]
slope2_kc18 <- coef(mod2_kc18)[2]

# calculate how much ammonium is in the samples
samples_kc18 <- data_kc18 %>%
  filter(sample != "standard") %>%
  mutate(nh4_conc = (mean_FLU - int2_kc18)/slope2_kc18)

# KCCA19 Second Beach ----
# Load the background fluorometry reading
data_kc19_bf <- read_csv("Data/Team_kelp/2022_08_22_KCCA19.csv") %>%
  filter(sample == "bf") %>% # just the BF sample
  transmute(bf = (FLU1 + FLU2 + FLU3)/3,
            site_code = site_code) 

# Load the standards and the sample data
data_kc19 <- read_csv("Data/Team_kelp/2022_08_22_KCCA19.csv") %>%
  as.data.frame() %>%
  left_join(data_kc19_bf, by = "site_code") %>%
  mutate(mean_FLU_unadj = (FLU1 + FLU2 + FLU3)/3,
         mean_FLU = mean_FLU_unadj - bf) %>%
  filter(sample != "bf")

# Build the standard curve for protocol 1!
standard_kc19 <- data_kc19 %>%
  filter(sample == "standard") %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * 200, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + 0.04, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L) #concentration
#of NH4 in seawater sample

#graph standard curve
ggplot(standard_kc19, aes(nh4_conc_final_umol_L, mean_FLU)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#model for protocol 1
mod_kc19 <- lm(mean_FLU ~ nh4_conc_final_umol_L, data = standard_kc19)
summary(mod_kc19)

#save coefficients
int_kc19 <- coef(mod_kc19)[1]
slope_kc19 <- coef(mod_kc19)[2]

# calculate ammonium in standards
standard_b_kc19 <- standard_kc19 %>%
  mutate(int = int_kc19, 
         slope = slope_kc19, 
         nh4_conc = abs(int/slope),
         true_nh4_conc = nh4_conc + nh4_conc_final_umol_L)

# Now move into protocol 2! Use the standards to calculate ammonium in samples
# Now BF-corrected FLU of standard curve against calculated NH4 concentration
mod2_kc19 <- lm(mean_FLU ~ true_nh4_conc, data = standard_b_kc19)
summary(mod2_kc19)

# pull coefficients from model
int2_kc19 <- coef(mod2_kc19)[1]
slope2_kc19 <- coef(mod2_kc19)[2]

# calculate how much ammonium is in the samples
samples_kc19 <- data_kc19 %>%
  filter(sample != "standard") %>%
  mutate(nh4_conc = (mean_FLU - int2_kc19)/slope2_kc19)

# KCCA21 Bordelais ----
# Load the background fluorometry reading
data_kc21_bf <- read_csv("Data/Team_kelp/2022_09_01_KCCA21.csv") %>%
  filter(sample == "bf") %>% # just the BF sample
  transmute(bf = (FLU1 + FLU2 + FLU3)/3,
            site_code = site_code) 

# Load the standards and the sample data
data_kc21 <- read_csv("Data/Team_kelp/2022_09_01_KCCA21.csv") %>%
  as.data.frame() %>%
  left_join(data_kc21_bf, by = "site_code") %>%
  mutate(mean_FLU_unadj = (FLU1 + FLU2 + FLU3)/3,
         mean_FLU = mean_FLU_unadj - bf) %>%
  filter(sample != "bf")

# Build the standard curve for protocol 1!
standard_kc21 <- data_kc21 %>%
  filter(sample == "standard") %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * 200, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + 0.04, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L) #concentration
#of NH4 in seawater sample

#graph standard curve
ggplot(standard_kc21, aes(nh4_conc_final_umol_L, mean_FLU)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#model for protocol 1
mod_kc21 <- lm(mean_FLU ~ nh4_conc_final_umol_L, data = standard_kc21)
summary(mod_kc21)

#save coefficients
int_kc21 <- coef(mod_kc21)[1]
slope_kc21 <- coef(mod_kc21)[2]

# calculate ammonium in standards
standard_b_kc21 <- standard_kc21 %>%
  mutate(int = int_kc21, 
         slope = slope_kc21, 
         nh4_conc = abs(int/slope),
         true_nh4_conc = nh4_conc + nh4_conc_final_umol_L)

# Now move into protocol 2! Use the standards to calculate ammonium in samples
# Now BF-corrected FLU of standard curve against calculated NH4 concentration
mod2_kc21 <- lm(mean_FLU ~ true_nh4_conc, data = standard_b_kc21)
summary(mod2_kc21)

# pull coefficients from model
int2_kc21 <- coef(mod2_kc21)[1]
slope2_kc21 <- coef(mod2_kc21)[2]

# calculate how much ammonium is in the samples
samples_kc21 <- data_kc21 %>%
  filter(sample != "standard") %>%
  mutate(nh4_conc = (mean_FLU - int2_kc21)/slope2_kc21)

# KCCA22 Taylor Rock ----
# Load the background fluorometry reading
data_kc22_bf <- read_csv("Data/Team_kelp/2022_09_05_KCCA22.csv") %>%
  filter(sample == "bf") %>% # just the BF sample
  transmute(bf = (FLU1 + FLU2 + FLU3)/3,
            site_code = site_code) 

# Load the standards and the sample data
data_kc22 <- read_csv("Data/Team_kelp/2022_09_05_KCCA22.csv") %>%
  as.data.frame() %>%
  left_join(data_kc22_bf, by = "site_code") %>%
  mutate(mean_FLU_unadj = (FLU1 + FLU2 + FLU3)/3,
         mean_FLU = mean_FLU_unadj - bf) %>%
  filter(sample != "bf")

# Build the standard curve for protocol 1!
standard_kc22 <- data_kc22 %>%
  filter(sample == "standard") %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * 200, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + 0.04, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L) #concentration
#of NH4 in seawater sample

#graph standard curve
ggplot(standard_kc22, aes(nh4_conc_final_umol_L, mean_FLU)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#model for protocol 1
mod_kc22 <- lm(mean_FLU ~ nh4_conc_final_umol_L, data = standard_kc22)
summary(mod_kc22)

#save coefficients
int_kc22 <- coef(mod_kc22)[1]
slope_kc22 <- coef(mod_kc22)[2]

# calculate ammonium in standards
standard_b_kc22 <- standard_kc22 %>%
  mutate(int = int_kc22, 
         slope = slope_kc22, 
         nh4_conc = abs(int/slope),
         true_nh4_conc = nh4_conc + nh4_conc_final_umol_L)

# Now move into protocol 2! Use the standards to calculate ammonium in samples
# Now BF-corrected FLU of standard curve against calculated NH4 concentration
mod2_kc22 <- lm(mean_FLU ~ true_nh4_conc, data = standard_b_kc22)
summary(mod2_kc22)

# pull coefficients from model
int2_kc22 <- coef(mod2_kc22)[1]
slope2_kc22 <- coef(mod2_kc22)[2]

# calculate how much ammonium is in the samples
samples_kc22 <- data_kc22 %>%
  filter(sample != "standard") %>%
  mutate(nh4_conc = (mean_FLU - int2_kc22)/slope2_kc22)

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
