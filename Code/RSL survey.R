# Code to analyze RLS pee samples
# May 11, 2021
# By: Em Lim

# Load packages and data ----

# Load packages
library(tidyverse)
library(visreg)
library(ggplot2)
theme_set(theme_bw())
library(RColorBrewer)
library(ggeffects)
source("Code/theme_black.R")
     
# Calculate NH4+ for April - May 2021 RLS samples -------
# Load bottle data
# Each row is an individual bottle

bottles1 <- read_csv("Data/RLS/2021_05_11_RLS_NH4_bottles.csv")
bottles2 <- read_csv("Data/RLS/2021_05_11_RLS_NH4_bottles2.csv")

# load fluorometry data
glow1 <- read_csv("Data/RLS/2021_05_18_fluorometry.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

glow2 <- read_csv("Data/RLS/2021_05_21_fluorometry.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

# load standard curve data
standard1 <- read_csv("Data/RLS/2021_05_18_standard_curve.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

standard2 <- read_csv("Data/RLS/2021_05_21_standard_curve.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

# # STANDARD CURVE-
#do the calculations for the standard curve
standard_1f <- standard1 %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * nh4_conc_og_umol, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + og_vol_L, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L, #concentration
         #of NH4 in seawater sample
         mean_FLU = (FLU1 + FLU2 + FLU3)/3) #mean FLU reading

standard_2f <- standard2 %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * nh4_conc_og_umol, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + og_vol_L, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L, #concentration
         #of NH4 in seawater sample
         mean_FLU = (FLU1 + FLU2 + FLU3)/3) #mean FLU reading


#linear model between the fluorometer reading and actual concentration of NH4
sc_mod1 <- lm(nh4_conc_final_umol_L ~ mean_FLU, data = standard_1f)
#check this summary for the R squared- the relationship should be TIGHT, so if it's lower than like 0.98, there might be some contamination
summary(sc_mod1)
#visualize curve to make sure it looks right
ggplot(standard_1f, aes(mean_FLU, nh4_conc_final_umol_L)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#extract intercept and slope to use later to determine sample concentrations
#based on fluorometry readings
int1 <- coef(sc_mod1)[1]
slope1 <- coef(sc_mod1)[2]

#linear model between the fluorometer reading and actual concentration of NH4
sc_mod2 <- lm(nh4_conc_final_umol_L ~ mean_FLU, data = standard_2f)
#check this summary for the R squared- the relationship should be TIGHT, so if it's lower than like 0.98, there might be some contamination
summary(sc_mod2)
#visualize curve to make sure it looks right
ggplot(standard_2f, aes(mean_FLU, nh4_conc_final_umol_L)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#extract intercept and slope to use later to determine sample concentrations
#based on fluorometry readings
int2 <- coef(sc_mod2)[1]
slope2 <- coef(sc_mod2)[2]

# Calculate matrix effects -

bottles1_a <- bottles1 %>%
  left_join(glow1, by = c("bottle", "site_ID")) %>%
  mutate(Fsm_zero = mean_FLU)

bottles2_a <- bottles2 %>%
  left_join(glow2, by = c("bottle", "site_ID")) %>%
  mutate(Fsm_zero = mean_FLU)

bottles1_spike <- bottles1_a %>%
  filter(sample == "matrix") %>%
  transmute(site_ID = site_ID,
            Fsm_spike = mean_FLU)

bottles2_spike <- bottles2_a %>%
  filter(sample == "matrix") %>%
  transmute(site_ID = site_ID,
            Fsm_spike = mean_FLU)

#idk
Fst_zero1 <- standard_1f$mean_FLU[standard_1f$nh4_vol_uL == "0"]
Fst_spike1 <- standard_1f$mean_FLU[standard_1f$nh4_vol_uL == "200"]

bottles1_b <- bottles1_a %>% 
  left_join(bottles1_spike, by = "site_ID") %>%
  mutate(
    Fst_zero = Fst_zero1,
    Fst_spike = Fst_spike1,
    matrix = 100 * ((Fst_spike1 - Fst_zero1 - (Fsm_spike - Fsm_zero))/
                      (Fst_spike1 - Fst_zero1)),
    matrix_est = 15,
    Fsm_cor_Holmes = Fsm_zero + (Fsm_zero * (matrix_est/100)),
    Fsm_cor_Taylor = (Fsm_zero / (1- matrix_est/100))
  ) %>%
  mutate(int = int1, #include values for the int and slope in for every column
         slope = slope1) %>% #to calculate the coversion to NH4 conc
  #those values come from our standard curve
  mutate(nh4_conc = ifelse(sample == "sample", (int1 + slope1 * Fsm_cor_Taylor),
                           (int1 + slope1 * (Fsm_cor_Taylor - Fst_spike1))
  ))

#repeat for 2
Fst_zero2 <- standard_2f$mean_FLU[standard_2f$nh4_vol_uL == "0"]
Fst_spike2 <- standard_2f$mean_FLU[standard_2f$nh4_vol_uL == "200"]

bottles2_b <- bottles2_a %>% 
  left_join(bottles2_spike, by = "site_ID") %>%
  mutate(
    Fst_zero = Fst_zero2,
    Fst_spike = Fst_spike2,
    matrix = 100 * ((Fst_spike2 - Fst_zero2 - (Fsm_spike - Fsm_zero))/
                      (Fst_spike2 - Fst_zero2)),
    matrix_est = 15,
    Fsm_cor_Holmes = Fsm_zero + (Fsm_zero * (matrix_est/100)),
    Fsm_cor_Taylor = (Fsm_zero / (1- matrix_est/100))
  ) %>%
  mutate(int = int2, #include values for the int and slope in for every column
         slope = slope2) %>% #to calculate the coversion to NH4 conc
  #those values come from our standard curve
  mutate(nh4_conc = ifelse(sample == "sample", (int2 + slope2 * Fsm_cor_Taylor),
                           (int2 + slope2 * (Fsm_cor_Taylor - Fst_spike2))
  ))


# merge two dfs
bottles_f <- rbind(bottles1_b, bottles2_b) %>%
  mutate(site = factor(site, levels = c("WouwerChannel", "EffinghamW", "RaymondKelpRock", "EffinghamArchipelago", "FaberS", "EussenRock", "BaeriaSouthN", "BaeriaNorthS", "BaeriaNorthN", "Ross", "WizardN", "WizardS", "Ohiat", "Kirby", "EdKingE", "EdKingSWPyramid", "Taylor", "DodgerChannel", "Kiixin", "Scotts", "DixonSW", "DixonInside"))) %>%
  select(-dilution)

#write_csv(bottles_f, "Output/Output_data/RSL_fluo_data.csv")


# Calculate NH4+ for the June 2021 RSL surveys --------
bottles_june <- read_csv("Data/RLS/2021_06_23_RLS_NH4_bottles.csv")

# load fluorometry data
glow_june <- read_csv("Data/RLS/2021_06_23_RSL_fluorometry.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

# load standard curve data
standard_june <- read_csv("Data/RLS/2021_06_23_standard_curve_data.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

#do the calculations for the standard curve
standard_june_f <- standard_june %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * nh4_conc_og_umol, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + og_vol_L, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L, #concentration
         #of NH4 in seawater sample
         mean_FLU = (FLU1 + FLU2 + FLU3)/3) #mean FLU reading

#graph standard curve
ggplot(standard_june_f, aes(mean_FLU, nh4_conc_final_umol_L)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#model
sc_mod_june <- lm(nh4_conc_final_umol_L ~ mean_FLU, data = standard_june_f)
summary(sc_mod_june)

#save coefficients
int_june <- coef(sc_mod_june)[1]
slope_june <- coef(sc_mod_june)[2]

# calculate matrix
bottles_june_depth <- bottles_june %>%
  left_join(glow_june, by = c("bottle")) %>%
  mutate(Fsm_zero = mean_FLU) %>%
  filter(depth > 0)

bottles_june_spike <- bottles_june_depth %>%
  filter(sample == "matrix") %>%
  transmute(site_ID = site_ID,
            Fsm_spike = mean_FLU)

bottles_june_zero <- bottles_june_depth %>%
  filter(syringe == "WP") %>%
  filter(sample != "matrix")

Fst_zero_june <- standard_june_f$mean_FLU[standard_june_f$nh4_vol_uL == "0"]
Fst_spike_june <- standard_june_f$mean_FLU[standard_june_f$nh4_vol_uL == "200"]

bottles_june_matrix <- bottles_june_zero %>% 
  left_join(bottles_june_spike, by = "site_ID")%>%
  mutate(
    Fst_zero = Fst_zero_june,
    Fst_spike = Fst_spike_june,
    matrix = 100 * ((Fst_spike_june - Fst_zero_june - (Fsm_spike - Fsm_zero))/
                      (Fst_spike_june - Fst_zero_june)),
    matrix_est = 15) %>%
  select(site_ID, Fsm_spike, Fst_zero, Fst_spike, matrix, matrix_est)

bottles_june_depth2 <- bottles_june_depth %>%
  left_join(bottles_june_matrix, by = "site_ID") %>%
  mutate(
    Fsm_cor_Holmes = Fsm_zero + (Fsm_zero * (matrix/100)),
    Fsm_cor_Taylor = (Fsm_zero / (1- matrix/100))
  ) %>%
  mutate(int = int_june, #include values for the int and slope in for every column
         slope = slope_june) %>% #to calculate the coversion to NH4 conc
  #those values come from our standard curve
  mutate(nh4_conc = int_june + slope_june * Fsm_cor_Taylor) %>%
  filter(sample != "matrix") %>%
  select(-sample_matrix)



ggplot(bottles_june_depth2, aes(site, nh4_conc)) + geom_boxplot()

# ALl surface samples are negative

# Calculate NH4+ for the July 2021 RSL surveys --------
bottles_july <- read_csv("Data/RLS/2021_07_20_RLS_NH4_bottles.csv")

# load fluorometry data
glow_july <- read_csv("Data/RLS/2021_07_20_RSL_fluorometry.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(flu_1, flu_2, flu_3)))

# load standard curve data
standard_july <- read_csv("Data/RLS/2021_07_20_standard.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

#do the calculations for the standard curve
standard_july_f <- standard_july %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * nh4_conc_og_umol, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + og_vol_L, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L, #concentration
         #of NH4 in seawater sample
         mean_FLU = (FLU1 + FLU2 + FLU3)/3) #mean FLU reading

#graph standard curve
ggplot(standard_july_f, aes(mean_FLU, nh4_conc_final_umol_L)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#model
sc_mod_july <- lm(nh4_conc_final_umol_L ~ mean_FLU, data = standard_july_f)
summary(sc_mod_july)

#save coefficients
int_july <- coef(sc_mod_july)[1]
slope_july <- coef(sc_mod_july)[2]

# calculate matrix
bottles_july2 <- bottles_july %>%
  left_join(glow_july, by = c("bottle")) %>%
  mutate(Fsm_zero = mean_FLU)

bottles_july_spike <- bottles_july2 %>%
  filter(sample == "matrix") %>%
  transmute(site_ID = site_ID,
            Fsm_spike = mean_FLU)

bottles_july_zero <- bottles_july2 %>%
  filter(syringe == "WP") %>%
  filter(sample != "matrix")

Fst_zero_july <- standard_july_f$mean_FLU[standard_july_f$nh4_vol_uL == "0"]
Fst_spike_july <- standard_july_f$mean_FLU[standard_july_f$nh4_vol_uL == "200"]

bottles_july_matrix <- bottles_july_zero %>% 
  left_join(bottles_july_spike, by = "site_ID")%>%
  mutate(
    Fst_zero = Fst_zero_july,
    Fst_spike = Fst_spike_july,
    matrix = 100 * ((Fst_spike_july - Fst_zero_july - (Fsm_spike - Fsm_zero))/
                      (Fst_spike_july - Fst_zero_july)),
    matrix_est = 15) %>%
  select(site_ID, Fsm_spike, Fst_zero, Fst_spike, matrix, matrix_est)

bottles_july3 <- bottles_july2 %>%
  left_join(bottles_july_matrix, by = "site_ID") %>%
  mutate(
    Fsm_cor_Holmes = Fsm_zero + (Fsm_zero * (matrix/100)),
    Fsm_cor_Taylor = (Fsm_zero / (1- matrix/100))
  ) %>%
  mutate(int = int_july, #include values for the int and slope in for every column
         slope = slope_july) %>% #to calculate the coversion to NH4 conc
  #those values come from our standard curve
  mutate(nh4_conc = int_july + slope_july * Fsm_cor_Taylor) %>%
  filter(sample != "matrix") %>%
  select(-sample_matrix)

ggplot(bottles_july3, aes(site, nh4_conc)) + geom_boxplot()

# # Calculate NH4+ for the April -May 2022 RLS samples -------
# Load bottle data
# Each row is an individual bottle
bottles3 <- read_csv("Data/RLS/2022_05_07_RLS_NH4_bottles.csv") %>%
  unite("site_ID_depth", c(site_ID,survey_depth), sep= "-", 
        remove = FALSE) # give each site a unique code based on the site code and survey depth 

bottles4 <- read_csv("Data/RLS/2022_05_07_RLS_NH4_bottles2.csv") %>%
  unite("site_ID_depth", c(site_ID,survey_depth), sep= "-", 
        remove = FALSE)

# load fluorometry data
#each row is a bottle
glow3 <- read_csv("Data/RLS/2022_05_17_fluorometry.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

glow4 <- read_csv("Data/RLS/2022_05_18_fluorometry.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

# load standard curve data
# one standard curve for each sampling period
standard3 <- read_csv("Data/RLS/2022_05_17_RLS_standard.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

standard4 <- read_csv("Data/RLS/2022_05_18_RLS_standard.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

#do the calculations for the standard curve
standard3_f <- standard3 %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * nh4_conc_og_umol, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + og_vol_L, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L, #concentration
         #of NH4 in seawater sample
         mean_FLU = (FLU1 + FLU2 + FLU3)/3,
         nh4_conc_adj = nh4_conc_final_umol_L + 0.3622925) # adjust the nh4 conc in the standard bottles because the DW had 0.3622925 uM NH4+ in it
# COME BACK HERE TO UPDATE SFU DW

standard4_f <- standard4 %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * nh4_conc_og_umol, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + og_vol_L, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L, #concentration
         #of NH4 in seawater sample
         mean_FLU = (FLU1 + FLU2 + FLU3)/3,
         nh4_conc_adj = nh4_conc_final_umol_L + 0.3622925) # adjust the nh4 conc in the standard bottles because the DW had 0.3622925 uM NH4+ in it
# COME BACK HERE TO UPDATE SFU DW 

#graph standard curves to visually inspect
ggplot(standard3_f, aes(mean_FLU, nh4_conc_adj)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

ggplot(standard4_f, aes(mean_FLU, nh4_conc_adj)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#model
sc_mod3 <- lm(nh4_conc_adj ~ mean_FLU, data = standard3_f)
summary(sc_mod3)

sc_mod4 <- lm(nh4_conc_adj ~ mean_FLU, data = standard4_f)
summary(sc_mod4)

#save coefficients (intercept and slope for NH4+ calculation later)
int3 <- coef(sc_mod3)[1]
slope3 <- coef(sc_mod3)[2]

int4 <- coef(sc_mod4)[1]
slope4 <- coef(sc_mod4)[2]

# Calculate matrix effects for first round of fluorometry-------

# calculate matrix as according to Taylor et al 2007
# % Matrix Effects = [(standard_spike - standard_zero) - (sample_spike - sample_zero)]/ (standard_spike - standard_zero) * 100

# join the bottle data with the fluorometer readings
bottles3_b <- bottles3 %>%
  left_join(glow3, by = c("bottle", "site_ID")) %>%
  mutate(Fsm_zero = mean_FLU)

# create a df of just the spike bottles, keep the site_ID_depth so we can join them later
bottles3_spike <- bottles3_b %>%
  filter(sample == "matrix") %>%
  transmute(site_ID_depth = site_ID_depth,
            Fsm_spike = mean_FLU)

# create a df of just the zero bottles (the non spiked bottle from the whirl paks)
bottles3_zero <- bottles3_b %>%
  filter(syringe == "WP") %>%
  filter(sample != "matrix")

# save the values of the 0 and 200 standards which we'll compare the whirl pak zero and spike to
Fst_zero3 <- standard3_f$mean_FLU[standard3_f$nh4_vol_uL == "0"]
Fst_spike3 <- standard3_f$mean_FLU[standard3_f$nh4_vol_uL == "200"]

# join the zero and spike bottles from the whirl pack back together by site_ID_depth
# then use the zero and spike samples + zero and spike standards to calculate the matrix for each sampling unit (site_ID_depth)
bottles3_matrix <- bottles3_zero %>% 
  left_join(bottles3_spike, by = "site_ID_depth") %>%
  mutate(
    Fst_zero = Fst_zero3,
    Fst_spike = Fst_spike3,
    matrix = 100 * ((Fst_spike - Fst_zero - (Fsm_spike - Fsm_zero))/
                      (Fst_spike - Fst_zero)),
    matrix_est = 20) %>%
  select(site_ID_depth, Fsm_spike, Fst_zero, Fst_spike, matrix, matrix_est) # just keep the relevant info for the matrix and the site_ID_depth to join this back up with the fluorometer readings for the rest of the samples

# now join this matrix info back up with all the samples, and use it to calculate the nh4 concentration
bottles3_c <- bottles3_b %>%
  left_join(bottles3_matrix, by = "site_ID_depth") %>%
  mutate(
    Fsm_cor_Holmes = Fsm_zero + (Fsm_zero * (matrix/100)),
    Fsm_cor_Taylor = (Fsm_zero / (1- matrix/100))
    ) %>%
  mutate(int = int3, #include values for the int and slope in for every column
         slope = slope3) %>% #to calculate the coversion to NH4 conc
  #those values come from our standard curve
  mutate(nh4_conc = int + slope * Fsm_cor_Taylor) %>%
  filter(sample != "matrix") %>%
  select(-sample)

# all of the samples are negative! this is terrible!
# what happens if we set the lowest nh4 value = 0, and add that to all of them???
min_nh4_3 <- - min(bottles3_c$nh4_conc) # make it positive so I can add it

bottles3_d <- bottles3_c %>%
  mutate(adj_nh4_conc = nh4_conc + min_nh4_3, # shift all the readings up using hte lowest as the 0
         fluo_date = "May_17")


# Calculate matrix effects for second round of fluorometry-------

# calculate matrix as according to Taylor et al 2007
# % Matrix Effects = [(standard_spike - standard_zero) - (sample_spike - sample_zero)]/ (standard_spike - standard_zero) * 100

# join the bottle data with the fluorometer readings
bottles4_b <- bottles4 %>%
  left_join(glow4, by = c("bottle", "site_ID")) %>%
  mutate(Fsm_zero = mean_FLU)

# create a df of just the spike bottles, keep the site_ID_depth so we can join them latet
bottles4_spike <- bottles4_b %>%
  filter(sample == "matrix") %>%
  transmute(site_ID_depth = site_ID_depth,
            Fsm_spike = mean_FLU)

# create a df of just the zero bottles (the non spiked bottle from the whirlpak)
bottles4_zero <- bottles4_b %>%
  filter(syringe == "WP") %>%
  filter(sample != "matrix")

# save the values of the 0 and 200 standards
Fst_zero4 <- standard4_f$mean_FLU[standard4_f$nh4_vol_uL == "0"]
Fst_spike4 <- standard4_f$mean_FLU[standard4_f$nh4_vol_uL == "200"]

# join the zero and spike bottles from the whirlpack back together by site_ID_depth
# then use the zero and spike samples + zero and spike standards to calculate the matrix for each sampling unit (site_ID_depth)
bottles4_matrix <- bottles4_zero %>% 
  left_join(bottles4_spike, by = "site_ID_depth") %>%
  mutate(
    Fst_zero = Fst_zero4,
    Fst_spike = Fst_spike4,
    matrix = 100 * ((Fst_spike - Fst_zero - (Fsm_spike - Fsm_zero))/
                      (Fst_spike - Fst_zero)),
    matrix_est = 20) %>%
  select(site_ID_depth, Fsm_spike, Fst_zero, Fst_spike, matrix, matrix_est) # just keep the relevant info for the matrix and the site_ID_depth to join this back up with the fluorometer readings for the rest of the samples

# now join this matrix info back up with all the samples, and use it to calculate the nh4 concentration
bottles4_c <- bottles4_b %>%
  left_join(bottles4_matrix, by = "site_ID_depth") %>%
  mutate(
    Fsm_cor_Holmes = Fsm_zero + (Fsm_zero * (matrix/100)),
    Fsm_cor_Taylor = (Fsm_zero / (1- matrix/100))
  ) %>%
  mutate(int = int4, #include values for the int and slope in for every column
         slope = slope4) %>% #to calculate the coversion to NH4 conc
  #those values come from our standard curve
  mutate(nh4_conc = int + slope * Fsm_cor_Taylor) %>%
  filter(sample != "matrix") %>%
  select(-sample)

# all of the samples are negative! this is terrible!
# what happens if we set the lowest nh4 value = 0, and add that to all of them???
min_nh4_4 <- - min(bottles4_c$nh4_conc) # make it positive so I can add it

bottles4_d <- bottles4_c %>%
  mutate(adj_nh4_conc = nh4_conc + min_nh4_4,
         fluo_date = "May_18")

# merge all 2022 samples
bottles2022 <- rbind(bottles3_d, bottles4_d) %>%
  select(site, site_ID, depth, temp_est, sal_est, date, period, matrix, nh4_conc, adj_nh4_conc)

  
# Look at the 2022 samples
ggplot(bottles2022, aes(adj_nh4_conc, site, fill = site)) +
  geom_boxplot(colour = "white") +
  theme_black() +
  labs(x= "Ammonium concentration (umol/L)", y = "Site") +
  theme(legend.position = "none") +
  facet_wrap(~fluo_date)

# They're all pretty low.........

# Smoosh three together???? ------
july_bottles <- bottles_july3 %>%
  select(site, site_ID, depth, temp_est, sal_est, date, period, matrix, nh4_conc) %>%
  mutate(nh4_conc = ifelse(nh4_conc <0, 0, nh4_conc),
         adj_nh4_conc = nh4_conc)

june_bottles <- bottles_june_depth2 %>% 
  select(site, site_ID, depth, temp_est, sal_est, date, period, matrix, nh4_conc) %>%
  filter(site != "WizardN") %>% # get rid of the weird negative ammonium Wizard site as it's probably a mistake 
  mutate(nh4_conc = ifelse(nh4_conc <0, 0, nh4_conc),
         adj_nh4_conc = nh4_conc)

may_bottles <- bottles_f %>%
  select(site, site_ID, depth, temp_est, sal_est, date, period, matrix, nh4_conc) %>%
  mutate(nh4_conc = ifelse(nh4_conc <0, 0, nh4_conc),
         adj_nh4_conc = nh4_conc)

# join 2021 and 2022
rls_data = rbind(may_bottles, june_bottles, july_bottles, bottles2022) 

# Which sites have the lowest and highest pee
rls_pee_summary <- rls_data %>%
  group_by(site) %>%
  summarize(pee = mean(adj_nh4_conc)) %>%
  arrange(desc(pee))

#What is the tide height at those sites
rls_tide_summary <- rls_data %>%
  group_by(site) %>%
  summarize(depth = mean(depth)) %>%
  arrange(desc(depth))


# Correlation analysis
may_bottles_cor <- may_bottles %>%
  filter(site != "Kirby") %>%
  filter(site != "Kiixin") %>%
  filter(site != "BaeriaNorthN")%>%
  filter(site != "BaeriaNorthS") %>%
  slice(-56)%>%
  slice(-52)

bottles2022_cor <- bottles2022 %>%
  filter(site != "AguilarPoint") %>%
  slice(4:57)

cor_data <- may_bottles_cor %>%
  left_join(bottles2022_cor, by = "site")

res <- cor.test(may_bottles_cor$adj_nh4_conc, bottles2022_cor$adj_nh4_conc, 
                method = "spearman")
summary(res)

pee <- lm(adj_nh4_conc ~ site * period, data = rls_data)


# Stats -------
# Does pee vary by site?
simple_model <- lm(nh4_conc ~ site_ID + period, data = rls_data)
summary(simple_model)

sum_stats_pee <- ggpredict(simple_model, terms = c("site_ID", "period")) %>% 
  #and then we'll just rename one of the columns so it's easier to plot
  rename(site_ID = x,
         nh4_conc = predicted,
         period = group)
#View(sum_stats_crabs)

# Graphing ----

# Dot and whisker?
ggplot() +
  geom_point(data = sum_stats_pee, 
             aes(y = site_ID, x = nh4_conc, colour = period),
             size = 4) +
  geom_errorbar(data = sum_stats_pee, 
                aes(y = site_ID,
                    x = nh4_conc,
                    colour = period,
                    # and you can decide which type of error to show here
                    # we're using 95% CI
                    xmin = conf.low,
                    xmax = conf.high),
                width = 0.2,
                size = 1.2)  +
  geom_point(data = rls_data, aes (y = site_ID, x = nh4_conc, colour = period), alpha = 0.5, height = 0, size = 2) +
  labs(x= "Ammonium concentration (umol/L)", y = "Site") +
  theme_black() +
  theme(legend.position="none") 

#ggsave("Output/Figures/RLS_pee_black.png", device = "png",
#       height = 9, width = 16, dpi = 400)

# Site vs pee
ggplot(rls_data, aes(adj_nh4_conc, site, colour = period)) +
  geom_boxplot() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15)) +
  labs(x= "Ammonium concentration (umol)", y = "Site") +
  theme_classic()

# black background
ggplot(rls_data, aes(nh4_conc, site, fill = site)) +
  geom_boxplot(colour = "white") +
  theme_black() +
  labs(x= "Ammonium concentration (umol/L)", y = "Site") +
  theme(legend.position = "none") +
  xlim(c(0, 5))

#ggsave("Output/Figures/RLS_sites_pee.png", device = "png",
#       height = 9, width = 16, dpi = 400)

# Visualize salinity
ggplot(bottles_f, aes(x = sal_est, y = site)) +
  geom_point() +
  labs(x = "Salinity Estimate", y = "Site")

# histogram
ggplot(bottles, aes(x = sal_est)) +
  geom_histogram(colour="black", fill="white") +
  geom_vline(aes(xintercept = 31.8),
             color="blue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept = 27.7),
             color="red", linetype="dashed", size=1)

# Sal vs pee
ggplot(bottles_f, aes(sal_est, nh4_conc)) + geom_point() +
  geom_smooth(method = lm)

# Visualize temperature
ggplot(bottles, aes(x = temp_est, y = site)) +
  geom_point() +
  labs(x = "Salinity Estimate", y = "Site")

ggplot(bottles, aes(x = temp_est)) +
  geom_histogram(colour="black", fill="white")

# stats -------------------------

pee <- lm(nh4_conc ~ site * depth, bottles_f)
anova(pee)
summary(pee)

visreg(pee, "depth", by = "site")
