# # Code to calculate ammonium in 2021 Spring and Summer RLS samples
# June 2, 2023
# Em Lim

# Load packages -----
library(tidyverse)
library(visreg)
library(ggplot2)

# 2021 -------
# Calculate NH4+ for April - May 2021 RLS samples -------
# Load bottle data
# Each row is an individual bottle

bottles1 <- read_csv("Data/RLS/2021/2021_05_11_RLS_NH4_bottles.csv")
bottles2 <- read_csv("Data/RLS/2021/2021_05_11_RLS_NH4_bottles2.csv")

# load fluorometry data
glow1 <- read_csv("Data/RLS/2021/2021_05_18_fluorometry.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

glow2 <- read_csv("Data/RLS/2021/2021_05_21_fluorometry.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

# load standard curve data
standard1 <- read_csv("Data/RLS/2021/2021_05_18_standard_curve.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

standard2 <- read_csv("Data/RLS/2021/2021_05_21_standard_curve.csv") %>%
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
    int = int1, #include values for the int and slope in for every column
    slope = slope1, #to calculate the conversion to NH4 conc
    matrix = 100 * ((Fst_spike1 - Fst_zero1 - (Fsm_spike - Fsm_zero))/
                      (Fst_spike1 - Fst_zero1)),
    matrix_est = 7.82,
    Fsm_cor_Holmes = Fsm_zero + (Fsm_zero * (matrix_est/100)),
    Fsm_cor_Taylor = (Fsm_zero / (1- matrix_est/100)),
    Fst_spike_cor = (Fst_spike / (1- matrix_est/100)), # correct for matrix
    nh4_conc = ifelse(sample == "sample", (int1 + slope1 * Fsm_cor_Taylor),
                           (int1 + slope1 * (Fsm_cor_Taylor - Fst_spike_cor))
  ))
# So I didn't use the matrix spikes to calculate the matrix effects, I just estimated them at 7.82% and used that for the 2021 bottles
# And then calculated what the matrix bottles would be if we hadn't spiked them by subtracting the 200 standard FLU from that bottle

#repeat for 2
Fst_zero2 <- standard_2f$mean_FLU[standard_2f$nh4_vol_uL == "0"]
Fst_spike2 <- standard_2f$mean_FLU[standard_2f$nh4_vol_uL == "200"]

bottles2_b <- bottles2_a %>% 
  left_join(bottles2_spike, by = "site_ID") %>%
  mutate(
    Fst_zero = Fst_zero2,
    Fst_spike = Fst_spike2,
    int = int2, #include values for the int and slope in for every column
    slope = slope2, 
    matrix = 100 * ((Fst_spike2 - Fst_zero2 - (Fsm_spike - Fsm_zero))/
                      (Fst_spike2 - Fst_zero2)),
    # Set estimated matrix = 7.82 for the sites where I messed up their matrix
    # But for DixonSW and DixonInside I did it right so keep those
    matrix_est = case_when(
      site == "DixonSW" & syringe == "WP1" & sample == "sample" ~ matrix, 
      site == "DixonInside" & syringe == "WP2" & sample == "sample" ~ matrix, 
      TRUE ~ as.numeric(7.82)),
    Fsm_cor_Holmes = Fsm_zero + (Fsm_zero * (matrix_est/100)),
    Fsm_cor_Taylor = (Fsm_zero / (1- matrix_est/100)),
    Fst_spike_cor = (Fst_spike / (1- matrix_est/100)), # correct for matrix
    nh4_conc = ifelse(sample == "sample", (int2 + slope2 * Fsm_cor_Taylor),
                           (int2 + slope2 * (Fsm_cor_Taylor - Fst_spike_cor))),
    cut_out = case_when(
      site == "DixonSW" & syringe == "WP1" & sample == "matrix" ~ "cut", 
      site == "DixonInside" & syringe == "WP2" & sample == "matrix" ~ "cut", 
      TRUE ~ as.character("keep"))) %>%
  filter(cut_out == "keep") %>% # manually filter out the two Dixon matrix bottles
  dplyr::select(-cut_out) # hide my manual shame

# merge two dfs
bottles_f_may2021 <- rbind(bottles1_b, bottles2_b) %>%
  mutate(site = factor(site, levels = c("WouwerChannel", "EffinghamW", "RaymondKelpRock", "EffinghamArchipelago", "FaberS", "EussenRock", "BaeriaSouthN", "BaeriaNorthS", "BaeriaNorthN", "Ross", "WizardN", "WizardS", "Ohiat", "Kirby", "EdKingE", "EdKingSWPyramid", "Taylor", "DodgerChannel", "Kiixin", "Scotts", "DixonSW", "DixonInside")),
         month = "May",
         year = "2021") %>%
  select(-dilution)

#write_csv(bottles_f_may2021, "Output/Output_data/RSL_fluo_data.csv")


# Calculate NH4+ for the June 2021 RSL surveys --------
bottles_june <- read_csv("Data/RLS/2021/2021_06_23_RLS_NH4_bottles.csv")

# load fluorometry data
glow_june <- read_csv("Data/RLS/2021/2021_06_23_RSL_fluorometry.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

# load standard curve data
standard_june <- read_csv("Data/RLS/2021/2021_06_23_standard_curve_data.csv") %>%
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
    matrix_est = 7.82) %>%
  select(site_ID, Fsm_spike, Fst_zero, Fst_spike, matrix, matrix_est)

bottles_june_2021 <- bottles_june_depth %>%
  left_join(bottles_june_matrix, by = "site_ID") %>%
  mutate(
    Fsm_cor_Holmes = Fsm_zero + (Fsm_zero * (matrix/100)),
    Fsm_cor_Taylor = (Fsm_zero / (1- matrix/100))
  ) %>%
  mutate(int = int_june, #include values for the int and slope in for every column
         slope = slope_june) %>% #to calculate the coversion to NH4 conc
  #those values come from our standard curve
  mutate(nh4_conc = int + slope * Fsm_cor_Taylor) %>%
  filter(sample != "matrix") %>%
  select(-sample_matrix) %>%
  mutate(month = "June",
         year = "2021")

ggplot(bottles_june_2021, aes(site, nh4_conc)) + geom_boxplot()

# ALl surface samples are negative

# Calculate NH4+ for the July 2021 RSL surveys --------
bottles_july <- read_csv("Data/RLS/2021/2021_07_20_RLS_NH4_bottles.csv")

# load fluorometry data
glow_july <- read_csv("Data/RLS/2021/2021_07_20_RSL_fluorometry.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(flu_1, flu_2, flu_3)))

# load standard curve data
standard_july <- read_csv("Data/RLS/2021/2021_07_20_standard.csv") %>%
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
    matrix_est = 7.82) %>%
  select(site_ID, Fsm_spike, Fst_zero, Fst_spike, matrix, matrix_est)

bottles_july_2021 <- bottles_july2 %>%
  left_join(bottles_july_matrix, by = "site_ID") %>%
  mutate(
    Fsm_cor_Holmes = Fsm_zero + (Fsm_zero * (matrix/100)),
    Fsm_cor_Taylor = (Fsm_zero / (1- matrix/100))
  ) %>%
  mutate(int = int_july, #include values for the int and slope in for every column
         slope = slope_july) %>% #to calculate the coversion to NH4 conc
  #those values come from our standard curve
  mutate(nh4_conc = int + slope * Fsm_cor_Taylor) %>%
  filter(sample != "matrix") %>%
  select(-sample_matrix) %>%
  mutate(month = "July",
         year = "2021")

ggplot(bottles_july_2021, aes(site, nh4_conc)) + geom_boxplot()

# Merge the 2021 samples -----
# First make the samples comparable
bottles_july_2021_merge <- bottles_july_2021 %>%
  select(site, site_ID, depth, temp_est, sal_est, date, month, year, matrix, nh4_conc)

bottles_june_2021_merge <- bottles_june_2021 %>%
  select(site, site_ID, depth, temp_est, sal_est, date, month, year, matrix, nh4_conc)
#filter(site != "WizardN") %>% # get rid of the weird negative ammonium Wizard site as it's probably a mistake

bottles_f_may2021_merge <- bottles_f_may2021 %>%
  select(site, site_ID, depth, temp_est, sal_est, date, month, year, matrix, nh4_conc)

# Now make one big dataframe
bottles2021 <- rbind(bottles_july_2021_merge, bottles_june_2021_merge, bottles_f_may2021_merge) %>%
  mutate(site = case_when(site_ID == "BMSC1" ~ "DodgerChannel",
                          site_ID == "BMSC2" ~ "Kirby",
                          site_ID == "BMSC3" ~ "Ohiat",
                          site_ID == "BMSC4" ~ "Kiixin",
                          site_ID == "BMSC5" ~ "Taylor",
                          site_ID == "BMSC6" ~ "BaeriaSouthN",
                          site_ID == "BMSC7" ~ "BaeriaNorthS",
                          site_ID == "BMSC8" ~ "BaeriaNorthN",
                          site_ID == "BMSC9" ~ "Scotts",
                          site_ID == "BMSC10" ~ "RossSlug",
                          site_ID == "BMSC11" ~ "WizardS",
                          site_ID == "BMSC12" ~ "WizardN",
                          site_ID == "BMSC13" ~ "EffinghamW",
                          site_ID == "BMSC14" ~ "EffinghamArchipelago",
                          site_ID == "BMSC15" ~ "RaymondKelpRock",
                          site_ID == "BMSC16" ~ "FaberS",
                          site_ID == "BMSC17" ~ "WouwerChannel",
                          site_ID == "BMSC18" ~ "EussenRock",
                          site_ID == "BMSC19" ~ "EdKingSWPyramid",
                          site_ID == "BMSC20" ~ "EdKingE",
                          site_ID == "BMSC21" ~ "DixonSW",
                          site_ID == "BMSC22" ~ "DixonInside",
                          site_ID == "BMSC23" ~ "AguilarPt",
                          site_ID == "BMSC24" ~ "SwissBoy",
                          site_ID == "BMSC25" ~ "GobyTown",
                          site_ID == "BMSC26" ~ "HosieSouth",
                          site_ID == "BMSC27" ~ "SanJoseNorth",
  ))


# Now write these files to a csv
write_csv(bottles2021, "Output/Output_data/RLS_nh4_2021.csv")


# 2022 ----
# first clear environment
rm(list=ls())

# # Calculate NH4+ for the April-May 2022 RLS samples -------
# Load bottle data
# Each row is an individual bottle
bottles3 <- read_csv("Data/RLS/2022/2022_05_07_RLS_NH4_bottles.csv") %>%
  unite("site_ID_depth", c(site_ID,survey_depth), sep= "-", 
        remove = FALSE) # give each site a unique code based on the site code and survey depth 

bottles4 <- read_csv("Data/RLS/2022/2022_05_07_RLS_NH4_bottles2.csv") %>%
  unite("site_ID_depth", c(site_ID,survey_depth), sep= "-", 
        remove = FALSE)

# load fluorometry data
#each row is a bottle
glow3 <- read_csv("Data/RLS/2022/2022_05_17_fluorometry.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

glow4 <- read_csv("Data/RLS/2022/2022_05_18_fluorometry.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

# load standard curve data
# one standard curve for each sampling period
standard3 <- read_csv("Data/RLS/2022/2022_05_17_RLS_standard.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

standard4 <- read_csv("Data/RLS/2022/2022_05_18_RLS_standard.csv") %>%
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
  mutate(adj_nh4_conc = nh4_conc + min_nh4_3, # shift all the readings up using the lowest as the 0
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

# merge all 2022 samples ----
bottles_may2022 <- rbind(bottles3_d, bottles4_d) %>%
  mutate(month = "May",
         year = "2022",
         nh4_conc = adj_nh4_conc) %>%
  select(site, site_ID, depth, temp_est, sal_est, date, month, year, matrix, nh4_conc) %>%
  mutate(site = case_when(site_ID == "BMSC1" ~ "DodgerChannel",
                          site_ID == "BMSC2" ~ "Kirby",
                          site_ID == "BMSC3" ~ "Ohiat",
                          site_ID == "BMSC4" ~ "Kiixin",
                          site_ID == "BMSC5" ~ "Taylor",
                          site_ID == "BMSC6" ~ "BaeriaSouthN",
                          site_ID == "BMSC7" ~ "BaeriaNorthS",
                          site_ID == "BMSC8" ~ "BaeriaNorthN",
                          site_ID == "BMSC9" ~ "Scotts",
                          site_ID == "BMSC10" ~ "RossSlug",
                          site_ID == "BMSC11" ~ "WizardS",
                          site_ID == "BMSC12" ~ "WizardN",
                          site_ID == "BMSC13" ~ "EffinghamW",
                          site_ID == "BMSC14" ~ "EffinghamArchipelago",
                          site_ID == "BMSC15" ~ "RaymondKelpRock",
                          site_ID == "BMSC16" ~ "FaberS",
                          site_ID == "BMSC17" ~ "WouwerChannel",
                          site_ID == "BMSC18" ~ "EussenRock",
                          site_ID == "BMSC19" ~ "EdKingSWPyramid",
                          site_ID == "BMSC20" ~ "EdKingE",
                          site_ID == "BMSC21" ~ "DixonSW",
                          site_ID == "BMSC22" ~ "DixonInside",
                          site_ID == "BMSC23" ~ "AguilarPt",
                          site_ID == "BMSC24" ~ "SwissBoy",
                          site_ID == "BMSC25" ~ "GobyTown",
                          site_ID == "BMSC26" ~ "HosieSouth",
                          site_ID == "BMSC27" ~ "SanJoseNorth",
  ))


# Now write these files to a csv
write_csv(bottles_may2022, "Output/Output_data/RLS_nh4_2022.csv")

# 2023 -----
# first clear environment
rm(list=ls())


# Source super cool function
source("Code/Functions.R")

# Now read allllllll the site data and let the magic function do all the work
samples_bmsc21 <- pee_calc("Data/RLS/2023/2023_05_08_BMSC21.csv")
samples_bmsc22 <- pee_calc("Data/RLS/2023/2023_05_08_BMSC22_b.csv")# drop 600 spike
samples_bmsc10 <- pee_calc("Data/RLS/2023/2023_05_09_BMSC10.csv")
samples_bmsc2 <- pee_calc("Data/RLS/2023/2023_05_09_BMSC2_b.csv")# drop 600 spike
samples_bmsc19 <- pee_calc("Data/RLS/2023/2023_05_10_BMSC19.csv")
samples_bmsc20 <- pee_calc("Data/RLS/2023/2023_05_10_BMSC20_b.csv")# drop 600 spike

# Samples read May 26
samples_bmsc9 <- pee_calc("Data/RLS/2023/2023_05_11_BMSC9.csv")
samples_bmsc23 <- pee_calc("Data/RLS/2023/2023_05_11_BMSC23.csv")
samples_bmsc4 <- pee_calc("Data/RLS/2023/2023_05_12_BMSC4.csv")
samples_bmsc3 <- pee_calc("Data/RLS/2023/2023_05_12_BMSC3.csv")
samples_bmsc1 <- pee_calc("Data/RLS/2023/2023_05_15_BMSC1.csv")
samples_bmsc5 <- pee_calc("Data/RLS/2023/2023_05_15_BMSC5.csv")

# Samples read May 27, 2023
samples_bmsc11 <- pee_calc("Data/RLS/2023/2023_05_16_BMSC11.csv")
samples_bmsc12 <- pee_calc("Data/RLS/2023/2023_05_16_BMSC12b.csv")
# used the BF from BMSC 11 because the bf sample was opa spiked earlier (for the bad one) and so when it was read the second time the bf wasn't a bf anymore

# Samples read May 30, 2023
samples_bmsc6 <- pee_calc("Data/RLS/2023/2023_05_17_BMSC6.csv")
samples_bmsc8 <- pee_calc("Data/RLS/2023/2023_05_17_BMSC8.csv")
samples_bmsc24 <- pee_calc("Data/RLS/2023/2023_05_18_BMSC24.csv")
samples_bmsc25 <- pee_calc("Data/RLS/2023/2023_05_18_BMSC25.csv")
samples_bmsc26 <- pee_calc("Data/RLS/2023/2023_05_19_BMSC26.csv")
samples_bmsc27 <- pee_calc("Data/RLS/2023/2023_05_19_BMSC27.csv")

# merge all samples
dfs = sapply(.GlobalEnv, is.data.frame) 

bottles2023 <- do.call(rbind, mget(names(dfs)[dfs])) %>%
  filter(nh4_conc < 3) %>% # remove the rave bottle 
  mutate(site_ID = site_code,
         temp_est = NA,
         sal_est = NA,
         date = collection_date,
         month = "May",
         year = "2023",
         matrix = NA) %>%
  select(site, site_ID, depth, temp_est, sal_est, date, month, year, matrix, nh4_conc)%>%
  mutate(site = case_when(site_ID == "BMSC1" ~ "DodgerChannel",
                          site_ID == "BMSC2" ~ "Kirby",
                          site_ID == "BMSC3" ~ "Ohiat",
                          site_ID == "BMSC4" ~ "Kiixin",
                          site_ID == "BMSC5" ~ "Taylor",
                          site_ID == "BMSC6" ~ "BaeriaSouthN",
                          site_ID == "BMSC7" ~ "BaeriaNorthS",
                          site_ID == "BMSC8" ~ "BaeriaNorthN",
                          site_ID == "BMSC9" ~ "Scotts",
                          site_ID == "BMSC10" ~ "RossSlug",
                          site_ID == "BMSC11" ~ "WizardS",
                          site_ID == "BMSC12" ~ "WizardN",
                          site_ID == "BMSC13" ~ "EffinghamW",
                          site_ID == "BMSC14" ~ "EffinghamArchipelago",
                          site_ID == "BMSC15" ~ "RaymondKelpRock",
                          site_ID == "BMSC16" ~ "FaberS",
                          site_ID == "BMSC17" ~ "WouwerChannel",
                          site_ID == "BMSC18" ~ "EussenRock",
                          site_ID == "BMSC19" ~ "EdKingSWPyramid",
                          site_ID == "BMSC20" ~ "EdKingE",
                          site_ID == "BMSC21" ~ "DixonSW",
                          site_ID == "BMSC22" ~ "DixonInside",
                          site_ID == "BMSC23" ~ "AguilarPt",
                          site_ID == "BMSC24" ~ "SwissBoy",
                          site_ID == "BMSC25" ~ "GobyTown",
                          site_ID == "BMSC26" ~ "HosieSouth",
                          site_ID == "BMSC27" ~ "SanJoseNorth",
  ))


rownames(bottles2023) <- NULL

# Now write these files to a csv
write_csv(bottles2023, "Output/Output_data/RLS_nh4_2023.csv")
