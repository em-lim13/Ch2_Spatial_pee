# # Code to calculate ammonium in 2021 Spring and Summer RLS samples
# June 2, 2023
# Em Lim

# Load packages -----
library(tidyverse)
library(visreg)
library(ggplot2)

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
write_csv(bottles2021, "Output/Output_data/RLS_nh4_2021")
