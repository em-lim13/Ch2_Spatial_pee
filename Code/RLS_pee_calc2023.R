# CODE to calculate ammonium in 2023 RLS samples
# May 24, 2023
# Em Lim

# Load packages -----
library(tidyverse)
library(visreg)

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
write_csv(bottles2023, "Output/Output_data/RLS_nh4_2023")
