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

rls2023 <- do.call(rbind, mget(names(dfs)[dfs])) %>%
  filter(nh4_conc < 3) # remove the rave bottle

# summary stats ----
rls_avg <- rls2023 %>%
  summarise(avg = mean(nh4_conc),
            max = max(nh4_conc))


# Exploring plots ----

ggplot(rls2023, aes(site, nh4_conc)) +
  geom_boxplot()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

