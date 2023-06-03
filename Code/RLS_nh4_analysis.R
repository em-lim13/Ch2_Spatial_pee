# Code to look at all the spatial data together
# June 2, 2023
# Em Lim

# Load packages and data ----
library(tidyverse)
library(visreg)
library(ggplot2)
theme_set(theme_bw())
library(PNWColors)
library(ggeffects)
source("Code/theme_black.R")

# Load data

# RLS nh4 from 2021
# This is the spring blitz, the june samples, and the July samples
# The may samples have estimated matrix effects, not great
# The June and July samples had a own matrix spike from each site but were still compared to DI
bottles2021 <- read_csv("Output/Output_data/RLS_nh4_2021")

# RLS from 2022
# These had a own matrix spike from each site but were still compared to DI
# They came out quite negative maybe because of SFU DI, and maybe because of temperature difference between the standards and samples
# I took the lowest nh4 reading and added it to everything else to "set" that sample to 0 and bump everything up
# Maybe an underestimation
bottles2022 <- read_csv("Output/Output_data/RLS_nh4_2022")

# RLS from 2023
# Did the full "proper" Taylor protocol with standard bottles + BF from each site
bottles2023 <- read_csv("Output/Output_data/RLS_nh4_2023")

# combine these three years into one!
rls_nh4 <- rbind(bottles2021, bottles2022, bottles2023) %>%
  mutate(year = as.factor(year)) %>%
  group_by(site) %>%
  mutate(nh4_avg = mean(nh4_conc)) %>%
  ungroup()

# Which sites have the lowest and highest pee
rls_pee_summary <- rls_nh4_all %>%
  group_by(site) %>%
  summarize(mean_nh4 = mean(nh4_conc)) %>%
  arrange(desc(mean_nh4)) 



# Plot these data??? ----
pal <- pnw_palette("Sailboat", 3)

ggplot(rls_nh4) +
  geom_point(aes(x = reorder(site, -nh4_avg), 
                 nh4_conc, colour = year, pch = month),
             size = 3, alpha = 0.75) +
  geom_point(aes(x = reorder(site, -nh4_avg), 
                 nh4_avg),
                 size = 3, colour = "black", pch = 20) +
  theme(axis.text.x = element_text(angle = 75, vjust = 0.5)) +
  labs(x = "Site", y = "NH4+ Concentration (umol/L)") +
  scale_colour_manual(values = pal)


# Summary stats -----

# What were the matrix effects like at the end of 2021 and in 2022?
rls_nh4 %>%
  filter(year == "2022") %>%
  summarise(matrix = mean(matrix))
# 14.6

rls_nh4 %>%
  filter(year == "2021") %>%
  filter(month != "May") %>%
  summarise(matrix = mean(matrix))
# 7.82
# Go back to the 2021 data and set the ME to 7.82

#What is the tide height at those sites
rls_tide_summary <- rls_data %>%
  group_by(site) %>%
  summarize(depth = mean(depth)) %>%
  arrange(desc(depth))


# Correlation analysis -----
# Make a DF where we only have overlap between sites
b_2021 <- rls_nh4 %>%
  filter(year == "2021") %>%
  transmute(site = site,
            nh4_2021 = nh4_avg) %>%
  unique()

b_2022 <- rls_nh4 %>%
  filter(year == "2022") %>%
  transmute(site = site,
            nh4_2022 = nh4_avg) %>%
  unique()

b_2023 <- rls_nh4 %>%
  filter(year == "2023") %>%
  transmute(site = site,
            nh4_2023 = nh4_avg) %>%
  unique()

# throw them together and drop the NAs
cor_data <- b_2021 %>%
  left_join(b_2022, by = "site") %>%
  left_join(b_2023, by = "site") %>%
  drop_na()

# OK so now do the correlations!

# 2021 vs 2022
cor.test(cor_data$nh4_2021, cor_data$nh4_2022, 
                method = "spearman")

# 2022 vs 2023
cor.test(cor_data$nh4_2022, cor_data$nh4_2023, 
         method = "spearman")

# 2023 vs 2021
cor.test(cor_data$nh4_2023, cor_data$nh4_2021, 
         method = "spearman")


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
ggplot(bottles_f_may2021, aes(x = sal_est, y = site)) +
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
ggplot(bottles_f_may2021, aes(sal_est, nh4_conc)) + geom_point() +
  geom_smooth(method = lm)

# Visualize temperature
ggplot(bottles, aes(x = temp_est, y = site)) +
  geom_point() +
  labs(x = "Salinity Estimate", y = "Site")

ggplot(bottles, aes(x = temp_est)) +
  geom_histogram(colour="black", fill="white")

# stats -------------------------
pee <- lm(nh4_conc ~ site * depth, bottles_f_may2021)
anova(pee)
summary(pee)

visreg(pee, "depth", by = "site")
