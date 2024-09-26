# Try using random forest??????
# Em Lim
# Jan 2024


# run one random forest with abundance data in the species cells, and one with the biomass of each species on each survey in each cell

# also remember to split inverts and fish
  # densities are different! the rectangles were different
  # maybe split M1 vs M2???

# should i convert my current abundance/biomass 

# load functions
source("Code/Functions.R")

# Load packages
library(tidyverse)
library(randomForest)
library(TMB)
library(glmmTMB)
library(ggeffects)

# try using rfsrc
library(randomForestSRC)

# Load data from RLS Blitz analysis
# rls_final is the df I've been analysing, with a row for each survey
rls <- read_csv("Output/Output_data/rls_final.csv", show_col_types = FALSE) %>%
  # rls_wider is just the species data. one column per species, one row per survey
  left_join(read_csv("Output/Output_data/rls_wider.csv", show_col_types = FALSE), by = "survey_id") %>%
  as.data.frame() %>%
  mutate(site_code = as.factor(site_code))

# add underscores to the stupid latin names
colnames(rls) <- colnames(rls) %>% 
  str_replace_all(" ","_")


# Load just the rls data so I can play with taxonomic level
rls_sp <- read_csv("Output/Output_data/rls_species.csv", show_col_types = FALSE) %>%
  select(-c(FID, 3:18, 29:33, 36:39)) %>% # cut extra columns
  clean_phylo_names() %>% # function to fix naming errors
  mutate(match = reporting_name == species_name) %>% 
  separate(species_name, c("genus", "species")) %>% # the additional pieced discarded just means it dropped the period for all the spp. species
  mutate(species = ifelse(species == "spp" | species == "shrimp", NA, species),
         genus = ifelse(genus == "Cottidae" |
                          genus == "Paguroidea" |
                          genus == "Unidentified" |
                          genus ==  "Brachyura" |
                          genus == "Actinopterygii" |
                          genus == "Embiotocidae" |
                          genus == "Percidae", NA, genus)
         )

# separate into m1 and m2
m1 <- rls_sp %>%
  filter(method == 1)

m2 <- rls_sp %>%
  filter(method == 2)

# use this to check the species with mismatching reporting vs scientific names
# false <- rls_sp %>%
#   filter(match == FALSE) %>%
#   select(-c(1:3, 11:13)) %>%
#   unique()

# try grouping by species and see how many we get
# omit indivs not identified to species? na.omit(species) ?
species <- rls_sp %>%
  group_by(genus, species) %>%
  summarise()
# 124 groups total
# 29 in m1, 107 m2

# try grouping by genus 
genus <- m2 %>%
  group_by(genus) %>%
  summarise(n_species = n_distinct(species))%>%
  drop_na(genus)
# 94 groups total
# 21 m1, 85 m2

# try grouping by family 
family <- m1 %>%
  group_by(family) %>%
  summarise(n_genus = n_distinct(genus),
            n_species = n_distinct(species)) %>%
  drop_na(family)
# 62 groups total
# 13 m1, 55 m2

# try grouping by order 
order <- m1 %>%
  group_by(order) %>%
  summarise(n_family = n_distinct(family),
            n_genus = n_distinct(genus),
            n_species = n_distinct(species))%>%
  drop_na(order)
# 26 groups total
# 5 groups m1, 23 m2 

# try grouping by class 
class <- rls_sp %>%
  group_by(class) %>%
  summarise(n_order = n_distinct(order),
            n_family = n_distinct(family),
            n_genus = n_distinct(genus),
            n_species = n_distinct(species))
# 12 groups

# So I'll probably go with family
# Am I supposed to make a model for each family???? -----

# pal for plotting
pal <- viridis::viridis(10)
pal3 <- c(pal[10], pal[8], pal[5])

# Top 4 M1 families -----
m1_fam_list <- m1 %>%
  group_by(family) %>%
  summarise(sum = sum(total)) %>%
  drop_na(family) %>%
  arrange(desc(sum)) %>%
  head(4) # Remember Cottidae is mostly an M2 fish

# create plots for m1
plot1 <- plot_sp_fun(m1, family = "Sebastidae", diagnose = TRUE)
plot2 <- plot_sp_fun(m1, family = "Hexagrammidae", diagnose = TRUE)
plot3 <- plot_sp_fun(m1, family = "Gobiidae", diagnose = TRUE)
plot4 <- plot_sp_fun(m1, family = "Embiotocidae", diagnose = TRUE)

m1_plot <- plot1 + plot2 + plot3 + plot4 & xlab(NULL) & ylab(NULL)

wrap_elements(panel = m1_plot) +
  labs(tag = expression(paste("Animal abundance/m"^2))) +
  theme(
    plot.tag = element_text(size = rel(2)),
    plot.tag.position = "bottom"
  )

# ggsave("Output/Figures/Kieran_m1_biomass.png", device = "png", height = 9, width = 16, dpi = 400)

# All fish -----
m1_fam_list <- rls_sp %>%
  group_by(family) %>%
  summarise(sum = sum(total)) %>%
  drop_na(family) %>%
  arrange(desc(sum))

# if i just do all fish
# Gobiidae, Cottidae, Hexagrammidae, Sebastidae

plot5 <- plot_sp_fun(rls_sp, family = "Gobiidae")
plot6 <- plot_sp_fun(rls_sp, family = "Cottidae")
plot7 <- plot_sp_fun(rls_sp, family = "Hexagrammidae")
plot8 <- plot_sp_fun(rls_sp, family = "Sebastidae")

fish_plot1 <- plot5 + plot6 + plot7 + plot8 & xlab(NULL) & ylab(NULL)

fish_plot <- wrap_elements(panel = fish_plot1) +
  labs(tag = expression(paste("Animal abundance/m"^2))) +
  theme(
    plot.tag = element_text(size = rel(2)),
    plot.tag.position = "bottom"
  )
print(fish_plot)

#ggsave("Output/Figures/Kieran_fish_biomass.png", device = "png", height = 9, width = 16, dpi = 400)


# can I do it for M2 -----
# Top 4 M1 families
m2_fam_list <- m2 %>%
  group_by(family) %>%
  summarise(sum = sum(total)) %>%
  drop_na(family) %>%
  arrange(desc(sum)) %>%
  head(4) # Remember Cottidae is mostly an M2 fish

# create plots
plot9 <- plot_sp_fun(m2, family = "Strongylocentrotidae", diagnose = TRUE)
plot10 <- plot_sp_fun(m2, family = "Turbinidae", diagnose = TRUE)
plot11 <- plot_sp_fun(m2, family = "Asterinidae", diagnose = TRUE)
plot12 <- plot_sp_fun(m2, family = "Gobiidae", diagnose = TRUE)
plot13 <- plot_sp_fun(m2, family = "Stichopodidae", diagnose = TRUE)
plot14 <- plot_sp_fun(m2, family = "Asteriidae", diagnose = TRUE)

m1_plot <- plot9 + plot10 + plot11 + plot12 + plot13 + plot14 & xlab(NULL) & ylab(NULL)


wrap_elements(panel = m1_plot) +
  labs(tag = expression(paste("Animal abundance/m"^2))) +
  theme(
    plot.tag = element_text(size = rel(2)),
    plot.tag.position = "bottom"
  )

#ggsave("Output/Figures/Kieran_m2_biomass.png", device = "png", height = 9, width = 16, dpi = 400)


# do it once manually
fam1 <- m1 %>% filter(family == "Gobiidae") %>%
  group_by(survey_id) %>%
  summarise(total = sum(total),
            weight_size_class_sum = sum(weight_size_class_sum)) %>%
  ungroup() %>%
  mutate(total_weight_g = weight_size_class_sum*1000)

fam <- rls %>%
  left_join(fam1, by = "survey_id") %>%
  replace_na(list(total = 0, total_weight_g = 0)) %>%
  mutate(total_scale = c(scale(total)))

# model?
mod_fam <- glmmTMB(nh4_avg ~ total_weight_g + tide_scale + 
                     (1|year) + (1|site_code), 
                   family = Gamma(link = 'log'),
                   data = fam)
print(summary(mod_fam))
print(plot(DHARMa::simulateResiduals(mod_fam)))

# plot predictions
predict_fam <- ggpredict(mod_fam, terms = c("total_weight_g", "tide_scale [-0.279]")) %>% 
  mutate(total_weight_g = x,
         tide_cat = factor(as.factor(case_when(group == "-1.067" ~ "Ebb",
                                               group == "-0.279" ~ "Slack",
                                               group == "1.066" ~ "Flood")),
                           levels = c("Ebb", "Slack", "Flood")))

# now plot these predictions
plot_fam <- plot_pred_fam(raw_data = fam,
            predict_data = predict_fam, 
            x_var = total_weight_g, y_var = nh4_avg, 
            pal = pal3)

print(plot_fam)

# Try to make a random forest I guess????

# let's make a simple one first
(RF.model = randomForest(nh4_avg ~ abundance_stand + tide_stand + shannon_stand + depth_avg_stand, data = rls))

plot(RF.model)
randomForest::importance(RF.model)
varImpPlot(RF.model)
partialPlot(RF.model, rls, "abundance_stand")

# try to get all species names so i can include them without typing each out?
species_cols <- paste(colnames(rls[30:160]), collapse = " + ")

# try to make the formula
my_formula <- as.formula(paste("nh4_avg ~ abundance_stand + tide_stand + shannon_stand + depth_avg_stand +", paste(colnames(rls[30:160]), collapse = " + ")))

my_formula <- as.formula(paste("nh4_avg ~ ", paste(colnames(rls[30:160]), collapse = " + ")))

# try to use that formula
(RF.model = randomForest(my_formula, data = rls))

# evaluate that?
plot(RF.model)
randomForest::importance(RF.model)
varImpPlot(RF.model)
partialPlot(RF.model, rls, "abundance_stand")



# I can also subset df so i can just call all columns
rls_reduced <- rls %>%
  select(c(nh4_avg, abundance_stand, tide_stand, shannon_stand, depth_avg_stand, 30:160))

(RF.model = randomForest(nh4_avg ~ ., data = rls_reduced))
