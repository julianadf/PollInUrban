rm(list=ls())
library(here)
library(tidyverse)
library(vegan)
library(lattice)
library(ggpubr)
library(effects)
library(visreg)
library(DHARMa)
library(MASS)
library(lme4)
library(car)
library(png)
library(cowplot)
library(visreg)
library(MuMIn)
library(glmmTMB)
library(sjPlot)
library(emmeans)
library(mgcv)
library(mgcViz)
library(corrplot)
library(Hmisc)

# Load data
landscape.raw <- read.csv2(here("data", "Landscape.csv" ), fileEncoding = "latin1")
transects.raw <- read.csv2(here("data", "Transect.csv" ), fileEncoding = "latin1")

# Create database for transect data and surrounding landscape
# landscape data
landscape.merge <- landscape.raw %>%
  mutate(City = as.factor(City)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Triplet = as.factor(Triplet)) %>%
  mutate(Category = as.factor(Category)) %>% 
  mutate(Cover_type = as.factor(Cover_type)) %>% 
  mutate(Cover_group = case_when(
    Cover_type %in% c("Restorable grassland", "Grazed", "Meadow", "Vegetated open land") ~ "Vegetated.open.land",
    Cover_type %in% c("Forest", "Forest on wetland","Temporarily non-forest","Temporarily non-forest on wetland") ~ "Forest",
    Cover_type %in% c("Sea and ocean", "Inland water") ~ "Water",
    Cover_type %in% c("Constructed not building", "Building", "Road.Railway") ~ "Urban",
    TRUE ~ as.character(Cover_type))) 

landscape <- landscape.merge %>% 
  group_by(City, Site, X, Y, Triplet, Category, Hab_Area, Cover_group) %>%
  summarise(
    Proportion_100m = sum(Proportion_100m, na.rm = TRUE)/100,
    Proportion_250m = sum(Proportion_250m, na.rm = TRUE)/100,
    Proportion_500m = sum(Proportion_500m, na.rm = TRUE)/100,
    Proportion_1km = sum(Proportion_1km, na.rm = TRUE)/100,
    .groups = "drop")
landscape

# per scale:
landscape_100 <- landscape %>%
  filter(Category != "", !is.na(Category)) %>% 
  group_by(City, Site, X, Y, Triplet, Category, Hab_Area, Cover_group) %>%
  summarise(Proportion_100m = sum(Proportion_100m), .groups = "drop") %>%
  pivot_wider(
    names_from = Cover_group,
    values_from = Proportion_100m,
    values_fill = 0)
landscape_100

landscape_250 <- landscape %>%
  filter(Category != "", !is.na(Category)) %>% 
  group_by(City, Site, X, Y, Triplet, Category, Hab_Area, Cover_group) %>%
  summarise(Proportion_250m = sum(Proportion_250m), .groups = "drop") %>%
  pivot_wider(
    names_from = Cover_group,
    values_from = Proportion_250m,
    values_fill = 0)
landscape_250

landscape_500 <- landscape %>%
  filter(Category != "", !is.na(Category)) %>% 
  group_by(City, Site, X, Y, Triplet, Category, Hab_Area, Cover_group) %>%
  summarise(Proportion_500m = sum(Proportion_500m), .groups = "drop") %>%
  pivot_wider(
    names_from = Cover_group,
    values_from = Proportion_500m,
    values_fill = 0)
landscape_500

landscape_1000 <- landscape %>%
  filter(Category != "", !is.na(Category)) %>% 
  group_by(City, Site, X, Y, Triplet, Category, Hab_Area, Cover_group) %>%
  summarise(Proportion_1km = sum(Proportion_1km), .groups = "drop") %>%
  pivot_wider(
    names_from = Cover_group,
    values_from = Proportion_1km,
    values_fill = 0)
landscape_1000

# Correlations ----
# 100 m:
ggplot(landscape, aes(x = Proportion_100m)) +
geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
  facet_wrap(~ Cover_group, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Distribution of land cover types (100 m buffer)",
    x = "Proportion",
    y = "Frequency")


corr.100 <- landscape_100 %>%
  dplyr::select("Arable", "Forest", "Non-vegetated open land", "Vegetated.open.land", "Urban", "Water", "Wetland")

cor_mat100 <- cor(corr.100, use = "pairwise.complete.obs")

# Spearman (rank-based, robust to outliers & non-normality)
spearman.100 <- rcorr(as.matrix(corr.100), type = "spearman")

# plot correlations
corrplot(spearman.100$r, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         p.mat = spearman.100$P, sig.level = 0.05, insig = "blank",
         title = "Spearman correlations 100 m (ns cells blanked)", mar = c(0,0,1,0))


# 250 m: 
ggplot(landscape, aes(x = Proportion_250m)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
  facet_wrap(~ Cover_group, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Distribution of land cover types (250 m buffer)",
    x = "Proportion",
    y = "Frequency")


corr.250 <- landscape_250 %>%
  dplyr::select("Arable", "Forest", "Non-vegetated open land", "Vegetated.open.land", "Urban", "Water", "Wetland")

cor_mat250 <- cor(corr.250, use = "pairwise.complete.obs")

# Spearman (rank-based, robust to outliers & non-normality)
spearman.250 <- rcorr(as.matrix(corr.250), type = "spearman")

# plot correlations
corrplot(spearman.250$r, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         p.mat = spearman.250$P, sig.level = 0.05, insig = "blank",
         title = "Spearman correlations 250 m (ns cells blanked)", mar = c(0,0,1,0))

# 500 m: 
ggplot(landscape, aes(x = Proportion_500m)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
  facet_wrap(~ Cover_group, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Distribution of land cover types (500 m buffer)",
    x = "Proportion",
    y = "Frequency")


corr.500 <- landscape_500 %>%
  dplyr::select("Arable", "Forest", "Non-vegetated open land", "Vegetated.open.land", "Urban", "Water", "Wetland")

cor_mat500 <- cor(corr.500, use = "pairwise.complete.obs")

# Spearman (rank-based, robust to outliers & non-normality)
spearman.500 <- rcorr(as.matrix(corr.500), type = "spearman")

# plot correlations
corrplot(spearman.500$r, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         p.mat = spearman.500$P, sig.level = 0.05, insig = "blank",
         title = "Spearman correlations 500 m (ns cells blanked)", mar = c(0,0,1,0))

# 1 km: 
ggplot(landscape, aes(x = Proportion_1km)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
  facet_wrap(~ Cover_group, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Distribution of land cover types (1 km buffer)",
    x = "Proportion",
    y = "Frequency")


corr.1000 <- landscape_1000 %>%
  dplyr::select("Arable", "Forest", "Non-vegetated open land", "Vegetated.open.land", "Urban", "Water", "Wetland")

cor_mat1000 <- cor(corr.1000, use = "pairwise.complete.obs")

# Spearman (rank-based, robust to outliers & non-normality)
spearman.1000 <- rcorr(as.matrix(corr.1000), type = "spearman")

# plot correlations
corrplot(spearman.1000$r, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         p.mat = spearman.1000$P, sig.level = 0.05, insig = "blank",
         title = "Spearman correlations 1 km (ns cells blanked)", mar = c(0,0,1,0))

# Test correlation between spatial scales
cor_by_cover <- landscape %>%
  group_split(Cover_group) %>%
  map(~ {
    vars <- .x %>%
      dplyr::select(Proportion_100m,
             Proportion_250m,
             Proportion_500m,
             Proportion_1km)
    
    rcorr(as.matrix(vars), type = "spearman")
  })

names(cor_by_cover) <- unique(landscape$Cover_group)

# forest
corrplot(
  cor_by_cover[["Forest"]]$r,
  method = "color",
  type = "upper",
  tl.col = "black",
  tl.srt = 45,
  addCoef.col = "black",
  p.mat = cor_by_cover[["Forest"]]$P,
  sig.level = 0.05,
  insig = "blank",
  title = "Spearman correlations - Forest",
  mar = c(0,0,1,0)
)

# arable
corrplot(
  cor_by_cover[["Arable"]]$r,
  method = "color",
  type = "upper",
  tl.col = "black",
  tl.srt = 45,
  addCoef.col = "black",
  p.mat = cor_by_cover[["Arable"]]$P,
  sig.level = 0.05,
  insig = "blank",
  title = "Spearman correlations - Arable",
  mar = c(0,0,1,0)
)

# Non-vegetated open land
corrplot(
  cor_by_cover[["Non-vegetated open land"]]$r,
  method = "color",
  type = "upper",
  tl.col = "black",
  tl.srt = 45,
  addCoef.col = "black",
  p.mat = cor_by_cover[["Non-vegetated open land"]]$P,
  sig.level = 0.05,
  insig = "blank",
  title = "Spearman correlations - Non-vegetated open land",
  mar = c(0,0,1,0)
)

# Urban
corrplot(
  cor_by_cover[["Urban"]]$r,
  method = "color",
  type = "upper",
  tl.col = "black",
  tl.srt = 45,
  addCoef.col = "black",
  p.mat = cor_by_cover[["Urban"]]$P,
  sig.level = 0.05,
  insig = "blank",
  title = "Spearman correlations - Urban",
  mar = c(0,0,1,0)
)

# Vegetated.open.land
corrplot(
  cor_by_cover[["Vegetated.open.land"]]$r,
  method = "color",
  type = "upper",
  tl.col = "black",
  tl.srt = 45,
  addCoef.col = "black",
  p.mat = cor_by_cover[["Vegetated.open.land"]]$P,
  sig.level = 0.05,
  insig = "blank",
  title = "Spearman correlations - Vegetated.open.land",
  mar = c(0,0,1,0)
)

# Water
corrplot(
  cor_by_cover[["Water"]]$r,
  method = "color",
  type = "upper",
  tl.col = "black",
  tl.srt = 45,
  addCoef.col = "black",
  p.mat = cor_by_cover[["Water"]]$P,
  sig.level = 0.05,
  insig = "blank",
  title = "Spearman correlations - Water",
  mar = c(0,0,1,0)
)

# Wetland
corrplot(
  cor_by_cover[["Wetland"]]$r,
  method = "color",
  type = "upper",
  tl.col = "black",
  tl.srt = 45,
  addCoef.col = "black",
  p.mat = cor_by_cover[["Wetland"]]$P,
  sig.level = 0.05,
  insig = "blank",
  title = "Spearman correlations - Wetland",
  mar = c(0,0,1,0)
)


# Create species richness database ----
# Community data
transects <- transects.raw %>% 
  #mutate(City = as.factor(City)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Triplet = as.factor(Triplet)) %>% 
  mutate(Category = as.factor(Category)) %>% 
  mutate(Visit = as.factor(Visit)) %>% 
  rename(Abund = Nr.obs) %>% 
  mutate(Order = as.factor(Order)) %>% 
  mutate(Pollinator.type = as.factor(Pollinator.type)) %>% 
  mutate(Insect.species = trimws(Insect.species)) %>% 
  rename(Taxon = Insect.species)
str(transects)

community.long.visit <- transects %>%
  group_by(City, Site, Triplet, Category, Visit, Taxon, Order) %>%
  summarise(Abund = sum(Abund, na.rm = TRUE), .groups = "drop")
community.long.visit

community.visit <- community.long.visit %>%
  filter(Taxon != "", !is.na(Taxon)) %>% # 2 empty rows, removed
  group_by(City, Site, Triplet, Category, Visit, Taxon) %>%
  summarise(Abundance = sum(Abund), .groups = "drop") %>%
  pivot_wider(
    names_from = Taxon,
    values_from = Abundance,
    values_fill = 0)
community.visit

#write.csv2(community.visit, "community.visit.csV", fileEncoding = "latin1")

richness.visit.db <- community.long.visit %>%
  group_by(City, Site, Triplet, Category, Visit) %>%
  summarise(
    richness = n_distinct(Taxon[!is.na(Taxon) & Taxon != ""]),
    abundance_total = sum(Abund),
    .groups = "drop"
  )
richness.visit.db

#databases by Site, visit and scale:

# 100 m:
rich.visit.100m <- richness.visit.db %>%
  left_join(landscape_100,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(rich.visit.100m, "rich.visit.100m.csv", fileEncoding = "latin1")

# 1 km:
rich.visit.1km <- richness.visit.db %>%
  left_join(landscape_1000,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(rich.visit.1km, "rich.visit.1km.csv", fileEncoding = "latin1")

# databases by Site and scale:

community <- community.long.visit %>%
  filter(Taxon != "", !is.na(Taxon)) %>% # 2 empty rows, removed
  group_by(City, Site, Triplet, Category, Taxon) %>%
  summarise(Abundance = sum(Abund), .groups = "drop") %>%
  pivot_wider(
    names_from = Taxon,
    values_from = Abundance,
    values_fill = 0)
community

#write.csv2(community, "community.csv", , fileEncoding = "latin1")

tot.richness.db <- community.long.visit %>%
  group_by(City, Site, Triplet, Category) %>%
  summarise(
    richness = n_distinct(Taxon[!is.na(Taxon) & Taxon != ""]),
    abundance_total = sum(Abund),
    .groups = "drop"
  )
tot.richness.db


write.csv2(tot.richness.db, "tot.richness.db.csv", fileEncoding = "latin1")

# Create databases for total richness per site and landscape scales
# 100 m:
tot.richness.100m <- tot.richness.db %>%
  left_join(landscape_100,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(tot.richness.100m, "tot.richness.100m.csv", fileEncoding = "latin1")


# 1 km:
tot.richness.1km <- tot.richness.db %>%
  left_join(landscape_1000,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(tot.richness.1km, "tot.richness.1km.csv", fileEncoding = "latin1")

# abundance - without honeybees (honeybee abundance as a separate variable) ----

