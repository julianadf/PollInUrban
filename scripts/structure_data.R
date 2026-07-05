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
      Cover_type %in% c("Restorable grassland", "Grazed", "Meadow") ~ "SNG",
      Cover_type %in% c("Forest", "Forest on wetland","Temporarily non-forest","Temporarily non-forest on wetland") ~ "Forest",
      Cover_type %in% c("Sea and ocean", "Inland water") ~ "Water",
      Cover_type %in% c("Constructed not building", "Building") ~ "Urban",
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

# Community data
transects <- transects.raw %>% 
  #mutate(City = as.factor(City)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Triplet = as.factor(Triplet)) %>% 
  mutate(Category = as.factor(Category)) %>% 
  mutate(Visit = as.factor(Visit)) %>% 
  rename(Abund = Nr.obs) %>% 
  mutate(Order = as.factor(Order)) %>% 
  mutate(Pollinator.type = as.factor(Pollinator.type))  %>%
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

#write_csv2(community.visit, "community.visit.csV")

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

# 250 m:
rich.visit.250m <- richness.visit.db %>%
  left_join(landscape_250,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(rich.visit.250m, "rich.visit.250m.csv",fileEncoding = "latin1")

# 500 m:
rich.visit.500m <- richness.visit.db %>%
  left_join(landscape_500,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(rich.visit.500m, "rich.visit.500m.csv", , fileEncoding = "latin1")

# 1 km:
rich.visit.1km <- richness.visit.db %>%
  left_join(landscape_1000,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(rich.visit.1km, "rich.visit.1km.csv", , fileEncoding = "latin1")

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

richness.pollinator.db <- community.long.visit %>%
  group_by(City, Site, Triplet, Category, Order) %>%
  summarise(
    richness = n_distinct(Taxon),
    abundance_total = sum(Abund),
    .groups = "drop") %>% 
  filter(Order != "")
richness.pollinator.db

#write.csv2(richness.pollinator.db, "richness.pollinator.db.csv", , fileEncoding = "latin1")
tot.richness.db <- community.long.visit %>%
  group_by(City, Site, Triplet, Category) %>%
  summarise(
    richness = n_distinct(Taxon[!is.na(Taxon) & Taxon != ""]),
    abundance_total = sum(Abund),
    .groups = "drop"
  )
tot.richness.db


#write.csv2(tot.richness.db, "tot.richness.db.csv", , fileEncoding = "latin1")

# Create databases for total richness per site and landscape scales
# 100 m:
tot.richness.100m <- tot.richness.db %>%
  left_join(landscape_100,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(tot.richness.100m, "tot.richness.100m.csv", fileEncoding = "latin1")

# 250 m:
tot.richness.250m <- tot.richness.db %>%
  left_join(landscape_250,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(tot.richness.250m, "tot.richness.250m.csv", fileEncoding = "latin1")


# 500 m:
tot.richness.500m <- tot.richness.db %>%
  left_join(landscape_500,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(tot.richness.500m, "tot.richness.500m.csv", fileEncoding = "latin1")

# 1 km:
tot.richness.1km <- tot.richness.db %>%
  left_join(landscape_1000,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(tot.richness.1km, "tot.richness.1km.csv", fileEncoding = "latin1")

# Informative figures and numbers
# Histograms for landcover types per scale
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
  dplyr::select("Arable", "Forest", "Non-vegetated open land", "Road.Railway", "SNG", "Urban", "Vegetated open land",
                "Water", "Wetland")

cor_mat100 <- cor(corr.100, use = "pairwise.complete.obs")

corrplot(cor_mat100,
  method = "color",
  type = "upper",
  order = "hclust",
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 45)

# Pearson (linear relationships)
pearson.100  <- rcorr(as.matrix(corr.100), type = "pearson")

# Spearman (rank-based, robust to outliers & non-normality)
spearman.100 <- rcorr(as.matrix(corr.100), type = "spearman")

# plot correlations
corrplot(pearson.100$r, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         p.mat = pearson.100$P, sig.level = 0.05, insig = "blank",
         title = "Pearson correlations (ns cells blanked)", mar = c(0,0,1,0))

corrplot(spearman.100$r, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         p.mat = spearman.100$P, sig.level = 0.05, insig = "blank",
         title = "Spearman correlations (ns cells blanked)", mar = c(0,0,1,0))


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
  dplyr::select("Arable", "Forest", "Non-vegetated open land", "Road.Railway", "SNG", "Urban", "Vegetated open land",
                "Water", "Wetland")

cor_mat250 <- cor(corr.250, use = "pairwise.complete.obs")

corrplot(
  cor_mat250,
  method = "color",
  type = "upper",
  order = "hclust",
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 45)

# Pearson (linear relationships)
pearson.250  <- rcorr(as.matrix(corr.250), type = "pearson")

# Spearman (rank-based, robust to outliers & non-normality)
spearman.250 <- rcorr(as.matrix(corr.250), type = "spearman")

# plot correlations
corrplot(pearson.250$r, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         p.mat = pearson.250$P, sig.level = 0.05, insig = "blank",
         title = "Pearson correlations (ns cells blanked)", mar = c(0,0,1,0))

corrplot(spearman.250$r, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         p.mat = spearman.250$P, sig.level = 0.05, insig = "blank",
         title = "Spearman correlations (ns cells blanked)", mar = c(0,0,1,0))


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
  dplyr::select("Arable", "Forest", "Non-vegetated open land", "Road.Railway", "SNG", "Urban", "Vegetated open land",
                "Water", "Wetland")

cor_mat500 <- cor(corr.500, use = "pairwise.complete.obs")

corrplot(
  cor_mat500,
  method = "color",
  type = "upper",
  order = "hclust",
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 45)

# Pearson (linear relationships)
pearson.500  <- rcorr(as.matrix(corr.500), type = "pearson")

# Spearman (rank-based, robust to outliers & non-normality)
spearman.500 <- rcorr(as.matrix(corr.500), type = "spearman")

# plot correlations
corrplot(pearson.500$r, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         p.mat = pearson.500$P, sig.level = 0.05, insig = "blank",
         title = "Pearson correlations (ns cells blanked)", mar = c(0,0,1,0))

corrplot(spearman.500$r, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         p.mat = spearman.500$P, sig.level = 0.05, insig = "blank",
         title = "Spearman correlations (ns cells blanked)", mar = c(0,0,1,0))


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
  dplyr::select("Arable", "Forest", "Non-vegetated open land", "Road.Railway", "SNG", "Urban", "Vegetated open land",
                "Water", "Wetland")

cor_mat1000 <- cor(corr.1000, use = "pairwise.complete.obs")

corrplot(
  cor_mat1000,
  method = "color",
  type = "upper",
  order = "hclust",
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 45)

# Pearson (linear relationships)
pearson.1000  <- rcorr(as.matrix(corr.1000), type = "pearson")

# Spearman (rank-based, robust to outliers & non-normality)
spearman.1000 <- rcorr(as.matrix(corr.1000), type = "spearman")

# plot correlations
corrplot(pearson.1000$r, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         p.mat = pearson.1000$P, sig.level = 0.05, insig = "blank",
         title = "Pearson correlations (ns cells blanked)", mar = c(0,0,1,0))

corrplot(spearman.1000$r, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         p.mat = spearman.1000$P, sig.level = 0.05, insig = "blank",
         title = "Spearman correlations (ns cells blanked)", mar = c(0,0,1,0))


# test correlations again after merging urban (buildings, constructed, roads) and removing other land cover types

sr.100 <- tot.richness.100m %>% 
  mutate(urban.road = Urban + Road.Railway)

sr.250 <- tot.richness.250m %>% 
  mutate(urban.road = Urban + Road.Railway)

sr.500 <- tot.richness.500m %>% 
  mutate(urban.road = Urban + Road.Railway)

sr.1 <- tot.richness.1km %>% 
  mutate(urban.road = Urban + Road.Railway)

# 100 m

corr.100 <- sr.100 %>%
  dplyr::select("Arable", "Forest", "SNG", "urban.road", "Vegetated open land")

# Spearman (rank-based, robust to outliers & non-normality)
spearman.100 <- rcorr(as.matrix(corr.100), type = "spearman")

corrplot(spearman.100$r, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         p.mat = spearman.100$P, sig.level = 0.05, insig = "blank",
         title = "100 m Spearman correlations (ns cells blanked)", mar = c(0,0,1,0))

# 250 m
corr.250 <- sr.250 %>%
  dplyr::select("Arable", "Forest", "SNG", "urban.road", "Vegetated open land")

# Spearman (rank-based, robust to outliers & non-normality)
spearman.250 <- rcorr(as.matrix(corr.250), type = "spearman")

corrplot(spearman.250$r, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         p.mat = spearman.250$P, sig.level = 0.05, insig = "blank",
         title = "250 m Spearman correlations (ns cells blanked)", mar = c(0,0,1,0))


# 500 m
corr.500 <- sr.500 %>%
  dplyr::select("Arable", "Forest", "SNG", "urban.road", "Vegetated open land")

# Spearman (rank-based, robust to outliers & non-normality)
spearman.500 <- rcorr(as.matrix(corr.500), type = "spearman")

corrplot(spearman.500$r, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         p.mat = spearman.500$P, sig.level = 0.05, insig = "blank",
         title = "500 m Spearman correlations (ns cells blanked)", mar = c(0,0,1,0))

# 1 km:

corr.1000 <- sr.1 %>%
  dplyr::select("Arable", "Forest", "SNG", "urban.road", "Vegetated open land")

# Spearman (rank-based, robust to outliers & non-normality)
spearman.1000 <- rcorr(as.matrix(corr.1000), type = "spearman")

corrplot(spearman.1000$r, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         p.mat = spearman.1000$P, sig.level = 0.05, insig = "blank",
         title = "1 km Spearman correlations (ns cells blanked)", mar = c(0,0,1,0))

# Abundance - remove honey bees ----

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

transects.nohoney <- transects %>% 
  filter(Taxon != "Apis mellifera")
transects.nohoney

community.long.visit.nohoney <- transects.nohoney %>%
  group_by(City, Site, Triplet, Category, Visit, Taxon, Order) %>%
  summarise(Abund = sum(Abund, na.rm = TRUE), .groups = "drop")
community.long.visit.nohoney

community.visit.nohoney <- community.long.visit.nohoney %>%
  filter(Taxon != "", !is.na(Taxon)) %>% # 2 empty rows, removed
  group_by(City, Site, Triplet, Category, Visit, Taxon) %>%
  summarise(Abundance = sum(Abund), .groups = "drop") %>%
  pivot_wider(
    names_from = Taxon,
    values_from = Abundance,
    values_fill = 0)
community.visit.nohoney

#write_csv2(community.visit, "community.visit.nohoney.csV")

richness.nohoney.visit.db <- community.long.visit.nohoney %>%
  group_by(City, Site, Triplet, Category, Visit) %>%
  summarise(
    richness = n_distinct(Taxon[!is.na(Taxon) & Taxon != ""]),
    abundance_total = sum(Abund),
    .groups = "drop"
  )
richness.nohoney.visit.db

#databases by Site, visit and scale:

# 100 m:
rich.visit.100m <- richness.nohoney.visit.db %>%
  left_join(landscape_100,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(rich.visit.100m, "rich.visit.100m.nohoney.csv", fileEncoding = "latin1")

# 250 m:
rich.visit.250m <- richness.nohoney.visit.db %>%
  left_join(landscape_250,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(rich.visit.250m, "rich.visit.250m.nohoney.csv",fileEncoding = "latin1")

# 500 m:
rich.visit.500m <- richness.nohoney.visit.db %>%
  left_join(landscape_500,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(rich.visit.500m, "rich.visit.500m.nohoney.csv", , fileEncoding = "latin1")

# 1 km:
rich.visit.1km <- richness.nohoney.visit.db %>%
  left_join(landscape_1000,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(rich.visit.1km, "rich.visit.1km.nohoney.csv", , fileEncoding = "latin1")

# databases by Site and scale:

community.nohoney <- community.long.visit.nohoney %>%
  filter(Taxon != "", !is.na(Taxon)) %>% # 2 empty rows, removed
  group_by(City, Site, Triplet, Category, Taxon) %>%
  summarise(Abundance = sum(Abund), .groups = "drop") %>%
  pivot_wider(
    names_from = Taxon,
    values_from = Abundance,
    values_fill = 0)
community.nohoney

richness.nohoney.pollinator.db <- community.long.visit.nohoney %>%
  group_by(City, Site, Triplet, Category, Order) %>%
  summarise(
    richness = n_distinct(Taxon),
    abundance_total = sum(Abund),
    .groups = "drop") %>% 
  filter(Order != "")
richness.nohoney.pollinator.db

#write.csv2(richness.pollinator.db, "richness.pollinator.db.csv", , fileEncoding = "latin1")
nohoney.tot.richness.db <- community.long.visit.nohoney %>%
  group_by(City, Site, Triplet, Category) %>%
  summarise(
    richness = n_distinct(Taxon[!is.na(Taxon) & Taxon != ""]),
    abundance_total = sum(Abund),
    .groups = "drop"
  )
nohoney.tot.richness.db



# Create databases for total richness per site and landscape scales
# 100 m:
tot.nohoney.richness.100m <- nohoney.tot.richness.db %>%
  left_join(landscape_100,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(tot.nohoney.richness.100m, "tot.nohoney.richness.100m.csv", fileEncoding = "latin1")

# 250 m:
tot.nohoney.richness.250m <- nohoney.tot.richness.db %>%
  left_join(landscape_250,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(tot.nohoney.richness.250m, "tot.nohoney.richness.250m.csv", fileEncoding = "latin1")


# 500 m:
tot.nohoney.richness.500m <- nohoney.tot.richness.db %>%
  left_join(landscape_500,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(tot.nohoney.richness.500m, "tot.nohoney.richness.500m.csv", fileEncoding = "latin1")

# 1 km:
tot.nohoney.richness.1km <- nohoney.tot.richness.db %>%
  left_join(landscape_1000,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(tot.nohoney.richness.1km, "tot.nohoney.richness.1km.csv", fileEncoding = "latin1")


# Abundance - solitary bees only ----

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

transects.wildbees <- transects %>% 
  filter(Taxon != "Apis mellifera" &
  Pollinator.type == "Bee")
  #str_count(Taxon, " ") >= 1)#keep only those with two words (Genus + species)
transects.wildbees

community.long.visit.wildbee <- transects.wildbees %>%
  group_by(City, Site, Triplet, Category, Visit, Taxon, Order) %>%
  summarise(Abund = sum(Abund, na.rm = TRUE), .groups = "drop")
community.long.visit.wildbee

community.visit.wildbee <- community.long.visit.wildbee %>%
  filter(Taxon != "", !is.na(Taxon)) %>% # 2 empty rows, removed
  group_by(City, Site, Triplet, Category, Visit, Taxon) %>%
  summarise(Abundance = sum(Abund), .groups = "drop") %>%
  pivot_wider(
    names_from = Taxon,
    values_from = Abundance,
    values_fill = 0)
community.visit.wildbee

#write.csv2(community.visit.wildbee, "community.visit.wildbee.csV", fileEncoding = "latin1")

richness.wildbee.visit.db <- community.long.visit.wildbee %>%
  group_by(City, Site, Triplet, Category, Visit) %>%
  summarise(
    richness = n_distinct(Taxon[!is.na(Taxon) & Taxon != ""]),
    abundance_total = sum(Abund),
    .groups = "drop"
  )
richness.wildbee.visit.db

#databases by Site, visit and scale:

# 100 m:
rich.visit.100m <- richness.wildbee.visit.db %>%
  left_join(landscape_100,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(rich.visit.100m, "rich.visit.100m.wildbee.csv", fileEncoding = "latin1")

# 250 m:
rich.visit.250m <- richness.wildbee.visit.db %>%
  left_join(landscape_250,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(rich.visit.250m, "rich.visit.250m.wildbee.csv",fileEncoding = "latin1")

# 500 m:
rich.visit.500m <- richness.wildbee.visit.db %>%
  left_join(landscape_500,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(rich.visit.500m, "rich.visit.500m.wildbee.csv", , fileEncoding = "latin1")

# 1 km:
rich.visit.1km <- richness.wildbee.visit.db %>%
  left_join(landscape_1000,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(rich.visit.1km, "rich.visit.1km.wildbee.csv", , fileEncoding = "latin1")

# databases by Site and scale:

community.wildbee <- community.long.visit.wildbee %>%
  filter(Taxon != "", !is.na(Taxon)) %>% # 2 empty rows, removed
  group_by(City, Site, Triplet, Category, Taxon) %>%
  summarise(Abundance = sum(Abund), .groups = "drop") %>%
  pivot_wider(
    names_from = Taxon,
    values_from = Abundance,
    values_fill = 0)
community.wildbee

richness.wildbee.pollinator.db <- community.long.visit.wildbee %>%
  group_by(City, Site, Triplet, Category, Order) %>%
  summarise(
    richness = n_distinct(Taxon),
    abundance_total = sum(Abund),
    .groups = "drop") %>% 
  filter(Order != "")
richness.wildbee.pollinator.db

#write.csv2(richness.wildbee.pollinator.db, "richness.wildbee.pollinator.csv", fileEncoding = "latin1")

wildbee.tot.richness.db <- community.long.visit.wildbee %>%
  group_by(City, Site, Triplet, Category) %>%
  summarise(
    richness = n_distinct(Taxon[!is.na(Taxon) & Taxon != ""]),
    abundance_total = sum(Abund),
    .groups = "drop"
  )
wildbee.tot.richness.db



# Create databases for total richness per site and landscape scales
# 100 m:
tot.wildbee.richness.100m <- wildbee.tot.richness.db %>%
  left_join(landscape_100,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(tot.wildbee.richness.100m, "tot.wildbee.richness.100m.csv", fileEncoding = "latin1")

# 250 m:
tot.wildbee.richness.250m <- wildbee.tot.richness.db %>%
  left_join(landscape_250,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(tot.wildbee.richness.250m, "tot.wildbee.richness.250m.csv", fileEncoding = "latin1")


# 500 m:
tot.wildbee.richness.500m <- wildbee.tot.richness.db %>%
  left_join(landscape_500,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(tot.wildbee.richness.500m, "tot.wildbee.richness.500m.csv", fileEncoding = "latin1")

# 1 km:
tot.wildbee.richness.1km <- wildbee.tot.richness.db %>%
  left_join(landscape_1000,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(tot.wildbee.richness.1km, "tot.wildbee.richness.1km.csv", fileEncoding = "latin1")


# Separate into pollinator groups (without honey  bees)----
transects <- transects.raw %>% 
  #mutate(City = as.factor(City)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Triplet = as.factor(Triplet)) %>% 
  mutate(Category = as.factor(Category)) %>% 
  mutate(Visit = as.factor(Visit)) %>% 
  rename(Abund = Nr.obs) %>% 
  mutate(Order = as.factor(Order), 
         Pollinator.type = ifelse(Pollinator.type == "Burnet moth", "Butterfly", Pollinator.type)) %>% 
  mutate(Pollinator.type = as.factor(Pollinator.type)) %>% 
  mutate(Insect.species = trimws(Insect.species)) %>% 
  rename(Taxon = Insect.species) 
str(transects)

transects.groups <- transects %>% 
  filter(Taxon != "Apis mellifera") 
#str_count(Taxon, " ") >= 1)#keep only those with two words (Genus + species)
transects.groups

community.long.visit.groups <- transects.groups %>%
  group_by(City, Site, Triplet, Category, Visit, Taxon, Pollinator.type) %>%
  summarise(Abund = sum(Abund, na.rm = TRUE), .groups = "drop")
community.long.visit.groups

community.visit.groups <- community.long.visit.groups %>%
  filter(Taxon != "", !is.na(Taxon)) %>% # 2 empty rows, removed
  group_by(City, Site, Triplet, Category, Visit, Pollinator.type) %>%
  summarise(Abundance = sum(Abund), .groups = "drop") %>%
  pivot_wider(
    names_from = Pollinator.type,
    values_from = Abundance,
    values_fill = 0)
community.visit.groups

#write.csv2(community.visit.wildbee, "community.visit.wildbee.csV", fileEncoding = "latin1")

abund.group.visit.db <- community.long.visit.groups %>%
  group_by(City, Site, Triplet, Category, Visit, Pollinator.type) %>%
  summarise(
    abundance_total = sum(Abund),
    .groups = "drop"
  )
abund.group.visit.db

abund.group.tot.db <- community.long.visit.groups %>%
  group_by(City, Site, Triplet, Category, Pollinator.type) %>%
  summarise(
    abundance_total = sum(Abund),
    .groups = "drop"
  )
abund.group.tot.db


#databases by Site, visit and scale:

# Create databases for total abundance per pollinator tyep per site and landscape scales
# 100 m:
tot.group.abund.100m <- abund.group.tot.db %>%
  left_join(landscape_100,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(tot.group.abund.100m, "tot.group.abund.100m.csv", fileEncoding = "latin1")

# 250 m:
tot.group.abund.250m <- abund.group.tot.db %>%
  left_join(landscape_250,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(tot.group.abund.250m, "tot.group.abund.250m.csv", fileEncoding = "latin1")


# 500 m:
tot.group.abund.500m <- abund.group.tot.db %>%
  left_join(landscape_500,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(tot.group.abund.500m, "tot.group.abund.500m.csv", fileEncoding = "latin1")

# 1 km:
tot.group.abund.1km <- abund.group.tot.db %>%
  left_join(landscape_1000,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(tot.group.abund.1km, "tot.group.abund.1km.csv", fileEncoding = "latin1")





abund_summary <- abund.group.visit.db %>%
  filter(
    !(Pollinator.type %in% c("Bumblebee", "Other", "Empty")) &
      !is.na(Pollinator.type)
  ) %>%
  group_by(Category, Visit, Pollinator.type) %>%
  summarise(
    mean_abund = mean(abundance_total),
    sd = sd(abundance_total),
    n = n(),
    se = sd / sqrt(n),
    .groups = "drop"
  )



ggplot(abund_summary,
       aes(x = Visit, y = mean_abund,
           color = Pollinator.type,
           group = Pollinator.type)) +
  
  geom_line() +
  geom_point(size = 2) +
  
  geom_errorbar(aes(
    ymin = mean_abund - se,
    ymax = mean_abund + se
  ), width = 0.1) +
  
  facet_wrap(~ Category) +
  
  theme_bw() +
  labs(
    title = "Seasonal changes in pollinator abundance (excluding bumblebees)",
    x = "Visit",
    y = "Mean abundance",
    color = "Pollinator group"
  )

# Create database with total abundance + honey bee abundance as covariate

abund.db <- transects %>%
  group_by(City, Site, Triplet, Category, Pollinator.type) %>%
  summarise(
    total_abundance = sum(Abund, na.rm = TRUE),
    honeybee_abundance = sum(
      Abund[Taxon == "Apis mellifera"],
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    wild_pollinator_abundance =
      total_abundance - honeybee_abundance
  )
head(abund.db)

# 100 m:
abund.100 <- abund.db %>%
  left_join(landscape_100,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(abund.100, "abund.db.100m.csv", fileEncoding = "latin1")

# 1 km:
abund.1 <- abund.db %>%
  left_join(landscape_1000,
            by = c("City", "Site", "Triplet", "Category"))

#write.csv2(abund.1, "abund.db.1km.csv", fileEncoding = "latin1")
