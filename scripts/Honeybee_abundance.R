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
library(performance)

# Abundance: all pollinators
# Load data
abund.100.raw <- read.csv2(here("data/databases", "abund.db.100m.csv" ), fileEncoding = "latin1")
abund.1.raw <- read.csv2(here("data/databases", "abund.db.1km.csv" ), fileEncoding = "latin1")

# Bumblebees
# structure data
# 100 m:-----
# Honey bees per site
honeybee_df <- abund.100.raw %>%
  group_by(Site) %>%
  summarise(honeybee = sum(honeybee_abundance, na.rm = TRUE), .groups = "drop")

# Keep all sites
all_sites <- abund.100.raw %>% 
  distinct(City, Site, Triplet, Category, X, Y, Hab_Area, Arable, Forest, Non.vegetated.open.land,
           Urban, Vegetated.open.land, Water,Wetland)

# Honeybee abundance per site
hb_abund <- abund.100.raw %>%
  filter(Pollinator.type == "Honey bee") %>% 
  group_by(City, Site, Triplet, Category, X, Y, Hab_Area, Arable, Forest, Non.vegetated.open.land,
           Urban, Vegetated.open.land, Water, Wetland) %>%
  summarise(abundance.hb = sum(total_abundance), .groups = "drop")

abund.100.hbs <- all_sites %>%
  left_join(hb_abund,
            by = c(
              "City", "Site", "Triplet", "Category",
              "X", "Y", "Hab_Area",
              "Arable", "Forest",
              "Non.vegetated.open.land",
              "Urban",
              "Vegetated.open.land",
              "Water",
              "Wetland")) %>%
  mutate(abundance.hb = replace_na(abundance.hb, 0))

# Add honey bee abundance
hbs.db.100 <- abund.100.hbs %>%
  left_join(honeybee_df, by = "Site")

hbs.db.100 <- hbs.db.100 |> 
  mutate(City = as.factor(City)) |> 
  mutate(Site = as.factor(Site)) |> 
  mutate(Triplet = as.factor(Triplet)) |> 
  mutate(Category = as.factor(Category)) 

# Landscape analysis - 100 m: ----
# 1. base model
base.hbs.100 <- glmmTMB(abundance.hb ~ Category + scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + 
                          scale(Forest) + 
                          (1 | City/Triplet),
                        family = nbinom2,
                        data = hbs.db.100)

summary(base.hbs.100)
mod_dharma1 <- base.hbs.100 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(base.hbs.100) 
AICc(base.hbs.100)
plot(allEffects(base.hbs.100))

res.100 <- simulateResiduals(base.hbs.100)
testSpatialAutocorrelation(
  res.100,
  x = tapply(hbs.db.100$X, hbs.db.100$Site, mean),
  y = tapply(hbs.db.100$Y, hbs.db.100$Site, mean)) 


# 2. Test interaction category * habitat area
cat.area <- glmmTMB(abundance.hb ~ Category * scale(Hab_Area) + 
                      scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                      (1 | City/Triplet),
                    family = nbinom2,
                    data = hbs.db.100)

summary(cat.area)
mod_dharma1 <- cat.area %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.area) 
AICc(cat.area)
plot(allEffects(cat.area))

res.100 <- simulateResiduals(cat.area)
testSpatialAutocorrelation(
  res.100,
  x = tapply(hbs.db.100$X, hbs.db.100$Site, mean),
  y = tapply(hbs.db.100$Y, hbs.db.100$Site, mean)) 

# 3. Test category  * urban
cat.urban <- glmmTMB(abundance.hb ~ Category * scale(Urban) +
                       scale(Hab_Area) + scale(Vegetated.open.land) + scale(Forest) + 
                       (1 | City/Triplet),
                     family = nbinom2,
                     data = hbs.db.100)

summary(cat.urban)
mod_dharma1 <- cat.urban %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.urban) 
AICc(cat.urban)
plot(allEffects(cat.urban))

res.100 <- simulateResiduals(cat.urban)
testSpatialAutocorrelation(
  res.100,
  x = tapply(hbs.db.100$X, hbs.db.100$Site, mean),
  y = tapply(hbs.db.100$Y, hbs.db.100$Site, mean)) 

# 4. Test interaction category * vegetated open land
cat.vol <- glmmTMB(abundance.hb ~ Category * scale(Vegetated.open.land) +
                     scale(Urban) + scale(Hab_Area) + scale(Forest) + 
                     (1 | City/Triplet),
                   family = nbinom2,
                   data = hbs.db.100)

summary(cat.vol)
mod_dharma1 <- cat.vol %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.vol) 
AICc(cat.vol)
plot(allEffects(cat.vol))

res.100 <- simulateResiduals(cat.vol)
testSpatialAutocorrelation(
  res.100,
  x = tapply(hbs.db.100$X, hbs.db.100$Site, mean),
  y = tapply(hbs.db.100$Y, hbs.db.100$Site, mean)) 

# 5. Test interaction category * forest
cat.forest <- glmmTMB(abundance.hb ~ Category * scale(Forest) +
                        scale(Vegetated.open.land) + scale(Urban) + scale(Hab_Area) +
                        (1 | City/Triplet),
                      family = nbinom2,
                      data = hbs.db.100)

summary(cat.forest)
mod_dharma1 <- cat.forest %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.forest) 
AICc(cat.forest)
plot(allEffects(cat.forest))

res.100 <- simulateResiduals(cat.forest)
testSpatialAutocorrelation(
  res.100,
  x = tapply(hbs.db.100$X, hbs.db.100$Site, mean),
  y = tapply(hbs.db.100$Y, hbs.db.100$Site, mean)) 

# 6. Test interaction habitat area * urban
area.urban <- glmmTMB(abundance.hb ~ Category + scale(Hab_Area) * scale(Urban) +
                        scale(Forest) + scale(Vegetated.open.land) +  
                        (1 | City/Triplet),
                      family = nbinom2,
                      data = hbs.db.100)

summary(area.urban)
mod_dharma1 <- area.urban %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.urban) 
AICc(area.urban)
plot(allEffects(area.urban))

res.100 <- simulateResiduals(area.urban)
testSpatialAutocorrelation(
  res.100,
  x = tapply(hbs.db.100$X, hbs.db.100$Site, mean),
  y = tapply(hbs.db.100$Y, hbs.db.100$Site, mean)) 


# 7.Test interaction habitat area * vegetated open land
area.vol <- glmmTMB(abundance.hb ~ Category + scale(Hab_Area) * scale(Vegetated.open.land) +
                      scale(Forest) + scale(Urban) +
                      (1 | City/Triplet),
                    family = nbinom2,
                    data = hbs.db.100)

summary(area.vol)
mod_dharma1 <- area.vol %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.vol) 
AICc(area.vol)
plot(allEffects(area.vol))

res.100 <- simulateResiduals(area.vol)
testSpatialAutocorrelation(
  res.100,
  x = tapply(hbs.db.100$X, hbs.db.100$Site, mean),
  y = tapply(hbs.db.100$Y, hbs.db.100$Site, mean)) 

# 8. Test interaction habitat area * forest
area.forest <- glmmTMB(abundance.hb ~ Category + scale(Hab_Area) * scale(Forest) +
                         scale(Vegetated.open.land) + scale(Urban) +
                         (1 | City/Triplet),
                       family = nbinom2,
                       data = hbs.db.100)

summary(area.forest)
mod_dharma1 <- area.forest %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.forest) 
AICc(area.forest)
plot(allEffects(area.forest))

res.100 <- simulateResiduals(area.forest)
testSpatialAutocorrelation(
  res.100,
  x = tapply(hbs.db.100$X, hbs.db.100$Site, mean),
  y = tapply(hbs.db.100$Y, hbs.db.100$Site, mean)) 

# Collect AICc values

mods <- list(
  base = base.hbs.100,
  cat.area = cat.area, 
  cat.urban = cat.urban,
  cat.open = cat.vol,
  cat.forest = cat.forest,
  area.urban = area.urban, 
  area.vol = area.vol, 
  area.forest = area.forest
)

map_dfr(mods, AIC, .id = "model")
map_dfr(mods, AICc, .id = "model")


aic_table <- map_dfr(
  mods,
  ~ data.frame(AICc = AICc(.x)),
  .id = "model"
)

aic_table <- aic_table %>%
  arrange(AICc) %>%
  mutate(delta_AICc = AICc - min(AICc))

aic_table

# 1 km:-----
# Honey bees per site
honeybee_df <- abund.1.raw %>%
  group_by(Site) %>%
  summarise(honeybee = sum(honeybee_abundance, na.rm = TRUE), .groups = "drop")

# Keep all sites
all_sites <- abund.1.raw %>% 
  distinct(City, Site, Triplet, Category, X, Y, Hab_Area, Arable, Forest, Non.vegetated.open.land,
           Urban, Vegetated.open.land, Water,Wetland)

# Honeybee abundance per site
hb_abund <- abund.1.raw %>%
  filter(Pollinator.type == "Honey bee") %>% 
  group_by(City, Site, Triplet, Category, X, Y, Hab_Area, Arable, Forest, Non.vegetated.open.land,
           Urban, Vegetated.open.land, Water, Wetland) %>%
  summarise(abundance.hb = sum(total_abundance), .groups = "drop")

abund.1.hbs <- all_sites %>%
  left_join(hb_abund,
            by = c(
              "City", "Site", "Triplet", "Category",
              "X", "Y", "Hab_Area",
              "Arable", "Forest",
              "Non.vegetated.open.land",
              "Urban",
              "Vegetated.open.land",
              "Water",
              "Wetland")) %>%
  mutate(abundance.hb = replace_na(abundance.hb, 0))

# Add honey bee abundance
hbs.db.1 <- abund.1.hbs %>%
  left_join(honeybee_df, by = "Site")

hbs.db.1 <- hbs.db.1 |> 
  mutate(City = as.factor(City)) |> 
  mutate(Site = as.factor(Site)) |> 
  mutate(Triplet = as.factor(Triplet)) |> 
  mutate(Category = as.factor(Category)) 

# Landscape analysis - 100 m: ----
# 1. base model
base.hbs.1 <- glmmTMB(abundance.hb ~ Category + scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + 
                          scale(Forest) + 
                          (1 | City/Triplet),
                        family = nbinom2,
                        data = hbs.db.1)

summary(base.hbs.1)
mod_dharma1 <- base.hbs.1 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(base.hbs.1) 
AICc(base.hbs.1)
plot(allEffects(base.hbs.1))

res.1 <- simulateResiduals(base.hbs.1)
testSpatialAutocorrelation(
  res.1,
  x = tapply(hbs.db.1$X, hbs.db.1$Site, mean),
  y = tapply(hbs.db.1$Y, hbs.db.1$Site, mean)) 


# 2. Test interaction category * habitat area
cat.area <- glmmTMB(abundance.hb ~ Category * scale(Hab_Area) + 
                      scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                      (1 | City/Triplet),
                    family = nbinom2,
                    data = hbs.db.1)

summary(cat.area)
mod_dharma1 <- cat.area %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.area) 
AICc(cat.area)
plot(allEffects(cat.area))

res.1 <- simulateResiduals(cat.area)
testSpatialAutocorrelation(
  res.1,
  x = tapply(hbs.db.1$X, hbs.db.1$Site, mean),
  y = tapply(hbs.db.1$Y, hbs.db.1$Site, mean)) 

# 3. Test category  * urban
cat.urban <- glmmTMB(abundance.hb ~ Category * scale(Urban) +
                       scale(Hab_Area) + scale(Vegetated.open.land) + scale(Forest) + 
                       (1 | City/Triplet),
                     family = nbinom2,
                     data = hbs.db.1)

summary(cat.urban)
mod_dharma1 <- cat.urban %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.urban) 
AICc(cat.urban)
plot(allEffects(cat.urban))

res.1 <- simulateResiduals(cat.urban)
testSpatialAutocorrelation(
  res.1,
  x = tapply(hbs.db.1$X, hbs.db.1$Site, mean),
  y = tapply(hbs.db.1$Y, hbs.db.1$Site, mean)) 

# 4. Test interaction category * vegetated open land
cat.vol <- glmmTMB(abundance.hb ~ Category * scale(Vegetated.open.land) +
                     scale(Urban) + scale(Hab_Area) + scale(Forest) + 
                     (1 | City/Triplet),
                   family = nbinom2,
                   data = hbs.db.1)

summary(cat.vol)
mod_dharma1 <- cat.vol %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.vol) 
AICc(cat.vol)
plot(allEffects(cat.vol))

res.1 <- simulateResiduals(cat.vol)
testSpatialAutocorrelation(
  res.1,
  x = tapply(hbs.db.1$X, hbs.db.1$Site, mean),
  y = tapply(hbs.db.1$Y, hbs.db.1$Site, mean)) 

# 5. Test interaction category * forest
cat.forest <- glmmTMB(abundance.hb ~ Category * scale(Forest) +
                        scale(Vegetated.open.land) + scale(Urban) + scale(Hab_Area) +
                        (1 | City/Triplet),
                      family = nbinom2,
                      data = hbs.db.1)

summary(cat.forest)
mod_dharma1 <- cat.forest %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.forest) 
AICc(cat.forest)
plot(allEffects(cat.forest))

res.1 <- simulateResiduals(cat.forest)
testSpatialAutocorrelation(
  res.1,
  x = tapply(hbs.db.1$X, hbs.db.1$Site, mean),
  y = tapply(hbs.db.1$Y, hbs.db.1$Site, mean)) 

# 6. Test interaction habitat area * urban
area.urban <- glmmTMB(abundance.hb ~ Category + scale(Hab_Area) * scale(Urban) +
                        scale(Forest) + scale(Vegetated.open.land) +  
                        (1 | City/Triplet),
                      family = nbinom2,
                      data = hbs.db.1)

summary(area.urban)
mod_dharma1 <- area.urban %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.urban) 
AICc(area.urban)
plot(allEffects(area.urban))

res.1 <- simulateResiduals(area.urban)
testSpatialAutocorrelation(
  res.1,
  x = tapply(hbs.db.1$X, hbs.db.1$Site, mean),
  y = tapply(hbs.db.1$Y, hbs.db.1$Site, mean)) 


# 7.Test interaction habitat area * vegetated open land
area.vol <- glmmTMB(abundance.hb ~ Category + scale(Hab_Area) * scale(Vegetated.open.land) +
                      scale(Forest) + scale(Urban) +
                      (1 | City/Triplet),
                    family = nbinom2,
                    data = hbs.db.1)

summary(area.vol)
mod_dharma1 <- area.vol %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.vol) 
AICc(area.vol)
plot(allEffects(area.vol))

res.1 <- simulateResiduals(area.vol)
testSpatialAutocorrelation(
  res.1,
  x = tapply(hbs.db.1$X, hbs.db.1$Site, mean),
  y = tapply(hbs.db.1$Y, hbs.db.1$Site, mean)) 

# 8. Test interaction habitat area * forest
area.forest <- glmmTMB(abundance.hb ~ Category + scale(Hab_Area) * scale(Forest) +
                         scale(Vegetated.open.land) + scale(Urban) +
                         (1 | City/Triplet),
                       family = nbinom2,
                       data = hbs.db.1)

summary(area.forest)
mod_dharma1 <- area.forest %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.forest) 
AICc(area.forest)
plot(allEffects(area.forest))

res.1 <- simulateResiduals(area.forest)
testSpatialAutocorrelation(
  res.1,
  x = tapply(hbs.db.1$X, hbs.db.1$Site, mean),
  y = tapply(hbs.db.1$Y, hbs.db.1$Site, mean)) 

# Collect AICc values

mods <- list(
  base = base.hbs.1,
  cat.area = cat.area, 
  cat.urban = cat.urban,
  cat.open = cat.vol,
  cat.forest = cat.forest,
  area.urban = area.urban, 
  area.vol = area.vol, 
  area.forest = area.forest
)

map_dfr(mods, AIC, .id = "model")
map_dfr(mods, AICc, .id = "model")


aic_table <- map_dfr(
  mods,
  ~ data.frame(AICc = AICc(.x)),
  .id = "model"
)

aic_table <- aic_table %>%
  arrange(AICc) %>%
  mutate(delta_AICc = AICc - min(AICc))

aic_table
