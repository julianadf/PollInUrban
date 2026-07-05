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

# Create database for transect data and surrounding landscape ----
landscape.100 <- read.csv2(here("data", "100m.csv" ), fileEncoding = "latin1")
landscape.1 <- read.csv2(here("data", "1km.csv" ), fileEncoding = "latin1")
generalists.specialists.raw <- read.csv2(here("data/databases", "transects.specialists.generalists.csv" ), fileEncoding = "latin1")

# Butterflies
all_sites <- generalists.specialists.raw %>%
  distinct(City, Site, Triplet, Category, honeybee.abund)

butterflies <- generalists.specialists.raw |> 
  filter(Pollinator.type == "Butterfly") |> 
  filter(Specialist != "", !is.na(Specialist)) |> 
  group_by(City, Site, Triplet, Category, Specialisation, honeybee.abund) %>%
  summarise(
    abundance.group = sum(Abund),
    .groups = "drop") |> 
  complete(nesting(City, Site, Triplet, Category, honeybee.abund),
           Specialisation = c("Oligophagous", "Polyphagous"),
           fill = list(abundance.group = 0))
butterflies


butterflies <- all_sites %>%
  crossing(Specialisation = c("Oligophagous", "Polyphagous")) %>%
  left_join(butterflies, by = c("City", "Site", "Triplet", "Category","Specialisation", "honeybee.abund")) %>%
  mutate(abundance.group = replace_na(abundance.group, 0))


# merge data
bf.100m <- butterflies %>%
  left_join(landscape.100, by = c("City", "Site", "Triplet", "Category"))

bf.1km <- butterflies %>%
  left_join(landscape.1, by = c("City", "Site", "Triplet", "Category"))

# Analyse 100 m:
base.butterfly.100m <- glmmTMB(abundance.group ~ Category + Specialisation + scale(honeybee.abund) + 
                                 scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                                 (1 | City/Triplet),
                               family = nbinom2,
                               data = bf.100m)

summary(base.butterfly.100m)
mod_dharma1 <- base.butterfly.100m %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(base.butterfly.100m) 
AICc(base.butterfly.100m) 
plot(allEffects(base.butterfly.100m))

res.100 <- simulateResiduals(base.butterfly.100m)
res.agg <- recalculateResiduals(res.100,  group = bf.100m$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bf.100m$X, bf.100m$Site, mean),
  y = tapply(bf.100m$Y, bf.100m$Site, mean))

# Category * specialisation
cat.butterfly.100m <- glmmTMB(abundance.group ~ Category * Specialisation + scale(honeybee.abund) + 
                                scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                                (1 | City/Triplet),
                              family = nbinom2,
                              data = bf.100m)

summary(cat.butterfly.100m)
mod_dharma1 <- cat.butterfly.100m %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.butterfly.100m) 
AICc(cat.butterfly.100m) 
plot(allEffects(cat.butterfly.100m))

res.100 <- simulateResiduals(cat.butterfly.100m)
res.agg <- recalculateResiduals(res.100,  group = bf.100m$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bf.100m$X, bf.100m$Site, mean),
  y = tapply(bf.100m$Y, bf.100m$Site, mean))


# specialist * honey bee abundance
honey.butterfly.100m <- glmmTMB(abundance.group ~ Category + Specialisation * scale(honeybee.abund) + 
                                  scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                                  (1 | City/Triplet),
                                family = nbinom2,
                                data = bf.100m)

summary(honey.butterfly.100m)
mod_dharma1 <- honey.butterfly.100m %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(honey.butterfly.100m) 
AICc(honey.butterfly.100m) 
plot(allEffects(honey.butterfly.100m))

res.100 <- simulateResiduals(honey.butterfly.100m)
res.agg <- recalculateResiduals(res.100,  group = bf.100m$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bf.100m$X, bf.100m$Site, mean),
  y = tapply(bf.100m$Y, bf.100m$Site, mean))


# specialist * habitat area
area.butterfly.100m <- glmmTMB(abundance.group ~ Category + Specialisation * scale(Hab_Area) + scale(honeybee.abund) + 
                                 scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                                 (1 | City/Triplet),
                               family = nbinom2,
                               data = bf.100m)

summary(area.butterfly.100m)
mod_dharma1 <- area.butterfly.100m %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.butterfly.100m) 
AICc(area.butterfly.100m) 
plot(allEffects(area.butterfly.100m))

res.100 <- simulateResiduals(area.butterfly.100m)
res.agg <- recalculateResiduals(res.100,  group = bf.100m$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bf.100m$X, bf.100m$Site, mean),
  y = tapply(bf.100m$Y, bf.100m$Site, mean))

# specialist * urban
urban.butterfly.100m <- glmmTMB(abundance.group ~ Category + Specialisation * scale(Urban) +
                                  scale(Hab_Area) + scale(honeybee.abund) + 
                                  scale(Vegetated.open.land) + scale(Forest) + 
                                  (1 | City/Triplet),
                                family = nbinom2,
                                data = bf.100m)

summary(urban.butterfly.100m)
mod_dharma1 <- urban.butterfly.100m %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(urban.butterfly.100m) 
AICc(urban.butterfly.100m) 
plot(allEffects(urban.butterfly.100m))

res.100 <- simulateResiduals(urban.butterfly.100m)
res.agg <- recalculateResiduals(res.100,  group = bf.100m$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bf.100m$X, bf.100m$Site, mean),
  y = tapply(bf.100m$Y, bf.100m$Site, mean))

# specialist * vegetated open land
vol.butterfly.100m <- glmmTMB(abundance.group ~ Category + Specialisation *  scale(Vegetated.open.land) +
                                scale(Urban) + scale(Hab_Area) + scale(honeybee.abund) + 
                                scale(Forest) + 
                                (1 | City/Triplet),
                              family = nbinom2,
                              data = bf.100m)

summary(vol.butterfly.100m)
mod_dharma1 <- vol.butterfly.100m %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(vol.butterfly.100m) 
AICc(vol.butterfly.100m) 
plot(allEffects(vol.butterfly.100m))

res.100 <- simulateResiduals(vol.butterfly.100m)
res.agg <- recalculateResiduals(res.100,  group = bf.100m$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bf.100m$X, bf.100m$Site, mean),
  y = tapply(bf.100m$Y, bf.100m$Site, mean))

# specialist * forest
forest.butterfly.100m <- glmmTMB(abundance.group ~ Category + Specialisation *  scale(Forest) +
                                   scale(Vegetated.open.land) + scale(Urban) + scale(Hab_Area) + scale(honeybee.abund) + 
                                   (1 | City/Triplet),
                                 family = nbinom2,
                                 data = bf.100m)

summary(forest.butterfly.100m)
mod_dharma1 <- forest.butterfly.100m %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(forest.butterfly.100m) 
AICc(forest.butterfly.100m) 
plot(allEffects(forest.butterfly.100m))

res.100 <- simulateResiduals(forest.butterfly.100m)
res.agg <- recalculateResiduals(res.100,  group = bf.100m$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bf.100m$X, bf.100m$Site, mean),
  y = tapply(bf.100m$Y, bf.100m$Site, mean))

# Collect AICc values

mods <- list(
  base = base.butterfly.100m,
  category = cat.butterfly.100m,
  honey = honey.butterfly.100m, 
  area = area.butterfly.100m, 
  urban = urban.butterfly.100m,
  open = vol.butterfly.100m,
  forest = forest.butterfly.100m
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

# Analyse 1 km:
base.butterfly.1km <- glmmTMB(abundance.group ~ Category + Specialisation + scale(honeybee.abund) + 
                                scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                                (1 | City/Triplet),
                              family = nbinom2,
                              data = bf.1km)

summary(base.butterfly.1km)
mod_dharma1 <- base.butterfly.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(base.butterfly.1km) 
AICc(base.butterfly.1km) 
plot(allEffects(base.butterfly.1km))

res.1 <- simulateResiduals(base.butterfly.1km)
res.agg <- recalculateResiduals(res.1,  group = bf.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bf.1km$X, bf.1km$Site, mean),
  y = tapply(bf.1km$Y, bf.1km$Site, mean))

# Category * specialisation
cat.butterfly.1km <- glmmTMB(abundance.group ~ Category * Specialisation + scale(honeybee.abund) + 
                               scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                               (1 | City/Triplet),
                             family = nbinom2,
                             data = bf.1km)

summary(cat.butterfly.1km)
mod_dharma1 <- cat.butterfly.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.butterfly.1km) 
AICc(cat.butterfly.1km) 
plot(allEffects(cat.butterfly.1km))

res.1 <- simulateResiduals(cat.butterfly.1km)
res.agg <- recalculateResiduals(res.1,  group = bf.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bf.1km$X, bf.1km$Site, mean),
  y = tapply(bf.1km$Y, bf.1km$Site, mean))


# specialist * honey bee abundance
honey.butterfly.1km <- glmmTMB(abundance.group ~ Category + Specialisation * scale(honeybee.abund) + 
                                 scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                                 (1 | City/Triplet),
                               family = nbinom2,
                               data = bf.1km)

summary(honey.butterfly.1km)
mod_dharma1 <- honey.butterfly.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(honey.butterfly.1km) 
AICc(honey.butterfly.1km) 
plot(allEffects(honey.butterfly.1km))

res.1 <- simulateResiduals(honey.butterfly.1km)
res.agg <- recalculateResiduals(res.1,  group = bf.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bf.1km$X, bf.1km$Site, mean),
  y = tapply(bf.1km$Y, bf.1km$Site, mean))


# specialist * habitat area
area.butterfly.1km <- glmmTMB(abundance.group ~ Category + Specialisation * scale(Hab_Area) + scale(honeybee.abund) + 
                                scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                                (1 | City/Triplet),
                              family = nbinom2,
                              data = bf.1km)

summary(area.butterfly.1km)
mod_dharma1 <- area.butterfly.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.butterfly.1km) 
AICc(area.butterfly.1km) 
plot(allEffects(area.butterfly.1km))

res.1 <- simulateResiduals(area.butterfly.1km)
res.agg <- recalculateResiduals(res.1,  group = bf.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bf.1km$X, bf.1km$Site, mean),
  y = tapply(bf.1km$Y, bf.1km$Site, mean))

# specialist * urban
urban.butterfly.1km <- glmmTMB(abundance.group ~ Category + Specialisation * scale(Urban) +
                                 scale(Hab_Area) + scale(honeybee.abund) + 
                                 scale(Vegetated.open.land) + scale(Forest) + 
                                 (1 | City/Triplet),
                               family = nbinom2,
                               data = bf.1km)

summary(urban.butterfly.1km)
mod_dharma1 <- urban.butterfly.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(urban.butterfly.1km) 
AICc(urban.butterfly.1km) 
plot(allEffects(urban.butterfly.1km))

res.1 <- simulateResiduals(urban.butterfly.1km)
res.agg <- recalculateResiduals(res.1,  group = bf.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bf.1km$X, bf.1km$Site, mean),
  y = tapply(bf.1km$Y, bf.1km$Site, mean))

# specialist * vegetated open land
vol.butterfly.1km <- glmmTMB(abundance.group ~ Category + Specialisation *  scale(Vegetated.open.land) +
                               scale(Urban) + scale(Hab_Area) + scale(honeybee.abund) + 
                               scale(Forest) + 
                               (1 | City/Triplet),
                             family = nbinom2,
                             data = bf.1km)

summary(vol.butterfly.1km)
mod_dharma1 <- vol.butterfly.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(vol.butterfly.1km) 
AICc(vol.butterfly.1km) 
plot(allEffects(vol.butterfly.1km))

res.1 <- simulateResiduals(vol.butterfly.1km)
res.agg <- recalculateResiduals(res.1,  group = bf.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bf.1km$X, bf.1km$Site, mean),
  y = tapply(bf.1km$Y, bf.1km$Site, mean))

# specialist * forest
forest.butterfly.1km <- glmmTMB(abundance.group ~ Category + Specialisation *  scale(Forest) +
                                  scale(Vegetated.open.land) + scale(Urban) + scale(Hab_Area) + scale(honeybee.abund) + 
                                  (1 | City/Triplet),
                                family = nbinom2,
                                data = bf.1km)

summary(forest.butterfly.1km)
mod_dharma1 <- forest.butterfly.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(forest.butterfly.1km) 
AICc(forest.butterfly.1km) 
plot(allEffects(forest.butterfly.1km))

res.1 <- simulateResiduals(forest.butterfly.1km)
res.agg <- recalculateResiduals(res.1,  group = bf.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bf.1km$X, bf.1km$Site, mean),
  y = tapply(bf.1km$Y, bf.1km$Site, mean))

# Collect AICc values

mods <- list(
  base = base.butterfly.1km,
  category = cat.butterfly.1km,
  honey = honey.butterfly.1km, 
  area = area.butterfly.1km, 
  urban = urban.butterfly.1km,
  open = vol.butterfly.1km,
  forest = forest.butterfly.1km
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


# Wild bees
all_sites <- generalists.specialists.raw %>%
  distinct(City, Site, Triplet, Category, honeybee.abund) 


bees <- generalists.specialists.raw |> 
  filter(Pollinator.type == "Bee") |> 
  filter(Specialist != "", !is.na(Specialist)) |> 
  filter(Specialisation != "Monolectic") |> #ignores the one monolectic observation
  group_by(City, Site, Triplet, Category, Specialisation, honeybee.abund) %>%
  summarise(
    abundance.group = sum(Abund),
    .groups = "drop") |> 
  complete(nesting(City, Site, Triplet, Category, honeybee.abund),
           Specialisation = c("Oligolectic", "Polylectic"), 
           fill = list(abundance.group = 0))
bees


bees <- all_sites %>%
  crossing(Specialisation = c("Oligolectic", "Polylectic")) %>%
  left_join(bees, by = c("City", "Site", "Triplet", "Category","Specialisation", "honeybee.abund")) %>%
  mutate(abundance.group = replace_na(abundance.group, 0))


# merge data
bees.100m <- bees %>%
  left_join(landscape.100, by = c("City", "Site", "Triplet", "Category"))

bees.1km <- bees %>%
  left_join(landscape.1, by = c("City", "Site", "Triplet", "Category"))

# Analyse 100 m:
base.bees.100m <- glmmTMB(abundance.group ~ Category + Specialisation + scale(honeybee.abund) + 
                            scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                            (1 | City/Triplet),
                          family = nbinom2,
                          data = bees.100m)

summary(base.bees.100m)
mod_dharma1 <- base.bees.100m %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(base.bees.100m) 
AICc(base.bees.100m) 
plot(allEffects(base.bees.100m))

res.100 <- simulateResiduals(base.bees.100m)
res.agg <- recalculateResiduals(res.100,  group = bees.100m$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.100m$X, bees.100m$Site, mean),
  y = tapply(bees.100m$Y, bees.100m$Site, mean))

# 1. Category * specialisation
cat.bees.100m <- glmmTMB(abundance.group ~ Category * Specialisation + scale(honeybee.abund) + 
                           scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                           (1 | City/Triplet),
                         family = nbinom2,
                         data = bees.100m)

summary(cat.bees.100m)
mod_dharma1 <- cat.bees.100m %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.bees.100m) 
AICc(cat.bees.100m) 
plot(allEffects(cat.bees.100m))

res.100 <- simulateResiduals(cat.bees.100m)
res.agg <- recalculateResiduals(res.100,  group = bees.100m$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.100m$X, bees.100m$Site, mean),
  y = tapply(bees.100m$Y, bees.100m$Site, mean))


# 2. specialist * honey bee abundance
honey.bees.100m <- glmmTMB(abundance.group ~ Category + Specialisation * scale(honeybee.abund) + 
                             scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                             (1 | City/Triplet),
                           family = nbinom2,
                           data = bees.100m)

summary(honey.bees.100m)
mod_dharma1 <- honey.bees.100m %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(honey.bees.100m) 
AICc(honey.bees.100m) 
plot(allEffects(honey.bees.100m))

res.100 <- simulateResiduals(honey.bees.100m)
res.agg <- recalculateResiduals(res.100,  group = bees.100m$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.100m$X, bees.100m$Site, mean),
  y = tapply(bees.100m$Y, bees.100m$Site, mean))

# 3. specialist * habitat area
area.bees.100m <- glmmTMB(abundance.group ~ Category + Specialisation * scale(Hab_Area) + scale(honeybee.abund) + 
                            scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                            (1 | City/Triplet),
                          family = nbinom2,
                          data = bees.100m)

summary(area.bees.100m)
mod_dharma1 <- area.bees.100m %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.bees.100m) 
AICc(area.bees.100m) 
plot(allEffects(area.bees.100m))

res.100 <- simulateResiduals(area.bees.100m)
res.agg <- recalculateResiduals(res.100,  group = bees.100m$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.100m$X, bees.100m$Site, mean),
  y = tapply(bees.100m$Y, bees.100m$Site, mean))

# 4- specialist * urban
urban.bees.100m <- glmmTMB(abundance.group ~ Category + Specialisation * scale(Urban) +
                             scale(Hab_Area) + scale(honeybee.abund) + 
                             scale(Vegetated.open.land) + scale(Forest) + 
                             (1 | City/Triplet),
                           family = nbinom2,
                           data = bees.100m)

summary(urban.bees.100m)
mod_dharma1 <- urban.bees.100m %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(urban.bees.100m) 
AICc(urban.bees.100m) 
plot(allEffects(urban.bees.100m))

res.100 <- simulateResiduals(urban.bees.100m)
res.agg <- recalculateResiduals(res.100,  group = bees.100m$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.100m$X, bees.100m$Site, mean),
  y = tapply(bees.100m$Y, bees.100m$Site, mean))

# 5. specialist * vegetated open land
vol.bees.100m <- glmmTMB(abundance.group ~ Category + Specialisation *  scale(Vegetated.open.land) +
                           scale(Urban) + scale(Hab_Area) + scale(honeybee.abund) + 
                           scale(Forest) + 
                           (1 | City/Triplet),
                         family = nbinom2,
                         data = bees.100m)

summary(vol.bees.100m)
mod_dharma1 <- vol.bees.100m %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(vol.bees.100m) 
AICc(vol.bees.100m) 
plot(allEffects(vol.bees.100m))

res.100 <- simulateResiduals(vol.bees.100m)
res.agg <- recalculateResiduals(res.100,  group = bees.100m$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.100m$X, bees.100m$Site, mean),
  y = tapply(bees.100m$Y, bees.100m$Site, mean))

# 6. specialist * forest
forest.bees.100m <- glmmTMB(abundance.group ~ Category + Specialisation *  scale(Forest) +
                              scale(Vegetated.open.land) + scale(Urban) + scale(Hab_Area) + scale(honeybee.abund) + 
                              (1 | City/Triplet),
                            family = nbinom2,
                            data = bees.100m)

summary(forest.bees.100m)
mod_dharma1 <- forest.bees.100m %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(forest.bees.100m) 
AICc(forest.bees.100m) 
plot(allEffects(forest.bees.100m))

res.100 <- simulateResiduals(forest.bees.100m)
res.agg <- recalculateResiduals(res.100,  group = bees.100m$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.100m$X, bees.100m$Site, mean),
  y = tapply(bees.100m$Y, bees.100m$Site, mean))

# Collect AICc values

mods <- list(
  base = base.bees.100m,
  category = cat.bees.100m,
  honey = honey.bees.100m, 
  area = area.bees.100m, 
  urban = urban.bees.100m,
  open = vol.bees.100m,
  forest = forest.bees.100m
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

# Analyse 1 km:
base.bees.1km <- glmmTMB(abundance.group ~ Category + Specialisation + scale(honeybee.abund) + 
                           scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                           (1 | City/Triplet),
                         family = nbinom2,
                         data = bees.1km)

summary(base.bees.1km)
mod_dharma1 <- base.bees.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(base.bees.1km) 
AICc(base.bees.1km) 
plot(allEffects(base.bees.1km))

res.1 <- simulateResiduals(base.bees.1km)
res.agg <- recalculateResiduals(res.1,  group = bees.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.1km$X, bees.1km$Site, mean),
  y = tapply(bees.1km$Y, bees.1km$Site, mean))

# Category * specialisation
cat.bees.1km <- glmmTMB(abundance.group ~ Category * Specialisation + scale(honeybee.abund) + 
                          scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                          (1 | City/Triplet),
                        family = nbinom2,
                        data = bees.1km)

summary(cat.bees.1km)
mod_dharma1 <- cat.bees.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.bees.1km) 
AICc(cat.bees.1km) 
plot(allEffects(cat.bees.1km))

res.1 <- simulateResiduals(cat.bees.1km)
res.agg <- recalculateResiduals(res.1,  group = bees.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.1km$X, bees.1km$Site, mean),
  y = tapply(bees.1km$Y, bees.1km$Site, mean))


# specialist * honey bee abundance
honey.bees.1km <- glmmTMB(abundance.group ~ Category + Specialisation * scale(honeybee.abund) + 
                            scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                            (1 | City/Triplet),
                          family = nbinom2,
                          data = bees.1km)

summary(honey.bees.1km)
mod_dharma1 <- honey.bees.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(honey.bees.1km) 
AICc(honey.bees.1km) 
plot(allEffects(honey.bees.1km))

res.1 <- simulateResiduals(honey.bees.1km)
res.agg <- recalculateResiduals(res.1,  group = bees.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.1km$X, bees.1km$Site, mean),
  y = tapply(bees.1km$Y, bees.1km$Site, mean))


# specialist * habitat area
area.bees.1km <- glmmTMB(abundance.group ~ Category + Specialisation * scale(Hab_Area) + scale(honeybee.abund) + 
                           scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                           (1 | City/Triplet),
                         family = nbinom2,
                         data = bees.1km)

summary(area.bees.1km)
mod_dharma1 <- area.bees.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.bees.1km) 
AICc(area.bees.1km) 
plot(allEffects(area.bees.1km))

res.1 <- simulateResiduals(area.bees.1km)
res.agg <- recalculateResiduals(res.1,  group = bees.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.1km$X, bees.1km$Site, mean),
  y = tapply(bees.1km$Y, bees.1km$Site, mean))

# specialist * urban
urban.bees.1km <- glmmTMB(abundance.group ~ Category + Specialisation * scale(Urban) +
                            scale(Hab_Area) + scale(honeybee.abund) + 
                            scale(Vegetated.open.land) + scale(Forest) + 
                            (1 | City/Triplet),
                          family = nbinom2,
                          data = bees.1km)

summary(urban.bees.1km)
mod_dharma1 <- urban.bees.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(urban.bees.1km) 
AICc(urban.bees.1km) 
plot(allEffects(urban.bees.1km))

res.1 <- simulateResiduals(urban.bees.1km)
res.agg <- recalculateResiduals(res.1,  group = bees.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.1km$X, bees.1km$Site, mean),
  y = tapply(bees.1km$Y, bees.1km$Site, mean))

# specialist * vegetated open land
vol.bees.1km <- glmmTMB(abundance.group ~ Category + Specialisation *  scale(Vegetated.open.land) +
                          scale(Urban) + scale(Hab_Area) + scale(honeybee.abund) + 
                          scale(Forest) + 
                          (1 | City/Triplet),
                        family = nbinom2,
                        data = bees.1km)

summary(vol.bees.1km)
mod_dharma1 <- vol.bees.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(vol.bees.1km) 
AICc(vol.bees.1km) 
plot(allEffects(vol.bees.1km))

res.1 <- simulateResiduals(vol.bees.1km)
res.agg <- recalculateResiduals(res.1,  group = bees.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.1km$X, bees.1km$Site, mean),
  y = tapply(bees.1km$Y, bees.1km$Site, mean))

# specialist * forest
forest.bees.1km <- glmmTMB(abundance.group ~ Category + Specialisation *  scale(Forest) +
                             scale(Vegetated.open.land) + scale(Urban) + scale(Hab_Area) + scale(honeybee.abund) + 
                             (1 | City/Triplet),
                           family = nbinom2,
                           data = bees.1km)

summary(forest.bees.1km)
mod_dharma1 <- forest.bees.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(forest.bees.1km) 
AICc(forest.bees.1km) 
plot(allEffects(forest.bees.1km))

res.1 <- simulateResiduals(forest.bees.1km)
res.agg <- recalculateResiduals(res.1,  group = bees.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.1km$X, bees.1km$Site, mean),
  y = tapply(bees.1km$Y, bees.1km$Site, mean))

# area * urban
area.urban.bees.1km <- glmmTMB(abundance.group ~ Category * Specialisation +  scale(Hab_Area) * scale(Urban) +
                                 scale(Forest) + scale(Vegetated.open.land) + scale(honeybee.abund) + 
                                 (1 | City/Triplet),
                               family = nbinom2,
                               data = bees.1km)

summary(area.urban.bees.1km)
mod_dharma1 <- area.urban.bees.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.urban.bees.1km) 
AICc(area.urban.bees.1km) 
plot(allEffects(area.urban.bees.1km))

res.1 <- simulateResiduals(area.urban.bees.1km)
res.agg <- recalculateResiduals(res.1,  group = bees.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.1km$X, bees.1km$Site, mean),
  y = tapply(bees.1km$Y, bees.1km$Site, mean))

# area * vegetated open land
area.vol.bees.1km <- glmmTMB(abundance.group ~ Category + Specialisation +  scale(Hab_Area) * scale(Vegetated.open.land) +
                               scale(Urban) + scale(Forest) + scale(honeybee.abund) + 
                               (1 | City/Triplet),
                             family = nbinom2,
                             data = bees.1km)

summary(area.vol.bees.1km)
mod_dharma1 <- area.vol.bees.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.vol.bees.1km) 
AICc(area.vol.bees.1km) 
plot(allEffects(area.vol.bees.1km))

res.1 <- simulateResiduals(area.vol.bees.1km)
res.agg <- recalculateResiduals(res.1,  group = bees.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.1km$X, bees.1km$Site, mean),
  y = tapply(bees.1km$Y, bees.1km$Site, mean))

# area * forest
area.forest.bees.1km <- glmmTMB(abundance.group ~ Category + Specialisation +  scale(Hab_Area) * scale(Forest) +
                                  scale(Vegetated.open.land) + scale(Urban) + scale(honeybee.abund) + 
                                  (1 | City/Triplet),
                                family = nbinom2,
                                data = bees.1km)

summary(area.forest.bees.1km)
mod_dharma1 <- area.forest.bees.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.forest.bees.1km) 
AICc(area.forest.bees.1km) 
plot(allEffects(area.forest.bees.1km))

res.1 <- simulateResiduals(area.forest.bees.1km)
res.agg <- recalculateResiduals(res.1,  group = bees.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.1km$X, bees.1km$Site, mean),
  y = tapply(bees.1km$Y, bees.1km$Site, mean))

# Collect AICc values

mods <- list(
  base = base.bees.1km,
  category = cat.bees.1km,
  honey = honey.bees.1km, 
  area = area.bees.1km, 
  urban = urban.bees.1km,
  open = vol.bees.1km,
  forest = forest.bees.1km, 
  area.urban = area.urban.bees.1km, 
  area.vol = area.vol.bees.1km, 
  area.forest = area.forest.bees.1km
  
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
