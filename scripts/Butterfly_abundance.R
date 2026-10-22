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

# Butterflies
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

# Butterfly abundance per site
bf_abund <- abund.100.raw %>%
  mutate(Pollinator.type = ifelse(Pollinator.type == "Burnet moth", "Butterfly", Pollinator.type)) %>%
  filter(Pollinator.type == "Butterfly") %>% 
  group_by(City, Site, Triplet, Category, X, Y, Hab_Area, Arable, Forest, Non.vegetated.open.land,
    Urban, Vegetated.open.land, Water, Wetland) %>%
  summarise(abundance.bf = sum(total_abundance), .groups = "drop")

# Join butterfly abundance to all sites
abund.100.bfs <- all_sites %>%
  left_join(bf_abund,
    by = c(
      "City", "Site", "Triplet", "Category",
      "X", "Y", "Hab_Area",
      "Arable", "Forest",
      "Non.vegetated.open.land",
      "Urban",
      "Vegetated.open.land",
      "Water",
      "Wetland")) %>%
  mutate(abundance.bf = replace_na(abundance.bf, 0))

# Add honey bee abundance
bfs.db.100 <- abund.100.bfs %>%
  left_join(honeybee_df, by = "Site")

bfs.db.100 <- bfs.db.100 |> 
  mutate(City = as.factor(City)) |> 
  mutate(Site = as.factor(Site)) |> 
  mutate(Triplet = as.factor(Triplet)) |> 
  mutate(Category = as.factor(Category)) 

# Landscape analysis - 100 m: ----
# 1. base model
base.bfs.100 <- glmmTMB(abundance.bf ~ Category + scale(honeybee) + scale(Hab_Area) + 
                           scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                           (1 | City/Triplet),
                         family = nbinom2,
                         data = bfs.db.100)

summary(base.bfs.100)
mod_dharma1 <- base.bfs.100 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(base.bfs.100) 
AICc(base.bfs.100)
plot(allEffects(base.bfs.100))

res.100 <- simulateResiduals(base.bfs.100)
res.agg <- recalculateResiduals(res.100,  group = bfs.db.100$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bfs.db.100$X, bfs.db.100$Site, mean),
  y = tapply(bfs.db.100$Y, bfs.db.100$Site, mean))


# 2. Test interaction category * honeybee abundance
cat.honey <- glmmTMB(abundance.bf ~ Category * scale(honeybee) + scale(Hab_Area) + 
                       scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                       (1 | City/Triplet),
                     family = nbinom2,
                     data = bfs.db.100)

summary(cat.honey)
mod_dharma1 <- cat.honey %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.honey) 
AICc(cat.honey)
plot(allEffects(cat.honey))

res.100 <- simulateResiduals(cat.honey)

testSpatialAutocorrelation(
  res.100,
  x = tapply(bfs.db.100$X, bfs.db.100$Site, mean),
  y = tapply(bfs.db.100$Y, bfs.db.100$Site, mean))

# 3. Test interaction category * habitat area
cat.area <- glmmTMB(abundance.bf ~ Category * scale(Hab_Area) + scale(honeybee) +   
                      scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                      (1 | City/Triplet),
                    family = nbinom2,
                    data = bfs.db.100)

summary(cat.area)
mod_dharma1 <- cat.area %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.area) 
AICc(cat.area)
plot(allEffects(cat.area))

res.100 <- simulateResiduals(cat.area)

testSpatialAutocorrelation(
  res.100,
  x = tapply(bfs.db.100$X, bfs.db.100$Site, mean),
  y = tapply(bfs.db.100$Y, bfs.db.100$Site, mean))

# 4. Test category  * urban
cat.urban <- glmmTMB(abundance.bf ~ Category * scale(Urban) +
                       scale(Hab_Area) + scale(honeybee) +   
                       scale(Vegetated.open.land) + scale(Forest) + 
                       (1 | City/Triplet),
                     family = nbinom2,
                     data = bfs.db.100)

summary(cat.urban)
mod_dharma1 <- cat.urban %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.urban) 
AICc(cat.urban)
plot(allEffects(cat.urban))

res.100 <- simulateResiduals(cat.urban)

testSpatialAutocorrelation(
  res.100,
  x = tapply(bfs.db.100$X, bfs.db.100$Site, mean),
  y = tapply(bfs.db.100$Y, bfs.db.100$Site, mean))

# 5. Test interaction category * vegetated open land
cat.vol <- glmmTMB(abundance.bf ~ Category * scale(Vegetated.open.land) +
                     scale(Urban) + scale(Hab_Area) + scale(honeybee) + scale(Forest) + 
                     (1 | City/Triplet),
                   family = nbinom2,
                   data = bfs.db.100)

summary(cat.vol)
mod_dharma1 <- cat.vol %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.vol) 
AICc(cat.vol)
plot(allEffects(cat.vol))

res.100 <- simulateResiduals(cat.vol)

testSpatialAutocorrelation(
  res.100,
  x = tapply(bfs.db.100$X, bfs.db.100$Site, mean),
  y = tapply(bfs.db.100$Y, bfs.db.100$Site, mean))

# 6. Test interaction category * forest
cat.forest <- glmmTMB(abundance.bf ~ Category * scale(Forest) +
                        scale(Vegetated.open.land) + scale(Urban) + scale(Hab_Area) + scale(honeybee) + 
                        (1 | City/Triplet),
                      family = nbinom2,
                      data = bfs.db.100)

summary(cat.forest)
mod_dharma1 <- cat.forest %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.forest) 
AICc(cat.forest)
plot(allEffects(cat.forest))

res.100 <- simulateResiduals(cat.forest)

testSpatialAutocorrelation(
  res.100,
  x = tapply(bfs.db.100$X, bfs.db.100$Site, mean),
  y = tapply(bfs.db.100$Y, bfs.db.100$Site, mean))

# 7. Test interaction habitat area * urban
area.urban <- glmmTMB(abundance.bf ~ Category + scale(Hab_Area) * scale(Urban) +
                        scale(Forest) + scale(Vegetated.open.land) +  scale(honeybee) + 
                        (1 | City/Triplet),
                      family = nbinom2,
                      data = bfs.db.100)

summary(area.urban)
mod_dharma1 <- area.urban %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.urban) 
AICc(area.urban)
plot(allEffects(area.urban))

res.100 <- simulateResiduals(area.urban)

testSpatialAutocorrelation(
  res.100,
  x = tapply(bfs.db.100$X, bfs.db.100$Site, mean),
  y = tapply(bfs.db.100$Y, bfs.db.100$Site, mean))


# 8.Test interaction habitat area * vegetated open land
area.vol <- glmmTMB(abundance.bf ~ Category + scale(Hab_Area) * scale(Vegetated.open.land) +
                      scale(Forest) + scale(honeybee) + scale(Urban) +
                      (1 | City/Triplet),
                    family = nbinom2,
                    data = bfs.db.100)

summary(area.vol)
mod_dharma1 <- area.vol %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.vol) 
AICc(area.vol)
plot(allEffects(area.vol))

res.100 <- simulateResiduals(area.vol)

testSpatialAutocorrelation(
  res.100,
  x = tapply(bfs.db.100$X, bfs.db.100$Site, mean),
  y = tapply(bfs.db.100$Y, bfs.db.100$Site, mean))

# 9. Test interaction habitat area * forest
area.forest <- glmmTMB(abundance.bf ~ Category + scale(Hab_Area) * scale(Forest) +
                         scale(Vegetated.open.land) + scale(honeybee) + scale(Urban) +
                         (1 | City/Triplet),
                       family = nbinom2,
                       data = bfs.db.100)

summary(area.forest)
mod_dharma1 <- area.forest %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.forest) 
AICc(area.forest)
plot(allEffects(area.forest))

res.100 <- simulateResiduals(area.forest)

testSpatialAutocorrelation(
  res.100,
  x = tapply(bfs.db.100$X, bfs.db.100$Site, mean),
  y = tapply(bfs.db.100$Y, bfs.db.100$Site, mean))

# 10. Null model
null.bfs.100 <- glmmTMB(abundance.bf ~ 1 + 
                          (1 | City/Triplet),
                        family = nbinom2,
                        data = bfs.db.100)

# Collect AICc values

mods <- list(
  base = base.bfs.100,
  cat.honey = cat.honey, 
  cat.area = cat.area, 
  cat.urban = cat.urban,
  cat.open = cat.vol,
  cat.forest = cat.forest,
  area.urban = area.urban, 
  area.vol = area.vol, 
  area.forest = area.forest, 
  null = null.bfs.100
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

# Butterfly abundance per site
bf_abund <- abund.1.raw %>%
  mutate(Pollinator.type = ifelse(Pollinator.type == "Burnet moth", "Butterfly", Pollinator.type)) %>%
  filter(Pollinator.type == "Butterfly") %>% 
  group_by(City, Site, Triplet, Category, X, Y, Hab_Area, Arable, Forest, Non.vegetated.open.land,
           Urban, Vegetated.open.land, Water, Wetland) %>%
  summarise(abundance.bf = sum(total_abundance), .groups = "drop")

# Join butterfly abundance to all sites
abund.1.bfs <- all_sites %>%
  left_join(bf_abund,
            by = c(
              "City", "Site", "Triplet", "Category",
              "X", "Y", "Hab_Area",
              "Arable", "Forest",
              "Non.vegetated.open.land",
              "Urban",
              "Vegetated.open.land",
              "Water",
              "Wetland")) %>%
  mutate(abundance.bf = replace_na(abundance.bf, 0))

# Add honey bee abundance
bfs.db.1 <- abund.1.bfs %>%
  left_join(honeybee_df, by = "Site")

bfs.db.1 <- bfs.db.1 |> 
  mutate(City = as.factor(City)) |> 
  mutate(Site = as.factor(Site)) |> 
  mutate(Triplet = as.factor(Triplet)) |> 
  mutate(Category = as.factor(Category)) 

# Landscape analysis - 1 km: ----
# 1. base model
base.bfs.1 <- glmmTMB(abundance.bf ~ Category + scale(honeybee) + scale(Hab_Area) + 
                          scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                          (1 | City/Triplet),
                        family = nbinom2,
                        data = bfs.db.1)

summary(base.bfs.1)
mod_dharma1 <- base.bfs.1 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(base.bfs.1) 
AICc(base.bfs.1)
plot(allEffects(base.bfs.1))

res.1 <- simulateResiduals(base.bfs.1)
testSpatialAutocorrelation(
  res.1,
  x = tapply(bfs.db.1$X, bfs.db.1$Site, mean),
  y = tapply(bfs.db.1$Y, bfs.db.1$Site, mean))


# 2. Test interaction category * honeybee abundance
cat.honey <- glmmTMB(abundance.bf ~ Category * scale(honeybee) + scale(Hab_Area) + 
                       scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                       (1 | City/Triplet),
                     family = nbinom2,
                     data = bfs.db.1)

summary(cat.honey)
mod_dharma1 <- cat.honey %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.honey) 
AICc(cat.honey)
plot(allEffects(cat.honey))

res.1 <- simulateResiduals(base.bfs.1)
testSpatialAutocorrelation(
  res.1,
  x = tapply(bfs.db.1$X, bfs.db.1$Site, mean),
  y = tapply(bfs.db.1$Y, bfs.db.1$Site, mean))

# 3. Test interaction category * habitat area
cat.area <- glmmTMB(abundance.bf ~ Category * scale(Hab_Area) + scale(honeybee) +   
                      scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                      (1 | City/Triplet),
                    family = nbinom2,
                    data = bfs.db.1)

summary(cat.area)
mod_dharma1 <- cat.area %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.area) 
AICc(cat.area)
plot(allEffects(cat.area))

res.1 <- simulateResiduals(cat.area)
testSpatialAutocorrelation(
  res.1,
  x = tapply(bfs.db.1$X, bfs.db.1$Site, mean),
  y = tapply(bfs.db.1$Y, bfs.db.1$Site, mean))

# 4. Test category  * urban
cat.urban <- glmmTMB(abundance.bf ~ Category * scale(Urban) +
                       scale(Hab_Area) + scale(honeybee) +   
                       scale(Vegetated.open.land) + scale(Forest) + 
                       (1 | City/Triplet),
                     family = nbinom2,
                     data = bfs.db.1)

summary(cat.urban)
mod_dharma1 <- cat.urban %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.urban) 
AICc(cat.urban)
plot(allEffects(cat.urban))

res.1 <- simulateResiduals(cat.urban)
testSpatialAutocorrelation(
  res.1,
  x = tapply(bfs.db.1$X, bfs.db.1$Site, mean),
  y = tapply(bfs.db.1$Y, bfs.db.1$Site, mean))

# 5. Test interaction category * vegetated open land
cat.vol <- glmmTMB(abundance.bf ~ Category * scale(Vegetated.open.land) +
                     scale(Urban) + scale(Hab_Area) + scale(honeybee) + scale(Forest) + 
                     (1 | City/Triplet),
                   family = nbinom2,
                   data = bfs.db.1)

summary(cat.vol)
mod_dharma1 <- cat.vol %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.vol) 
AICc(cat.vol)
plot(allEffects(cat.vol))

res.1 <- simulateResiduals(cat.vol)
testSpatialAutocorrelation(
  res.1,
  x = tapply(bfs.db.1$X, bfs.db.1$Site, mean),
  y = tapply(bfs.db.1$Y, bfs.db.1$Site, mean))

# 6. Test interaction category * forest
cat.forest <- glmmTMB(abundance.bf ~ Category * scale(Forest) +
                        scale(Vegetated.open.land) + scale(Urban) + scale(Hab_Area) + scale(honeybee) + 
                        (1 | City/Triplet),
                      family = nbinom2,
                      data = bfs.db.1)

summary(cat.forest)
mod_dharma1 <- cat.forest %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.forest) 
AICc(cat.forest)
plot(allEffects(cat.forest))

res.1 <- simulateResiduals(cat.forest)
testSpatialAutocorrelation(
  res.1,
  x = tapply(bfs.db.1$X, bfs.db.1$Site, mean),
  y = tapply(bfs.db.1$Y, bfs.db.1$Site, mean))

# 7. Test interaction habitat area * urban
area.urban <- glmmTMB(abundance.bf ~ Category + scale(Hab_Area) * scale(Urban) +
                        scale(Forest) + scale(Vegetated.open.land) +  scale(honeybee) + 
                        (1 | City/Triplet),
                      family = nbinom2,
                      data = bfs.db.1)

summary(area.urban)
mod_dharma1 <- area.urban %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.urban) 
AICc(area.urban)
plot(allEffects(area.urban))

res.1 <- simulateResiduals(area.urban)
testSpatialAutocorrelation(
  res.1,
  x = tapply(bfs.db.1$X, bfs.db.1$Site, mean),
  y = tapply(bfs.db.1$Y, bfs.db.1$Site, mean))


# 8.Test interaction habitat area * vegetated open land
area.vol <- glmmTMB(abundance.bf ~ Category + scale(Hab_Area) * scale(Vegetated.open.land) +
                      scale(Forest) + scale(honeybee) + scale(Urban) +
                      (1 | City/Triplet),
                    family = nbinom2,
                    data = bfs.db.1)

summary(area.vol)
mod_dharma1 <- area.vol %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.vol) 
AICc(area.vol)
plot(allEffects(area.vol))

res.1 <- simulateResiduals(area.vol)
testSpatialAutocorrelation(
  res.1,
  x = tapply(bfs.db.1$X, bfs.db.1$Site, mean),
  y = tapply(bfs.db.1$Y, bfs.db.1$Site, mean))

# 9. Test interaction habitat area * forest
area.forest <- glmmTMB(abundance.bf ~ Category + scale(Hab_Area) * scale(Forest) +
                         scale(Vegetated.open.land) + scale(honeybee) + scale(Urban) +
                         (1 | City/Triplet),
                       family = nbinom2,
                       data = bfs.db.1)

summary(area.forest)
mod_dharma1 <- area.forest %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.forest) 
AICc(area.forest)
plot(allEffects(area.forest))

res.1 <- simulateResiduals(area.forest)
testSpatialAutocorrelation(
  res.1,
  x = tapply(bfs.db.1$X, bfs.db.1$Site, mean),
  y = tapply(bfs.db.1$Y, bfs.db.1$Site, mean))

# 10. Null model
null.bfs.1 <- glmmTMB(abundance.bf ~ 1 + 
                        (1 | City/Triplet),
                      family = nbinom2,
                      data = bfs.db.1)

# Collect AICc values

mods <- list(
  base = base.bfs.1,
  cat.honey = cat.honey, 
  cat.area = cat.area, 
  cat.urban = cat.urban,
  cat.open = cat.vol,
  cat.forest = cat.forest,
  area.urban = area.urban, 
  area.vol = area.vol, 
  area.forest = area.forest, 
  null = null.bfs.1
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
