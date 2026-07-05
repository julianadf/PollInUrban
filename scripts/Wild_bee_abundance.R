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

# Wild bees
# structure data
# 100 m:-----
honeybee_df <- abund.100.raw %>%
  group_by(Site) %>%
  summarise(honeybee = sum(honeybee_abundance, na.rm = TRUE), .groups = "drop")


abund.100.bees <- abund.100.raw %>%
  filter(Pollinator.type == "Bee") |> 
  group_by(City, Site, Triplet, Category, X, Y, Hab_Area, 
           Arable, Forest, Non.vegetated.open.land, Urban, Vegetated.open.land, Water, Wetland) 
abund.100.bees

bees.db.100 <- left_join(honeybee_df, abund.100.bees, by = c("Site"))
bees.db.100 <- bees.db.100 |> 
  mutate(City = as.factor(City)) |> 
  mutate(Site = as.factor(Site)) |> 
  mutate(Triplet = as.factor(Triplet)) |> 
  mutate(Category = as.factor(Category)) 

# Landscape analysis - pollinator groups 100 m: ----
# 1. base model
base.bees.100 <- glmmTMB(total_abundance ~ Category + scale(honeybee) + scale(Hab_Area) + 
                            scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                            (1 | City/Triplet),
                          family = nbinom2,
                          data = bees.db.100)

summary(base.bees.100)
mod_dharma1 <- base.bees.100 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(base.bees.100) 
AICc(base.bees.100)
plot(allEffects(base.bees.100))

res.100 <- simulateResiduals(base.bees.100)
testSpatialAutocorrelation(
  res.100,
  x = tapply(bees.db.100$X, bees.db.100$Site, mean),
  y = tapply(bees.db.100$Y, bees.db.100$Site, mean))


# 2. Test interaction category * honeybee abundance
cat.honey <- glmmTMB(total_abundance ~ Category * scale(honeybee) + scale(Hab_Area) + 
                              scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                              (1 | City/Triplet),
                            family = nbinom2,
                            data = bees.db.100)

summary(cat.honey)
mod_dharma1 <- cat.honey %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.honey) 
AICc(cat.honey)
plot(allEffects(cat.honey))

res.100 <- simulateResiduals(cat.honey)
res.agg <- recalculateResiduals(res.100,  group = bees.db.100$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.db.100$X, bees.db.100$Site, mean),
  y = tapply(bees.db.100$Y, bees.db.100$Site, mean)) # marginal autocorrelation!!!

# 3. Test interaction category * habitat area
cat.area <- glmmTMB(total_abundance ~ Category * scale(Hab_Area) + scale(honeybee) +   
                             scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                             (1 | City/Triplet),
                           family = nbinom2,
                           data = bees.db.100)

summary(cat.area)
mod_dharma1 <- cat.area %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.area) 
AICc(cat.area)
plot(allEffects(cat.area))

res.100 <- simulateResiduals(cat.area)
res.agg <- recalculateResiduals(res.100,  group = bees.db.100$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.db.100$X, bees.db.100$Site, mean),
  y = tapply(bees.db.100$Y, bees.db.100$Site, mean)) # marginal autocorrelation!!!

# 4. Test category  * urban
cat.urban <- glmmTMB(total_abundance ~ Category * scale(Urban) +
                              scale(Hab_Area) + scale(honeybee) +   
                              scale(Vegetated.open.land) + scale(Forest) + 
                              (1 | City/Triplet),
                            family = nbinom2,
                            data = bees.db.100)

summary(cat.urban)
mod_dharma1 <- cat.urban %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.urban) 
AICc(cat.urban)
plot(allEffects(cat.urban))

res.100 <- simulateResiduals(cat.urban)
res.agg <- recalculateResiduals(res.100,  group = bees.db.100$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.db.100$X, bees.db.100$Site, mean),
  y = tapply(bees.db.100$Y, bees.db.100$Site, mean))

# 5. Test interaction category * vegetated open land
cat.vol <- glmmTMB(total_abundance ~ Category * scale(Vegetated.open.land) +
                            scale(Urban) + scale(Hab_Area) + scale(honeybee) + scale(Forest) + 
                            (1 | City/Triplet),
                          family = nbinom2,
                          data = bees.db.100)

summary(cat.vol)
mod_dharma1 <- cat.vol %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.vol) 
AICc(cat.vol)
plot(allEffects(cat.vol))

res.100 <- simulateResiduals(cat.vol)
res.agg <- recalculateResiduals(res.100,  group = bees.db.100$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.db.100$X, bees.db.100$Site, mean),
  y = tapply(bees.db.100$Y, bees.db.100$Site, mean)) # autocorrelation!!!

# 6. Test interaction category * forest
cat.forest <- glmmTMB(total_abundance ~ Category * scale(Forest) +
                               scale(Vegetated.open.land) + scale(Urban) + scale(Hab_Area) + scale(honeybee) + 
                               (1 | City/Triplet),
                             family = nbinom2,
                             data = bees.db.100)

summary(cat.forest)
mod_dharma1 <- cat.forest %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.forest) 
AICc(cat.forest)
plot(allEffects(cat.forest))

res.100 <- simulateResiduals(cat.forest)
res.agg <- recalculateResiduals(res.100,  group = bees.db.100$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.db.100$X, bees.db.100$Site, mean),
  y = tapply(bees.db.100$Y, bees.db.100$Site, mean))

# 7. Test interaction habitat area * urban
area.urban <- glmmTMB(total_abundance ~ Category + scale(Hab_Area) * scale(Urban) +
                        scale(Forest) + scale(Vegetated.open.land) +  scale(honeybee) + 
                        (1 | City/Triplet),
                      family = nbinom2,
                      data = bees.db.100)

summary(area.urban)
mod_dharma1 <- area.urban %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.urban) 
AICc(area.urban)
plot(allEffects(area.urban))

res.100 <- simulateResiduals(area.urban)
res.agg <- recalculateResiduals(res.100,  group = bees.db.100$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.db.100$X, bees.db.100$Site, mean),
  y = tapply(bees.db.100$Y, bees.db.100$Site, mean))


# 8.Test interaction habitat area * vegetated open land
area.vol <- glmmTMB(total_abundance ~ Category + scale(Hab_Area) * scale(Vegetated.open.land) +
                      scale(Forest) + scale(honeybee) + scale(Urban) +
                      (1 | City/Triplet),
                    family = nbinom2,
                    data = bees.db.100)

summary(area.vol)
mod_dharma1 <- area.vol %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.vol) 
AICc(area.vol)
plot(allEffects(area.vol))

res.100 <- simulateResiduals(area.vol)
res.agg <- recalculateResiduals(res.100,  group = bees.db.100$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.db.100$X, bees.db.100$Site, mean),
  y = tapply(bees.db.100$Y, bees.db.100$Site, mean))

# 9. Test interaction habitat area * forest
area.forest <- glmmTMB(total_abundance ~ Category + scale(Hab_Area) * scale(Forest) +
                         scale(Vegetated.open.land) + scale(honeybee) + scale(Urban) +
                         (1 | City/Triplet),
                       family = nbinom2,
                       data = bees.db.100)

summary(area.forest)
mod_dharma1 <- area.forest %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.forest) 
AICc(area.forest)
plot(allEffects(area.forest))

res.100 <- simulateResiduals(area.forest)
res.agg <- recalculateResiduals(res.100,  group = bees.db.100$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.db.100$X, bees.db.100$Site, mean),
  y = tapply(bees.db.100$Y, bees.db.100$Site, mean))

# Collect AICc values

mods <- list(
  base = base.bees.100,
  cat.honey = cat.honey, 
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
# structure data
honeybee_df <- abund.1.raw %>%
  group_by(Site) %>%
  summarise(honeybee = sum(honeybee_abundance, na.rm = TRUE), .groups = "drop")


abund.1.bees <- abund.1.raw %>%
  filter(Pollinator.type == "Bee") |> 
  group_by(City, Site, Triplet, Category, X, Y, Hab_Area, 
           Arable, Forest, Non.vegetated.open.land, Urban, Vegetated.open.land, Water, Wetland) 
abund.1.bees

bees.db.1 <- left_join(honeybee_df, abund.1.bees, by = c("Site"))
bees.db.1 <- bees.db.1 |> 
  mutate(City = as.factor(City)) |> 
  mutate(Site = as.factor(Site)) |> 
  mutate(Triplet = as.factor(Triplet)) |> 
  mutate(Category = as.factor(Category)) 

# 1. base model
base.bees.1 <- glmmTMB(total_abundance ~ Category + scale(honeybee) + scale(Hab_Area) + 
                           scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                           (1 | City/Triplet),
                         family = nbinom2,
                         data = bees.db.1)

summary(base.bees.1)
mod_dharma1 <- base.bees.1 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(base.bees.1) 
AICc(base.bees.1)
plot(allEffects(base.bees.1))

res.1 <- simulateResiduals(base.bees.1)
res.agg <- recalculateResiduals(res.1,  group = bees.db.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.db.1$X, bees.db.1$Site, mean),
  y = tapply(bees.db.1$Y, bees.db.1$Site, mean))


# 2. Test interaction category * honeybee abundance
cat.honey <- glmmTMB(total_abundance ~ Category * scale(honeybee) + scale(Hab_Area) + 
                       scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                       (1 | City/Triplet),
                     family = nbinom2,
                     data = bees.db.1)

summary(cat.honey)
mod_dharma1 <- cat.honey %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.honey) 
AICc(cat.honey)
plot(allEffects(cat.honey))

res.1 <- simulateResiduals(cat.honey)
res.agg <- recalculateResiduals(res.1,  group = bees.db.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.db.1$X, bees.db.1$Site, mean),
  y = tapply(bees.db.1$Y, bees.db.1$Site, mean))

# 3. Test interaction category * habitat area
cat.area <- glmmTMB(total_abundance ~ Category * scale(Hab_Area) + scale(honeybee) +   
                      scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                      (1 | City/Triplet),
                    family = nbinom2,
                    data = bees.db.1)

summary(cat.area)
mod_dharma1 <- cat.area %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.area) 
AICc(cat.area)
plot(allEffects(cat.area))

res.1 <- simulateResiduals(cat.area)
res.agg <- recalculateResiduals(res.1,  group = bees.db.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.db.1$X, bees.db.1$Site, mean),
  y = tapply(bees.db.1$Y, bees.db.1$Site, mean))

# 4. Test category  * urban
cat.urban <- glmmTMB(total_abundance ~ Category * scale(Urban) +
                       scale(Hab_Area) + scale(honeybee) +   
                       scale(Vegetated.open.land) + scale(Forest) + 
                       (1 | City/Triplet),
                     family = nbinom2,
                     data = bees.db.1)

summary(cat.urban)
mod_dharma1 <- cat.urban %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.urban) 
AICc(cat.urban)
plot(allEffects(cat.urban))

res.1 <- simulateResiduals(cat.urban)
res.agg <- recalculateResiduals(res.1,  group = bees.db.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.db.1$X, bees.db.1$Site, mean),
  y = tapply(bees.db.1$Y, bees.db.1$Site, mean))

# 5. Test interaction category * vegetated open land
cat.vol <- glmmTMB(total_abundance ~ Category * scale(Vegetated.open.land) +
                     scale(Urban) + scale(Hab_Area) + scale(honeybee) + scale(Forest) + 
                     (1 | City/Triplet),
                   family = nbinom2,
                   data = bees.db.1)

summary(cat.vol)
mod_dharma1 <- cat.vol %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.vol) 
AICc(cat.vol)
plot(allEffects(cat.vol))

res.1 <- simulateResiduals(cat.vol)
res.agg <- recalculateResiduals(res.1,  group = bees.db.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.db.1$X, bees.db.1$Site, mean),
  y = tapply(bees.db.1$Y, bees.db.1$Site, mean))

# 6. Test interaction category * forest
cat.forest <- glmmTMB(total_abundance ~ Category * scale(Forest) +
                        scale(Vegetated.open.land) + scale(Urban) + scale(Hab_Area) + scale(honeybee) + 
                        (1 | City/Triplet),
                      family = nbinom2,
                      data = bees.db.1)

summary(cat.forest)
mod_dharma1 <- cat.forest %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.forest) 
AICc(cat.forest)
plot(allEffects(cat.forest))

res.1 <- simulateResiduals(cat.forest)
res.agg <- recalculateResiduals(res.1,  group = bees.db.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.db.1$X, bees.db.1$Site, mean),
  y = tapply(bees.db.1$Y, bees.db.1$Site, mean))

# 7. Test interaction habitat area * urban
area.urban <- glmmTMB(total_abundance ~ Category + scale(Hab_Area) * scale(Urban) +
                        scale(Forest) + scale(Vegetated.open.land) +  scale(honeybee) + 
                        (1 | City/Triplet),
                      family = nbinom2,
                      data = bees.db.1)

summary(area.urban)
mod_dharma1 <- area.urban %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.urban) 
AICc(area.urban)
plot(allEffects(area.urban))

res.1 <- simulateResiduals(area.urban)
res.agg <- recalculateResiduals(res.1,  group = bees.db.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.db.1$X, bees.db.1$Site, mean),
  y = tapply(bees.db.1$Y, bees.db.1$Site, mean))


# 8.Test interaction habitat area * vegetated open land
area.vol <- glmmTMB(total_abundance ~ Category + scale(Hab_Area) * scale(Vegetated.open.land) +
                      scale(Forest) + scale(honeybee) + scale(Urban) +
                      (1 | City/Triplet),
                    family = nbinom2,
                    data = bees.db.1)

summary(area.vol)
mod_dharma1 <- area.vol %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.vol) 
AICc(area.vol)
plot(allEffects(area.vol))

res.1 <- simulateResiduals(area.vol)
res.agg <- recalculateResiduals(res.1,  group = bees.db.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.db.1$X, bees.db.1$Site, mean),
  y = tapply(bees.db.1$Y, bees.db.1$Site, mean))

# 9. Test interaction habitat area * forest
area.forest <- glmmTMB(total_abundance ~ Category + scale(Hab_Area) * scale(Forest) +
                         scale(Vegetated.open.land) + scale(honeybee) + scale(Urban) +
                         (1 | City/Triplet),
                       family = nbinom2,
                       data = bees.db.1)

summary(area.forest)
mod_dharma1 <- area.forest %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.forest) 
AICc(area.forest)
plot(allEffects(area.forest))

res.1 <- simulateResiduals(area.forest)
res.agg <- recalculateResiduals(res.1,  group = bees.db.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(bees.db.1$X, bees.db.1$Site, mean),
  y = tapply(bees.db.1$Y, bees.db.1$Site, mean))

# Collect AICc values

mods <- list(
  base = base.bees.1,
  #cat.honey = cat.honey, # did not converge
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
