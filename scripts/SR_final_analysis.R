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
library(purrr)


# Load data
sr.100.raw <- read.csv2(here("data/databases", "tot.richness.100m.csv" ), fileEncoding = "latin1")
sr.1.raw <- read.csv2(here("data/databases", "tot.richness.1km.csv" ), fileEncoding = "latin1")

# structure data
sr.100 <- sr.100.raw %>% 
  mutate(City = as.factor(City)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Triplet = as.factor(Triplet)) %>% 
  mutate(Category = as.factor(Category))

sr.1 <- sr.1.raw %>% 
  mutate(City = as.factor(City)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Triplet = as.factor(Triplet)) %>% 
  mutate(Category = as.factor(Category))

# Compare models with 2 random effects (City, and City/Triplet)

# 100 m: 
mod.city.100 <- glmmTMB(richness ~ Category + scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) +
    (1 | City),
  family = poisson, data = sr.100)

summary(mod.city.100)
mod_dharma1 <- mod.city.100 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.city.100) 
AICc(mod.city.100) 

mod.nested.100 <- glmmTMB(richness ~ Category + scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) +
    (1 | City/Triplet),
  family = poisson,
  data = sr.100)

summary(mod.nested.100)
mod_dharma1 <- mod.nested.100 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.nested.100) 
AICc(mod.nested.100) 

anova(mod.city.100, mod.nested.100)
AICc(mod.city.100, mod.nested.100)

res.100 <- simulateResiduals(mod.nested.100)

testSpatialAutocorrelation(
  res.100,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))

# 1 km:
mod.city <- glmmTMB(richness ~ Category + scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) +
                      (1 | City),
                    family = poisson, data = sr.1)

summary(mod.city)
mod_dharma1 <- mod.city %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.city) 
AICc(mod.city) 

mod.nested.1 <- glmmTMB(richness ~ Category + scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) +
                        (1 | City/Triplet),
                      family = poisson,
                      data = sr.1)

summary(mod.nested)
mod_dharma1 <- mod.nested %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.nested) 
AICc(mod.nested) 

anova(mod.city, mod.nested)
AICc(mod.city, mod.nested)

res.1000 <- simulateResiduals(mod.nested)

testSpatialAutocorrelation(
  res.1000,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))




# Final analysis ----

# 100 m: test intractions one by one ----
# base model
mod.base <- glmmTMB(richness ~ Category + scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) +
                          (1 | City/Triplet),
                        family = poisson,
                        data = sr.100)

summary(mod.base)
mod_dharma1 <- mod.base %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.base) 
AICc(mod.base) 
plot(allEffects(mod.base))
visreg(mod.base, scale = "response")

res.100 <- simulateResiduals(mod.base)

testSpatialAutocorrelation(
  res.100,
  x = tapply(sr.100$X, sr.100$Site, mean),
  y = tapply(sr.100$Y, sr.100$Site, mean))


# Category*area
mod.area.cat <- glmmTMB(richness ~ Category * scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) +
                        (1 | City/Triplet),
                      family = poisson,
                      data = sr.100)

summary(mod.area.cat)
mod_dharma1 <- mod.area.cat %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.area.cat) 
AICc(mod.area.cat) 
plot(allEffects(mod.area.cat))
visreg(mod.area.cat, scale = "response")

res.100 <- simulateResiduals(mod.area.cat)

testSpatialAutocorrelation(
  res.100,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))

anova(mod.nested.100, mod.area.cat)

# Urban
mod.urban <- glmmTMB(richness ~ Category * scale(Urban) + scale(Hab_Area)  + scale(Vegetated.open.land) + scale(Forest) +
                          (1 | City/Triplet),
                        family = poisson,
                        data = sr.100)

summary(mod.urban)
mod_dharma1 <- mod.urban %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.urban) 
AICc(mod.urban) 
plot(allEffects(mod.urban))
visreg(mod.urban, scale = "response")

res.100 <- simulateResiduals(mod.urban)

testSpatialAutocorrelation(
  res.100,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))

anova(mod.nested.100, mod.urban)

# Vegetated open land
mod.vol <- glmmTMB(richness ~ Category * scale(Vegetated.open.land) + scale(Hab_Area)  + scale(Urban) + scale(Forest) +
                       (1 | City/Triplet),
                     family = poisson,
                     data = sr.100)

summary(mod.vol)
mod_dharma1 <- mod.vol %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.vol) 
AICc(mod.vol) 
plot(allEffects(mod.vol))
visreg(mod.vol, scale = "response")

res.100 <- simulateResiduals(mod.vol)

testSpatialAutocorrelation(
  res.100,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))

# Forest
mod.forest <- glmmTMB(richness ~ Category * scale(Forest) + scale(Hab_Area)  + scale(Urban) + scale(Vegetated.open.land) +
                     (1 | City/Triplet),
                   family = poisson,
                   data = sr.100)

summary(mod.forest)
mod_dharma1 <- mod.forest %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.forest) 
AICc(mod.forest) 
plot(allEffects(mod.forest))
visreg(mod.forest, scale = "response")

res.100 <- simulateResiduals(mod.forest)

testSpatialAutocorrelation(
  res.100,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))

# Habitat area * Urban cover
mod.area.urban <- glmmTMB(richness ~ Category +  scale(Hab_Area) * scale(Urban) + scale(Vegetated.open.land) + scale(Forest) +
                        (1 | City/Triplet),
                      family = poisson,
                      data = sr.100)

summary(mod.area.urban)
mod_dharma1 <- mod.area.urban %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.area.urban) 
AICc(mod.area.urban) 
plot(allEffects(mod.area.urban))
visreg(mod.area.urban, scale = "response")

res.100 <- simulateResiduals(mod.area.urban)

testSpatialAutocorrelation(
  res.100,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))

# Habitat area * Vegetated open land
mod.area.vol <- glmmTMB(richness ~ Category +  scale(Hab_Area) * scale(Vegetated.open.land) + scale(Urban) + scale(Forest) +
                            (1 | City/Triplet),
                          family = poisson,
                          data = sr.100)

summary(mod.area.vol)
mod_dharma1 <- mod.area.vol %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.area.vol) 
AICc(mod.area.vol) 
plot(allEffects(mod.area.vol))
visreg(mod.area.vol, scale = "response")

res.100 <- simulateResiduals(mod.area.vol)

testSpatialAutocorrelation(
  res.100,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))

# Habitat area * Forest
mod.area.forest <- glmmTMB(richness ~ Category +  scale(Hab_Area) * scale(Forest) + scale(Vegetated.open.land) + scale(Urban) +
                          (1 | City/Triplet),
                        family = poisson,
                        data = sr.100)

summary(mod.area.forest)
mod_dharma1 <- mod.area.forest %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.area.forest) 
AICc(mod.area.forest) 
plot(allEffects(mod.area.forest))
visreg(mod.area.forest, scale = "response")

res.100 <- simulateResiduals(mod.area.forest)

testSpatialAutocorrelation(
  res.100,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))

# Collect AICc values

mods <- list(
  base = mod.base,
  urban = mod.urban,
  open = mod.vol,
  forest = mod.forest,
  area = mod.area.cat, 
  area.urban = mod.area.urban, 
  area.vol = mod.area.vol, 
  area.forest = mod.area.forest
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



# 1 km: test intractions one by one ----
# base
mod.base.1 <- glmmTMB(richness ~ Category + scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) +
                          (1 | City/Triplet),
                        family = poisson,
                        data = sr.1)

summary(mod.base.1)
mod_dharma1 <- mod.base.1 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.base.1) 
AICc(mod.base.1) 
plot(allEffects(mod.base.1))
visreg(mod.base.1, scale = "response")

res.1 <- simulateResiduals(mod.base.1)

testSpatialAutocorrelation(
  res.1,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))

anova(mod.nested.1, mod.area.cat)

# Category*area
mod.area.cat <- glmmTMB(richness ~ Category * scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) +
                          (1 | City/Triplet),
                        family = poisson,
                        data = sr.1)

summary(mod.area.cat)
mod_dharma1 <- mod.area.cat %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.area.cat) 
AICc(mod.area.cat) 
plot(allEffects(mod.area.cat))
visreg(mod.area.cat, scale = "response")

res.1 <- simulateResiduals(mod.area.cat)

testSpatialAutocorrelation(
  res.1,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))

anova(mod.nested.1, mod.area.cat)

# Urban
mod.urban <- glmmTMB(richness ~ Category * scale(Urban) + scale(Hab_Area)  + scale(Vegetated.open.land) + scale(Forest) +
                       (1 | City/Triplet),
                     family = poisson,
                     data = sr.1)

summary(mod.urban)
mod_dharma1 <- mod.urban %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.urban) 
AICc(mod.urban) 
plot(allEffects(mod.urban))
visreg(mod.urban, scale = "response")

res.1 <- simulateResiduals(mod.urban)

testSpatialAutocorrelation(
  res.1,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))

anova(mod.nested.1, mod.urban)

# Vegetated open land
mod.vol <- glmmTMB(richness ~ Category * scale(Vegetated.open.land) + scale(Hab_Area)  + scale(Urban) + scale(Forest) +
                     (1 | City/Triplet),
                   family = poisson,
                   data = sr.1)

summary(mod.vol)
mod_dharma1 <- mod.vol %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.vol) 
AICc(mod.vol) 
plot(allEffects(mod.vol))
visreg(mod.vol, scale = "response")

res.1 <- simulateResiduals(mod.vol)

testSpatialAutocorrelation(
  res.1,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))

# Forest
mod.forest <- glmmTMB(richness ~ Category * scale(Forest) + scale(Hab_Area)  + scale(Urban) + scale(Vegetated.open.land) +
                        (1 | City/Triplet),
                      family = poisson,
                      data = sr.1)

summary(mod.forest)
mod_dharma1 <- mod.forest %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.forest) 
AICc(mod.forest) 
plot(allEffects(mod.forest))
visreg(mod.forest, scale = "response")

res.1 <- simulateResiduals(mod.forest)

testSpatialAutocorrelation(
  res.1,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))

# Habitat area * Urban cover
mod.area.urban <- glmmTMB(richness ~ Category +  scale(Hab_Area) * scale(Urban) + scale(Vegetated.open.land) + scale(Forest) +
                            (1 | City/Triplet),
                          family = poisson,
                          data = sr.1)

summary(mod.area.urban)
mod_dharma1 <- mod.area.urban %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.area.urban) 
AICc(mod.area.urban) 
plot(allEffects(mod.area.urban))
visreg(mod.area.urban, scale = "response")

res.1 <- simulateResiduals(mod.area.urban)

testSpatialAutocorrelation(
  res.1,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))

# Habitat area * Vegetated open land
mod.area.vol <- glmmTMB(richness ~ Category +  scale(Hab_Area) * scale(Vegetated.open.land) + scale(Urban) + scale(Forest) +
                          (1 | City/Triplet),
                        family = poisson,
                        data = sr.1)

summary(mod.area.vol)
mod_dharma1 <- mod.area.vol %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.area.vol) 
AICc(mod.area.vol) 
plot(allEffects(mod.area.vol))
visreg(mod.area.vol, scale = "response")

res.1 <- simulateResiduals(mod.area.vol)

testSpatialAutocorrelation(
  res.100,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))

# Habitat area * Forest
mod.area.forest <- glmmTMB(richness ~ Category +  scale(Hab_Area) * scale(Forest) + scale(Vegetated.open.land) + scale(Urban) +
                             (1 | City/Triplet),
                           family = poisson,
                           data = sr.1)

summary(mod.area.forest)
mod_dharma1 <- mod.area.forest %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.area.forest) 
AICc(mod.area.forest) 
plot(allEffects(mod.area.forest))
visreg(mod.area.forest, scale = "response")

res.1 <- simulateResiduals(mod.area.forest)

testSpatialAutocorrelation(
  res.1,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))

# Collect AICc values

mods <- list(
  base = mod.base.1,
  urban = mod.urban,
  open = mod.vol,
  forest = mod.forest,
  area = mod.area.cat, 
  area.urban = mod.area.urban, 
  area.vol = mod.area.vol, 
  area.forest = mod.area.forest
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
  mutate(delta_AICc = AICc - min(AIC))

aic_table

