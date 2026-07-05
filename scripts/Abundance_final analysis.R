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

# Analysis without the pollinator groups has been removed from the analysis, but the code remains 
#in "scripts/Abundance_analysis_honeybees.R"

# Pollinator groups
# structure data
# 100 m:-----
honeybee_df <- abund.100.raw %>%
  group_by(Site) %>%
  summarise(honeybee = sum(honeybee_abundance, na.rm = TRUE), .groups = "drop")


abund.100.groups <- abund.100.raw %>%
  mutate(Pollinator.type = ifelse(Pollinator.type == "Burnet moth", "Butterfly",Pollinator.type)) %>%
  filter(Pollinator.type %in% c("Bee", "Bumblebee", "Butterfly", "Hoverfly")) |> 
  group_by(City, Site, Triplet, Category, Pollinator.type, X, Y, Hab_Area, 
           Arable, Forest, Non.vegetated.open.land, Urban, Vegetated.open.land, Water, Wetland) |> 
  summarise(wild.pollinator.abund = sum(total_abundance), .groups = "drop_last")
abund.100.groups

group.db.100 <- left_join(honeybee_df, abund.100.groups, by = c("Site"))
group.db.100<- group.db.100 |> 
  mutate(City = as.factor(City)) |> 
  mutate(Site = as.factor(Site)) |> 
  mutate(Triplet = as.factor(Triplet)) |> 
  mutate(Category = as.factor(Category)) |>
  mutate(Pollinator.type = as.factor(Pollinator.type))


# Landscape analysis - pollinator groups 100 m: ----
# 1. base model
base.group.100 <- glmmTMB(wild.pollinator.abund ~ Category + Pollinator.type + scale(honeybee) + scale(Hab_Area) + 
                            scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                            (1 | City/Triplet),
                          family = nbinom2,
                          data = group.db.100)

summary(base.group.100)
mod_dharma1 <- base.group.100 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(base.group.100) 
AICc(base.group.100)
plot(allEffects(base.group.100))

res.100 <- simulateResiduals(base.group.100)
res.agg <- recalculateResiduals(res.100,  group = group.db.100$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(group.db.100$X, group.db.100$Site, mean),
  y = tapply(group.db.100$Y, group.db.100$Site, mean))


# 2. Test interaction Category * pollinator group
cat.pollinator <- glmmTMB(wild.pollinator.abund ~ Category * Pollinator.type + scale(honeybee) + scale(Hab_Area) + 
                        scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                        (1 | City/Triplet),
                      family = nbinom2,
                      data = group.db.100)

summary(cat.pollinator)
mod_dharma1 <- cat.pollinator %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.pollinator) 
AICc(cat.pollinator)
plot(allEffects(cat.pollinator))

res.100 <- simulateResiduals(cat.pollinator)
res.agg <- recalculateResiduals(res.100,  group = group.db.100$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(group.db.100$X, group.db.100$Site, mean),
  y = tapply(group.db.100$Y, group.db.100$Site, mean))

# 3. Test interaction Pollinator group * honeybee abundance
pollinator.honey <- glmmTMB(wild.pollinator.abund ~ Category + Pollinator.type * scale(honeybee) + scale(Hab_Area) + 
                              scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                              (1 | City/Triplet),
                            family = nbinom2,
                            data = group.db.100)

summary(pollinator.honey)
mod_dharma1 <- pollinator.honey %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(pollinator.honey) 
AICc(pollinator.honey)
plot(allEffects(pollinator.honey))

res.100 <- simulateResiduals(pollinator.honey)
res.agg <- recalculateResiduals(res.100,  group = group.db.100$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(group.db.100$X, group.db.100$Site, mean),
  y = tapply(group.db.100$Y, group.db.100$Site, mean))

# 4. Test interaction Pollinator group * habitat area
pollinator.area <- glmmTMB(wild.pollinator.abund ~ Category + Pollinator.type * scale(Hab_Area) +scale(honeybee) +   
                              scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                              (1 | City/Triplet),
                            family = nbinom2,
                            data = group.db.100)

summary(pollinator.area)
mod_dharma1 <- pollinator.area %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(pollinator.area) 
AICc(pollinator.area)
plot(allEffects(pollinator.area))

res.100 <- simulateResiduals(pollinator.area)
res.agg <- recalculateResiduals(res.100,  group = group.db.100$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(group.db.100$X, group.db.100$Site, mean),
  y = tapply(group.db.100$Y, group.db.100$Site, mean))

# 5. Test interaction Pollinator group * urban
pollinator.urban <- glmmTMB(wild.pollinator.abund ~ Category + Pollinator.type * scale(Urban) +
                             scale(Hab_Area) + scale(honeybee) +   
                             scale(Vegetated.open.land) + scale(Forest) + 
                             (1 | City/Triplet),
                           family = nbinom2,
                           data = group.db.100)

summary(pollinator.urban)
mod_dharma1 <- pollinator.urban %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(pollinator.urban) 
AICc(pollinator.urban)
plot(allEffects(pollinator.urban))

res.100 <- simulateResiduals(pollinator.urban)
res.agg <- recalculateResiduals(res.100,  group = group.db.100$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(group.db.100$X, group.db.100$Site, mean),
  y = tapply(group.db.100$Y, group.db.100$Site, mean))

# 6. Test interaction Pollinator group * vegetated open land
pollinator.vol <- glmmTMB(wild.pollinator.abund ~ Category + Pollinator.type * scale(Vegetated.open.land) +
                              scale(Urban) + scale(Hab_Area) + scale(honeybee) + scale(Forest) + 
                              (1 | City/Triplet),
                            family = nbinom2,
                            data = group.db.100)

summary(pollinator.vol)
mod_dharma1 <- pollinator.vol %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(pollinator.vol) 
AICc(pollinator.vol)
plot(allEffects(pollinator.vol))

res.100 <- simulateResiduals(pollinator.vol)
res.agg <- recalculateResiduals(res.100,  group = group.db.100$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(group.db.100$X, group.db.100$Site, mean),
  y = tapply(group.db.100$Y, group.db.100$Site, mean))

# 7. Test interaction Pollinator group * forest
pollinator.forest <- glmmTMB(wild.pollinator.abund ~ Category + Pollinator.type * scale(Forest) +
                            scale(Vegetated.open.land) + scale(Urban) + scale(Hab_Area) + scale(honeybee) + 
                            (1 | City/Triplet),
                          family = nbinom2,
                          data = group.db.100)

summary(pollinator.forest)
mod_dharma1 <- pollinator.forest %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(pollinator.forest) 
AICc(pollinator.forest)
plot(allEffects(pollinator.forest))

res.100 <- simulateResiduals(pollinator.forest)
res.agg <- recalculateResiduals(res.100,  group = group.db.100$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(group.db.100$X, group.db.100$Site, mean),
  y = tapply(group.db.100$Y, group.db.100$Site, mean))

# 8. Test interaction habitat area * urban
area.urban <- glmmTMB(wild.pollinator.abund ~ Category + Pollinator.type + scale(Hab_Area) * scale(Urban) +
                        scale(Forest) + scale(Vegetated.open.land) +  scale(honeybee) + 
                               (1 | City/Triplet),
                             family = nbinom2,
                             data = group.db.100)

summary(area.urban)
mod_dharma1 <- area.urban %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.urban) 
AICc(area.urban)
plot(allEffects(area.urban))

res.100 <- simulateResiduals(area.urban)
res.agg <- recalculateResiduals(res.100,  group = group.db.100$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(group.db.100$X, group.db.100$Site, mean),
  y = tapply(group.db.100$Y, group.db.100$Site, mean))


# 9.Test interaction habitat area * vegetated open land
area.vol <- glmmTMB(wild.pollinator.abund ~ Category + Pollinator.type + scale(Hab_Area) * scale(Vegetated.open.land) +
                        scale(Forest) + scale(honeybee) + scale(Urban) +
                        (1 | City/Triplet),
                      family = nbinom2,
                      data = group.db.100)

summary(area.vol)
mod_dharma1 <- area.vol %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.vol) 
AICc(area.vol)
plot(allEffects(area.vol))

res.100 <- simulateResiduals(area.vol)
res.agg <- recalculateResiduals(res.100,  group = group.db.100$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(group.db.100$X, group.db.100$Site, mean),
  y = tapply(group.db.100$Y, group.db.100$Site, mean))

# 10. Test interaction habitat area * forest
area.forest <- glmmTMB(wild.pollinator.abund ~ Category + Pollinator.type + scale(Hab_Area) * scale(Forest) +
                      scale(Vegetated.open.land) + scale(honeybee) + scale(Urban) +
                      (1 | City/Triplet),
                    family = nbinom2,
                    data = group.db.100)

summary(area.forest)
mod_dharma1 <- area.forest %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.forest) 
AICc(area.forest)
plot(allEffects(area.forest))

res.100 <- simulateResiduals(area.forest)
res.agg <- recalculateResiduals(res.100,  group = group.db.100$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(group.db.100$X, group.db.100$Site, mean),
  y = tapply(group.db.100$Y, group.db.100$Site, mean))

# Collect AICc values

mods <- list(
  base = base.group.100,
  cat.group = cat.pollinator, 
  group.honey = pollinator.honey, 
  area.group = pollinator.area, 
  urban = pollinator.urban,
  open = pollinator.vol,
  forest = pollinator.forest,
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
honeybee_df <- abund.1.raw %>%
  group_by(Site) %>%
  summarise(honeybee = sum(honeybee_abundance, na.rm = TRUE), .groups = "drop")


abund.1.groups <- abund.1.raw %>%
  mutate(Pollinator.type = ifelse(Pollinator.type == "Burnet moth", "Butterfly",Pollinator.type)) %>%
  filter(Pollinator.type %in% c("Bee", "Bumblebee", "Butterfly", "Hoverfly")) |> 
  group_by(City, Site, Triplet, Category, Pollinator.type, X, Y, Hab_Area, 
           Arable, Forest, Non.vegetated.open.land, Urban, Vegetated.open.land, Water, Wetland) |> 
  summarise(wild.pollinator.abund = sum(total_abundance), .groups = "drop_last")
abund.1.groups

group.db.1 <- left_join(honeybee_df, abund.1.groups, by = c("Site"))

# Landscape analysis - pollinator groups 100 m: ----
# 1. base model
base.group.1 <- glmmTMB(wild.pollinator.abund ~ Category + Pollinator.type + scale(honeybee) + scale(Hab_Area) + 
                            scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                            (1 | City/Triplet),
                          family = nbinom2,
                          data = group.db.1)

summary(base.group.1)
mod_dharma1 <- base.group.1 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(base.group.1) 
AICc(base.group.1)
plot(allEffects(base.group.1))

res.1 <- simulateResiduals(base.group.1)
res.agg <- recalculateResiduals(res.1,  group = group.db.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(group.db.1$X, group.db.1$Site, mean),
  y = tapply(group.db.1$Y, group.db.1$Site, mean))


# 2. Test interaction Category * pollinator group
cat.pollinator <- glmmTMB(wild.pollinator.abund ~ Category * Pollinator.type + scale(honeybee) + scale(Hab_Area) + 
                            scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                            (1 | City/Triplet),
                          family = nbinom2,
                          data = group.db.1)

summary(cat.pollinator)
mod_dharma1 <- cat.pollinator %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.pollinator) 
AICc(cat.pollinator)
plot(allEffects(cat.pollinator))

res.1 <- simulateResiduals(cat.pollinator)
res.agg <- recalculateResiduals(res.1,  group = group.db.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(group.db.1$X, group.db.1$Site, mean),
  y = tapply(group.db.1$Y, group.db.1$Site, mean))

# 3. Test interaction Pollinator group * honeybee abundance
pollinator.honey <- glmmTMB(wild.pollinator.abund ~ Category + Pollinator.type * scale(honeybee) + scale(Hab_Area) + 
                              scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                              (1 | City/Triplet),
                            family = nbinom2,
                            data = group.db.1)

summary(pollinator.honey)
mod_dharma1 <- pollinator.honey %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(pollinator.honey) 
AICc(pollinator.honey)
plot(allEffects(pollinator.honey))

res.1 <- simulateResiduals(pollinator.honey)
res.agg <- recalculateResiduals(res.1,  group = group.db.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(group.db.1$X, group.db.1$Site, mean),
  y = tapply(group.db.1$Y, group.db.1$Site, mean))

# 4. Test interaction Pollinator group * habitat area
pollinator.area <- glmmTMB(wild.pollinator.abund ~ Category + Pollinator.type * scale(Hab_Area) +scale(honeybee) +   
                             scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                             (1 | City/Triplet),
                           family = nbinom2,
                           data = group.db.1)

summary(pollinator.area)
mod_dharma1 <- pollinator.area %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(pollinator.area) 
AICc(pollinator.area)
plot(allEffects(pollinator.area))

res.1 <- simulateResiduals(pollinator.area)
res.agg <- recalculateResiduals(res.1,  group = group.db.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(group.db.1$X, group.db.1$Site, mean),
  y = tapply(group.db.1$Y, group.db.1$Site, mean))

# 5. Test interaction Pollinator group * urban
pollinator.urban <- glmmTMB(wild.pollinator.abund ~ Category + Pollinator.type * scale(Urban) +
                              scale(Hab_Area) + scale(honeybee) +   
                              scale(Vegetated.open.land) + scale(Forest) + 
                              (1 | City/Triplet),
                            family = nbinom2,
                            data = group.db.1)

summary(pollinator.urban)
mod_dharma1 <- pollinator.urban %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(pollinator.urban) 
AICc(pollinator.urban)
plot(allEffects(pollinator.urban))

res.1 <- simulateResiduals(pollinator.urban)
res.agg <- recalculateResiduals(res.1,  group = group.db.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(group.db.1$X, group.db.1$Site, mean),
  y = tapply(group.db.1$Y, group.db.1$Site, mean))

# 6. Test interaction Pollinator group * vegetated open land
pollinator.vol <- glmmTMB(wild.pollinator.abund ~ Category + Pollinator.type * scale(Vegetated.open.land) +
                            scale(Urban) + scale(Hab_Area) + scale(honeybee) + scale(Forest) + 
                            (1 | City/Triplet),
                          family = nbinom2,
                          data = group.db.1)

summary(pollinator.vol)
mod_dharma1 <- pollinator.vol %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(pollinator.vol) 
AICc(pollinator.vol)
plot(allEffects(pollinator.vol))

res.1 <- simulateResiduals(pollinator.vol)
res.agg <- recalculateResiduals(res.1,  group = group.db.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(group.db.1$X, group.db.1$Site, mean),
  y = tapply(group.db.1$Y, group.db.1$Site, mean))

# 7. Test interaction Pollinator group * forest
pollinator.forest <- glmmTMB(wild.pollinator.abund ~ Category + Pollinator.type * scale(Forest) +
                               scale(Vegetated.open.land) + scale(Urban) + scale(Hab_Area) + scale(honeybee) + 
                               (1 | City/Triplet),
                             family = nbinom2,
                             data = group.db.1)

summary(pollinator.forest)
mod_dharma1 <- pollinator.forest %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(pollinator.forest) 
AICc(pollinator.forest)
plot(allEffects(pollinator.forest))

res.1 <- simulateResiduals(pollinator.forest)
res.agg <- recalculateResiduals(res.1,  group = group.db.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(group.db.1$X, group.db.1$Site, mean),
  y = tapply(group.db.1$Y, group.db.1$Site, mean))

# 8. Test interaction habitat area * urban
area.urban <- glmmTMB(wild.pollinator.abund ~ Category + Pollinator.type + scale(Hab_Area) * scale(Urban) +
                        scale(Forest) + scale(Vegetated.open.land) +  scale(honeybee) + 
                        (1 | City/Triplet),
                      family = nbinom2,
                      data = group.db.1)

summary(area.urban)
mod_dharma1 <- area.urban %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.urban) 
AICc(area.urban)
plot(allEffects(area.urban))

res.1 <- simulateResiduals(area.urban)
res.agg <- recalculateResiduals(res.1,  group = group.db.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(group.db.1$X, group.db.1$Site, mean),
  y = tapply(group.db.1$Y, group.db.1$Site, mean))


# 9.Test interaction habitat area * vegetated open land
area.vol <- glmmTMB(wild.pollinator.abund ~ Category + Pollinator.type + scale(Hab_Area) * scale(Vegetated.open.land) +
                      scale(Forest) + scale(honeybee) + scale(Urban) +
                      (1 | City/Triplet),
                    family = nbinom2,
                    data = group.db.1)

summary(area.vol)
mod_dharma1 <- area.vol %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.vol) 
AICc(area.vol)
plot(allEffects(area.vol))

res.1 <- simulateResiduals(area.vol)
res.agg <- recalculateResiduals(res.1,  group = group.db.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(group.db.1$X, group.db.1$Site, mean),
  y = tapply(group.db.1$Y, group.db.1$Site, mean))

# 10. Test interaction habitat area * forest
area.forest <- glmmTMB(wild.pollinator.abund ~ Category + Pollinator.type + scale(Hab_Area) * scale(Forest) +
                         scale(Vegetated.open.land) + scale(honeybee) + scale(Urban) +
                         (1 | City/Triplet),
                       family = nbinom2,
                       data = group.db.1)

summary(area.forest)
mod_dharma1 <- area.forest %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.forest) 
AICc(area.forest)
plot(allEffects(area.forest))

res.1 <- simulateResiduals(area.forest)
res.agg <- recalculateResiduals(res.1,  group = group.db.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(group.db.1$X, group.db.1$Site, mean),
  y = tapply(group.db.1$Y, group.db.1$Site, mean))

# Collect AICc values

mods <- list(
  base = base.group.1,
  cat.group = cat.pollinator, 
  group.honey = pollinator.honey, 
  area.group = pollinator.area, 
  urban = pollinator.urban,
  group.open = pollinator.vol,
  forest = pollinator.forest,
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





