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

# sum all pollinators
# 100 m:
abund.100.total <- abund.100.raw %>%
  mutate(City = as.factor(City)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Triplet = as.factor(Triplet)) %>% 
  mutate(Category = as.factor(Category)) %>% 
  group_by(City, Site, Triplet, Category, X, Y, Hab_Area, 
           Arable, Forest, Non.vegetated.open.land, Urban, Vegetated.open.land, Water, Wetland) %>%
  summarise(total.wild.abund = sum(wild_pollinator_abundance, na.rm = TRUE),
    honeybee.abund = sum(honeybee_abundance, na.rm = TRUE),
    .groups = "drop")
abund.100.total

# 1 km:
abund.1.total <- abund.1.raw %>%
  group_by(City, Site, Triplet, Category, X, Y, Hab_Area,
    Arable, Forest, Non.vegetated.open.land, Urban, Vegetated.open.land, Water, Wetland) %>%
  summarise(total.wild.abund = sum(wild_pollinator_abundance, na.rm = TRUE),
    honeybee.abund = sum(honeybee_abundance, na.rm = TRUE),
    .groups = "drop")
abund.1.total


# Landscape analysis - 100 m:  ----
# Base model:
mod.nested.100 <- glmmTMB(total.wild.abund ~ Category + scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + 
                            scale(Forest) + scale(honeybee.abund) +
                            (1 | City/Triplet),
                          family = nbinom2,
                          data = abund.100.total)

summary(mod.nested.100)
mod_dharma1 <- mod.nested.100 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.nested.100) 
AICc(mod.nested.100) 

res.100 <- simulateResiduals(mod.nested.100)

testSpatialAutocorrelation(
  res.100,
  x = tapply(abund.100.total$X, abund.100.total$Site, mean),
  y = tapply(abund.100.total$Y, abund.100.total$Site, mean))


# 100 m: test intractions one by one ----
# Category*area
mod.area.cat <- glmmTMB(total.wild.abund ~ Category * scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + 
                          scale(Forest) + scale(honeybee.abund) +
                          (1 | City/Triplet),
                        family = nbinom2,
                        data = abund.100.total)

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
  x = tapply(abund.100.total$X, abund.100.total$Site, mean),
  y = tapply(abund.100.total$Y, abund.100.total$Site, mean))


# Urban
mod.urban <- glmmTMB(total.wild.abund ~ Category  * scale(Urban) + scale(Hab_Area) + scale(Vegetated.open.land) + 
                       scale(Forest) + scale(honeybee.abund) +
                       (1 | City/Triplet),
                     family = nbinom2,
                     data = abund.100.total)

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
  x = tapply(abund.100.total$X, abund.100.total$Site, mean),
  y = tapply(abund.100.total$Y, abund.100.total$Site, mean))

# Vegetated open land
mod.vol <- glmmTMB(total.wild.abund ~ Category  * scale(Vegetated.open.land) + scale(Urban) + scale(Hab_Area) + 
                     scale(Forest) + scale(honeybee.abund) +
                     (1 | City/Triplet),
                   family = nbinom2,
                   data = abund.100.total)

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
  x = tapply(abund.100.total$X, abund.100.total$Site, mean),
  y = tapply(abund.100.total$Y, abund.100.total$Site, mean))

# Forest
mod.forest <- glmmTMB(total.wild.abund ~ Category  * scale(Forest) + scale(Urban) + scale(Hab_Area) + 
                        scale(Vegetated.open.land) + scale(honeybee.abund) +
                        (1 | City/Triplet),
                      family = nbinom2,
                      data = abund.100.total)

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
  x = tapply(abund.100.total$X, abund.100.total$Site, mean),
  y = tapply(abund.100.total$Y, abund.100.total$Site, mean))

# Honey bee
mod.honey <- glmmTMB(total.wild.abund ~ Category  * scale(honeybee.abund) + scale(Forest) + scale(Urban) + scale(Hab_Area) + 
                        scale(Vegetated.open.land)  +
                        (1 | City/Triplet),
                      family = nbinom2,
                      data = abund.100.total)

summary(mod.honey)
mod_dharma1 <- mod.honey %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.honey) 
AICc(mod.honey) 
plot(allEffects(mod.honey))
visreg(mod.honey, scale = "response")

res.100 <- simulateResiduals(mod.honey)

testSpatialAutocorrelation(
  res.100,
  x = tapply(abund.100.total$X, abund.100.total$Site, mean),
  y = tapply(abund.100.total$Y, abund.100.total$Site, mean))

# Habitat area * Urban cover
mod.area.urban <- glmmTMB(total.wild.abund ~ Category + scale(Hab_Area) * scale(Urban) + scale(Vegetated.open.land) + 
                            scale(Forest) + scale(honeybee.abund) +
                            (1 | City/Triplet),
                          family = nbinom2,
                          data = abund.100.total)

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
  x = tapply(abund.100.total$X, abund.100.total$Site, mean),
  y = tapply(abund.100.total$Y, abund.100.total$Site, mean))

# Habitat area * Vegetated open land
mod.area.vol <- glmmTMB(total.wild.abund ~ Category + scale(Hab_Area) * scale(Vegetated.open.land) + scale(Urban) +
                          scale(Forest) + scale(honeybee.abund) +
                          (1 | City/Triplet),
                        family = nbinom2,
                        data = abund.100.total)

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
  x = tapply(abund.100.total$X, abund.100.total$Site, mean),
  y = tapply(abund.100.total$Y, abund.100.total$Site, mean))

# Habitat area * Forest
mod.area.forest <- glmmTMB(total.wild.abund ~Category + scale(Hab_Area) * scale(Forest) + scale(Vegetated.open.land) + 
                             scale(Urban) + scale(honeybee.abund) +
                             (1 | City/Triplet),
                           family = nbinom2,
                           data = abund.100.total)

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
  x = tapply(abund.100.total$X, abund.100.total$Site, mean),
  y = tapply(abund.100.total$Y, abund.100.total$Site, mean))

# Collect AICc values

mods <- list(
  base = mod.nested.100,
  urban = mod.urban,
  open = mod.vol,
  forest = mod.forest,
  honeybee = mod.honey,
  area = mod.area.cat, 
  area.urban = mod.area.urban, 
  area.vol = mod.area.vol, 
  area.forest = mod.area.forest
)

map_dfr(mods, AIC, .id = "model")
map_dfr(mods, AICc, .id = "model")


aic_table <- map_dfr(
  mods,
  ~ data.frame(AIC = AICc(.x)),
  .id = "model"
)

aic_table <- aic_table %>%
  arrange(AIC) %>%
  mutate(delta_AIC = AIC - min(AIC))

aic_table

# Landscape analysis - 1 km:  ----
# Base model:
mod.nested.1 <- glmmTMB(total.wild.abund ~ Category + scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + 
                            scale(Forest) + scale(honeybee.abund) +
                            (1 | City/Triplet),
                          family = nbinom2,
                          data = abund.1.total)

summary(mod.nested.1)
mod_dharma1 <- mod.nested.1 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.nested.1) 
AICc(mod.nested.1) 
plot(allEffects(mod.nested.1))

res.1 <- simulateResiduals(mod.nested.1)

testSpatialAutocorrelation(
  res.1,
  x = tapply(abund.1.total$X, abund.1.total$Site, mean),
  y = tapply(abund.1.total$Y, abund.1.total$Site, mean))


# 1 km: test interactions one by one ----
# Category*area
mod.area.cat <- glmmTMB(total.wild.abund ~ Category * scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + 
                          scale(Forest) + scale(honeybee.abund) +
                          (1 | City/Triplet),
                        family = nbinom2,
                        data = abund.1.total)

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
  x = tapply(abund.1.total$X, abund.1.total$Site, mean),
  y = tapply(abund.1.total$Y, abund.1.total$Site, mean))


# Urban
mod.urban <- glmmTMB(total.wild.abund ~ Category  * scale(Urban) + scale(Hab_Area) + scale(Vegetated.open.land) + 
                       scale(Forest) + scale(honeybee.abund) +
                       (1 | City/Triplet),
                     family = nbinom2,
                     data = abund.1.total)

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
  x = tapply(abund.1.total$X, abund.1.total$Site, mean),
  y = tapply(abund.1.total$Y, abund.1.total$Site, mean))

# Vegetated open land
mod.vol <- glmmTMB(total.wild.abund ~ Category  * scale(Vegetated.open.land) + scale(Urban) + scale(Hab_Area) + 
                     scale(Forest) + scale(honeybee.abund) +
                     (1 | City/Triplet),
                   family = nbinom2,
                   data = abund.1.total)

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
  x = tapply(abund.1.total$X, abund.1.total$Site, mean),
  y = tapply(abund.1.total$Y, abund.1.total$Site, mean))

# Forest
mod.forest <- glmmTMB(total.wild.abund ~ Category  * scale(Forest) + scale(Urban) + scale(Hab_Area) + 
                        scale(Vegetated.open.land) + scale(honeybee.abund) +
                        (1 | City/Triplet),
                      family = nbinom2,
                      data = abund.1.total)

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
  x = tapply(abund.1.total$X, abund.1.total$Site, mean),
  y = tapply(abund.1.total$Y, abund.1.total$Site, mean))

# Honey bee
mod.honey <- glmmTMB(total.wild.abund ~ Category  * scale(honeybee.abund) + scale(Forest) + scale(Urban) + scale(Hab_Area) + 
                       scale(Vegetated.open.land)  +
                       (1 | City/Triplet),
                     family = nbinom2,
                     data = abund.1.total)

summary(mod.honey)
mod_dharma1 <- mod.honey %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.honey) 
AICc(mod.honey) 
plot(allEffects(mod.honey))
visreg(mod.honey, scale = "response")

res.1 <- simulateResiduals(mod.honey)

testSpatialAutocorrelation(
  res.1,
  x = tapply(abund.1.total$X, abund.1.total$Site, mean),
  y = tapply(abund.1.total$Y, abund.1.total$Site, mean))

# Habitat area * Urban cover
mod.area.urban <- glmmTMB(total.wild.abund ~ Category + scale(Hab_Area) * scale(Urban) + scale(Vegetated.open.land) + 
                            scale(Forest) + scale(honeybee.abund) +
                            (1 | City/Triplet),
                          family = nbinom2,
                          data = abund.1.total)

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
  x = tapply(abund.1.total$X, abund.1.total$Site, mean),
  y = tapply(abund.1.total$Y, abund.1.total$Site, mean))

# Habitat area * Vegetated open land
mod.area.vol <- glmmTMB(total.wild.abund ~ Category + scale(Hab_Area) * scale(Vegetated.open.land) + scale(Urban) +
                          scale(Forest) + scale(honeybee.abund) +
                          (1 | City/Triplet),
                        family = nbinom2,
                        data = abund.1.total)

summary(mod.area.vol)
mod_dharma1 <- mod.area.vol %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(mod.area.vol) 
AICc(mod.area.vol) 
plot(allEffects(mod.area.vol))
visreg(mod.area.vol, scale = "response")

res.1 <- simulateResiduals(mod.area.vol)

testSpatialAutocorrelation(
  res.1,
  x = tapply(abund.1.total$X, abund.1.total$Site, mean),
  y = tapply(abund.1.total$Y, abund.1.total$Site, mean))

# Habitat area * Forest
mod.area.forest <- glmmTMB(total.wild.abund ~Category + scale(Hab_Area) * scale(Forest) + scale(Vegetated.open.land) + 
                             scale(Urban) + scale(honeybee.abund) +
                             (1 | City/Triplet),
                           family = nbinom2,
                           data = abund.1.total)

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
  x = tapply(abund.1.total$X, abund.1.total$Site, mean),
  y = tapply(abund.1.total$Y, abund.1.total$Site, mean))

# Collect AICc values

mods <- list(
  base = mod.nested.1,
  urban = mod.urban,
  open = mod.vol,
  forest = mod.forest,
  honeybee = mod.honey,
  area = mod.area.cat, 
  area.urban = mod.area.urban, 
  area.vol = mod.area.vol, 
  area.forest = mod.area.forest
)

map_dfr(mods, AIC, .id = "model")
map_dfr(mods, AICc, .id = "model")


aic_table <- map_dfr(
  mods,
  ~ data.frame(AIC = AICc(.x)),
  .id = "model"
)

aic_table <- aic_table %>%
  arrange(AIC) %>%
  mutate(delta_AIC = AIC - min(AIC))

aic_table

# plot stuff

abund_long <- abund.1.total %>%
  dplyr::select(Category,
         total.wild.abund,
         honeybee.abund) %>%
  pivot_longer(
    cols = c(total.wild.abund, honeybee.abund),
    names_to = "Group",
    values_to = "Abundance"
  ) %>%
  mutate(
    Group = dplyr::recode(
      Group,
      total.wild.abund = "Wild pollinators",
      honeybee.abund = "Honey bees"
    )
  )

ggplot(abund_long,
       aes(x = Category,
           y = Abundance,
           fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  theme_bw() +
  labs(
    x = "Habitat type",
    y = "Abundance",
    fill = ""
  ) +
  theme(
    axis.text.x = element_text(angle = 15, hjust = 1)
  )


ggplot(abund_long,
       aes(x = Category,
           y = Abundance,
           fill = Category)) +
  geom_boxplot() +
  facet_wrap(~ Group, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Habitat type",
    y = "Abundance"
  ) +
  theme(
    axis.text.x = element_text(angle = 15, hjust = 1),
    legend.position = "none"
  )

# 

# Pollinator groups ----
# structure data
# 100 m:
honeybee_df <- abund.100.raw %>%
  group_by(Site) %>%
  summarise(honeybee = sum(honeybee_abundance, na.rm = TRUE), .groups = "drop")


abund.100.groups <- abund.100.raw %>%
  mutate(Pollinator.type = ifelse(Pollinator.type == "Burnet moth", "Butterfly",Pollinator.type)) %>%
  filter(Pollinator.type %in% c("Bee", "Bumblebee", "Butterfly", "Hoverfly")) |> 
  group_by(City, Site, Triplet, Category, Pollinator.type, X, Y, Hab_Area, 
           Arable, Forest, Non.vegetated.open.land, Urban, Vegetated.open.land, Water, Wetland) |> 
  summarise(wild.pollinator.abund = sum(total_abundance))
abund.100.groups
  
group.db.100 <- left_join(honeybee_df, abund.100.groups, by = c("Site"))

# Landscape analysis - pollinator groups 100 m: ----
# base model
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


# Test interaction Category * pollinator group
pollinator <- glmmTMB(wild.pollinator.abund ~ Category * Pollinator.type + scale(honeybee) + scale(Hab_Area) + 
                            scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                            (1 | City/Triplet),
                          family = nbinom2,
                          data = group.db.100)

summary(pollinator)
mod_dharma1 <- pollinator %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(pollinator) 
AICc(pollinator)
plot(allEffects(pollinator))

# Test interaction Pollinator group * honeybee abundance
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

# 1 km: 
honeybee_df <- abund.1.raw %>%
  group_by(Site) %>%
  summarise(honeybee = sum(honeybee_abundance, na.rm = TRUE), .groups = "drop")


abund.1.groups <- abund.1.raw %>%
  mutate(Pollinator.type = ifelse(Pollinator.type == "Burnet moth", "Butterfly",Pollinator.type)) %>%
  filter(Pollinator.type %in% c("Bee", "Bumblebee", "Butterfly", "Hoverfly")) |> 
  group_by(City, Site, Triplet, Category, Pollinator.type, X, Y, Hab_Area, 
           Arable, Forest, Non.vegetated.open.land, Urban, Vegetated.open.land, Water, Wetland) |> 
  summarise(wild.pollinator.abund = sum(total_abundance))
abund.1.groups

group.db.1 <- left_join(honeybee_df, abund.1.groups, by = c("Site"))

# Landscape analysis - pollinator groups 1 km: ----
# base model
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


# Test interaction Category * pollinator group
pollinator <- glmmTMB(wild.pollinator.abund ~ Category * Pollinator.type + scale(honeybee) + scale(Hab_Area) + 
                        scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                        (1 | City/Triplet),
                      family = nbinom2,
                      data = group.db.1)

summary(pollinator)
mod_dharma1 <- pollinator %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(pollinator) 
AICc(pollinator)
plot(allEffects(pollinator))

# Test interaction Pollinator group * honeybee abundance
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
