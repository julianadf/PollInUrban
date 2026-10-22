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
library(ggimage)
library(ggeffects)

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
  filter(Pollinator.type %in% c("Wild bee", "Bumblebee", "Butterfly", "Hoverfly")) |> 
  group_by(City, Site, Triplet, Category, Pollinator.type, X, Y, Hab_Area, 
           Arable, Forest, Non.vegetated.open.land, Urban, Vegetated.open.land, Water, Wetland) |> 
  summarise(wild.pollinator.abund = sum(total_abundance), .groups = "drop_last")
abund.100.groups

group.db.100 <- left_join(honeybee_df, abund.100.groups, by = c("Site"))
group.db.100 <- group.db.100 |> 
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


# 2. FINAL::: Test interaction Category * pollinator group
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


# urban
fig.urban <- visreg(
  pollinator.honey,
  "Urban",
  scale = "response",
  gg = TRUE
) +
  labs(
    x = "Urban cover (100 m)",
    y = "Wild pollinator abundance"
  ) +
  theme_bw()

fig.urban

fig.urban + theme(axis.text.x = element_text(size=14),
               axis.text.y= element_text(size=17),
               axis.title.x = element_text(size=16),
               axis.title.y = element_text(size=17),
               strip.text = element_text(size=20))


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


# 11. Null model
null.group.100 <- glmmTMB(wild.pollinator.abund ~ 1 + 
                            (1 | City/Triplet),
                          family = nbinom2,
                          data = group.db.100)

summary(null.group.100)
mod_dharma1 <- null.group.100 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(null.group.100) 
AICc(null.group.100)
plot(allEffects(null.group.100))

res.100 <- simulateResiduals(null.group.100)
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
  area.forest = area.forest, 
  null = null.group.100
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
  filter(Pollinator.type %in% c("Wild bee", "Bumblebee", "Butterfly", "Hoverfly")) |> 
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

# 6. Final model: Test interaction Pollinator group * vegetated open land
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

# Figures

# Post-hoc: pollinator abundance per habitat type
Pvol.em <- emmeans(pollinator.vol, ~ Category , type="response") # probabilities
Pvol.emms <- emmeans(pollinator.vol, "Category") # estimated marginal means
# Effect size:
print(eff_size(Pvol.em, sigma=sigma(pollinator.vol), edf=df.residual(pollinator.vol))) #
eff_size(pairs(Pvol.em), sigma(pollinator.vol), df.residual(pollinator.vol), method = "identity") # another way of writing the above
# 
pairs(Pvol.emms)
plot(Pvol.emms, comparisons=TRUE)

# Adam Flöhr:
em.Pvol <- as.data.frame(Pvol.em)
ggplot(em.Pvol, aes(Category, response)) +
  geom_point() +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1) +
  ylab("Wild pollinator abundance") +
  xlab(" ") +
  #coord_cartesian(ylim = c(0,1)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

em <- multcomp::cld(Pvol.em, Letters = letters)
em <- as.data.frame(em)
gg.abund <- ggplot(em, aes(Category, response)) +
  geom_point() +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) +#, width = 0.1) +
  #coord_cartesian(ylim = c(0,1.1)) +
  geom_text(aes(y = 34, label = .group)) +
  xlab(" ") +
  ylab("Pollinator abundance")+
  theme_bw() +
  theme(axis.text.x = element_text(size=14),
        axis.text.y= element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16)
  )
gg.abund

# Post-hoc: Abundance pollinator type
pollinType.em <- emmeans(pollinator.vol, ~  Pollinator.type*Vegetated.open.land, type="response") # probabilities
pollinType.emms <- emmeans(pollinator.vol, "Pollinator.type") # estimated marginal means
# Effect size:
print(eff_size(pollinType.em, sigma=sigma(pollinator.vol), edf=df.residual(pollinator.vol))) #
eff_size(pairs(pollinType.em), sigma(pollinator.vol), df.residual(pollinator.vol), method = "identity") # another way of writing the above
# 
pairs(pollinType.emms)
plot(pollinType.emms, comparisons=TRUE)

# Adam Flöhr:
em.Ptype <- as.data.frame(pollinType.em)
ggplot(em.Ptype, aes(Pollinator.type, response)) +
  geom_point() +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1) +
  #ylab("Wild pollinator abundance") +
  xlab(" ") +
  #coord_cartesian(ylim = c(0,1)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

em <- multcomp::cld(pollinType.em, Letters = letters)
em <- as.data.frame(em)
gg.Ptype <- ggplot(em, aes(Pollinator.type, response)) +
  geom_point() +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) +#, width = 0.1) +
  #coord_cartesian(ylim = c(0,1.1)) +
  geom_text(aes(y = 71, label = .group)) +
  xlab(" ") +
  ylab("Pollinator abundance")+
  theme_bw() +
  theme(axis.text.x = element_text(size=14),
        axis.text.y= element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16)
  )
gg.Ptype


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

# 11. Null model
null.group.1 <- glmmTMB(wild.pollinator.abund ~ 1 +
                          (1 | City/Triplet),
                        family = nbinom2,
                        data = group.db.1)

summary(null.group.1)
mod_dharma1 <- null.group.1 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(null.group.1) 
AICc(null.group.1)
plot(allEffects(null.group.1))

res.1 <- simulateResiduals(null.group.1)
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
  area.forest = area.forest, 
  null = null.group.1
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


# Make nice figure for Erik -----
# Load data
abund <- read.csv2(here("data/databases", "abund.db.100m.csv"), fileEncoding = "latin1")

# Prepare data
abund <- abund.100.raw %>%
  mutate(Pollinator.type = dplyr::recode(Pollinator.type, "Burnet moth" = "Butterfly")) %>%
  filter(Pollinator.type %in%  c("Wild bee", "Honey bee", "Bumblebee", "Butterfly", "Hoverfly"))

# ------------------------------------------
# Pollinator-group means ± 95% CI
# ------------------------------------------

poll.data <- abund %>%
  group_by(Category, Pollinator.type) %>%
  summarise(
    n = n(),
    mean_abund = mean(total_abundance),
    se = sd(total_abundance) / sqrt(n),
    lower95 = pmax(0, mean_abund - 1.96 * se),
    upper95 = mean_abund + 1.96 * se,
    .groups = "drop"
  )

# ------------------------------------------
# Total abundance per site
# ------------------------------------------

total.summary <- abund %>%
  group_by(Site, Category) %>%
  summarise(
    total_abund = sum(total_abundance),
    .groups = "drop"
  ) %>%
  group_by(Category) %>%
  summarise(
    n = n(),
    mean_total = mean(total_abund),
    se = sd(total_abund) / sqrt(n),
    lower95 = pmax(0, mean_total - 1.96 * se),
    upper95 = mean_total + 1.96 * se,
    .groups = "drop"
  )

# ------------------------------------------
# Plot
# ------------------------------------------

poll.cols <- c(
  "Wild bee" = "#F2C14E",
  "Honey bee" = "#E76F51",
  "Bumblebee" = "#5DA5C9",
  "Butterfly" = "#8E63CE",
  "Hoverfly" = "#74C476"
)

ggplot() +
  
  # White bars (total abundance)
  geom_col(
    data = total.summary,
    aes(x = Category, y = mean_total),
    fill = "white",
    colour = "black",
    width = 0.9
  ) +
  
  geom_errorbar(
    data = total.summary,
    aes(
      x = Category,
      ymin = mean_total,
      ymax = upper95
    ),
    width = 0.12,
    linewidth = 0.7
  ) +
  
  # Coloured-bar confidence intervals
  geom_errorbar(
    data = poll.data,
    aes(
      x = Category,
      ymin = mean_abund,
      ymax = upper95,
      group = Pollinator.type
    ),
    position = position_dodge(width = 0.75),
    width = 0.09,
    linewidth = 0.5
  ) +
  
  # Coloured bars
  geom_col(
    data = poll.data,
    aes(
      x = Category,
      y = mean_abund,
      fill = Pollinator.type
    ),
    position = position_dodge(width = 0.75),
    width = 0.75
  ) +
  
  scale_fill_manual(values = poll.cols) +
  
  labs(
    x = "Site type",
    y = "Abundance",
    fill = "Pollinator type:"
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Figures for publication: ----
# Final model (100 m):
cat.pollinator <- glmmTMB(wild.pollinator.abund ~ Category * Pollinator.type + scale(honeybee) + scale(Hab_Area) + 
                            scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                            (1 | City/Triplet),
                          family = nbinom2,
                          data = group.db.100)

summary(cat.pollinator)
car::Anova(cat.pollinator, type = 3)
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


# Post-hoc: pollinator abundance per habitat type
abund.em <- emmeans(cat.pollinator, ~ Pollinator.type * Category, type="response") # probabilities
abund.emms <- emmeans(cat.pollinator, ~ Category | Pollinator.type, type = "response")# estimated marginal means
# Effect size:
print(eff_size(abund.em, sigma=sigma(cat.pollinator), edf=df.residual(cat.pollinator))) #
eff_size(pairs(abund.em), sigma(cat.pollinator), df.residual(cat.pollinator), method = "identity") # another way of writing the above
# 
pairs(abund.emms)
plot(abund.emms, comparisons=TRUE)

# Adam Flöhr:
em.abund <- as.data.frame(abund.em)
ggplot(em.abund, aes(Category, response, by= )) +
  facet_wrap("Pollinator.type", scales = "free_y") +
  geom_point() +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1) +
  ylab("Wild pollinator abundance") +
  xlab(" ") +
  #coord_cartesian(ylim = c(0,1)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

em <- multcomp::cld(abund.emms,Letters = letters)
em <- as.data.frame(em)
Th <- as.data.frame(em)

em$label_y <- em$asymp.UCL * 1.19

gg.abund100m <- ggplot(em, aes(Category, response, colour = Category)) +
  facet_wrap("Pollinator.type", scales = "free_y", 
             labeller = labeller(Pollinator.type = c("Bumblebee" = "Bumblebees", 
                                                     "Butterfly" = "Butterflies",
                                                     "Hoverfly" = "Hoverflies", 
                                                     "Wild bee" = "Wild bees"))) +
  geom_point(size = 3) +
  scale_color_manual(values = c("DPC" = "#440154", "Semi-natural grassland" = "#21908C", "Sown meadow" = "#f46d43")) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1) +
  #coord_cartesian(ylim = c(0,1.1)) +
  geom_text(aes(y = label_y, label = .group), colour = "#000000") +
  xlab(" ") +
  ylab("Pollinator abundance") +
  #scale_x_discrete(labels = c("DPC" = "Designed Plant\nCommunity", "Semi-natural grassland" = "Semi-natural\ngrassland", "Sown meadow" = "Sown\nmeadow")) +
  scale_x_discrete(labels = c("DPC" = "Designed", "Semi-natural grassland" = "Semi-natural", "Sown meadow" = "Sown")) +
  theme_bw() +
  theme(axis.text.x = element_text(size=14, angle = 0, hjust = 0.5, vjust = 1), #, face = "bold"),
        axis.text.y= element_text(size=13), #, face = "bold"),
        #axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=14), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 12),
        legend.position = "none", 
        panel.grid = element_blank())
gg.abund100m

ggsave("Fig2.png", gg.abund100m, width = 24, height = 16, units = "cm",  dpi = 600)

# effect of urban cover (average effect of urban cover on pollinator abundance = Urban
# effect marginalized over habitat type and pollinator group)
urban.emm <- emmeans(cat.pollinator, ~ Urban,  at = list(Urban = seq(min(group.db.100$Urban), max(group.db.100$Urban), length = 50)),
  type = "response")
urban.emm <- as.data.frame(urban.emm)

urban.100m <- ggplot(urban.emm, aes(Urban, response)) + geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL),
    fill = "#4C72B0", alpha = 0.2) +
  #geom_point(data = group.db.100, aes(Urban, wild.pollinator.abund, color=Pollinator.type), alpha = 0.5, size = 2) +
  geom_line(linewidth = 1, colour = "#4C72B0") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 16)) +
  labs(x = "Urban cover (100 m)", y = "Pollinator abundance")
urban.100m

#ggsave("Figure3.png", urban.100m, width = 18, height = 12, units = "cm",  dpi = 600)

# Final model (1 km):interaction Pollinator group * vegetated open land ----
pollinator.vol <- glmmTMB(wild.pollinator.abund ~ Category + Pollinator.type * scale(Vegetated.open.land) +
                            scale(Urban) + scale(Hab_Area) + scale(honeybee) + scale(Forest) + 
                            (1 | City/Triplet),
                          family = nbinom2,
                          data = group.db.1)

summary(pollinator.vol)
car::Anova(pollinator.vol, type = 3)
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

# Figure pollinator type * VOL
vol.emm <- emmeans(pollinator.vol, ~ Pollinator.type * Vegetated.open.land, 
                   at = list(Vegetated.open.land = seq(min(group.db.1$Vegetated.open.land), 
                                                       max(group.db.1$Vegetated.open.land),
                                                       length.out = 100)),
                   type = "response")

vol.emm <- as.data.frame(vol.emm)

gg.vol1km <- ggplot(vol.emm, aes(x = Vegetated.open.land, y = response)) +
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), fill = "#4C72B0", alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1, colour = "#4C72B0") +
  facet_wrap(~ Pollinator.type, scales = "free_y") +
  #geom_point(data = group.db.1, aes(Vegetated.open.land, wild.pollinator.abund), color = "#4C72B0", alpha = 0.5, size = 2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.text = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 16)) +
  labs(x = "Vegetated open land cover (1 km)", y = "Pollinator abundance")

gg.vol1km

#ggsave("FigS5.png", gg.vol1km, width = 18, height = 12, units = "cm",  dpi = 600)













#Add silhouettes
# icons <- data.frame(Pollinator.type = c("Bumblebee", "Butterfly", "Hoverfly", "Wild bee"),
#                     image = c("silhouettes/bumblebee.png", "silhouettes/butterfly.png", "silhouettes/hoverfly.png", "silhouettes/wild bee.png"),
#                     Category = "Sown meadow",
#                     response = c(184, 27, 43, 35))
# 
# 
# gg.abund100m <- ggplot(em, aes(Category, response, colour = Category)) +
#   facet_wrap("Pollinator.type", scales = "free_y", 
#              labeller = labeller(Pollinator.type = c("Bumblebee" = "Bumblebees", 
#                                                      "Butterfly" = "Butterflies",
#                                                      "Hoverfly" = "Hoverflies", 
#                                                      "Wild bee" = "Wild bees"))) +
#   geom_point(size = 3) +
#   scale_color_manual(values = c("DPC" = "#440154", "Semi-natural grassland" = "#21908C", "Sown meadow" = "#f46d43")) +
#   geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1) +
#   #coord_cartesian(ylim = c(0,1.1)) +
#   geom_text(aes(y = label_y, label = .group), colour = "#000000") +
#   xlab(" ") +
#   ylab("Pollinator abundance") +
#   scale_x_discrete(labels = c("DPC" = "DPC", "Semi-natural grassland" = "Semi-natural\ngrassland", "Sown meadow" = "Sown\nmeadow")) +
#   theme_bw() +
#   geom_image(data = icons, aes(x = Category, y = response, image = image), size = 0.19, inherit.aes = FALSE) +
#   theme(axis.text.x = element_text(size=13, angle = 0, hjust = 0.5, vjust = 1, face = "bold"),
#         axis.text.y= element_text(size=13, face = "bold"),
#         #axis.title.x = element_text(size=16),
#         axis.title.y = element_text(size=16), 
#         strip.background = element_blank(), 
#         strip.text = element_text(size = 12),
#         legend.position = "none", 
#         panel.grid = element_blank())
# gg.abund100m

