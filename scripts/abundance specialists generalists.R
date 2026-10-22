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



# Site-level honey bee abundance
honeybee_df <- generalists.specialists.raw %>%
  group_by(City, Site, Triplet, Category) %>%
  summarise(
    honeybee.abund = max(honeybee.abund, na.rm = TRUE),
    .groups = "drop"
  )

# Site-level landscape data
landscape.100 <- landscape.100 %>%
  distinct(City, Site, X, Y, Triplet, Category, Hab_Area, Urban, Forest, Vegetated.open.land)

# Specialist/generalist abundances
sg_abund <- generalists.specialists.raw %>%
  filter(Specialist %in% c("Yes", "No")) %>%
  group_by(City, Site, Triplet, Category, Specialist) %>%
  summarise(abundance.group = sum(Abund, na.rm = TRUE),.groups = "drop")

# Complete all sites × specialist categories: 100 m
specialist.100m <- landscape.100 %>%
  dplyr::select(City, Site, X, Y, Triplet, Category, Hab_Area,
         Urban, Forest, Vegetated.open.land) %>%
  crossing(Specialist = c("Yes", "No")) %>%
  left_join(
    sg_abund,
    by = c(
      "City", "Site", "Triplet",
      "Category", "Specialist")) %>%
  mutate(abundance.group = replace_na(abundance.group, 0)) %>%
  left_join(honeybee_df, by = c("City", "Site", "Triplet", "Category"))

str(specialist.100m)


# Complete all sites × specialist categories: 1 km
specialist.1km <- landscape.1 %>%
  dplyr::select(City, Site, X, Y, Triplet, Category, Hab_Area,
                Urban, Forest, Vegetated.open.land) %>%
  crossing(Specialist = c("Yes", "No")) %>%
  left_join(
    sg_abund,
    by = c(
      "City", "Site", "Triplet",
      "Category", "Specialist")) %>%
  mutate(abundance.group = replace_na(abundance.group, 0)) %>%
  left_join(honeybee_df, by = c("City", "Site", "Triplet", "Category"))

str(specialist.1km)

# Analyse
# 100 m: generalist yes/no ----
base.specialist.100m <- glmmTMB(abundance.group ~ Category + Specialist + scale(honeybee.abund) + 
                                  scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                                  (1 | City/Triplet),
                        family = nbinom2,
                        data = specialist.100m)

summary(base.specialist.100m)
mod_dharma1 <- base.specialist.100m %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(base.specialist.100m) 
AICc(base.specialist.100m) 
plot(allEffects(base.specialist.100m))

res.100 <- simulateResiduals(base.specialist.100m)
res.agg <- recalculateResiduals(res.100,  group = specialist.100m$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(specialist.100m$X, specialist.100m$Site, mean),
  y = tapply(specialist.100m$Y, specialist.100m$Site, mean))

# Test interactions one by one
# Best model 100 m: Category * specialist ----
cat.specialist.100m <- glmmTMB(abundance.group ~ Category * Specialist + scale(honeybee.abund) + 
                                  scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                                  (1 | City/Triplet),
                                family = nbinom2,
                                data = specialist.100m)

summary(cat.specialist.100m)
car::Anova(cat.specialist.100m, type = 3)
mod_dharma1 <- cat.specialist.100m %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.specialist.100m) 
AICc(cat.specialist.100m) 
plot(allEffects(cat.specialist.100m))

res.100 <- simulateResiduals(cat.specialist.100m)
res.agg <- recalculateResiduals(res.100,  group = specialist.100m$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(specialist.100m$X, specialist.100m$Site, mean),
  y = tapply(specialist.100m$Y, specialist.100m$Site, mean))


# specialist * honey bee abundance
honey.specialist.100m <- glmmTMB(abundance.group ~ Category + Specialist * scale(honeybee.abund) + 
                                 scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                                 (1 | City/Triplet),
                               family = nbinom2,
                               data = specialist.100m)

summary(honey.specialist.100m)
mod_dharma1 <- honey.specialist.100m %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(honey.specialist.100m) 
AICc(honey.specialist.100m) 
plot(allEffects(honey.specialist.100m))

res.100 <- simulateResiduals(honey.specialist.100m)
res.agg <- recalculateResiduals(res.100,  group = specialist.100m$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(specialist.100m$X, specialist.100m$Site, mean),
  y = tapply(specialist.100m$Y, specialist.100m$Site, mean))

# specialist * habitat area
area.specialist.100m <- glmmTMB(abundance.group ~ Category + Specialist * scale(Hab_Area)+
                                   scale(honeybee.abund) + 
                                    + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                                   (1 | City/Triplet),
                                 family = nbinom2,
                                 data = specialist.100m)

summary(area.specialist.100m)
mod_dharma1 <- area.specialist.100m %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.specialist.100m) 
AICc(area.specialist.100m) 
plot(allEffects(area.specialist.100m))

res.100 <- simulateResiduals(area.specialist.100m)
res.agg <- recalculateResiduals(res.100,  group = specialist.100m$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(specialist.100m$X, specialist.100m$Site, mean),
  y = tapply(specialist.100m$Y, specialist.100m$Site, mean))

# specialist * urban
urban.specialist.100m <- glmmTMB(abundance.group ~ Category + Specialist * scale(Urban) +
                                  scale(honeybee.abund) + scale(Hab_Area)+
                                  + scale(Vegetated.open.land) + scale(Forest) + 
                                  (1 | City/Triplet),
                                family = nbinom2,
                                data = specialist.100m)

summary(urban.specialist.100m)
mod_dharma1 <- urban.specialist.100m %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(urban.specialist.100m) 
AICc(urban.specialist.100m) 
plot(allEffects(urban.specialist.100m))

res.100 <- simulateResiduals(urban.specialist.100m)
res.agg <- recalculateResiduals(res.100,  group = specialist.100m$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(specialist.100m$X, specialist.100m$Site, mean),
  y = tapply(specialist.100m$Y, specialist.100m$Site, mean))

# specialist * vegetated open land
vol.specialist.100m <- glmmTMB(abundance.group ~ Category + Specialist * scale(Vegetated.open.land) +
                                   scale(Urban) + scale(honeybee.abund) + scale(Hab_Area)+ scale(Forest) + 
                                   (1 | City/Triplet),
                                 family = nbinom2,
                                 data = specialist.100m)

summary(vol.specialist.100m)
mod_dharma1 <- vol.specialist.100m %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(vol.specialist.100m) 
AICc(vol.specialist.100m) 
plot(allEffects(vol.specialist.100m))

res.100 <- simulateResiduals(vol.specialist.100m)
res.agg <- recalculateResiduals(res.100,  group = specialist.100m$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(specialist.100m$X, specialist.100m$Site, mean),
  y = tapply(specialist.100m$Y, specialist.100m$Site, mean))

# specialist * forest
forest.specialist.100m <- glmmTMB(abundance.group ~ Category + Specialist * scale(Forest) +
                                 scale(Urban) + scale(honeybee.abund) + scale(Hab_Area)+ scale(Vegetated.open.land) + 
                                 (1 | City/Triplet),
                               family = nbinom2,
                               data = specialist.100m)

summary(forest.specialist.100m)
mod_dharma1 <- forest.specialist.100m %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(forest.specialist.100m) 
AICc(forest.specialist.100m) 
plot(allEffects(forest.specialist.100m))

res.100 <- simulateResiduals(forest.specialist.100m)
res.agg <- recalculateResiduals(res.100,  group = specialist.100m$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(specialist.100m$X, specialist.100m$Site, mean),
  y = tapply(specialist.100m$Y, specialist.100m$Site, mean))


# Null model
null.specialist.100m <- glmmTMB(abundance.group ~ 1 + 
                                  (1 | City/Triplet),
                                family = nbinom2,
                                data = specialist.100m)

# Collect AICc values

mods <- list(
  base = base.specialist.100m,
  category = cat.specialist.100m,
  honey = honey.specialist.100m, 
  area = area.specialist.100m, 
  urban = urban.specialist.100m,
  open = vol.specialist.100m,
  forest = forest.specialist.100m, 
  null = null.specialist.100m
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


# 1 km: generalist yes/no ----
base.specialist.1km <- glmmTMB(abundance.group ~ Category + Specialist + scale(honeybee.abund) + 
                                  scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                                  (1 | City/Triplet),
                                family = nbinom2,
                                data = specialist.1km)

summary(base.specialist.1km)
mod_dharma1 <- base.specialist.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(base.specialist.1km) 
AICc(base.specialist.1km) 
plot(allEffects(base.specialist.1km))

res.1 <- simulateResiduals(base.specialist.1km)
res.agg <- recalculateResiduals(res.1,  group = specialist.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(specialist.1km$X, specialist.1km$Site, mean),
  y = tapply(specialist.1km$Y, specialist.1km$Site, mean))

# Test interactions one by one
# Category * specialist
cat.specialist.1km <- glmmTMB(abundance.group ~ Category * Specialist + scale(honeybee.abund) + 
                                 scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                                 (1 | City/Triplet),
                               family = nbinom2,
                               data = specialist.1km)

summary(cat.specialist.1km)
car::Anova(cat.specialist.1km, type = 3)
mod_dharma1 <- cat.specialist.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(cat.specialist.1km) 
AICc(cat.specialist.1km) 
plot(allEffects(cat.specialist.1km))

res.1 <- simulateResiduals(cat.specialist.1km)
res.agg <- recalculateResiduals(res.1,  group = specialist.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(specialist.1km$X, specialist.1km$Site, mean),
  y = tapply(specialist.1km$Y, specialist.1km$Site, mean))


# Plot
# Post-hoc: pollinator abundance per habitat type
special.em <- emmeans(cat.specialist.1km, ~ Specialist * Category, type="response") # probabilities
special.emms <- emmeans(cat.specialist.1km, ~ Category | Specialist, type = "response")# estimated marginal means
# Effect size:
print(eff_size(special.em, sigma=sigma(cat.specialist.1km), edf=df.residual(cat.specialist.1km))) #
eff_size(pairs(special.em), sigma(cat.specialist.1km), df.residual(cat.specialist.1km), method = "identity") # another way of writing the above
# 
pairs(special.emms)
plot(special.emms, comparisons=TRUE)

# Adam Flöhr:
em.special <- as.data.frame(special.em)
ggplot(em.special, aes(Category, response, by= )) +
  facet_wrap("Specialist", scales = "free_y") +
  geom_point() +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1) +
  ylab("Abundance") +
  xlab(" ") +
  #coord_cartesian(ylim = c(0,1)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

em <- multcomp::cld(special.emms,Letters = letters)
em <- as.data.frame(em)
Th <- as.data.frame(em)

em$label_y <- em$asymp.UCL * 1.1

gg.special1km <- ggplot(em, aes(Category, response, colour = Category)) +
  facet_wrap("Specialist", scales = "free_y", 
             labeller = labeller(Specialist = c("No" = "Generalist", "Yes" = "Specialist"))) + 
  geom_point(size = 3) +
  scale_color_manual(values = c("DPC" = "#440154", "Semi-natural grassland" = "#21908C", "Sown meadow" = "#f46d43")) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1) +
  #coord_cartesian(ylim = c(0,1.1)) +
  geom_text(aes(y = label_y, label = .group), colour = "#000000") +
  xlab(" ") +
  ylab("Abundance") +
  scale_x_discrete(labels = c("DPC" = "Designed", "Semi-natural grassland" = "Semi-natural", "Sown meadow" = "Sown")) +
  theme_bw() +
  theme(axis.text.x = element_text(size=13, angle = 0, hjust = 0.5, vjust = 1),
        axis.text.y= element_text(size=13),
        #axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 12),
        legend.position = "none", 
        panel.grid = element_blank())
gg.special1km

ggsave("Fig4.png", gg.special1km, width = 28, height = 12, units = "cm",  dpi = 600)

# specialist * honey bee abundance
honey.specialist.1km <- glmmTMB(abundance.group ~ Category + Specialist * scale(honeybee.abund) + 
                                   scale(Hab_Area) + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                                   (1 | City/Triplet),
                                 family = nbinom2,
                                 data = specialist.1km)

summary(honey.specialist.1km)
mod_dharma1 <- honey.specialist.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(honey.specialist.1km) 
AICc(honey.specialist.1km) 
plot(allEffects(honey.specialist.1km))

res.1 <- simulateResiduals(honey.specialist.1km)
res.agg <- recalculateResiduals(res.1,  group = specialist.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(specialist.1km$X, specialist.1km$Site, mean),
  y = tapply(specialist.1km$Y, specialist.1km$Site, mean))

# specialist * habitat area
area.specialist.1km <- glmmTMB(abundance.group ~ Category + Specialist * scale(Hab_Area)+
                                  scale(honeybee.abund) + 
                                  + scale(Urban) + scale(Vegetated.open.land) + scale(Forest) + 
                                  (1 | City/Triplet),
                                family = nbinom2,
                                data = specialist.1km)

summary(area.specialist.1km)
mod_dharma1 <- area.specialist.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(area.specialist.1km) 
AICc(area.specialist.1km) 
plot(allEffects(area.specialist.1km))

res.1 <- simulateResiduals(area.specialist.1km)
res.agg <- recalculateResiduals(res.1,  group = specialist.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(specialist.1km$X, specialist.1km$Site, mean),
  y = tapply(specialist.1km$Y, specialist.1km$Site, mean))

# specialist * urban
urban.specialist.1km <- glmmTMB(abundance.group ~ Category + Specialist * scale(Urban) +
                                   scale(honeybee.abund) + scale(Hab_Area)+
                                   + scale(Vegetated.open.land) + scale(Forest) + 
                                   (1 | City/Triplet),
                                 family = nbinom2,
                                 data = specialist.1km)

summary(urban.specialist.1km)
mod_dharma1 <- urban.specialist.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(urban.specialist.1km) 
AICc(urban.specialist.1km) 
plot(allEffects(urban.specialist.1km))

res.1 <- simulateResiduals(urban.specialist.1km)
res.agg <- recalculateResiduals(res.1,  group = specialist.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(specialist.1km$X, specialist.1km$Site, mean),
  y = tapply(specialist.1km$Y, specialist.1km$Site, mean))

# specialist * vegetated open land
vol.specialist.1km <- glmmTMB(abundance.group ~ Category + Specialist * scale(Vegetated.open.land) +
                                 scale(Urban) + scale(honeybee.abund) + scale(Hab_Area)+ scale(Forest) + 
                                 (1 | City/Triplet),
                               family = nbinom2,
                               data = specialist.1km)

summary(vol.specialist.1km)
mod_dharma1 <- vol.specialist.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(vol.specialist.1km) 
AICc(vol.specialist.1km) 
plot(allEffects(vol.specialist.1km))

res.1 <- simulateResiduals(vol.specialist.1km)
res.agg <- recalculateResiduals(res.1,  group = specialist.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(specialist.1km$X, specialist.1km$Site, mean),
  y = tapply(specialist.1km$Y, specialist.1km$Site, mean))

# specialist * forest
forest.specialist.1km <- glmmTMB(abundance.group ~ Category + Specialist * scale(Forest) +
                                    scale(Urban) + scale(honeybee.abund) + scale(Hab_Area)+ scale(Vegetated.open.land) + 
                                    (1 | City/Triplet),
                                  family = nbinom2,
                                  data = specialist.1km)

summary(forest.specialist.1km)
mod_dharma1 <- forest.specialist.1km %>% simulateResiduals(n=1000)
plot(mod_dharma1)
check_collinearity(forest.specialist.1km) 
AICc(forest.specialist.1km) 
plot(allEffects(forest.specialist.1km))

res.1 <- simulateResiduals(forest.specialist.1km)
res.agg <- recalculateResiduals(res.1,  group = specialist.1km$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.agg,
  x = tapply(specialist.1km$X, specialist.1km$Site, mean),
  y = tapply(specialist.1km$Y, specialist.1km$Site, mean))


# Null model
null.specialist.1km <- glmmTMB(abundance.group ~ 1 + 
                                 (1 | City/Triplet),
                               family = nbinom2,
                               data = specialist.1km)

# Collect AICc values

mods <- list(
  base = base.specialist.1km,
  category = cat.specialist.1km,
  honey = honey.specialist.1km, 
  area = area.specialist.1km, 
  urban = urban.specialist.1km,
  open = vol.specialist.1km,
  forest = forest.specialist.1km, 
  null = null.specialist.1km
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



# Create database ----
# Load data
transects.raw <- read.csv2(here("data", "Transect.csv" ), fileEncoding = "latin1")
bee.traits.raw <- read.csv2(here("data", "Bees Skåne traits.csv" ), fileEncoding = "latin1")
butterfly.traits.raw <- read.csv2(here("data", "Butterflies Skåne_filtered_traits.csv" ), fileEncoding = "latin1")

# Select the trait variables and match species in our study
bee.traits <- bee.traits.raw |> 
  dplyr::select("Species", "Lecticism..mono...oligo...poly..", "Pollinator.group") |> 
  rename(Lecticism = Lecticism..mono...oligo...poly..) 
bee.traits

butterfly.traits <- butterfly.traits.raw |> 
  dplyr::select("Species", "Host_spec_genus", "Pollinator.group") |> 
  rename(Lecticism = Host_spec_genus) 
butterfly.traits


traits <- bind_rows(bee.traits, butterfly.traits)

transects <- transects.raw %>% 
  mutate(City = as.factor(City)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Triplet = as.factor(Triplet)) %>% 
  mutate(Category = as.factor(Category)) %>% 
  mutate(Visit = as.factor(Visit)) %>% 
  rename(Abund = Nr.obs) %>% 
  mutate(Order = as.factor(Order)) %>% 
  mutate(Insect.species = trimws(Insect.species)) %>% 
  rename(Taxon = Insect.species)
str(transects)

transects.bf.bees <- transects |> 
  filter(Pollinator.type %in% c("Bee", "Bumblebee", "Butterfly", "Burnet moth")) |> 
  mutate(Pollinator.type = ifelse(Pollinator.type == "Burnet moth", "Butterfly", Pollinator.type)) 
transects.bf.bees


transects_traits <- transects.bf.bees %>%
  left_join(traits, by = c("Taxon" = "Species"))


setdiff(
  levels(transects$Site),
  unique(transects_traits$Site[!is.na(transects_traits$Lecticism)])
)


sum(!is.na(transects_traits$Lecticism)) # 661 non-matches

transects_traits %>%
  filter(is.na(Lecticism)) %>%
  distinct(Taxon) %>%
  arrange(Taxon)






