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


# Load data
sr.100.raw <- read.csv2(here("data/databases", "tot.richness.100m.csv" ), fileEncoding = "latin1")
sr.250.raw <- read.csv2(here("data/databases", "tot.richness.250m.csv" ), fileEncoding = "latin1")
sr.500.raw <- read.csv2(here("data/databases", "tot.richness.500m.csv" ), fileEncoding = "latin1")
sr.1.raw <- read.csv2(here("data/databases", "tot.richness.1km.csv" ), fileEncoding = "latin1")

# Merge urban category into one consisting of Buildings (51), Constructed, not building (52), road/railway (53)
sr.100 <- sr.100.raw %>% 
  mutate(urban.road = Urban + Road.Railway)

sr.250 <- sr.250.raw %>% 
  mutate(urban.road = Urban + Road.Railway)

sr.500 <- sr.500.raw %>% 
  mutate(urban.road = Urban + Road.Railway)

sr.1 <- sr.1.raw %>% 
  mutate(urban.road = Urban + Road.Railway)


# 100 m ----
mod.100 <- glmmTMB(richness ~ Category + scale(Area) * scale(Urban) + scale(Vegetated.open.land) + scale(SNG) +
                      (1 | City),
                    family = poisson, data = sr.100)


summary(mod.100)
mod_dharma1 <- mod.100 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(mod.100) 
#check_collinearity(mod.100) 
plot(allEffects(mod.100))
visreg(mod.100, scale = "response")

res.100 <- simulateResiduals(mod.100)

#res_agg <- recalculateResiduals(res,  group = sr.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.100,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))

# 250 m ----
mod.250 <- glmmTMB(richness ~ Category + scale(Area) * scale(Urban) + scale(Vegetated.open.land) + scale(SNG) +
                     (1 | City),
                   family = poisson, data = sr.250)


summary(mod.250)
mod_dharma1 <- mod.250 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(mod.250) 
#check_collinearity(mod.250) 
plot(allEffects(mod.250))
visreg(mod.250, scale = "response")

res.250 <- simulateResiduals(mod.250)

#res_agg <- recalculateResiduals(res,  group = sr.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.250,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))


# 500 m ----
mod.500 <- glmmTMB(richness ~ Category + scale(Area) * scale(Urban) + scale(Vegetated.open.land) + scale(SNG) +
                     (1 | City),
                   family = poisson, data = sr.500)


summary(mod.500)
mod_dharma1 <- mod.500 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(mod.500) 
#check_collinearity(mod.250) 
plot(allEffects(mod.500))
visreg(mod.500, scale = "response")

res.500 <- simulateResiduals(mod.500)

#res_agg <- recalculateResiduals(res,  group = sr.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.500,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))


# 1 km ----
mod.1000 <- glmmTMB(richness ~ Category + scale(Area) * scale(Urban) + scale(Vegetated.open.land) + scale(SNG) +
                          (1 | City),
                        family = poisson, data = sr.1)


summary(mod.1000)
mod_dharma1 <- mod.1000 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(mod.1000) 
#check_collinearity(mod.100) 
plot(allEffects(mod.1000))
visreg(mod.100, scale = "response")

res.1000 <- simulateResiduals(mod.1000)

#res_agg <- recalculateResiduals(res,  group = sr.1$Site) # used when I have several observations per site. This aggregates residuals per location

testSpatialAutocorrelation(
  res.1000,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))





# Test the effect of city on the model: can remove
m_with_city <- glmmTMB(
  richness ~ Category + Urban + SNG + (1 | City/Triplet/Site),
  family = nbinom2,
  data = rich.100m)

m_no_city <- glmmTMB(
  richness ~ Category + Urban + SNG + (1 | Triplet/Site),
  family = nbinom2,
  data = rich.100m)

AIC(m_with_city, m_no_city)

# Deal with the spatial autocorrelation
