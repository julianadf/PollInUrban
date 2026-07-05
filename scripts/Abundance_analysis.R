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


# Abundance: with honey bees ----
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


# 100 m
mod.100 <- glmmTMB(abundance_total ~ Category + scale(Area) * scale(Urban) + scale(Vegetated.open.land) + scale(SNG) +
                     (1 | City),
                   family = nbinom2, data = sr.100)


summary(mod.100)
mod_dharma1 <- mod.100 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(mod.100) 
#check_collinearity(mod.100) 
plot(allEffects(mod.100))
visreg(mod.100, scale = "response")

res.100 <- simulateResiduals(mod.100)

testSpatialAutocorrelation(
  res.100,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))

# 250 m 
mod.250 <- glmmTMB(abundance_total ~ Category + scale(Area) * scale(Urban) + scale(Vegetated.open.land) + scale(SNG) +
                     (1 | City),
                   family = nbinom2, data = sr.250)


summary(mod.250)
mod_dharma1 <- mod.250 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(mod.250) 
#check_collinearity(mod.250) 
plot(allEffects(mod.250))
visreg(mod.250, scale = "response")

res.250 <- simulateResiduals(mod.250)

testSpatialAutocorrelation(
  res.250,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))


# 500 m 
mod.500 <- glmmTMB(abundance_total ~ Category + scale(Area) * scale(Urban) + scale(Vegetated.open.land) + scale(SNG) +
                     (1 | City),
                   family = nbinom2, data = sr.500)


summary(mod.500)
mod_dharma1 <- mod.500 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(mod.500) 
#check_collinearity(mod.250) 
plot(allEffects(mod.500))
visreg(mod.500, scale = "response")

res.500 <- simulateResiduals(mod.500)

testSpatialAutocorrelation(
  res.500,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))


# 1 km 
mod.1000 <- glmmTMB(abundance_total ~ Category + scale(Area) * scale(Urban) + scale(Vegetated.open.land) + scale(SNG) +
                      (1 | City),
                    family = nbinom2, data = sr.1)


summary(mod.1000)
mod_dharma1 <- mod.1000 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(mod.1000) 
#check_collinearity(mod.100) 
plot(allEffects(mod.1000))
visreg(mod.100, scale = "response")

res.1000 <- simulateResiduals(mod.1000)

testSpatialAutocorrelation(
  res.1000,
  x = tapply(sr.1$X, sr.1$Site, mean),
  y = tapply(sr.1$Y, sr.1$Site, mean))


# Abundance: without honey bees ----
# Load data
sr.100.raw <- read.csv2(here("data/databases", "tot.nohoney.richness.100m.csv" ), fileEncoding = "latin1")
sr.250.raw <- read.csv2(here("data/databases", "tot.nohoney.richness.250m.csv" ), fileEncoding = "latin1")
sr.500.raw <- read.csv2(here("data/databases", "tot.nohoney.richness.500m.csv" ), fileEncoding = "latin1")
sr.1.raw <- read.csv2(here("data/databases", "tot.nohoney.richness.1km.csv" ), fileEncoding = "latin1")

# Merge urban category into one consisting of Buildings (51), Constructed, not building (52), road/railway (53)
sr.100.nohoney <- sr.100.raw %>% 
  mutate(urban.road = Urban + Road.Railway)

sr.250.nohoney <- sr.250.raw %>% 
  mutate(urban.road = Urban + Road.Railway)

sr.500.nohoney <- sr.500.raw %>% 
  mutate(urban.road = Urban + Road.Railway)

sr.1.nohoney <- sr.1.raw %>% 
  mutate(urban.road = Urban + Road.Railway)

# 100 m
mod.100 <- glmmTMB(abundance_total ~ Category + scale(Area) + scale(Urban) + scale(Vegetated.open.land) + scale(SNG) +
                     (1 | City),
                   family = nbinom2, data = sr.100.nohoney)


summary(mod.100)
mod_dharma1 <- mod.100 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(mod.100) 
#check_collinearity(mod.100) 
plot(allEffects(mod.100))
visreg(mod.100, scale = "response")

res.100 <- simulateResiduals(mod.100)

testSpatialAutocorrelation(
  res.100,
  x = tapply(sr.100.nohoney$X, sr.100.nohoney$Site, mean),
  y = tapply(sr.100.nohoney$Y, sr.100.nohoney$Site, mean))

# 250 m 
mod.250 <- glmmTMB(abundance_total ~ Category + scale(Area) * scale(Urban) + scale(Vegetated.open.land) + scale(SNG) +
                     (1 | City),
                   family = nbinom2, data = sr.250.nohoney)


summary(mod.250)
mod_dharma1 <- mod.250 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(mod.250) 
#check_collinearity(mod.250) 
plot(allEffects(mod.250))
visreg(mod.250, scale = "response")

res.250 <- simulateResiduals(mod.250)

testSpatialAutocorrelation(
  res.250,
  x = tapply(sr.250.nohoney$X, sr.250.nohoney$Site, mean),
  y = tapply(sr.250.nohoney$Y, sr.250.nohoney$Site, mean))


# 500 m 
mod.500 <- glmmTMB(abundance_total ~ Category + scale(Area) * scale(Urban) + scale(Vegetated.open.land) + scale(SNG) +
                     (1 | City),
                   family = nbinom2, data = sr.500.nohoney)


summary(mod.500)
mod_dharma1 <- mod.500 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(mod.500) 
#check_collinearity(mod.250) 
plot(allEffects(mod.500))
visreg(mod.500, scale = "response")

res.500 <- simulateResiduals(mod.500)

testSpatialAutocorrelation(
  res.500,
  x = tapply(sr.500.nohoney$X, sr.500.nohoney$Site, mean),
  y = tapply(sr.500.nohoney$Y, sr.500.nohoney$Site, mean))


# 1 km 
mod.1000 <- glmmTMB(abundance_total ~ Category + scale(Area) + scale(Urban) + scale(Vegetated.open.land) + scale(SNG) +
                      (1 | City),
                    family = nbinom2, data = sr.1.nohoney)


summary(mod.1000)
mod_dharma1 <- mod.1000 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(mod.1000) 
#check_collinearity(mod.100) 
plot(allEffects(mod.1000))
visreg(mod.100, scale = "response")

res.1000 <- simulateResiduals(mod.1000)

testSpatialAutocorrelation(
  res.1000,
  x = tapply(sr.1.nohoney$X, sr.1.nohoney$Site, mean),
  y = tapply(sr.1.nohoney$Y, sr.1.nohoney$Site, mean))

# Abundance: wild bees only ----
# Load data
sr.100.raw <- read.csv2(here("data/databases", "tot.nohoney.richness.100m.csv" ), fileEncoding = "latin1")
sr.250.raw <- read.csv2(here("data/databases", "tot.nohoney.richness.250m.csv" ), fileEncoding = "latin1")
sr.500.raw <- read.csv2(here("data/databases", "tot.nohoney.richness.500m.csv" ), fileEncoding = "latin1")
sr.1.raw <- read.csv2(here("data/databases", "tot.nohoney.richness.1km.csv" ), fileEncoding = "latin1")

# Merge urban category into one consisting of Buildings (51), Constructed, not building (52), road/railway (53)
sr.100.nohoney <- sr.100.raw %>% 
  mutate(urban.road = Urban + Road.Railway)

sr.250.nohoney <- sr.250.raw %>% 
  mutate(urban.road = Urban + Road.Railway)

sr.500.nohoney <- sr.500.raw %>% 
  mutate(urban.road = Urban + Road.Railway)

sr.1.nohoney <- sr.1.raw %>% 
  mutate(urban.road = Urban + Road.Railway)

# 100 m
mod.100 <- glmmTMB(abundance_total ~ Category + scale(Area) + scale(Urban) + scale(Vegetated.open.land) + scale(SNG) +
                     (1 | City),
                   family = nbinom2, data = sr.100.nohoney)


summary(mod.100)
mod_dharma1 <- mod.100 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(mod.100) 
#check_collinearity(mod.100) 
plot(allEffects(mod.100))
visreg(mod.100, scale = "response")

res.100 <- simulateResiduals(mod.100)

testSpatialAutocorrelation(
  res.100,
  x = tapply(sr.100.nohoney$X, sr.100.nohoney$Site, mean),
  y = tapply(sr.100.nohoney$Y, sr.100.nohoney$Site, mean))

# 250 m 
mod.250 <- glmmTMB(abundance_total ~ Category + scale(Area) * scale(Urban) + scale(Vegetated.open.land) + scale(SNG) +
                     (1 | City),
                   family = nbinom2, data = sr.250.nohoney)


summary(mod.250)
mod_dharma1 <- mod.250 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(mod.250) 
#check_collinearity(mod.250) 
plot(allEffects(mod.250))
visreg(mod.250, scale = "response")

res.250 <- simulateResiduals(mod.250)

testSpatialAutocorrelation(
  res.250,
  x = tapply(sr.250.nohoney$X, sr.250.nohoney$Site, mean),
  y = tapply(sr.250.nohoney$Y, sr.250.nohoney$Site, mean))


# 500 m 
mod.500 <- glmmTMB(abundance_total ~ Category + scale(Area) * scale(Urban) + scale(Vegetated.open.land) + scale(SNG) +
                     (1 | City),
                   family = nbinom2, data = sr.500.nohoney)


summary(mod.500)
mod_dharma1 <- mod.500 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(mod.500) 
#check_collinearity(mod.250) 
plot(allEffects(mod.500))
visreg(mod.500, scale = "response")

res.500 <- simulateResiduals(mod.500)

testSpatialAutocorrelation(
  res.500,
  x = tapply(sr.500.nohoney$X, sr.500.nohoney$Site, mean),
  y = tapply(sr.500.nohoney$Y, sr.500.nohoney$Site, mean))


# 1 km 
mod.1000 <- glmmTMB(abundance_total ~ Category + scale(Area) + scale(Urban) + scale(Vegetated.open.land) + scale(SNG) +
                      (1 | City),
                    family = nbinom2, data = sr.1.nohoney)


summary(mod.1000)
mod_dharma1 <- mod.1000 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(mod.1000) 
#check_collinearity(mod.100) 
plot(allEffects(mod.1000))
visreg(mod.100, scale = "response")

res.1000 <- simulateResiduals(mod.1000)

testSpatialAutocorrelation(
  res.1000,
  x = tapply(sr.1.nohoney$X, sr.1.nohoney$Site, mean),
  y = tapply(sr.1.nohoney$Y, sr.1.nohoney$Site, mean))


# Abundance: group wise (without honey bees and ignoring group "other")
# Load data
ab.100.raw <- read.csv2(here("data/databases", "tot.group.abund.100m.csv" ), fileEncoding = "latin1")
ab.250.raw <- read.csv2(here("data/databases", "tot.group.abund.250m.csv" ), fileEncoding = "latin1")
ab.500.raw <- read.csv2(here("data/databases", "tot.group.abund.500m.csv" ), fileEncoding = "latin1")
ab.1.raw <- read.csv2(here("data/databases", "tot.group.abund.1km.csv" ), fileEncoding = "latin1")

# take away group "other" & merge urban category into one consisting of Buildings (51), Constructed, not building (52), road/railway (53)
ab.100.groups <- ab.100.raw %>% 
  mutate(urban.road = Urban + Road.Railway) %>% 
  filter(Pollinator.type != "Other") %>% 
  filter(Pollinator.type != "Empty")

ab.250.groups <- ab.250.raw %>% 
  mutate(urban.road = Urban + Road.Railway) %>% 
  filter(Pollinator.type != "Other")%>% 
  filter(Pollinator.type != "Empty")

ab.500.groups <- ab.500.raw %>% 
  mutate(urban.road = Urban + Road.Railway) %>% 
  filter(Pollinator.type != "Other") %>% 
  filter(Pollinator.type != "Empty")

ab.1.groups <- ab.1.raw %>% 
  mutate(urban.road = Urban + Road.Railway) %>% 
  filter(Pollinator.type != "Other")%>% 
  filter(Pollinator.type != "Empty")

# 100 m
mod.100 <- glmmTMB(abundance_total ~ Pollinator.type * Category + scale(Area) + scale(Urban) + scale(Vegetated.open.land) + scale(SNG) +
                     (1 | City),
                   family = nbinom2, data = ab.100.groups)


summary(mod.100)
mod_dharma1 <- mod.100 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(mod.100) 
#check_collinearity(mod.100) 
plot(allEffects(mod.100))
visreg(mod.100, scale = "response")

res.100 <- simulateResiduals(mod.100)

testSpatialAutocorrelation(
  res.100,
  x = tapply(sr.100.nohoney$X, sr.100.nohoney$Site, mean),
  y = tapply(sr.100.nohoney$Y, sr.100.nohoney$Site, mean))

# 250 m 
mod.250 <- glmmTMB(abundance_total ~ Category + scale(Area) * scale(Urban) + scale(Vegetated.open.land) + scale(SNG) +
                     (1 | City),
                   family = nbinom2, data = ab.250.groups)


summary(mod.250)
mod_dharma1 <- mod.250 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(mod.250) 
#check_collinearity(mod.250) 
plot(allEffects(mod.250))
visreg(mod.250, scale = "response")

res.250 <- simulateResiduals(mod.250)

testSpatialAutocorrelation(
  res.250,
  x = tapply(sr.250.nohoney$X, sr.250.nohoney$Site, mean),
  y = tapply(sr.250.nohoney$Y, sr.250.nohoney$Site, mean))


# 500 m 
mod.500 <- glmmTMB(abundance_total ~ Category + scale(Area) * scale(Urban) + scale(Vegetated.open.land) + scale(SNG) +
                     (1 | City),
                   family = nbinom2, data = ab.500.groups)


summary(mod.500)
mod_dharma1 <- mod.500 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(mod.500) 
#check_collinearity(mod.250) 
plot(allEffects(mod.500))
visreg(mod.500, scale = "response")

res.500 <- simulateResiduals(mod.500)

testSpatialAutocorrelation(
  res.500,
  x = tapply(sr.500.nohoney$X, sr.500.nohoney$Site, mean),
  y = tapply(sr.500.nohoney$Y, sr.500.nohoney$Site, mean))


# 1 km 
mod.1000 <- glmmTMB(abundance_total ~ Category + scale(Area) + scale(Urban) + scale(Vegetated.open.land) + scale(SNG) +
                      (1 | City),
                    family = nbinom2, data = ab.1.groups)


summary(mod.1000)
mod_dharma1 <- mod.1000 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(mod.1000) 
#check_collinearity(mod.100) 
plot(allEffects(mod.1000))
visreg(mod.100, scale = "response")

res.1000 <- simulateResiduals(mod.1000)

testSpatialAutocorrelation(
  res.1000,
  x = tapply(sr.1.nohoney$X, sr.1.nohoney$Site, mean),
  y = tapply(sr.1.nohoney$Y, sr.1.nohoney$Site, mean))

