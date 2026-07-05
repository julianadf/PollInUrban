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
library(betapart)

# Load data
community <- read.csv2(here("data/databases", "community.nohoney.csv"), fileEncoding = "latin1")

comm.matrix <- community %>%
  dplyr::select(-City, -Site, -Triplet, -Category)


meta <- community %>%
  dplyr::select(Site, Triplet, Category)

rownames(comm.matrix) <- community$Site

set.seed(123)

nmds <- metaMDS(comm.matrix, distance = "bray", k = 2,trymax = 100)
nmds$stress #0.2524

# Need to remove rare species due to high stress. There are 196 species in total. Remove species seen 
# in less than 3 sites (96)
comm.matrix <- comm.matrix[, colSums(comm.matrix > 0) >= 3]

set.seed(123)

nmds <- metaMDS(comm.matrix, distance = "bray", k = 2,trymax = 100)
nmds$stress #0.27756

# Still high stress. Transform data
comm.hel <- decostand(comm.matrix, method = "hellinger")

set.seed(123)

nmds <- metaMDS(comm.hel, distance = "bray", k = 2, trymax = 100)
nmds$stress #0.2432

scores.comm <- as.data.frame(scores(nmds, display = "sites"))
scores.comm$Site <- rownames(comm.matrix)

scores.comm <- scores.comm %>%
  left_join(meta, by = "Site")


ggplot(scores.comm, aes(NMDS1, NMDS2, color = Category)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  theme_bw()


ggplot(scores.comm, aes(NMDS1, NMDS2, color = Category)) +
  geom_point(size = 3) +
  geom_line(aes(group = Triplet), alpha = 0.4) +
  theme_bw()

