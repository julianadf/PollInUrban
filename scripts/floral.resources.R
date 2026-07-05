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
library(corrplot)
library(Hmisc)

# Load data
floral.raw <- read.csv2(here("data", "Floral.resources.csv" ), fileEncoding = "latin1")

cover_lookup <- c(
  "0-1" = 0.5,
  "1-5" = 3,
  "5-25" = 15,
  "25-50" = 37.5,
  "50-75" = 62.5,
  "75-100" = 87.5
)

floral.site <- floral.raw %>%
  mutate(
    Cover = cover_lookup[Blossom.cover.percent]
  ) %>%
  group_by(City, Site, Triplet, Category, Visit) %>%
  summarise(
    Floral_cover = sum(Cover, na.rm = TRUE),
    .groups = "drop"
  )


floral.summary <- floral.site %>%
  group_by(Category, Visit) %>%
  summarise(
    mean_cover = mean(Floral_cover),
    sd = sd(Floral_cover),
    n = n(),
    se = sd / sqrt(n),
    .groups = "drop"
  )


ggplot(floral.summary,
       aes(x = factor(Visit),
           y = mean_cover,
           colour = Category,
           group = Category)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean_cover - se,
        ymax = mean_cover + se),
    width = 0.15
  ) +
  theme_bw() +
  labs(
    x = "Visit",
    y = "Total floral cover",
    colour = "Habitat type",
    title = "Seasonal dynamics of floral resources"
  )


ggplot(floral.site,
       aes(x = factor(Visit),
           y = Floral_cover,
           fill = Category)) +
  geom_boxplot() +
  facet_wrap(~ Category) +
  theme_bw() +
  labs(
    x = "Visit",
    y = "Total floral cover",
    title = "Floral resources across the season"
  )



floral.richness <- floral.raw %>%
  group_by(City, Site, Triplet, Category, Visit) %>%
  summarise(
    plant_richness = n_distinct(Plant.species),
    .groups = "drop"
  )


ggplot(floral.richness,
       aes(x = factor(Visit),
           y = plant_richness,
           colour = Category,
           group = Category)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  theme_bw() +
  labs(
    y = "Plant species richness",
    x = "Visit"
  )


# Standardise efter transect length


floral_transect <- floral.raw %>%
  mutate(
    Cover = cover_lookup[Blossom.cover.percent]
  ) %>%
  group_by(City, Site, Triplet, Category, Visit, Transect) %>%
  summarise(
    Transect_length = first(Transect.length.m),
    Floral_cover = sum(Cover, na.rm = TRUE),
    Floral_cover_m = Floral_cover / Transect_length,
    .groups = "drop"
  )


floral_site <- floral_transect %>%
  group_by(City, Site, Triplet, Category, Visit) %>%
  summarise(
    Floral_cover_m = mean(Floral_cover_m),
    .groups = "drop"
  )



floral_summary <- floral_site %>%
  group_by(Category, Visit) %>%
  summarise(
    mean_cover = mean(Floral_cover_m),
    se = sd(Floral_cover_m) / sqrt(n()),
    .groups = "drop"
  )


ggplot(
  floral_summary,
  aes(x = factor(Visit),
      y = mean_cover,
      colour = Category,
      group = Category)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean_cover - se,
        ymax = mean_cover + se),
    width = 0.15
  ) +
  theme_bw() +
  labs(
    x = "Visit",
    y = "Floral cover per metre transect",
    colour = "Habitat type",
    title = "Seasonal dynamics of floral resources"
  )



ggplot(
  floral_site,
  aes(x = factor(Visit),
      y = Floral_cover_m,
      fill = Category)
) +
  geom_boxplot() +
  facet_wrap(~ Category) +
  theme_bw() +
  labs(
    x = "Visit",
    y = "Floral cover per metre transect",
    title = "Floral resource availability through the season"
  )


