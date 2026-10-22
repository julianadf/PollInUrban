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

## Create database 
# Load data
transects.raw <- read.csv2(here("data", "Transect.csv" ), fileEncoding = "latin1")
bee.redl.raw <- read.csv2(here("data", "Bees Skåne redlist.csv" ), fileEncoding = "latin1")
butterfly.redl.raw <- read.csv2(here("data", "Butterflies Skåne_filtered_traits_redlisted.csv" ), fileEncoding = "latin1")

# Select the trait variables and match species in our study
bee.red <- bee.redl.raw |> 
  dplyr::select("Species", "RL2025", "Pollinator.group")
bee.red

butterfly.red <- butterfly.redl.raw |> 
  dplyr::select("Species", "RL2025", "Pollinator.group")
butterfly.red


red.list <- bind_rows(bee.red, butterfly.red)

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

transects.bb <- transects |> 
  filter(Pollinator.type %in% c("Bee", "Bumblebee", "Butterfly", "Burnet moth")) |> 
  mutate(Pollinator.type = ifelse(Pollinator.type == "Burnet moth", "Butterfly", Pollinator.type)) 
transects.bf.bees

unmatched_species <- transects.bb %>%
distinct(Taxon) %>%
anti_join(red.list,by = c("Taxon" = "Species"))
unmatched_species

# Unmatched species:
# 1  Bombus lucorum coll.
# 2     Anthophora retusa = LC
# 3         Osmia bicolor = LC
# 4                Nomada
# 5              Wild bee
# 6            Anthophora
# 7               Andrena
# 8            Halictinae
# 9            Halictidae
# 10               Bombus
# 11             Pieridae
# 12             Leptidea
# 13            Sphecodes
# 14         Lasioglossum
# 15               Pieris
# 16            Seladonia
# 17                Osmia
# 18           Lycaenidae
# 19  Glaucopsyche alexis = LC
# 20         Polyommatini
# 21             Argynnis
# 22          Nymphalidae
# 23          Coenonympha
# 24              Hylaeus
# 25         Megachilidae

# No need to complete.

# Check species that match
matched_species <- transects.bb %>%
distinct(Taxon) %>% 
  semi_join(red.list, by = c("Taxon" = "Species"))
matched_species
nrow(matched_species)

#create database
transects.redlist <- transects.bb %>%
  left_join(red.list, by = c("Taxon" = "Species"))

redlisted.only <- transects.redlist %>%
  filter(!is.na(RL2025),RL2025 != "")

redlisted.only %>%
  summarise(
    observations = n(),
    total_abundance = sum(Abund),
    n_species = n_distinct(Taxon),
    n_sites = n_distinct(Site)
  )


redlisted.only %>%
  group_by(Taxon) %>%
  summarise(
    abundance = sum(Abund),
    sites = n_distinct(Site)
  ) %>%
  arrange(desc(abundance))

