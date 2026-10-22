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
library(forcats)

# Load data
transects.raw <- read.csv2(here("data", "Transect.csv" ), fileEncoding = "latin1")

# Make bar plot for separate pollinator groups
species.abund <- transects.raw %>%
  filter(Taxon.category == "Art") %>%
  mutate(Pollinator.type = ifelse(Pollinator.type == "Burnet moth", "Butterfly",Pollinator.type)) %>%
  filter(Pollinator.type %in% c("Wild bee", "Bumblebee", "Butterfly", "Hoverfly")) %>%
  group_by(Category, Pollinator.type, Insect.species) %>%
  summarise(
    abundance = sum(Nr.obs, na.rm = TRUE),
    .groups = "drop"
  )
species.abund 

# barplot bumblebees
bbees <- species.abund |> filter(Pollinator.type == "Bumblebee") |> 
  mutate(Insect.species = fct_reorder(Insect.species, abundance, .fun = sum))

bumble <- readPNG("silhouettes/bumblebee.png")

abund.bbees <- ggplot(bbees, aes(x = Insect.species, y = abundance, fill = Category)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "DPC" = "#440154",
      "Semi-natural grassland" = "#21908C",
      "Sown meadow" = "#f46d43")) +
  labs(x = "Species", y = "Abundance", fill = "Habitat type") +
  theme_bw() +
  theme(panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.text.y = element_text(, face = "italic"),
    axis.title = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "top")

abund.bbees <- ggdraw(abund.bbees) + draw_image(bumble, x = 0.93, y = 0.04, width = 0.1, height = 0.5, hjust = 1, vjust = 0)
abund.bbees

ggsave("FigSbbees.png", abund.bbees, width = 20, height = 15, units = "cm",  dpi = 600)
  
# barplot butterflies
bflies <- species.abund |> filter(Pollinator.type == "Butterfly") |> 
  mutate(Insect.species = fct_reorder(Insect.species, abundance, .fun = sum))

butterfly <- readPNG("silhouettes/butterfly.png")

abund.bflies <- ggplot(bflies, aes(x = Insect.species, y = abundance, fill = Category)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "DPC" = "#440154",
      "Semi-natural grassland" = "#21908C",
      "Sown meadow" = "#f46d43")) +
  labs(x = "Species", y = "Abundance", fill = "Habitat type") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.text.y = element_text(, face = "italic"),
    axis.title = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "top")  


abund.bflies <- ggdraw(abund.bflies) + draw_image(butterfly, x = 0.99, y = 0.12, width = 0.2, height = 0.2, hjust = 1, vjust = 0)
abund.bflies

ggsave("FigSbflies.png", abund.bflies, width = 20, height = 15, units = "cm",  dpi = 600)

# barplot hoverflies
hflies <- species.abund |> filter(Pollinator.type == "Hoverfly") |> 
  mutate(Insect.species = fct_reorder(Insect.species, abundance, .fun = sum))

hoverfly <- readPNG("silhouettes/hoverfly.png")

abund.hflies <- ggplot(hflies, aes(x = Insect.species, y = abundance, fill = Category)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "DPC" = "#440154",
      "Semi-natural grassland" = "#21908C",
      "Sown meadow" = "#f46d43")) +
  labs(x = "Species", y = "Abundance", fill = "Habitat type") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.text.y = element_text(, face = "italic"),
    axis.title = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "top")  

abund.hflies <- ggdraw(abund.hflies) + draw_image(hoverfly, x = 0.93, y = 0.10, width = 0.15, height = 0.15, hjust = 1, vjust = 0)
abund.hflies

ggsave("FigShflies.png", abund.hflies, width = 20, height = 22, units = "cm",  dpi = 600)


# barplot wild bees
wbees <- species.abund |> filter(Pollinator.type == "Wild bee") |> 
  mutate(Insect.species = fct_reorder(Insect.species, abundance, .fun = sum))

wildbee <- readPNG("silhouettes/wild bee.png")

abund.wbees <- ggplot(wbees, aes(x = Insect.species, y = abundance, fill = Category)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "DPC" = "#440154",
      "Semi-natural grassland" = "#21908C",
      "Sown meadow" = "#f46d43")) +
  labs(x = "Species", y = "Abundance", fill = "Habitat type") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.text.y = element_text(, face = "italic"),
    axis.title = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "top")  

abund.wbees <- ggdraw(abund.wbees) + draw_image(wildbee, x = 0.93, y = 0.10, width = 0.15, height = 0.15, hjust = 1, vjust = 0)
abund.wbees


ggsave("FigSwbees.png", abund.wbees, width = 20, height = 22, units = "cm",  dpi = 600)

# Species accumulation curves -----
accum_data <- transects.raw %>%
  filter(Taxon.category == "Art") %>%
  mutate(Pollinator.type = ifelse(Pollinator.type == "Burnet moth", "Butterfly",Pollinator.type)) %>%
  filter(Pollinator.type %in% c("Wild bee", "Bumblebee", "Butterfly", "Hoverfly")) %>%
  group_by(Site, Category, Pollinator.type, Insect.species) %>%
  summarise(Abundance = sum(Nr.obs), .groups = "drop") %>%
  mutate(Presence = 1)
par(mfrow = c(2, 2))

# bumblebees
bb_data <- accum_data %>%
  filter(Pollinator.type == "Bumblebee")

accum_list <- list()

for(cat in unique(bb_data$Category)) {
  
  dat <- bb_data %>%
    filter(Category == cat)
  
  comm <- dat %>%
    dplyr::select(Site, Insect.species, Presence) %>%
    pivot_wider(
      names_from = Insect.species,
      values_from = Presence,
      values_fill = 0
    )
  
  rownames(comm) <- comm$Site
  comm <- comm[, -1]
  
  acc <- specaccum(comm, method = "random")
  
  accum_list[[cat]] <- data.frame(
    Sites = seq_along(acc$richness),
    Richness = acc$richness,
    SD = acc$sd,
    Category = cat
  )
}

accum_bbees <- bind_rows(accum_list)

# plot that
accum.bbees <- ggplot(accum_bbees, aes(Sites, Richness, colour = Category, fill = Category)) +
  facet_wrap(~Category)+
  geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = c("DPC" = "#440154",
      "Semi-natural grassland" = "#21908C",
      "Sown meadow" = "#f46d43")) +
  scale_fill_manual(values = c("DPC" = "#440154",
      "Semi-natural grassland" = "#21908C",
      "Sown meadow" = "#f46d43")) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(), 
        strip.text = element_text(size = 12)) +
  geom_image(data = icons, aes(x = Sites, y = Richness, image = image), size = 0.19, inherit.aes = FALSE) +
  labs(x = "Number of sites", y = "Accumulated species", colour = "Habitat type", fill = "Habitat type")
accum.bbees

ggsave("accumbbees.png", accum.bbees, width = 20, height = 14, units = "cm",  dpi = 600)


# butterflies
bf_data <- accum_data %>%
  filter(Pollinator.type == "Butterfly")

accum_list <- list()

for(cat in unique(bf_data$Category)) {
  
  dat <- bf_data %>%
    filter(Category == cat)
  
  comm <- dat %>%
    dplyr::select(Site, Insect.species, Presence) %>%
    pivot_wider(
      names_from = Insect.species,
      values_from = Presence,
      values_fill = 0
    )
  
  rownames(comm) <- comm$Site
  comm <- comm[, -1]
  
  acc <- specaccum(comm, method = "random")
  
  accum_list[[cat]] <- data.frame(
    Sites = seq_along(acc$richness),
    Richness = acc$richness,
    SD = acc$sd,
    Category = cat
  )
}

accum_bflies <- bind_rows(accum_list)

# plot that
accum.bflies <- ggplot(accum_bflies, aes(Sites, Richness, colour = Category, fill = Category)) +
  facet_wrap(~Category)+
  geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = c("DPC" = "#440154",
                                 "Semi-natural grassland" = "#21908C",
                                 "Sown meadow" = "#f46d43")) +
  scale_fill_manual(values = c("DPC" = "#440154",
                               "Semi-natural grassland" = "#21908C",
                               "Sown meadow" = "#f46d43")) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_text(size = 12)) +
  labs(x = "Number of sites", y = "Accumulated species", colour = "Habitat type", fill = "Habitat type")

ggsave("accumbflies.png", accum.bflies, width = 20, height = 14, units = "cm",  dpi = 600)


# hoverflies
hf_data <- accum_data %>%
  filter(Pollinator.type == "Hoverfly")

accum_list <- list()

for(cat in unique(hf_data$Category)) {
  
  dat <- hf_data %>%
    filter(Category == cat)
  
  comm <- dat %>%
    dplyr::select(Site, Insect.species, Presence) %>%
    pivot_wider(
      names_from = Insect.species,
      values_from = Presence,
      values_fill = 0
    )
  
  rownames(comm) <- comm$Site
  comm <- comm[, -1]
  
  acc <- specaccum(comm, method = "random")
  
  accum_list[[cat]] <- data.frame(
    Sites = seq_along(acc$richness),
    Richness = acc$richness,
    SD = acc$sd,
    Category = cat
  )
}

accum_hflies <- bind_rows(accum_list)

# plot that
accum.hflies <- ggplot(accum_hflies, aes(Sites, Richness, colour = Category, fill = Category)) +
  facet_wrap(~Category)+
  geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = c("DPC" = "#440154",
                                 "Semi-natural grassland" = "#21908C",
                                 "Sown meadow" = "#f46d43")) +
  scale_fill_manual(values = c("DPC" = "#440154",
                               "Semi-natural grassland" = "#21908C",
                               "Sown meadow" = "#f46d43")) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_text(size = 12)) +
  labs(x = "Number of sites", y = "Accumulated species", colour = "Habitat type", fill = "Habitat type")

ggsave("accumhflies.png", accum.hflies, width = 20, height = 14, units = "cm",  dpi = 600)

# wild bees
wb_data <- accum_data %>%
  filter(Pollinator.type == "Wild bee")

accum_list <- list()

for(cat in unique(wb_data$Category)) {
  
  dat <- wb_data %>%
    filter(Category == cat)
  
  comm <- dat %>%
    dplyr::select(Site, Insect.species, Presence) %>%
    pivot_wider(
      names_from = Insect.species,
      values_from = Presence,
      values_fill = 0
    )
  
  rownames(comm) <- comm$Site
  comm <- comm[, -1]
  
  acc <- specaccum(comm, method = "random")
  
  accum_list[[cat]] <- data.frame(
    Sites = seq_along(acc$richness),
    Richness = acc$richness,
    SD = acc$sd,
    Category = cat
  )
}

accum_wbees <- bind_rows(accum_list)

# plot that
accum.wbees <- ggplot(accum_wbees, aes(Sites, Richness, colour = Category, fill = Category)) +
  facet_wrap(~Category)+
  geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = c("DPC" = "#440154",
                                 "Semi-natural grassland" = "#21908C",
                                 "Sown meadow" = "#f46d43")) +
  scale_fill_manual(values = c("DPC" = "#440154",
                               "Semi-natural grassland" = "#21908C",
                               "Sown meadow" = "#f46d43")) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_text(size = 12)) +
  labs(x = "Number of sites", y = "Accumulated species", colour = "Habitat type", fill = "Habitat type")

ggsave("accumwbees.png", accum.wbees, width = 20, height = 14, units = "cm",  dpi = 600)
