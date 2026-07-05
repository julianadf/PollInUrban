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
community.visit <- read.csv2(here("data/databases/temporal beta div", "community.visit.csv"), fileEncoding = "latin1")

df.all <- read.csv2(here("data/databases/temporal beta div", "community.visit.csv"), fileEncoding = "latin1")

# Species richness per visit per habitat type
rich_per_visit <- community.visit %>%
  mutate(
    richness = rowSums(across(-(City:Visit), ~ . > 0))
  )


rich_summary <- rich_per_visit %>%
  group_by(Category, Visit) %>%
  summarise(
    mean_richness = mean(richness),
    sd = sd(richness),
    n = n(),
    se = sd / sqrt(n),
    .groups = "drop"
  )

ggplot(rich_per_visit, aes(x = factor(Visit), y = richness, fill = Category)) +
  geom_boxplot() +
  theme_bw() +
  labs(
    x = "Visit",
    y = "Species richness",
    title = "Species richness per visit by habitat"
  )


m_rich <- glmmTMB(richness ~ Category * Visit + (1 | Site), family = nbinom2,
  data = rich_per_visit)
summary(m_rich)
mod_dharma1 <- m_rich %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(m_rich) 
#check_collinearity(global.model) 
plot(allEffects(m_rich))

emm <- emmeans(m_rich, ~ Category | Visit)
pairs(emm)


ggplot(rich_per_visit, aes(x = Category, y = richness, fill = Category)) +
  geom_boxplot() +
  facet_wrap(~ Visit) +
  stat_compare_means(method = "kruskal.test") +
  theme_bw()

# plot abundances per species
comm_long <- community.visit %>%
  pivot_longer(
    cols = -(City:Visit),  # keep metadata columns
    names_to = "Species",
    values_to = "Abundance")

comm_long <- comm_long %>%
  filter(Abundance > 0)


comm_sum <- comm_long %>%
  group_by(Species, Category, Visit) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop")


ggplot(comm_sum, aes(x = Species, y = total_abundance, fill = Category)) +
  geom_col(position = "stack") +
  facet_wrap(~ Visit, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
  labs(title = "Species abundances per habitat across visits",
    x = "Species",
    y = "Total abundance")
  

common_species <- comm_sum %>%
  group_by(Species) %>%
  summarise(total = sum(total_abundance)) %>%
  filter(total > 10) %>%   # adjust threshold
  pull(Species)

comm_sum_filt <- comm_sum %>%
  mutate(Species = reorder(Species, total_abundance, sum)) %>% 
  filter(Species %in% common_species)

ggplot(comm_sum_filt, aes(x = Species, y = total_abundance, fill = Category)) +
  geom_col() +
  facet_wrap(~ Visit, scales = "free_y") +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "top") +
  labs(title="Species more than 10 occurrences only (49 species of 196)")

                                     

# calculate beta-div and plot
cat.species <- df.all[,5:ncol(df.all)]
occurrences <- df.all[,6:ncol(df.all)]
category <- df.all$Category
visit <- df.all$Visit


# Beta diversity: abundance based
beta.abund <- beta.pair.abund(occurrences, index.family = "bray")
# Treatment
bd.cat <-betadisper(beta.abund[[3]], category)
bd.cat <-betadisper(beta.abund[[3]], category)
plot(bd.cat, label=FALSE)
# Visit
bd.visit <-betadisper(beta.abund[[3]], visit)
plot(bd.visit)

# Compare species compositions across time
visit1 <- df.all %>% 
  filter(Visit == "1")
visit1.comm <- visit1[,6:ncol(visit1)]
#convert to presence-absence
pa.visit1 <- presabs <- ifelse(visit1.comm>0,1,0)


visit2 <- df.all %>% 
  filter(Visit == "2")
visit2.comm <- visit2[,6:ncol(visit2)]
#convert to presence-absence
pa.visit2 <- presabs <- ifelse(visit2.comm>0,1,0)

visit3 <- df.all %>% 
  filter(Visit == "3")
visit3.comm <- visit3[,6:ncol(visit3)]
#convert to presence-absence
pa.visit3 <- presabs <- ifelse(visit3.comm>0,1,0)

visit4 <- df.all %>% 
  filter(Visit == "4")
visit4.comm <- visit4[,6:ncol(visit4)]
#convert to presence-absence
pa.visit4 <- presabs <- ifelse(visit4.comm>0,1,0)


# Using all data:

beta.temporal1.2 <- beta.temp(pa.visit1, pa.visit2, index.family = "sorensen")
# beta.sim = turnover
# beta.sor = sorensen dissimilarity
# beta.sne = nestedness

#write.csv2(beta.temporal1.2, "beta.temp.v1.v2.csv", fileEncoding = "latin1")

beta.temporal2.3 <- beta.temp(pa.visit2, pa.visit3, index.family = "sorensen")
#write.csv2(beta.temporal2.3, "beta.temp.v2.v3.csv", fileEncoding = "latin1")

beta.temporal3.4 <- beta.temp(pa.visit3, pa.visit4, index.family = "sorensen")
#write.csv2(beta.temporal3.4, "beta.temp.v3.v4.csv", fileEncoding = "latin1")

beta.temporal1.3 <- beta.temp(pa.visit1, pa.visit3, index.family = "sorensen")
#write.csv2(beta.temporal1.3, "beta.temp.v1.v3.csv", fileEncoding = "latin1")

beta.temporal1.4 <- beta.temp(pa.visit1, pa.visit4, index.family = "sorensen")
#write.csv2(beta.temporal1.4, "beta.temp.v1.v4.csv", fileEncoding = "latin1")

beta.temporal2.4 <- beta.temp(pa.visit2, pa.visit4, index.family = "sorensen")
#write.csv2(beta.temporal2.4, "beta.temp.v2.v4.csv", fileEncoding = "latin1")


# plot and test:
beta.visit <- read.csv2(here("data/databases/temporal beta div", "beta.temporal.between.visits.csv"))


beta.div <- ggplot(beta.visit, aes(x=Category, y=value, fill=Beta.type)) + geom_boxplot() +
  facet_wrap(vars(Visit)) +
  ggtitle("Temporal beta diversity (incidence based: Sorensen)")
beta.div

# transform data and remove missing values
n <- nrow(beta.visit)
beta.visit$value2 <- (beta.visit$value * (n - 1) + 0.5) / n 

beta.visit <- beta.visit %>%
  filter(!is.na(value))


# test turnover
beta.visit <- beta.visit %>% 
  mutate(Category = as.factor(Category)) %>% 
  mutate(Visit = as.factor(Visit))

beta.turnover <- beta.visit %>% 
  filter(Beta.type %in% "beta.sim")

beta.turnover$value2 <- pmin(pmax(beta.turnover$value2, 0.001), 0.999) # because i get all kinds of issues with the diagnostics

turnover.mod <- glmmTMB(value2 ~ Category * Visit + (1 | Site), family = beta_family(), data = beta.turnover)

summary(turnover.mod)
mod_dharma1 <- turnover.mod %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(turnover.mod) 
#check_collinearity(global.model) 
plot(allEffects(turnover.mod))


# test sorensen dissimilarity
beta.diss <- beta.visit %>% 
  filter(Beta.type %in% "beta.sor")

diss.mod <- glmmTMB(value2 ~ Category * Visit + (1 | Site), family = beta_family(), data = beta.diss)
summary(diss.mod)
mod_dharma1 <- diss.mod %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(diss.mod) 
#check_collinearity(global.model) 
plot(allEffects(diss.mod))


# test nestedness
beta.nest <- beta.visit %>% 
  filter(Beta.type %in% "beta.sne")

nested.mod <- glmmTMB(value2 ~ Category * Visit + (1 | Site), family = beta_family(), data = beta.nest)
summary(nested.mod)
mod_dharma1 <- nested.mod %>% simulateResiduals(n=1000)
plot(mod_dharma1)
AICc(nested.mod) 
#check_collinearity(global.model) 
plot(allEffects(nested.mod))
