# Install/Load packages ------------------------------------
packages = c("adespatial", "tidyverse", "vegan", "readxl", "ggtext", 
             "patchwork", "cowplot", "ggrepel", "FD", "ade4", "ape",
             "picante", "betareg", "lmtest", "rcompanion", "multcompView", "emmeans")
#-----------------------------------------------------------
lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
    } 
    library(x, character.only = TRUE)
  }
)
# Importing data --------------------------------------------
data <- read_excel("raw_data.xlsx", sheet="abundance")
abund <- data %>%
  select(Micfur:Rhilal)
species <- read_excel("raw_data.xlsx", sheet = "species")
env <- read_excel("raw_data.xlsx", sheet = "environmental") %>% 
  bind_cols(Season = as.factor(data$Season), 
            transect = factor(data$transect, 
                              levels = c("2m","5m", "7m", "10m", "15m", "16m")))

# Measuring Beta diversity contributions ---------------------
BDt <- beta.div(abund, method="hellinger",sqrt.D=FALSE,samp=TRUE,
                nperm=2,save.D=FALSE)
env$LCBD <- BDt$LCBD

# Species contribution to beta diversity ---------------------
# Data frame set up
species$scbd <- BDt$SCBD
species$species <- as.factor(species$species)
species$guilds <- factor(species$guilds, ordered = TRUE, 
                         levels = c('Demersal','Benthopelagic', 'Reef-associated', 
                                    'Pelagic-neritic', 'Pelagic-oceanic'))

# SCBD values summary
SCBD_summary = species %>% 
  group_by(guilds) %>% 
  summarise(scdb_mean = mean(scbd),
            scdb_sd = sd(scbd),
            scdb_sum = sum(scbd)) 

# Boxplot 
SCBD_bxp <- ggplot(species) +
  geom_boxplot(aes(y = guilds, x = scbd, color = guilds, fill = guilds, alpha = 0.3),
               coef = 30, width = 0.3) +
  geom_jitter(aes(y = as.numeric(guilds)-0.35, x = scbd, color = guilds, fill = guilds, alpha = 0.3),
              height = 0.1, shape = 21, color = 'black', size = 3) +
  scale_x_continuous(breaks = seq(0, 0.12, by = 0.01)) +
  xlab("SCBD") +
  ylab(NULL) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")

SCBD_bxp
#ggsave('SCBD_bxp.png', SCBD_bxp, units = 'cm', width = 30, height = 25)

# Barplot
tx = 20
SCBD_barp <- ggplot(species, aes(x = reorder(species, scbd), y = scbd, fill = guilds)) + 
  geom_bar(stat = "identity", alpha=0.6) +
  ylab("SCBD") + 
  xlab(NULL) +
  theme(legend.position="bottom") + 
  guides(fill = guide_legend(title = "")) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,0.12),
                     breaks = seq(0, 0.12, by = 0.01)) +
  theme_classic() +
  coord_flip(clip = "off") +
  theme(axis.text.y = element_text(face = "italic"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size = tx),
        legend.text = element_text(size = tx-5),
        axis.text.x = element_text(size = tx-8),
        axis.title.x = element_text(size = tx),
        legend.position="bottom")
SCBD_barp

SCBD_plt = SCBD_barp + inset_element(SCBD_bxp,
                                     left = 0.05, bottom = 0.01,
                                     right = 0.97, top = 0.6)
SCBD_plt

ggsave('SCBD.png', SCBD_plt, units = 'cm', width = 30, height = 35)

# Kruskal-Wallis test
kruskal.test(scbd ~ guilds, data = species)

# Beta Regression Predictors Set Up
species$relative_abundance <- species$total_abundance / sum(species$total_abundance)
Y <- abund
Y[Y>0] <- 1
species$n_sites <- colSums(Y)

# Beta Regression Results
mod_scbd <- betareg(scbd ~ n_sites + relative_abundance ,
                    data = species)
car::vif(mod_scbd)
summary(mod_scbd)

# Beta Regression Plots
Sites_plt <- ggplot(data = species) + 
  labs(x = "Occurrence (Number of Sites)", 
       y = "SCBD") +
geom_point(aes(x = n_sites, y = scbd),
             size = 4, shape = 21, alpha = 0.6,
           fill = "darkgreen") +
geom_smooth(aes(x = n_sites, y = scbd),
              method = 'lm', se = TRUE, color = "black") +
  scale_x_continuous(breaks = seq(0, 25, by = 5)) + 
  theme_classic(base_size = 15)+
  ylim(-0.005, 0.18)

Sites_plt

ggsave('Sites_plt.png', Sites_plt, units = 'cm', width = 30, height = 25)

Rabundance_plt <- ggplot(data = species) + 
  labs(x = "Relative Total Abundance", 
       y = "SCBD") +
  geom_point(aes(x = relative_abundance, y = scbd),
             size = 4, shape = 21, alpha = 0.6,
             fill = "darkgreen") +
  geom_smooth(aes(x = relative_abundance, y = scbd),
              method = 'lm', se = TRUE, color = "black") +
  theme_classic(base_size = 15)+
  ylim(-0.005, 0.18)

Rabundance_plt

ggsave('Rabundance_plt.png', Rabundance_plt, units = 'cm', width = 30, height = 25)

SCBD_mod <- plot_grid(Rabundance_plt, Sites_plt, ncol=1, nrow= 2,
                      labels = "AUTO", label_size = 12, hjust = -4.8)
SCBD_mod

ggsave('SCBD_mod.png', SCBD_mod, units = 'cm', width = 20, height = 25)

# Local contribution to beta diversity --------------------------
# Beta Regression Predictors Set Up
env$Rich <- specnumber(abund)
env$abundance <- apply(abund, 1, sum)
env$relative_abundance <- env$abundance / (sum(env$abundance))


species$common_length <- as.numeric(species$common_length)
trait_categ <- cbind.data.frame(
  guilds = species$guilds,
  genus = species$genus,
  family = species$family,
  order = species$order,
  class = species$class
)
trait_num <- cbind.data.frame(
  common_length = species$common_length,
  relative_abundance = species$relative_abundance
)
rownames(species) <- species$initials
rownames(trait_categ) <- rownames(species)
rownames(trait_num) <- rownames(species)
ktab_list <- ktab.list.df(list(trait_categ, trait_num))
dist_mist <- dist.ktab(ktab_list, type = c("N", "Q"))
dist_mist <- as.matrix(dist_mist)
env$fric <- dbFD(dist_mist, abund)$FRic

# Beta Regression Results
mod1 <- betareg(LCBD ~ fric + Rich + relative_abundance,
                data = env)
car::vif(mod1)
summary(mod1)

mod2 <- betareg(LCBD ~ BT + BS + Phi, data = env)
car::vif(mod2)
summary(mod2)

env$transect <- factor(env$transect, levels = c("2m", "5m", "7m",
                                                "10m", "15m", "16m"),
                       labels = c("2 m", "5 m", "7 m", "10 m",
                                  "15 m", "16 m"))

# Plots
LCBD_Phi <- ggplot(data = env) + 
  labs(x = "Mean sedimenter diameter (phi)", 
       y = "LCBD") +
  geom_point(aes(x = Phi, y = LCBD, fill = factor(transect)),
             size = 4, shape = 21, alpha = 0.7) +
  geom_smooth(aes(x = Phi, y = LCBD),
              method = 'lm', se = TRUE, color = "black") +
  scale_x_continuous(breaks = seq(0, 6, by = .2)) +
  scale_y_continuous(breaks = seq(0, 0.1, by = .01)) +
  guides(fill = guide_legend(title = "Depth (m)")) +
  theme_classic(base_size = 15) +
  theme(legend.position="none")

LCBD_BT <- ggplot(data = env) + 
  labs(x = "Bottom temperature", 
       y = "LCBD") +
  geom_point(aes(x = BT, y = LCBD, fill = factor(transect)),
             size = 4, shape = 21, alpha = 0.7) +
  geom_smooth(aes(x = BT, y = LCBD),
              method = 'lm', se = TRUE, color = "black") +
  scale_x_continuous(breaks = seq(20, 30, by = 1)) +
  scale_y_continuous(breaks = seq(0, 0.1, by = .01)) +
  guides(fill = guide_legend(title = "Depth (m)")) +
  theme_classic(base_size = 15) +
  theme(legend.position="none")

LCBD_BS <-  ggplot(data = env ) + 
  labs(x = "Bottom salinity", 
       y = "LCBD") +
  geom_point(aes(x = BS, y = LCBD, fill = factor(transect)),
             size = 4, shape = 21, alpha = 0.7) +
  geom_smooth(aes(x = BS, y = LCBD),
              method = 'lm', se = TRUE, color = "black") +
  scale_x_continuous(breaks = seq(29, 40, by = 1)) +
  scale_y_continuous(breaks = seq(0, 0.1, by = .01)) +
  guides(fill = guide_legend(title = "Depth")) +
  theme_classic(base_size = 15) +
  theme(legend.position="bottom")

LCBD <- plot_grid(LCBD_Phi, LCBD_BT, LCBD_BS, ncol=1, nrow=3, 
                  labels = "AUTO", label_size = 12, hjust = -4.5)
LCBD

ggsave('LCBD.png', LCBD, units = 'cm', width = 15, height = 40)

LCBD_Fric <- ggplot(data = env) + 
  labs(x = "Functional Richness", 
       y = "LCBD") +
  geom_point(aes(x = fric, y = LCBD, fill = factor(transect)),
             size = 4, shape = 21, alpha = 0.7) +
  geom_smooth(aes(x = fric, y = LCBD),
              method = 'lm', se = TRUE, color = "black") +
scale_y_continuous(breaks = seq(0, 0.1, by = .01))+
  guides(fill = guide_legend(title = "Depth (m)")) +
  theme_classic(base_size = 15) +
  theme(legend.position="none")

LCBD_Rich <- ggplot(data = env) + 
  labs(x = "Richness", 
       y = "LCBD") +
  geom_point(aes(x = Rich, y = LCBD, fill = factor(transect)),
             size = 4, shape = 21, alpha = 0.7) +
  geom_smooth(aes(x = Rich, y = LCBD),
              method = 'lm', se = TRUE, color = "black") +
  scale_y_continuous(breaks = seq(0, 0.1, by = .01))+
  scale_x_continuous(breaks = seq(0,40, by=10))+
  guides(fill = guide_legend(title = "Depth (m)"))+
  theme_classic(base_size = 15)+
  theme(legend.position="bottom")

LCBD_abundance <- ggplot(data = env) + 
  labs(x = "Relative Abundance", 
       y = "LCBD") +
  geom_point(aes(x = relative_abundance, y = LCBD, fill = factor(transect)),
             size = 4, shape = 21, alpha = 0.7) +
  geom_smooth(aes(x = relative_abundance, y = LCBD),
              method = 'lm', se = TRUE, color = "black") +
  scale_y_continuous(breaks = seq(0, 0.1, by = .01))+
  xlim(0,0.25)+
  guides(fill = guide_legend(title = "Depth (m)"))+
  theme_classic(base_size = 15)+
  theme(legend.position="none") 

LCBD_Community <- plot_grid(LCBD_abundance, LCBD_Fric, LCBD_Rich, ncol=1, nrow=3, 
                  labels = "AUTO", label_size = 12, hjust = -4.5)
LCBD_Community

ggsave('LCBD_Community.png', LCBD_Community, units = 'cm', width = 15, height = 30)

# Redundancy Analysis -----------------------------------------
sc = 3
rda_sp = rda(decostand(abund, method = 'hellinger') ~ BT + BS + Phi, data = env)
rda_sc = scores(rda_sp, scaling = sc)

biplot <- ggplot() +
  geom_segment(data = data.frame(rda_sc$biplot),
               mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(4, "mm")),
               color = 'red') +
  geom_text(data = data.frame(rda_sc$biplot),
            aes(x = RDA1 * 1.05, y = RDA2 * 1.05,
                label = rownames(data.frame(rda_sc$biplot))),
            color = "red") +
  labs(x = 'RDA 1', y = 'RDA 2') +
  theme_bw(base_size = 15)

biplot

# RDA species
rda_sp <- biplot +
  geom_point(data = data.frame(rda_sc$species, species),
             mapping = aes(x = RDA1, y = RDA2, fill = guilds, size = scbd),
             shape = 22, alpha = 0.7) +
  geom_label_repel(data = data.frame(rda_sc$species, species) %>%
                     mutate(should_be_labeled=ifelse(scbd > 0.033,TRUE, FALSE)) %>%
                     filter(should_be_labeled==TRUE),
                   aes(x = RDA1, y = RDA2, label=species, fill=guilds, alpha = 0.9),
                   show.legend = FALSE, size = 3,
                   box.padding = unit(0.5, "lines")) +
  guides(size = FALSE, alpha = FALSE,
         fill = guide_legend(title = "Groups", 
                             override.aes = list(size = 3))) + 
  labs(x = "RDA1 (18.01%)", y = "RDA2 (0.05%)")

rda_sp  

# RDA sites
rda_sites <- biplot +
  geom_point(data = data.frame(rda_sc$sites, env),
             mapping = aes(x = RDA1, y = RDA2, fill = factor(transect), size = LCBD),
             shape = 21, alpha = 0.7) +
  guides(size = FALSE,
         fill = guide_legend(title = "Depth",
                             override.aes = list(size = 3)))+ 
  labs(x = "RDA1 (18.01%)", y = "RDA2 (0.05%)")


rda_plt <- rda_sp / rda_sites
rda_plt

ggsave('RDA.png', rda_plt, units = 'cm', width = 22, height = 30)
