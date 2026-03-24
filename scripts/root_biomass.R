library(tidyverse)
library(ggplot2)
library(dunn.test)
library(dplyr)

root_data <- read.csv(file = "data_raw/branch_root_data.csv") %>%
  mutate(treatment = factor(treatment,
                            levels = c("control", "heatwave", "extended"),
                            labels = c("Control", "Heat wave", "Extended season")))
      str(root_data)
      root_data$thin_white_branched <- as.numeric(root_data$thin_white_branched)
      root_data$thin_white_unbranched <- as.numeric(root_data$thin_white_unbranched)
      root_data$thick_white_unbranched <- as.numeric(root_data$thin_white_unbranched)
      root_data$thin_beige_branched <- as.numeric(root_data$thin_beige_branched)
      root_data$thin_beige_unbranched <- as.numeric(root_data$thin_beige_unbranched)
      root_data$thick_beige_unbranched <- as.numeric(root_data$thick_beige_unbranched)
      root_data$thin_brown_branched <- as.numeric(root_data$thin_brown_branched )
      root_data$thin_brown_unbranched <- as.numeric(root_data$thin_brown_unbranched )
      root_data$thick_brown_unbranched <- as.numeric(root_data$thick_brown_unbranched )
      root_data$thin_black_branched  <- as.numeric(root_data$thin_black_branched)
      root_data$thin_black_unbranched  <- as.numeric(root_data$thin_black_unbranched)
      root_data$thick_black_unbranched  <- as.numeric(root_data$thick_black_unbranched)

root_morphology <- root_data %>%
  pivot_longer(cols = shrub_biomass:graminoid_biomass,
               names_to = "functional_type_biomass",
               values_to = "pft_biomass") %>%
  filter(pft_biomass != "NA")

# root biomass ~ functional type by treatment graph
ggplot(data = root_morphology, aes(x = treatment, y = pft_biomass, fill = treatment)) +
  geom_boxplot(alpha = 0.8) +
  facet_grid(col = vars(functional_type_biomass),
             labeller = labeller(functional_type_biomass = c("graminoid_biomass" = "Graminoid",
                                                             "shrub_biomass" = "Shrub"))) + 
  labs(x = "Treatment",
       y = "Root biomass (g)") +
  geom_jitter(width = 0.12, alpha = 0.4, size = 1) +
  scale_fill_manual(values = c("Control" = "#6c6563",
                               "Heat wave" = "#b56d5e",
                               "Extended season" = "#bbbc81")) +
  theme_classic() +
  theme(legend.position = c(0.12, 0.87),
          legend.title = element_blank(),
          legend.text = element_text(size = 9), 
        strip.background = element_rect(colour = "black", fill = "#EDEDED"),
          strip.text.x = element_text(size = 12, colour = "black"), 
        axis.text = element_text(size = 9),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 13, margin = margin(t = 12)),
        axis.title.y = element_text(size = 13, margin = margin(r = 12)),
        panel.spacing = unit(1, "lines"),
          panel.border = element_rect(color = "black", fill = NA, size = 1))

# root biomass ~ functional type by treatment analysis
root_morphology %>%
  group_by(functional_type_biomass) %>%
  summarise(p_value = kruskal.test(pft_biomass ~ treatment)$p.value)

graminoid <-  root_morphology %>%
  filter(functional_type_biomass == "graminoid_biomass") 
dunn.test(root_morphology$pft_biomass, root_morphology$treatment, method = "bonferroni")

shrub <-  root_morphology %>%
  filter(functional_type_biomass == "shrub_biomass") 
dunn.test(root_morphology$pft_biomass, root_morphology$treatment, method = "bonferroni")

