library(tidyverse)
library(ggplot2)
library(dplyr)
library(lmerTest)
library(lme4)

root_data <- read.csv(file = "data_raw/branch_root_data.csv") %>%
  mutate(treatment = factor(treatment,
                            levels = c("control", "heatwave", "extended"),
                            labels = c("Control", "Heat wave", "Extended season")))

    str(root_data)
    root_data$thin_white_branched <- as.numeric(root_data$thin_white_branched)
    root_data$thin_white_unbranched <- as.numeric(root_data$thin_white_unbranched)
    root_data$thick_white_unbranched <- as.numeric(root_data$thick_white_unbranched)
    root_data$thin_beige_branched <- as.numeric(root_data$thin_beige_branched)
    root_data$thin_beige_unbranched <- as.numeric(root_data$thin_beige_unbranched)
    root_data$thick_beige_unbranched <- as.numeric(root_data$thick_beige_unbranched)
    root_data$thin_brown_branched <- as.numeric(root_data$thin_brown_branched )
    root_data$thin_brown_unbranched <- as.numeric(root_data$thin_brown_unbranched )
    root_data$thick_brown_unbranched <- as.numeric(root_data$thick_brown_unbranched )
    root_data$thin_black_branched  <- as.numeric(root_data$thin_black_branched)
    root_data$thin_black_unbranched  <- as.numeric(root_data$thin_black_unbranched)
    root_data$thick_black_unbranched  <- as.numeric(root_data$thick_black_unbranched)
    
co2_flux <- read_csv(file = "data_raw/shifted_co2_data.csv")

co2_flux_shifted <- co2_flux %>%
  mutate(measurement.week = ifelse(Treatment == "extended",
                                   measurement.week - 3,
                                   measurement.week),
         Treatment = factor(Treatment,
                            levels = c("Control", "Heatwave", "Extended"),
                            labels = c("Control", "Heat wave", "Extended season")))


# calculated co2 flux for each plant
avg_season_flux <- co2_flux_shifted %>%
  group_by(plant_id, Treatment) %>%
  summarise(total_co2_flux = sum(seasonal_co2_flux, na.rm = TRUE))


# combining co2 flux with root data frames
flux_root_combine <- left_join(avg_season_flux, root_data, by = 'plant_id') %>%
  filter(Treatment != "NA",
         total_root_biomass != "NA")


# co2 flux ~ root biomass by treatment graph
ggplot(data = flux_root_combine, aes(x = total_root_biomass, y = (total_co2_flux*86400), color = Treatment, fill = Treatment)) +
  geom_smooth(method = lm, alpha = 0.2, size = 1) +
  geom_point(alpha = 0.4) +
  labs(x = "Root biomass (g)",
       y = "Cumulative CO₂ respiration (µmol m⁻²)") +
  scale_color_manual(values = c("Control" = "#6c6563",
                               "Heat wave" = "#b56d5e",
                               "Extended season" = "#bbbc81")) +
  scale_fill_manual(values = c("Control" = "#6c6563",
                               "Heat wave" = "#b56d5e",
                               "Extended season" = "#bbbc81")) +
  theme_classic() +
  theme(legend.position = c(0.13, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 9),
        axis.title.x = element_text(size = 13, margin = margin(t = 12)),
        axis.title.y = element_text(size = 13, margin = margin(r = 12))) +
  annotate("text",
           x = 0.025,
           y = 0.8e7,
           label = "Treatment: p < 0.001, Biomass: p < 0.01",
           hjust = 0)


# co2 flux ~ root biomass by treatment analysis
flux_root_combine_anova <- lm(total_co2_flux ~ treatment*total_root_biomass, data = flux_root_combine)
anova(flux_root_combine_anova)
plot(flux_root_combine_anova)
