library(tidyverse)
library(ggplot2)
library(dunn.test)
library(dplyr)
library(data.table)
library(brms)
library(patchwork)

R.version.string
Sys.setenv(PATH = paste0("C:\\rtools43\\usr\\bin;", Sys.getenv("C:\rtools43")))
options(buildtools.check = function(action) TRUE )


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

root_morphology <- root_data %>%
  pivot_longer(cols = shrub_biomass:graminoid_biomass,
               names_to = "functional_type_biomass",
               values_to = "pft_biomass") %>%
  filter(pft_biomass != "NA")

# fixing duplicates
root_morphology <- data.table(root_morphology)
root_morphology[,pft_biomass := as.numeric(pft_biomass)]
root_morphology <- root_morphology[!is.na(pft_biomass),]
root_morphology <- root_morphology[!duplicated(root_morphology),]


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


# GRAMINOID bayesian analysis 
prior_root_biomass <- set_prior("normal(-0.03,0.05)" , class = "Intercept")
prior_root_biomass <- c(prior_root_biomass,set_prior("normal(0,0.03)" , class = "b" , coef  = "treatmentHeatwave"))
prior_root_biomass <- c(prior_root_biomass,set_prior("normal(0,0.03)" , class = "b" , coef  = "treatmentExtendedseason"))
prior_root_biomass <- c(prior_root_biomass,set_prior("normal(0.01,0.5)" , class = "sigma"))

bayesian_model_gram <- brm( pft_biomass |  trunc(lb = 0) ~   treatment , # |  trunc(lb=0) 
                       data = graminoid,
                       iter = 5000,
                       warmup = 1000,
                       cores = 3,
                       chains = 3,
                       prior = prior_root_biomass ,
                       control = list(adapt_delta = 0.99),
                       family = gaussian(),
                       # threads = threading(3),
                       init = 0)

summary(bayesian_model_gram, prob = 0.9)
plot(bayesian_model_gram)
pp_check(bayesian_model_gram)
  # no difference between treatments for graminoid roots



# SHRUB bayesian analysis
prior_root_biomass_shrub <- set_prior("normal(-0.05,0.15)" , class = "Intercept")
prior_root_biomass_shrub <- c(prior_root_biomass_shrub,set_prior("normal(0,0.09)" , class = "b" , coef  = "treatmentHeatwave"))
prior_root_biomass_shrub <- c(prior_root_biomass_shrub,set_prior("normal(0,0.09)" , class = "b" , coef  = "treatmentExtendedseason"))
prior_root_biomass_shrub <- c(prior_root_biomass_shrub,set_prior("normal(0.02,0.05)" , class = "sigma"))

bayesian_model_shrub <- brm( pft_biomass |  trunc(lb=0)   ~  treatment , # |  trunc(lb=0)  # model formula
                       data = shrub, # dataset
                       iter = 5000, # number of smapling iteration
                       warmup = 1000, # discarded iterations at the start
                       cores = 3, # a core compute a chain, 3 time faster
                       chains = 3, # number of independant models that we want to converge
                       prior = prior_root_biomass_shrub ,
                       control = list(adapt_delta = 0.99),
                       family = gaussian(), # distribution
                       #threads = threading(3), # even faster
                       init = 0) # more stable sampling

summary(bayesian_model_shrub, prob = 0.9)
plot(bayesian_model_shrub)
pp_check(bayesian_model_shrub)

### SHRUB statistical reporting with the fit models ###
new_data_shrub <- data.table(treatment = sort(unique(shrub$treatment)))
preds_shrub <- cbind(new_data_shrub, fitted(bayesian_model_shrub, newdata = new_data_shrub, probs = c(0.05, 0.95)))

preds_full_shrub <- data.table(fitted(bayesian_model_shrub,
                                      newdata = new_data_shrub,
                                      summary = FALSE))

colnames(preds_full_shrub) <- as.character(new_data_shrub$treatment)

preds_delta_shrub <- preds_full_shrub[, .(`Heat wave` = `Heat wave` - Control, `Extended season` = `Extended season` - Control)]

# plot set up
preds_full_shrub <- melt(preds_full_shrub, measure.vars = names(preds_full_shrub), variable.name = "treatment", value.name = "Estimate")

preds_delta_shrub <- melt(preds_delta_shrub, measure.vars = names(preds_delta_shrub), variable.name = "treatment", value.name = "Estimate")

preds_full_shrub$treatment <- as.factor(preds_full_shrub$treatment)

# distribution of the mean plot
plot_shrub1 <- ggplot(preds_shrub, aes(x = treatment, y = Estimate)) +
  geom_violin(data = preds_full_shrub, trim = FALSE, alpha = 0.8, aes(fill = treatment)) +
  geom_pointrange(aes(ymin = Q5, ymax = Q95), color = "white", alpha = 0.7) +
  geom_point(data = shrub, aes(y = pft_biomass, x = treatment), alpha = 0.6, size = 2) +
  theme_classic() +
  scale_fill_manual(values = c("Control" = "#6c6563",
                               "Heat wave" = "#b56d5e",
                               "Extended season" = "#bbbc81")) +
  labs(x = "Treatment",
       y = "Estimated root biomass") +
  theme(legend.position = "none",
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 15, margin = margin(t = 5)),
        axis.title.y = element_text(size = 15, margin = margin(r = 1)))
plot_shrub1
# how much of the delta control - treatment is positive? 
# if more than 90% i'm confident in treatment effect
preds_delta[,1-mean(Estimate < 0), by = treatment]
preds_delta[,mean(Estimate), by = treatment]

preds_delta$treatment <- factor(preds_delta$treatment,
                                levels = c("Heat wave", "Extended season", "Control"))

plot_shrub2 <- ggplot(preds_delta, aes(x = Estimate, fill = treatment)) + 
                geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "black") +
                geom_density(alpha = 0.5) +
                theme_classic()+  
                scale_fill_manual(values = c("Control" = "#6c6563",
                                             "Heat wave" = "#b56d5e",
                                             "Extended season" = "#bbbc81")) +
                labs(x = "Estimated difference in root biomass",
                     y = "Probability density") +
                facet_wrap(~treatment,scales = "free_y",nrow = 3) +
                theme(legend.position = "none",
                      legend.title = element_blank(),
                      legend.text = element_text(size = 10),
                      strip.text = element_blank(),
                      strip.background = element_blank(), 
                      panel.spacing = unit(1, "lines"),
                      axis.text = element_text(size = 10),
                      axis.title.x = element_text(size = 15, margin = margin(t = 5)),
                      axis.title.y = element_text(size = 15, margin = margin(r = 1))) +
                geom_segment(data = subset(preds_delta, treatment == "Extended season"),
                             aes(x = 0.03, y = 25, xend = 0.15, yend = 25),
                             arrow = arrow(type = "open", length = unit(0.2, "inches")),
                             color = "#6c6563",
                             linewidth = 2.5) +
                geom_text(data = subset(preds_delta, treatment == "Heat wave"),
                          aes(x = 0.02, y = 70, label = "33% chance shrub root\nbiomass increases"),
                          size = 5.5,
                          hjust = 0,
                          fontface = "plain") +
                geom_text(data = subset(preds_delta, treatment == "Extended season"),
                          aes(x = 0.03, y = 40, label = "93% chance shrub root\nbiomass increases"),
                          size = 5.5,
                          hjust = 0,
                          fontface = "plain") 
plot_shrub2



### GRAMINOID statistical reporting with the fit models ###
new_data_gram <- data.table(treatment = sort(unique(graminoid$treatment)))
preds_gram <- cbind(new_data_gram, fitted(bayesian_model_gram,
                                  newdata = new_data_gram, probs = c(0.05, 0.95)))

preds_full_gram <- data.table(fitted(bayesian_model_gram,
                                newdata = new_data_gram, summary = F))
colnames(preds_full_gram) <- as.character(new_data_gram$treatment)

preds_delta_gram <- preds_full_gram[,.(`Heat wave` = `Heat wave` - Control ,
                             `Extended season` =`Extended season` - Control)]

# plot set up
preds_full_gram <- melt(preds_full_gram, measure.vars = names(preds_full_gram), variable.name = "treatment", value.name = "Estimate")
preds_full_gram$treatment <- as.factor(preds_full_gram$treatment)

preds_delta_gram <- melt(preds_delta_gram, measure.vars = names(preds_delta_gram), variable.name = "treatment", value.name = "Estimate")

# distribution of the mean plot
plot_gram1 <- ggplot(preds_gram, aes(x = treatment, y = Estimate)) +
              geom_violin(data = preds_full_gram, trim = FALSE, alpha = 0.8, aes(fill = treatment)) +
              geom_pointrange(aes(ymin = Q5, ymax = Q95), color = "white", alpha = 0.7) + 
              geom_point(data = graminoid, aes(y = pft_biomass, x = treatment), alpha = 0.6, size = 2) +
              theme_classic() +
              scale_fill_manual(values = c("Control" = "#6c6563",
                                           "Heat wave" = "#b56d5e",
                                           "Extended season" = "#bbbc81")) +
              labs(x = "Treatment",
                   y = "Estimated root biomass") +
              theme(legend.position = "none",
                    legend.title = element_blank(),
                    legend.text = element_text(size = 10),
                    strip.text = element_blank(),
                    panel.spacing = unit(1, "lines"),
                    axis.text = element_text(size = 10),
                    axis.title.x = element_text(size = 15, margin = margin(t = 5)),
                    axis.title.y = element_text(size = 15, margin = margin(r = 1))) 
plot_gram1

# how much of the delta control - treatment is positive? 
# if more than 90% i'm confident in treatment effect
preds_delta_gram[,1-mean(Estimate < 0), by = treatment]
preds_delta_gram[,mean(Estimate), by = treatment]

plot_gram2 <- ggplot(preds_delta_gram, aes(x = Estimate, fill = treatment)) + 
              geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "black") +
              geom_density(alpha = 0.5) +
              theme_classic()+  
              scale_fill_manual(values = c("Control" = "#6c6563",
                                           "Heat wave" = "#b56d5e",
                                           "Extended season" = "#bbbc81"))+
              labs(x = "Estimated difference in root biomass",
                   y = "Probability density") +
              facet_wrap(~treatment, scales = "free_y", nrow = 3, strip.position = "right") +
              theme(legend.position = "none",
                    legend.title = element_blank(),
                    legend.text = element_text(size = 10),
                    strip.text = element_text(angle = 90, size = 12),
                      strip.background = element_rect(fill = "white", color = "black"),
                    panel.spacing = unit(1, "lines"),
                    axis.text = element_text(size = 10),
                    axis.title.x = element_text(size = 15, margin = margin(t = 5)),
                    axis.title.y = element_text(size = 15, margin = margin(r = 1))) +
              geom_text(data = subset(preds_delta, treatment == "Heat wave"),
                        aes(x = 0.008, y = 200, label = "55% chance graminoid\nroot biomass increases"),
                        size = 5.5,
                        hjust = 0,
                        fontface = "plain") +
              geom_text(data = subset(preds_delta, treatment == "Extended season"),
                        aes(x = 0.008, y = 225, label = "48% chance graminoid\nroot biomass increases"),
                        size = 5.5,
                        hjust = 0,
                        fontface = "plain") 
              
  
plot_gram2

### combining plots
plot_shrub1 <- plot_shrub1 + ggtitle(NULL)
plot_gram1 <- plot_gram1 + ggtitle(NULL)

(plot_shrub1 | plot_gram1) +
  plot_layout(axis_titles = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "none",
        plot.tag = element_text(face = "bold", size = 20))


plot_shrub2 <- plot_shrub2 + ggtitle(NULL)
plot_gram2 <- plot_gram2 + ggtitle(NULL)

(plot_shrub2 | plot_gram2) +
  plot_layout(axis_titles = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "none",
        plot.tag = element_text(face = "bold", size = 20))

