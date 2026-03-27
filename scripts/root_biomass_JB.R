#### packages ####
library(tidyverse)
library(ggplot2)
library(dunn.test)
library(dplyr)
library(data.table)
library(brms)

#### Data loading ####

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


root_morphology <- data.table(root_morphology)
root_morphology[,pft_biomass := as.numeric(pft_biomass)]
root_morphology <- root_morphology[!is.na(pft_biomass),]
root_morphology <- root_morphology[!duplicated(root_morphology),]

# root biomass ~ functional type by treatment graph
ggplot(data = root_morphology, aes(x = treatment, y = pft_biomass, fill = treatment)) +
  geom_boxplot(alpha = 0.8,outliers = F) +
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
dunn.test(graminoid$pft_biomass, graminoid$treatment, method = "bonferroni")

shrub <-  root_morphology %>%
  filter(functional_type_biomass == "shrub_biomass") 
dunn.test(shrub$pft_biomass, shrub$treatment, method = "bonferroni")

#### Fitting the models ####
library(brms)

#### setting the priors
## the mean of the distribution has to be negative so that the tail end is in the potive realm
prior_root_biomass <- set_prior("normal(-0.03,0.05)" , class = "Intercept")
## my prior belief is no effect of treatment on biomass, and that the maximum effect 2* standard deviation of the normal distribution
prior_root_biomass <- c(prior_root_biomass,set_prior("normal(0,0.03)" , class = "b" , coef  = "treatmentHeatwave"))
prior_root_biomass <- c(prior_root_biomass,set_prior("normal(0,0.03)" , class = "b" , coef  = "treatmentExtendedseason"))
## MY prior belief is that the the vaiability around the mean biomass is between 0 and 0.01 + 0.05*2 , it is 
## probably less variable in real life
prior_root_biomass <- c(prior_root_biomass,set_prior("normal(0.01,0.05)" , class = "sigma"))

## truncated normal distribution, could also be exponential distribution
hist(graminoid$pft_biomass)
hist(shrub$pft_biomass)

bayesian_model_gram <- brm( pft_biomass |  trunc(lb=0)   ~  treatment ,  # model formula, trunc is to force the distributino to be above 0
                       data = graminoid, # dataset
                       iter = 5000, # number of sampling iteration
                       warmup = 1000, # discarded iterations at the start
                       cores = 3, # a core compute a chain, 3 time faster
                       chains = 3, # number of independent models that we want to converge
                       prior = prior_root_biomass , # my prior knowledge
                       control = list(adapt_delta = 0.99), # increase the time for computation, but allows for a better estimation of the posterior
                       family = gaussian(), # distribution
                       #threads = threading(3), # even faster
                       init = 0) # more stable sampling
summary(bayesian_model_gram,prob = 0.9)
plot(bayesian_model_gram) ##convergence OK
pp_check(bayesian_model_gram) ## fitting the distributino somewhat ok

## prior for the shrubs,the values are bigger because we found heavier roots during the experiment 
prior_root_biomass <- set_prior("normal(-0.05,0.15)" , class = "Intercept")
prior_root_biomass <- c(prior_root_biomass,set_prior("normal(0,0.09)" , class = "b" , coef  = "treatmentHeatwave"))
prior_root_biomass <- c(prior_root_biomass,set_prior("normal(0,0.09)" , class = "b" , coef  = "treatmentExtendedseason"))
prior_root_biomass <- c(prior_root_biomass,set_prior("normal(0.02,0.05)" , class = "sigma"))

bayesian_model_shrub <- brm( pft_biomass |  trunc(lb=0)   ~  treatment , ## model formula
                       data = shrub, # dataset
                       iter = 5000, # number of sampling iteration
                       warmup = 1000, # discarded iterations at the start
                       cores = 3, # a core compute a chain, 3 time faster
                       chains = 3, # number of independent models that we want to converge
                       prior = prior_root_biomass ,
                       control = list(adapt_delta = 0.99),
                       family = gaussian(), # distribution
                       #threads = threading(3), # even faster
                       init = 0) # more stable sampling
s
summary(bayesian_model_shrub,prob = 0.9)
plot(bayesian_model_shrub) ##convergence OK
pp_check(bayesian_model_shrub) ## fitting the distribution somewhat ok

#### statistical reporting with the fit models ####

### this is for the shrub, repeat this process for the gram !
## I create a dummy data frame with the 3 treatments
new_data <-  data.table(treatment = sort(unique(shrub$treatment)))
## I gather the distribution the root biomass mean value (posterior) with the fitted function
## I can report the mean (Estimate) and the 90 % credible interval
preds <- cbind(new_data,fitted(bayesian_model_shrub,# change here to the gram model
                               newdata = new_data,probs = c(0.05,0.95)))

## I get all the draws (full distribution, summary = F), usefull to compute a new distributino I want to study
preds_full <- data.table(fitted(bayesian_model_shrub, # change here to the gram model
                                newdata = new_data,summary = F))
colnames(preds_full) <- as.character(new_data$treatment)

## the new distributino I want to study is the difference between a treatment and the control
preds_delta <- preds_full[,.(`Heat wave` = `Heat wave` - Control ,
                             `Extended season` =`Extended season` - Control)]

## steps to make ggplot happy
preds_full <- melt(preds_full,variable.name = "treatment",value.name = "Estimate")
preds_full$treatment <- as.factor(preds_full$treatment)

preds_delta <- melt(preds_delta,variable.name = "treatment",value.name = "Estimate")

## the plots
# distribution of the mean 
ggplot(preds,aes(x = treatment,y = Estimate))+
  geom_violin(data = preds_full,aes(fill = treatment))+
  geom_pointrange(aes(ymin =  Q5, ymax = Q95),color="white")+
  geom_point(data = shrub,aes(y = pft_biomass,x = treatment))+
  theme_classic()+  
  scale_fill_manual(values = c("Control" = "#6c6563",
                             "Heat wave" = "#b56d5e",
                              "Extended season" = "#bbbc81")) 

## how much of the delta control - treatment is positive ?
## if more than 90 %, I'm confident in my treatment effect
preds_delta[,1-mean(Estimate<0),by = treatment]
preds_delta[,mean(Estimate),by = treatment]

## distribution of the delta control - treatment
ggplot(preds_delta,aes(x = Estimate,fill = treatment ))+
  geom_vline(xintercept = 0,lty = 2)+
  geom_density(alpha = 0.5)+
  theme_classic()+  
  scale_fill_manual(values = c("Control" = "#6c6563",
                               "Heat wave" = "#b56d5e",
                               "Extended season" = "#bbbc81"))+
  facet_wrap(~treatment,scales = "free_y",nrow = 3)

## distribution of the delta control - treatment, I like this less
ggplot(preds_delta,aes(x = treatment, y = Estimate ,fill = treatment ))+
  geom_hline(yintercept = 0,lty = 2)+
  geom_violin(alpha = 0.95)+
  theme_classic()+  
  scale_fill_manual(values = c("Control" = "#6c6563",
                               "Heat wave" = "#b56d5e",
                               "Extended season" = "#bbbc81"))



# https://distribution-explorer.github.io/

