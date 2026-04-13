# Author:  Giacomo Bignardi
# Project: Letter to James et al. — Shared and unique taste in acoustic preferences
# Purpose: Fit a GLMM on repeated trials and decompose variance into shared
#          (stimulus-driven) and unique (individual) aesthetic taste components

# Load necessary packages (auto-installs if missing) ####
required_packages <- c("dplyr", "lme4", "readr")
installed <- rownames(installed.packages())
to_install <- required_packages[!required_packages %in% installed]
if (length(to_install)) install.packages(to_install)
lapply(required_packages, library, character.only = TRUE)

# Set seed
set.seed(42)

# Load data ####
# data obtained from James, L. S., Woolley, S. C., Sakata, J. T., Hilton, C. B., Ryan, M. J., & Mehr, S. A. (2026). Humans share acoustic preferences with other animals. Science, 391(6791), 1246-1249.
# https://github.com/themusiclab/animal-sounds
url <- "https://raw.githubusercontent.com/themusiclab/animal-sounds/main/data/animal-sounds-data.csv"
animal_sounds <- read_csv(url)

# Load function ####
source("vca.fun.R")

# Tidy data ####
# Ensure grouping variables are treated as factors
animal_sounds <- animal_sounds %>%
  mutate(user_id = as.factor(user_id),
         StimID  = as.factor(StimID))

# Name of species
animal_sounds %>%
    group_by(Category) %>%
    summarise(Species = paste(unique(Species), collapse = ", "))

# Tetain only stimuli seen more than once per person
# (i.e., repeated-exposure trials, needed to estimate within-person consistency)
# Use exact same sampling strategy as James et al. for (stimulus) intra-rater reliability 
repeated_trials <- animal_sounds %>%
  mutate(person_stim_id = paste(user_id, StimID)) %>%
  group_by(person_stim_id) %>%
   filter(n() > 1) %>%
  ungroup()

# Number of observation per individuals
repeated_trials %>%
  group_by(user_id) %>%
  summarise(n_stim = n_distinct(StimID)) %>% 
  reframe(
    min(n_stim),
    max(n_stim),
    median(n_stim)
  )
  
repeated_trials %>%
  group_by(user_id) %>%
  summarise(n_stim = n_distinct(StimID)) %>% 
  filter(n_stim > 1) %>% 
  distinct(user_id) %>% 
  nrow()

# Number of observation per trait
repeated_trials %>%
  filter(!is.na(correct)) %>%
  count(StimID) %>% 
  reframe(
    min(n),
    max(n),
    median(n)
  )

# Fit intercept-only GLMM ####
# Random effects:
#   (1 | user_id)          — individual differences in overall agreement
#   (1 | StimID)           — stimulus differences in average agreement (shared)
#   (1 | user_id:StimID)   — person-by-stimulus interaction (unique)
glmm_fit <- glmer(
  correct ~ 1 + (1 | user_id) + (1 | StimID) + (1 | user_id:StimID),
  family = binomial("logit"),
  data   = repeated_trials
)
summary(glmm_fit)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: binomial  ( logit )
# Formula: correct ~ 1 + (1 | user_id) + (1 | StimID) + (1 | user_id:StimID)
# Data: repeated_trials
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# 8380.5    8407.6   -4186.3    8372.5      6415 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -1.8674 -0.7743  0.4482  0.7102  1.5746 
# 
# Random effects:
#   Groups         Name        Variance Std.Dev.
# user_id:StimID (Intercept) 0.89366  0.9453  
# user_id        (Intercept) 0.09357  0.3059  
# StimID         (Intercept) 0.38860  0.6234  
# Number of obs: 6419, groups:  user_id:StimID, 3074; user_id, 1812; StimID, 110
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  0.15484    0.07788   1.988   0.0468 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# # Obtain 95% CI ####
# NOTE THIS MAY TAKE A MIN OR SO
ci_profile <- confint(glmm_fit, method = "profile", parm = "theta_", level = 0.95)

# Combine summary and CIs into one clean dataframe
vc_df <- as.data.frame(lme4::VarCorr(glmm_fit))[, c("grp", "vcov", "sdcor")]
result <- data.frame(
  group    = vc_df$grp,
  variance = vc_df$vcov,
  sd       = vc_df$sdcor,
  ci_lower_sd  = ci_profile[, 1],
  ci_upper_sd  = ci_profile[, 2],
  ci_lower_var = ci_profile[, 1]^2,
  ci_upper_var = ci_profile[, 2]^2,
  row.names = NULL
)

# Obtain variance components####
vpc <- VCA(glmm_fit,
           subject_label  = "user_id",
           stimulus_label = "StimID",
           ci             = FALSE)

# Results####
print(result, digits = 2)
# group          variance    sd ci_lower_sd ci_upper_sd ci_lower_var ci_upper_var
# user_id:StimID    0.894  0.95        0.79        1.10         0.63         1.20
#        user_id    0.094  0.31        0.00        0.51         0.00         0.26
#         StimID    0.389  0.62        0.50        0.77         0.25         0.60

print(vpc, digits = 1)
# VPC total
# 1            VPC_Stimulus  0.08
# 2          VPC_Individual  0.02
# 3 VPC_Stimulus*Individual  0.19
# 4            VPC_Residual  0.71
# 5                      b1  0.70 # THIS IS THE ESTIMATE OF INTEREST (UNIQUE)
# 6                      b2  0.72
# 7                     mb1  0.30 # THIS IS THE ESTIMATE OF INTEREST (SHARED)
# 8                     mb2  0.28

# Sensitivity analysis####
## Intra-rater reliability ####
intra_rater <- repeated_trials%>%
  group_by(StimID, user_id, Species, Category) %>%
  filter(n() == 2) %>%
  summarise(agreement = as.integer(mean(correct) != 0.5), .groups = "drop") %>%
  group_by(user_id, Species, Category) %>% # note that here we average across individuals, not stimuli
  summarise(agreement = mean(agreement), .groups = "drop") 
low_agreement_ids <- intra_rater %>%
  filter(agreement == 0) %>%
  pull(user_id)
print(length(unique(low_agreement_ids)))
# There are quite a few individuals who disagree completely with themselves ("unreliable judges")

# Ratio of unreliable individuals
length(unique(low_agreement_ids))/length(unique(repeated_trials$user_id))

# Re-fit intercept-only GLMM to only reliable raters
rel_glmm_fit <- glmer(
  correct ~ 1 + (1 | user_id) + (1 | StimID) + (1 | user_id:StimID),
  family = binomial("logit"),
  data   = subset(repeated_trials, !(repeated_trials$user_id %in% low_agreement_ids))
)
summary(rel_glmm_fit)

# Obtain variance components
rel_vpc <- VCA(rel_glmm_fit,
           subject_label  = "user_id",
           stimulus_label = "StimID",
           ci             = FALSE)
print(rel_vpc, digits = 1)

## Animal strenght preferences ####
# Test results across levels of animal preferences as in James et al.
# Re-fit intercept-only GLMM to AnimalStrength>=66
strong_glmm_fit <- glmer(
  correct ~ 1 + (1 | user_id) + (1 | StimID) + (1 | user_id:StimID),
  family = binomial("logit"),
  data   = subset(repeated_trials,AnimalStrength>=66)
)
summary(strong_glmm_fit)

strongest_glmm_fit <- glmer(
  correct ~ 1 + (1 | user_id) + (1 | StimID) + (1 | user_id:StimID),
  family = binomial("logit"),
  data   = subset(repeated_trials,AnimalStrength>=75)
)
summary(strongest_glmm_fit)

# Obtain variance components
vpc_strong <- VCA(strong_glmm_fit,
                  subject_label  = "user_id",
                  stimulus_label = "StimID",
                  ci             = FALSE)

vpc_strongest <- VCA(strongest_glmm_fit,
                     subject_label  = "user_id",
                     stimulus_label = "StimID",
                     ci             = FALSE)

print(vpc_strong, digits = 1)
print(vpc_strongest, digits = 1)
