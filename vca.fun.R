# Author: Giacomo Bignardi
# Date: 2026-03-24
# Adapted from https://github.com/giacomobignardi/empirical-aesthetics-VCA

# ------------------------------------------------------------------------------
# VCA()
# Computes Variance Components Analysis.
# Supports LMM (lmer) and GLMM (glmer).
#
# Args:
#   model         : a fitted lmer or glmer model object
#   subject_label : character string matching the subject random effect grouping
#                   factor as it appears in the model (e.g. "Sub" | "user_id")
#   stimulus_label: character string matching the stimulus random effect grouping
#                   factor as it appears in the model (e.g. "Obj" | "stim_id")
#   ci            : if TRUE, returns a named numeric vector for
#                   bootstrapping; if FALSE, returns a tidy data.frame with
#                   Variance Partitioning Coefficient (VPC) and Beholder Indices (bi)
# References: 
# Martinez, J. E., Funk, F., & Todorov, A. (2020). Quantifying idiosyncratic and shared contributions to judgment. Behavior Research Methods, 52(4), 1428-1444.
# Hönekopp, J. (2006). Once more: is beauty in the eye of the beholder? Relative contributions of private and shared taste to judgments of facial attractiveness. Journal of Experimental Psychology: Human Perception and Performance, 32(2), 199.
# ------------------------------------------------------------------------------

VCA <- function(model, subject_label, stimulus_label, ci = FALSE) {
  
  ## Detect model type ####
  is_lmm  <- lme4::isLMM(model)
  is_glmm <- lme4::isGLMM(model)
  if (!is_lmm && !is_glmm) stop("Model must be a fitted lmer or glmer object.")
  
  ## Extract variance components via VarCorr (works for lmer & glmer) ####
  vc_df <- as.data.frame(lme4::VarCorr(model))
  
  get_var <- function(grp) {
    row <- vc_df[vc_df$grp == grp & is.na(vc_df$var2), ]
    if (nrow(row) == 0) stop(paste("Random effect group not found:", grp))
    row$vcov
  }
  
  # Build interaction label: try both orderings in case of label mismatch
  interaction_label <- if (
    any(vc_df$grp == paste0(subject_label, ":", stimulus_label))
  ) {
    paste0(subject_label, ":", stimulus_label)
  } else if (
    any(vc_df$grp == paste0(stimulus_label, ":", subject_label))
  ) {
    paste0(stimulus_label, ":", subject_label)
  } else {
    stop(paste(
      "Interaction random effect not found for labels:",
      subject_label, "and", stimulus_label
    ))
  }
  
  var_sub     <- get_var(subject_label)
  var_obj     <- get_var(stimulus_label)
  var_sub_obj <- get_var(interaction_label)
  
  ## 3. Residual variance ####
  # LMM    : sigma^2 from the model
  # GLMM   : distribution-specific latent-scale approximation
  #   logit  link -> pi^2 / 3   (logistic distribution variance)
  #   probit link -> 1          (standard normal variance)
  #   log    link -> 0          (Poisson; no residual variance on latent scale)
  var_residual <- if (is_lmm) {
    lme4::sigma(model)^2
  } else {
    fam <- family(model)
    switch(fam$link,
           logit  = pi^2 / 3,
           probit = 1,
           log    = 0,
           stop(paste("Unsupported GLMM link function:", fam$link))
    )
  }
  
  ## Partition variance ####
  total_variance      <- var_sub + var_obj + var_sub_obj + var_residual
  repeatable_variance <- var_sub + var_obj + var_sub_obj
  
  ## Variance components over total variance ####
  vpc_obj      <- var_obj      / total_variance
  vpc_sub      <- var_sub      / total_variance
  vpc_sub_obj  <- var_sub_obj  / total_variance
  vpc_residual <- var_residual / total_variance
  
  ## Summary indices (Honekoepp / modified Beholder) ####
  b1       <- var_sub_obj / (var_sub_obj + var_obj)                             # Beholder index 1
  b2       <- (var_sub_obj + var_sub) / (var_sub_obj + var_sub + var_obj )      # Beholder index 2
  mb1       <- 1 - b1                                                           # Modified Beholder index 1
  mb2       <- 1 - b2                                                           # Modified Beholder index 2
  
  # Output ####
  if (!ci) {
    data.frame(
      VPC = c(
        "VPC_Stimulus", "VPC_Individual", "VPC_Stimulus*Individual", "VPC_Residual",
        "b1", "b2", "mb1", "mb2"
      ),
      total = c(
        vpc_obj, vpc_sub, vpc_sub_obj, vpc_residual,
        b1, b2, mb1, mb2
      )
    )
  } else {
    # Named vector for downstream use in bootstrapping
    c(
      vpc_obj, vpc_sub, vpc_sub_obj, vpc_residual,
      b1, b2, mb1, mb2
    )
  }
}