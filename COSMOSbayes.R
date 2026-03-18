# ==============================================================================
# Project: Bayesian Interpretation for a Large-scale Trial using Mechanistically 
#          Informed Priors: Effects of Cocoa Extract on Cardiovascular Disease
# Author: Rikuta Hamaya
# 
# Description: This script conducts a Bayesian survival analysis to evaluate 
# the effect of cocoa extract on cardiovascular disease (CVD). It calculates 
# mechanistically informed priors based on SBP, FMD, and CRP pathways, fits 
# Bayesian Weibull Proportional Hazards models, and generates summary tables, 
# risk differences, and posterior likelihood visualizations.
# ==============================================================================

# Load necessary libraries
library(tidyverse)
library(haven)
library(survHE)
library(rstan)
library(bayesplot)
library(tidybayes)
library(scales)
library(MASS)
library(survminer)

# Set working directory (Update this to your local path before running)
# setwd("/path/to/your/directory")


# 1. PRIORS GENERATION ----------------------------------------------------
# These calculations are based on previous literature.
# Priors were defined by three mechanisms: CF-BP-CVD, CF-FMD-CVD, and CF-CRP-CVD
# Please see the methods section of the manuscript for exact derivations.

set.seed(42)          # Set seed for reproducibility
n_sim <- 100000       # Number of Monte Carlo samples for prior distributions

# --- 1.1 BP pathway (SBP reduction) ---
# PriorSBP-CVD: HR per 5 mmHg reduction = 0.91 (0.88, 0.95)
HR_bp5   <- 0.91
L_bp5    <- 0.88
U_bp5    <- 0.95
logHR_bp5 <- log(HR_bp5)
logL_bp5  <- log(L_bp5)
logU_bp5  <- log(U_bp5)
SE_logHR_bp5 <- ((logHR_bp5 - logL_bp5)/1.96 + (logU_bp5 - logHR_bp5)/1.96)/2

# Per 1 mmHg reduction (note: reduction is negative change -> logHR negative)
logHR_bp_per1 <- logHR_bp5 / -5          # ≈ -0.01886
SE_bp_per1    <- SE_logHR_bp5 / 5        # ≈ 0.00392

# Cocoa -> SBP: mean -2.37 mmHg (CI -4.30 to -0.44)
cocoa_bp_mean  <- -2.37
cocoa_bp_lower <- -4.30                  # larger reduction
cocoa_bp_upper <- -0.44                  # smaller reduction

Est_BP   <- logHR_bp_per1 * cocoa_bp_mean
H_est_BP <- logHR_bp_per1 * cocoa_bp_upper
L_est_BP <- logHR_bp_per1 * cocoa_bp_lower

cat("BP pathway (log-HR):\n")
cat("Est_BP   =", round(Est_BP, 3), "\n")
cat("H_est_BP =", round(H_est_BP, 3), "\n")
cat("L_est_BP =", round(L_est_BP, 3), "\n")

# Approximate SE_BP propagating uncertainty
SE_BP <- ((H_est_BP + 1.96*SE_bp_per1 - Est_BP)/1.96 + 
            (Est_BP - (L_est_BP - 1.96*SE_bp_per1))/1.96)/2

cat("SE_BP ≈", round(SE_BP, 5), "\n")
cat("HR point =", round(exp(Est_BP), 3), "\n")
cat("95% CI  =", round(exp(Est_BP - 1.96*SE_BP), 3), "-", 
    round(exp(Est_BP + 1.96*SE_BP), 3), "\n\n")

# --- 1.2 FMD pathway ---
HR_fmd   <- 0.88
L_fmd    <- 0.84
U_fmd    <- 0.91
logHR_fmd <- log(HR_fmd)
logL_fmd  <- log(L_fmd)
logU_fmd  <- log(U_fmd)
SE_fmd    <- ((logHR_fmd - logL_fmd)/1.96 + (logU_fmd - logHR_fmd)/1.96)/2

# Cocoa -> FMD: mean +1.34% (CI 1.00 to 1.68)
cocoa_fmd_mean  <- 1.34
cocoa_fmd_lower <- 1.00
cocoa_fmd_upper <- 1.68

Est_FMD   <- logHR_fmd * cocoa_fmd_mean
H_est_FMD <- logHR_fmd * cocoa_fmd_lower
L_est_FMD <- logHR_fmd * cocoa_fmd_upper

cat("FMD pathway (log-HR):\n")
cat("Est_FMD   =", round(Est_FMD, 3), "\n")
cat("H_est_FMD =", round(H_est_FMD, 3), "\n")
cat("L_est_FMD =", round(L_est_FMD, 3), "\n")

SE_FMD <- ((H_est_FMD + 1.96*SE_fmd - Est_FMD)/1.96 + 
             (Est_FMD - (L_est_FMD - 1.96*SE_fmd))/1.96)/2

cat("SE_FMD ≈", round(SE_FMD, 5), "\n")
cat("HR point =", round(exp(Est_FMD), 3), "\n")
cat("95% CI  =", round(exp(Est_FMD - 1.96*SE_FMD), 3), "-", 
    round(exp(Est_FMD + 1.96*SE_FMD), 3), "\n\n")

# --- 1.3 CRP pathway ---
HR_crp   <- 0.86
L_crp    <- 0.75
U_crp    <- 0.99
logHR_crp <- log(HR_crp)
logL_crp  <- log(L_crp)
logU_crp  <- log(U_crp)

reduction_crp <- 1.85
logHR_per1_crp <- logHR_crp / -reduction_crp     # per 1 mg/L reduction
SE_per1_crp    <- (((logHR_crp - logL_crp)/1.96 + (logU_crp - logHR_crp)/1.96)/2) / reduction_crp

# Cocoa -> CRP: mean -0.98 mg/L (CI -1.69 to -0.27)
cocoa_crp_mean  <- -0.978   
cocoa_crp_lower <- -1.687
cocoa_crp_upper <- -0.269

Est_CRP   <- logHR_per1_crp * cocoa_crp_mean
H_est_CRP <- logHR_per1_crp * cocoa_crp_upper
L_est_CRP <- logHR_per1_crp * cocoa_crp_lower

cat("CRP pathway (log-HR):\n")
cat("Est_CRP   =", round(Est_CRP, 3), "\n")
cat("H_est_CRP =", round(H_est_CRP, 3), "\n")
cat("L_est_CRP =", round(L_est_CRP, 3), "\n")

SE_CRP <- ((H_est_CRP + 1.96*SE_per1_crp - Est_CRP)/1.96 + 
             (Est_CRP - (L_est_CRP - 1.96*SE_per1_crp))/1.96)/2

cat("SE_CRP ≈", round(SE_CRP, 5), "\n")
cat("HR point =", round(exp(Est_CRP), 3), "\n")
cat("95% CI  =", round(exp(Est_CRP - 1.96*SE_CRP), 3), "-", 
    round(exp(Est_CRP + 1.96*SE_CRP), 3), "\n\n")

# --- 1.4 Sample the pathway-specific effects (log-HR scale) ---
BP_samples  <- rnorm(n_sim, Est_BP,  SE_BP)
FMD_samples <- rnorm(n_sim, Est_FMD, SE_FMD)
CRP_samples <- rnorm(n_sim, Est_CRP, SE_CRP)

cat("Sampled pathway means (should ≈ Est_...):\n")
cat("BP  mean =", round(mean(BP_samples), 4), " | sd =", round(sd(BP_samples), 4), "\n")
cat("FMD mean =", round(mean(FMD_samples), 4), " | sd =", round(sd(FMD_samples), 4), "\n")
cat("CRP mean =", round(mean(CRP_samples), 4), " | sd =", round(sd(CRP_samples), 4), "\n\n")

# --- 1.5 Construct uncalibrated priors for each scenario (log-HR scale) ---

# Scenario 1: solely via SBP (BP)
prior_1 <- BP_samples

# Scenario 2: solely via CRP
prior_2 <- CRP_samples

# Scenario 3: independently via SBP + CRP (additive, no shared adjustment)
prior_3 <- BP_samples + CRP_samples

# Scenario 4: SBP + CRP with 50% shared effects
# Reduces mean benefit to reflect redundant overlap
prior_4 <- (1 / (1 + 0.5 * (2 - 1))) * (BP_samples + CRP_samples)

# Scenario 5: all three pathways with 50% shared effects
prior_5 <- (1 / (1 + 0.5 * (3 - 1))) * (BP_samples + CRP_samples + FMD_samples)

# Skeptical: centered at HR=1 (log-HR=0)
prior_skeptical <- rnorm(n_sim, mean = 0, sd = 0.1)

# Pessimistic: centered at HR=1.1 (log-HR≈0.0953)
prior_pessimistic <- rnorm(n_sim, mean = log(1.1), sd = 0.1)

# List of all priors
priors_uncalib <- list(
  "SBP only" = prior_1,
  "CRP only" = prior_2,
  "SBP+CRP independent" = prior_3,
  "SBP+CRP 50% shared" = prior_4,
  "SBP+CRP+FMD 50% shared" = prior_5,
  "Skeptical (HR=1)" = prior_skeptical,
  "Pessimistic (HR=1.1)" = prior_pessimistic
)

# Summarize uncalibrated priors
for (name in names(priors_uncalib)) {
  p <- priors_uncalib[[name]]
  mean_p <- mean(p)
  sd_p <- sd(p)
  ess_p <- 4 / sd_p^2
  cat(name, ":\n")
  cat("  Mean log-HR =", round(mean_p, 4), "\n")
  cat("  SD          =", round(sd_p, 4), "\n")
  cat("  ESS approx  =", round(ess_p), "\n")
  cat("  HR point    =", round(exp(mean_p), 3), "\n")
  cat("  95% HR CI   =", round(exp(mean_p - 1.96*sd_p), 3), "-", 
      round(exp(mean_p + 1.96*sd_p), 3), "\n\n")
}

# --- 1.6 Calibrate each prior to target ESS 30% ---
target_ess <- 866 * 0.3
sd_target  <- sqrt(4 / target_ess) 

priors_calib <- list()
for (name in names(priors_uncalib)) {
  p_uncalib <- priors_uncalib[[name]]
  mean_p <- mean(p_uncalib)
  sd_uncalib <- sd(p_uncalib)
  inflation_factor <- sd_target / sd_uncalib
  p_calib <- mean_p + inflation_factor * (p_uncalib - mean_p)
  priors_calib[[name]] <- p_calib
}

set.seed(42)
n_sim <- length(priors_calib[[1]])
non_informative <- rnorm(n_sim, mean = 0, sd = 10)  # sd=10 -> extremely wide

# Combine into a final named list
all_priors_calib <- c(
  list("0 Non-informative" = non_informative),
  priors_calib
)

# Function to summarize and export priors
summarize_prior <- function(samples, scenario_name) {
  mean_loghr <- mean(samples)
  sd_loghr   <- sd(samples)
  
  hr_point   <- exp(mean_loghr)
  hr_lower   <- exp(mean_loghr - 1.96 * sd_loghr)
  hr_upper   <- exp(mean_loghr + 1.96 * sd_loghr)
  
  if (sd_loghr > 5) {
    ci_str <- "(0, 3e8)"
    hr_str <- "1.00"
  } else {
    ci_str <- sprintf("(%.2f, %.2f)", round(hr_lower, 2), round(hr_upper, 2))
    hr_str <- sprintf("%.2f %s", round(hr_point, 2), ci_str)
  }
  
  tibble(
    Scenario = scenario_name,
    `Mean prior HR for CVD (95% CI)` = hr_str,
    Rational = NA_character_
  )
}

prior_table <- map_dfr(
  names(all_priors_calib),
  ~ summarize_prior(all_priors_calib[[.x]], .x)
)

prior_table$Rational=
  c(
    "No prior information",
    "Primarily through BP-lowering effects",
    "Primarily through anti-inflammatory effects",
    "Independently through BP-lowering and anti-inflammatory effects",
    "Through BP-lowering and anti-inflammatory effects (50% shared)",
    "Through BP-lowering, anti-inflammatory, and endothelial improvement (50% shared)",
    "Belief in no effects (HR=1)",
    "Belief in harmful effects (HR=1.1)"
  )

print(prior_table, width = Inf)
# write.csv(prior_table,"prior_summary_table.csv", row.names = FALSE)


# 2. DATA PREPARATION -----------------------------------------------------

# -------------------------------------------------------------------------
# NOTE FOR REPRODUCIBILITY:
# The original dataset from the trial is not publicly available. 
# To run the code below, you must provide a dataframe named `main_use` 
# containing the following required columns:
#
# t2mistrrevcvddthoth : Numeric. Follow-up time to CVD event or censoring.
# mistrrevcvddthoth   : Binary (1/0). Event indicator (1 = event, 0 = censored).
# interventionC       : Binary (1/0). Cocoa extract intervention arm.
# interventionM       : Binary (1/0). Multivitamin intervention arm.
# agerand             : Numeric. Age at randomization.
# men                 : Binary (1/0). Sex indicator (1 = Male, 0 = Female).
# whi                 : Binary (1/0). WHI cohort indicator.
# -------------------------------------------------------------------------

# Example placeholder for loading local data
# main_use <- read.csv("path/to/your/dataset.csv")

# Ensure time-to-event is strictly positive for survival analysis
# main_use <- main_use %>% filter(t2mistrrevcvddthoth > 0)


# 3. FREQUENTIST ANALYSIS -------------------------------------------------
# Generate standard baseline models for comparison

# Cox Proportional Hazards Model
# coxph(Surv(t2mistrrevcvddthoth, mistrrevcvddthoth) ~ interventionC + agerand + men + interventionM + whi, data = main_use) %>% summary()

# Weibull Model
# survreg(Surv(t2mistrrevcvddthoth, mistrrevcvddthoth) ~ interventionC + agerand + men + interventionM + whi, dist="weibull", data=main_use) %>% summary()


# 4. BAYESIAN ANALYSIS (survHE) -------------------------------------------
# We fit Bayesian Weibull Proportional Hazards models corresponding to each prior scenario.

# Define the base formula
form <- Surv(t2mistrrevcvddthoth, mistrrevcvddthoth) ~ interventionC + agerand + men + interventionM + whi

# 4.1 Non-informative Prior
CVD_ni = fit.models(formula = form, data = main_use, distr = c("weibullPH"), method = "hmc",
                    priors = list(
                      exp = list(mu_beta = rep(0,6), sigma_beta= rep(10,6)),
                      weibullPH = list(a_alpha=0.1, b_alpha=0.1, mu_beta = rep(0,6), sigma_beta= rep(10,6))
                    ),
                    chains=2, iter=3000, seed=1233)

# 4.2 SBP Only Prior
CVD_sbp = fit.models(formula = form, data = main_use, distr = c("weibullPH"), method = "hmc",
                     priors = list(
                       weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                        mu_beta = c(0, mean(all_priors_calib$`SBP only`), 0,0,0,0), 
                                        sigma_beta= c(10, sd(all_priors_calib$`SBP only`), 10,10,10,10))
                     ),
                     chains=2, iter=3000, seed=1233)

# 4.3 CRP Only Prior
CVD_crp = fit.models(formula = form, data = main_use, distr = c("weibullPH"), method = "hmc",
                     priors = list(
                       weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                        mu_beta = c(0, mean(all_priors_calib$`CRP only`), 0,0,0,0), 
                                        sigma_beta= c(10, sd(all_priors_calib$`CRP only`), 10,10,10,10))
                     ),
                     chains=2, iter=3000, seed=1233)

# 4.4 SBP + CRP Independent Prior
CVD_sbpcrp = fit.models(formula = form, data = main_use, distr = c("weibullPH"), method = "hmc",
                        priors = list(
                          weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                           mu_beta = c(0, mean(all_priors_calib$`SBP+CRP independent`), 0,0,0,0), 
                                           sigma_beta= c(10, sd(all_priors_calib$`SBP+CRP independent`), 10,10,10,10))
                        ),
                        chains=2, iter=3000, seed=1233)

# 4.5 SBP + CRP (50% Shared) Prior
CVD_sbpcrp50 = fit.models(formula = form, data = main_use, distr = c("weibullPH"), method = "hmc",
                          priors = list(
                            weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                             mu_beta = c(0, mean(all_priors_calib$`SBP+CRP 50% shared`), 0,0,0,0), 
                                             sigma_beta= c(10, sd(all_priors_calib$`SBP+CRP 50% shared`), 10,10,10,10))
                          ),
                          chains=2, iter=3000, seed=1233)

# 4.6 SBP + CRP + FMD (50% Shared) Prior
CVD_sbpcrpfmd = fit.models(formula = form, data = main_use, distr = c("weibullPH"), method = "hmc",
                           priors = list(
                             weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                              mu_beta = c(0, mean(all_priors_calib$`SBP+CRP+FMD 50% shared`), 0,0,0,0), 
                                              sigma_beta= c(10, sd(all_priors_calib$`SBP+CRP+FMD 50% shared`), 10,10,10,10))
                           ),
                           chains=2, iter=3000, seed=1233)

# 4.7 Skeptical Prior
CVD_skep = fit.models(formula = form, data = main_use, distr = c("weibullPH"), method = "hmc",
                      priors = list(
                        weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                         mu_beta = c(0, mean(all_priors_calib$`Skeptical (HR=1)`), 0,0,0,0), 
                                         sigma_beta= c(10, sd(all_priors_calib$`Skeptical (HR=1)`), 10,10,10,10))
                      ),
                      chains=2, iter=3000, seed=1233)

# 4.8 Pessimistic Prior
CVD_pess = fit.models(formula = form, data = main_use, distr = c("weibullPH"), method = "hmc",
                      priors = list(
                        weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                         mu_beta = c(0, mean(all_priors_calib$`Pessimistic (HR=1.1)`), 0,0,0,0), 
                                         sigma_beta= c(10, sd(all_priors_calib$`Pessimistic (HR=1.1)`), 10,10,10,10))
                      ),
                      chains=2, iter=3000, seed=1233)

# 4.9 Model for K-M curve (Intervention only)
cvd9km = fit.models(formula = Surv(t2mistrrevcvddthoth,mistrrevcvddthoth) ~ as.factor(interventionC),
                    data = main_use, distr = c("weibullPH"), method = "hmc",
                    priors = list(
                      weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                       mu_beta = c(0, mean(all_priors_calib$`SBP+CRP+FMD 50% shared`)), 
                                       sigma_beta= c(10, sd(all_priors_calib$`SBP+CRP+FMD 50% shared`)))
                    ),
                    chains=2, iter=3000, seed=1233)

# Group all CVD model fits into a list
CVD <- grep("CVD_", names(.GlobalEnv), value=TRUE)
CVD_list <- do.call("list", mget(CVD))


# 5. MAIN RESULTS (TABLE 3) -----------------------------------------------

# Helper function to extract Hazard Ratios and 95% Credible Intervals
comb = function(res){
  temp = res$models$`Weibull (PH)`
  n4 = temp@.MISC$summary$quan[2,5] # Median HR
  n5 = temp@.MISC$summary$quan[2,1] # Lower CI
  n6 = temp@.MISC$summary$quan[2,9] # Upper CI
  
  print(paste0('Weibull: ',
               round(exp(n4),2)," (",
               round(exp(n5),2),", ",
               round(exp(n6),2),")"))
}

# Print estimates for each model scenario
for (i in 1:length(CVD_list)){
  print(names(CVD_list)[[i]])
  comb(CVD_list[[i]])
}


# 6. RISK DIFFERENCE CALCULATION ------------------------------------------
# Calculates absolute risk difference (RD) at a specific evaluation time.

t_eval <- 365*4 # 4 years follow-up

RDcalc = function(cvdfit){
  temp = cvdfit$models$`Weibull (PH)`
  post <- rstan::extract(temp)
  alpha <- post$alpha
  beta <- post$beta
  
  # Evaluate covariates at their means (Adjust based on your actual data structure)
  # intercept, interventionC, agerand, men, interventionM, whi
  X_control <- c(1, 0, mean(main_use$agerand), mean(main_use$men), 0, 0) 
  X_treat   <- c(1, 1, mean(main_use$agerand), mean(main_use$men), 0, 0)
  
  # Calculate linear predictors
  lp_control <- as.vector(beta %*% X_control)
  lp_treat   <- as.vector(beta %*% X_treat)
  
  S_control <- exp(- (t_eval ^ alpha) * exp(lp_control))
  S_treat   <- exp(- (t_eval ^ alpha) * exp(lp_treat))
  
  risk_control <- 1 - S_control
  risk_treat <- 1 - S_treat
  RD <- risk_treat - risk_control
  
  return(c(mean(RD), quantile(RD, c(0.025, 0.975))))
}

# Apply to all models and build data frame
rd_results <- lapply(CVD_list, RDcalc)
rd_df <- do.call(rbind, rd_results) %>% as.data.frame()
colnames(rd_df) <- c("RD_mean", "RD_lower", "RD_upper")

rd_df$est = paste0(round(rd_df$RD_mean*100,2), " (", round(rd_df$RD_lower*100,2), ", ", round(rd_df$RD_upper*100,2), ")")
# write.csv(rd_df,"rd_df.csv")


# 7. POSTERIOR LIKELIHOOD# ==============================================================================
# Project: Bayesian Interpretation for a Large-scale Trial using Mechanistically 
#          Informed Priors: Effects of Cocoa Extract on Cardiovascular Disease
# Author: Rikuta Hamaya
# 
# Description: This script conducts a Bayesian survival analysis to evaluate 
# the effect of cocoa extract on cardiovascular disease (CVD). It calculates 
# mechanistically informed priors based on SBP, FMD, and CRP pathways, fits 
# Bayesian Weibull Proportional Hazards models, and generates summary tables, 
# risk differences, and posterior likelihood visualizations.
# ==============================================================================

# Load necessary libraries
library(tidyverse)
library(haven)
library(survHE)
library(rstan)
library(bayesplot)
library(tidybayes)
library(scales)
library(MASS)
library(survminer)

# Set working directory (Update this to your local path before running)
# setwd("/path/to/your/directory")


# 1. PRIORS GENERATION ----------------------------------------------------
# These calculations are based on previous literature.
# Priors were defined by three mechanisms: CF-BP-CVD, CF-FMD-CVD, and CF-CRP-CVD
# Please see the methods section of the manuscript for exact derivations.

set.seed(42)          # Set seed for reproducibility
n_sim <- 100000       # Number of Monte Carlo samples for prior distributions

# --- 1.1 BP pathway (SBP reduction) ---
# PriorSBP-CVD: HR per 5 mmHg reduction = 0.91 (0.88, 0.95)
HR_bp5   <- 0.91
L_bp5    <- 0.88
U_bp5    <- 0.95
logHR_bp5 <- log(HR_bp5)
logL_bp5  <- log(L_bp5)
logU_bp5  <- log(U_bp5)
SE_logHR_bp5 <- ((logHR_bp5 - logL_bp5)/1.96 + (logU_bp5 - logHR_bp5)/1.96)/2

# Per 1 mmHg reduction (note: reduction is negative change -> logHR negative)
logHR_bp_per1 <- logHR_bp5 / -5          # ≈ -0.01886
SE_bp_per1    <- SE_logHR_bp5 / 5        # ≈ 0.00392

# Cocoa -> SBP: mean -2.37 mmHg (CI -4.30 to -0.44)
cocoa_bp_mean  <- -2.37
cocoa_bp_lower <- -4.30                  # larger reduction
cocoa_bp_upper <- -0.44                  # smaller reduction

Est_BP   <- logHR_bp_per1 * cocoa_bp_mean
H_est_BP <- logHR_bp_per1 * cocoa_bp_upper
L_est_BP <- logHR_bp_per1 * cocoa_bp_lower

cat("BP pathway (log-HR):\n")
cat("Est_BP   =", round(Est_BP, 3), "\n")
cat("H_est_BP =", round(H_est_BP, 3), "\n")
cat("L_est_BP =", round(L_est_BP, 3), "\n")

# Approximate SE_BP propagating uncertainty
SE_BP <- ((H_est_BP + 1.96*SE_bp_per1 - Est_BP)/1.96 + 
            (Est_BP - (L_est_BP - 1.96*SE_bp_per1))/1.96)/2

cat("SE_BP ≈", round(SE_BP, 5), "\n")
cat("HR point =", round(exp(Est_BP), 3), "\n")
cat("95% CI  =", round(exp(Est_BP - 1.96*SE_BP), 3), "-", 
    round(exp(Est_BP + 1.96*SE_BP), 3), "\n\n")

# --- 1.2 FMD pathway ---
HR_fmd   <- 0.88
L_fmd    <- 0.84
U_fmd    <- 0.91
logHR_fmd <- log(HR_fmd)
logL_fmd  <- log(L_fmd)
logU_fmd  <- log(U_fmd)
SE_fmd    <- ((logHR_fmd - logL_fmd)/1.96 + (logU_fmd - logHR_fmd)/1.96)/2

# Cocoa -> FMD: mean +1.34% (CI 1.00 to 1.68)
cocoa_fmd_mean  <- 1.34
cocoa_fmd_lower <- 1.00
cocoa_fmd_upper <- 1.68

Est_FMD   <- logHR_fmd * cocoa_fmd_mean
H_est_FMD <- logHR_fmd * cocoa_fmd_lower
L_est_FMD <- logHR_fmd * cocoa_fmd_upper

cat("FMD pathway (log-HR):\n")
cat("Est_FMD   =", round(Est_FMD, 3), "\n")
cat("H_est_FMD =", round(H_est_FMD, 3), "\n")
cat("L_est_FMD =", round(L_est_FMD, 3), "\n")

SE_FMD <- ((H_est_FMD + 1.96*SE_fmd - Est_FMD)/1.96 + 
             (Est_FMD - (L_est_FMD - 1.96*SE_fmd))/1.96)/2

cat("SE_FMD ≈", round(SE_FMD, 5), "\n")
cat("HR point =", round(exp(Est_FMD), 3), "\n")
cat("95% CI  =", round(exp(Est_FMD - 1.96*SE_FMD), 3), "-", 
    round(exp(Est_FMD + 1.96*SE_FMD), 3), "\n\n")

# --- 1.3 CRP pathway ---
HR_crp   <- 0.86
L_crp    <- 0.75
U_crp    <- 0.99
logHR_crp <- log(HR_crp)
logL_crp  <- log(L_crp)
logU_crp  <- log(U_crp)

reduction_crp <- 1.85
logHR_per1_crp <- logHR_crp / -reduction_crp     # per 1 mg/L reduction
SE_per1_crp    <- (((logHR_crp - logL_crp)/1.96 + (logU_crp - logHR_crp)/1.96)/2) / reduction_crp

# Cocoa -> CRP: mean -0.98 mg/L (CI -1.69 to -0.27)
cocoa_crp_mean  <- -0.978   
cocoa_crp_lower <- -1.687
cocoa_crp_upper <- -0.269

Est_CRP   <- logHR_per1_crp * cocoa_crp_mean
H_est_CRP <- logHR_per1_crp * cocoa_crp_upper
L_est_CRP <- logHR_per1_crp * cocoa_crp_lower

cat("CRP pathway (log-HR):\n")
cat("Est_CRP   =", round(Est_CRP, 3), "\n")
cat("H_est_CRP =", round(H_est_CRP, 3), "\n")
cat("L_est_CRP =", round(L_est_CRP, 3), "\n")

SE_CRP <- ((H_est_CRP + 1.96*SE_per1_crp - Est_CRP)/1.96 + 
             (Est_CRP - (L_est_CRP - 1.96*SE_per1_crp))/1.96)/2

cat("SE_CRP ≈", round(SE_CRP, 5), "\n")
cat("HR point =", round(exp(Est_CRP), 3), "\n")
cat("95% CI  =", round(exp(Est_CRP - 1.96*SE_CRP), 3), "-", 
    round(exp(Est_CRP + 1.96*SE_CRP), 3), "\n\n")

# --- 1.4 Sample the pathway-specific effects (log-HR scale) ---
BP_samples  <- rnorm(n_sim, Est_BP,  SE_BP)
FMD_samples <- rnorm(n_sim, Est_FMD, SE_FMD)
CRP_samples <- rnorm(n_sim, Est_CRP, SE_CRP)

cat("Sampled pathway means (should ≈ Est_...):\n")
cat("BP  mean =", round(mean(BP_samples), 4), " | sd =", round(sd(BP_samples), 4), "\n")
cat("FMD mean =", round(mean(FMD_samples), 4), " | sd =", round(sd(FMD_samples), 4), "\n")
cat("CRP mean =", round(mean(CRP_samples), 4), " | sd =", round(sd(CRP_samples), 4), "\n\n")

# --- 1.5 Construct uncalibrated priors for each scenario (log-HR scale) ---

# Scenario 1: solely via SBP (BP)
prior_1 <- BP_samples

# Scenario 2: solely via CRP
prior_2 <- CRP_samples

# Scenario 3: independently via SBP + CRP (additive, no shared adjustment)
prior_3 <- BP_samples + CRP_samples

# Scenario 4: SBP + CRP with 50% shared effects
# Reduces mean benefit to reflect redundant overlap
prior_4 <- (1 / (1 + 0.5 * (2 - 1))) * (BP_samples + CRP_samples)

# Scenario 5: all three pathways with 50% shared effects
prior_5 <- (1 / (1 + 0.5 * (3 - 1))) * (BP_samples + CRP_samples + FMD_samples)

# Skeptical: centered at HR=1 (log-HR=0)
prior_skeptical <- rnorm(n_sim, mean = 0, sd = 0.1)

# Pessimistic: centered at HR=1.1 (log-HR≈0.0953)
prior_pessimistic <- rnorm(n_sim, mean = log(1.1), sd = 0.1)

# List of all priors
priors_uncalib <- list(
  "SBP only" = prior_1,
  "CRP only" = prior_2,
  "SBP+CRP independent" = prior_3,
  "SBP+CRP 50% shared" = prior_4,
  "SBP+CRP+FMD 50% shared" = prior_5,
  "Skeptical (HR=1)" = prior_skeptical,
  "Pessimistic (HR=1.1)" = prior_pessimistic
)

# Summarize uncalibrated priors
for (name in names(priors_uncalib)) {
  p <- priors_uncalib[[name]]
  mean_p <- mean(p)
  sd_p <- sd(p)
  ess_p <- 4 / sd_p^2
  cat(name, ":\n")
  cat("  Mean log-HR =", round(mean_p, 4), "\n")
  cat("  SD          =", round(sd_p, 4), "\n")
  cat("  ESS approx  =", round(ess_p), "\n")
  cat("  HR point    =", round(exp(mean_p), 3), "\n")
  cat("  95% HR CI   =", round(exp(mean_p - 1.96*sd_p), 3), "-", 
      round(exp(mean_p + 1.96*sd_p), 3), "\n\n")
}

# --- 1.6 Calibrate each prior to target ESS 30% ---
target_ess <- 866 * 0.3
sd_target  <- sqrt(4 / target_ess) 

priors_calib <- list()
for (name in names(priors_uncalib)) {
  p_uncalib <- priors_uncalib[[name]]
  mean_p <- mean(p_uncalib)
  sd_uncalib <- sd(p_uncalib)
  inflation_factor <- sd_target / sd_uncalib
  p_calib <- mean_p + inflation_factor * (p_uncalib - mean_p)
  priors_calib[[name]] <- p_calib
}

set.seed(42)
n_sim <- length(priors_calib[[1]])
non_informative <- rnorm(n_sim, mean = 0, sd = 10)  # sd=10 -> extremely wide

# Combine into a final named list
all_priors_calib <- c(
  list("0 Non-informative" = non_informative),
  priors_calib
)

# Function to summarize and export priors
summarize_prior <- function(samples, scenario_name) {
  mean_loghr <- mean(samples)
  sd_loghr   <- sd(samples)
  
  hr_point   <- exp(mean_loghr)
  hr_lower   <- exp(mean_loghr - 1.96 * sd_loghr)
  hr_upper   <- exp(mean_loghr + 1.96 * sd_loghr)
  
  if (sd_loghr > 5) {
    ci_str <- "(0, 3e8)"
    hr_str <- "1.00"
  } else {
    ci_str <- sprintf("(%.2f, %.2f)", round(hr_lower, 2), round(hr_upper, 2))
    hr_str <- sprintf("%.2f %s", round(hr_point, 2), ci_str)
  }
  
  tibble(
    Scenario = scenario_name,
    `Mean prior HR for CVD (95% CI)` = hr_str,
    Rational = NA_character_
  )
}

prior_table <- map_dfr(
  names(all_priors_calib),
  ~ summarize_prior(all_priors_calib[[.x]], .x)
)

prior_table$Rational=
  c(
    "No prior information",
    "Primarily through BP-lowering effects",
    "Primarily through anti-inflammatory effects",
    "Independently through BP-lowering and anti-inflammatory effects",
    "Through BP-lowering and anti-inflammatory effects (50% shared)",
    "Through BP-lowering, anti-inflammatory, and endothelial improvement (50% shared)",
    "Belief in no effects (HR=1)",
    "Belief in harmful effects (HR=1.1)"
  )

print(prior_table, width = Inf)
# write.csv(prior_table,"prior_summary_table.csv", row.names = FALSE)


# 2. DATA PREPARATION -----------------------------------------------------

# -------------------------------------------------------------------------
# NOTE FOR REPRODUCIBILITY:
# The original dataset from the trial is not publicly available. 
# To run the code below, you must provide a dataframe named `main_use` 
# containing the following required columns:
#
# t2mistrrevcvddthoth : Numeric. Follow-up time to CVD event or censoring.
# mistrrevcvddthoth   : Binary (1/0). Event indicator (1 = event, 0 = censored).
# interventionC       : Binary (1/0). Cocoa extract intervention arm.
# interventionM       : Binary (1/0). Multivitamin intervention arm.
# agerand             : Numeric. Age at randomization.
# men                 : Binary (1/0). Sex indicator (1 = Male, 0 = Female).
# whi                 : Binary (1/0). WHI cohort indicator.
# -------------------------------------------------------------------------

# Example placeholder for loading local data
# main_use <- read.csv("path/to/your/dataset.csv")

# Ensure time-to-event is strictly positive for survival analysis
# main_use <- main_use %>% filter(t2mistrrevcvddthoth > 0)


# 3. FREQUENTIST ANALYSIS -------------------------------------------------
# Generate standard baseline models for comparison

# Cox Proportional Hazards Model
# coxph(Surv(t2mistrrevcvddthoth, mistrrevcvddthoth) ~ interventionC + agerand + men + interventionM + whi, data = main_use) %>% summary()

# Weibull Model
# survreg(Surv(t2mistrrevcvddthoth, mistrrevcvddthoth) ~ interventionC + agerand + men + interventionM + whi, dist="weibull", data=main_use) %>% summary()


# 4. BAYESIAN ANALYSIS (survHE) -------------------------------------------
# We fit Bayesian Weibull Proportional Hazards models corresponding to each prior scenario.

# Define the base formula
form <- Surv(t2mistrrevcvddthoth, mistrrevcvddthoth) ~ interventionC + agerand + men + interventionM + whi

# 4.1 Non-informative Prior
CVD_ni = fit.models(formula = form, data = main_use, distr = c("weibullPH"), method = "hmc",
                    priors = list(
                      exp = list(mu_beta = rep(0,6), sigma_beta= rep(10,6)),
                      weibullPH = list(a_alpha=0.1, b_alpha=0.1, mu_beta = rep(0,6), sigma_beta= rep(10,6))
                    ),
                    chains=2, iter=3000, seed=1233)

# 4.2 SBP Only Prior
CVD_sbp = fit.models(formula = form, data = main_use, distr = c("weibullPH"), method = "hmc",
                     priors = list(
                       weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                        mu_beta = c(0, mean(all_priors_calib$`SBP only`), 0,0,0,0), 
                                        sigma_beta= c(10, sd(all_priors_calib$`SBP only`), 10,10,10,10))
                     ),
                     chains=2, iter=3000, seed=1233)

# 4.3 CRP Only Prior
CVD_crp = fit.models(formula = form, data = main_use, distr = c("weibullPH"), method = "hmc",
                     priors = list(
                       weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                        mu_beta = c(0, mean(all_priors_calib$`CRP only`), 0,0,0,0), 
                                        sigma_beta= c(10, sd(all_priors_calib$`CRP only`), 10,10,10,10))
                     ),
                     chains=2, iter=3000, seed=1233)

# 4.4 SBP + CRP Independent Prior
CVD_sbpcrp = fit.models(formula = form, data = main_use, distr = c("weibullPH"), method = "hmc",
                        priors = list(
                          weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                           mu_beta = c(0, mean(all_priors_calib$`SBP+CRP independent`), 0,0,0,0), 
                                           sigma_beta= c(10, sd(all_priors_calib$`SBP+CRP independent`), 10,10,10,10))
                        ),
                        chains=2, iter=3000, seed=1233)

# 4.5 SBP + CRP (50% Shared) Prior
CVD_sbpcrp50 = fit.models(formula = form, data = main_use, distr = c("weibullPH"), method = "hmc",
                          priors = list(
                            weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                             mu_beta = c(0, mean(all_priors_calib$`SBP+CRP 50% shared`), 0,0,0,0), 
                                             sigma_beta= c(10, sd(all_priors_calib$`SBP+CRP 50% shared`), 10,10,10,10))
                          ),
                          chains=2, iter=3000, seed=1233)

# 4.6 SBP + CRP + FMD (50% Shared) Prior
CVD_sbpcrpfmd = fit.models(formula = form, data = main_use, distr = c("weibullPH"), method = "hmc",
                           priors = list(
                             weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                              mu_beta = c(0, mean(all_priors_calib$`SBP+CRP+FMD 50% shared`), 0,0,0,0), 
                                              sigma_beta= c(10, sd(all_priors_calib$`SBP+CRP+FMD 50% shared`), 10,10,10,10))
                           ),
                           chains=2, iter=3000, seed=1233)

# 4.7 Skeptical Prior
CVD_skep = fit.models(formula = form, data = main_use, distr = c("weibullPH"), method = "hmc",
                      priors = list(
                        weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                         mu_beta = c(0, mean(all_priors_calib$`Skeptical (HR=1)`), 0,0,0,0), 
                                         sigma_beta= c(10, sd(all_priors_calib$`Skeptical (HR=1)`), 10,10,10,10))
                      ),
                      chains=2, iter=3000, seed=1233)

# 4.8 Pessimistic Prior
CVD_pess = fit.models(formula = form, data = main_use, distr = c("weibullPH"), method = "hmc",
                      priors = list(
                        weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                         mu_beta = c(0, mean(all_priors_calib$`Pessimistic (HR=1.1)`), 0,0,0,0), 
                                         sigma_beta= c(10, sd(all_priors_calib$`Pessimistic (HR=1.1)`), 10,10,10,10))
                      ),
                      chains=2, iter=3000, seed=1233)

# 4.9 Model for K-M curve (Intervention only)
cvd9km = fit.models(formula = Surv(t2mistrrevcvddthoth,mistrrevcvddthoth) ~ as.factor(interventionC),
                    data = main_use, distr = c("weibullPH"), method = "hmc",
                    priors = list(
                      weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                       mu_beta = c(0, mean(all_priors_calib$`SBP+CRP+FMD 50% shared`)), 
                                       sigma_beta= c(10, sd(all_priors_calib$`SBP+CRP+FMD 50% shared`)))
                    ),
                    chains=2, iter=3000, seed=1233)

# Group all CVD model fits into a list
CVD <- grep("CVD_", names(.GlobalEnv), value=TRUE)
CVD_list <- do.call("list", mget(CVD))


# 5. MAIN RESULTS (TABLE 3) -----------------------------------------------

# Helper function to extract Hazard Ratios and 95% Credible Intervals
comb = function(res){
  temp = res$models$`Weibull (PH)`
  n4 = temp@.MISC$summary$quan[2,5] # Median HR
  n5 = temp@.MISC$summary$quan[2,1] # Lower CI
  n6 = temp@.MISC$summary$quan[2,9] # Upper CI
  
  print(paste0('Weibull: ',
               round(exp(n4),2)," (",
               round(exp(n5),2),", ",
               round(exp(n6),2),")"))
}

# Print estimates for each model scenario
for (i in 1:length(CVD_list)){
  print(names(CVD_list)[[i]])
  comb(CVD_list[[i]])
}


# 6. RISK DIFFERENCE CALCULATION ------------------------------------------
# Calculates absolute risk difference (RD) at a specific evaluation time.

t_eval <- 365*4 # 4 years follow-up

RDcalc = function(cvdfit){
  temp = cvdfit$models$`Weibull (PH)`
  post <- rstan::extract(temp)
  alpha <- post$alpha
  beta <- post$beta
  
  # Evaluate covariates at their means (Adjust based on your actual data structure)
  # intercept, interventionC, agerand, men, interventionM, whi
  X_control <- c(1, 0, mean(main_use$agerand), mean(main_use$men), 0, 0) 
  X_treat   <- c(1, 1, mean(main_use$agerand), mean(main_use$men), 0, 0)
  
  # Calculate linear predictors
  lp_control <- as.vector(beta %*% X_control)
  lp_treat   <- as.vector(beta %*% X_treat)
  
  S_control <- exp(- (t_eval ^ alpha) * exp(lp_control))
  S_treat   <- exp(- (t_eval ^ alpha) * exp(lp_treat))
  
  risk_control <- 1 - S_control
  risk_treat <- 1 - S_treat
  RD <- risk_treat - risk_control
  
  return(c(mean(RD), quantile(RD, c(0.025, 0.975))))
}

# Apply to all models and build data frame
rd_results <- lapply(CVD_list, RDcalc)
rd_df <- do.call(rbind, rd_results) %>% as.data.frame()
colnames(rd_df) <- c("RD_mean", "RD_lower", "RD_upper")

rd_df$est = paste0(round(rd_df$RD_mean*100,2), " (", round(rd_df$RD_lower*100,2), ", ", round(rd_df$RD_upper*100,2), ")")
# write.csv(rd_df,"rd_df.csv")


# 7. POSTERIOR LIKELIHOOD --------------------------------------
# Evaluates the probability that HR is less than 1.0, 0.9, and 0.8 

pl = replicate(7, rep(NA, length(all_priors_calib))) %>% data.frame()

# Prior probabilities
for (i in 1:length(all_priors_calib)){
  pl[i,1] = names(all_priors_calib)[[i]]
  AA = all_priors_calib[[i]]
  vec = AA %>% exp()
  pl[i,2:4] = c(mean(vec<1), mean(vec<0.9), mean(vec<0.8))*100
}

# Posterior probabilities
desired_order <- c("CVD_ni","CVD_sbp","CVD_crp","CVD_sbpcrp","CVD_sbpcrp50",
                   "CVD_sbpcrpfmd","CVD_skep","CVD_pess")
CVD_list <- CVD_list[desired_order]

for (i in 1:length(CVD_list)){
  AA = CVD_list[[i]]
  vec = as.matrix(AA$models$`Weibull (PH)`)[,2] %>% exp()
  pl[i,5:7] = c(mean(vec<1), mean(vec<0.9), mean(vec<0.8), mean(vec<0.7))*100 # Adjusted mapping
}

colnames(pl)=c("Prior","pre1","pre0.9","pre0.8","post1","post0.9","post0.8")
pl$Prior[pl$Prior=="0 Non-informative"] = "Non-informative"

plnew_long2 <- pl %>%
  pivot_longer(
    cols = c("pre1", "pre0.9", "pre0.8", "post1", "post0.9", "post0.8"),
    names_to = "prob_type",
    values_to = "prob_val"
  ) %>%
  mutate(
    Group = ifelse(grepl("^pre", prob_type), "A: Prior probability", "B: Posterior probability"),
    Threshold = case_when(
      prob_type %in% c("pre1", "post1") ~ "HR<1.0",
      prob_type %in% c("pre0.9", "post0.9") ~ "HR<0.9",
      TRUE ~ "HR<0.8"
    )
  )

plnew_long2$scenario <- factor(plnew_long2$Prior, levels = rev(unique(plnew_long2$Prior)))
plnew_long2$Threshold = factor(plnew_long2$Threshold, levels = c("HR<1.0", "HR<0.9", "HR<0.8"))

# Generate Tile Plot
ggplot(plnew_long2, aes(x = Threshold, y = scenario, fill = prob_val)) +
  geom_tile() +
  geom_text(aes(label = ifelse(prob_val >99.9, ">99.9%", sprintf("%.1f", prob_val))), size = 3) +
  facet_wrap(~ Group) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(x = NULL, y = "Scenario", fill = "Probability") +
  theme_minimal() +
  theme(panel.grid = element_blank())


# 8. KAPLAN-MEIER PLOTS ----------------------------------------

main_use$interventionC = as.factor(main_use$interventionC)
km = survfit(Surv(t2mistrrevcvddthoth, mistrrevcvddthoth) ~ interventionC, data = main_use)

km_plot <- ggsurvplot(
  km,
  data = main_use,
  linetype = "interventionC",
  palette = c("blue", "blue"),
  censor = FALSE,
  ylim = c(0.90, 1),
  xlab = "Follow-up (days)",
  ylab = "Survival",
  legend.title = "Intervention Group",
  legend.labs = c("Placebo", "Active Cocoa"),
  ggtheme = theme_classic(),
  risk.table = FALSE
)

# Custom plot adjustments
km_plot$plot <- km_plot$plot +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_x_continuous(breaks = c(500, 1000, 1500)) +
  theme(
    legend.position = "bottom",         
    legend.key = element_blank(),       
    legend.title = element_text(size = 12), 
    legend.text = element_text(size = 10)    
  )

print(km_plot)


# 9. MODEL DIAGNOSTICS ----------------------------
# Evaluates MCMC chain mixing and autocorrelation for a representative model.

rstan::traceplot(CVD_sbpcrpfmd$models$`Weibull (PH)`)
bayesplot::mcmc_acf(as.matrix(CVD_sbpcrpfmd$models$`Weibull (PH)`))
--------------------------------------
# Evaluates the probability that HR is less than 1.0, 0.9, and 0.8 

pl = replicate(7, rep(NA, length(all_priors_calib))) %>% data.frame()

# Prior probabilities
for (i in 1:length(all_priors_calib)){
  pl[i,1] = names(all_priors_calib)[[i]]
  AA = all_priors_calib[[i]]
  vec = AA %>% exp()
  pl[i,2:4] = c(mean(vec<1), mean(vec<0.9), mean(vec<0.8))*100
}

# Posterior probabilities
desired_order <- c("CVD_ni","CVD_sbp","CVD_crp","CVD_sbpcrp","CVD_sbpcrp50",
                   "CVD_sbpcrpfmd","CVD_skep","CVD_pess")
CVD_list <- CVD_list[desired_order]

for (i in 1:length(CVD_list)){
  AA = CVD_list[[i]]
  vec = as.matrix(AA$models$`Weibull (PH)`)[,2] %>% exp()
  pl[i,5:7] = c(mean(vec<1), mean(vec<0.9), mean(vec<0.8), mean(vec<0.7))*100 # Adjusted mapping
}

colnames(pl)=c("Prior","pre1","pre0.9","pre0.8","post1","post0.9","post0.8")
pl$Prior[pl$Prior=="0 Non-informative"] = "Non-informative"

plnew_long2 <- pl %>%
  pivot_longer(
    cols = c("pre1", "pre0.9", "pre0.8", "post1", "post0.9", "post0.8"),
    names_to = "prob_type",
    values_to = "prob_val"
  ) %>%
  mutate(
    Group = ifelse(grepl("^pre", prob_type), "A: Prior probability", "B: Posterior probability"),
    Threshold = case_when(
      prob_type %in% c("pre1", "post1") ~ "HR<1.0",
      prob_type %in% c("pre0.9", "post0.9") ~ "HR<0.9",
      TRUE ~ "HR<0.8"
    )
  )

plnew_long2$scenario <- factor(plnew_long2$Prior, levels = rev(unique(plnew_long2$Prior)))
plnew_long2$Threshold = factor(plnew_long2$Threshold, levels = c("HR<1.0", "HR<0.9", "HR<0.8"))

# Generate Tile Plot
ggplot(plnew_long2, aes(x = Threshold, y = scenario, fill = prob_val)) +
  geom_tile() +
  geom_text(aes(label = ifelse(prob_val >99.9, ">99.9%", sprintf("%.1f", prob_val))), size = 3) +
  facet_wrap(~ Group) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(x = NULL, y = "Scenario", fill = "Probability") +
  theme_minimal() +
  theme(panel.grid = element_blank())


# 8. KAPLAN-MEIER PLOTS ----------------------------------------

main_use$interventionC = as.factor(main_use$interventionC)
km = survfit(Surv(t2mistrrevcvddthoth, mistrrevcvddthoth) ~ interventionC, data = main_use)

km_plot <- ggsurvplot(
  km,
  data = main_use,
  linetype = "interventionC",
  palette = c("blue", "blue"),
  censor = FALSE,
  ylim = c(0.90, 1),
  xlab = "Follow-up (days)",
  ylab = "Survival",
  legend.title = "Intervention Group",
  legend.labs = c("Placebo", "Active Cocoa"),
  ggtheme = theme_classic(),
  risk.table = FALSE
)

# Custom plot adjustments
km_plot$plot <- km_plot$plot +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_x_continuous(breaks = c(500, 1000, 1500)) +
  theme(
    legend.position = "bottom",         
    legend.key = element_blank(),       
    legend.title = element_text(size = 12), 
    legend.text = element_text(size = 10)    
  )

print(km_plot)


# 9. MODEL DIAGNOSTICS ----------------------------
# Evaluates MCMC chain mixing and autocorrelation for a representative model.

rstan::traceplot(CVD_sbpcrpfmd$models$`Weibull (PH)`)
bayesplot::mcmc_acf(as.matrix(CVD_sbpcrpfmd$models$`Weibull (PH)`))