### cntrl+find "MODIFY AS NEEDED"
seed_number <- 1
set.seed(seed_number)

library(cmdstanr) # used to be 'rstan' package
library(loo) # waic
library(bayesplot) # mcmc_trace
library(ggplot2)
library(dplyr)
library(tibble)
library(beepr)

thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript (shell script)
    return(dirname(normalizePath(sub(needle, "", cmdArgs[match]))))
    #return(getwd())
  } else if(Sys.getenv("RSTUDIO")=="1") {
    # Rstudio
    return(dirname(rstudioapi::getSourceEditorContext()$path))
    #script_name <- tools::file_path_sans_ext(basename(rstudioapi::getSourceEditorContext()$path))
  }
}

code_dir <- thisFile()
setwd(code_dir)

parent_source <- dirname(code_dir)
base_name <- basename(code_dir)
sub_base_length <- nchar(base_name)-nchar("code")
main_base_name <- substr(base_name,
                         1, nchar(base_name)-sub_base_length)
sub_base_name <- substr(base_name,
                        nchar(base_name)-sub_base_length+1, nchar(base_name))
output_dir <- paste0(parent_source,"/","outputs",sub_base_name)
data_dir <- paste0(parent_source,"/","data")

source(paste0("functions",sub_base_name,".R"))

### MODIFY AS NEEDED
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(ncores)) {
  ncores <- parallelly::availableCores()
}

### MODIFY AS NEEDED
stan_file_name <- paste0("DPTcorr_lessREs_TjointNC_IB",sub_base_name,".stan")
output_file_name <- paste0("DPTcorr_lessREs_TjointNC_IB",sub_base_name,".RDS")





# Data
dat <- read.csv(paste0(data_dir,"/","mkim_proj2.csv"))
n_rows_dat <- nrow(dat)
ID <- dat$ID
N <- length(unique(ID))
ages <- dat$TIME
min_age <- min(ages)
max_age <- max(ages)

obs_indxs_6MWD <- which(!is.na(dat[, "SMWT"]))
n_obs_6MWD <- length(obs_indxs_6MWD)
y_obs_6MWD <- dat[obs_indxs_6MWD, "SMWT"]
ages_6MWD <- ages[obs_indxs_6MWD]
ID_6MWD <- ID[obs_indxs_6MWD] # min, max, length(unique(...)) same as 'ID'

obs_indxs_BFLH <- which(!is.na(dat[, "BFLH"]))
n_obs_BFLH <- length(obs_indxs_BFLH)
y_obs_BFLH <- dat[obs_indxs_BFLH, "BFLH"]
min_obs_BFLH <- min(y_obs_BFLH)
max_obs_BFLH <- max(y_obs_BFLH)
ages_BFLH <- ages[obs_indxs_BFLH]
ID_BFLH <- ID[obs_indxs_BFLH] # min, max, length(unique(...)) same as 'ID'

obs_indxs_GRA <- which(!is.na(dat[, "GRA"]))
n_obs_GRA <- length(obs_indxs_GRA)
y_obs_GRA <- dat[obs_indxs_GRA, "GRA"]
min_obs_GRA <- min(y_obs_GRA)
max_obs_GRA <- max(y_obs_GRA)
ages_GRA <- ages[obs_indxs_GRA]
ID_GRA <- ID[obs_indxs_GRA] # min, max, length(unique(...)) same as 'ID'

obs_indxs_MG <- which(!is.na(dat[, "MG"]))
n_obs_MG <- length(obs_indxs_MG)
y_obs_MG <- dat[obs_indxs_MG, "MG"]
min_obs_MG <- min(y_obs_MG)
max_obs_MG <- max(y_obs_MG)
ages_MG <- ages[obs_indxs_MG]
ID_MG <- ID[obs_indxs_MG] # min, max, length(unique(...)) same as 'ID'

obs_indxs_SOL <- which(!is.na(dat[, "SOL"]))
n_obs_SOL <- length(obs_indxs_SOL)
y_obs_SOL <- dat[obs_indxs_SOL, "SOL"]
min_obs_SOL <- min(y_obs_SOL)
max_obs_SOL <- max(y_obs_SOL)
ages_SOL <- ages[obs_indxs_SOL]
ID_SOL <- ID[obs_indxs_SOL] # min, max, length(unique(...)) same as 'ID'

obs_indxs_VL <- which(!is.na(dat[, "VL"]))
n_obs_VL <- length(obs_indxs_VL)
y_obs_VL <- dat[obs_indxs_VL, "VL"]
min_obs_VL <- min(y_obs_VL)
max_obs_VL <- max(y_obs_VL)
ages_VL <- ages[obs_indxs_VL]
ID_VL <- ID[obs_indxs_VL] # min, max, length(unique(...)) same as 'ID'

stan_dat <- list(N=N, min_age=min_age, max_age=max_age,
                 
                 n_obs_6MWD=n_obs_6MWD,
                 y_obs_6MWD=y_obs_6MWD,
                 ages_6MWD=ages_6MWD, ID_6MWD=ID_6MWD,
                 
                 n_obs_BFLH=n_obs_BFLH,
                 min_obs_BFLH=min_obs_BFLH, max_obs_BFLH=max_obs_BFLH,
                 y_obs_BFLH=y_obs_BFLH,
                 ages_BFLH=ages_BFLH, ID_BFLH=ID_BFLH,
                 
                 n_obs_GRA=n_obs_GRA,
                 min_obs_GRA=min_obs_GRA, max_obs_GRA=max_obs_GRA,
                 y_obs_GRA=y_obs_GRA,
                 ages_GRA=ages_GRA, ID_GRA=ID_GRA,
                 
                 n_obs_MG=n_obs_MG,
                 min_obs_MG=min_obs_MG, max_obs_MG=max_obs_MG,
                 y_obs_MG=y_obs_MG,
                 ages_MG=ages_MG, ID_MG=ID_MG,
                 
                 n_obs_SOL=n_obs_SOL,
                 min_obs_SOL=min_obs_SOL, max_obs_SOL=max_obs_SOL,
                 y_obs_SOL=y_obs_SOL,
                 ages_SOL=ages_SOL, ID_SOL=ID_SOL,
                 
                 n_obs_VL=n_obs_VL,
                 min_obs_VL=min_obs_VL, max_obs_VL=max_obs_VL,
                 y_obs_VL=y_obs_VL,
                 ages_VL=ages_VL, ID_VL=ID_VL
)





# PAPER ESTIMATES + RSEs for F.E.'s
param_types <- c("DP", "DPT", "gamma", "S0")
FE_names_per_traj <- paste0("bar_", param_types)
traj_names <- c("6MWD", "BFLH", "GRA", "MG", "SOL", "VL")
n_FEs_per_traj <- length(FE_names_per_traj)
n_trajs <- length(traj_names)
# estimates + RSEs of F.E.'s in univariate models
univ_estimates <- matrix(data=c(NA, 34.62, 20.68, 18.82, 27.48, 35.55,
                                12.33, 9.91, 12.97, 12.03, 12.78, 10.22,
                                56.16, 10.8, 12.13, 10.84, 8.2, 8.53,
                                365.6, 41.84, 38.63, 40.47, 40.8, 40.27
),
nrow=n_FEs_per_traj, ncol=n_trajs, byrow=TRUE,
dimnames=list(FE_names_per_traj, traj_names))
univ_RSEs <- matrix(data=c(NA, 6.23, 35.5, 15.8, 12.6, 3.90,
                           2.53, 2.19, 8.18, 4.59, 4.59, 2.60,
                           12.4, 9.23, 31.0, 20.8, 13.5, 5.92,
                           1.88, 1.88, 1.61, 1.40, 1.50, 1.21
),
nrow=n_FEs_per_traj, ncol=n_trajs, byrow=TRUE,
dimnames=list(paste0(FE_names_per_traj, "(RSE)"), traj_names))
# estimates + RSEs of F.E.'s in pairwise joint models w/o covariates
# (1st column (6MWD) contains averages)
joint_no_covariates_estimates <- 
  matrix(data=c(NA, 40.07, 28.11, 29.71, 34.40, 36.69,
                12.34, 10.18, 14.32, 13.74, 13.72, 10.49,
                54.12, 7.94, 9.66, 6.31, 4.49, 8.42,
                365.76, 39.94, 38.45, 38.95, 38.12, 40.31
  ),
  nrow=n_FEs_per_traj, ncol=n_trajs, byrow=TRUE,
  dimnames=list(FE_names_per_traj, traj_names))
joint_no_covariates_RSEs <- 
  matrix(data=c(NA, 3.41, 10.60, 8.89, 6.32, 4.65,
                2.52, 2.63, 3.79, 5.25, 3.81, 2.79,
                11.34, 5.12, 13.00, 12.6, 11.80, 5.79,
                1.83, 1.58, 1.45, 1.59, 1.73, 1.54
  ),
  nrow=n_FEs_per_traj, ncol=n_trajs, byrow=TRUE,
  dimnames=list(paste0(FE_names_per_traj, "(RSE)"), traj_names))





### prior specification for F.E.'s
# "univariate", "joint_no_covariates", "joint_yes_covariates"
FE_priors <- FE_priors_func(
  N,
  ### MODIFY AS NEEDED:
  joint_no_covariates_estimates, joint_no_covariates_RSEs,
  
  FE_names_per_traj, traj_names,
  n_FEs_per_traj, n_trajs
)
(FE_Nprior_means <- FE_priors$mean_estimates)
(FE_Nprior_SDs <- FE_priors$SD_estimates)
(FE_logNprior_mus <- FE_priors$logN_mu_estimates)
(FE_logNprior_sigmas <- FE_priors$logN_sigma_estimates) 





nc <- 4
ni <- 10000
nw <- floor(ni/4)
noi <- ni-nw

### MODIFY AS NEEDED
max_td <- 10 # 10 is default for 'max_treedepth' (12 or 15 -- bigger search)
ad <- 0.8 # 0.8 is default for 'adapt_delta' (0.99 or 0.999 -- dec. div.)

stan_model <- cmdstan_model(file.path(code_dir, stan_file_name))

start_time <- Sys.time()
stan_fit <- stan_model$sample(
  data = stan_dat,
  seed = seed_number,
  chains = nc, parallel_chains = ncores,
  iter_warmup = nw, iter_sampling = ni - nw,
  adapt_delta = ad, max_treedepth = max_td
  # init = init_fun
)
end_time <- Sys.time()
stan_time <- round(as.numeric(difftime(end_time, start_time, units='hours')),2)

trial_name <- paste0("seed=", seed_number,
                     ", nc=", nc, ", ni=", ni, ", nw=", nw,
                     ", adapt_delta=", ad, ", max_treedepth=", max_td)
stan_diagnostics <- stan_fit$sampler_diagnostics()
divergences_total <- divergences_total <- sum(stan_diagnostics[, , "divergent__"])
treedepth <- stan_diagnostics[, , "treedepth__"]
hits <- sum(treedepth == max_td)

stan_time_and_more <- 
  paste0(stan_file_name, " (", trial_name, "): ", stan_time, " hours, ",
         divergences_total, " divergences, ",
         hits, " transitions hit max_treedepth=", max_td)
setwd(output_dir)
write(stan_time_and_more, file="stan_time.txt", append=TRUE)
write("", file="stan_time.txt", append=TRUE)
setwd(code_dir)
beep(3)

stan_fit$save_object(file = file.path(output_dir, output_file_name))

stan_posterior <- stan_fit$draws(format = "df")
stan_posterior_transformed <- stan_posterior %>% mutate(
  # MODIFY AS NEEDED (ADD/DELETE FIXED EFFECTS IF NEEDED)
  bar_DPT_6MWD = exp(log_bar_DPT_6MWD),
  bar_gamma_6MWD = exp(log_bar_gamma_6MWD),
  bar_S0_6MWD = exp(log_bar_S0_6MWD),
  
  bar_DP_BFLH = exp(log_bar_DP_BFLH),
  bar_DPT_BFLH = exp(log_bar_DPT_BFLH),
  bar_gamma_BFLH = exp(log_bar_gamma_BFLH),
  bar_S0_BFLH = exp(log_bar_S0_BFLH),
  
  bar_DP_GRA = exp(log_bar_DP_GRA),
  bar_DPT_GRA = exp(log_bar_DPT_GRA),
  bar_gamma_GRA = exp(log_bar_gamma_GRA),
  bar_S0_GRA = exp(log_bar_S0_GRA),
  
  bar_DP_MG = exp(log_bar_DP_MG),
  bar_DPT_MG = exp(log_bar_DPT_MG),
  bar_gamma_MG = exp(log_bar_gamma_MG),
  bar_S0_MG = exp(log_bar_S0_MG),
  
  bar_DP_SOL = exp(log_bar_DP_SOL),
  bar_DPT_SOL = exp(log_bar_DPT_SOL),
  bar_gamma_SOL = exp(log_bar_gamma_SOL),
  bar_S0_SOL = exp(log_bar_S0_SOL),
  
  bar_DP_VL = exp(log_bar_DP_VL),
  bar_DPT_VL = exp(log_bar_DPT_VL),
  bar_gamma_VL = exp(log_bar_gamma_VL),
  bar_S0_VL = exp(log_bar_S0_VL)
)





# WAIC
log_lik_variables <- c("log_lik_6MWD", "log_lik_BFLH", "log_lik_GRA", 
                       "log_lik_MG", "log_lik_SOL", "log_lik_VL")
log_lik_matrix <- stan_fit$draws(variables = log_lik_variables,
                                 format="matrix")
waic_result <- waic(log_lik_matrix)
chain_id_vector <- rep(1:nc, each = noi)
r_eff <- relative_eff(exp(log_lik_matrix), chain_id = chain_id_vector)
loo_result <- loo(log_lik_matrix, r_eff = r_eff,
                  cores = ncores)

setwd(output_dir)
sink("waic_loo.txt")
print("--- WAIC ---")
print(waic_result)
cat("\n\n--- LOO ---\n")
print(loo_result)
sink()

saveRDS(list(waic_result=waic_result,
             loo_result=loo_result),
        file="waic_loo.RDS")
setwd(code_dir)





### MODIFY AS NEEDED
eta <- 1  # from LKJ(eta)
s_prior_sd <- 1
omega_prior_sd <- 1
corr_RE_names <- # MAKE SURE TO LIST THIS IN ORDER OF s[1]:s[JA]!
  c("DPT_6MWD", "DPT_BFLH", "DPT_GRA", "DPT_MG", "DPT_SOL", "DPT_VL")

JA <- length(corr_RE_names)
s_strings <- paste0("s[", 1:JA, "]")
R_pairs <- combn(1:JA, 2)
R_strings <- apply(R_pairs, 2, function(x) paste0("R[", x[1], ",", x[2], "]"))

all_column_names <- colnames(stan_posterior)
omega_strings <- all_column_names[grep("^omega_", all_column_names)]
uncorr_RE_names <- sub("^omega_", "", omega_strings)
JB <- length(uncorr_RE_names)
# I have 'eta_strings' here so it has same ordering as 'omega_strings'
eta_strings <- paste0("eta_", uncorr_RE_names)

# check
all_expected_re_params <- c(s_strings, R_strings, omega_strings)
(is_valid <- all(all_expected_re_params %in% all_column_names))
stopifnot(
  "Model Mismatch: The parameter names generated in R (s_strings, R_strings, omega_strings) do not all exist as columns in the Stan fit output. Check `corr_RE_names` for typos or check for changes in the Stan file." = is_valid
)
cat("Success: R script RE definitions are consistent with the fitted Stan model.\n")

omitted_FE_names <- c()





# POSTERIOR ANALYSIS NOTES
# 1. posterior median of parameter shouldn't be approx. 0
# 1b. posterior CRI shouldn't contain 0
# 2. posterior S.D. of parameter shouldn't be very very low (or very very high?),
#    but relatively small is good (<0.1 is good)
# 3. high RSE (particularly > 40%) is bad
# 4. high shrinkage (particularly > 70%) is bad
#    (shrinkage issues can persist if omega is too flexible (weak prior),
#    or data is weak, or random effect varies little across subjects)
# 5. ess_bulk =  ESS for estimating posterior center
# 6. ess_tail = ESS for estimating posterior quantiles


FE_log_strings <- all_column_names[grep("^log_bar_", all_column_names)]
FE_log_draws <- stan_fit$draws(variables = FE_log_strings, format = "df")
# Apply exp() to all columns except metadata
FE_draws <- FE_log_draws %>%
  mutate(across(-c(.chain, .iteration, .draw), exp))
# Remove metadata columns
FE_summary <- FE_draws %>%
  select(-c(.chain, .iteration, .draw)) %>%
  summarise(across(everything(),
                   list(median = median,
                        mean = mean,
                        sd = sd,
                        q5 = ~quantile(., 0.05),
                        q95 = ~quantile(., 0.95)),
                   .names = "{.col}::{.fn}")) %>%
  tidyr::pivot_longer(everything(),
                      names_to = c("variable", ".value"),
                      names_sep = "::")

FE_cvg_dx <- stan_fit$summary(variables = FE_log_strings) %>%
  select(variable, rhat, ess_bulk, ess_tail)

FE_dx <- FE_summary %>%
  # both 'FE_summary' and 'FE_cvg_dx' have "log_bar_" prefix under 'variable'
  left_join(FE_cvg_dx, by = "variable") %>%
  # remove "log_bar_" prefix
  mutate(variable = sub("^log_bar_", "", variable)) %>%
  select(variable, rhat, ess_bulk, ess_tail, median, mean, sd, q5, q95) %>%
  rename(
    CRI_lower = q5,
    CRI_upper = q95
  ) %>%
  mutate(
    # same as rounding to 3 decimal places
    RSE = round((sd / abs(mean))* 100, 1),
    across(
      .cols = c(median, mean, sd, CRI_lower, CRI_upper),
      .fns = ~ round(.x, 3)
    ),
    across(
      .cols = c(ess_bulk, ess_tail),
      .fns = ~ round(.x, 0)
    )
  ) %>%
  column_to_rownames("variable")

setwd(output_dir)
write.csv(FE_dx, "FE_dx.csv", row.names=TRUE)
setwd(code_dir)



cor_summary <- stan_fit$summary(variables = R_strings)
cor_dx <- cor_summary %>%
  select(variable, rhat, ess_bulk, ess_tail, median, mean, sd, q5, q95) %>%
  rename(
    CRI_lower = q5,
    CRI_upper = q95
  ) %>%
  mutate(
    # same as rounding to 3 decimal places
    RSE = round((sd / abs(mean))* 100, 1),
    across(
      .cols = c(median, mean, sd, CRI_lower, CRI_upper),
      .fns = ~ round(.x, 3)
    ),
    across(
      .cols = c(ess_bulk, ess_tail),
      .fns = ~ round(.x, 0)
    )
  ) %>%
  column_to_rownames("variable")
setwd(output_dir) 
write.csv(cor_dx, "cor_dx.csv", row.names = TRUE)
setwd(code_dir)

# If posterior matches prior, parameter is not informed by data
# If posterior is tighter than prior, parameter is informed by data
setwd(output_dir)
dir.create("cor_plots")
setwd(paste0(output_dir,"/cor_plots"))
for (corr_indx in 1:length(R_strings)) {
  param_name <- R_strings[corr_indx]
  file_name <- paste0(param_name, ".pdf")
  
  pdf(file_name, width = 6, height = 5)
  
  hist(stan_posterior[[param_name]],
       breaks = 30,
       main = paste0("Posterior of ", param_name, " vs. Prior=LKJ(",eta,")"),
       xlab = param_name,
       probability = TRUE)
  curve(dbeta((x + 1) / 2, eta, eta) / 2,
        from = -1, to = 1, col = "red", lwd = 2, add = TRUE)
  mtext("(Posterior samples should be narrower than (or shifted from) the prior)", line = 0.5)
  dev.off()
}
setwd(code_dir)



s_summary <- stan_fit$summary(variables = s_strings)
shrinkage_list_s <- list()
for (ja in 1:JA) {
  s_string <- s_strings[ja]
  eta_means <- apply(
    stan_posterior[, paste0("eta_mat[", ja, ",", 1:N, "]")],
    2, mean)
  s_mean <- s_summary %>%
    filter(variable == s_string) %>%
    pull(mean)
  shrinkage <- 1 - sd(eta_means)/s_mean
  shrinkage_list_s[[s_string]] <- tibble(
    variable = s_string,
    shrinkage = shrinkage
  )
}
shrinkage_df_s <- bind_rows(shrinkage_list_s)
s_dx <- s_summary %>%
  # Safely join the shrinkage data. 'left_join' ensures rows are matched by 'variable' name.
  left_join(shrinkage_df_s, by = "variable") %>%
  select(variable, rhat, ess_bulk, ess_tail, median, mean, sd, q5, q95, shrinkage) %>%
  rename(
    CRI_lower = q5,
    CRI_upper = q95
  ) %>%
  mutate(
    # same as rounding to 3 decimal places
    RSE = round((sd / abs(mean))* 100, 1),
    across(
      .cols = c(median, mean, sd, CRI_lower, CRI_upper, shrinkage),
      .fns = ~ round(.x, 3)
    ),
    across(
      .cols = c(ess_bulk, ess_tail),
      .fns = ~ round(.x, 0)
    )
  ) %>%
  column_to_rownames("variable")
setwd(output_dir)
write.csv(s_dx, "s_dx.csv", row.names = TRUE)
setwd(code_dir)

# If posterior matches prior, parameter is not informed by data
# If posterior is tighter than prior, parameter is informed by data
setwd(output_dir)
dir.create("s_plots")
setwd(paste0(output_dir,"/s_plots"))
for (ja in 1:JA) {
  param_name <- s_strings[ja]
  file_name <- paste0(param_name, ".pdf")
  
  pdf(file_name, width = 6, height = 5)
  
  hist(stan_posterior[[param_name]],
       breaks = 30,
       main = paste0("Posterior of ", param_name, " vs. Prior=N+(0,",s_prior_sd,")"),
       xlab = param_name,
       probability = TRUE,
       xlim = c(0, max(stan_posterior[[param_name]]) * 1.2))
  curve(2 * dnorm(x, 0, s_prior_sd),
        from = 0, col = "red", lwd = 2, add = TRUE)
  mtext("(Posterior samples should be narrower than (or shifted from) the prior)", line = 0.5)
  dev.off()
}
setwd(code_dir)



# Shrinkage for subject R.E.'s
omega_summary <- stan_fit$summary(variables = omega_strings)
shrinkage_list_omega <- list()
for (jb in 1:JB) {
  omega_string <- omega_strings[jb]
  eta_means <- apply(
    stan_posterior[, paste0(eta_strings[jb],"[", 1:N, "]")],
    2, mean)
  omega_mean <- omega_summary %>%
    filter(variable == omega_string) %>%
    pull(mean)
  shrinkage <- 1 - sd(eta_means)/omega_mean
  shrinkage_list_omega[[omega_string]] <- tibble(
    variable = omega_string, 
    shrinkage = shrinkage
  )
}
shrinkage_df_omega <- bind_rows(shrinkage_list_omega)
omega_dx <- omega_summary %>%
  # Safely join the shrinkage data. 'left_join' ensures rows are matched by 'variable' name.
  left_join(shrinkage_df_omega, by = "variable") %>%
  select(variable, rhat, ess_bulk, ess_tail, median, mean, sd, q5, q95, shrinkage) %>%
  rename(
    CRI_lower = q5,
    CRI_upper = q95
  ) %>%
  mutate(
    # same as rounding to 3 decimal places
    RSE = round((sd / abs(mean))* 100, 1),
    across(
      .cols = c(median, mean, sd, CRI_lower, CRI_upper, shrinkage),
      .fns = ~ round(.x, 3)
    ),
    across(
      .cols = c(ess_bulk, ess_tail),
      .fns = ~ round(.x, 0)
    )
  ) %>%
  column_to_rownames("variable")
setwd(output_dir)
write.csv(omega_dx, "omega_dx.csv", row.names = TRUE)
setwd(code_dir)

# If posterior matches prior, parameter is not informed by data
# If posterior is tighter than prior, parameter is informed by data
setwd(output_dir)
base_name_omega1 <- "omega_plots"
#base_name_omega2 <- "N_meds_eta_plots(shrinkage)"
base_name_omega2 <- "N_means_eta_plots(shrinkage)"
base_name_omega3 <- "N_SDs_eta_plots"
dir.create(base_name_omega1)
dir.create(base_name_omega2)
dir.create(base_name_omega3)
for (jb in 1:JB) {
  param_name <- omega_strings[jb]
  
  omega_mean <- omega_summary %>%
    filter(variable == param_name) %>%
    pull(mean)
  
  file_name <- paste0(base_name_omega1, "/", param_name, ".pdf")
  pdf(file_name, width = 6, height = 5)
  hist(stan_posterior[[param_name]],
       breaks = 30,
       main = paste0("Posterior of ", param_name, " vs. Prior=N+(0,",omega_prior_sd,")"),
       sub = "(If posterior is hugging 0, then data isn't identifying subject variation)",
       xlab = param_name,
       probability = TRUE,
       xlim = c(0, max(stan_posterior[[param_name]]) * 1.2))
  curve(2 * dnorm(x, 0, omega_prior_sd),
        from = 0, col = "red", lwd = 2, add = TRUE)
  mtext("(Posterior samples should be narrower than (or shifted from) the prior)", line = 0.5)
  dev.off()
  
  # ideally peaks around 0.3-0.6; <0 or >1 means some identifiability issue
  # clustering around 0 means close to group mean => high shrinkage
  eta_means <- apply(stan_posterior[, paste0(eta_strings[jb],"[", 1:N, "]")], 2, mean)
  file_name <- paste0(base_name_omega2, "/", eta_strings[jb], ".pdf")
  pdf(file_name, width = 6, height = 5)
  hist(eta_means,
       freq = FALSE,
       breaks = 30,
       main = paste0("N post. means of subject R.E.'s vs. N(0, omega_mean)"),
       sub="(3. wide=underest. omega or overdispersed eta_{i}'s?)",
       xlab = paste0("Posterior means of ",eta_strings[jb],"[i]'s"))
  curve(dnorm(x, mean = 0, sd = omega_mean),
        col = "red", lwd = 2, add = TRUE)
  mtext("(Hist: 1. match=good, 2. narrow=high shrinkage?)", line = 0.5)
  dev.off()
  
  eta_sds <- apply(stan_posterior[, paste0(eta_strings[jb],"[", 1:N, "]")], 2, sd)
  file_name <- paste0(base_name_omega3, "/", eta_strings[jb], ".pdf")
  pdf(file_name, width = 6, height = 5)
  hist(eta_sds, 
       breaks = 30,
       main = paste0("N post. S.D.'s of subject R.E.'s"),
       sub = paste0("(hopefully most S.D.'s are small and similar)"),
       xlab = paste0("Posterior S.D.'s of ",eta_strings[jb],"[i]'s"))
  mtext("(heteregeneous or large S.D.'s may imply weak identifiability of R.E.'s)", line = 0.5)
  dev.off()
}
setwd(code_dir)



# TRACEPLOTS
theme_set(theme_bw())
setwd(output_dir)
dir.create("traceplots")
setwd(paste0(output_dir,"/traceplots"))

trace_regex_pars_vec <- 
  c("^tau_", # residual variance 
    "^log_bar_DP_", "^log_bar_DPT_", "^log_bar_gamma_", "^log_bar_S0_", # F.E.'s
    "^omega_DP_", "^omega_gamma_", "^omega_S0_", # group S.D.'s of non-correlated R.E.'s
    "^s\\[" # group S.D.'s of correlated R.E.'s
    #"^R\\[" # group corr's of correlated R.E.'s
  )
trace_filenames_vec <- 
  paste0(c(
    "tau_trace",
    "log_bar_DP_trace", "log_bar_DPT_trace", "log_bar_gamma_trace", "log_bar_S0_trace",
    "omega_DP_trace", "omega_gamma_trace", "omega_S0_trace",
    "s_trace"
    #"R_trace"
  ),".pdf")

trace_plot_width=8
trace_plot_height=6
for (trace_indx in 1:length(trace_filenames_vec)) {
  print(paste0("plotting ",trace_filenames_vec[trace_indx]," ..."))
  
  trace_plot_temp <-
    mcmc_trace(stan_posterior,
               regex_pars = trace_regex_pars_vec[trace_indx])
  if (trace_filenames_vec[trace_indx]=="R_trace.pdf") {
    trace_plot_width=10
    trace_plot_height=8
  }
  ggsave(filename = trace_filenames_vec[trace_indx],
         plot = trace_plot_temp,
         width=trace_plot_width, height=trace_plot_height)
}
setwd(code_dir)



# PAIRS PLOTS
# NOTE: red x's are divergences (where sampler crashed)
# NOTE: correlation with lp__ (oval plot) means bad sampling
setwd(output_dir)
dir.create("pairs_plots")
setwd(paste0(output_dir,"/pairs_plots"))

np <- nuts_params(stan_fit)

# All F.E.'s within 6MWD or within a muscle group
# NOTE: check for oval plots (collinearity between two params)
# -- expected for our model b/c trajectory parameters
# naturally related to each other by defn.
pairs_regex_pars_vec <- 
  paste0("^log_bar_.*_", traj_names )
pairs_filenames_vec <- 
  paste0("pairs_log(fixed)_",traj_names,".pdf")

for (pplot_indx in 1:length(pairs_filenames_vec)) {
  print(paste0("plotting ",pairs_filenames_vec[pplot_indx]," ..."))
  
  pplot_temp <- mcmc_pairs(
    stan_posterior,
    regex_pars = pairs_regex_pars_vec[pplot_indx],
    pars = "lp__",
    np=np
  )
  ggsave(pairs_filenames_vec[pplot_indx], pplot_temp,
         width=10, height=8)
}

# F.E.'s vs. (corresponding) correlated R.E.'s
# NOTE: check for funnels btwn F.E. and S.D. of R.E.'s,
# which means poor identification (hard to tell apart)
for (ja in 1:JA) {
  FE_string_temp <- paste0("log_bar_",corr_RE_names[ja])
  RE_string_temp <- s_strings[ja]
  filename_temp <- paste0("pairs-",
                          FE_string_temp,"-",
                          RE_string_temp,".pdf")
  
  print(paste0("plotting ",filename_temp," ..."))
  
  pplot_temp <- mcmc_pairs(
    stan_posterior,
    pars = c(FE_string_temp,
             RE_string_temp,
             "lp__"),
    np = np
  )
  ggsave(filename_temp,
         pplot_temp,
         width = 8, height = 6)
}

# F.E.'s vs. (corresponding) non-correlated R.E.'s
# NOTE: check for funnels btwn F.E. and S.D. of R.E.'s,
# which means poor identification (hard to tell apart)
for (jb in 1:JB) {
  FE_string_temp <- paste0("log_bar_",uncorr_RE_names[jb])
  RE_string_temp <- omega_strings[jb]
  filename_temp <- paste0("pairs-",
                          FE_string_temp,"-",
                          RE_string_temp,".pdf")
  print(paste0("plotting ",filename_temp," ..."))
  
  pplot_temp <- mcmc_pairs(
    stan_posterior,
    pars = c(FE_string_temp,
             RE_string_temp,
             "lp__"),
    np = np
  )
  ggsave(filename_temp,
         pplot_temp,
         width = 8, height = 6)
}
setwd(code_dir)



beep(8)