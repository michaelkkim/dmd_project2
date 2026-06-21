### cntrl+find "MODIFY AS NEEDED"
set.seed(1)

library(rstan)
library(loo) # waic
library(Ryacas0)
library(Ryacas)
library(readr)

library(foreach)
library(doParallel)
library(beepr)

library(zoo) # zoo
library(compositions) # acomp, clr, clrInv

library(RColorBrewer) # setting colors
library(ggplot2) # plotting
library(dplyr) # data filtering

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





dat <- read.csv(paste0(data_dir,"/","mkim_proj2.csv"))
n_rows_dat <- nrow(dat)
ID <- dat$ID
unique_ID <- unique(ID) # unique_ID is same aw index! (convenient for us)
N <- length(unique_ID)
muscle_names <- c("BFLH", "GRA", "MG", "SOL", "VL")
L <- length(muscle_names)
muscle_colors <- c("brown", "red", "magenta", "blue", "cyan")

# missing data summary (multiply by 'n_rows_dat' to get counts):
missing_perc1_func(dat[,c("SMWT",muscle_names)])
# distribution of subject-wise measurement times
Js <- as.numeric(table(dat[,"ID"]))
Js_dist_func(Js, N)

ages_dat <- dat$TIME
min_age_dat <- min(ages_dat)
max_age_dat <- max(ages_dat)

SMWT_dat <- dat$SMWT
min_SMWT_dat <- min(SMWT_dat, na.rm=TRUE) # 10
max_SMWT_dat <- max(SMWT_dat, na.rm=TRUE) # 543



# Yoon et al.:
# - "Firstly, the models cannot be extrapolated to predict the
# trajectories of 6MWD and MRI-T2 for individuals under the age of 5,
# because the trajectory of 6MWD should increase during the muscle
# development period, typically between the ages of 2 to 5."
# - "The inclusion criteria of the study are as follow; ambulatory and
# non-ambulatory males (ages 5-30 years ate baseline testing) ..."
min_age_search <- 5
max_age_search <- 30










### MODIFY AS NEEDED
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(ncores) || ncores < 1) {
  ncores <- parallelly::availableCores() - 1
}
### MODIFY AS NEEDED
output_file_name <- paste0("DPTcorr_lessREs_TjointNC_IB",sub_base_name,".RDS")
joint_model_stan_results <- 
  readRDS(paste0(output_dir,"/",output_file_name))



### Combine chains of STAN posterior?
joint_model_draws_3Darray <- joint_model_stan_results$draws(format = "draws_array")
### MODIFY AS NEEDED (chains)
#chain1_matrix <- joint_model_draws_3Darray[, 1, , drop = TRUE]
#joint_model_draws <- chain1_matrix
joint_model_draws <- joint_model_stan_results$draws(format = "df")
joint_model_mat <- as.matrix(joint_model_draws)
joint_model_colnames <- colnames(joint_model_mat)
post_means <- colMeans(joint_model_mat)
noi <- nrow(joint_model_mat)
noi_target <- 500 ### MODIFY AS NEEDED
thin_interval <- floor(noi / noi_target)
iters <- seq(from = 1, to = noi, by = thin_interval)
noi_target <- length(iters) # should be same as above 'noi_target'

min_SMWT_grid <- 1
max_SMWT_grid <- 601
grid_tick <- 1 ### MODIFY AS NEEDED
SMWT_grid <- seq(min_SMWT_grid, max_SMWT_grid, grid_tick)
grid_size <- length(SMWT_grid)

year_diff <- 1



### MODIFY AS NEEDED
grid_val_target <- 396 # change to where there is a sudden jump or drop
# useful when 'grid_tick' changes
target_grid_indx <- 
  which.min(abs(SMWT_grid - grid_val_target))
#mean_bounded_ages_old <- ### MODIFY AS NEEDED
#  read.csv("/Users/michaelkim/Desktop/Research/weights_6MWD/outputs1/J1_seed1_trials/1grid_7.5kthin/mean_bounded_ages.csv")
subjects_of_interest <- unique_ID ### MODIFY AS NEEDED
#subjects_of_interest <- as.numeric(which(!is.na( ### MODIFY AS NEEDED
#  mean_bounded_ages_old[, paste0("SMWT_",grid_val_target)])))
N_subset <- length(subjects_of_interest)
subjects_dropped <- unique_ID[-subjects_of_interest]
N_dropped <- length(subjects_dropped)



if (N_dropped > 0) {
  blue_palette <- colorRampPalette(c("lightblue", "navy"))(N_subset)
  red_palette <- colorRampPalette(c("mistyrose", "darkred"))(N_dropped)
  subject_colors <- rep(NA, N)
  subject_colors[subjects_of_interest] <- blue_palette
  subject_colors[subjects_dropped] <- red_palette
} else { # N_dropped ==0
  blue_palette <- colorRampPalette(c("lightblue", "navy"))(N_subset)
  subject_colors <- blue_palette
}
#dat$cohort_color <- subject_colors[dat$ID]



### Parameter matrices
log_S0_6MWD_mat <- joint_model_mat[, paste0("log_S0_6MWD[", 1:N, "]")]
log_DPT_6MWD_mat <- joint_model_mat[, paste0("log_DPT_6MWD[", 1:N, "]")]
log_gamma_6MWD_mat <- matrix(rep(joint_model_mat[, "log_bar_gamma_6MWD"], N), ncol=N)
muscle_param_list <- list()
param_prefixes <- c("log_S0", "log_DP", "log_DPT", "log_gamma")
for (m in muscle_names) {
  muscle_param_list[[m]] <- list()
  
  for (p_prefix in param_prefixes) {
    muscle_mat_temp <- matrix(NA, nrow=noi, ncol=N)
    
    FE_string <- paste0(sub("_", "_bar_", p_prefix), "_", m)
    if (!FE_string %in% joint_model_colnames) {
      stop("Could not find required FE parameter: ", FE_string)
    }
    
    for (i in 1:N) {
      RE_string <- paste0(p_prefix, "_", m, "[", i, "]")
      if (RE_string %in% joint_model_colnames) {
        muscle_mat_temp[, i] <- joint_model_mat[, RE_string]
      } else {
        muscle_mat_temp[, i] <- joint_model_mat[, FE_string]
      }
    }
    
    muscle_param_list[[m]][[p_prefix]] <- muscle_mat_temp
  }
}










### Plot ages vs 6MWD (data)
setwd(output_dir)
png(file=paste0("6MWD_traj_dat.png"),
    width=6.5, height=3.5, units="in", res=800, pointsize=6.9)
par(mfrow=c(1,1))
plot(ages_dat, SMWT_dat,
     type = "n", # no points yet
     xlab = "Age",
     ylab = "6MWD",
     main = "6MWD vs. Age (Data)",
     xlim = range(ages_dat, na.rm=TRUE),
     ylim = range(SMWT_dat, na.rm=TRUE)
)
for (i in 1:N) {
  subject_data <- subset(dat, ID == i)
  lines(subject_data$TIME,
        subject_data$SMWT,
        col = subject_colors[i],
        lwd = 1)
}
if (N_dropped > 0) {
  legend("topright", 
         legend = c(paste0("Slower Progressors (N=", N_subset, ")"), 
                    paste0("Faster Progressors (N=", N_dropped, ")")),
         col = c("navy", "darkred"),
         lty = 1, lwd = c(1, 1), bty = "n")
}
dev.off()
setwd(code_dir)



### Plot ages vs 6MWD (model) with credible bands
setwd(output_dir)
png(file=paste0("6MWD_pop_traj.png"),
    width=6.5, height=3.5, units="in", res=800, pointsize=6.9)
age_grid <- seq(min_age_search, max_age_search, length.out=200)
plot(1, type="n",
     xlim = range(age_grid),
     ylim = c(0, max(dat$SMWT, na.rm=TRUE) * 1.1),
     xlab = "Age",
     ylab = "6MWD",
     main = "6MWD vs. Age (Model)"
)
# (make) posterior mean population trajectory + pointwise 95% posterior credible bands
draws_6MWD_pop <- t(apply(joint_model_mat, 1, function(x) {
  get_SMWT_trajectory(age_grid,
                      exp(x["log_bar_S0_6MWD"]),
                      x["log_bar_DPT_6MWD"],
                      exp(x["log_bar_gamma_6MWD"]))
}))

predicted_6MWD_pop <- colMeans(draws_6MWD_pop, na.rm=TRUE)
band_6MWD_pop <- apply(draws_6MWD_pop, 2, quantile,
                       probs = c(0.025, 0.975), na.rm=TRUE)

# individual trajectories
#for (i in 1:N) {
  #S0_6MWD_i <- exp(mean(log_S0_6MWD_mat[, i]))
  #log_DPT_6MWD_i <- mean(log_DPT_6MWD_mat[, i])
  #gamma_6MWD_i <- exp(mean(log_gamma_6MWD_mat[, i]))
  
  #predicted_6MWD_i <- 
    #get_SMWT_trajectory(age_grid,
                        #S0_6MWD_i, log_DPT_6MWD_i, gamma_6MWD_i)
  
  #line_col <- adjustcolor(subject_colors[i], alpha.f = 0.15) 
  
  #lines(age_grid,
        #predicted_6MWD_i,
        #col = line_col,
        #lwd=1)
#}

# (plot) posterior mean population trajectory + pointwise 95% posterior credible bands
polygon(c(age_grid, rev(age_grid)),
        c(band_6MWD_pop[1, ], rev(band_6MWD_pop[2, ])),
        col = adjustcolor("grey", alpha.f = 0.35),
        border = NA)
lines(age_grid, predicted_6MWD_pop, col="black", lwd=1.5)

if (N_dropped > 0) {
  legend("topright",
         legend = c("Population Mean", "95% Credible Band",
                    paste0("Slower Progressors (N=", N_subset, ")"), 
                    paste0("Faster Progressors (N=", N_dropped, ")")),
         col = c("black", adjustcolor("gray", alpha.f = 0.35), "navy", "darkred"),
         lty = c(1, NA, 1, 1),
         lwd = c(1.5, NA, 1, 1),
         pch = c(NA, 15, NA, NA),
         pt.cex = c(NA, 2, NA, NA),
         bty = "n")
} else {
  legend("topright", 
         legend = c("Population Mean", "95% Credible Band"),
         col = c("black", adjustcolor("gray", alpha.f = 0.35)),
         lty = c(1, NA),
         lwd = c(1.5, NA),
         pch = c(NA, 15),
         pt.cex = c(NA, 2),
         bty = "n")
}

dev.off()
setwd(code_dir)



### Histogram of 6MWD (data)
setwd(output_dir)
png(file=paste0("6MWD_hist_dat.png"),
    width=6.5, height=3.5, units="in", res=800, pointsize=6.9)
par(mfrow=c(1,1))
hist(SMWT_dat,
     main=paste0("Histogram of SMWT_dat (n_rows_dat=",n_rows_dat," observations)"))
dev.off()
setwd(code_dir)



### Histogram of max_6MWD (data)
max_6MWD_dat <- sapply(1:N, function(i) 
  max(dat[which(ID==i), "SMWT"], na.rm=TRUE))
setwd(output_dir)
png(file=paste0("max6MWD_hist_dat.png"),
    width=6.5, height=3.5, units="in", res=800, pointsize=6.9)
par(mfrow=c(1,1))
hist(max_6MWD_dat,
     main=paste0("Histogram of max_6MWD_dat (N=",N," subjects)"))
dev.off()
setwd(code_dir)



### Boxplots of max_6MWD by group (data)
if (N_dropped > 0) {
  max_6MWD_dat_soi <- max_6MWD_dat[subjects_of_interest]
  max_6MWD_dat_dropped <- max_6MWD_dat[subjects_dropped]
  max_6MWD_dat_list <- list("Slower" = max_6MWD_dat_soi,
                            "Faster" = max_6MWD_dat_dropped)
  setwd(output_dir)
  png(file=paste0("max_6MWD_box_dat.png"),
      width=6.5, height=3.5, units="in", res=800, pointsize=6.9)
  boxplot(max_6MWD_dat_list,
          main = "Maximum 6MWD by Group (Data)",
          ylab = "Max_6MWD",
          xlab = "Group",
          col = c("deepskyblue3", "indianred2"))
  dev.off()
  setwd(code_dir) 
}


### Boxplots of max_6MWD by group (model)
if (N_dropped > 0) {
  setwd(output_dir)
  max_6MWD_model <- sapply(unique_ID, function(i) {
    S0_6MWD_i <- exp(mean(log_S0_6MWD_mat[, i]))
    log_DPT_6MWD_i <- mean(log_DPT_6MWD_mat[, i])
    gamma_6MWD_i <- exp(mean(log_gamma_6MWD_mat[, i]))
    
    # Assumes your "get_SMWT_trajectory" function exists from previous script parts
    predicted_6MWD_i <- get_SMWT_trajectory(age_grid,
                                            S0_6MWD_i, log_DPT_6MWD_i, gamma_6MWD_i)
    return(max(predicted_6MWD_i))
  })
  
  max_6MWD_model_soi <- max_6MWD_model[subjects_of_interest]
  max_6MWD_model_dropped <- max_6MWD_model[subjects_dropped]
  max_6MWD_model_list <- list("Slower" = max_6MWD_model_soi,
                              "Faster" = max_6MWD_model_dropped)
  png(file=paste0("max_6MWD_box_model.png"),
      width=6.5, height=3.5, units="in", res=800, pointsize=6.9)
  boxplot(max_6MWD_model_list,
          main = "Maximum 6MWD by Group (Model)",
          ylab = "Max_6MWD",
          xlab = "Group",
          col = c("deepskyblue3", "indianred2"))
  dev.off()
  setwd(code_dir) 
}





### Plot ages vs (original or one-year) MRI-T2 (data, model)
# I don't think I can do one-year MRI-T2 for data
setwd(output_dir)
dir.create("muscle_plots")
setwd(paste0(output_dir,"/","muscle_plots"))
for (l in 1:L) {
  m_name <- muscle_names[l]
  muscle_dat <- dat[, m_name]
  png(file=paste0(m_name, "_traj_dat.png"), 
      width=6.5, height=3.5, units="in", res=800, pointsize=6.9)
  par(mfrow=c(1,1))
  plot(ages_dat, muscle_dat,
       type = "n", # no points yet
       xlab = "Age",
       ylab = paste0(m_name, " (MRI-T2)"),
       main = paste0(m_name, " vs. Age (Data)"),
       xlim = range(ages_dat, na.rm=TRUE),
       ylim = range(muscle_dat, na.rm=TRUE)
  )
  for (i in 1:N) {
    subject_data <- subset(dat, ID == i)
    lines(subject_data$TIME,
          subject_data[, m_name],
          col = subject_colors[i],
          lwd = 1)
  }
  if (N_dropped > 0) {
    legend("topright", 
           legend = c(paste0("Slower Progressors (N=", N_subset, ")"), 
                      paste0("Faster Progressors (N=", N_dropped, ")")),
           col = c("navy", "darkred"),
           lty = 1, lwd = c(1, 1), bty = "n")
  }
  dev.off()
}
  


for (l in 1:L) {
  m_name <- muscle_names[l]
  muscle_dat <- dat[, m_name]
  png(file=paste0(m_name,"_pop_traj.png"),
      width=6.5, height=3.5, units="in", res=800, pointsize=6.9)
  plot(1, type="n",
       xlim = range(age_grid),
       ylim = c(min(muscle_dat, na.rm=TRUE),
                max(muscle_dat, na.rm=TRUE) * 1.1),
       xlab = "Age",
       ylab = paste0(m_name, " (MRI-T2)"),
       main = paste0(m_name, " vs. Age (Model)")
  )
  
  # (make) posterior mean population trajectory + pointwise 95% posterior credible bands
  draws_muscle_pop <- t(apply(joint_model_mat, 1, function(x) {
    get_muscle_trajectory(age_grid,
                          x[paste0("log_bar_S0_", m_name)],
                          x[paste0("log_bar_DP_", m_name)],
                          x[paste0("log_bar_DPT_", m_name)],
                          exp(x[paste0("log_bar_gamma_", m_name)]))
  }))
  
  predicted_muscle_pop <- colMeans(draws_muscle_pop, na.rm=TRUE)
  band_muscle_pop <- apply(draws_muscle_pop, 2, quantile,
                           probs = c(0.025, 0.975), na.rm=TRUE)
  
  # individual trajectories
  #for (i in 1:N) {
    #log_S0_T2_i <- mean(muscle_param_list[[m_name]]$log_S0[, i])
    #log_DP_T2_i <- mean(muscle_param_list[[m_name]]$log_DP[, i])
    #log_DPT_T2_i <- mean(muscle_param_list[[m_name]]$log_DPT[, i])
    #gamma_T2_i <- exp(mean(muscle_param_list[[m_name]]$log_gamma[, i]))
    #predicted_muscle_i <- 
      #get_muscle_trajectory(age_grid,
                            #log_S0_T2_i, log_DP_T2_i, log_DPT_T2_i, gamma_T2_i)
    #line_col <- adjustcolor(subject_colors[i], alpha.f = 0.15)
    #lines(age_grid,
          #predicted_muscle_i,
          #col = line_col,
          #lwd=1)
  #}
  
  # (plot) posterior mean population trajectory + pointwise 95% posterior credible bands
  polygon(c(age_grid, rev(age_grid)),
          c(band_muscle_pop[1, ], rev(band_muscle_pop[2, ])),
          col = adjustcolor("gray", alpha.f = 0.35),
          border = NA)
  lines(age_grid, predicted_muscle_pop, col="black", lwd=1.5)
  
  if (N_dropped > 0) {
    legend("topright",
           legend = c("Population Mean", "95% Credible Band",
                      paste0("Slower Progressors (N=", N_subset, ")"), 
                      paste0("Faster Progressors (N=", N_dropped, ")")),
           col = c("black", adjustcolor("gray", alpha.f = 0.35), "navy", "darkred"),
           lty = c(1, NA, 1, 1),
           lwd = c(1.5, NA, 1, 1),
           pch = c(NA, 15, NA, NA),
           pt.cex = c(NA, 2, NA, NA),
           bty = "n")
  } else {
    legend("topright", 
           legend = c("Population Mean", "95% Credible Band"),
           col = c("black", adjustcolor("gray", alpha.f = 0.35)),
           lty = c(1, NA),
           lwd = c(1.5, NA),
           pch = c(NA, 15),
           pt.cex = c(NA, 2),
           bty = "n")
  }
  
  dev.off()
}

  

for (l in 1:L) {
  m_name <- muscle_names[l]
  png(file=paste0(year_diff,"year_",m_name,"_traj_model.png"),
      width=6.5, height=3.5, units="in", res=800, pointsize=6.9)
  plot(1, type="n",
       xlim = range(age_grid),
       ylim = c(0, 20),
       xlab = "Initial Age",
       ylab = paste0(m_name, " (", year_diff, "-Year Change in MRI-T2)"),
       main = paste0(year_diff, "-Year Change in ", m_name, " vs. Initial Age (Model)")
  )
  for (i in 1:N) {
    # technically I should maybe get the posterior mean/median of
    # all trajectory values (slightly more computation time/complexity)
    # but using posterior mean/median of parameters then getting trajectory
    # should be sufficient visual summary/representation
    log_S0_T2_i <- mean(muscle_param_list[[m_name]]$log_S0[, i])
    log_DP_T2_i <- mean(muscle_param_list[[m_name]]$log_DP[, i])
    log_DPT_T2_i <- mean(muscle_param_list[[m_name]]$log_DPT[, i])
    gamma_T2_i <- exp(mean(muscle_param_list[[m_name]]$log_gamma[, i]))
    predicted_muscle_i <- 
      get_muscle_trajectory(age_grid,
                            log_S0_T2_i, log_DP_T2_i, log_DPT_T2_i, gamma_T2_i)
    year_later_predicted_muscle_i <-
      get_muscle_trajectory(age_grid + year_diff,
                            log_S0_T2_i, log_DP_T2_i, log_DPT_T2_i, gamma_T2_i)
    line_col <- adjustcolor(subject_colors[i], alpha.f = 0.15)
    lines(age_grid,
          year_later_predicted_muscle_i - predicted_muscle_i,
          col = line_col,
          lwd=1)
  }
  
  # posterior mean of fixed effects
  log_S0_T2_pop <- mean(joint_model_mat[, paste0("log_bar_S0_",m_name)])
  log_DP_T2_pop <- mean(joint_model_mat[, paste0("log_bar_DP_", m_name)])
  log_DPT_T2_pop <- mean(joint_model_mat[, paste0("log_bar_DPT_", m_name)])
  gamma_T2_pop <- exp(mean(joint_model_mat[, paste0("log_bar_gamma_", m_name)]))
  
  predicted_muscle_pop <- 
    get_muscle_trajectory(age_grid,
                          log_S0_T2_pop, log_DP_T2_pop, log_DPT_T2_pop, gamma_T2_pop)
  year_later_predicted_muscle_pop <-
    get_muscle_trajectory(age_grid + year_diff,
                          log_S0_T2_pop, log_DP_T2_pop, log_DPT_T2_pop, gamma_T2_pop)
  
  lines(age_grid, 
        year_later_predicted_muscle_pop - predicted_muscle_pop, 
        col = "black", 
        lwd=1.5)
  
  if (N_dropped > 0) {
    legend("topright",
           legend = c("Population Mean",
                      paste0("Slower Progressors (N=", N_subset, ")"), 
                      paste0("Faster Progressors (N=", N_dropped, ")")),
           col = c("black", "navy", "darkred"),
           lty = 1, lwd = c(1.5, 1, 1), bty = "n")
  } else {
    legend("topright", 
           legend = c("Population Mean", "Individual Fits"),
           col = c("black", "lightblue"),
           lty = 1, lwd = c(1.5, 1), bty = "n") 
  }
  dev.off() 
}
setwd(code_dir)





##### Combined plot of ages vs 6MWD (model); ages vs MRI-T2 (model)
### repeated stuff
age_grid <- seq(min_age_search, max_age_search, length.out=200)
# 6MWD (make) posterior mean population trajectory
draws_6MWD_pop <- t(apply(joint_model_mat, 1, function(x) {
  get_SMWT_trajectory(age_grid,
                      exp(x["log_bar_S0_6MWD"]),
                      x["log_bar_DPT_6MWD"],
                      exp(x["log_bar_gamma_6MWD"]))
}))
predicted_6MWD_pop <- colMeans(draws_6MWD_pop, na.rm=TRUE)
### slightly modified repeated stuff
# MRI-T2 (make) posterior mean population trajectory
predicted_muscle_pop_mat <- matrix(NA, nrow=length(age_grid), ncol=L)
colnames(predicted_muscle_pop_mat) <- muscle_names
for (l in 1:L) {
  m_name <- muscle_names[l]
  
  draws_muscle_pop <- t(apply(joint_model_mat, 1, function(x) {
    get_muscle_trajectory(age_grid,
                          x[paste0("log_bar_S0_", m_name)],
                          x[paste0("log_bar_DP_", m_name)],
                          x[paste0("log_bar_DPT_", m_name)],
                          exp(x[paste0("log_bar_gamma_", m_name)]))
  }))
  
  predicted_muscle_pop_mat[, l] <- colMeans(draws_muscle_pop, na.rm = TRUE)
}
### new stuff
setwd(output_dir)
cairo_pdf(file = "combined_pop_traj.pdf",
          width=6.5, height=3.5, pointsize=6.9)
par(mfrow = c(1, 2),
    mar = c(4.2, 4.2, 2.5, 1))
# MRI-T2 (plot) posterior mean population trajectory
matplot(age_grid,
        predicted_muscle_pop_mat,
        type = "l",
        lty = 1,
        lwd = 2,
        col = muscle_colors,
        xlab = "Age",
        #ylab = "MRI_T2",
        #main = "Posterior mean population trajectories of MRI_T2 of leg muscles"
        ylab = expression(MRI[T2]),
        main = expression("Leg muscle " *  MRI[T2] * " posterior mean population trajectories")
        )
legend("bottomright",
       legend = muscle_names,
       col = muscle_colors,
       lty = 1,
       lwd = 2,
       bty = "n",
       cex = 0.8)
# 6MWD (plot) posterior mean population trajectory
matplot(age_grid,
     predicted_6MWD_pop,
     type = "l",
     lty = 1,
     lwd = 2,
     col = "black",
     xlab = "Age",
     #ylab = "6MWD",
     #main = "Posterior mean population trajectory of 6MWD"
     ylab = expression("6MWD"),
     main = expression("6MWD posterior mean population trajectory")
     )
dev.off()
setwd(code_dir)










### (initialize) Step 2.1 of grid algorithm
array_dim1 <- c(N, grid_size, noi_target)
array_names1 <- list(paste0("ID_",1:N),
                     paste0("SMWT_",SMWT_grid),
                     paste0("iter_",iters))
# memory issues (also currently not used later after algorithm finishes)
bounded_ages <- all_ages <- error_search <-
  array(data=NA, dim=array_dim1, dimnames=array_names1)


### (initialize) Step 2.2 of grid algorithm
matrix_names1 <- list(paste0("ID_",1:N),
                      paste0("SMWT_",SMWT_grid))
mean_bounded_ages <- matrix(data=NA, nrow=N, ncol=grid_size,
                            dimnames=matrix_names1)

### (initialize) Steps 2.5-2.6 of grid algorithm
optim_method <- "L-BFGS-B" # Nelder-Mead, BFGS, CG, L-BFGS-B
grad_bool <- TRUE # for optimization (doesn't matter unless using CG)
upper_bool <- TRUE # for optimization (doesn't matter but why not)
matrix_names2 <- list(paste0("SMWT_",SMWT_grid),
                      muscle_names)
w_mat <- matrix(data=NA, nrow=grid_size, ncol=L,
                dimnames=matrix_names2)

array_dim3 <- c(grid_size, noi_target, L)
array_names3 <- list(paste0("SMWT_",SMWT_grid),
                     paste0("iter_",iters),
                     muscle_names)
w_it_array <- array(data=NA, dim=array_dim3, dimnames=array_names3)



### 6MWD Grid Algo
matrix_init <- matrix(data=NA, nrow=N, ncol=grid_size)

cl <- makeCluster(ncores)
registerDoParallel(cl)

start_time <- Sys.time()

results_list <- foreach(it_indx = 1:noi_target, 
                        .errorhandling = "pass"
) %dopar%
{
  it <- iters[it_indx]
  
  # diagnostics (save in bigger matrix or array)
  bounded_ages_temp <- matrix_init
  error_search_temp <- matrix_init
  
  # save in bigger matrix or array
  w_it_mat_temp <- matrix(data=NA, nrow=grid_size, ncol=L)
  
  # not saved/returned but calculations
  one_year_muscles_it_temp <- array(data=NA, dim=c(N, grid_size, L))

  for (grid_indx in 1:grid_size) {
    s <- SMWT_grid[grid_indx]
    
    for (i in subjects_of_interest) {
      
      ### Step 1.1 of grid algorithm
      S0_F_temp <- exp(log_S0_6MWD_mat[it, i])
      log_DPT_F_temp <- log_DPT_6MWD_mat[it, i]
      gamma_F_temp <- exp(log_gamma_6MWD_mat[it, i])
      
      # see if any age becomes less than 'min_age_search' or above 'max_age_search'
      tryCatch({
        age_F_temp <- uniroot(
          f = stable_SMWT_function0,
          target_SMWT = s,
          S0_F = S0_F_temp,
          log_DPT_F = log_DPT_F_temp,
          gamma_F = gamma_F_temp,
          interval = c(min_age_search, max_age_search), 
          extendInt = "yes" )$root 
        # note: there is more info other than root
        
        if ( (age_F_temp < min_age_search) || (age_F_temp > max_age_search) ) {
          error_search_temp[i, grid_indx] <- age_F_temp
        }
        if ( (age_F_temp >= min_age_search) && (age_F_temp <= max_age_search) ) {
          bounded_ages_temp[i, grid_indx] <- age_F_temp
        }
        
      }, error=function(e){
        error_search_temp[i, grid_indx] <- -1
      }
      )
      
      ### Step 1.2 of grid algorithm
      if (is.finite(bounded_ages_temp[i, grid_indx])) {
        for (l in 1:L) {
          m_name <- muscle_names[l]
          
          # these could be either the FE or RE value
          # (depending on if RE was omitted)
          log_S0_T2_temp <- muscle_param_list[[m_name]]$log_S0[it, i]
          log_DP_T2_temp <- muscle_param_list[[m_name]]$log_DP[it, i]
          log_DPT_T2_temp <- muscle_param_list[[m_name]]$log_DPT[it, i]
          gamma_T2_temp <- exp(muscle_param_list[[m_name]]$log_gamma[it, i])
          
          # iterative
          current_T2_it_temp <- stable_muscle_function(
            bounded_ages_temp[i, grid_indx],
            log_S0_T2_temp, log_DP_T2_temp, log_DPT_T2_temp, gamma_T2_temp)
          year_later_T2_it_temp <- stable_muscle_function(
            bounded_ages_temp[i, grid_indx] + year_diff,
            log_S0_T2_temp, log_DP_T2_temp, log_DPT_T2_temp, gamma_T2_temp)
          one_year_muscles_it_temp[i, grid_indx, l] <-
            year_later_T2_it_temp - current_T2_it_temp
        }
      }
      
    }
    
    ### Steps 1.3-1.4 of grid algorithm
    w_it_temp <- rep(NA, L)
    
    valid_subjects_it <- 
      subjects_of_interest[is.finite(
        bounded_ages_temp[subjects_of_interest, grid_indx])]
    n_valid_subjects_it <- length(valid_subjects_it)
    # don't compute weights if too few data
    if ( (n_valid_subjects_it >= floor(0.1*N_subset)) 
         && (n_valid_subjects_it >=2) )  {
      subset_it_results <- 
        one_year_muscles_it_temp[valid_subjects_it, grid_indx, ]
      # b/c 'one_year_muscles_it_temp' could contain NaN, Inf, or -Inf
      # if 'stable_muscle_function' has numerical instability
      # (at least two subjects with complete data?)
      if (sum(complete.cases(subset_it_results)) >= 2) {
        mu_tilde_it_temp <- colMeans(subset_it_results, na.rm=TRUE)
        Sigma_tilde_it_temp <- cov(subset_it_results, use="pairwise.complete.obs")
        # (covariance matrix doesn't contain NA's?)
        if (!any(is.na(Sigma_tilde_it_temp))) {
          # (covariance matrix is numerically stable and non-singular?)
          if (rcond(Sigma_tilde_it_temp) > .Machine$double.eps) {
            w_stuffs_it_temp <- 
              muscle_optim_func(mu_tilde_it_temp, Sigma_tilde_it_temp,
                                optim_method, grad_bool, upper_bool)
            w_it_temp <- w_stuffs_it_temp$w
          }
        }
      }
    }
    w_it_mat_temp[grid_indx, ] <- w_it_temp
  }
  
  return(list(
    # diagnostics (save in bigger matrix or array)
    bounded_ages_slice = bounded_ages_temp,
    error_search_slice = error_search_temp,
    
    # save in bigger matrix or array
    weights_it_slice = w_it_mat_temp
    # MAKE SURE THIS DOESN'T END WITH COMMA!!!!!
  ))
}

end_time <- Sys.time()
stopCluster(cl)

algo_time <- round(as.numeric(difftime(end_time, start_time, units='hours')),2)
trial_name <- paste0("noi=", noi, ", noi_target=", noi_target,
                     ", grid_tick=", grid_tick,
                     ", N_subset=", N_subset, ", N_dropped=", N_dropped)
algo_time_and_more <-
  paste0(output_file_name, " (", trial_name, "): ", algo_time, " hours")
setwd(output_dir)
write(algo_time_and_more, file="algo_time.txt", append=TRUE)
setwd(code_dir)
beep(3)





# loop through the results list and fill them slice by slice
for (it_indx in 1:noi_target) {
  results_temp <- results_list[[it_indx]]
  
  # Check if the worker for this grid_index produced an error
  if (!inherits(results_temp, "error")) {
    
    # diagnostics (saved in bigger matrix or array)
    bounded_ages[, , it_indx] <- results_temp$bounded_ages_slice
    error_search[, , it_indx] <- results_temp$error_search_slice
    
    # saved in bigger matrix or array
    w_it_array[, it_indx , ] <- results_temp$weights_it_slice
  } 
}










cl2 <- makeCluster(ncores)
registerDoParallel(cl2)

start_time2 <- Sys.time()

# (novel) point estimator of weights
results_list2 <- foreach(grid_indx = 1:grid_size, 
                        .errorhandling = "pass"
) %dopar%
{
  s <- SMWT_grid[grid_indx]
  
  # save in bigger matrix or array
  mean_bounded_ages_s <- rep(NA, N)
  
  # not saved/returned but calculations
  one_year_muscles_s <- array(data=NA, dim=c(N, noi_target, L))
  mean_one_year_muscles_s <- matrix(data=NA, nrow=N, ncol=L)
  
  for (i in subjects_of_interest) {
    ### Step 2.2
    n_successful_iters <- sum(!is.na(bounded_ages[i, grid_indx, ]))
    if (n_successful_iters >= floor(0.1*noi_target)) {
      # note: mean(..., na.rm=TRUE) produces 'NaN' when all values = NA
      mean_bounded_ages_s[i] <- 
        mean(bounded_ages[i, grid_indx, ], na.rm=TRUE) ### MODIFY AS NEEDED
    }
    
    ### Step 2.3 of grid algorithm (Dr. Daniels said used 'bounded_ages')
    if (is.finite(mean_bounded_ages_s[i])) {
      for (l in 1:L) {
        m_name <- muscle_names[l]
        
        for (it_indx in 1:noi_target) {
          it <- iters[it_indx]
          
          # these could be either the FE or RE value
          # (depending on if RE was omitted)
          log_S0_T2_temp <- muscle_param_list[[m_name]]$log_S0[it, i]
          log_DP_T2_temp <- muscle_param_list[[m_name]]$log_DP[it, i]
          log_DPT_T2_temp <- muscle_param_list[[m_name]]$log_DPT[it, i]
          gamma_T2_temp <- exp(muscle_param_list[[m_name]]$log_gamma[it, i])
          
          # pt. est.
          current_T2_temp <- stable_muscle_function(
            mean_bounded_ages_s[i],
            log_S0_T2_temp, log_DP_T2_temp, log_DPT_T2_temp, gamma_T2_temp)
          year_later_T2_temp <- stable_muscle_function(
            mean_bounded_ages_s[i] + year_diff,
            log_S0_T2_temp, log_DP_T2_temp, log_DPT_T2_temp, gamma_T2_temp)
          
          one_year_muscles_s[i, it_indx, l] <- 
            year_later_T2_temp - current_T2_temp
        }
        
        ### Step 2.4
        # note that averaging over the trajectories is more robust
        # than just using posterior means/medians of muscle parameters
        # and getting a single trajectory value with the 
        # 'mean_bounded_ages_s'
        mean_one_year_muscles_s[i,l] <- 
          mean(one_year_muscles_s[i,,l], na.rm=TRUE)
      }
    }
    
  }
  
  ### Steps 2.5-2.6 of grid algorithm
  # if only like two individuals have age for a certain grid value,
  # we prob don't want to use it (threshold: 10% data at least)
  # ('is.na' detects 'NaN's as well 
  # when bounded_ages[i, grid_indx, ] contains only NA values
  # across all iterations)
  
  # 'subjects_of_interest' could be all subjects, 'unique_ID' == 1:N
  # 'N_subset' could be 'N'
  w_temp <- rep(NA, L)
  
  valid_subset_indices <- 
    subjects_of_interest[is.finite(mean_bounded_ages_s[subjects_of_interest])]
  n_valid_in_subset <- length(valid_subset_indices)
  # don't compute weights if too few data
  if ( (n_valid_in_subset >= floor(0.1 * N_subset)) 
       && (n_valid_in_subset >=2) ) {
    subset_results <- 
      mean_one_year_muscles_s[valid_subset_indices, ]
    # b/c 'mean_one_year_muscles_s' could contain NaN, Inf, or -Inf
    # if 'stable_muscle_function' has numerical instability
    # (at least two subjects with complete data?)
    if (sum(complete.cases(subset_results)) >= 2) {
      mu_tilde_temp <- colMeans(subset_results, na.rm=TRUE)
      Sigma_tilde_temp <- cov(subset_results, use="pairwise.complete.obs")
      # (covariance matrix doesn't contain NA's?)
      if (!any(is.na(Sigma_tilde_temp))) {
        # (covariance matrix is numerically stable and non-singular?)
        if (rcond(Sigma_tilde_temp) > .Machine$double.eps) {
          w_stuffs_temp <- muscle_optim_func(mu_tilde_temp, Sigma_tilde_temp,
                                             optim_method, grad_bool, upper_bool)
          #w_orig_temp <- w_stuffs_temp$w_orig
          w_temp <- w_stuffs_temp$w
        }
      }
    }
  }
  
  return(list(
    weights_pt_est_slice = w_temp,
    mean_bounded_ages_slice = mean_bounded_ages_s
    #mean_one_year_muscles_slice = mean_one_year_muscles_s
    
    # MAKE SURE THIS DOESN'T END WITH COMMA!!!!!
  ))
  
}

end_time2 <- Sys.time()
stopCluster(cl2)

algo_time2 <- round(as.numeric(difftime(end_time2, start_time2, units='hours')),2)
algo_time3 <- algo_time + algo_time2
algo_time2_and_time3 <-
  paste0(algo_time2, " hours; total = ", algo_time+algo_time2, " hours.")
setwd(output_dir)
write(algo_time2_and_time3, file="algo_time.txt", append=TRUE)
write("", file="algo_time.txt", append=TRUE)
setwd(code_dir)
beep(5)





# loop through the results list and fill them slice by slice
for (grid_indx in 1:grid_size) {
  results_temp2 <- results_list2[[grid_indx]]
  
  # Check if the worker for this grid_index produced an error
  if (!inherits(results_temp2, "error")) {
    # saved in bigger matrix or array
    w_mat[grid_indx, ] <- results_temp2$weights_pt_est_slice
    mean_bounded_ages[, grid_indx] <- results_temp2$mean_bounded_ages_slice
  }
}



setwd(output_dir)
write.csv(mean_bounded_ages, "mean_bounded_ages.csv")
setwd(code_dir)










### (simple diagnostics) Step 1.1
if (grid_tick >=10) {
  dim(which(((error_search < min_age_search) & (error_search > -1)), arr.ind=TRUE))
  dim(which(error_search > max_age_search, arr.ind=TRUE))
  dim(which(error_search == -1, arr.ind=TRUE))
  length(error_search)
}





### (summary of errors) Step 1.1
summary_vars <- c("*N_no_age_solution",
                  "*min_SMWT_no_age_solution",
                  "*BOOL_no_age_solution_after",
                  "*N_age_solution_after",
                  "N_smaller_than_min_age_search",
                  "N_larger_than_max_age_search")
error_summary <- array(data=NA, dim=c(N, length(summary_vars), noi_target))
dimnames(error_summary) <- list(paste0("ID_",1:N),
                                summary_vars,
                                paste0("iter_",iters))
for (it_indx in 1:noi_target) {
  
  for (i in 1:N) {
    error_search_it_i <- error_search[i,,it_indx]
    
    no_age_indxs <- which(error_search_it_i==-1)
    
    if (length(no_age_indxs) > 0) {
      n_no_age <- length(no_age_indxs)
      
      first_no_age_indx <- no_age_indxs[1]
      SMWT_first_no_age <- SMWT_grid[first_no_age_indx]
      
      slice_after_first_no_age <- error_search_it_i[first_no_age_indx:grid_size]
      
      n_age_soln_after <- length(slice_after_first_no_age) - 
        sum(slice_after_first_no_age == -1, na.rm = TRUE)
      
      if (n_age_soln_after==0) {
        NO_age_soln_after <- TRUE
      } else {
        NO_age_soln_after <- FALSE
      }
      
      error_summary[i,1:4,it_indx] <-
        c(n_no_age, SMWT_first_no_age, NO_age_soln_after, n_age_soln_after)
    }
    
    smaller_than_min_age_search_indxs <- 
      which((error_search_it_i < min_age_search) & (error_search_it_i > -1))
    if (length(smaller_than_min_age_search_indxs) > 0) {
      error_summary[i,5,it_indx] <- length(smaller_than_min_age_search_indxs)
    }
    
    larger_than_max_age_search_indxs <- 
      which(error_search_it_i > max_age_search)
    if (length(larger_than_max_age_search_indxs) > 0) {
      error_summary[i,6,it_indx] <- length(larger_than_max_age_search_indxs)
    }
    
    
  }
}

it_indx_interest <- noi_target
setwd(output_dir)
write.csv(error_summary[,,it_indx_interest], 
          paste0("error_summ_iter",iters[it_indx_interest],".csv"))
setwd(code_dir)





### (PLOTS i.e. Sanity check) Step 2.2
setwd(output_dir)
png(file=paste0("mean_bounded_ages.png"),
    width=6.5, height=3.5, units="in", res=800, pointsize=6.9)
par(mfrow=c(1,1))
plot(x=1, y=1, type="n", 
     xlim=c(min(mean_bounded_ages, na.rm=TRUE),
            max(mean_bounded_ages, na.rm=TRUE)),
     ylim=c(min(SMWT_grid), max(SMWT_grid)), 
     xlab=paste0("Age"),
     ylab=paste0("6MWD grid value (sequence of ",grid_tick," from ",min_SMWT_grid," to ",max_SMWT_grid,")"),
     main=paste0("Solved age (within the (",min_age_search,",",max_age_search,") bound)"),
     #sub=paste0()
)
for (i in 1:N) {
  nonNA_indxs_temp <- 
    as.numeric(which(!is.na(mean_bounded_ages[i,]), arr.ind=TRUE))
  mean_bounded_ages_i <- 
    as.numeric(mean_bounded_ages[i, nonNA_indxs_temp])
  lines(mean_bounded_ages_i, SMWT_grid[nonNA_indxs_temp], col=i)
}
dev.off()
setwd(code_dir)










# set up for plots
w_mat_ours <- w_mat
w_mat_mean <- 
  apply(w_it_array, c(1, 3), mean, na.rm = TRUE)

w_mat_choice <- w_mat_mean

if (all(w_mat_choice==w_mat_ours, na.rm=TRUE)) {
  post_est_label <- "Posterior estimate of"
} else if (all(w_mat_choice==w_mat_mean, na.rm=TRUE)) {
  post_est_label <- "Posterior mean"
}

notNA_rows <- which(complete.cases(w_mat_choice))
notNA_rows <- notNA_rows[-length(notNA_rows)] # sometimes plot looks funky for largest 6MWD index
n_grid_points_plot <- length(notNA_rows)
x_values <- 1:n_grid_points_plot

w_mat_clean <- w_mat_choice[notNA_rows, ]
SMWT_grid_clean <- SMWT_grid[notNA_rows]
w_it_array_clean <- w_it_array[notNA_rows, , ]





### Calculate Posterior Probabilities of P(Weight_A > Weight_B)
prob_matrix_list <- vector("list", n_grid_points_plot)
names(prob_matrix_list) <- paste0("SMWT_", SMWT_grid_clean)

for (grid_indx in 1:n_grid_points_plot) {
  prob_matrix <- matrix(NA, nrow = L, ncol = L,
                        dimnames = list(muscle_names, muscle_names))
  for (l1 in 1:L) {
    for (l2 in 1:L) {
      if (l1==l2) {
        next # this is to keep diagonal as NAs (cleaner)
      }
      draws_l1 <- w_it_array_clean[grid_indx, , l1]
      draws_l2 <- w_it_array_clean[grid_indx, , l2]
      
      # posterior P(row_muscle > col_muscle)
      prob_l1_gt_l2 <- mean(draws_l1 > draws_l2)
      
      prob_matrix[l1, l2] <- prob_l1_gt_l2
    }
  }
  
  prob_matrix_list[[grid_indx]] <- prob_matrix
}
setwd(output_dir)
dir.create("prob_row_gt_col")
setwd(paste0(output_dir,"/","prob_row_gt_col"))
for (grid_indx in 1:n_grid_points_plot) {
  write.csv(round(prob_matrix_list[[grid_indx]], 3),
            paste0("SMWT",SMWT_grid[grid_indx],".csv"))
}
setwd(code_dir)


# Average of Posterior Probabilities P(Weight_A > Weight_B)
avg_prob_matrix <- 
  calculate_average_matrix(prob_matrix_list, 1:n_grid_points_plot)
indices_le396 <- which(SMWT_grid_clean <= grid_val_target)
indices_gt396 <- which(SMWT_grid_clean > grid_val_target)
avg_prob_matrix_le396 <-
  calculate_average_matrix(prob_matrix_list, indices_le396)
avg_prob_matrix_gt396 <- 
  calculate_average_matrix(prob_matrix_list, indices_gt396)

setwd(output_dir)
write.csv(round(avg_prob_matrix,3),
          "prob_row_gt_col_avg.csv")
write.csv(round(avg_prob_matrix_le396,3), 
          paste0("prob_row_gt_col_avg_le",grid_val_target,".csv"))
write.csv(round(avg_prob_matrix_gt396,3),
          paste0("prob_row_gt_col_avg_gt",grid_val_target,".csv"))
setwd(code_dir)


# score = rowMeans of avg_prob_matrix
overall_scores <- rowMeans(avg_prob_matrix, na.rm = TRUE)
le396_scores <- rowMeans(avg_prob_matrix_le396, na.rm = TRUE)
gt396_scores <- rowMeans(avg_prob_matrix_gt396, na.rm = TRUE)
score_matrix <- rbind(overall_scores, le396_scores, gt396_scores)
rownames(score_matrix) <- c("Overall", paste0("<=SMWT", grid_val_target),
                            paste0(">SMWT", grid_val_target))
setwd(output_dir)
write.csv(round(score_matrix,3), "muscle_scores.csv")
setwd(code_dir)










# Posterior densities of weights for each iteration
setwd(output_dir)
dir.create(paste0(output_dir,"/","w_dist_it"))
setwd(paste0(output_dir,"/","w_dist_it"))
w_it_array_clean_df_long <- as.data.frame.table(w_it_array_clean)
colnames(w_it_array_clean_df_long) <- 
  c("grid_index", "iteration", "muscle", "weight")
# convert 'grid_index' column to numeric
w_it_array_clean_df_long$SMWT <- 
  SMWT_grid_clean[as.numeric(w_it_array_clean_df_long$grid_index)]

for (grid_indx in 1:n_grid_points_plot) {
  s_clean <- SMWT_grid_clean[grid_indx]
  df_subset <- w_it_array_clean_df_long %>%
    filter(SMWT == s_clean) %>%
    filter(!is.na(weight))
  
  if (nrow(df_subset) == 0) {
    next
  }
  
  p <- ggplot(df_subset, aes(x = weight)) +
    
    geom_density(
      aes(fill = muscle),
      bounds = c(0, 1),
      alpha = 0.8) +
    
    facet_wrap(~ muscle, scales = "free_y") +
    
    labs(
      title = paste0("Posterior Densities of Muscle Weights at 6MWD = ", s_clean),
      x = "Weight",
      y = "Density"
    ) +
    
    scale_fill_manual(values = muscle_colors) +
    
    theme_bw(base_size = 7) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold") # bold facet titles
    )
    
  # The filename includes the 6MWD value, so each plot is unique.
  ggsave(
    filename = paste0("w_dens_6MWD", s_clean, ".png"),
    plot = p,
    width=6.5, height=3.5, units="in",dpi=800
  )
}
setwd(code_dir)










# CSV of weight posterior means and credible intervals (for every 10 multiple of 6MWD)
the_multiple <- 10
SMWT_grid_report <- c(1, SMWT_grid_clean[SMWT_grid_clean %% 10 == 0])

w_mean_and_ci_df <- data.frame()
for (grid_indx in 1:length(SMWT_grid_report)) {
  s_report <- SMWT_grid_report[grid_indx]
  
  for (l in 1:L) {
    muscle_name <- muscle_names[l]
    
    w_temp <- w_it_array_clean[grid_indx, , l]
    
    w_post_mean_temp <- mean(w_temp, na.rm = TRUE)
    w_cri_temp <- quantile(w_temp, probs = c(0.025, 0.975), na.rm = TRUE)
    
    # Append to the results data frame
    w_mean_and_ci_df <- rbind(w_mean_and_ci_df, data.frame(
      SMWT = s_report,
      Muscle = muscle_name,
      Weight = round(w_post_mean_temp, 3),
      CI_Lower = round(w_cri_temp[1], 3),
      CI_Upper = round(w_cri_temp[2], 3)
    ))
  }
}

setwd(output_dir)
write.csv(w_mean_and_ci_df, paste0("w_mean_and_ci_by_",the_multiple,".csv"),
          row.names = FALSE)
setwd(code_dir)








### one year T2 weight vs 6MWD plots
grid_sample_size <- colSums(!is.na(mean_bounded_ages))
grid_sample_size_clean <- as.numeric(grid_sample_size[notNA_rows])
# note that weights weren't calculated for a grid index in the algorithm
# if the grid's sample size is less than 11 (10% of overall sample size)
print(which(grid_sample_size_clean < floor(0.1*N)))

muscle_symbols <- c(1, 2, 15, 4, 5)


if (all(w_mat_clean == w_mat_mean[notNA_rows, ])) {
  w_plot_name <- paste0(year_diff,"yearW_mean_vs_6MWD")
} else {
  w_plot_name <- paste0(year_diff,"yearW_novel_vs_6MWD")
}
if (N_dropped > 0) {
  w_plot_name <- paste0("NSubset_",w_plot_name)
}

setwd(output_dir)
#png(file=paste0(w_plot_name,".png"), 
#    width=6.5, height=3.5, units="in", res=800, pointsize=6.9)
cairo_pdf(file = paste0(w_plot_name, ".pdf"),
          width=6.5, height=3.5, pointsize=6.9)

par(mfrow = c(1, 1),
    mar = c(5.7, 4.9, 3.6, 1.0),
    cex.main = 0.95)

plot(1, type = "n",  # no points
     xlim = c(n_grid_points_plot, 1),
     #ylim = c(min(w_mat_clean), max(w_mat_clean) * 1.05),
     ylim = c(0,1),
     xlab = "",  # custom x-axis later
     #ylab = "Weight",
     ylab = "",
     xaxt = 'n', # custom x-axis later
     yaxt = 'n', # custom y-axis later
     main =
     )

# Credible Intervals (plot first as background?)
# transparent lines don't help
for (l in 1:L) {
  lower_bound <- apply(w_it_array_clean[,,l], 1, quantile, probs=0.025, na.rm=TRUE)
  upper_bound <- apply(w_it_array_clean[,,l], 1, quantile, probs=0.975, na.rm=TRUE)
  
  rgb_vals <- col2rgb(muscle_colors[l]) / 255
  
  arrows(x0 = x_values,
         y0 = lower_bound,
         x1 = x_values,
         y1 = upper_bound,
         code = 3, # draw bars at end of line
         angle = 90, # bars perpendicular to line
         length = 0.03,
         # more transparent vertical lines:
         col = rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], alpha = 0.3),
         lwd = 0.5) # width
  
  # more opaque endpoints
  segments(x0 = x_values - 0.02, x1 = x_values + 0.02,
           y0 = lower_bound, y1 = lower_bound,
           col = rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], alpha = 0.8),
           lwd = 1)
  segments(x0 = x_values - 0.02, x1 = x_values + 0.02,
           y0 = upper_bound, y1 = upper_bound,
           col = rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], alpha = 0.8),
           lwd = 1)
  
  # Credible Bands
  #polygon(c(x_values, rev(x_values)), c(lower_bound, rev(upper_bound)),
  #        col = rgb(t(col2rgb(l))/255, alpha=0.2), border = NA)
  
}

#matlines(x = x_values, y = w_mat_clean, 
#         lty = 1, col = 1:L, lwd = 1.5)
#matpoints(x = x_values, y = w_mat_clean,
#          pch = 1:L, col = 1:L, cex = 0.8)
for (l in 1:L) {
  points(x = x_values, 
         y = w_mat_clean[, l], ### MODIFY AS NEEDED
         pch = muscle_symbols[l],  # custom symbols
         col = muscle_colors[l],   # custom colors
         bg = muscle_colors[l],    # fill in color for some symbols
         cex = 0.8)
}

if (grid_tick < 5) {
  #x_axis_indices <- seq(1, n_grid_points_plot, by = 5)
  x_axis_indices <- c(1, which(SMWT_grid_clean %% 10 == 0))
} else {
  x_axis_indices <- 1:n_grid_points_plot
}
x_axis_labels <- paste0(SMWT_grid_clean[x_axis_indices], 
                        " (", grid_sample_size_clean[x_axis_indices], ")")
axis(side = 1,
     at = x_values[x_axis_indices],
     labels = x_axis_labels,
     cex.axis = 0.75,
     las = 2)

axis(side = 2,
     cex.axis = 0.85,
     las = 1)

# axis labels and title
mtext("6MWD (sample size)",
      side = 1,
      line = 3.9,
      cex = 1.2,
      las = 1)

mtext("Weight",
      side = 2,
      line = 2.8,
      cex = 1.2,
      las = 0)

title(main = bquote(.(post_est_label) ~ " weights for 1-year " ~ MRI[T2] ~
                      " changes in lower extremity muscles vs. 6MWD"),
      line = 2,
      cex.main = 1.2)

legend("top", 
       legend = muscle_names, 
       pch = muscle_symbols,
       col = muscle_colors,
       pt.bg = muscle_colors, # fill in colors for some symbols
       #lty = 1,
       #lwd=1.5,
       horiz = TRUE,
       cex = 0.9,
       pt.cex = 1.1,
       xpd = TRUE,
       inset = c(0, -0.07),
       bty = "n"
)

dev.off()
setwd(code_dir)










# 2 x 3 figure for each of the 5 plots:
setwd(output_dir)
cairo_pdf(file = paste0(w_plot_name, "_2x3.pdf"),
          #width=6.5, height=3.5, pointsize=6.9)
          width=12, height=8, pointsize=8)
par(mfrow = c(2, 3),
    mar = c(5.4, 4.2, 3.8, 1.2),
    oma = c(4.2, 4.2, 4.2, 0.8),
    cex.axis = 0.95,
    cex.main = 1.65)
for (l in 1:L) {
  
  plot(1, type = "n",  # no points
       xlim = c(n_grid_points_plot, 1),
       #ylim = c(min(w_mat_clean), max(w_mat_clean) * 1.05),
       ylim = c(0,1),
       xlab = "",  # custom x-axis later
       #ylab = "Weight",
       ylab = "",
       xaxt = 'n', # custom x-axis later
       yaxt = "n", # custom y-axis later
       main = paste0(muscle_names[l], " vs. 6MWD")
       #main = paste0(post_est_label,
                     #" weights for 1-year MRI-T2 changes in ",muscle_names[l]," vs. 6MWD")
  )
  
  lower_bound <- apply(w_it_array_clean[,,l], 1, quantile, probs=0.025, na.rm=TRUE)
  upper_bound <- apply(w_it_array_clean[,,l], 1, quantile, probs=0.975, na.rm=TRUE)
  
  rgb_vals <- col2rgb(muscle_colors[l]) / 255
  
  arrows(x0 = x_values,
         y0 = lower_bound,
         x1 = x_values,
         y1 = upper_bound,
         code = 3, # draw bars at end of line
         angle = 90, # bars perpendicular to line
         length = 0.03,
         # more transparent vertical lines:
         col = rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], alpha = 0.3),
         lwd = 0.5) # width
  
  # more opaque endpoints
  segments(x0 = x_values - 0.02, x1 = x_values + 0.02,
           y0 = lower_bound, y1 = lower_bound,
           col = rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], alpha = 0.8),
           lwd = 1)
  segments(x0 = x_values - 0.02, x1 = x_values + 0.02,
           y0 = upper_bound, y1 = upper_bound,
           col = rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], alpha = 0.8),
           lwd = 1)
  
  # Credible Bands
  #polygon(c(x_values, rev(x_values)), c(lower_bound, rev(upper_bound)),
  #        col = rgb(t(col2rgb(l))/255, alpha=0.2), border = NA)
  
  points(
    x = x_values,
    y = w_mat_clean[, l],
    pch = muscle_symbols[l],
    col = muscle_colors[l],
    bg = muscle_colors[l],
    cex = 1.15
  )
  
  if (grid_tick < 5) {
    x_axis_indices <- c(1, which(SMWT_grid_clean %% 10 == 0))
  } else {
    x_axis_indices <- 1:n_grid_points_plot

  }
  
  x_axis_labels <- paste0(
    SMWT_grid_clean[x_axis_indices],
    " (", grid_sample_size_clean[x_axis_indices], ")"
  )
  
  axis(side = 1,
       at = x_values[x_axis_indices],
       labels = x_axis_labels,
       cex.axis = 1.25,
       las = 2
       )
  
  axis(side = 2,
       cex.axis = 1.25,
       las = 1)
}
plot.new()
legend("center",
       legend = muscle_names,
       pch = muscle_symbols,
       col = muscle_colors,
       pt.bg = muscle_colors,
       bty = "n",
       cex = 1.65,
       pt.cex = 2.0)

mtext("6MWD (sample size)",
      outer = TRUE,
      side = 1,
      line = 2.4,
      cex = 1.85,
      las = 1)

mtext("Weight",
      outer = TRUE,
      side = 2,
      line = 0.8,
      cex = 1.85)

mtext(bquote(.(post_est_label) ~ " weights for 1-year " ~ MRI[T2] ~
               " changes in lower extremity muscles vs. 6MWD"),
      outer = TRUE,
      side = 3,
      line = 0.3,
      cex = 1.9)

dev.off()
setwd(code_dir)

beep(1)