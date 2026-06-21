time_txt <- function(mult_chains, ff_dat_name, burn_in_prop, nchains, ncores,
                     start_time, end_time, id="") {
  if (id!="") {
    id <- mult_chains[[1]]$id
  }
  orig_iters <- mult_chains[[1]]$iters
  burn_in_number <- orig_iters*burn_in_prop
  r_updates <- mult_chains[[1]]$r_updates
  
  r_jump_dist_name <- mult_chains[[1]]$r_jump_dist_name
  L <- mult_chains[[1]]$L
  N <- mult_chains[[1]]$N
  J_max <- mult_chains[[1]]$J_max
  
  nu <- mult_chains[[1]]$nu
  kappa <- mult_chains[[1]]$kappa
  
  runtime <- paste0(difftime(end_time, start_time, units='hours'))
  
  DA_MH_time <- paste0("id=",id,", ",ff_dat_name," > L=",L, " > r_JUMP_DIST=",r_jump_dist_name,
                       " (N=",N,", J_max=",J_max,", Q=",Q,")",
                       "; chains=",nchains,
                       ", cores=", ncores,
                       "; orig_iters=",orig_iters,
                       ", burn-in=",burn_in_number,
                       ", r_updates=[",paste0(r_updates,collapse="|"),"]",
                       "; nu=[",paste0(nu,collapse="|"),"]",
                       "; kappa=[",paste0(kappa,collapse="|"),"]",
                       "; ",runtime, " hours")
  
  write(DA_MH_time, file=paste0("time.txt"), append=TRUE)
  write("", file=paste0("time.txt"), append=TRUE)
}



stan_dat_single_prep <- function(dat_single) {
  min_y_obs <- min(dat_single, na.rm=TRUE)
  max_y_obs <- max(dat_single, na.rm=TRUE)
  
  y_obs <- dat_single
  y_obs_bool <- 1*(!is.na(y_obs))
  # replace NAs with sentinel value '-999' b/c STAN doesn't accept NAs:
  y_obs[is.na(y_obs)] <- -999
  
  return(list(
    min_y_obs=min_y_obs, max_y_obs=max_y_obs,
    y_obs=y_obs, y_obs_bool=y_obs_bool))
}



stan_dat_pair_prep <- function(dat_single1, dat_single2) {
  min_y_obs1 <- min(dat_single1, na.rm=TRUE)
  max_y_obs1 <- max(dat_single1, na.rm=TRUE)
  
  y_obs1 <- dat_single1
  y_obs_bool1 <- 1*(!is.na(y_obs1))
  # replace NAs with dummy value '0' b/c STAN doesn't like NAs:
  y_obs1[is.na(y_obs1)] <- 0
  
  
  min_y_obs2 <- min(dat_single2, na.rm=TRUE)
  max_y_obs2 <- max(dat_single2, na.rm=TRUE)
  
  y_obs2 <- dat_single2
  y_obs_bool2 <- 1*(!is.na(y_obs2))
  # replace NAs with dummy value '0' b/c STAN doesn't like NAs:
  y_obs2[is.na(y_obs2)] <- 0
  
  return(list(
    min_y_obs1=min_y_obs1, max_y_obs1=max_y_obs1, y_obs1=y_obs1, y_obs_bool1=y_obs_bool1,
    min_y_obs2=min_y_obs2, max_y_obs2=max_y_obs2, y_obs2=y_obs2, y_obs_bool2=y_obs_bool2
    ))
}



FE_priors_func <- function(N,
                           estimates, RSEs,
                           FE_names_per_traj, traj_names,
                           n_FEs_per_traj, n_trajs) {
  SD_estimates <- matrix(
    nrow=n_FEs_per_traj, ncol=n_trajs,
    dimnames=list(paste0(FE_names_per_traj, "(SD estimate)"), traj_names)
    )
  
  logN_mu_estimates <- matrix(
    nrow=n_FEs_per_traj, ncol=n_trajs,
    dimnames=list(paste0(FE_names_per_traj, "(logN_mu estimate)"), traj_names)
    )
  
  logN_sigma_estimates <- matrix(
    nrow=n_FEs_per_traj, ncol=n_trajs,
    dimnames=list(paste0(FE_names_per_traj, "(logN_sigma estimate)"), traj_names)
    )
  
  for (row_indx in 1:n_FEs_per_traj) {
    for (col_indx in 1:n_trajs) {
      target_RSE <- RSEs[row_indx, col_indx]
      target_mean <- estimates[row_indx, col_indx]
      target_SE <- target_RSE * target_mean / 100 
      target_SD <- target_SE * 3
      SD_estimates[row_indx, col_indx] <- round(target_SD, 2)
      
      v <- (target_SD / target_mean) ^ 2
      target_logN_sigma <- sqrt(log(v + 1))
      target_logN_mu <- log(target_mean) - (target_logN_sigma ^ 2 / 2)
      logN_mu_estimates[row_indx, col_indx] <- round(target_logN_mu, 2)
      logN_sigma_estimates[row_indx, col_indx] <- round(target_logN_sigma, 2)
    }
  }
  
  return(list(mean_estimates=estimates, # input
              SD_estimates=SD_estimates, # output
              logN_mu_estimates=logN_mu_estimates, # output
              logN_sigma_estimates=logN_sigma_estimates) # output
         )
}



missing_perc1_func <- function(y_mat_Lcols_amb) { # for: balanced1_sim_func, unbalanced1_sim_func
  N.J <- nrow(y_mat_Lcols_amb)
  L <- ncol(y_mat_Lcols_amb)
  
  missing_perc_by_muscle <- colSums(is.na(y_mat_Lcols_amb))/N.J # vector
  
  mat_NAs_for_each_row <- as.matrix(as.data.frame(
    table(rowSums(is.na(y_mat_Lcols_amb)))))
  class(mat_NAs_for_each_row) <- "numeric"
  for (l_indx in 0:L) {
    if (!l_indx %in% mat_NAs_for_each_row[,"Var1"]) {
      mat_NAs_for_each_row <- rbind(mat_NAs_for_each_row, c(l_indx,0))
    }
  }
  mat_NAs_for_each_row <- mat_NAs_for_each_row[
    order(mat_NAs_for_each_row[,"Var1"]),]
  mat_NAs_for_each_row <- cbind(mat_NAs_for_each_row,
                                mat_NAs_for_each_row[,"Freq"]/N.J)
  colnames(mat_NAs_for_each_row) <- c("NAs_for_each_row","Freq","%")
  missing_perc_by_rowNAs <- mat_NAs_for_each_row[
    mat_NAs_for_each_row[,"NAs_for_each_row"]!=0,"%"] # vector
  
  return(list(missing_perc_by_muscle=missing_perc_by_muscle,
              missing_perc_by_rowNAs=missing_perc_by_rowNAs))
}

missing_perc2_func <- function(y_list, L) { # for: unbalanced2_sim_func
  N <- length(y_list) # N == nrow(y_sim_full) == length(y_sim_full_unbal) == length(y_sim_wNA_unbal)
  ps <- sapply(1:N, function(i) length(y_list[[i]]))
  Js <- ps/L
  p_max <- max(ps)
  J_max <- max(Js)
  
  y_mat_Lcols_amb <- matrix(unlist(y_list), ncol=L, byrow=TRUE)
  
  unique_Js <- sort(unique(Js))
  
  missing_perc_by_muscle <- matrix(nrow=J_max, ncol=L)
  missing_perc_by_rowNAs <- matrix(nrow=J_max, ncol=L)
  for (j_indx in unique_Js) {
    p_indx <- j_indx*L
    y_list_temp <- Filter(function(x) length(x)==p_indx,
                          y_list)
    y_mat_Lcols_temp <- matrix(unlist(y_list_temp), ncol=L, byrow=TRUE)
    N.J_temp <- nrow(y_mat_Lcols_temp)
    
    missing_perc_by_muscle[j_indx,] <- colSums(is.na(y_mat_Lcols_temp)/N.J_temp)
    
    mat_NAs_for_each_row_temp <- as.matrix(as.data.frame(
      table(rowSums(is.na(y_mat_Lcols_temp)))))
    class(mat_NAs_for_each_row_temp) <- "numeric"
    for (l_indx in 0:L) {
      if (!l_indx %in% mat_NAs_for_each_row_temp[,"Var1"]) {
        mat_NAs_for_each_row_temp <- rbind(mat_NAs_for_each_row_temp, c(l_indx,0))
      }
    }
    mat_NAs_for_each_row_temp <- mat_NAs_for_each_row_temp[
      order(mat_NAs_for_each_row_temp[,"Var1"]),]
    mat_NAs_for_each_row_temp <- cbind(mat_NAs_for_each_row_temp,
                                       mat_NAs_for_each_row_temp[,"Freq"]/N.J_temp)
    colnames(mat_NAs_for_each_row_temp) <- c("NAs_for_each_row","Freq","%")
    missing_perc_by_rowNAs_temp <- mat_NAs_for_each_row_temp[
      mat_NAs_for_each_row_temp[,"NAs_for_each_row"]!=0,"%"]
    missing_perc_by_rowNAs[j_indx,] <- missing_perc_by_rowNAs_temp
  }
  
  return(list(missing_perc_by_muscle=missing_perc_by_muscle,
              missing_perc_by_rowNAs=missing_perc_by_rowNAs))
}

Js_dist_func <- function(Js, N) { # for: unbalanced1_sim_func, unbalanced2_sim_func
  mat_Js_dist <- as.matrix(as.data.frame(table(Js)))
  class(mat_Js_dist) <- "numeric"
  #for (j_indx in 1:J_max) {
  #  if (!j_indx %in% mat_Js_dist[,"Js"]) {
  #    mat_Js_dist <- rbind(mat_Js_dist, c(j_indx,0))
  #  }
  #}
  #mat_Js_dist <- mat_Js_dist[order(mat_Js_dist[,"Js"]),]
  mat_Js_dist <- cbind(mat_Js_dist, mat_Js_dist[,"Freq"]/N)
  colnames(mat_Js_dist) <- c("Js", "Freq", "%")
  perc_Js_dist <- mat_Js_dist[,"%"]
  
  return(list(mat_Js_dist=mat_Js_dist,
              perc_Js_dist=perc_Js_dist))
}



get_param_value <- function(RE_prefix, muscle_name, i, it, stan_mat, stan_mat_colnames) {
  RE_string <- paste0(RE_prefix, "_", muscle_name, "[", i, "]")
  
  # case 1
  if (RE_string %in% stan_mat_colnames) {
    return(as.numeric(stan_mat[it, RE_string]))
  }
  
  # case 2 ('else' functionality cuz of return above)
  FE_prefix <- sub("_", "_bar_", RE_prefix)
  FE_string <- paste0(FE_prefix, "_", muscle_name)
  if (FE_string %in% stan_mat_colnames) {
    return(as.numeric(stan_mat[it, FE_string]))
  }
  
  # case 3 ('else' functionality cuz of return above)s
  warning(paste("Neither", RE_string, "nor", FE_string, "found in Stan matrix."))
  return(NA)
}



stable_log_sum_exp <- function(x_vec) {
  mx <- max(x_vec)
  # If the max is -Inf, all elements are -Inf, so the result is -Inf
  if (mx == -Inf) {
    return(-Inf)
  }
  # log(sum(exp(x_vec))) = mx + log(sum(exp(x_vec - mx)))
  # Prevents overflox since largest value inside the exp() is 0
  return(mx + log(sum(exp(x_vec - mx))))
}

SMWT_function <- function(age, S0_F, DPT_F, gamma_F) {
  # numerically unstable (could overflow and produce NaN)
  S0_F * ( 1 - age ^ gamma_F / ( DPT_F ^ gamma_F + age ^ gamma_F ) )
}
SMWT_function0 <- function(age, target_SMWT, S0_F, DPT_F, gamma_F) {
  # numerically unstable (could overflow and produce NaN)
  target_SMWT - S0_F * ( 1 - age ^ gamma_F / ( DPT_F ^ gamma_F + age ^ gamma_F ) )
}
stable_SMWT_function0 <- function(age, target_SMWT, S0_F, log_DPT_F, gamma_F) {
  # NOT VECTORIZED (better for algo. due to performance)
  log_numer <- gamma_F * log(age)
  log_denomA <- gamma_F * log_DPT_F
  log_denomB <- log_numer
  log_denom <- stable_log_sum_exp(c(log_denomA, log_denomB))
  log_fraction <- log_numer - log_denom
  
  predicted_SMWT <- S0_F * (1 - exp(log_fraction))
  # return difference for 'uniroot'
  return(predicted_SMWT - target_SMWT) # ordering shouldn't matter
}
get_SMWT_trajectory <- function(age, S0_F, log_DPT_F, gamma_F) {
  # VECTORIZED (better for plotting, i.e. vector 'age')
  n_ages <- length(age)
  
  log_numer <- gamma_F * log(age)
  log_denomA <- rep(gamma_F * log_DPT_F, n_ages)
  log_denomB <- log_numer
  
  row_maxes <- pmax(log_denomA, log_denomB)
  log_denom <- row_maxes + log(
    exp(log_denomA - row_maxes) + exp(log_denomB - row_maxes)
    )
  
  log_fraction <- log_numer - log_denom
  
  predicted_SMWT <- S0_F * (1 - exp(log_fraction))
  return(predicted_SMWT)
}

muscle_function <- function(age, S0_T2, DP_T2, DPT_T2, gamma_T2) {
  S0_T2 + ( DP_T2 * age ^ gamma_T2 / ( DPT_T2 ^ gamma_T2 + age ^ gamma_T2 ) )
}
stable_muscle_function <- function(age, log_S0_T2, log_DP_T2, log_DPT_T2, gamma_T2) {
  # NOT VECTORIZED (better for algo. due to performance)
  log_numer <- gamma_T2 * log(age)
  log_denomA <- gamma_T2 * log_DPT_T2
  log_denomB <- log_numer
  log_denom <- stable_log_sum_exp(c(log_denomA, log_denomB))
  log_fraction <- log_numer - log_denom
  
  log_term1 <- log_S0_T2
  log_term2 <- log_DP_T2 + log_fraction
  log_predicted <- stable_log_sum_exp(c(log_term1, log_term2))
  predicted <- exp(log_predicted)
  return(predicted)
}
get_muscle_trajectory <- function(age, log_S0_T2, log_DP_T2, log_DPT_T2, gamma_T2) {
  # VECTORIZED (better for plotting, i.e. vector 'age')
  n_ages <- length(age)
  
  log_numer <- gamma_T2 * log(age)
  log_denomA <- rep(gamma_T2 * log_DPT_T2, n_ages)
  log_denomB <- log_numer
  
  row_maxes <- pmax(log_denomA, log_denomB)
  log_denom <- row_maxes + log(
    exp(log_denomA - row_maxes) + exp(log_denomB - row_maxes)
    )
  log_fraction <- log_numer - log_denom
  
  log_term1 <- rep(log_S0_T2, n_ages)
  log_term2 <- rep(log_DP_T2, n_ages) + log_fraction
  
  row_maxes2 <- pmax(log_term1, log_term2)
  log_predicted <- row_maxes2 + log(
    exp(log_term1 - row_maxes2) + exp(log_term2 - row_maxes2)
  )
  predicted <- exp(log_predicted)
  return(predicted)
}



calculate_average_matrix <- function(prob_matrix_list, indices) {
  prob_matrices_temp <- prob_matrix_list[indices]
  
  # Sum matrices in the list; treat any NAs as 0
  summed_matrix <- Reduce(function(m1, m2)
    m1 + ifelse(is.na(m2), 0, m2),
    prob_matrices_temp,
    init = matrix(0, L, L))
  
  # Divide by the number of matrices to get the average
  avg_matrix <- summed_matrix / length(prob_matrices_temp)
  
  diag(avg_matrix) <- NA
  return(avg_matrix)
}



SRM <- function(w,mu,Sigma) {
  return( c(t(w)%*%mu) / c(sqrt(t(w)%*%Sigma%*%w)) )
}
SRM_grad <- function(w,mu,Sigma) { # confirmed okay-looking in 'SRM_not_convex2.R'
  return( c((t(w)%*%Sigma%*%w)^{-0.5})*mu
          -c(t(w)%*%mu)*c((t(w)%*%Sigma%*%w)^{-1.5})*(Sigma%*%w) )
}
softmax <- function(w) {
  w1 <- w - max(w)
  exp(w1)/sum(exp(w1))
}
muscle_optim_func <- function(muscle_mu, muscle_Sigma,
                              optim_method="L-BFGS-B", grad_bool=TRUE, upper_bool=TRUE) {
  L <- nrow(muscle_Sigma)
  if (optim_method=="L-BFGS-B") {
    lower_optim_bound <- rep(0,L)
    upper_optim_bound <- rep(1,L)
  } else {
    lower_optim_bound <- -Inf
    upper_optim_bound <- Inf
  }
  
  if (grad_bool==TRUE) {
    if (upper_bool==TRUE) {
      w_summ <- optim(par=rep(1/L,L), fn=SRM, gr=SRM_grad,
                      mu=muscle_mu, Sigma=muscle_Sigma,
                      method=optim_method,
                      lower=lower_optim_bound, upper=upper_optim_bound,
                      control=list(fnscale=-1))
    } else {
      w_summ <- optim(par=rep(1/L,L), fn=SRM, gr=SRM_grad,
                      mu=muscle_mu, Sigma=muscle_Sigma,
                      method=optim_method,
                      lower=lower_optim_bound, 
                      control=list(fnscale=-1))
    }
  } else {
    if (upper_bool==TRUE) {
      w_summ <- optim(par=rep(1/L,L), fn=SRM,
                      mu=muscle_mu, Sigma=muscle_Sigma,
                      method=optim_method,
                      lower=lower_optim_bound, upper=upper_optim_bound,
                      control=list(fnscale=-1))
    } else {
      w_summ <- optim(par=rep(1/L,L), fn=SRM,
                      mu=muscle_mu, Sigma=muscle_Sigma,
                      method=optim_method,
                      lower=lower_optim_bound, 
                      control=list(fnscale=-1))
    }
  }
  
  w_orig <- w_summ$par
  if (optim_method=="L-BFGS-B") {
    w <- abs(w_orig)/sum(abs(w_orig))
  } else {
    w <- softmax(w_orig)
  }
  #SRMval <- SRM(w, muscle_mu, muscle_Sigma)
  
  return(list(w_orig=w_orig,
              w=w#,
              #SRMval=SRMval
  ))
}