data {
  int N; // number of subjects
  real min_age;
  real max_age;

  int n_obs_6MWD;
  vector[n_obs_6MWD] y_obs_6MWD; // OMIT BOUNDS BECAUSE CENSORING
  array[n_obs_6MWD] int<lower=1, upper=N> ID_6MWD;
  vector<lower=min_age, upper=max_age> [n_obs_6MWD] ages_6MWD;

  int n_obs_BFLH;
  real min_obs_BFLH;
  real max_obs_BFLH;
  vector<lower=min_obs_BFLH, upper=max_obs_BFLH> [n_obs_BFLH] y_obs_BFLH;
  array[n_obs_BFLH] int<lower=1, upper=N> ID_BFLH;
  vector<lower=min_age, upper=max_age> [n_obs_BFLH] ages_BFLH;

  int n_obs_GRA;
  real min_obs_GRA;
  real max_obs_GRA;
  vector<lower=min_obs_GRA, upper=max_obs_GRA> [n_obs_GRA] y_obs_GRA;
  array[n_obs_GRA] int<lower=1, upper=N> ID_GRA;
  vector<lower=min_age, upper=max_age> [n_obs_GRA] ages_GRA;

  int n_obs_MG;
  real min_obs_MG;
  real max_obs_MG;
  vector<lower=min_obs_MG, upper=max_obs_MG> [n_obs_MG] y_obs_MG;
  array[n_obs_MG] int<lower=1, upper=N> ID_MG;
  vector<lower=min_age, upper=max_age> [n_obs_MG] ages_MG;

  int n_obs_SOL;
  real min_obs_SOL;
  real max_obs_SOL;
  vector<lower=min_obs_SOL, upper=max_obs_SOL> [n_obs_SOL] y_obs_SOL;
  array[n_obs_SOL] int<lower=1, upper=N> ID_SOL;
  vector<lower=min_age, upper=max_age> [n_obs_SOL] ages_SOL;

  int n_obs_VL;
  real min_obs_VL;
  real max_obs_VL;
  vector<lower=min_obs_VL, upper=max_obs_VL> [n_obs_VL] y_obs_VL;
  array[n_obs_VL] int<lower=1, upper=N> ID_VL;
  vector<lower=min_age, upper=max_age> [n_obs_VL] ages_VL;
}

parameters {
  // DP_6MWD has no FE or RE (set as 1 in paper) -- FE include => 7242/30000 div, 22500/30000 max_td hit; RE include => 30000/30000 max_td hit
  // DP_GRA has no RE (missing for all pairwise joint models with or without covariates) -- 77.97% shrinkage, 55.20% RSE

  // in order to construct mixed-error
  real<lower=0> tau_const_6MWD;
  real<lower=0> tau_prop_6MWD;

  // in order to construct proportional error
  real<lower=0> tau_BFLH;
  real<lower=0> tau_GRA;
  real<lower=0> tau_MG;
  real<lower=0> tau_SOL;
  real<lower=0> tau_VL;


  
  real log_bar_DPT_6MWD; // fixed eff
  real log_bar_gamma_6MWD; // fixed eff (min=0, max=infty)
  real log_bar_S0_6MWD; // fixed eff (extrapolated value at age 0) (prob no upper limit for 6MWD)

  real<lower=0> omega_S0_6MWD; // s.d. for BSV random effect

  vector[N] eta_S0_6MWD_raw; // BSV random effect



  real log_bar_DP_BFLH; // fixed eff (min_increase=0; max_increase=max(muscle)-min(muscle))
  real log_bar_DPT_BFLH; // fixed eff
  real log_bar_gamma_BFLH; // fixed eff (min=0, max=infty)
  real log_bar_S0_BFLH; // fixed eff (extrapolated value at age 0) (prob no upper limit for 6MWD)

  real<lower=0> omega_DP_BFLH; // s.d. for BSV random effect
  real<lower=0> omega_S0_BFLH; // s.d. for BSV random effect

  vector[N] eta_DP_BFLH_raw; // BSV random effect
  vector[N] eta_S0_BFLH_raw; // BSV random effect



  real log_bar_DP_GRA; // fixed eff (min_increase=0; max_increase=max(muscle)-min(muscle))
  real log_bar_DPT_GRA; // fixed eff
  real log_bar_gamma_GRA; // fixed eff (min=0, max=infty)
  real log_bar_S0_GRA; // fixed eff (extrapolated value at age 0) (prob no upper limit for 6MWD)

  real<lower=0> omega_S0_GRA; // s.d. for BSV random effect

  vector[N] eta_S0_GRA_raw; // BSV random effect



  real log_bar_DP_MG; // fixed eff (min_increase=0; max_increase=max(muscle)-min(muscle))
  real log_bar_DPT_MG; // fixed eff
  real log_bar_gamma_MG; // fixed eff (min=0, max=infty)
  real log_bar_S0_MG; // fixed eff (extrapolated value at age 0) (prob no upper limit for 6MWD)

  real<lower=0> omega_DP_MG; // s.d. for BSV random effect
  real<lower=0> omega_gamma_MG; // s.d. for BSV random effect
  real<lower=0> omega_S0_MG; // s.d. for BSV random effect

  vector[N] eta_DP_MG_raw; // BSV random effect
  vector[N] eta_gamma_MG_raw; // BSV random effect
  vector[N] eta_S0_MG_raw; // BSV random effect



  real log_bar_DP_SOL; // fixed eff (min_increase=0; max_increase=max(muscle)-min(muscle))
  real log_bar_DPT_SOL; // fixed eff
  real log_bar_gamma_SOL; // fixed eff (min=0, max=infty)
  real log_bar_S0_SOL; // fixed eff (extrapolated value at age 0) (prob no upper limit for 6MWD)

  real<lower=0> omega_DP_SOL; // s.d. for BSV random effect
  real<lower=0> omega_gamma_SOL; // s.d. for BSV random effect
  real<lower=0> omega_S0_SOL; // s.d. for BSV random effect

  vector[N] eta_DP_SOL_raw; // BSV random effect
  vector[N] eta_gamma_SOL_raw; // BSV random effect
  vector[N] eta_S0_SOL_raw; // BSV random effect



  real log_bar_DP_VL; // fixed eff (min_increase=0; max_increase=max(muscle)-min(muscle))
  real log_bar_DPT_VL; // fixed eff
  real log_bar_gamma_VL; // fixed eff (min=0, max=infty)
  real log_bar_S0_VL; // fixed eff (extrapolated value at age 0) (prob no upper limit for 6MWD)

  real<lower=0> omega_DP_VL; // s.d. for BSV random effect
  real<lower=0> omega_gamma_VL; // s.d. for BSV random effect
  real<lower=0> omega_S0_VL; // s.d. for BSV random effect

  vector[N] eta_DP_VL_raw; // BSV random effect
  vector[N] eta_gamma_VL_raw; // BSV random effect
  vector[N] eta_S0_VL_raw; // BSV random effect



  vector<lower=0>[6] s;
  cholesky_factor_corr[6] R_chol;
  matrix[6, N] z;
}

transformed parameters {
  vector[N] eta_S0_6MWD = omega_S0_6MWD * eta_S0_6MWD_raw; // BSV random effect

  vector[N] eta_DP_BFLH = omega_DP_BFLH * eta_DP_BFLH_raw; // BSV random effect
  vector[N] eta_S0_BFLH = omega_S0_BFLH * eta_S0_BFLH_raw; // BSV random effect

  vector[N] eta_S0_GRA = omega_S0_GRA * eta_S0_GRA_raw; // BSV random effect

  vector[N] eta_DP_MG = omega_DP_MG * eta_DP_MG_raw; // BSV random effect
  vector[N] eta_gamma_MG = omega_gamma_MG * eta_gamma_MG_raw; // BSV random effect
  vector[N] eta_S0_MG = omega_S0_MG * eta_S0_MG_raw; // BSV random effect

  vector[N] eta_DP_SOL = omega_DP_SOL * eta_DP_SOL_raw; // BSV random effect
  vector[N] eta_gamma_SOL = omega_gamma_SOL * eta_gamma_SOL_raw; // BSV random effect
  vector[N] eta_S0_SOL = omega_S0_SOL * eta_S0_SOL_raw; // BSV random effect

  vector[N] eta_DP_VL = omega_DP_VL * eta_DP_VL_raw; // BSV random effect
  vector[N] eta_gamma_VL = omega_gamma_VL * eta_gamma_VL_raw; // BSV random effect
  vector[N] eta_S0_VL = omega_S0_VL * eta_S0_VL_raw; // BSV random effect


  corr_matrix[6] R = R_chol * R_chol';
  matrix[6, N] eta_mat = diag_pre_multiply(s, R_chol) * z; // eta_mat ~ N(0, diag(s) * R_chol * R_chol' * diag(s))

  
  vector[N] log_S0_6MWD = log_bar_S0_6MWD + eta_S0_6MWD;
  vector[N] log_DPT_6MWD = log_bar_DPT_6MWD + (eta_mat[1,])';

  vector[N] log_S0_BFLH = log_bar_S0_BFLH + eta_S0_BFLH;
  vector[N] log_DP_BFLH = log_bar_DP_BFLH + eta_DP_BFLH;
  vector[N] log_DPT_BFLH = log_bar_DPT_BFLH + (eta_mat[2,])';

  vector[N] log_S0_GRA = log_bar_S0_GRA + eta_S0_GRA;
  vector[N] log_DPT_GRA = log_bar_DPT_GRA + (eta_mat[3,])';

  vector[N] log_S0_MG = log_bar_S0_MG + eta_S0_MG;
  vector[N] log_DP_MG = log_bar_DP_MG + eta_DP_MG;
  vector[N] log_DPT_MG = log_bar_DPT_MG + (eta_mat[4,])';
  vector[N] log_gamma_MG = log_bar_gamma_MG + eta_gamma_MG;

  vector[N] log_S0_SOL = log_bar_S0_SOL + eta_S0_SOL;
  vector[N] log_DP_SOL = log_bar_DP_SOL + eta_DP_SOL;
  vector[N] log_DPT_SOL = log_bar_DPT_SOL + (eta_mat[5,])';
  vector[N] log_gamma_SOL = log_bar_gamma_SOL + eta_gamma_SOL;

  vector[N] log_S0_VL = log_bar_S0_VL + eta_S0_VL;
  vector[N] log_DP_VL = log_bar_DP_VL + eta_DP_VL;
  vector[N] log_DPT_VL = log_bar_DPT_VL + (eta_mat[6,])';
  vector[N] log_gamma_VL = log_bar_gamma_VL + eta_gamma_VL;



  //vector[n_obs_6MWD] y_pred_6MWD = exp(log_bar_S0_6MWD) .* ( 1 - ages_6MWD .^ exp(log_gamma_6MWD_long) ./ ( exp(log_DPT_6MWD_long) .^ exp(log_gamma_6MWD_long) + ages_6MWD .^ exp(log_gamma_6MWD_long) ) );
  vector[n_obs_6MWD] y_pred_6MWD;
  { // local block; all temporary variables
    vector[n_obs_6MWD] log_S0_6MWD_long = log_S0_6MWD[ID_6MWD];
    vector[n_obs_6MWD] log_DPT_6MWD_long = log_DPT_6MWD[ID_6MWD];
    real gamma_6MWD_fixed = exp(log_bar_gamma_6MWD);
    
    // Safely calculate the fraction N/D = ages^gamma / (DPT^gamma + ages^gamma) on the log scale
    vector[n_obs_6MWD] log_numer = gamma_6MWD_fixed .* log(ages_6MWD);
    vector[n_obs_6MWD] log_denomA = gamma_6MWD_fixed .* log_DPT_6MWD_long;
    vector[n_obs_6MWD] log_denomB = log_numer;
    vector[n_obs_6MWD] log_denom = log_sum_exp(log_denomA, log_denomB);
    vector[n_obs_6MWD] log_fraction = log_numer - log_denom;

    // exponentiating at end = more numerically stable
    y_pred_6MWD = exp(log_S0_6MWD_long) .* (1 - exp(log_fraction));
  }


  //vector[n_obs_BFLH] y_pred_BFLH = exp(log_S0_BFLH_long) + ( exp(log_DP_BFLH_long) .* ages_BFLH .^ exp(log_bar_gamma_BFLH) ./ ( exp(log_DPT_BFLH_long) .^ exp(log_bar_gamma_BFLH) + ages_BFLH .^ exp(log_bar_gamma_BFLH) ) );
  vector[n_obs_BFLH] log_y_pred_BFLH;
  { // local block; all temporary variables
    vector[n_obs_BFLH] log_S0_BFLH_long = log_S0_BFLH[ID_BFLH];
    vector[n_obs_BFLH] log_DP_BFLH_long = log_DP_BFLH[ID_BFLH];
    vector[n_obs_BFLH] log_DPT_BFLH_long = log_DPT_BFLH[ID_BFLH];
    real gamma_BFLH_fixed = exp(log_bar_gamma_BFLH);
    
    // Safely calculate the fraction N/D = ages^gamma / (DPT^gamma + ages^gamma) on the log scale
    vector[n_obs_BFLH] log_numer = gamma_BFLH_fixed .* log(ages_BFLH);
    vector[n_obs_BFLH] log_denomA = gamma_BFLH_fixed .* log_DPT_BFLH_long;
    vector[n_obs_BFLH] log_denomB = log_numer;
    vector[n_obs_BFLH] log_denom = log_sum_exp(log_denomA, log_denomB);
    vector[n_obs_BFLH] log_fraction = log_numer - log_denom;

    // exponentiating at end = more numerically stable
    vector[n_obs_BFLH] log_term1 = log_S0_BFLH_long;
    vector[n_obs_BFLH] log_term2 = log_DP_BFLH_long + log_fraction;
    vector[n_obs_BFLH] log_y_pred = log_sum_exp(log_term1, log_term2);

    log_y_pred_BFLH = log_y_pred;
  }


  //vector[n_obs_GRA] y_pred_GRA = exp(log_S0_GRA_long) + ( exp(log_bar_DP_GRA) .* ages_GRA .^ exp(log_bar_gamma_GRA) ./ ( exp(log_DPT_GRA_long) .^ exp(log_bar_gamma_GRA) + ages_GRA .^ exp(log_bar_gamma_GRA) ) );
  vector[n_obs_GRA] log_y_pred_GRA;
  { // local block; all temporary variables
    vector[n_obs_GRA] log_S0_GRA_long = log_S0_GRA[ID_GRA];
    vector[n_obs_GRA] log_DPT_GRA_long = log_DPT_GRA[ID_GRA];
    real gamma_GRA_fixed = exp(log_bar_gamma_GRA);
    
    // Safely calculate the fraction N/D = ages^gamma / (DPT^gamma + ages^gamma) on the log scale
    vector[n_obs_GRA] log_numer = gamma_GRA_fixed .* log(ages_GRA);
    vector[n_obs_GRA] log_denomA = gamma_GRA_fixed .* log_DPT_GRA_long;
    vector[n_obs_GRA] log_denomB = log_numer;
    vector[n_obs_GRA] log_denom = log_sum_exp(log_denomA, log_denomB);
    vector[n_obs_GRA] log_fraction = log_numer - log_denom;
    
    // exponentiating at end = more numerically stable
    vector[n_obs_GRA] log_term1 = log_S0_GRA_long;
    vector[n_obs_GRA] log_term2 = log_bar_DP_GRA + log_fraction;
    vector[n_obs_GRA] log_y_pred = log_sum_exp(log_term1, log_term2);

    log_y_pred_GRA = log_y_pred;
  }


  //vector[n_obs_MG] y_pred_MG = exp(log_S0_MG_long) + ( exp(log_DP_MG_long) .* ages_MG .^ exp(log_gamma_MG_long) ./ ( exp(log_DPT_MG_long) .^ exp(log_gamma_MG_long) + ages_MG .^ exp(log_gamma_MG_long) ) );
  vector[n_obs_MG] log_y_pred_MG;
  { // local block; all temporary variables
    vector[n_obs_MG] log_S0_MG_long = log_S0_MG[ID_MG];
    vector[n_obs_MG] log_DP_MG_long = log_DP_MG[ID_MG];
    vector[n_obs_MG] log_DPT_MG_long = log_DPT_MG[ID_MG];
    vector[n_obs_MG] gamma_MG_long = exp(log_gamma_MG[ID_MG]);

    // Safely calculate the fraction N/D = ages^gamma / (DPT^gamma + ages^gamma) on the log scale
    vector[n_obs_MG] log_numer = gamma_MG_long .* log(ages_MG);
    vector[n_obs_MG] log_denomA = gamma_MG_long .* log_DPT_MG_long;
    vector[n_obs_MG] log_denomB = log_numer;
    vector[n_obs_MG] log_denom = log_sum_exp(log_denomA, log_denomB);
    vector[n_obs_MG] log_fraction = log_numer - log_denom;
    
    // exponentiating at end = more numerically stable
    vector[n_obs_MG] log_term1 = log_S0_MG_long;
    vector[n_obs_MG] log_term2 = log_DP_MG_long + log_fraction;
    vector[n_obs_MG] log_y_pred = log_sum_exp(log_term1, log_term2);
    
    log_y_pred_MG = log_y_pred;
  }


  //vector[n_obs_SOL] y_pred_SOL = exp(log_bar_S0_SOL) + ( exp(log_DP_SOL_long) .* ages_SOL .^ exp(log_gamma_SOL_long) ./ ( exp(log_DPT_SOL_long) .^ exp(log_gamma_SOL_long) + ages_SOL .^ exp(log_gamma_SOL_long) ) );
  vector[n_obs_SOL] log_y_pred_SOL;
  { // local block; all temporary variables
    vector[n_obs_SOL] log_S0_SOL_long = log_S0_SOL[ID_SOL];
    vector[n_obs_SOL] log_DP_SOL_long = log_DP_SOL[ID_SOL];
    vector[n_obs_SOL] log_DPT_SOL_long = log_DPT_SOL[ID_SOL];
    vector[n_obs_SOL] gamma_SOL_long = exp(log_gamma_SOL[ID_SOL]);

    // Safely calculate the fraction N/D = ages^gamma / (DPT^gamma + ages^gamma) on the log scale
    vector[n_obs_SOL] log_numer = gamma_SOL_long .* log(ages_SOL);
    vector[n_obs_SOL] log_denomA = gamma_SOL_long .* log_DPT_SOL_long;
    vector[n_obs_SOL] log_denomB = log_numer;
    vector[n_obs_SOL] log_denom = log_sum_exp(log_denomA, log_denomB);
    vector[n_obs_SOL] log_fraction = log_numer - log_denom;

    // exponentiating at end = more numerically stable
    vector[n_obs_SOL] log_term1 = log_S0_SOL_long;
    vector[n_obs_SOL] log_term2 = log_DP_SOL_long + log_fraction;
    vector[n_obs_SOL] log_y_pred = log_sum_exp(log_term1, log_term2);

    log_y_pred_SOL = log_y_pred;
  }


  //vector[n_obs_VL] y_pred_VL = exp(log_S0_VL_long) + ( exp(log_bar_DP_VL) .* ages_VL .^ exp(log_gamma_VL_long) ./ ( exp(log_DPT_VL_long) .^ exp(log_gamma_VL_long) + ages_VL .^ exp(log_gamma_VL_long) ) );
  vector[n_obs_VL] log_y_pred_VL;
  { // local block; all temporary variables
    vector[n_obs_VL] log_S0_VL_long = log_S0_VL[ID_VL];
    vector[n_obs_VL] log_DP_VL_long = log_DP_VL[ID_VL];
    vector[n_obs_VL] log_DPT_VL_long = log_DPT_VL[ID_VL];
    vector[n_obs_VL] gamma_VL_long = exp(log_gamma_VL[ID_VL]);
    
    // Safely calculate the fraction N/D = ages^gamma / (DPT^gamma + ages^gamma) on the log scale
    vector[n_obs_VL] log_numer = gamma_VL_long .* log(ages_VL);
    vector[n_obs_VL] log_denomA = gamma_VL_long .* log_DPT_VL_long;
    vector[n_obs_VL] log_denomB = log_numer;
    vector[n_obs_VL] log_denom = log_sum_exp(log_denomA, log_denomB);
    vector[n_obs_VL] log_fraction = log_numer - log_denom;
    
    // exponentiating at end = more numerically stable
    vector[n_obs_VL] log_term1 = log_S0_VL_long;
    vector[n_obs_VL] log_term2 = log_DP_VL_long + log_fraction;
    vector[n_obs_VL] log_y_pred = log_sum_exp(log_term1, log_term2);

    log_y_pred_VL = log_y_pred;
  }
}

model {
  tau_const_6MWD ~ normal(0, 5);
  tau_prop_6MWD ~ normal(0, 1);
  
  tau_BFLH ~ normal(0, 1);
  tau_GRA ~ normal(0, 1);
  tau_MG ~ normal(0, 1);
  tau_SOL ~ normal(0, 1);
  tau_VL ~ normal(0, 1);

  to_vector(z) ~ std_normal();

  R_chol ~ lkj_corr_cholesky(1);

  s[1] ~ normal(0, 1);
  s[2] ~ normal(0, 1);
  s[3] ~ normal(0, 1);
  s[4] ~ normal(0, 1);
  s[5] ~ normal(0, 1);
  s[6] ~ normal(0, 1);

  omega_S0_6MWD ~ normal(0,1);

  omega_DP_BFLH ~ normal(0, 1);
  omega_S0_BFLH ~ normal(0, 1);

  omega_S0_GRA ~ normal(0, 1);

  omega_DP_MG ~ normal(0, 1);
  omega_gamma_MG ~ normal(0, 1);
  omega_S0_MG ~ normal(0, 1);

  omega_DP_SOL ~ normal(0, 1);
  omega_gamma_SOL ~ normal(0, 1);
  omega_S0_SOL ~ normal(0, 1);

  omega_DP_VL ~ normal(0, 1);
  omega_gamma_VL ~ normal(0, 1);
  omega_S0_VL ~ normal(0, 1);


  log_bar_DPT_6MWD ~ normal(2.51, 0.08);
  log_bar_gamma_6MWD ~ normal(3.94, 0.33);
  log_bar_S0_6MWD ~ normal(5.90, 0.05);

  eta_S0_6MWD_raw ~ std_normal();


  log_bar_DP_BFLH ~ normal(3.69, 0.10);
  log_bar_DPT_BFLH ~ normal(2.32, 0.08);
  log_bar_gamma_BFLH ~ normal(2.06, 0.15);
  log_bar_S0_BFLH ~ normal(3.69, 0.05);

  eta_DP_BFLH_raw ~ std_normal();
  eta_S0_BFLH_raw ~ std_normal();


  log_bar_DP_GRA ~ normal(3.29, 0.31);
  log_bar_DPT_GRA ~ normal(2.66, 0.11);
  log_bar_gamma_GRA ~ normal(2.20, 0.38);
  log_bar_S0_GRA ~ normal(3.65, 0.04);

  eta_S0_GRA_raw ~ std_normal();


  log_bar_DP_MG ~ normal(3.36, 0.26);
  log_bar_DPT_MG ~ normal(2.61, 0.16);
  log_bar_gamma_MG ~ normal(1.78, 0.37);
  log_bar_S0_MG ~ normal(3.66, 0.05);

  eta_DP_MG_raw ~ std_normal();
  eta_gamma_MG_raw ~ std_normal();
  eta_S0_MG_raw ~ std_normal();
  

  log_bar_DP_SOL ~ normal(3.52, 0.19);
  log_bar_DPT_SOL ~ normal(2.61, 0.11);
  log_bar_gamma_SOL ~ normal(1.44, 0.34);
  log_bar_S0_SOL ~ normal(3.64, 0.05);

  eta_DP_SOL_raw ~ std_normal();
  eta_gamma_SOL_raw ~ std_normal();
  eta_S0_SOL_raw ~ std_normal();
  

  log_bar_DP_VL ~ normal(3.59, 0.14);
  log_bar_DPT_VL ~ normal(2.35, 0.08);
  log_bar_gamma_VL ~ normal(2.12, 0.17);
  log_bar_S0_VL ~ normal(3.70, 0.05);

  eta_DP_VL_raw ~ std_normal();
  eta_gamma_VL_raw ~ std_normal();
  eta_S0_VL_raw ~ std_normal();

  for (i in 1:n_obs_6MWD) {
    // mixed-error s.d.
    real current_sd_6MWD = sqrt(square(y_pred_6MWD[i] * tau_prop_6MWD) + square(tau_const_6MWD));
    if (y_obs_6MWD[i] <= 10.0) {
      // Left-Censored data likelihood for values at the floor (obs values=10)
      target += normal_lcdf(10.0 | y_pred_6MWD[i], current_sd_6MWD);
    } else {
      // Regular observation likelihood
      target += normal_lpdf(y_obs_6MWD[i] | y_pred_6MWD[i], current_sd_6MWD)
              - normal_lccdf(0 | y_pred_6MWD[i], current_sd_6MWD);
    }
  }

  y_obs_BFLH ~ lognormal(log_y_pred_BFLH, tau_BFLH);
  y_obs_GRA ~ lognormal(log_y_pred_GRA, tau_GRA);
  y_obs_MG ~ lognormal(log_y_pred_MG, tau_MG);
  y_obs_SOL ~ lognormal(log_y_pred_SOL, tau_SOL);
  y_obs_VL ~ lognormal(log_y_pred_VL, tau_VL);
}

generated quantities {
  vector[n_obs_6MWD] log_lik_6MWD;
  vector[n_obs_BFLH] log_lik_BFLH;
  vector[n_obs_GRA] log_lik_GRA;
  vector[n_obs_MG] log_lik_MG;
  vector[n_obs_SOL] log_lik_SOL;
  vector[n_obs_VL] log_lik_VL;

  // Use explicit loops to calculate the pointwise log-likelihood for each outcome.
  // type-mismatch error when I do something like:
  //log_lik_6MWD = normal_lpdf(y_obs_6MWD | y_pred_6MWD, y_pred_6MWD .* tau_6MWD) - normal_lccdf(0 | y_pred_6MWD, y_pred_6MWD .* tau_6MWD);

  for (i in 1:n_obs_6MWD) {
    real current_sd_6MWD = sqrt(square(y_pred_6MWD[i] * tau_prop_6MWD) + square(tau_const_6MWD));
    if (y_obs_6MWD[i] <= 10.0) {
      log_lik_6MWD[i] = normal_lcdf(10.0 | y_pred_6MWD[i], current_sd_6MWD);
    } else {
      // For WAIC calculations, we need to be consistent with muscle outcomes being truncated at zero
      // when calculating the probability of observed data point given the model
      log_lik_6MWD[i] = normal_lpdf(y_obs_6MWD[i] | y_pred_6MWD[i], current_sd_6MWD)
                        - normal_lccdf(0 | y_pred_6MWD[i], current_sd_6MWD); 
    }
  }

  for (i in 1:n_obs_BFLH) {
    log_lik_BFLH[i] = lognormal_lpdf(y_obs_BFLH[i] | log_y_pred_BFLH[i], tau_BFLH);
  }

  for (i in 1:n_obs_GRA) {
    log_lik_GRA[i] = lognormal_lpdf(y_obs_GRA[i] | log_y_pred_GRA[i], tau_GRA);
  }

  for (i in 1:n_obs_MG) {
    log_lik_MG[i] = lognormal_lpdf(y_obs_MG[i] | log_y_pred_MG[i], tau_MG);
  }

  for (i in 1:n_obs_SOL) {
    log_lik_SOL[i] = lognormal_lpdf(y_obs_SOL[i] | log_y_pred_SOL[i], tau_SOL);
  }

  for (i in 1:n_obs_VL) {
    log_lik_VL[i] = lognormal_lpdf(y_obs_VL[i] | log_y_pred_VL[i], tau_VL);
  }
}
