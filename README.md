# dmd_project2
Code, posterior summaries, model diagnostics, plots for the paper: `A novel approach to optimize composite biomarkers for continuous measures of clinical function: Application to muscle imaging biomarkers in Duchenne muscular dystrophy`

### `code1/...`
- `stan_fit_all1.R` executes Hamiltonian Monte Carlo (HMC) with No-U-Turn Sampler (NUTS) implemented in `Stan` software for our joint nonlinear mixed-effects (NLME) model of DMD data (incomplete longitudinal trajectories of 6MWD and MRI-T2 of five leg muscles) -- see Section 3.1 of paper; also contains code for `stan_results` outputs folder
- `DPTcorr_lessREs_TjointNC_IB1.stan` specifies our joint NLME model in `Stan` language
- `functions1.R` contains various helper functions
- `grid_it_weights_joint.R` executes main algorithm for computing posterior samples of weights of one-year MRI-T2 changes of five leg muscles across a fine grid of 6MWD values -- see Section 4.1 of paper; also contains code for `main_plots` outputs folder

### `outputs1/...`
- `stan_results`: Posterior summaries and model diagnostics of our joint nonlinear mixed-effects (NLME) model implemented in `Stan` software
- `main_plots`: Trajectory plots (data, model), histograms (data), solved ages of weights algorithm, average posterior dominance probabilities of muscle weights, main weights plot (one-year MRI-T2 changes of five leg muscles across a fine grid of 6MWD values)
