# ==== Setup ====
M <- 25
seed <- 1234
set.seed(seed)

nb0 <- 1000  # burn-in
ns0 <- 5000  # posterior samples

prior_configs <- c("noninformative", "moderate", "informative")
results_jump <- list()

# ==== Step 1: Generate datasets ====
sim_datasets_jump <- vector("list", M)
for (m in 1:M) {
     sim_datasets_jump[[m]] <- simulate_gbm_with_jumps(seed = seed + m)
}

# ==== Step 2: Fit jump model for each dataset and prior ====
for (prior in prior_configs) {
     mu_list <- vector("list", M)
     sig_list <- vector("list", M)
     
     for (m in 1:M) {
          dat <- sim_datasets_jump[[m]]
          y <- log(dat$price)
          
          # Prior setup inside mcmc2
          fit <- mcmc2_with_priors(y = y, nsim = ns0, nburnin = nb0, prior_config = prior)
          
          mu_list[[m]] <- fit$mu[-(1:nb0)]
          sig_list[[m]] <- fit$sig2[-(1:nb0)]
          
          cat("Prior:", prior, "| Dataset:", m, "\n")
     }
     
     results_jump[[prior]] <- list(mu = mu_list, sigma = sig_list)
}

# ==== Step 4: Evaluate performance ====
eval_results_jump <- list()

for (prior in prior_configs) {
     mu_list <- results_jump[[prior]]$mu
     sig_list <- results_jump[[prior]]$sigma
     
     rmse_values <- numeric(M)
     mae_values <- numeric(M)
     
     for (m in 1:M) {
          true_values <- log(sim_datasets_jump[[m]]$price)[-1]
          preds <- generate_gbm_predictions(th = mu_list[[m]], sig2 = sig_list[[m]], y0 = true_values[1], nsteps = 252)
          eval <- evaluate_predictions(preds, true_values)
          
          rmse_values[m] <- eval["RMSE"]
          mae_values[m] <- eval["MAE"]
     }
     
     eval_results_jump[[prior]] <- list(RMSE = mean(rmse_values), MAE = mean(mae_values))
}

# ==== Print Results ====
print("Evaluation results for GBM with jumps:")
print(eval_results_jump)
