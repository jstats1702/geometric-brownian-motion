source("~/Dropbox/PAPERS/projects/brownian/helper functions jumps.R")

# ==== Setup ====
M <- 1000
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
       tryCatch({
         dat <- sim_datasets_jump[[m]]
         y   <- log(dat$price)
         
         # Prior setup inside mcmc2
         fit <- mcmc2_with_priors(y = y, nsim = ns0, nburnin = nb0, prior_config = prior)
         
         mu_list[[m]]  <- fit$mu[-(1:nb0)]
         sig_list[[m]] <- fit$sig2[-(1:nb0)]
         
         cat("Prior:", prior, "| Dataset:", m, "– OK\n")
       }, error = function(e) {
         message(sprintf("⚠ Error with prior=%s, dataset=%d: %s", prior, m, e$message))
         # record a placeholder so your list indices stay aligned
         mu_list[[m]]  <- NA
         sig_list[[m]] <- NA
       })
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
       tryCatch({
         true_values <- log(sim_datasets_jump[[m]]$price)[-1]
         
         preds <- generate_gbm_predictions(
           th     = mu_list[[m]], 
           sig2   = sig_list[[m]], 
           y0     = true_values[1], 
           nsteps = 252
         )
         
         eval <- evaluate_predictions(preds, true_values)
         
         rmse_values[m] <- eval["RMSE"]
         mae_values[m]  <- eval["MAE"]
         
         cat("Dataset:", m, "– Evaluation OK\n")
       }, error = function(e) {
         message(sprintf("⚠ Error during evaluation of dataset %d: %s", m, e$message))
         rmse_values[m] <- NA
         mae_values[m]  <- NA
       })
     }
    
     eval_results_jump[[prior]] <- list(RMSE = mean(rmse_values, na.rm = T), MAE = mean(mae_values, na.rm = T))
}

# ==== Print Results ====
print("Evaluation results for GBM with jumps:")
print(eval_results_jump)
