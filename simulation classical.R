setwd("~/Dropbox/PAPERS/projects/brownian")
source("~/Dropbox/PAPERS/projects/brownian/helper functions classical.R")

pckgs<-c("moments","MASS","pscl","mnormt")

lapply(pckgs,require,character.only=T)

M <- 1000
seed <- 1234

# Step 1: Generate M datasets
sim_datasets <- generate_multiple_gbm_datasets(M)

# Step 2: Fit models for each prior and dataset
prior_configs <- c("noninformative", "moderate", "informative")

results <- list()

# MCMC settings
nb0 <- 1000  # burn-in
ns0 <- 5000  # number of post-burn-in samples

for (prior in prior_configs) {
     mu_list <- vector("list", M)
     sig_list <- vector("list", M)
     
     for (m in 1:M) {
          # Fit the model for each dataset with the current prior
          dat <- sim_datasets[[m]]
          y <- log(dat$price)
          t <- dat$time
          
          fit <- mcmc1_with_priors(y = y, t = t, nsim = ns0, nburnin = nb0, prior_config = prior)
          
          # Store the posterior samples for mu and sigma
          mu_list[[m]] <- fit$mu[-(1:nb0)]
          sig_list[[m]] <- fit$sig2[-(1:nb0)]
          
          cat("Prior:", prior, "| Dataset:", m, "\n")
     }
     
     # Store results for each prior
     results[[prior]] <- list(mu = mu_list, sigma = sig_list)
}

# Step 3: Evaluate predictive performance (RMSE and MAE)
evaluate_predictions <- function(predictions, true_values) {
     pred_mean <- colMeans(predictions)
     rmse <- sqrt(mean((pred_mean - true_values)^2))
     mae <- mean(abs(pred_mean - true_values))
     return(c(RMSE = rmse, MAE = mae))
}

# Evaluate predictions and summarize
eval_results <- list()

for (prior in prior_configs) {
     mu_list <- results[[prior]]$mu
     sig_list <- results[[prior]]$sigma
     
     rmse_values <- numeric(M)
     mae_values <- numeric(M)
     
     for (m in 1:M) {
          # Get the true values for this dataset (use log-prices as true values)
          true_values <- log(sim_datasets[[m]]$price)[-1]
          
          # Generate predictions for this dataset
          preds <- generate_gbm_predictions(th = mu_list[[m]], sig2 = sig_list[[m]], y0 = true_values[1], nsteps = 252)
          
          # Evaluate prediction performance
          eval <- evaluate_predictions(preds, true_values)
          
          # Store RMSE and MAE
          rmse_values[m] <- eval["RMSE"]
          mae_values[m] <- eval["MAE"]
     }
     
     # Store evaluation results for the current prior
     eval_results[[prior]] <- list(RMSE = mean(rmse_values), MAE = mean(mae_values))
}

# Print the results
print(eval_results)
