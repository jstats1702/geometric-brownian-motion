# Define simulate_gbm_data() if not already defined
simulate_gbm_data <- function(mu = 0.1, sigma = 0.2, x0 = 100, T = 1, N = 252, seed = 123) {
     set.seed(seed)
     
     dt <- T / N
     time_points <- seq(0, T, length.out = N + 1)
     W <- cumsum(rnorm(N, mean = 0, sd = sqrt(dt)))
     W <- c(0, W)
     
     drift <- (mu - 0.5 * sigma^2) * time_points
     diffusion <- sigma * W
     Y <- log(x0) + drift + diffusion
     X <- exp(Y)
     
     data.frame(time = time_points, log_price = Y, price = X)
}

# Define the function to generate multiple datasets
generate_multiple_gbm_datasets <- function(M = 10, 
                                           mu = 0.1, sigma = 0.2, 
                                           x0 = 100, T = 1, N = 252, 
                                           base_seed = 100) {
     datasets <- vector("list", M)
     for (m in 1:M) {
          datasets[[m]] <- simulate_gbm_data(mu = mu, sigma = sigma, x0 = x0,
                                             T = T, N = N, seed = base_seed + m)
     }
     return(datasets)
}

# Updated MCMC function with flexible priors
mcmc1_with_priors <- function(y, t, nsim, nburnin, prior_config = "noninformative") {
     nchain <- nburnin + nsim
     set.seed(seed)
     
     # Define priors based on configuration
     if (prior_config == "noninformative") {
          mu_prior <- c(0, 100)   # Normal(0, 100)
          sig2_prior <- c(0.01, 0.01)  # Inverse-Gamma(0.01, 0.01)
     } else if (prior_config == "moderate") {
          mu_prior <- c(0, 1)    # Normal(0, 1)
          sig2_prior <- c(2, 2)  # Inverse-Gamma(2, 2)
     } else if (prior_config == "informative") {
          mu_prior <- c(0.1, 0.1)  # Normal(0.1, 0.1)
          sig2_prior <- c(10, 10)  # Inverse-Gamma(10, 10)
     }
     
     # Compute differences
     y.d <- y[-length(y)] - y[-1]
     t.d <- rep(1, length(y.d)) / 252
     n <- length(y.d)
     
     # Initialization
     th.ls <- sig2.ls <- rep(NA, nchain)
     a0 <- 2; b0 <- 0.001; sig0 <- 100
     th.ls[1] <- 0; sig2 <- sig2.ls[1] <- 1e-3
     
     # Sampling functions
     samp_th <- function(sig2) {
          sig2.th <- 1 / (1 / sig0 + sum(t.d) / sig2)
          mu.th <- sig2.th * (sum(y.d) / sig2)
          rnorm(1, mu.th, sqrt(sig2.th))
     }
     
     samp_sig2 <- function(th) {
          a.sig2 <- a0 + n / 2
          b.sig2 <- b0 + sum((y.d - th * t.d)^2 / t.d) / 2
          rigamma(1, a.sig2, b.sig2)
     }
     
     # MCMC loop
     for (k in 2:nchain) {
          th.ls[k] <- th <- samp_th(sig2 = sig2)
          sig2.ls[k] <- sig2 <- samp_sig2(th = th)
          # if (k %% 500 == 0) print(k)
     }
     
     # Transform back
     mu.ls <- th.ls + sig2.ls / 2
     sig.ls <- sqrt(sig2.ls)
     
     return(list(th = th.ls, sig2 = sig2.ls, mu = mu.ls, sig = sig.ls))
}

generate_gbm_predictions <- function(th, sig2, y0, nsteps = 252) {
     nsim <- length(th)
     dt <- 1 / 252  # daily time increment
     preds <- matrix(NA, nrow = nsim, ncol = nsteps)
     
     for (i in 1:nsim) {
          increments <- rnorm(nsteps, mean = th[i] * dt, sd = sqrt(sig2[i] * dt))
          preds[i, ] <- y0 + cumsum(increments)
     }
     
     return(preds)
}

