# ==== Function: Simulate GBM with Jumps ====
simulate_gbm_with_jumps <- function(mu = 0.1, sigma = 0.2, 
                                    lambda = 0.02, mu_z = -0.05, sigma_z = 0.1,
                                    x0 = 100, T = 1, N = 252, seed = 123) {
     set.seed(seed)
     dt <- T / N
     time_points <- seq(0, T, length.out = N + 1)
     
     W <- cumsum(rnorm(N, mean = 0, sd = sqrt(dt)))
     W <- c(0, W)
     
     J <- rpois(N, lambda * dt)
     Z <- rnorm(N, mean = mu_z, sd = sigma_z)
     jumps <- cumsum(J * Z)
     
     drift <- (mu - 0.5 * sigma^2) * time_points
     diffusion <- sigma * W
     
     Y <- log(x0) + drift + diffusion + c(0, jumps)
     X <- exp(Y)
     
     data.frame(time = time_points, log_price = Y, price = X)
}

mcmc2_with_priors <- function(y, nsim, nburnin, prior_config = "noninformative") {
     nchain <- nburnin + nsim
     set.seed(seed)
     
     y.d <- rev(y)[-1] - rev(y)[-length(y)]
     n <- length(y.d)
     t.d <- rep(1, n) / 252
     
     # Prior configurations
     if (prior_config == "noninformative") {
          prior <- list(sig0 = 100, a0 = 0.01, b0 = 0.01, a_lam = 1, b_lam = 1)
     } else if (prior_config == "moderate") {
          prior <- list(sig0 = 10, a0 = 2, b0 = 2, a_lam = 2, b_lam = 2)
     } else if (prior_config == "informative") {
          prior <- list(sig0 = 1, a0 = 10, b0 = 10, a_lam = 5, b_lam = 1)
     }
     
     th.ls <- sig2.ls <- lamt.ls <- muz.ls <- sigz2.ls <- rep(NA, nchain)
     jt.ls <- zt.ls <- matrix(NA, nchain, n)
     th <- th.ls[1] <- muz <- muz.ls[1] <- 0
     sig2 <- sig2.ls[1] <- lamt <- lamt.ls[1] <- sigz2 <- sigz2.ls[1] <- 1e-3
     jt <- jt.ls[1,] <- zt <- zt.ls[1,] <- 0
     
     # Sampling steps
     for (k in 2:nchain) {
          # z_t
          sig2_zt <- 1 / (1 / sigz2 + jt^2 / sig2 / t.d)
          mu_zt <- sig2_zt * (muz / sigz2 + jt * (y.d - th * t.d) / sig2 / t.d)
          zt <- rnorm(n, mu_zt, sqrt(sig2_zt))
          zt.ls[k,] <- zt
          
          # j_t
          p1 <- lamt * dnorm(y.d, th * t.d + zt, sqrt(sig2 * t.d))
          p0 <- (1 - lamt) * dnorm(y.d, th * t.d, sqrt(sig2 * t.d))
          jt <- rbinom(n, 1, p1 / (p1 + p0))
          jt.ls[k,] <- jt
          
          # lambda
          lamt <- rbeta(1, prior$a_lam + sum(jt), prior$b_lam + sum(1 - jt))
          lamt.ls[k] <- lamt
          
          # theta
          sig2_th <- 1 / (1 / prior$sig0 + sum(t.d) / sig2)
          mu_th <- sig2_th * sum(y.d - jt * zt) / sig2
          th <- rnorm(1, mu_th, sqrt(sig2_th))
          th.ls[k] <- th
          
          # sigma^2
          a_sig2 <- prior$a0 + n / 2
          b_sig2 <- prior$b0 + sum((y.d - th * t.d - jt * zt)^2 / t.d) / 2
          sig2 <- rigamma(1, a_sig2, b_sig2)
          sig2.ls[k] <- sig2
          
          # mu_z
          sig2_muz <- 1 / (1 / prior$sig0 + n / sigz2)
          mu_muz <- sig2_muz * sum(zt) / sigz2
          muz <- rnorm(1, mu_muz, sqrt(sig2_muz))
          muz.ls[k] <- muz
          
          # sigma_z^2
          a_sigz2 <- prior$a0 + n / 2
          b_sigz2 <- prior$b0 + sum((zt - muz)^2) / 2
          sigz2 <- rigamma(1, a_sigz2, b_sigz2)
          sigz2.ls[k] <- sigz2
          
          # if (k %% 500 == 0) cat("Iteration:", k, "\n")
     }
     
     mu.ls <- th.ls + sig2.ls / 2
     sig.ls <- sqrt(sig2.ls)
     lam.ls <- lamt.ls / (1 / 252)
     
     return(list(th = th.ls, sig2 = sig2.ls, mu = mu.ls, sig = sig.ls,
                 lamt = lamt.ls, lam = lam.ls, jt = jt.ls, zt = zt.ls,
                 muz = muz.ls, sigz2 = sigz2.ls))
}

# ==== Step 3: Generate Predictions ====
generate_gbm_predictions <- function(th, sig2, y0, nsteps = 252) {
     nsim <- length(th)
     dt <- 1 / 252
     preds <- matrix(NA, nrow = nsim, ncol = nsteps)
     
     for (i in 1:nsim) {
          increments <- rnorm(nsteps, mean = th[i] * dt, sd = sqrt(sig2[i] * dt))
          preds[i, ] <- y0 + cumsum(increments)
     }
     
     return(preds)
}

evaluate_predictions <- function(predictions, true_values) {
     pred_mean <- colMeans(predictions)
     rmse <- sqrt(mean((pred_mean - true_values)^2))
     mae <- mean(abs(pred_mean - true_values))
     return(c(RMSE = rmse, MAE = mae))
}

