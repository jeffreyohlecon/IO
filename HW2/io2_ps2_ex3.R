data <- read.csv("ps2_ex3.csv")

negLogLik <- function(params, data) {
  beta  <- params[1]
  phi   <- params[2]
  delta <- params[3]
  
  # Compute the lower and upper thresholds for each market.
  lower_bound <- pnorm(phi + delta * log(data$n) - data$x * beta)
  upper_bound <- pnorm(phi + delta * log(data$n + 1) - data$x * beta)
  
  # The probability of observing n entrants is the difference.
  prob <- upper_bound - lower_bound
  
  # To avoid log(0), replace any very small probabilities with a small number.
  prob[prob < 1e-10] <- 1e-10
  
  # Compute the negative log-likelihood over all markets
  nll <- -sum(log(prob))
  return(nll)
}

init_params <- c(beta = 1, phi = 1, delta = 1)

optim_results <- optim(par = init_params, 
                       fn = negLogLik, 
                       data = data, 
                       method = "BFGS", 
                       control = list(maxit = 10000))

cat("Estimated Parameters:\n")
cat("beta =", optim_results$par[1], "\n")
cat("phi  =", optim_results$par[2], "\n")
cat("delta =", optim_results$par[3], "\n")