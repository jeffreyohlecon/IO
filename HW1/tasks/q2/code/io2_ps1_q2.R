library(dplyr)
library(AER)
library(tidyr)
library(xtable)

# Read in the dataset
df <- read.csv("../input/ps1_ex2.csv")
names(df) <- c("market","product","s","p","x","z")

#----------------------------------------------------------------------------#
# 1. Prepare data & estimate (alpha, beta) via IV
#    Model:  ln(s_j/s_0) = - alpha * p_j + beta * x_j
#            => y_j = - alpha * p_j + beta * x_j, with y_j = ln(s_j) - ln(s_0)
#            p_j instrumented by z_j
#----------------------------------------------------------------------------#
# Compute outside share per market: s0 = 1 - sum(s_j for j=1-6)
df <- df %>%
  group_by(market) %>%
  mutate(s0 = 1 - sum(s)) %>%
  ungroup()

# log-share ratio y = ln(s_j) - ln(s0)
df <- df %>%
  mutate(
    ln_s  = log(s),
    ln_s0 = log(s0),
    y     = ln_s - ln_s0
  )

# Run IV
iv_model <- ivreg(formula = y ~ 0 + p + x | 0 + z + x, data = df)
iv_summary <- summary(iv_model)

# extract coefficients
a_p <- coef(iv_model)["p"]   # this is negative of alpha
a_x <- coef(iv_model)["x"]   # this is beta
alpha_hat <- -a_p
beta_hat  <- a_x

# Create summary table
iv_table <- data.frame(
  Parameter = c("Alpha", "Beta"),
  Estimate = c(alpha_hat, beta_hat),
  StdError = c(iv_summary$coefficients["p", "Std. Error"],
               iv_summary$coefficients["x", "Std. Error"])
)
iv_latex_table <- xtable(iv_table, 
                         caption = "IV Regression Results: Estimating Alpha and Beta",
                         label = "tab:iv_results")
print(iv_latex_table, type = "latex", include.rownames = FALSE, file = "../output/iv_results.tex")

#----------------------------------------------------------------------------#
# 2. Compute own- and cross-price elasticities
#    e_{j,k} = (∂ s_j / ∂ p_k) * (p_k / s_j)
#      = - alpha * p_j (1 - s_j)   if j=k
#      = + alpha * p_k * s_k       if j != k
#----------------------------------------------------------------------------#
J <- 6
market_list <- unique(df$market)
elasticities_list <- list()

for(m in market_list) {
  
  df_m <- df %>% filter(market == m) %>% arrange(product)
  p_vec <- df_m$p
  s_vec <- df_m$s
  
  E_m <- matrix(0, nrow=J, ncol=J)
  for(j in 1:J){
    for(k in 1:J){
      if(j == k) {
        E_m[j,k] <- - alpha_hat * p_vec[j] * (1 - s_vec[j])
      } else {
        E_m[j,k] <- alpha_hat * p_vec[k] * s_vec[k]
      }
    }
  }
  elasticities_list[[as.character(m)]] <- E_m
}

# Average across markets
all_mats <- array(unlist(elasticities_list), dim = c(J,J,length(elasticities_list)))
E_avg <- apply(all_mats, c(1,2), mean)

cat("\n-- Average Own- and Cross-Price Elasticities (6x6) --\n")
print(round(E_avg, 4))

# Make table
dimnames(E_avg) <- list(
  paste0("Product",1:6),
  paste0("Product",1:6)
)

E_table <- xtable(E_avg, 
                  caption="Average Own- and Cross-Price Elasticities (6x6)",
                  label="tab:elasticities")
print(E_table, type="latex", include.rownames=TRUE, file="../output/elasticities.tex")

#----------------------------------------------------------------------------#
# 3. Recover Marginal Costs: Single-Product Nash-Bertrand
#    p_j - c_j = 1 / [ alpha (1 - s_j) ] => c_j = p_j - 1/(alpha (1 - s_j))
#----------------------------------------------------------------------------#
df <- df %>%
  mutate(mc = p - 1/(alpha_hat * (1 - s)))

avg_mc <- df %>%
  group_by(product) %>%
  summarise(
    avg_mc    = mean(mc),
    avg_price = mean(p),
    avg_share = mean(s)
  ) %>%
  ungroup()

cat("\n-- Average marginal costs, prices, shares by product --\n")
print(avg_mc)

mc_table <- xtable(avg_mc, 
                   caption="Average Marginal Costs, Prices, and Shares by Product",
                   label="tab:mc")
print(mc_table, type="latex", include.rownames=FALSE, file = "../output/prices_and_shares.tex")

#----------------------------------------------------------------------------#
# 4. Counterfactual: Remove product 1, solve new equilibrium among 2-6
#    We'll define a function to do iterative best responses
#----------------------------------------------------------------------------#
solve_bertrand_singleprod <- function(c_vec, alpha, beta, xi_vec, x_vec,
                                      tol=1e-9, max_iter=2000) {
  # c_vec:  cost for each product (length 5)
  # alpha, beta: scalars
  # xi_vec: unobserved "quality" (length 5)
  # x_vec:  observed characteristic (length 5)
  # returns a numeric vector p_vec (equilibrium prices)
  
  p_current <- c_vec + 0.5
  
  for(iter in seq_len(max_iter)) {
    delta_current <- -alpha * p_current + beta * x_vec + xi_vec
    exp_delta     <- exp(delta_current)
    sum_exp       <- 1 + sum(exp_delta)
    s_current     <- exp_delta / sum_exp  # logit share
    
    p_next <- c_vec + 1/(alpha*(1 - s_current))
    if(max(abs(p_next - p_current)) < tol) {
      return(p_next)
    }
    p_current <- p_next
  }
  
  warning("Did not converge")
  return(p_current)
}

# For the CF, need xi_j from baseline:
#   ln(s_j / s0) = - alpha p_j + beta x_j + xi_j
# => xi_j = ln(s_j / s0) + alpha p_j - beta x_j
df_cf <- data.frame()

for(m in market_list) {
  # the full baseline outside share for market m
  s0_full <- unique(df$s0[df$market == m])
  
  # subset: only products 2-6
  df_m <- df %>% filter(market == m, product != 1) %>% arrange(product)
  
  df_m <- df_m %>%
    mutate(
      xi = log(s / s0_full) + alpha_hat * p - beta_hat * x
    )
  
  c_vec  <- df_m$mc
  x_vec  <- df_m$x
  xi_vec <- df_m$xi
  
  # Solve new eq. prices
  p_cf <- solve_bertrand_singleprod(c_vec, alpha_hat, beta_hat, xi_vec, x_vec)
  
  # Compute shares
  delta_cf    <- -alpha_hat * p_cf + beta_hat * x_vec + xi_vec
  exp_delta_cf<- exp(delta_cf)
  sum_exp_cf  <- 1 + sum(exp_delta_cf)
  s_cf        <- exp_delta_cf / sum_exp_cf
  s0_cf       <- 1 / sum_exp_cf
  
  df_m$price_cf <- p_cf
  df_m$share_cf <- s_cf
  df_m$market   <- m
  
  df_cf <- rbind(df_cf, df_m)
}

# Summarize
avg_cf <- df_cf %>%
  group_by(product) %>%
  summarise(
    price_cf_avg = mean(price_cf),
    share_cf_avg = mean(share_cf)
  ) %>%
  ungroup()

cat("\n-- Counterfactual average prices (product 2-6) --\n")
print(avg_cf)

cf_table <- xtable(avg_cf,
                   caption="Counterfactual Average Prices and Shares (products 2-6)",
                   label="tab:cf")
print(cf_table, type="latex", include.rownames=FALSE, file="../output/prices_and_shares_cf.tex")

#----------------------------------------------------------------------------#
# 5. Profit changes & Consumer Surplus changes
#----------------------------------------------------------------------------#
# (a) Baseline profits: pi_j = (p_j - mc_j)* s_j   (per-capita or per-consumer measure)
df <- df %>%
  mutate(profit = (p - mc)*s)

# (b) Counterfactual profits
df_cf <- df_cf %>%
  mutate(profit_cf = (price_cf - mc)*share_cf)

# Merge
df_base <- df %>% select(market, product, profit)
df_cf_sub <- df_cf %>% select(market, product, profit_cf)

df_compare <- full_join(df_base, df_cf_sub, by=c("market","product")) %>%
  mutate(
    profit_cf = if_else(is.na(profit_cf), 0, profit_cf),  # product 1 -> 0 in CF
    delta_profit = profit_cf - profit
  )

profit_summary <- df_compare %>%
  group_by(product) %>%
  summarise(
    avg_profit_base = mean(profit),
    avg_profit_cf   = mean(profit_cf),
    avg_delta_profit= mean(delta_profit)
  ) %>%
  ungroup()

cat("\n-- Profit changes by product (average across markets) --\n")
print(profit_summary)

profit_table <- xtable(profit_summary,
                       caption="Average Baseline vs. Counterfactual Profits by Product",
                       label="tab:profit")
print(profit_table, type="latex", include.rownames=FALSE, file="../output/profits.tex")

# (c) Consumer Surplus changes
#     CS = (1/alpha)*ln(1 + sum_j exp(delta_j))
# with delta_j^0 = log(s_j/s0_j) in baseline, and similarly in CF

cs_function <- function(alpha, delta_vec) {
  (1/alpha)*log(1 + sum(exp(delta_vec)))
}

cs_changes <- data.frame()

for(m in market_list) {
  
  # baseline: all 6 products
  df_m_base <- df %>% filter(market==m)
  s0_m_base <- 1 - sum(df_m_base$s)
  # delta0_j = ln(s_j / s0_m)
  delta0 <- log(df_m_base$s / s0_m_base)
  
  cs_base_m <- cs_function(alpha_hat, delta0)
  
  # CF: no product 1
  df_m_cf <- df_cf %>% filter(market==m)
  s0_m_cf <- 1 - sum(df_m_cf$share_cf)
  delta1  <- log(df_m_cf$share_cf / s0_m_cf)
  
  cs_cf_m <- cs_function(alpha_hat, delta1)
  
  cs_changes <- rbind(cs_changes,
                      data.frame(market   = m,
                                 cs_base  = cs_base_m,
                                 cs_cf    = cs_cf_m,
                                 delta_cs = cs_cf_m - cs_base_m))
}

avg_delta_cs <- mean(cs_changes$delta_cs)

cat(sprintf("\nAverage change in consumer surplus (per capita) = %.6f\n", avg_delta_cs))

# Table
cs_summary_table <- data.frame(
  "Mean CS Change" = round(avg_delta_cs,4)
)
cs_xt <- xtable(cs_summary_table,
                caption="Average Change in Consumer Surplus (per capita)",
                label="tab:cs")
print(cs_xt, type="latex", include.rownames=FALSE, file="../output/consumer_surplus.tex")
