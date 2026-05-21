# #############################################################
#
#  Threat 2. Residual and differential confounding across sites
#
# #############################################################

rm(list = ls())
objects()
ls()

############################################################
# Figure 2a: Site-specific effect heterogeneity plot
############################################################

set.seed(321)
library(tidyverse)

# Parameters
n_site <- 10
n_per_site <- 500
beta_true <- 0.40

sim_list <- vector("list", n_site)

for (s in 1:n_site) {
  
  # Unmeasured confounder distribution differs by site
  u_mean <- runif(1, -1.5, 1.5)
  u <- rnorm(n_per_site, mean = u_mean, sd = 1)
  
  # Site-specific confounding strength varies
  gamma_x <- runif(1, -1.0, 1.0)   # U -> X
  gamma_y <- runif(1, -1.2, 1.2)   # U -> Y
  
  # Exposure
  x <- gamma_x * u + rnorm(n_per_site, 0, 1)
  
  # Outcome
  y <- beta_true * x + gamma_y * u + rnorm(n_per_site, 0, 1)
  
  sim_list[[s]] <- tibble(
    site = paste0("Site ", s),
    y = y,
    x = x
  )
}

dat_all <- bind_rows(sim_list)

# Site-specific naive regressions
site_results <- dat_all %>%
  group_by(site) %>%
  group_modify(~{
    fit <- lm(y ~ x, data = .x)
    sm <- summary(fit)$coef
    
    tibble(
      beta_hat = sm["x","Estimate"],
      se = sm["x","Std. Error"],
      lcl = beta_hat - 1.96*se,
      ucl = beta_hat + 1.96*se
    )
  }) %>%
  ungroup()

# Plot
p2a <- ggplot(site_results,
                        aes(x = reorder(site, beta_hat), y = beta_hat)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.12) +
  geom_hline(yintercept = beta_true,
             linetype = "dashed",
             color = "red",
             linewidth = 0.8) +
  geom_hline(yintercept = 0,
             linetype = "dotted",
             color = "grey40") +
  coord_flip() +
  labs(
    x = "Site",
    y = "Estimated Exposure Effect",
    title = "Site-specific effect estimates show marked heterogeneity despite a common true causal effect"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9)
  )

print(p2a)


############################################################
# Figure 2b: Simulated bias amplification plot
############################################################

set.seed(456)
library(tidyverse)

# function to simulate data for a given level of imbalance
sim_bias <- function(imbalance = 0, n_site = 6, n_per_site = 400,
                     beta_true = 0.4) {
  sim_list <- vector("list", n_site)
  for (s in 1:n_site) {
    # site-specific mean of unmeasured confounder increases with imbalance
    u_mean <- (s - 1) * imbalance  # bigger imbalance => bigger spread across sites
    u <- rnorm(n_per_site, mean = u_mean, sd = 0.5)
    
    # exposure depends on u
    x <- 0.6 * u + rnorm(n_per_site, 0, 1)
    # outcome depends on exposure and u
    y <- beta_true * x + 0.8 * u + rnorm(n_per_site, 0, 1)
    
    sim_list[[s]] <- tibble(site = factor(s),
                            y = y,
                            x = x,
                            u = u)
  }
  dat <- bind_rows(sim_list)
  
  # naive pooled model, not adjusting for u
  fit_naive <- lm(y ~ x, data = dat)
  beta_naive <- coef(summary(fit_naive))["x", "Estimate"]
  
  tibble(imbalance = imbalance,
         beta_hat = beta_naive)
}

# try several imbalance levels
imbalance_vals <- seq(0, 0.6, by = 0.1)
bias_res <- map_dfr(imbalance_vals, sim_bias)

p2b <- bias_res %>%
  ggplot(aes(x = imbalance, y = beta_hat)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 2) +
  geom_hline(yintercept = 0.4, linetype = "dashed", color = "red") +
  scale_x_continuous(limits = c(0, 0.61), breaks = seq(0, 0.7, by = 0.1)) +
  scale_y_continuous(limits = c(0.35, 0.9), breaks = seq(0, 1, by = 0.1)) +
  labs(
    #title = "Bias increases as cross-site imbalance in an unmeasured confounder grows",
    x = "Cross-site imbalance in unmeasured confounder",
    y = "Naive pooled estimate of exposure effect"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8)
  )

print(p2b)




