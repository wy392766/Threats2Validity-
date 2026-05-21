# #############################################################
#
#  Threat 4. Causal identification in the presence of rich covariates
#
# #############################################################

rm(list = ls())
objects()
ls() 

############################################################
# Figure 4a: Bias-vs-adjustment plot
############################################################

library(tidyverse)

res <- tibble(
  model = factor(
    c(
      "Model 1: Crude\nY ~ A",
      "Model 2: + C1\nY ~ A + C1",
      "Model 3: + C1 + C2\nY ~ A + C1 + C2",
      "Model 4: Large set\nY ~ A + C1 + C2 + C3 + C4 + C5"
    ),
    levels = c(
      "Model 1: Crude\nY ~ A",
      "Model 2: + C1\nY ~ A + C1",
      "Model 3: + C1 + C2\nY ~ A + C1 + C2",
      "Model 4: Large set\nY ~ A + C1 + C2 + C3 + C4 + C5"
    )
  ),
  est = c(2.8, 2.45, 2.35, 2.35),
  lcl = c(2.5, 2.05, 1.65, 1.2),
  ucl = c(3.1, 2.8, 3.1, 3.4)
)

true_effect <- 2.05

fig4a <- ggplot(res, aes(x = model, y = est)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.10, linewidth = 0.6) +
  geom_line(aes(group = 1), linewidth = 0.5, alpha = 0.6) +
  geom_hline(
    yintercept = true_effect,
    linetype = "dashed",
    color = "red",
    linewidth = 0.8
  ) +
  scale_y_continuous(
    limits = c(0.5, 3.7),
    breaks = seq(0.5, 3.5, by = 0.5)
  ) +
  labs(
    x = "Adjustment set",
    y = "Estimated effect of A on Y"
  ) +
  coord_cartesian(ylim = c(1.1, 3.65)) +
  theme_minimal(base_size = 11) +
  theme(
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8, angle = 20, hjust = 0.7),
    axis.text.y = element_text(size = 8),
    panel.grid.minor = element_blank()
  )

print(fig4a)



# #############################################################
#  Figure 4B. Exposure misclassification produces persistent 
#             bias despite increasing sample size
# #############################################################

library(tidyverse)
set.seed(2025)

# Simulation parameters
n_sims <- 300
sample_sizes <- c(
  seq(200, 2000, by = 200),
  seq(3000, 10000, by = 1000),
  25000, 50000, 100000
)

# True parameters
true_beta <- 0.80

# Stronger exposure misclassification to illustrate bias clearly
sens <- 0.60
spec <- 0.70

simulate_one <- function(n, true_beta, sens, spec) {
  
  # True binary exposure
  X_true <- rbinom(n, 1, 0.5)
  
  # Outcome depends on true exposure
  Y <- 1 + true_beta * X_true + rnorm(n, 0, 1)
  
  # Observed exposure is misclassified
  X_obs <- ifelse(
    X_true == 1,
    rbinom(n, 1, sens),
    rbinom(n, 1, 1 - spec)
  )
  
  fit_true <- lm(Y ~ X_true)
  fit_obs  <- lm(Y ~ X_obs)
  
  tibble(
    n = n,
    beta_true_exposure = coef(fit_true)[2],
    beta_misclassified = coef(fit_obs)[2]
  )
}

results <- map_dfr(sample_sizes, function(n) {
  map_dfr(1:n_sims, ~ simulate_one(n, true_beta, sens, spec))
})

summary_df <- results %>%
  group_by(n) %>%
  summarise(
    mean_true = mean(beta_true_exposure),
    low_true = quantile(beta_true_exposure, 0.025),
    high_true = quantile(beta_true_exposure, 0.975),
    
    mean_obs = mean(beta_misclassified),
    low_obs = quantile(beta_misclassified, 0.025),
    high_obs = quantile(beta_misclassified, 0.975),
    .groups = "drop"
  )

fig4b <- ggplot(summary_df, aes(x = n)) +
  geom_ribbon(
    aes(ymin = low_true, ymax = high_true, fill = "True exposure"),
    alpha = 0.18
  ) +
  geom_line(
    aes(y = mean_true, color = "True exposure"),
    linewidth = 1.1
  ) +
  geom_ribbon(
    aes(ymin = low_obs, ymax = high_obs, fill = "Misclassified exposure"),
    alpha = 0.18
  ) +
  geom_line(
    aes(y = mean_obs, color = "Misclassified exposure"),
    linewidth = 1.1,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = true_beta,
    linetype = "dotted",
    color = "black",
    linewidth = 0.8
  ) +
  scale_x_continuous(
    trans = "log10",
    breaks = c(200, 500, 1000, 5000, 10000, 50000, 100000),
    labels = scales::comma
  ) +
  labs(
    x = "Sample size (log10 scale)",
    y = expression("Estimated effect " * (hat(beta))),
    color = "Model",
    fill = "Model"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(fig4b)

















