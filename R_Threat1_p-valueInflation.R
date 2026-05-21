# #############################################################
#
#  Threat 1. p-value inflation and trivial effects
#
# #############################################################

rm(list = ls())
objects()
ls()

############################################################
# Figure 1a: P-value vs sample size (fixed small effect)
############################################################

library(tidyverse)

set.seed(12345)

# Simulation parameters
true_d <- 0.15          # small standardized mean difference
sigma  <- 1
ns <- seq(100, 9000, by = 200)
B <- 300               # repeated simulations per sample size

# Function to simulate one two-sample t-test
get_p_for_n <- function(n, d, sigma = 1) {
  n1 <- n2 <- n / 2
  x1 <- rnorm(n1, mean = 0, sd = sigma)
  x2 <- rnorm(n2, mean = d * sigma, sd = sigma)
  t.test(x1, x2, var.equal = TRUE)$p.value
}

# Repeated simulations
df_sim <- expand_grid(n = ns, rep = 1:B) %>%
  mutate(
    p_value = map_dbl(n, get_p_for_n, d = true_d, sigma = sigma)
  )

# Median p-value across simulations
df_summary <- df_sim %>%
  group_by(n) %>%
  summarise(
    p_median = median(p_value),
    p_q25 = quantile(p_value, 0.25),
    p_q75 = quantile(p_value, 0.75),
    .groups = "drop"
  )

# Theoretical expected p-value using expected t-statistic
df_theory <- tibble(n = ns) %>%
  mutate(
    n1 = n / 2,
    n2 = n / 2,
    expected_t = true_d / sqrt(1 / n1 + 1 / n2),
    df = n - 2,
    p_theory = 2 * pt(-abs(expected_t), df = df)
  )

# Plot
fig1a <- ggplot() +
  geom_point(
    data = df_sim,
    aes(x = n, y = p_value),
    size = 0.5,
    alpha = 0.12,
    color = "steelblue"
  ) +
  geom_line(
    data = df_summary,
    aes(x = n, y = p_median),
    color = "darkblue",
    linewidth = 1.1
  ) +
  geom_line(
    data = df_theory,
    aes(x = n, y = p_theory),
    color = "black",
    linetype = "dashed",
    linewidth = 0.9
  ) +
  geom_hline(
    yintercept = 0.05,
    linetype = "dashed",
    color = "red",
    linewidth = 0.7
  ) +
  scale_y_log10(
    limits = c(1e-10, 1),
    breaks = 10^(-10:0)
  ) +
  labs(
    x = "Total sample size",
    y = "P-value (log10 scale)",
    caption = "Light points show repeated simulations; blue line shows median simulated p-value; black dashed line shows theoretical expected p-value."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    plot.caption = element_text(size = 8, hjust = 0)
  )

print(fig1a)


############################################################
# Figure 1b: Effect size vs -log10(p) scatter
############################################################

library(tidyverse)

set.seed(456)
n_obs <- 5000
n_tests <- 200

sim_tests <- map_dfr(1:n_tests, function(j) {
  
  beta_true <- rnorm(1, mean = 0, sd = 0.10)
  x <- rnorm(n_obs, mean = 0, sd = 1)
  error <- rnorm(n_obs, mean = 0, sd = 1)
  y <- beta_true * x + error
  
  fit <- lm(y ~ x)
  sm <- summary(fit)$coef
  
  tibble(
    test_id = j,
    beta_true = beta_true,
    beta_est = sm["x", "Estimate"],
    se = sm["x", "Std. Error"],
    p_value = sm["x", "Pr(>|t|)"],
    neglog10p = -log10(p_value)
  )
})

fig1b <- sim_tests %>%
  ggplot(aes(x = beta_est, y = neglog10p)) +
  geom_point(alpha = 0.70, size = 1.4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  labs(
    x = "Estimated effect size (beta)",
    y = expression(-log[10](p)),
    caption = "Each point represents one simulated regression test. True beta values were drawn from N(0, 0.10²)."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    plot.caption = element_text(size = 8, hjust = 0)
  )

print(fig1b)


############################################################
# Figure 1c: Two-panel comparison small-N vs large-N
############################################################

library(tidyverse)
set.seed(789)

simulate_reg <- function(n, p = 10) {
  
  # Unmeasured confounder
  U <- rnorm(n)
  # Predictors
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  x4 <- rnorm(n)
  
  # x5 and x6 are correlated with unmeasured U,
  # creating spurious/confounded associations
  x5 <- 0.40 * U + rnorm(n, 0, sqrt(1 - 0.40^2))
  x6 <- 0.30 * U + rnorm(n, 0, sqrt(1 - 0.30^2))
  
  # Pure null predictors
  x7  <- rnorm(n)
  x8  <- rnorm(n)
  x9  <- rnorm(n)
  x10 <- rnorm(n)
  
  dat <- tibble(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
  
  # Outcome model:
  # x1 and x2 have meaningful true effects
  # x3 and x4 have tiny true effects
  # U affects outcome but is omitted from the fitted model
  y <- 0.30 * x1 +
    0.20 * x2 +
    0.02 * x3 -
    0.015 * x4 +
    0.25 * U +
    rnorm(n, 0, 1)
  
  dat$y <- y
  
  fit <- lm(y ~ ., data = dat)
  sm <- summary(fit)$coef
  
  tibble(
    predictor = rownames(sm)[-1],
    estimate = sm[-1, "Estimate"],
    p_value = sm[-1, "Pr(>|t|)"]
  )
}

predictor_labels <- tibble(
  predictor = paste0("x", 1:10),
  role = c(
    "True meaningful effect",
    "True meaningful effect",
    "Trivial true effect",
    "Trivial true effect",
    "Spurious/confounded",
    "Spurious/confounded",
    "Null",
    "Null",
    "Null",
    "Null"
  )
)

res_small <- simulate_reg(n = 1000) %>%
  mutate(sample = "n = 1,000")
res_large <- simulate_reg(n = 1000000) %>%
  mutate(sample = "n = 1,000,000")
res_both <- bind_rows(res_small, res_large) %>%
  left_join(predictor_labels, by = "predictor") %>%
  mutate(
    neglog10p = -log10(p_value),
    significant = if_else(p_value < 0.05, "p < 0.05", "Not significant"),
    predictor_role = paste0(predictor, " (", role, ")")
  )

fig1c <- res_both %>%
  ggplot(
    aes(
      x = reorder(predictor_role, neglog10p),
      y = neglog10p,
      fill = role
    )
  ) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ sample, ncol = 2) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  labs(
    x = "Predictor",
    y = expression(-log[10](p)),
    fill = "Predictor type",
    caption = "Dashed horizontal line corresponds to p = 0.05. x5 and x6 are non-causal but associated with an omitted confounder."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    strip.text = element_text(size = 9, face = "bold"),
    plot.caption = element_text(size = 8, hjust = 0)
  )

print(fig1c)



