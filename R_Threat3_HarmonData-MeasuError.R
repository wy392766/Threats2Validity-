# #############################################################
#
#  Threat 3. Harmonization and measurement error
#
# #############################################################

rm(list = ls())
objects()
ls()

#############################################################
# Figure 3a: Bland–Altman or scatter-difference plot
#############################################################

set.seed(123)
library(tidyverse)

# simulate "true" site-specific measure (e.g., BMI from site lab)
n <- 300
true_val <- rnorm(n, mean = 25, sd = 4)

# site-specific measure (more precise)
site_measure <- true_val + rnorm(n, 0, 0.8)

# harmonized measure (coarser, maybe rounded or from questionnaire)
harmon_measure <- true_val + rnorm(n, 0, 1.5)

bland_dat <- tibble(
  mean_val = (site_measure + harmon_measure) / 2,
  diff_val = site_measure - harmon_measure
)

mean_diff <- mean(bland_dat$diff_val)
sd_diff   <- sd(bland_dat$diff_val)
loa_upper <- mean_diff + 1.96 * sd_diff
loa_lower <- mean_diff - 1.96 * sd_diff

p3a <- bland_dat %>%
  ggplot(aes(x = mean_val, y = diff_val)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = mean_diff, color = "blue", linetype = "solid") +
  geom_hline(yintercept = loa_upper, color = "red", linetype = "dashed") +
  geom_hline(yintercept = loa_lower, color = "red", linetype = "dashed") +
  labs(
    #title = "Bland–Altman plot: site-specific vs harmonized measure",
    x = "Mean of two measures",
    y = "Site-specific – Harmonized"
  ) +
  theme_minimal(base_size = 10)

print(p3a)


# ############################################################
# Figure 3b: Regression attenuation plot
# ############################################################

library(tidyverse)

set.seed(2026)
n <- 800

# Latent true measurement
true_value <- rnorm(n, mean = 25, sd = 4.5)

# Site-specific measurement: relatively precise
site_specific <- true_value + rnorm(n, mean = 0, sd = 0.8)

# Harmonized measurement: substantial harmonization-related error
harmonized <- true_value + rnorm(n, mean = 0, sd = 8)

# Outcome depends on the latent true value
beta_true <- 0.7
outcome <- 5 + beta_true * true_value + rnorm(n, mean = 0, sd = 3)

dat <- tibble(
  id = 1:n,
  true_value = true_value,
  site_specific = site_specific,
  harmonized = harmonized,
  outcome = outcome
)

# Regression models
mod_true <- lm(outcome ~ true_value, data = dat)
mod_site <- lm(outcome ~ site_specific, data = dat)
mod_harm <- lm(outcome ~ harmonized, data = dat)

coef_true <- coef(mod_true)[2]
coef_site <- coef(mod_site)[2]
coef_harm <- coef(mod_harm)[2]

r2_true <- summary(mod_true)$r.squared
r2_site <- summary(mod_site)$r.squared
r2_harm <- summary(mod_harm)$r.squared

plot_dat <- dat %>%
  select(outcome, true_value, site_specific, harmonized) %>%
  pivot_longer(
    cols = c(true_value, site_specific, harmonized),
    names_to = "measure_type",
    values_to = "x"
  ) %>%
  mutate(
    measure_type = factor(
      measure_type,
      levels = c("true_value", "site_specific", "harmonized"),
      labels = c("True value", "Site-specific", "Harmonized")
    )
  )

p3b <- ggplot(plot_dat, aes(x = x, y = outcome)) +
  geom_point(
    shape = 16,
    size = 1.1,
    alpha = 0.18,
    color = "gray35"
  ) +
  geom_smooth(
    data = subset(plot_dat, measure_type == "True value"),
    method = "lm",
    se = FALSE,
    linewidth = 0.9,
    color = "black"
  ) +
  geom_smooth(
    data = subset(plot_dat, measure_type == "Site-specific"),
    method = "lm",
    se = FALSE,
    linewidth = 0.8,
    color = "blue"
  ) +
  geom_smooth(
    data = subset(plot_dat, measure_type == "Harmonized"),
    method = "lm",
    se = FALSE,
    linewidth = 0.8,
    linetype = "dashed",
    color = "red3"
  ) +
  annotate(
    "text",
    x = min(plot_dat$x) +15,
    y = max(dat$outcome) - 0.5,
    label = paste0("True slope = ", round(coef_true, 2),
                   ", R² = ", round(r2_true, 2)),
    hjust = 0,
    size = 3.2,
    color = "black"
  ) +
  annotate(
    "text",
    x = min(plot_dat$x) +15,
    y = max(dat$outcome) - 1.4,
    label = paste0("Site-specific slope = ", round(coef_site, 2),
                   ", R² = ", round(r2_site, 2)),
    hjust = 0,
    size = 3.2,
    color = "blue"
  ) +
  annotate(
    "text",
    x = min(plot_dat$x) +15,
    y = max(dat$outcome) - 2.3,
    label = paste0("Harmonized slope = ", round(coef_harm, 2),
                   ", R² = ", round(r2_harm, 2)),
    hjust = 0,
    size = 3.2,
    color = "red3"
  ) +
  labs(
    x = "Measurement value",
    y = "Outcome"
  ) +
  coord_cartesian(xlim = c(10, 41)) +
    theme_minimal(base_size = 11) +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

print(p3b)


############################################################
# Figure 3c: Attenuation plot from simulation
############################################################

set.seed(202)
library(tidyverse)

# true data-generating: Y = 2*X_true + error
n <- 2000
X_true <- rnorm(n, 0, 1)
eps <- rnorm(n, 0, 1)
Y <- 2 * X_true + eps  # true slope = 2

# try increasing levels of measurement error in X
err_sd_vals <- seq(0, 2, by = 0.2)

att_res <- map_dfr(err_sd_vals, function(esd) {
  X_obs <- X_true + rnorm(n, 0, esd)
  fit <- lm(Y ~ X_obs)
  tibble(
    meas_error_sd = esd,
    est_beta = coef(fit)[2]
  )
})

p3c <- att_res %>%
  ggplot(aes(x = meas_error_sd, y = est_beta)) +
  geom_line(color = "steelblue", linewidth = 0.75) +
  geom_point(color = "steelblue", size = 2) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  scale_x_continuous(limits = c(0, 2.0), breaks = seq(0, 2.0, by = 0.2)) +
  scale_y_continuous(limits = c(0.2, 2.0), breaks = seq(0.2, 2.0, by = 0.2)) +
  labs(
    #title = "Attenuation of effect with increasing measurement error",
    x = "SD of measurement error added to exposure",
    y = "Estimated slope (Y ~ X_observed)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  )

print(p3c)

