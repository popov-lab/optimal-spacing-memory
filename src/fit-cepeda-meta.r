library(ggplot2)
library(dplyr)
library(here)

optimal_sac <- function(y, delta, d) {
  denum <- delta^(-1 / (d + 1)) * (1 + y)^(d / (1 + d)) - 1
  y / denum - 1
}

cepeda2006 <- read.csv(here("data/cepeda2006_optimal_gaps.csv"))

SECONDS_PER_DAY <- 24 * 60 * 60
MAX_DELAY <- 1000 * SECONDS_PER_DAY
MIN_DELAY <- 5
test_lag <- 10^seq(log10(MIN_DELAY), log10(MAX_DELAY), length = 1000)

sac_pred <- data.frame(
  test_delay_days = test_lag / SECONDS_PER_DAY,
  optimal_gap_days = optimal_sac(test_lag, 0.5, 0.1) / SECONDS_PER_DAY
)

obj <- function(par) {
  delta <- 1 / (1 + exp(-par["delta_"]))
  d <- 1 / (1 + exp(-par["d_"]))
  pred <- optimal_sac(cepeda2006$test_delay_days * SECONDS_PER_DAY, delta, d)
  sum((log(pred) - log(cepeda2006$optimal_gap_days * SECONDS_PER_DAY))^2)
}

optim(c(delta_ = 0, d_ = 0), obj)

ggplot(cepeda2006, aes(test_delay_days, optimal_gap_days)) +
  geom_point() +
  scale_x_log10(name = "Test delay (days)", breaks = breaks, labels = breaks) +
  scale_y_log10(name = "Study gap (days)", breaks = breaks, labels = breaks) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  geom_line(data = sac_pred, color = "red") +
  theme_bw()
