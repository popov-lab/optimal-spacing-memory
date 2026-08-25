library(ggplot2)
library(purrr)
library(dplyr)

optimal_power <- function(y, delta, d) {
  denum <- delta^(-1 / (d + 1)) * (1+y)^(d / (1 + d)) - 1
  y/denum - 1
}

construct_curves <- function(y, delta, d) {
  data.frame(
    delta = delta,
    d = d,
    y = y,
    x = optimal_power(y + 1, delta, d)
  )
}

SECONDS_PER_DAY <- 24 * 60 * 60
MAX_DELAY <- 1000 * SECONDS_PER_DAY
MIN_DELAY <- 1
test_lag <- 10^seq(log10(MIN_DELAY), log10(MAX_DELAY), length = 1000)

conditions <- expand.grid(delta = c(0.25, 0.5, 0.8, 1), d = c(0.01, 0.25, 0.5, 0.8))
curves_data <- pmap(conditions, \(delta, d) construct_curves(test_lag, delta, d)) |> dplyr::bind_rows()

breaks <- c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000)

curves_data |> 
  filter(y > 1/delta - 1)  |> 
  ggplot(aes(y / SECONDS_PER_DAY, x / SECONDS_PER_DAY, color = as.factor(delta))) +
  geom_line() +
  scale_x_log10(name = "Test delay (days)", breaks = breaks, labels = breaks) +
  scale_y_log10(name = "Study gap (days)", breaks = breaks, labels = breaks) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  ggtitle("Panels: decay rate (d)") +
  scale_color_discrete("Learning rate (delta)") +
  facet_wrap(~d, scales = "free") +
  theme(legend.position = "bottom")

ggsave("figures/f2_spacing_parameters_log_log.png", f2, width = 8, height = 8)



