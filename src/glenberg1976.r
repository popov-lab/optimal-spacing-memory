library(ggplot2)
library(dplyr)
library(here)
library(patchwork)
library(tidyr)

sac <- function(x, y, delta, d) {
  delta * ((1 + y)^(-d) + (1 + x + y)^(-d) - delta * (1 + x)^(-d) * (1 + y)^(-d))
}

glenberg_design <- expand.grid(lag = c(0, 1, 4, 8, 20, 40), ri = c(2, 8, 32, 64))

glenberg_design$sac <- sac(
  glenberg_design$lag * 3, glenberg_design$ri * 3,
  delta = 0.7, d = 0.1
)


ggplot(glenberg_design, aes(lag, sac, color = as.factor(ri))) +
  geom_point() +
  geom_line()
