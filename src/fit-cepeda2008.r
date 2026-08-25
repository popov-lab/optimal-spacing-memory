library(ggplot2)
library(dplyr)
library(here)
library(patchwork)
library(tidyr)

cepeda <- read.csv(here("data/cepeda_spacing_recall.csv"))
cepeda$acc <- cepeda$recall_pct / 100
cepeda2008 <- filter(cepeda, experiment == "Cepeda et al. (2008)")

cepeda2008$isi_days[cepeda2008$isi_days == 0] <- 0.00256

sac <- function(x, y, delta, d, tau) {
  f <- function(t) (1 + t / tau)^(-d)
  delta * (f(y) + f(x + y) - delta * f(x) * f(y))
}

expit <- function(s, mu = 0, sigma = 1) {
  1 / (1 + exp((mu - s) / sigma))
}

logit <- function(p) {
  log(p) - log(1 - p)
}

obj <- function(par) {
  delta <- expit(par["delta_"])
  d <- expit(par["d_"])
  sigma <- exp(par["sigma_"])
  mu <- par["mu"]
  tau <- exp(par["tau_"])
  s <- sac(
    cepeda2008$isi_days, cepeda2008$ri_days,
    delta, d, tau
  )
  sum((expit(s, mu, sigma) - cepeda2008$acc)^2)
}

fit <- optim(c(delta_ = logit(0.4013), d_ = logit(0.1367), mu = 0.2465, sigma_ = log(0.0267), tau_ = log(0.0318)), obj,
  control = list(maxit = 1e6)
)
fit
pars <- list(
  delta = expit(fit$par["delta_"]),
  d = expit(fit$par["d_"]),
  mu = fit$par["mu"],
  sigma = exp(fit$par["sigma_"]),
  tau = exp(fit$par["tau_"])
)
cepeda2008 |>
  mutate(sac = sac(
    isi_days,
    ri_days,
    delta = pars$delta,
    d = pars$d,
    tau = pars$tau
  )) |>
  mutate(sac_acc = logistic(sac, pars$mu, pars$sigma)) |>
  pivot_longer(c(acc, sac_acc)) |>
  ggplot(aes(isi_days, value, color = as.factor(ri_days))) +
  geom_point() +
  geom_line() +
  facet_wrap(~name, scales = "free") +
  theme(legend.position = "bottom")

cepeda2008 |>
  mutate(sac = sac(
    isi_days,
    ri_days,
    delta = logistic(fit$par["delta_"]),
    d = logistic(fit$par["d_"]),
    tau = exp(fit$par["tau_"])
  )) |>
  mutate(sac_acc = logistic(sac, fit$par["mu"], exp(fit$par["sigma_"]))) |>
  ggplot(aes(isi_days, acc, color = as.factor(ri_days))) +
  geom_point() +
  geom_line(aes(y = sac_acc)) +
  theme(legend.position = "bottom")
