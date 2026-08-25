library(ggplot2)
library(dplyr)
library(here)
library(patchwork)
library(tidyr)

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

transform_parameters <- function(par) {
  list(
    delta = expit(par[["delta_"]]),
    d = expit(par[["d_"]]),
    mu = par[["mu"]],
    sigma = exp(par[["sigma_"]]),
    tau = exp(par[["tau_"]])
  )
}

obj <- function(par) {
  pars <- transform_parameters(par)
  s <- sac(
    cepeda2008$isi_days, cepeda2008$ri_days,
    pars$delta, pars$d, pars$tau
  )
  sqrt(mean((expit(s, pars$mu, pars$sigma) - cepeda2008$acc)^2))
}

cepeda <- read.csv(here("data/cepeda_spacing_recall.csv"))
cepeda$acc <- cepeda$recall_pct / 100
cepeda2008 <- filter(cepeda, experiment == "Cepeda et al. (2008)")

cepeda2008$isi_days[cepeda2008$isi_days == 0] <- 0.00256

inits <- c(delta_ = logit(0.40), d_ = logit(0.13), mu = 0.24, sigma_ = log(0.0269), tau_ = log(0.033))

N_RANDOM_STARTS <- 49L

random_inits <- replicate(
  N_RANDOM_STARTS,
  c(
    delta_ = logit(runif(1, 0.05, 0.95)),
    d_ = logit(runif(1, 0.05, 0.95)),
    mu = runif(1, -0.25, 1.25),
    sigma_ = runif(1, log(0.005), log(0.5)),
    tau_ = runif(1, log(0.001), log(100))
  ),
  simplify = FALSE
)

starts <- c(list(inits), random_inits)
fits <- lapply(starts, function(start) {
  tryCatch(
    optim(start, obj, control = list(maxit = 1e6)),
    error = function(error) NULL
  )
})

finite_fits <- vapply(fits, function(fit) {
  !is.null(fit) && is.finite(fit$value)
}, logical(1))

if (!any(finite_fits)) {
  stop("All optimization starts failed")
}

converged_fits <- finite_fits & vapply(fits, function(fit) {
  !is.null(fit) && fit$convergence == 0
}, logical(1))

candidate_fits <- if (any(converged_fits)) converged_fits else finite_fits
candidate_indices <- which(candidate_fits)
best_index <- candidate_indices[
  which.min(vapply(fits[candidate_indices], function(fit) fit$value, numeric(1)))
]

fit <- fits[[best_index]]
fit$n_starts <- length(starts)
fit$n_successful_starts <- sum(finite_fits)
fit$start <- starts[[best_index]]

pars <- transform_parameters(fit$par)

cepeda2008 |>
  mutate(sac = sac(
    isi_days,
    ri_days,
    delta = pars$delta,
    d = pars$d,
    tau = pars$tau
  )) |>
  mutate(sac_acc = expit(sac, pars$mu, pars$sigma)) |>
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
    delta = pars$delta,
    d = pars$d,
    tau = pars$tau
  )) |>
  mutate(sac_acc = expit(sac, pars$mu, pars$sigma)) |>
  ggplot(aes(isi_days, acc, color = as.factor(ri_days))) +
  geom_point() +
  geom_line(aes(y = sac_acc)) +
  theme(legend.position = "bottom")
