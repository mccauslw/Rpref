# The function create_alpha_prior precomputes information associated with the prior
# distribution of alpha, including lookup tables used for Metropolis-Hastings updates
# of the conditional posterior distribution of alpha given pi, described in Appendix B.
# Inputs are:
#  - n_fact, the value n!, where n is the number of elements in the choice universe.
#  - a, the shape parameter of the Gamma prior distribution of alpha.
#  - b, the rate parameter of the Gamma prior distribution of alpha.
#  - h, the distance between grid points in a grid of values of p_bar
#  - eps, the probability defining two extreme prior quantiles of alpha.
# Output is a structure with the fields:
#  - a, the shape parameter of the Gamma prior distribution of alpha.
#  - b, the rate parameter of the Gamma prior distribution of alpha.
#  - alpha_mode, the mode of alpha|pi as a function of p_bar on a grid of p_bar values.
#  - p_grid, a grid of values of p_bar
#  - psi_diff, the value psi(1+alpha) - psi(1+alpha/n!) evaluated on the grid of
#              alpha_mode values.
#  - p_min, the minimum value of p_bar on the grid.
#  - h, the distance between grid points on the grid of values of p_bar.
#  - funcs, evaluations of prior cdf and pdf functions on a grid of values of alpha
create_alpha_prior <- function(n_fact, a, b, h = 0.05, eps = 1e-7) {
  # Prior quantiles eps and 1-eps of alpha
  alpha_min = qgamma(eps, a, b)
  alpha_max = qgamma(eps, a, b, lower.tail=F)
  p_min <- -grad(alpha_min, 0, a, b, n_fact)
  p_max <- -grad(alpha_max, 0, a, b, n_fact)
  p_grid = seq(p_min, p_max, by=h)
  n_grid = length(p_grid)
  alpha_mode = rep(NA, n_grid)
  psi_diff = rep(NA, n_grid)
  for (i in seq(n_grid)) {
    res <- uniroot(grad, c(0,2000), p_grid[i], a, b, n_fact)
    alpha_mode[i] = res$root
    psi_diff[i] = digamma(1 + alpha_mode[i]) - digamma(1+alpha_mode[i]/n_fact)
  }

  # Compute prior density of alpha on a grid
  n_grid = 200
  grid <- seq(0, qgamma(0.995, a, b), length.out=n_grid)
  funcs <- list(cdf = list(x=grid, func=pgamma(grid, a, b), nse=rep(0, 200)),
                pdf = list(x=grid, func=dgamma(grid, a, b), nse=rep(0, 200)))
  list(a=a, b=b, alpha_mode = alpha_mode, p_grid = p_grid,
       psi_diff = psi_diff, p_min = p_min, h = h, funcs = funcs)
}

# The function grad evaluates the 1st derivative of log f(alpha|pi) with respect
# to alpha.
# Inputs are:
#  - alpha, the current value of alpha
#  - p_bar, a sufficient statistic for pi, equal to the arithmetic mean of the log
#           preference probabiliites.
#  - a, the shape parameter of the Gamma prior distribution of alpha.
#  - b, the rate parameter of the Gamma prior distribution of alpha.
#  - n_fact, the value n!, where n is the number of elements in the choice universe.
# Output is the value of the 1st derivative
grad <- function(alpha, p_bar, a, b, n_fact) {
  (a + n_fact - 2)/alpha - (b - p_bar) +
    digamma(1+alpha) - digamma(1+alpha/n_fact)
}
