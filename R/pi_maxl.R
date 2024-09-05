# Functions associated with the random preference (RP) log likelihood
#
#' Evaluate random preference (RP) log likelihood
#'
#' Evaluate RP log likelihood as a function of pi,
#' the vector of preference probabilities, for given data N.
#'
#' @param pi vector of preference probabilities of length n!
#' @param u list with precomputed information about choice universes
#'      of size n, created using [create_universe()]
#' @param N matrix containing observed choice counts in a choice experiment
#'
#' @return The value of the RP log likelihood at the given value of pi
#' @export
#'
#' @examples
#' library(RanCh)
#' n <- 5
#' n_fact = factorial(n)
#' pi <- rep(1/n_fact, n_fact) # Uniform distribution over preference orders
#' u <- create_universe(n)
#' N <- vectorize_counts(u, RanCh::MMS_2019_counts[1, , ])
#' pi_ln_like <- compute_pi_ln_like(pi, u, N)
#'
#' @inherit create_universe author references
#'
compute_pi_ln_like <- function(pi, u, N) {
  P_Ax <- u$pi_to_P %*% pi
  P_Ax[N$by_Ax == 0] <- 1
  ln_like <- sum(N$by_Ax * log(P_Ax))
}

#' Evaluate gradient of random preference (RP) log likelihood
#'
#' Evaluate gradient, with repect to pi, of RP log likelihood,
#' as a function of pi, for given data N.
#'
#' @inheritParams compute_pi_ln_like
#'
#' @return Gradient of RP log likelihood, a vector of size n!.
#' @export
#'
#' @examples
#' library(RanCh)
#' n <- 5
#' n_fact = factorial(n)
#' pi <- rep(1/n_fact, n_fact) # Uniform distribution over preference orders
#' u <- create_universe(n)
#' N <- vectorize_counts(u, RanCh::MMS_2019_counts[1, , ])
#' pi_score <- compute_pi_score(pi, u, N)
#'
#' @inherit create_universe author references
#'
compute_pi_score <- function(pi, u, N) {
  P_Ax <- as.vector(u$pi_to_P %*% pi)
  P_Ax[N$by_Ax == 0] <- 1
  score <- colSums((N$by_Ax/P_Ax) * u$pi_to_P)
}

#' Maximize the RP log likelihood for given data N
#'
#' Using an algorithm for convex/concave programming on probability simplices,
#' iteratively find the maximum value of the RP log likelihood.
#' The basic algorithm is described in Chok and Vasil (2023) and a modification
#' to the algorithm to incorporate the momentum idea from the Machine Learning
#' literature is described in McCausland (unpub.)
#' While the RP log likelihood is concave in pi, for choice universes with more
#' than 3 elements, it is not strictly concave, pi is not identified and the
#' returned value of pi will depend on the starting value.
#'
#' @inheritParams compute_pi_ln_like
#' @param init_pi starting value for vector pi that defaults to the discrete uniform distriubtion.
#' @param delta parameter of a halt condition. If the range of log likelihood gradient values
#          associated with positive elements of pi is less than delta, and the gradient
#          values associated with elements of pi equal to zero are all less than those
#          associated with positive elements, then the iterations halt.
#' @param epsilon lowest value of pi that is considered distinct from zero.
#' @param n_max_iter maximum number of iterations
#' @param beta a momemtum parameter
#' @param debug a flag that should normally be set to FALSE. If TRUE, some diagnostic
#' information is displayed to the screen.
#' @param doplot a flag that should normally be set to FALSE. If TRUE, a diagnostic
#' plots is displayed.
#'
#' @return A list with the elements
#' \describe{
#'   \item{pi}{final value of the probability vector pi}
#'   \item{ln_maxl}{final value of the RP log likelihood.}
#'   \item{score}{final value of the RP log likelihood gradient.}
#'   \item{z}{final search direction, a decaying linear function of previous score values.}
#'   \item{S}{logical vector of length n! indicating which elements of pi are greater than one.}
#'   \item{n_iter}{number of realized iterations.}
#'   \item{score_range}{difference between maximum and minimum gradient values
#'   associated with positive elements of pi.}
#'   \item{ln_maxl_diff}{difference of RP log likelihood values between the
#'   2nd last and last iterations.}
#' }
#' @export
#'
#' @examples
#' library(RanCh)
#' n <- 5
#' u <- create_universe(n)
#' N <- vectorize_counts(u, RanCh::MMS_2019_counts[1, , ])
#' pi_ln_maxl <- compute_pi_ln_maxl(u, N)
#'
#' @inherit create_universe author references
#' @references
#' Chok, J. and G. M. Vasil (2023). Convex Optimization over a Probability Simplex, arXiv:2305.09046
#'
compute_pi_ln_maxl <- function(u, N, init_pi = rep(1/u$n_orders, u$n_orders),
                               delta = 0.0001, epsilon = 1e-8, n_max_iter = 15000,
                               beta = 0.95,
                               debug = FALSE, doplot = FALSE) {
  eta <- rep(NA, 4)
  iota <- rep(1, u$n_orders)
  v <- rep(NA, 4)
  prev_val <- -Inf
  n_iter <- 0
  curr_pi <- init_pi
  curr_val <- compute_pi_ln_like(curr_pi, u, N)
  curr_score <- compute_pi_score(curr_pi, u, N)
  curr_z <- curr_score
  S <- rep(TRUE, u$n_orders); Q <- !S
  while ((max(curr_score[S]) - min(curr_score[S]) > delta) &&
         n_iter < n_max_iter) {

    if (doplot) {
      plot(curr_pi, curr_score, xlim=c(0,0.08), ylim=c(1000, 1100))
      readline("Next plot ")
    }
    i_min <- which.min(curr_z[S])
    eta_max <- 1 / (sum(curr_pi * curr_z) - min(curr_z[S]))
    d <- curr_pi * (curr_z - sum(curr_pi * curr_z))

    phi_inv = 0.5 * (sqrt(5)-1)
    eta[1] <- 0
    eta[2] <- eta_max
    eta[3] <- (1 - phi_inv) * eta_max
    eta[4] <- phi_inv * eta_max
    v[1] <- curr_val
    for (i in 3:4) v[i] <- compute_pi_ln_like(curr_pi + eta[i] * d, u, N)
    pi_eta_max <- curr_pi + eta_max * d
    pi_eta_max[pi_eta_max < 0] = 0
    v_eta_max <- compute_pi_ln_like(pi_eta_max, u, N)
    v[2] <- v_eta_max

    # Line search
    i <- 1
    while (i <= 10 || v[1] > v[3]) {
      if (v[3] < v[4]) {
        eta[1] <- eta[3]; v[1] <- v[3]
        eta[3] <- eta[4]; v[3] <- v[4]
        eta[4] <- (1 - phi_inv) * eta[1] + phi_inv * eta[2]
        v[4] <- compute_pi_ln_like(curr_pi + eta[4] * d, u, N)
      }
      else {
        eta[2] <- eta[4]; v[2] <- v[4]
        eta[4] <- eta[3]; v[4] <- v[3]
        eta[3] <- (1 - phi_inv) * eta[2] + phi_inv * eta[1]
        v[3] <- compute_pi_ln_like(curr_pi + eta[3] * d, u, N)
      }
      i <- i+1
    }

    prev_val <- curr_val
    j <- which.max(v)
    curr_pi <- curr_pi + eta[j] * d
    if ((j == 2) && v[2] == v_eta_max) curr_pi <- pi_eta_max

    S <- (curr_pi >= epsilon)
    Q <- !S

    if (debug)
      print(c(j, eta[j], eta_max, v[j] - curr_val, v_eta_max - curr_val, sum(S)))

    curr_pi[Q] <- 0
    curr_pi[S] = curr_pi[S]/sum(curr_pi[S])
    curr_val <- compute_pi_ln_like(curr_pi, u, N)
    curr_score <- compute_pi_score(curr_pi, u, N)
    curr_z <- beta * curr_z + curr_score
    n_iter <- n_iter + 1
  }
  list(pi = curr_pi, ln_maxl = curr_val,
       score = curr_score, z = curr_z, S = S, n_iter = n_iter,
       score_range = max(curr_score[S]) - min(curr_score[S]),
       ln_maxl_diff = curr_val - prev_val)
}

# This is obsolete code using MCMC simulation to get a value of pi with high
# posterior probability, for use as a starting value.
# WJM: consider deleting.
sim_MCMC <- function(u, M, N, alpha, gamma_p = NULL) {
  if (is.null(gamma_p))
    gamma_p = stats::rgamma(u$n_orders, alpha/u$n_orders)
  big_bl_n = factorial(u$n-1)
  small_bl_n = factorial(u$n-2)

  sum_gamma_p = rep.int(NA, M)
  sum_ln_gamma_p = rep.int(NA, M)
  ln_like = rep.int(NA, M)

  sum_gamma_p = sum(gamma_p)
  pi = gamma_p / sum_gamma_p
  ln_like_den = compute_pi_ln_like(pi, u, N)

  sum_gamma_p[1] = sum_gamma_p; print(sum_gamma_p[1])
  sum_ln_gamma_p[1] = sum(log(gamma_p)); print(sum_ln_gamma_p[1])
  ln_like[1] = ln_like_den; print(ln_like[1])

  for (m in seq(2:M)) {
    for (i in seq(1)) {
      for (big_bl in seq(0, u$n-1)) {
        gamma_p_star = gamma_p
        gamma_p_star[(big_bl * big_bl_n + 1):
                       ((big_bl + 1) * big_bl_n)] =
          stats::rgamma(big_bl_n, alpha/u$n_orders)
        sum_gamma_p_star = sum(gamma_p_star)
        pi_star = gamma_p_star / sum_gamma_p_star
        ln_like_num = compute_pi_ln_like(pi_star, u, N)
        if (stats::runif(1) < exp(ln_like_num - ln_like_den)) {
          pi = pi_star
          ln_like_den = ln_like_num
          gamma_p = gamma_p_star
        }
      }
      for (small_bl in seq(0, u$n * (u$n-1) - 1)) {
        gamma_p_star = gamma_p
        gamma_p_star[(small_bl * small_bl_n + 1):
                       ((small_bl + 1) * small_bl_n)] =
          stats::rgamma(small_bl_n, alpha/u$n_orders)
        sum_gamma_p_star = sum(gamma_p_star)
        pi_star = gamma_p_star / sum_gamma_p_star
        ln_like_num = compute_pi_ln_like(pi_star, u, N)
        if (stats::runif(1) < exp(ln_like_num - ln_like_den)) {
          pi = pi_star
          ln_like_den = ln_like_num
          gamma_p = gamma_p_star
        }
      }
    }
    ln_like[m] = ln_like_den
    sum_gamma_p[m] = sum(gamma_p)
    sum_ln_gamma_p[m] = -sum(log(gamma_p))
  }
  list(sum_gamma_p = sum_gamma_p, sum_ln_gamma_p = sum_ln_gamma_p,
       ln_like = ln_like, final_pi = pi)
}
