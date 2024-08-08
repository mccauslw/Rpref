#' Compute max-min log likelihood value
#'
#' Compute the max-min log likelihood value, which is both the Random Choice
#' (RC) log likelihood evaluated at the RC model where all choice distributions
#' are discrete uniform, and the Random Preference (RP) log likelihood evaluated
#' at the RP model where all preferences are equally likely. See Section 4.1
#' of the reference below.
#'
#' @template param-u
#' @template param-N
#'
#' @return Max-min log likelihood value
#'
#' @export
#'
#' @examples
#' library(RanCh)
#' n <- 5
#' u <- create_universe(n)
#' N <- vectorize_counts(u, RanCh::MMS_2019_counts[1, , ])
#' compute_max_min_ln_marl(u, N)
#'
#' @template SMC_reference
#' @template author
compute_max_min_ln_marl <- function(u, N) {
  ln_marl <- 0
  for (A in 1:u$n_subsets) {
    if (u$A_table[A, 'nA'] > 1) {
      ln_marl <- ln_marl - N$by_A[A] * log(u$A_table[A, 'nA'])
    }
  }
  ln_marl
}

#' Compute Random Choice (RC) model P maximising the RC likelihood and the
#' maximum value of the RC log likelihood.
#'
#' @template param-u
#' @template param-N
#'
#' @return A list with the following elements
#' \describe{
#'  \item{ln_maxl}{Maximum value of RC log likelihood}
#'  \item{P_mle}{RC model P maximizing the RC log likelihood}
#' }
#'
#' @export
#'
#' @examples
#' library(RanCh)
#' n <- 5
#' u <- create_universe(n)
#' N <- vectorize_counts(u, RanCh::MMS_2019_counts[1, , ])
#' compute_max_min_ln_marl(u, N)
#'
#' @template SMC_reference
#' @template author
compute_P_ln_maxl <- function(u, N) {
  P_mle <- rep(0, u$n_probs)
  ln_maxl <- 0
  for (A in 1:u$n_subsets) {
    if (u$A_table[A, 'nA'] > 1) {
      Ax <- u$A_table[A, 'Ax']
      nA <- u$A_table[A, 'nA']
      N_Ax <- N$by_Ax[Ax:(Ax+nA-1)]
      N_A <- N$by_A[A]
      P_Ax <- N_Ax / N_A
      P_mle[Ax:(Ax+nA-1)] <- P_Ax
      ln_maxl <- ln_maxl + sum(N_Ax[N_Ax!=0]*log(P_Ax[N_Ax!=0]))
    }
  }
  names(P_mle) <- rownames(u$Ax_table)
  list(ln_maxl = ln_maxl, P_mle = P_mle)
}

#' Compute the log marginal likelihood for the uniform-P model
#'
#' Compute the log marginal likelihood for the Bayesian Random Choice model
#' where choice probability vectors are mutually independent and each is
#' uniformly distributed on its corresponding probability simplex.
#'
#' @template param-u
#' @template param-N
#'
#' @return The value of the log marginal likelihood
#'
#' @export
#'
#' @examples
#' library(RanCh)
#' n <- 5
#' u <- create_universe(n)
#' N <- vectorize_counts(u, RanCh::MMS_2019_counts[1, , ])
#' compute_uniform_P_ln_marl(u, N)
#'
#' @template SMC_reference
#' @template author
compute_uniform_P_ln_marl <- function(u, N) {
  alpha_Ax <- matrix(1, nrow=nrow(u$Ax_table), ncol=1)
  rownames(alpha_Ax) <- rownames(u$Ax_table)
  ln_Pr_RC_by_A <- compute_ln_Pr_by_A(u, 'RC', alpha_Ax, N)
  ln_marl <- sum(ln_Pr_RC_by_A)
}

#' Use SMC to simulate the target distribution alpha|N of the Dirichlet RC model
#'
#' Use SMC to simulate the target distribution alpha|N, the posterior distribution
#' alpha|N of the scalar parameter alpha of the Dirichlet RC model.
#'
#' @template param-u
#' @template param-J
#' @template param-M
#' @template param-alpha_prior
#' @template param-N
#'
#' @return A list with the following elements
#' \describe{
#'    \item{alpha}{A M/J x J matrix with a sample from the target distribution}
#'    \item{marl_stats}{Simulation statistics for the marginal likelihood}
#'    \item{n_plus}{Number of non-zero counts in data}
#'    \item{theta}{Parameter of Beta-Gamma importance distribution}
#' }
#' @export
#'
#' @importFrom extraDistr rbetapr dbetapr
#'
#' @examples
#' library(RanCh)
#' n <- 5
#' u <- create_universe(n)
#' alpha_prior <- create_alpha_prior(n, 4, 0.1)
#' N <- vectorize_counts(u, RanCh::MMS_2019_counts[1, , ])
#' RC_sim <- run_RC_sim(u, 20, 1000, alpha_prior, N)
run_RC_sim <- function(u, J, M, alpha_prior, N) {

  # Compute alpha proposal distribution based on f(alpha) * Pr[N|alpha,lambda=0]
  params <- compute_proposal_params(u, alpha_prior, N)

  # Draw alpha from gamma proposal distribution
  alpha <- extraDistr::rbetapr(M, params[1], params[2], params[3])
  alpha_Ax <- compute_alpha_Ax(u, alpha)
  ln_Pr_RC_by_A <- compute_ln_Pr_by_A(u, 'RC', alpha_Ax, N)

  # Compute IS (importance sampling) log weights w
  ln_num <- stats::dgamma(alpha, alpha_prior$a, alpha_prior$b, log = TRUE)
  ln_num <- ln_num + compute_ln_like(u, 0.0, ln_Pr_RC_by_A, NULL)
  ln_den <- extraDistr::dbetapr(alpha, params[1], params[2], params[3], log = TRUE)
  w <- ln_num - ln_den # log weights

  # Subtract an offset to the log weights to avoid underflow, and incorporate
  # offset into the cumulative log marginal likelihood
  ln_marl_offset <- max(w) # Compute an offset to avoid underflow
  gr_cum_ln_marl <- rep(ln_marl_offset, J)
  w <- w - ln_marl_offset
  W <- exp(w)

  # Compute IS marginal likelihood for lambda = 0, taking into account offset
  C_stage_stats <- weights_to_C_stage_stats(W, ln_marl_offset, 0, gr_cum_ln_marl)

  # Initial selection step to get an unweighted sample
  W = matrix(W, nrow=M/J, ncol=J)
  for (j in 1:J) {
    probabilities = W[,j] / sum(W[,j])
    selection = sample.int(M/J, M/J, replace = TRUE, prob = probabilities)
    alpha[(1+(j-1)*M/J):(j*M/J)] <- alpha[(selection+(j-1)*(M/J))]
  }

  n_plus = sum(N$by_Ax > 0)
  list(alpha=alpha, marl_stats = C_stage_stats, n_plus = n_plus, theta = params)
}

# Compute sample statistics and standard errors for a sample targeting
# a distribution.
# The sample x is a vector of length M, with J mutually independent
# subvectors, each of size M/J.
# The inputs are
# - x, the sample in the form of a vector of length M
# - J, the number of independent subvectors (i.e. the number of groups of particles)
# - p, a vector of probabilities for which to compute quantiles
# The output is a list with fields
# - mu, sample mean
# - std, sample standard deviation
# - nse, numerical standard error for the sample mean
# - rne, relative numerical efficiency for the sample mean
# - p, the vector of probabilities given as input
# - q, the vector of corresponding quantiles
# - q_nse, a vector of numerical standard errors for the quantiles
ind_groups_stats <- function(x, J, p) {
  M <- length(x)
  s2 <- stats::var(x); std = sqrt(s2)
  # x_by_group is a (M/J) x J matrix, x organized as J independent columns, with
  # each column generated by a group of particles.
  x_by_group <- matrix(x, nrow=M/J, ncol=J)
  mu_by_group <- colMeans(x_by_group)    # J-vector of group means
  mu <- mean(mu_by_group)                # Global mean
  nse2 <- stats::var(mu_by_group)/J; nse <- sqrt(nse2)
  q_by_group <- apply(x_by_group, 2, stats::quantile, probs = p)
  q <- rowMeans(q_by_group)
  q_nse2 = apply(q_by_group, 1, stats::var)/J; q_nse <- sqrt(q_nse2)
  list(mu = mu, std = std, nse = nse, rne = (s2/M)/nse2,
       p = p, q = q, q_nse = q_nse)
}

# Compute marginal likelihood statistics for the vector of importance weights
# from a single value of lambda
weights_to_C_stage_stats <- function(W, cum_ln_marl, cum_ln_marl_nse2,
                                     gr_cum_ln_marl) {
  M <- length(W); J <- length(gr_cum_ln_marl)

  # To compute globally, not by group
  W_var <- stats::var(W)      # Global variance of weights
  ESS <- sum(W)^2 / sum(W^2)  # Effective sample size of weights

  # Organize weights by group, compute group means
  W <- matrix(W, nrow=M/J, ncol=J)      # Weights organized by groups

  # To compute by group, j=1,...,J
  gr_ESS <- colSums(W)^2 / colSums(W^2) # ESS
  gr_marl <- colMeans(W)             # Marginal likelihood factor
  # Cumulative log marginal likelihood
  gr_cum_ln_marl <- gr_cum_ln_marl + log(gr_marl)

  # Final computations
  marl <- mean(gr_marl)              # Marginal likelihood increment
  marl_nse2 <- stats::var(gr_marl)/J # Numerical variance of mlike increment
  marl_nse <- sqrt(marl_nse2)        # Numerical standard error of mlike increment
  marl_rne <- (W_var/M) / marl_nse2  # Relative numerical efficiency of ml
  ln_marl <- log(marl)               # Log marginal likelihood increment
  ln_marl_nse <- marl_nse / marl     # Numerical standard error of log mlike increment
  list(gr_ESS = gr_ESS, gr_marl = gr_marl, gr_cum_ln_marl = gr_cum_ln_marl,
       ESS = ESS, W_var = W_var,
       marl = marl, marl_nse = marl_nse, marl_rne = marl_rne,
       ln_marl = ln_marl, ln_marl_nse = ln_marl_nse,
       cum_ln_marl = cum_ln_marl + ln_marl,
       cum_ln_marl_nse2 = cum_ln_marl_nse2 + ln_marl_nse^2)
}

# Compute log multinomial density for sequence, rather than counts
lnf_seq_multinom <- function(p, N_x) {
  sum(N_x[N_x != 0] * log(p[N_x != 0]))
}

# Compute expectation of maximum likelihood P at maximum likelihood probabilities
compute_EP_ln_maxl <- function(u, N, alpha, M) {
  P_mle <- rep(0, u$n_probs)
  Eln_maxl <- 0
  for (A in 1:u$n_subsets) {
    if (u$A_table[A, 'nA'] > 1) {
      Ax <- u$A_table[A, 'Ax']
      nA <- u$A_table[A, 'nA']
      N_Ax <- N$by_Ax[Ax:(Ax+nA-1)]
      N_A <- N$by_A[A]
      al_post = N_Ax + rep(alpha/nA, nA)
      P = bayesm::rdirichlet(M, N_Ax)
      lnf = apply(P, MARGIN=1, lnf_seq_multinom, N_Ax)
      Eln_maxl <- Eln_maxl + mean(lnf)
    }
  }
  Eln_maxl
}

# Aggregated preference gamma weights to choice gamma weights
compute_gamma_Ax <- function(u, gamma) {
  gamma_Ax <- u$pi_to_P %*% gamma
  rownames(gamma_Ax) <- u$Ax_strings
  gamma_Ax
}

# Compute choice alpha parameters from scalar alpha parameter
compute_alpha_Ax <- function(u, alpha) {
  M <- ncol(alpha)
  card_Ax <- u$A_table[u$Ax_table[, 'A'], 'nA']
  alpha_Ax <- outer(card_Ax, alpha, FUN = function(x,y){y/x})
  rownames(alpha_Ax) <- rownames(u$Ax_table)
  alpha_Ax
}

log_prior_gamma <- function(u, gamma_p, alpha_p) {
  colSums(dgamma(gamma_p, alpha_p, log = TRUE))
}

# Compute one of the following log probabilities, by choice set
#
#   (1) ln_Pr_RC_by_A, n_subsets x M
#     - gives, for set A (row) and particle m (col), the probability of choice
#       count vector N_A given alpha_m and lambda = 0
#   (2) ln_Pr_RP_by_A, n_subsets x M
#     - gives, for set A (row) and particle m (col), the probability of choice
#       count vector N_A given gamma_p (col m) and lambda = 1
#
compute_ln_Pr_by_A <- function(u, type, weight_Ax, N) {
  ln_Pr_by_A <- matrix(0, nrow=u$n_subsets, ncol=ncol(weight_Ax))
  rownames(ln_Pr_by_A) <- u$A_strings

  if (type == 'RP')
    gamma_total <- weight_Ax[1,] + weight_Ax[2,]

  for (A in 1:u$n_subsets) {
    if (u$A_table[A, 'nA'] > 1) {
      Ax <- u$A_table[A, 'Ax']
      nA <- u$A_table[A, 'nA']
      N_Ax <- N$by_Ax[Ax:(Ax+nA-1)]
      N_A <- N$by_A[A]
      if (type == 'RP') { # weight_Ax is gamma_Ax
        ga_Ax <- weight_Ax[Ax:(Ax+nA-1), ]
        ln_Pr_by_A[A, ] <- colSums(N_Ax*log(ga_Ax)) - N_A*log(gamma_total)
      }
      if (type == 'RC') { # weight_Ax is alpha_Ax
        al_Ax <- as.matrix(weight_Ax[Ax:(Ax+nA-1), ])
        ln_Pr_by_A[A, ] <-
          (lgamma(colSums(al_Ax)) - colSums(lgamma(al_Ax))
           - lgamma(colSums(al_Ax) + N_A) + colSums(lgamma(al_Ax + N_Ax)))
      }
    }
  }
  ln_Pr_by_A
}

# If lambda is the scalar zero, ln_Pr_RP_by_A is not referenced
compute_ln_like <- function(u, lambda, ln_Pr_RC_by_A, ln_Pr_RP_by_A) {
  if (identical(lambda, 0.00)) {
    ln_like <- colSums(ln_Pr_RC_by_A)
  } else {
    n_lambda <- length(lambda)
    ln_like <- matrix(0, nrow=ncol(ln_Pr_RP_by_A), ncol=n_lambda)
    for (i in 1:n_lambda) {
      ln_like[, i] <-
        colSums(log(
          lambda[i] * exp(ln_Pr_RP_by_A) + (1-lambda[i]) * exp(ln_Pr_RC_by_A)
        ))
    }
  }
  ln_like
}

AR_gamma <- function(gamma, alpha, phi) {
  # gamma is length M,
  # alpha is length M,
  # phi is scalar
  M <- length(alpha)
  be <- stats::rbeta(M, phi*alpha, (1-phi)*alpha)
  gamma <- be * gamma + stats::rgamma(M, (1-phi)*alpha)
  gamma
}

update_alpha <- function(u, N, alpha_prior, alpha, gamma_p, lambda, ln_Pr_RP_by_A,
                         n_reps = 2) {
  M <- length(alpha)
  # Metropolis-Hastings proposal distribution for target distribution
  # alpha|pi is based on the close approximation gamma(x) = 1/x for
  # small values of x. The proposal distribution is Ga(a_post, b_post)

  # First compute a_post, b_post and draw alpha_star
  g_bar <- pmax(colMeans(log(gamma_p)), log(.Machine$double.xmin))
  G0 <- colSums(gamma_p)
  p_bar <- g_bar - log(G0)
  n_fact = u$n_orders

  alpha_mode <- alpha_prior$alpha_mode[(p_bar - alpha_prior$p_min)/alpha_prior$h]
  a_delta <- -3
  a_bar <- alpha_prior$a + n_fact - 1 + a_delta
  b_delta <- -alpha_prior$psi_diff[(p_bar - alpha_prior$p_min)/alpha_prior$h] +
    a_delta/alpha_mode
  b_bar <- alpha_prior$b - p_bar + b_delta

  alpha_Ax <- compute_alpha_Ax(u, alpha)
  ln_Pr_RC_by_A <- compute_ln_Pr_by_A(u, 'ind', alpha_Ax, N)
  ln_ll <- compute_ln_like(u, lambda, ln_Pr_RC_by_A, ln_Pr_RP_by_A)
  ln_f_over_g <- lgamma(alpha + 1) -
    n_fact * lgamma(alpha/n_fact + 1) -
    a_delta * log(alpha) + b_delta * alpha

  aPr <- 0
  for (i_rep in seq(n_reps)) {
    # Next compute ln Pr[N|alpha_star, lambda = 0], ln Pr[N|alpha, lambda = 0]
    alpha_star <- stats::rgamma(M, a_bar, b_bar)
    alpha_Ax_star <- compute_alpha_Ax(u, alpha_star)
    ln_Pr_RC_by_A_star <- compute_ln_Pr_by_A(u, 'ind', alpha_Ax_star, N)
    ln_ll_star <- compute_ln_like(u, lambda, ln_Pr_RC_by_A_star, ln_Pr_RP_by_A)

    # Next evaluate proposal density g and target density f at alpha and alpha_star
    ln_f_over_g_star <- lgamma(alpha_star + 1) -
      n_fact * lgamma(alpha_star/n_fact + 1) -
      a_delta * log(alpha_star) + b_delta * alpha_star

    # Finally, evaluate log Hasting ratio, accept proposals with probability
    # min(1, H) and return updated alpha.
    ln_H <- ln_ll_star - ln_ll + ln_f_over_g_star - ln_f_over_g
    H = pmin(1, exp(ln_H))
    accept <- stats::runif(M) <= H
    aPr <- aPr + mean(H)
    alpha[accept] <- alpha_star[accept]
    ln_f_over_g[accept] <- ln_f_over_g_star[accept]
    ln_ll[accept] <- ln_ll_star[accept]
  }
  G = stats::rgamma(M, alpha)
  gamma_p <- scale(gamma_p, center=F, scale=G0/G)
  aPr <- aPr / n_reps
  list(alpha = alpha, gamma_p = gamma_p, aPr = aPr, mu = mean(alpha))
}
