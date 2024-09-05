#' Compute max-min log likelihood value
#'
#' Compute the max-min log likelihood value, which is both the Random Choice
#' (RC) log likelihood evaluated at the RC model where all choice distributions
#' are discrete uniform, and the Random Preference (RP) log likelihood evaluated
#' at the RP model where all preferences are equally likely. See Section 4.1
#' of the reference below.
#'
#' @inheritParams compute_pi_ln_like
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
#' @inherit create_universe author references
#'
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
#' @inheritParams compute_pi_ln_like
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
#' @inherit create_universe author references
#'
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
#' @inheritParams compute_pi_ln_like
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
#' @inherit create_universe author references
#'
compute_uniform_P_ln_marl <- function(u, N) {
  alpha_Ax <- matrix(1, nrow=nrow(u$Ax_table), ncol=1)
  rownames(alpha_Ax) <- rownames(u$Ax_table)
  ln_Pr_RC_by_A <- compute_ln_Pr_by_A(u, 'RC', alpha_Ax, N)
  ln_marl <- sum(ln_Pr_RC_by_A)
}

#' Use SMC to simulate the posterior distribution alpha|N of the Dirichlet RC model
#'
#' Use SMC to simulate the target distribution alpha|N, the posterior distribution
#' alpha|N of the scalar parameter alpha of the Dirichlet RC model.
#'
#' @inheritParams compute_pi_ln_like
#' @param J number of independent groups of particles
#' @param M umber of particles, must be multiple of J
#' @param alpha_prior list with information on the prior distribution of alpha,
#'                    created using [create_alpha_prior()]
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
#' J <- 20
#' M <- 1000
#' RC_sim <- run_RC_sim(u, J, M, alpha_prior, N)
#'
#' @inherit create_universe author references
#'
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

#' Use SMC to simulate posterior distributions of Dirichlet RC and hybrid models
#'
#'
#' @inheritParams run_RC_sim
#' @param lambda_values Vector of indices of the hybrid models
#' @param cycle_schedule Schedule of parameters used for the various SMC cycles
#'
#' @return List with the following elements
#' \describe{
#'    \item{alpha}{A M/J x J matrix with a sample of alpha from the final target distribution}
#'    \item{gamma}{?WJM? sample of gamma from the final target distribution}
#'    \item{lambda_stats}{?WJM? simulation statistics}
#'    \item{cycle_stats}{?WJM? simulation statistics}
#'    \item{big_aPr}{Acceptance probabilities for big blocks}
#'    \item{sm_aPr}{Acceptance probabilities for small blocks}
#'    \item{alpha_aPr}{Acceptance probabilities for alpha proposals}
#'    \item{alpha_mu}{?WJM?}
#' }
#'
#' @export
#'
#' @importFrom stats rgamma runif
#'
#' @examples
#' library(RanCh)
#' n <- 5
#' u <- create_universe(n)
#' alpha_prior <- create_alpha_prior(n, 4, 0.1)
#' N <- vectorize_counts(u, RanCh::MMS_2019_counts[1, , ])
#' J <- 20
#' M <- 1000
#' RP_sim <- run_RP_sim(u, J, M, alpha_prior, N, lambda_values, cycle_schedule)
#'
#' @inherit create_universe author references
#'
run_RP_sim <- function(u, J, M, alpha_prior, N, lambda_values, cycle_schedule) {

  lambda_aggregate_names =
    c('ESS', 'W_var', 'marl', 'marl_nse', 'marl_rne',
      'ln_marl', 'ln_marl_nse', 'cum_ln_marl', 'cum_ln_marl_nse2')
  n_lambda_values = length(lambda_values)
  n_cycles = nrow(cycle_schedule)

  # Table to hold marginal likelihood statistics
  lambda_stats <- list(
    gr_ESS = matrix(nrow = J, ncol = n_lambda_values),
    gr_marl = matrix(nrow = J, ncol = n_lambda_values),
    gr_cum_ln_marl = matrix(nrow = J, ncol = n_lambda_values),
    aggregates =
      matrix(nrow = n_lambda_values, ncol = 1 + length(lambda_aggregate_names),
             dimnames = list(NULL, c('lambda', lambda_aggregate_names)))
  )
  lambda_stats$aggregates[, 'lambda'] <- lambda_values

  # Table for cycle parameters pertaining to gamma update
  n_bl <- c(u$n, u$n*(u$n-1)); bl_len <- c(factorial(u$n-1), factorial(u$n-2))
  bl_data <- list(big = list(), sm = list())
  for (blt in 1:length(n_bl)) {
    bl_data[[blt]]$n_bl = n_bl[blt]
    bl_data[[blt]]$bl_len = bl_len[blt]
    bl_data[[blt]]$indices <- seq(n_bl[blt])
    bl_data[[blt]]$start <- (seq(n_bl[blt]) - 1) * bl_len[blt] + 1
    bl_data[[blt]]$end <- seq(n_bl[blt]) * bl_len[blt]
    bl_data[[blt]]$aPr <- matrix(NA, nrow = n_cycles, ncol = n_bl[blt])
  }

  # alpha information by cycle
  alpha_mu <- rep(NA, n_cycles)
  alpha_aPr <- rep(NA, n_cycles)

  # Start with simulation based on lambda = 0
  RC_sim <- run_RC_sim(u, J, M, alpha_prior, N)
  alpha <- RC_sim$alpha
  gr_cum_ln_marl <- RC_sim$marl_stats$gr_cum_ln_marl
  cum_ln_marl <- RC_sim$marl_stats$cum_ln_marl
  cum_ln_marl_nse2 <- RC_sim$marl_stats$cum_ln_marl_nse2

  # Compute quantities depending on alpha
  alpha_p <- outer(rep.int(1/u$n_orders, u$n_orders), alpha)
  rownames(alpha_p) <- u$order_strings
  alpha_Ax <- compute_alpha_Ax(u, alpha)
  rownames(alpha_Ax) <- u$Ax_strings
  ln_Pr_RC_by_A <- compute_ln_Pr_by_A(u, 'RC', alpha_Ax, N)

  # Draw gamma_p from gamma_p|alpha, same as gamma_p|alpha, N, lambda = 0
  # Columns are iid, the rows are independent.
  # An element i,j is distributed as Gamma(alpha_i, 1), where alpha_i
  # is the i'th element of alpha_p.
  gamma_p <- matrix(stats::rgamma(u$n_orders * M, alpha_p),
                    nrow=u$n_orders, ncol=M)
  rownames(gamma_p) <- u$order_strings
  gamma_Ax <- compute_gamma_Ax(u, gamma_p)
  ln_Pr_rp_by_A <- compute_ln_Pr_by_A(u, 'RP', gamma_Ax, N)

  # Initialize loop over cycles
  old_ll <- compute_ln_like(u, 0.0, ln_Pr_RC_by_A, NULL)
  first_lambda_index = 1

  # Loop over S-C-M cycles
  for (cycle_index in 1:n_cycles) {

    n_sweeps <- c(cycle_schedule[[cycle_index, 'n_big_sweeps']],
                  cycle_schedule[[cycle_index, 'n_sm_sweeps']])
    phi_sweeps <- c(cycle_schedule[[cycle_index, 'phi_big_sweeps']],
                    cycle_schedule[[cycle_index, 'phi_sm_sweeps']])
    last_lambda_index <- cycle_schedule[[cycle_index, 'lambda_break_indices']]
    lambda_index_range <- seq(first_lambda_index, last_lambda_index)
    cycle_lambda_values <- lambda_values[lambda_index_range]
    n_lambda_index <- length(lambda_index_range)
    last_lambda <- cycle_lambda_values[n_lambda_index]

    # C stage: compute importance weights
    #####################################

    ll <- compute_ln_like(u, cycle_lambda_values, ln_Pr_RC_by_A, ln_Pr_rp_by_A)
    w <- ll - as.vector(old_ll)

    for (lambda_i in lambda_index_range) {
      single_lambda_stats <-
        weights_to_C_stage_stats(exp(w[, lambda_i - first_lambda_index + 1]),
                                 cum_ln_marl, cum_ln_marl_nse2, gr_cum_ln_marl)
      for (name in lambda_aggregate_names)
        lambda_stats$aggregates[lambda_i, name] <- single_lambda_stats[[name]]
      lambda_stats$gr_ESS[, lambda_i] <- single_lambda_stats$gr_ESS
      lambda_stats$gr_marl[, lambda_i] <- single_lambda_stats$gr_marl
      lambda_stats$gr_cum_ln_marl[, lambda_i] <- single_lambda_stats$gr_cum_ln_marl
    }

    # Values at last lambda value of cycle
    gr_cum_ln_marl <- single_lambda_stats$gr_cum_ln_marl
    cum_ln_marl <- single_lambda_stats$cum_ln_marl
    cum_ln_marl_nse2 <- single_lambda_stats$cum_ln_marl_nse2

    # Selection (resampling) stage
    ##############################

    # Resample alpha and gamma using final weights, group by group
    W = matrix(exp(w[, n_lambda_index]), nrow=M/J, ncol=J)
    for (j in 1:J) {
      # Old multinomial sampling
      #probabilities = W[,j] / sum(W[,j])
      #if (any(is.nan(probabilities))) {
      #  print(cycle_index)
      #  print(summary(alpha))
      #  print(summary(colSums(gamma)))
      #}
      #selection <- sample.int(M/J, M/J, replace=TRUE, prob=probabilities)

      # New residual sampling
      En <- W[,j] / mean(W[,j])
      indices = seq(M/J)
      sure_n = floor(En)
      resid_n = En - sure_n
      sure_M = sum(sure_n)
      resid_M = M/J - sure_M
      selection <- c(rep(seq(M/J), sure_n),
                     sample.int(M/J, resid_M, replace=TRUE, prob=resid_n))

      gamma_p[,(1+(j-1)*M/J):(j*M/J)] <- gamma_p[,(selection+(j-1)*(M/J))]
      alpha[(1+(j-1)*M/J):(j*M/J)] <- alpha[(selection+(j-1)*(M/J))]
    }
    alpha_p = outer(rep.int(1/u$n_orders, u$n_orders), alpha)

    # Mutation stage
    ################

    # Instead of recomputing gamma_Ax, alpha_Ax, ln_Pr_rp_by_A,
    # construct at resampling step
    gamma_Ax <- compute_gamma_Ax(u, gamma_p)
    alpha_Ax <- compute_alpha_Ax(u, alpha)
    ln_Pr_RC_by_A <- compute_ln_Pr_by_A(u, 'RC', alpha_Ax, N)
    ln_Pr_rp_by_A <- compute_ln_Pr_by_A(u, 'RP', gamma_Ax, N)
    den <- compute_ln_like(u, last_lambda, ln_Pr_RC_by_A, ln_Pr_rp_by_A)

    for (blt in 1:length(bl_data)) {
      bld <- bl_data[[blt]]
      bl_data[[blt]]$aPr[cycle_index, ] <- 0
      phi <- phi_sweeps[blt]
      for (i_sweep in 1:n_sweeps[blt]) {
        for (i_bl in 1:(bld$n_bl)) {

          # Redraw a random sample of gamma vectors
          small_gamma_indices <- seq(bld$start[i_bl], bld$end[i_bl])
          small_pi_to_P <- u$pi_to_P[, small_gamma_indices]
          small_alpha_p <- alpha_p[small_gamma_indices, ]
          small_alpha_total <- colSums(small_alpha_p)
          small_gamma_p <- gamma_p[small_gamma_indices, ]

          small_gamma_total <- colSums(small_gamma_p)
          AR_small_gamma_total <-
            AR_gamma(small_gamma_total, small_alpha_total, phi)
          small_gamma_p_star <- matrix(rgamma(bld$bl_len * M, small_alpha_p),
                                       nrow=bld$bl_len, ncol=M)
          small_gamma_total_star <-
            pmax(colSums(small_gamma_p_star), .Machine$double.xmin)
          scale_value <- small_gamma_total_star/AR_small_gamma_total
          small_gamma_p_star <-
            scale(small_gamma_p_star, center=F, scale=scale_value)
          small_gamma_p_star <- pmax(small_gamma_p_star, .Machine$double.xmin)

          # Recompute likelihood values for proposal, accept or reject
          gamma_Ax_diff <- small_pi_to_P %*% (small_gamma_p_star-small_gamma_p)
          gamma_Ax_star <- pmax(gamma_Ax + gamma_Ax_diff, .Machine$double.xmin)
          ln_Pr_rp_by_A_star <- compute_ln_Pr_by_A(u, 'RP', gamma_Ax_star, N)
          num <- compute_ln_like(u, last_lambda, ln_Pr_RC_by_A,
                                 ln_Pr_rp_by_A_star)
          H <- pmin(1, exp(num-den))
          accept <- (runif(M) < H)
          gamma_p[small_gamma_indices, accept] <- small_gamma_p_star[, accept]
          gamma_Ax <- compute_gamma_Ax(u, gamma_p)
          den[accept] <- num[accept]
          bl_data[[blt]]$aPr[cycle_index, i_bl] <-
            bl_data[[blt]]$aPr[cycle_index, i_bl] + mean(H)
        }
      }
      bl_data[[blt]]$aPr[cycle_index, ] <-
        bl_data[[blt]]$aPr[cycle_index, ] / n_sweeps[blt]
    }
    ln_Pr_rp_by_A <- compute_ln_Pr_by_A(u, 'RP', gamma_Ax, N)

    # Update alpha, update likelihood for next cycle
    res <- update_alpha(u, N, alpha_prior, alpha, gamma_p, last_lambda, ln_Pr_rp_by_A)
    alpha <- res$alpha
    gamma_p <- res$gamma_p
    alpha_aPr[cycle_index] = res$aPr
    alpha_mu[cycle_index] = res$mu
    alpha_Ax <- compute_alpha_Ax(u, alpha)
    alpha_p <- outer(rep.int(1/u$n_orders, u$n_orders), alpha)
    ln_Pr_RC_by_A <- compute_ln_Pr_by_A(u, 'RC', alpha_Ax, N)
    old_ll <- compute_ln_like(u, last_lambda, ln_Pr_RC_by_A, ln_Pr_rp_by_A)

    first_lambda_index = last_lambda_index + 1
  }

  list(alpha = alpha, gamma = gamma_p,
       lambda_stats = lambda_stats,
       cycle_stats = lambda_stats$aggregates[cycle_schedule$lambda_break_indices,],
       big_aPr = bl_data[[1]]$aPr, sm_aPr = bl_data[[2]]$aPr,
       alpha_aPr = alpha_aPr, alpha_mu = alpha_mu)
}

#' Compute sample statistics and standard errors for a collection of SMC samples
#'
#' Given a collection of J independent Sequential Monte Carlo (SMC) samples,
#' compute mean, standard deviation, numerical standard error for the mean,
#' relative numerical efficiency for the mean, specified quantiles and their
#' numerical standard errors.
#'
#' @param x vector of length M, with J mutually independent subvectors of length
#'  M/J, each subvector arising from a SMC simulation
#' @inheritParams run_RC_sim
#' @param p vector of probabilities for which to compute quantiles
#'
#' @return A list with the elements
#' \describe{
#'   \item{mu}{sample mean}
#'   \item{std}{sample standard deviation (measure of posterior uncertainty)}
#'   \item{nse}{numerical standard error for the sample mean (measure of
#'   simulation noise)}
#'   \item{rne}{relative numerical efficiency for the sample mean}
#'   \item{p}{same as input of that name}
#'   \item{q}{sample quantiles corresponding to p}
#'   \item{q_nse}{numerical standard errors for the quantiles}
#' }
#' @export
#'
#' @examples
#' library(RanCh)
#' n <- 5
#' u <- create_universe(n)
#' alpha_prior <- create_alpha_prior(n, 4, 0.1)
#' N <- vectorize_counts(u, RanCh::MMS_2019_counts[1, , ])
#' J <- 20
#' M <- 1000
#' RC_sim <- run_RC_sim(u, J, M, alpha_prior, N)
#' ind_groups_stats(RC_sim$alpha)
#'
#' @inherit create_universe author references
#'
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

#' Compute pdf and cdf on the grid x_grid for a collection of SMC samples
#'
#' Given a collection of J independent Sequential Monte Carlo (SMC) samples,
#' compute the sample probability density function (pdf) and cumulative
#' distribution function (cdf), with numerical standard errors, at a grid
#' of points of evaluation
#'
#' @inheritParams ind_groups_stats
#' @param x_grid vector of points of evaluation
#'
#' @return A list with elements pdf and cdf, each with elements
#' \describe{
#'   \item{x}{vector of argument values}
#'   \item{func}{vector of function values}
#'   \item{nse}{vector of numerical standard errors for the function values}
#' }
#' @export
#'
#' @examples
compute_pdf_cdf_on_grid <- function(x, J, x_grid) {
  if (length(x_grid) == 1) {
    x_grid = seq(min(x), max(x), len=x_grid)
  }
  M <- length(x)
  Kp1 <- length(x_grid); K = Kp1 - 1
  x_by_group <- matrix(x, nrow=M/J, ncol=J)
  pdf_by_group <- matrix(nrow=K, ncol=J)
  cdf_by_group <- matrix(nrow=Kp1, ncol=J)
  for (j in 1:J) {
    h <- hist(x_by_group[,j], breaks = x_grid, plot=F)
    pdf_by_group[,j] <- h$density
    cdf_by_group[,j] <- c(0, cumsum(h$counts/sum(h$counts)))
  }
  pdf <- rowMeans(pdf_by_group)
  cdf <- rowMeans(cdf_by_group)
  pdf_nse <- sqrt(rowMeans((pdf_by_group-pdf)^2)/J)
  cdf_nse <- sqrt(rowMeans((cdf_by_group-cdf)^2)/J)
  list(pdf = list(x=h$mids, func=pdf, nse = pdf_nse),
       cdf = list(x=x_grid, func=cdf, nse = cdf_nse))
}

#' Compute the posterior pdf and cdf of binary choice probabilities for the
#' Dirichlet Random Preference model
#'
#' For each binary choice probability, compute its posterior pdf and cdf for the
#' Dirichlet RP model, on the grid p_grid of binary choice probabilities,
#' given a collection of SMC samples that target the posterior distribution
#' of gamma in the Dirichlet RP model
#'
#' @inheritParams ind_groups_stats
#' @inheritParams compute_pi_ln_like
#' @param gamma_p matrix, M by n!, where each column is a vector of
#' unnormalized preference probabilities
#' @param p_grid vector of choice probability values
#'
#' @return A list with elements correponding to binary subsets of choice
#' objects. For a universe of size n=5, there are n choose 2, or 10 such
#' subsets. Each element of the list is a list organized in the same way
#' as the return value of [compute_pdf_cdf_on_grid()].
#' @export
#'
#' @examples
#'
#' @inherit create_universe author references
#'
compute_RP_binp_funcs <- function(u, gamma_p, J, N, p_grid) {
  M <- ncol(gamma_p)
  n_grid_pts = length(p_grid)
  Ax_binaries <- u$A_table[u$A_table[,'nA']==2, 'Ax']
  n_binary <- length(Ax_binaries)
  binp_funcs <- vector("list", n_binary)

  # bin_probs is a n_binary x M matrix of binary choice probabilities
  bin_probs <- t(t(u$pi_to_P[Ax_binaries,] %*% gamma_p)/colSums(gamma_p))
  for (i_bin in seq(n_binary)) {
    binp_funcs[[i_bin]] <- compute_pdf_cdf_on_grid(bin_probs[i_bin,], J, p_grid)
  }
  binp_funcs
}

#' Compute the posterior pdf and cdf of binary choice probabilities for the
#' Dirichlet Random Choice model
#'
#' For each binary choice probability, compute its posterior pdf and cdf for the
#' Dirichlet RC model, on the grid p_grid of binary choice probabilities,
#' given a collection of SMC samples that target the posterior distribution
#' of alpha in the Dirichlet RC model
#'
#' @inheritParams compute_RP_binp_funcs
#' @param alpha vector of length M, sample of values of alpha
#'
#' @inherit compute_RP_binp_funcs return
#' @export
#'
#' @examples
#'
#' @inherit create_universe author references
#'
compute_RC_binp_funcs <- function(u, alpha, J, N, p_grid) {
  n_grid <- length(p_grid)
  M <- length(alpha)
  Ax_binaries <- u$A_table[u$A_table[,'nA']==2, 'Ax']
  n_binary <- length(Ax_binaries)
  n1 <- N$by_Ax[Ax_binaries]      # first object count for n_binary choice sets
  n0 <- N$by_Ax[Ax_binaries + 1]  # second object count for n_binary choice sets

  binp_funcs <- vector("list", n_binary)

  for (i_bin in seq(n_binary)) {
    a1v <- n1[i_bin] + 0.5 * alpha  # Vector of M Dirichlet weights, 1st object
    a2v <- n0[i_bin] + 0.5 * alpha  # Vector of M Dirichlet weights, 2nd object
    a1 <- matrix(a1v, nrow = M/J, ncol = J)
    a2 <- matrix(a2v, nrow = M/J, ncol = J)

    cdf <- matrix(0, nrow = n_grid, ncol = J)
    pdf <- matrix(0, nrow = n_grid, ncol = J)
    for (j in seq(J)) {
      for (m in seq(M/J)) {
        pdf <- pdf + dbeta(p_grid, a1[m, j], a2[m, j])
        cdf <- cdf + pbeta(p_grid, a1[m, j], a2[m, j])
      }
      pdf = pdf/(M/J); cdf <- cdf/(M/J)
    }
    pdf_mean = rowMeans(pdf); cdf_mean = rowMeans(cdf)
    binp_funcs[[i_bin]] <- list(
      pdf = list(x=p_grid, func=pdf_mean, nse=rowMeans((pdf - pdf_mean)^2/J)),
      cdf = list(x=p_grid, func=cdf_mean, nse=rowMeans((cdf - cdf_mean)^2/J))
    )
  }
  binp_funcs
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
  ln_Pr_RC_by_A <- compute_ln_Pr_by_A(u, 'RC', alpha_Ax, N)
  ln_ll <- compute_ln_like(u, lambda, ln_Pr_RC_by_A, ln_Pr_RP_by_A)
  ln_f_over_g <- lgamma(alpha + 1) -
    n_fact * lgamma(alpha/n_fact + 1) -
    a_delta * log(alpha) + b_delta * alpha

  aPr <- 0
  for (i_rep in seq(n_reps)) {
    # Next compute ln Pr[N|alpha_star, lambda = 0], ln Pr[N|alpha, lambda = 0]
    alpha_star <- stats::rgamma(M, a_bar, b_bar)
    alpha_Ax_star <- compute_alpha_Ax(u, alpha_star)
    ln_Pr_RC_by_A_star <- compute_ln_Pr_by_A(u, 'RC', alpha_Ax_star, N)
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
  G <- stats::rgamma(M, alpha)
  gamma_p <- scale(gamma_p, center=F, scale=G0/G)
  aPr <- aPr / n_reps
  list(alpha = alpha, gamma_p = gamma_p, aPr = aPr, mu = mean(alpha))
}
