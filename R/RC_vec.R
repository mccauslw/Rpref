#' Vectorize a choice count matrix
#'
#' Create a list with two vectors of choice counts
#'
#' @inheritParams compute_pi_ln_like
#' @param choice_counts Matrix of choice counts indexed by choice set (row) and
#' choice object (column)
#'
#' @return List with the following elements
#' \describe{
#'   \item{by_Ax}{Vector of choice counts, indexed by Ax, a flat index coding
#'   both choice set and choice object}
#'   \item{by_A}{Vector of choice counts, indexed by choice set
#'   (aggregates by_Ax)}
#' }
#' @export
#'
#' @examples
#' n <- 5
#' u <- create_universe(n)
#' N <- vectorize_counts(u, RanCh::MMS_2019_counts[1, , ])
vectorize_counts <- function(u, choice_counts) {
  N <- list(by_Ax = rep.int(0, u$n_probs),
            by_A = rep.int(0, u$n_subsets))
  for (Ax in 1:u$n_probs) {
    A <- u$Ax_table[Ax,'A']
    x <- u$Ax_table[Ax,'x']
    N$by_Ax[Ax] <- choice_counts[A, x]
    N$by_A[A] <- N$by_A[A] + choice_counts[A, x]
  }
  names(N$by_Ax) <- u$Ax_strings
  names(N$by_A) <- u$A_strings
  N
}
