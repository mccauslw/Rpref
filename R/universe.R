#' Precompute tables and matrices depending only on n, the number of objects
#'
#' For a choice universe of size n, precompute a list of frequently used
#' quantities, including tables and convenient row names. The result is a
#' required input for many functions in this package.
#'
#' @param n Number of elements in choice universe
#'
#' @return A list with the following elements
#'\describe{
#'   \item{n}{Same as input with that name.}
#'   \item{n_orders}{Number of preference orders, equal to n!}
#'   \item{n_subsets}{Number of non-empty subsets of choice universe}
#'   \item{n_probs}{Total number of choice probabilities, over all doubleton and larger subsets}
#'   \item{orders}{List of preference orders, with each order a numeric n-vector}
#'   \item{order_strings}{List of strings (e.g. 'abc', 'acb') representing preference orders}
#'   \item{inv_orders}{List of inverse (in the sense of permutation inverse) of preference orders}
#'   \item{singletons}{Vector of singletons, with sets coded as bit strings}
#'   \item{Ax_table}{n_probs by 2 matrix. Each row gives a unique pair of choice set and element}
#'   \item{Ax_strings}{Vector of strings naming rows of Ax_table; the string 'ab;a', for example,
#'   identifies the pair consisting of choice set \{a,b\} and choice object a}
#'   \item{A_table}{n_subsets by 2 matrix. Each row gives the number of choice probabilities associated
#'   with the choice subset and the index of the first row of Ax_table where these choice probabilities
#'   are stored.}
#'   \item{A_strings}{Vector of strings naming rows of A_table; the strings 'a' and
#'   'ab', for example, identify the subsets \{a\} and \{a,b\}}
#'   \item{pi_to_P}{Matrix giving choice probabilities when left multiplying preference probabilities}
#'   \item{pi_to_P_logical}{Matrix, same as pi_to_P but with TRUE/FALSE instead of 1/0}
#' }
#'
#' @export
#'
#' @importFrom combinat combn permn
#' @importFrom Matrix invertPerm
#'
#' @examples
#' create_universe(3)      # Example where output is not excessively long
#' u = create_universe(5)  # Example corresponding to datasets in reference, above.
#'
#' @template SMC_reference
#' @template author
create_universe <- function(n) {
  n_objects <- n               # Number of objects in universe
  n_orders <- factorial(n)     # Number of orders (permutations) over objects
  orders <- combinat::permn(n) # List of all orders (permutations)

  # Put permutations in lexicographic order
  order_strings <- sapply(orders, function(x){paste(letters[x], collapse='')})
  p <- order(order_strings)
  order_strings <- order_strings[p]
  orders <- orders[p]
  inv_orders <- lapply(orders, Matrix::invertPerm)  # Compute inverse permutations

  n_subsets <- 2^n-1           # Number of non-empty subsets
  singletons <- 2^(0:(n-1))    # Singleton sets in binary notation
  n_probs <- n*(2^(n-1)-1)     # Number of probabilities in an RCS

  # Initialize tables and matrices
  pi_to_P_logical <- matrix(nrow=n_probs, ncol=n_orders)
  A_table <- matrix(0, nrow=n_subsets, ncol=2)
  colnames(A_table) <- c('nA', 'Ax')
  Ax_table <- matrix(0, nrow=n_probs, ncol=2)
  colnames(Ax_table) <- c('A', 'x')

  # Fill in tables and matrices
  A_strings <- vector('character', n_subsets)
  Ax_strings <- vector('character', n_probs)
  Ax_index <- 1
  for (set_size in 1:n) {
    combin_table <- combinat::combn(n, set_size, simplify=FALSE)
    for (combin_index in 1:length(combin_table)) {
      # A is a subset of {1,...,n} of size set_size
      A_as_vector <- combin_table[[combin_index]]
      A_index <- sum(singletons[A_as_vector])
      A_strings[A_index] <- paste(letters[A_as_vector], collapse='')
      if (set_size > 1) {
        A_table[A_index, 'nA'] <- set_size
        A_table[A_index, 'Ax'] <- Ax_index
        for (x in A_as_vector) {
          # x is an element of A
          x_as_string <- paste(letters[x])
          Ax_table[Ax_index, 'A'] <- A_index
          Ax_table[Ax_index, 'x'] <- x
          Ax_strings[Ax_index] <- paste(A_strings[A_index], x_as_string, sep=';')
          for (order_index in 1:n_orders) {
            inv_order <- inv_orders[[order_index]]
            pi_to_P_logical[Ax_index, order_index] <-
              (inv_order[x]==min(inv_order[A_as_vector]))
          }
          Ax_index <- Ax_index + 1
        }
      }
    }
  }
  pi_to_P <- pi_to_P_logical
  pi_to_P[pi_to_P] <- 1
  rownames(Ax_table) <- Ax_strings
  rownames(A_table) <- A_strings
  dimnames(pi_to_P) <- list(Ax_strings, order_strings)
  dimnames(pi_to_P_logical) <- list(Ax_strings, order_strings)

  list(n=n, n_orders=n_orders, n_subsets=n_subsets, n_probs=n_probs,
       orders=orders, inv_orders=inv_orders, singletons=singletons,
       Ax_table=Ax_table, A_table=A_table,
       Ax_strings=Ax_strings, A_strings=A_strings, order_strings=order_strings,
       pi_to_P=pi_to_P, pi_to_P_logical=pi_to_P_logical)
}

