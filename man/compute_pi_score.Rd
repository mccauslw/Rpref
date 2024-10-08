% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pi_maxl.R
\name{compute_pi_score}
\alias{compute_pi_score}
\title{Evaluate gradient of random preference (RP) log likelihood}
\usage{
compute_pi_score(pi, u, N)
}
\arguments{
\item{pi}{Vector of preference probabilities of length n!}

\item{u}{List with precomputed information about choice universes
of size n, created using \code{\link[=create_universe]{create_universe()}}}

\item{N}{Matrix containing observed choice counts in a choice experiment.}
}
\value{
Gradient of RP log likelihood, a vector of size n!.
}
\description{
Evaluate gradient, with repect to pi, of RP log likelihood,
as a function of pi, for given data N.
}
\examples{
library(RanCh)
n <- 5
n_fact = factorial(n)
pi <- rep(1/n_fact, n_fact) # Uniform distribution over preference orders
u <- create_universe(n)
N <- vectorize_counts(u, RanCh::MMS_2019_counts[1, , ])
pi_score <- compute_pi_score(pi, u, N)

}
\references{
McCausland, W. (2024). Sequential Monte Carlo for Random Prefernces. Unpublished manuscript.
}
\author{
William McCausland, \email{william.j.mccausland@umontreal.ca}
}
