% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SUBSET_gibbs.R
\name{SUBSET_gibbs}
\alias{SUBSET_gibbs}
\title{SUBspace Shrinkage via Exponential Tilting for estimated
projection matrix}
\usage{
SUBSET_gibbs(
  Mean0,
  Sigma0,
  n_draws,
  pi0_samples,
  P_phi = c("power", "geometric")[1],
  prior_phi,
  v_sequence = seq(0.25, 5, by = 0.25),
  Z_approximation = list(df = 10 + 1, phi_range = NULL),
  verbose = TRUE,
  cl
)
}
\arguments{
\item{Mean0}{vector. The posterior mean corresponding
to the base prior (i.e., untilted prior).}

\item{Sigma0}{matrix. The posterior covariance
matrix corresponding to the base prior. NOTE: The
variable names are extracted from Sigma0, so use
colnames(Sigma0) = ...}

\item{n_draws}{positive integer of how many posterior
samples are desired.}

\item{pi0_samples}{Samples of the model parameters
from the \emph{prior} (not including phi)}

\item{P_phi}{Either (1) a function taking in a scalar and returning
a projection matrix for the model parameters; or (2) "power"
in which case the projection matrix will correspond to
span((1,1,...)',(1^x,2^x,...)'); or (3) "geometric" in
which case the projection matrix will corrspond to
span((1,1,...)',(x^1,x^2,...)').}

\item{prior_phi}{list with named argument "r".  prior_phi$r(n)
should provide n random draws from the prior on phi.}

\item{v_sequence}{vector of positive values. See details.}

\item{Z_approximation}{If set to be NULL, Z will be
evaluated "exactly" using the pi0_samples.  Otherwise a list
must be provided in order to approximate Z(phi) via
natural cubic splines.  If the latter, the list must
have named arguments df, the degrees of freedom (see
ns()), and phi_range, a vector of length 2 giving the
lower and upper bound of \eqn{\phi}.  phi_range is
optional, and if left missing will be obtained from
the range of prior_phi(1000).}

\item{verbose}{logical. Should any output be provided?}

\item{cl}{(optional) cluster from the parallel package.}
}
\value{
Object of class "subset_gibbs".
\itemize{
\item samples: 3-dim array.  First index is the samples,
second index is the variable, and third index is the
value of v (the exponential tilting parameter)
\item acc_rate: numeric giving the acceptance rate for phi
\item Mean0: the posterior mean without tilting
\item Sigma0: the posterior covariance without tilting
\item pi0_samples: Samples of the model parameters
from the \emph{prior} (not including phi)
\item P_phi: a function taking in a scalar and returning
a projection matrix for the model parameters
\item prior_phi: list with named argument "r".
\item v_sequence: vector of exponential tilting parameters.
}
}
\description{
SUBSET_gibbs() will compute perform a 2-block Gibbs sampler,
where the blocks correspond to the primary model parameters
and the projection matrix parameter (only one parameter
supported at this time).  The Gibbs algorithm relies on
the asymptotic normal approximation to the conditional posterior
using an exponential tilted prior that shrinks the
parameters towards a linear subspace for a given projection
matrix.
}
\details{
SUBSET takes a base prior \eqn{\pi_0} and tilts it towards
a linear subspace via the multiplicative term
\deqn{e^{-\frac{v}{2}\theta'(I-P)\theta}.}
The user must perform the initial Bayesian inference for
this base prior, and compute the posterior mean and covariance.
\emph{It is imperative that the necessary conditions for an asymptotic
normal approximation of the posterior hold for the base prior.}
For a projection matrix \eqn{P} parameterized by \eqn{\phi},
SUBSET_gibbs() will obtain joint draws of the model parameters and
\eqn{\phi}. \code{v_sequence} is a user provided positive
number or vector of positive numbers that control the shrinkage:
higher values imply more shrinkage of the posterior while \eqn{v=0}
corresponds to the base prior with no shrinkage.
}
