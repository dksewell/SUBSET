% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Proj.R
\name{Proj}
\alias{Proj}
\title{Create a projection matrix}
\usage{
Proj(x)
}
\arguments{
\item{x}{Either a matrix/numeric vector providing the column span of
the projection matrix P, or else a list of such matrices/vecors.  If
the latter, a block-diagonal projection matrix will be constructed.}
}
\value{
A projection matrix.
}
\description{
Create a projection matrix
}
\examples{
Proj(rep(1,3))
Proj(cbind(1,1:3))
Proj(list(diag(2),cbind(1,1:3),rep(1,4)))
}
