\name{mrarma}
\alias{mrarma}
\title{
ARMA(p,q) Errors in Marginal Regression Models
}
\description{
Sets ARMA(p,q) errors in marginal regression models.
}
\usage{
mrarma(p, q)
}
\arguments{
  \item{p}{
order of the autoregressive component.
}
  \item{q}{
order of the moving average component.
}
}
\value{
An object of class \code{\link{cor.mr}} representing a correlation matrix with ARMA(p,q) structure.
}
\references{
Masarotto, G. and Varin, C. (2010). Gaussian dependence models for non-Gaussian marginal regression. Preprint.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{mr}}
}
\keyword{regression}
\keyword{nonlinear}
