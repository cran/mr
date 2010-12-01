\name{mrse}
\alias{mrse}
\title{
Variance-Covariance Matrix for Fitted Marginal Regression Models
}
\description{
Returns the standard errors vector of a fitted marginal regression model.
}
\usage{
mrse(x, ...)
}
\arguments{
  \item{x}{
a fitted marginal regression model of class \code{\link{mr}}.
}
  \item{\dots}{
further arguments passed to function \code{\link{mrparvar}}.
}
}
\value{
A vector containing the estimated standard errors of the parameter of a fitted \code{\link{mr}} model.
}
\references{
Masarotto, G. and Varin. C. (2010). Gaussian dependence models for non-Gaussian marginal regression. Preprint. 
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{mr}}, \code{\link{mrparvar}}
}
\examples{
data(polio)
## ARMA(2,1)      
m <- mr(y=y.polio,x=x.polio,marginal=mrnb(),vcov=mrarma(2,1),options=list(seed=71271, nrep=100))
m
## various type of standard errors
cbind(mrse(m,"he"),mrse(m,"v"),mrse(m,"s"),mrse(m,"hac"))
}
\keyword{regression}
\keyword{nonliner}
