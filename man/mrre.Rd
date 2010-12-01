\name{mrre}
\alias{mrre}
\title{
Random-Effects Type Errors in Marginal Regression Models
}
\description{
Sets "random-effects type" errors in marginal regression models.
}
\usage{
mrre(z)
}
\arguments{
\item{z}{
	random effects design matrix. 	
}
}
\value{
An object of class \code{\link{cor.mr}} representing a correlation matrix with "random-effects type" structure. See the "Salamander" data example in function \code{\link{mr}}.
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