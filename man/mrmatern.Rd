\name{mrmatern}
\alias{mrmatern}
\title{
Matern-type Errors in Marginal Regression Models
}
\description{
Sets a Matern correlation matrix for the errors in marginal regression models. 
}
\usage{
mrmatern(D, k=0.5)
}
\arguments{
  \item{D}{
matrix with values of the distances between pairs of data locations.
}
  \item{k}{
value of the smoothness parameter of the Matern correlation class.
}
}
\value{
An object of class \code{\link{cor.mr}} representing a Matern correlation matrix.
}
\details{
The correlation matrix is inherited from function \code{\link{matern}} of the geoR package.
}
\references{
Masarotto, G. and Varin, C. (2010). Gaussian dependence models for non-Gaussian marginal regression. Preprint.

Ribeiro Jr., P.J. and Diggle, P.J. (2001) geoR: A package for geostatistical analysis. \emph{R-NEWS} \bold{1}, 15--18. 
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{matern}}, \code{\link{mr}}
}
\keyword{regression}
\keyword{nonlinear}