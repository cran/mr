\name{estfun.mr}
\alias{estfun.mr}
\alias{bread.mr}
\title{
Bread and Salami for Sandwich in Marginal Regression Models
}
\description{
Extracts empirical estimating functions (aka salami) and the inverse of (minus) the Hessian (aka bread) of a fitted marginal regression model for computation of the sandwich variance-covariance matrix (Zeleis, 2006). 
}
\usage{
estfun.mr(x, ...)
bread.mr(x, ...)
}
\arguments{

  \item{x}{
a fitted marginal regression model of class \code{\link{mr}}.
}
  \item{\dots}{
further arguments passed to methods.
}
}
\details{
These functions are passed to \code{\link{sandwich}} (Zeleis, 2006).
}
\value{
\code{estfun.mr}: a matrix containing the empirical estimating functions. 

\code{bread.mr}: a matrix containing the inverse of (minus) the Hessian. 
}
\references{
Zeleis, A. (2006). Object-oriented computation of sandwich estimators. \emph{Journal of Statistical Software} \bold{16}, issue 9.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{bread}}, \code{\link{estfun}}, \code{\link{mr}}, \code{\link{sandwich}}
}
\keyword{regression}
\keyword{nonlinear}
