\name{mrparvar}
\alias{mrparvar}
\title{
Variance-Covariance Matrix for Fitted Marginal Regression Models
}
\description{
Returns the variance-covariance matrix of the parameters of a fitted marginal regression model.
}
\usage{
mrparvar(x, type = c("hessian", "sandwich", "vscore", "cluster", "hac"), ...)
}
\arguments{
  \item{x}{
a fitted marginal regression model of class \code{\link{mr}}.
}
  \item{type}{
type of variance-covariance matrix. This can be one of the following:
\tabular{ll}{
\code{hessian} \tab inverse of the observed Fisher information (default). \cr
\code{sandwich}\tab sandwich variance matrix. \cr
\code{vscore} \tab cross-product of scores. \cr
\code{cluster} \tab sandwich variance matrix for cluster/longitudinal data. \cr
\code{hac} \tab heteroskedasticity and autocorrelation consistent estimator for time series of Andrews (1991).
}
}
  \item{\dots}{
further arguments passed to function \code{\link{vcovHAC}}.
}
}
\details{
\code{sandwich} and \code{vscore} are based on the predictive decomposition of the joint density. \code{cluster} uses the decomposition of the data in independent subject-specific blocks. \code{hac} is appropriate for time series and uses the implementation in package \code{\link{sandwich}} (Zeileis, 2006).
}
\value{
A matrix containing the estimated covariances between the parameter estimates of a fitted \code{\link{mr}} model.
}
\references{
Andrews, D.W.K. (1991). Heteroskedasticity and autocorrelation consistent covariance matrix estimation. \emph{Econometrica}, \bold{59}, 817--858.

Masarotto, G. and Varin. C. (2010). Gaussian dependence models for non-Gaussian marginal regression. Preprint. 

Zeleis, A. (2006). Object-oriented computation of sandwich estimators. \emph{Journal of Statistical Software} \bold{16}, issue 9.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{mr}}, \code{\link{sandwich}}, \code{\link{vcovHAC}}
}
\examples{
data(polio)
## working independence       
m <- mr(y=y.polio,x=x.polio,marginal=mrnb(),vcov=mrind())
m
## HAC variance-covariance matrix
round(mrparvar(m, "hac"), 2)
}
\keyword{regression}
\keyword{nonliner}
