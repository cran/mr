\name{mrprof}
\alias{mrprof}
\title{
Profile likelihood for Marginal Regression Models
}
\description{
Computes the profile likelihood for a fitted marginal regression model.
}
\usage{
mrprof(m, which, low, up, npoints = 10, display = TRUE, alpha = 0.05)
}
\arguments{
  \item{m}{
a fitted marginal regression model of class \code{\link{mr}}.
}
  \item{which}{
the index of the parameter which should be profiled. 
}
  \item{low}{
the lower limit used in computation of the profile likelihood.
}
  \item{up}{
the upper limit used in computation of the profile likelihood.
}
  \item{npoints}{
number of points used in computation of the profile likelihood. Default is \code{10}.
}
  \item{display}{
should the profile likelihood be displayed or not? default is \code{TRUE}. 
}
  \item{alpha}{
the significance level, default is \code{0.05}. 
}
}
\details{
If display requested, then the profile likelihood is smoothed by cubic spline interpolation.
}
\value{
A list with the following components:
\item{points}{points at which the profile likelihood is evaluated.}
\item{profile}{values of the profile likelihood.}
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{mr}}
}
\examples{
data(polio)
## marginal negative binomial model with ARMA(2,1) errors
m <- mr(y=y.polio,x=x.polio,marginal=mrnb(),vcov=mrarma(2,1),options=list(seed=71271, nrep=100))
prof <- mrprof(m, which=2, low=-9, up=1)
prof
}
\keyword{regression}
\keyword{nonlinear}
