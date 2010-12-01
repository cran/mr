\name{mroptions}
\alias{mroptions}
\title{
Setting Options for Fitting Marginal Regression Models
}
\description{
Sets options which affect fitting of marginal regression models.
}
\usage{
mroptions(seed = round(runif(1, 1, 1e+05)), nrep = c(100, 1000), no.se = FALSE, opt = mropt)
}
\arguments{
  \item{seed}{
seed of the pseudorandom generator used in the importance sampling algorithm  for likelihood approximation in case of noncontinuous responses.
}
  \item{nrep}{
Monte Carlo size of the importance sampling algorithm  for likelihood approximation in case of noncontinuous responses.
}
  \item{no.se}{
logical. Should standard errors be computed and returned or not?
}
  \item{opt}{
    the function used to optimize the likelihood. The default uses
    \code{nlminb}. See \code{mr:::mropt}.
}
}
\value{
A list containing the options.
}
\references{
Masarotto, G. and Varin. C. (2010). Gaussian dependence models for non-Gaussian marginal regression. Preprint.
}
\author{
Guido Masarotto and Cristiano Varin. 
}
\seealso{
\code{\link{mr}}
}
\keyword{regression}
\keyword{nonlinear}
