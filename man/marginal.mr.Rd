\name{marginal.mr}
\title{Marginal Component in Marginal Regression Models}
\alias{marginal.mr-class}
\alias{marginal.mr}
\description{
  Class of marginal components available in the \code{mr} library. 
}
\value{
At the moment, the following are implemented:
\tabular{ll}{ 
\code{\link{mrgs}} \tab normal margins. \cr 
\code{\link{mrbn}} \tab binomial margins. \cr
\code{\link{mrps}} \tab Poisson margins. \cr 
\code{\link{mrnb}} \tab negative binomial margins. \cr
\code{\link{mrsn}} \tab skew-normal margins.  
}
}
\references{
 Masarotto, G. and Varin, C. (2010). Gaussian dependence models for non-Gaussian marginal regression. Preprint.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{\code{\link{mr}},
  \code{\link{mrgs}},
  \code{\link{mrbn}},
  \code{\link{mrps}}, 
  \code{\link{mrnb}},
  \code{\link{mrsn}}
}
\keyword{regression}
\keyword{nonlinear}
