\name{cor.mr}
\title{Dependence Component in Marginal Regression Models}
\alias{cor.mr-class}
\alias{cor.mr}
\description{
  Class of correlation matrix of the errors available in the \code{mr} library. 
}
\value{
At the moment, the following are implemented:
\tabular{ll}{ 
\code{\link{mrind}} \tab working independence. \cr
\code{\link{mrarma}} \tab ARMA(p,q).  \cr 
\code{\link{mrcluster}} \tab longitudinal/clustered data.  \cr
\code{\link{mrmatern}} \tab  Matern spatial correlation. \cr
\code{\link{mrre}} \tab random-effects type. \cr 
}
}
\references{
 Masarotto, G. and Varin, C. (2010). Gaussian dependence models for non-Gaussian marginal regression. Preprint.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{\code{\link{mr}},
  \code{\link{mrind}},
  \code{\link{mrarma}},
  \code{\link{mrcluster}}, 
  \code{\link{mrmatern}},
  \code{\link{mrre}}
}
\keyword{regression}
\keyword{nonlinear}
