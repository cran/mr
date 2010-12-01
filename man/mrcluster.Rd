\name{mrcluster}
\alias{mrcluster}
\title{
Longitudinal and Clustered Data Errors in Marginal Regression Models
}
\description{
Sets longitudinal/clustered data errors in marginal regression models.
}
\usage{
mrcluster(id, type = c("ar1", "ma1", "exch", "unstr"))
}
\arguments{
  \item{id}{
subject id.
}
  \item{type}{
a character string specifying the correlation structure. At the moment, the following are implemented: 
\tabular{ll}{
\code{ar1} \tab autoregressive of order 1. \cr
\code{ma1} \tab moving average of order 1. \cr
\code{exch} \tab echangeable. \cr
\code{unstr} \tab unstructured.
}
}
}
\value{
An object of class \code{\link{cor.mr}} representing a correlation matrix for longitudinal or clustered data.
}
\details{
The correlation matrices are inherited from the \code{\link{nlme}} package.
}
\references{
Masarotto, G. and Varin, C. (2010). Gaussian dependence models for non-Gaussian marginal regression. Preprint.

Pinheiro, J.C., and Bates, D.M. (2000). \emph{Mixed-Effects Models in S and S-PLUS}. Springer.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{mr}}, \code{\link{nlme}}
}
\keyword{regression}
\keyword{nonlinear}