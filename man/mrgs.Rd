\name{mrgs}
\alias{mrgs}
\alias{mrbn}
\alias{mrps}
\alias{mrnb}
\alias{mrsn}
\title{
Sets the Marginal Component in Marginal Regression Models 
}
\description{
These functions set the marginal component in marginal regression models.

At the moment, the following are implemented:\tabular{ll}{ 
\code{mrgs} \tab normal margins. \cr 
\code{mrbn} \tab binomial margins. \cr
\code{mrps} \tab Poisson margins. \cr 
\code{mrnb} \tab negative binomial margins. \cr
\code{mrsn} \tab skew-normal margins (Azzalini, 1985).
}
}
\usage{
mrgs(link = "identity")
mrbn(link = "logit")
mrps(link = "log")
mrnb(link = "log")
mrsn(link = "identity")
}

\arguments{
  \item{link}{
a specification for the model link function. See \code{\link{family}} for the special case of generalized linear models. 
}
}
\details{
Binomial margins implemented in \code{mrbn} require that y is the form of a two-column matrix with the columns giving the numbers of successes and failures.
Negative binomial margins implemented in \code{mrnb} are parameterized such that \eqn{var(Y)=E(Y)+k E(Y)^2}. 
Skew-normal margins implemented in \code{mrsn} are inherited from package \code{\link{SN}}. 
}
\value{
An object of class \code{\link{marginal.mr}} representing the marginal component. 
}
\references{
Azzalini, A. (1985). A class of distributions which includes the normal ones. \emph{Scandinavian Journal of Statistics} \bold{12}, 171--178.
 
Masarotto, G. and Varin, C. (2010). Gaussian dependence models for non-Gaussian marginal regression. Preprint. 
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{mr}}, \code{\link{dsn}}
}
\examples{
data(polio)
## marginal negative binomial model with ARMA(2,1) errors
m <- mr(y=y.polio,x=x.polio,marginal=mrnb(),vcov=mrarma(2,1),options=list(seed=71271, nrep=300))
m
}
\keyword{regression}
\keyword{nonlinear}
