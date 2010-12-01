\name{mr}
\alias{mr}
\alias{mr-class}
\title{
Fitting Marginal Regression Models by Maximum (Simulated) Likelihood
}
\description{
Fits marginal regression models described in Masarotto and Varin (2010) by maximum (simulated) likelihood.
}
\usage{
mr(x = rep(1, NROW(y)), y, marginal, vcov, start, fixed, options = mroptions())
}
\arguments{
  \item{x}{
model matrix.
}
  \item{y}{
response vector. 
}
\item{marginal}{
an object of class \code{\link{marginal.mr}} specifying the marginal part of the model. 
}
  \item{vcov}{
an object of class \code{\link{cor.mr}} representing the correlation matrix of the errors. 
}
  \item{start}{
optional numeric vector with starting values for the model parameters.
}
  \item{fixed}{
optional numeric vector of the same length as the total number of parameters. If supplied, only \code{NA} entries in fixed will be varied. 
}
\item{options}{
list of options passed to function \code{\link{mroptions}}.
}
}
\details{
Computation of the exact likelihood is possible only for continuous responses, otherwise the likelihood function is approximated by importance sampling. See Masarotto and Varin (2010) for details. Functions \code{\link{coefficients}}, \code{\link{logLik}}, \code{\link{mrparvar}}, \code{\link{mrse}} and \code{\link{residuals.mr}} can be used to extract various useful features of the value returned by \code{\link{mr}}.
}
\value{
An object of class \code{"mr"} with the following components:
\item{estimate}{the vector of parameter estimates.}
\item{maximum}{the maximum (simulated) likelihood.}
\item{vcov}{the variance-covariance matrix of parameter estimates.}
\item{hessian}{(minus) the Hessian at the maximum likelihood estimate.}
\item{jac}{the Jacobian at the maximum likelihood estimate.}
\item{y}{the y vector used.}
\item{x}{the model matrix used.}
\item{n}{number of observations.}
\item{call}{the matched call.}
\item{not.na}{vector of binary indicators of missing observations.}
\item{fixed}{numeric vector indicating which parameters are constants.}
\item{ibeta}{indices of marginal parameters.}
\item{igamma}{indices of dependence parameters.}
\item{options}{fitting options used, see \code{\link{mroptions}}.}
}
\references{
Masarotto, G. and Varin, C. (2010). Gaussian dependence models for non-Gaussian marginal regression. Preprint. 
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{cor.mr}}, \code{\link{marginal.mr}}, \code{\link{mroptions}}
}
\examples{
## Warning: Likelihood is approximated using only a limited number of Monte Carlo replications.
## Polio data. Marginal negative binomial model with ARMA(2,1) errors.
data(polio)
mr(x=x.polio,y=y.polio,marginal=mrnb(),vcov=mrarma(2,1),options=list(seed=71271, nrep=100))
## Salamanders data. Analysis of the first experiment.
data(salamanders)
attach(salamanders)
mr(x=cbind(1,x[,1]+x[,3],x[,1]+x[,2],x[,1]),y=cbind(y1,1-y1),marginal=mrbn(),vcov=mrre(list(zf,zm)),
   options=list(seed=71271,nrep=1000))
detach(salamanders) 
}
\keyword{regression}
\keyword{nonlinear}
