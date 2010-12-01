\name{residuals.mr}
\alias{residuals.mr}
\title{Quantile Residuals from Fitted Marginal Regression Models}
\usage{
\method{residuals}{mr}(object, type=c("conditional","marginal"),
                       method=c("random","interval","mid"),...)
}
\arguments{
  \item{object}{an object of class \code{mr}, typically the result of
    a call to \code{\link{mr}}.}
  \item{type }{the type of quantile residuals which should be returned.
    The alternatives are: \code{"conditional"} (default) and \code{"marginal"}.}
  \item{method }{different methods available for quantile residuals in case of noncontinuous responses:
    \code{"random"} for randomized quantile residuals (default), 
    \code{"interval"} for interval quantile residuals as defined in Zucchini and MacDonald (2009) and 
    \code{"mid"} for mid interval qantile residuals as defined in Zucchini and MacDonald (2009).}
\item{...}{further arguments passed to or from other methods.}
}
\description{
  Computes various type of quantile residuals for validation of a fitted marginal regression model described by Masarotto and Varin (2010).
}
\details{
 Quantile residuals are defined in Dunn and Smyth (1996). Two different types are available:
\tabular{ll}{
\code{conditional} \tab quantile residuals that account for the dependence. \cr
\code{marginal} \tab quantile residuals that do not account for the dependence.
}
Conditional quantile residuals are normal quantiles of Rosenblatt (1952) transformations and they are appropriate for validation of the marginal regression models discussed by Masarotto and Varin (2010). If the responses are noncontinuous, then the conditional quantile residuals are not well defined. This difficulty is overcame by randomized quantile residuals available through option \code{method="random"}. Alternatively, Zucchini and MacDonald (2009) suggest the use of interval quantile residuals (\code{method="interval"}) and mid interval quantile residuals (\code{method="mid"}). Note that randomized quantile residuals are realizations of independent standard normal variables under model conditions, while interval and mid quantile residuals are not.
}
\seealso{
  \code{\link{mr}}
}
\references{
Dunn, P.K. and Smyth, G.K. (1996). Randomized quantile residuals. \emph{Journal of Computational and Graphical Statistics} \bold{5}, 236--244.

Masarotto, G. and Varin, C. (2010). Gaussian dependence models for non-Gaussian marginal regression. Preprint.

Rosenblatt, M. (1952). Remarks on a multivariate transformation. \emph{The Annals of Mathematical Statistics} \bold{23}, 470--472.

Zucchini, W. and MacDonald, I.L. (2009). \emph{Hidden Markov Models for Time Series}. Chapman and Hall/CRC.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\note{
It is appropriate to inspect several sets of randomized quantile residuals before to take a decision about the model. See Masarotto and Varin (2010) for more details.
}
\examples{
data(polio)
## marginal negative binomial model with ARMA(2,1) errors
m <- mr(y=y.polio,x=x.polio,marginal=mrnb(),vcov=mrarma(2,1),options=list(seed=71271, nrep=100))
## randomized quantile residuals
res <- residuals(m)
## normal probability plot
qqnorm(res)
qqline(res)
## autocorrelation plot
acf(res)
}
\keyword{models}
\keyword{regression}

