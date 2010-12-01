## Marginals

## Marginal must have class marginal.mr and the following elements:
## - start(y,x): compute initial estimates ignoring correlations
## - dp(y,x,beta): evaluate [d,p]
## - q(u,y,x,beta): evaluate quantiles
## - npar(x): number of parameters
## - type: response type (integer or numeric)

# Gaussian
mrgs <- function(link = "identity" ) {
    fm <- gaussian( substitute( link ) )
    ans <- list()
    ans$start <- function(y, x) {
        m <- glm.fit( x , y, family=fm )
        ans <- c( coef( m ) , log(sqrt( mean( residuals( m )^2 ) ) ))
        names(ans) <- c(paste("beta",1:length(coef(m)),sep=""),"log.sigma")
        ans
    }
    ans$npar <- function(x) NCOL(x)+1
    ans$dp <- function(y, x, beta) {
        nb <- length(beta)
        mu <- fm$linkinv( x %*% beta[-nb] )
        sd <- exp( beta[nb] )
        cbind( dnorm( y , mu , sd ) , pnorm( y , mu , sd ) )
    }
    ans$sim <- function( u, y, x, beta) {
        nb <- length(beta)
        mu <- fm$linkinv( x %*% beta[-nb] )
        sd <- exp( beta[nb] )
        qnorm( u , mu , sd )
    }
    ans$type <- "numeric"
    class(ans) <- c( "marginal.mr")
    ans
}

# Binomial
# y is [#success,#failure]
mrbn <- function(link = "logit") {
    fm <- binomial( substitute( link ) )
    ans <- list()
    ans$start <- function(y, x) {
        res <- coef(glm.fit( x , y, family=fm ))
        names(res) <- paste("beta",1:length(res),sep="")
        res
    }
    ans$npar <- function(x) NCOL(x)
    ans$dp <- function(y, x, beta) {
        mu <- fm$linkinv( x %*% beta )
        cbind(dbinom( y[,1] , y[,1]+y[,2], mu ) ,
              pbinom( y[,1] , y[,1]+y[,2] , mu ) )
    }
    ans$sim <- function(u, y, x, beta) {
        mu <- fm$linkinv( x %*% beta )
        cbind(qbinom( u, y[,1]+y[,2] , mu ), y[,2])
    }
    ans$type <- "integer"
    class(ans) <- c( "marginal.mr")
    ans
}


# Poisson
mrps <- function(link = "log") {
    fm <- poisson( substitute( link ) )
    ans <- list()
    ans$start <- function(y, x) {
        res <- coef(glm.fit( x , y, family=fm ))
        names(res) <- paste("beta",1:length(res),sep="")
        res
    }
    ans$npar <- function(x) NCOL(x)
    ans$dp <- function(y, x, beta) {
        mu <- fm$linkinv( x %*% beta )
        cbind( dpois( y , mu ) , ppois( y , mu ) )
    }
    ans$sim <- function(u, y, x, beta) {
        mu <- fm$linkinv( x %*% beta )
        qpois( u , mu )
    }
    ans$type <- "integer"
    class(ans) <- c( "marginal.mr")
    ans
}


# Negative binomial
# var(y) = E(y) + k*E(y)^2 (k>0)
mrnb <- function(link = "log" ) {
    fm <- poisson( substitute( link ) )
    h <- 2
    h2 <- 2 - h
    eps <- .Machine$double.eps
    ans <- list()
    ans$start <- function(y, x) {
        m <- glm.fit( x , y, family=fm )
        mu <- fitted(m)
        res <- c( coef( m ) , max(0.01, mean(((y-mu)^2-mu)/mu^h)))
        names(res) <- c(paste("beta",1:length(coef(m)),sep=""),"k")
        res
    }
    ans$npar <- function(x) NCOL(x)+1
    ans$dp <- function(y, x, beta) {
        nb <- length(beta)
        if (beta[nb]<=0) rep(NA,NROW(y))
        mu <- fm$linkinv( x %*% beta[-nb] )
        size <- pmax(eps, mu^h2 / beta[nb])
        cbind( dnbinom( y , mu=mu , size=size) , pnbinom( y , mu=mu , size=size) )
    }
    ans$sim <- function( u, y, x, beta) {
        nb <- length(beta)
        if (beta[nb]<=0) rep(NA,NROW(y))
        mu <- fm$linkinv( x %*% beta[-nb] )
        size <- mu^h2 / beta[nb]
        qnbinom( u , mu=mu, size=size)
    }
    ans$type <- "integer"
    class(ans) <- c( "marginal.mr")
    ans
}



## marginale skew-normal
mrsn <- function( link="identity" ) {
    if ( !require(sn) ) stop("This requires package sn")  
    fm <- gaussian( substitute( link ) ) ##  ;-)
    eps <- sqrt(.Machine$double.eps)
    ans <- list()
    ans$start <- function(y, x) {
        m <- sn.mle( x, y , plot.it=FALSE)
        ## "centered" sn parameterization with further transformation
        ## to keep unconstrained optimization
        ## Starting point is taken from sn.ml
        qrX <- qr(x)
        r <- qr.resid(qrX,y)
        s <- sqrt(sum(qr.resid(qrX, y)^2)/length(y))
        gamma1 <- sum(qr.resid(qrX, y)^3)/(length(y) * s^3)
        if (abs(gamma1) > 0.99 ) gamma1 <- sign(gamma1) * 0.99
        ans <- c(qr.coef(qrX, y), log(s), log((gamma1+0.995)/(0.995-gamma1)) )
        names(ans) <- c(paste("beta",1:NCOL(x),sep=""),"log.scale", "logit.shape")
        ans
    }
    ans$npar <- function(x) NCOL(x)+2
    ans$dp <- function(y, x, beta ) {
        p <- NCOL(x)
        n <- length(y)
        ## back transformation
        beta[p+1] <- exp(beta[p+1])
        beta[p+2] <- exp(beta[p+2])
        beta[p+2] <- 0.995*(beta[p+2]-1)/(beta[p+2]+1)
        ## now go back to "direct" parameterization
        beta <- cp.to.dp( beta ) 
        mu <- fm$linkinv( x %*% beta[1:p] )
        scale <- beta[p+1]
        alpha <- beta[p+2] 
        cbind(dsn( y, mu, scale, alpha ) ,
              psn( y, mu, scale, alpha) )
    }
    ans$sim <- function( u, y, x, beta){
        p <- NCOL(x)
        ## back transformation
        beta[p+1] <- exp(beta[p+1])
        beta[p+2] <- exp(beta[p+2])
        beta[p+2] <- 0.995*(beta[p+2]-1)/(beta[p+2]+1)
        ## now go back to "direct" parameterization
        beta <- cp.to.dp( beta ) 
        mu <- fm$linkinv( x %*% beta[1:p] )
        scale <- beta[p+1]
        alpha <- beta[p+2] 
        qsn( u, mu, scale, alpha )
    }
    ans$type <- "numeric"
    class(ans) <- "marginal.mr"
    ans
}
