# vcov must have class vcov.mr and the following elements:
# - npar: number of parameters. Wow. Really?
# - start(marginal, betastart): computes initial estimates
# - chol(gamma,not.na): returns the cholesky factor of the vcov matrix (only for the not.na obs)
#                       returns NULL if gamma is outside parameter space

# independent
mrind <- function() {
    ans <- list()
    ans$npar <- 0
    ans$start <- function() double(0)
    ans$chol <- function(gamma , not.na) diag(rep(1,sum(not.na)))
    ans$independent <- TRUE
    class( ans ) <- "cor.mr"
    ans
}

# arma
mrarma<- function( p , q ) {
    start <- rep( 0 , p+q )
    names(start) <- c(if ( p ) paste("ar",1:p,sep="") else NULL ,
                      if ( q ) paste("ma",1:q,sep="") else NULL )
    iar <- if ( p ) 1:p else NULL
    ima <- if ( q ) (p+1):(p+q) else NULL
    ans <- list()
    ans$npar <- length(start)
    ans$start <- function() start
    ans$chol <- function( gamma , not.na ) {
        if ( ( p && any(Mod(polyroot(c(1,-gamma[iar])))<1.01) ) ||
             ( q && any(Mod(polyroot(c(1, gamma[ima])))<1.01) ) )
            return( NULL )
        n <- length(not.na)
        rho <- ARMAacf(gamma[iar],gamma[ima],n-1)
        r <- seq(1,n)[not.na]
        chol(outer( r , r , function(i,j) rho[1+abs(i-j)] ))
    }
    class( ans ) <- "cor.mr"
    ans
}

## assume that it is not possible that all the observations inside a cluster
## can be missing
mrcluster <- function(id, type=c("ar1", "ma1", "exch", "unstr")) { 
    type <- match.arg(type)
    ng <- 1:length(unique(id))
    if (!(length(ng)>1)) stop("only one strata")
    ans <- list(type=type,id=id)
    if(type=="unstr")
        ans$npar <- choose(max(table(id)), 2)
    else
        ans$npar <- 1
    start <- rep(0, ans$npar)
    names(start) <- paste("gamma", 1:ans$npar, sep="")
    data <- data.frame(id=id)
    fn <- switch(type, "ar1"=function(g) nlme::corAR1(g, form= ~1|id),
                 "ma1"=function(g) nlme::corARMA(g, form= ~1|id, p=0, q=1),
                 "exch"=function(g) nlme::corCompSymm(g, form= ~1|id),
                 "unstr"=function(g) nlme::corSymm(g, form= ~1|id))
    ans$start <- function() start
    ans$chol <- function(gamma, not.na) {
        q <- try(nlme::corMatrix(nlme::Initialize(fn(gamma),data=data)),silent=TRUE)
        if (inherits(q,"try-error")) return(NULL)
        g <- split(not.na,id)
        q <- try(lapply(ng,function(i) chol(q[[i]][g[[i]],g[[i]]])),silent=TRUE)
        if (inherits(q,"try-error") ) NULL else q
    }
    class( ans ) <- "cor.mr"
    ans
}


## V = I + sum_i (gamma_i R_i) and D = diag(diag(V))
## if (gammas.are.cors) omega=I+V-D
## otherwise omega = D^(-1/2)(I+V)D^(-1/2) 
mrvsum<- function(R, gammas.are.cors=TRUE) {
    if (!is.list(R)) stop("first argument must be a list")
    n <- NROW(R[[1]])
    if (any(sapply(R,function(x) dim(x)!=n))) stop("Matrices with different dimension")
    R <- matrix(sapply(R,function(x) if(gammas.are.cors) x-diag(diag(x)) else x),n*n)
    O <- as.double(diag(n))
    eps <- sqrt(.Machine$double.eps)
    start <- rep( 0 , NCOL(R) )
    names(start) <- paste("gamma", 1:length(start), sep="")
    ans <- list()
    ans$npar <- length(start)
    ans$start <- function() start
    ans$chol <- function( gamma , not.na ) {
        g <- matrix(rowSums(sweep(R,2,gamma,"*"))+O,n)
        if (!gammas.are.cors) {
            d <- diag(g)
            if (any(d<=eps)) return(NULL)
            d <- sqrt(d)
            g <- g / outer(d,d)
        }
        g <- try(chol(g[not.na,not.na]),silent=TRUE)
        if (inherits(g,"try-error") ) NULL else g
    }
    class( ans ) <- "cor.mr"
    ans
}

mrre <- function(z) {
    z.good <- FALSE
    if (is.matrix(z)) {
        z.good <- TRUE
        R <- list(z %*% t(z))
    } else if (is.list(z)) {
        if (all(sapply(z,is.matrix)) || all(sapply(z,function(zi) NROW(zi)==NROW(z[[1]])))) {
            z.good <- TRUE
            R <- lapply(z,function(zi) zi %*% t(zi) )
        }
    }
    if (!z.good ) stop("first argument must be a matrix or a list of matrices with the same number of rows")
    mrvsum(R)
}

## added by Sammy 2010-11-22
mrmatern <- function(D, k=0.5){
  ans <- list()
  ans$npar <- 1
  start <-median(D) 
  names(start) <- c("gamma")
  ans$start <- function( marginal, beta ) start
  ans$chol <- function( gamma, not.na ){
    if( any(gamma<=0) ) return( NULL )
    S <- geoR:::matern(D, gamma, k)
    q <- try(chol(S),silent=TRUE)
    if( inherits(q,"try-error") ) NULL else q
  }
  class( ans ) <- "cor.mr"
  ans
}
