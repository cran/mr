# Return a function which computes an approximation of log(c(p(y[1]),p(y[2]|y[1]),...))
mrloglik <- function(m, irep, magic=NA) {
    n <- m$n
    not.na <- m$not.na
    y <- m$y[not.na,]
    x <- m$x[not.na,,drop=FALSE]
    is.int <- m$marginal$type == "integer"
    ind.lik <- !is.null(m$vcov$independent) 
    dp <- m$marginal$dp
    ibeta <- m$ibeta
    igamma <- m$igamma
    vchol <- m$vcov$chol
    if (missing(irep)) irep <- max(m$options$nrep)
    seed <- m$options$seed
    id <- if(is.null(m$vcov$id)) rep(1,n) else m$vcov$id
    lstrata <- rle(id)$length
    ipar <- c(n,irep,length(lstrata),lstrata)
    Cfun <- if (is.int) "csampler" else "bsolver"
    magic <- rep(magic,m$n)
    theta <- m$fixed
    ifree <- is.na(theta)
    cache <- new.env()
    function( theta.free ) {
        theta[ifree] <- theta.free
        beta <- theta[ibeta]
        if (!identical(cache$beta,beta)) {
            dp <- dp(y,x,beta)
            assign("beta",beta,envir=cache)
            assign("dp",dp,envir=cache)
        } else {
            dp <- get("dp",envir=cache)
        }
        if ( !all(is.finite(dp)) ) return( magic )
        if ( ind.lik ) return( log(dp[,1]) )
        gamma <- theta[igamma]
        if (!identical(cache$gamma,gamma)) {
            q <- vchol( gamma , not.na )
            assign("gamma",gamma,envir=cache)
            assign("q",q,envir=cache)
        } else {
            q <- get("q",envir=cache)
        }
        if ( is.null(q) ) return ( magic )
        if (is.int) set.seed(seed)
        lk <- .C(Cfun, as.integer(ipar), as.double(unlist(q)), as.double(dp),
                 lk=double(n), NAOK=TRUE,dup=FALSE,package="mr")$lk
        if ( all(is.finite(lk)) ) lk else magic
    }    
} 



# Workhorse which is called from  mr  and mrprof
mrtruefit <- function(x) {
    ## saving/restoring the random seed
    if ( x$marginal$type == "integer" ) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            seed.keep <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
            on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
        }
    }
    ifree <- is.na( x$fixed )
    theta <- x$fixed
    theta[ifree] <- x$estimate[ifree]
    big <- -sqrt(.Machine$double.xmax)
    for ( irep in x$options$nrep ) {
        log.lik <- mrloglik(x,irep,big)
        ans <- x$options$opt( theta[ifree] , log.lik )
        theta[ifree] <- ans$estimate
    }
    names(theta) <- names(x$estimate)
    x$estimate <- theta
    x$maximum <- ans$maximum
    x
}

# Add/replace an estimate of jacobian and hessian to an mr object
# Hessian is approximated by finite differences
# in a rotated space in which the true hessian is near
# the identity matrix. The rotation is based on
# jacobian crossprod
mrjhess <- function(x, options=x$options, only.jac = FALSE) {
    if ( !inherits( x , "mr" ) ) stop("First argument must be a mr object")
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        seed.keep <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
    }
    theta <- x$estimate
    ifree <- is.na( x$fixed )
    theta.free <- theta[ifree]
    log.lik <- mrloglik(x)
    eps <- .Machine$double.eps^(1/4)
    relStep <- 0.1
    maxtry <- 10
    delta <- ifelse(abs(theta.free)<1, eps, eps*theta.free)
    di <- function(i,delta) {
        x1 <- x2 <- theta.free
        x1[i] <- x1[i] - delta[i]
        x2[i] <- x2[i] + delta[i]
        (log.lik(x2)-log.lik(x1))/(2*delta[i])
    }
    while (1) {
        x$jac <- sapply(seq_along(theta.free),di,delta)
        if( all(is.finite(x$jac)) ) break
        delta <- delta/2
        maxtry <- maxtry - 1
        if (maxtry<0) stop("impossible to compute a finite jacobian")
    }
    if (!only.jac) {
        a <- svd(x$jac)
        a$d <- pmax(a$d,sqrt(.Machine$double.eps)*a$d[1])
        x$hessian <- nlme:::fdHess(rep(0,length(theta.free)),
                                   function(tx) sum(log.lik(theta.free+a$v%*%(tx/a$d))),
                                   minAbsPar=1,.relStep=relStep)$Hessian
        x$hessian = (x$hessian+t(x$hessian))/2
        x$hessian <- a$v%*%(outer(a$d,a$d)*(x$hessian))%*%t(a$v)
    }
    x$options$no.se <- FALSE
    x
}


mrgradcheck <- function(x, options=x$options) {
    j <- mrjhess(x,options,only.jac=TRUE)$jac
    g <- colSums(j)
    sum(g*solve(crossprod(j),g))
}


mropt <- function(start,loglik) {
    ans <- nlminb(start,function(x) -sum(loglik(x)))
    if(ans$convergence) warning(paste("nlminb exits with code",ans$convergence))
    list(estimate=ans$par,maximum=ans$objective)
}


mroptions <- function(seed=round(runif(1,1,100000)), nrep=c(100,1000), no.se=FALSE, opt=mropt) {
    list(seed=seed,nrep=nrep,no.se=no.se,opt=opt)
}
    
mr <- function(x=rep(1,NROW(y)) , y, marginal, vcov, start, fixed, options=mroptions()) {
    ## arguments check
    if ( !inherits( marginal , "marginal.mr" ) ) stop("Unknown marginal")
    if ( !inherits( vcov , "cor.mr" ) ) stop("Unknown cor")
    ## mr object
    nb <- marginal$npar(x)
    ng <- vcov$npar
    not.na <- apply(cbind(y,x),1,function(z) !any(is.na(z)))
    m <- structure(list(y=as.matrix(y), x=as.matrix(x),
                        not.na=not.na, n=sum(not.na), marginal=marginal, vcov=vcov,
                        ibeta=1:nb, igamma=if (ng) (nb+1):(nb+ng) else NULL,
                        call=match.call()), class="mr")
    m$estimate <- if ( missing(start) ) c(marginal$start(m$y,m$x),vcov$start()) else start
    if (length(m$estimate) != nb+ng) stop("mismatch in the number of initial parameters")
    if ( missing(fixed) ) {
        m$fixed <- rep( NA , length(m$estimate) )
    } else {
        if (length(fixed) != length(m$estimate) ) stop("fixed has a wrong length")
        m$fixed <- fixed
    }
    m$options <- do.call(mroptions,options)
    # compute estimate
    m <- mrtruefit(m)
    # and s.e. and return
    if (m$options$no.se) m else mrjhess(m)
}


coef.mr <- function(object,...) object$estimate

logLik.mr <- function(object,...) {
    ans <- -object$maximum
    attr(ans,"df") <- sum( is.na(object$fixed) )
    class(ans) <- "logLik"
    ans
}

mrpinv <- function(h,tol) {
    h <- svd(h)
    idx <- h$d > sqrt(.Machine$double.eps)*h$d[1]
    ans <- h$u[,idx,drop=FALSE]%*%( (1/h$d[idx])*t(h$u[,idx,drop=FALSE]))
    attr(ans,"df") <- sum(idx)
    ans
}

estfun.mr <- function(x,...) x$jac
bread.mr <- function(x,...) mrpinv(x$hessian)*x$n

mrparvar <- function(x, type=c("hessian","sandwich","vscore","cluster","hac"),...) {
    if ( !inherits( x , "mr" ) ) stop("First argument must be a mr object")
    type <- match.arg(type)
    if (type=="sandwich") {
        if (is.null(x$hessian) || is.null(x$jac)) x <- mrjhess(x)
        h <- mrpinv(x$hessian)
        v <- h %*% crossprod(x$jac) %*% h
    } else if (type=="hessian")  {
        if (is.null(x$hessian) || is.null(x$jac)) x <- mrjhess(x)
        v <- mrpinv(x$hessian)
    } else if (type=="vscore")  {
        if (is.null(x$jac)) x <- mrjhess(x,only.jac=TRUE)
        v <- mrpinv(crossprod(x$jac))
    } else if (type=="cluster") {
        if (is.null(x$hessian) || is.null(x$jac)) x <- mrjhess(x)
        if (is.null(x$vcov$id)) stop("no cluster found")
        h <- mrpinv(x$hessian)
        v <- sapply(split(1:x$n,x$vcov$id),function(i) colSums(x$jac[i,])) 
        v <- h %*% (v %*% t(v) ) %*% h
    } else {
        if (!require(sandwich)) stop("HAC requires package sandwich")
        if (is.null(x$hessian) || is.null(x$jac)) x <- mrjhess(x)
        v <- vcovHAC(x,...)
    }
    ans <- matrix(0,length(x$estimate),length(x$estimate))
    ifree <- is.na(x$fixed)
    ans[ifree,ifree] <- v
    colnames(ans) <- rownames(ans) <- names(x$estimate)
    ans
}

mrse <- function(x,...) sqrt(diag(mrparvar(x,...)))

print.mr <- function (x, digits = max(3, getOption("digits") - 3), k = 2 ,...)
{
    cat("\nCall:", deparse(x$call, width.cutoff = 75), "", sep = "\n")
    cat("Parameters:\n")
    par <- x$estimate
    se <- if (x$options$no.se) NULL else mrse(x,...)
    xpar <- round( rbind( as.numeric(par) , s.e.=se ) , digits=digits )
    colnames(xpar) <- names(par)
    print.default(xpar, print.gap = 2)
    cat("\nlog likelihood = ", format(round(-x$maximum, 2)),
        ",  aic = ", format(round(AIC(x,k=k), 2)), "\n", sep = "")
    invisible(x)
}


residuals.mr <- function (object, type=c("conditional","marginal"),
                          method=c("random","interval","mid"),...) {
    type <- match.arg(type)
    method <- match.arg(method)
    is.int <- object$marginal$type == "integer"
    cond <- is.null(object$vcov$independent) && (type=="conditional") 
    if (is.int) {
        ##before saving/setting seed to obtain different randomization 
        if (method=="random") u <- runif(object$n) 
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            seed.keep <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
            on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
        }
        set.seed(object$options$seed)
    }
    a <- object$marginal$dp(object$y[object$not.na,],
                            object$x[object$not.na,,drop=FALSE],
                            object$estimate[object$ibeta])
    if (cond) {
        q <- object$vcov$chol(object$estimate[object$igamma],object$not.na)
        id <- if(is.null(object$vcov$id)) rep(1,object$n) else object$vcov$id
        lstrata <- rle(id)$length
        ipar <- c(object$n,max(object$options$nrep),length(lstrata),lstrata)
        a <- matrix(.C(if (is.int) "csampler" else "bsolver",
                       as.integer(ipar), as.double(unlist(q)), a=as.double(a),
                       double(object$n), NAOK=TRUE,dup=FALSE,package="mr")$a, object$n)
    }
    if (is.int) {
        if (method=="random") {
            ans <- rep(NA,sum(object$not.na))
            ans[object$not.na] <- qnorm(a[,2]-u*a[,1])            
        } else if (method=="interval") {
            ans <- matrix(NA,sum(object$not.na),2)
            ans[object$not.na,] <- qnorm(cbind(a[,2],a[,2]-a[,1]))
        } else {
            ans <- rep(NA,sum(object$not.na))
            ans[object$not.na] <- qnorm(a[,2]-a[,1]/2)            
        }
    } else {
            ans <- rep(NA,sum(object$not.na))
            ans[object$not.na] <- if (cond) a[,2] else qnorm(a[,2])            
    }
    ans
}


mrplotint <- function(r) {
    r <- as.ts(r)
    tm <- time(r)
    plot(NA,NA,xlim=range(tm),xlab="",ylim=range(r),ylab="pseudo residuals")
    for (i in 1:length(tm)) lines(c(tm[i],tm[i]),r[i,])
}

mrprof <- function(m , which , low , up , npoints = 10 , display = TRUE , alpha = 0.05 ) {
    if ( !inherits( m , "mr" ) ) stop("first argument must be a mr object")
    points <- seq( low , up , length = npoints )
    prof <- function(x,which,m) {
        m$fixed[which] <- x
        mr:::mrtruefit(m)$maximum
    }
    loglik <- sapply(points,prof,which,m)
    points <- c( points , m$par[which] )
    ord <- order( points)
    points <- points[ord]
    loglik <- -c(loglik,m$maximum)[ord]
    if ( display ) {
        npoints <- seq( low , up , length = 200 )
        plot( npoints , splinefun( points, loglik )(npoints) , type = "l" ,
             xlab = names(m$estimate)[which], ylab = "log-likelihood profile")
        grid()
        abline( h = max(loglik) - qchisq( 1 - alpha , 1 )/2 , lty = "dashed" )
    }
    invisible(list(points=points,profile=loglik))
}

mrhausman <- function(m, method=c("one-step","full"),
                      nboot=if (method=="one-step") 1000 else 200,
                      nrep=200) {
    if ( !inherits( m , "mr" ) ) stop("first argument must be a mr object")
    method <- match.arg(method)
    m$options$nrep <- nrep
    m$options$no.se <- TRUE
    ## ind lik
    mind <- mr(m$x,m$y,m$marg,mrind(),fixed=m$fixed[m$ibeta],options=m$options)
    ## differences in the marginal parameters
    betaml <- coef(m)[m$ibeta]
    betaind <- mind$estimate
    delta <- betaind - betaml
    ## simulator
    R <- m$vcov$chol( coef(m)[m$igamma] , m$not.na )
    if (!is.list(R)) R <- list(R)
    R <- lapply(R,t)
    y <- m$y[m$not.na,]
    x <- m$x[m$not.na,,drop=FALSE]
    sim <- function() {
        u <- pnorm(unlist(lapply(R, function(x) x %*% rnorm(NROW(x)))))
        m$marginal$sim(u,y,x,betaml)
    }
    ## bootstrap estimate of the the variance of the difference
    vdelta <- var(if (method=="one-step") dh1(m,mind,sim,nboot) else dhf(m,mind,sim,nboot))
    ## return
    se <- sqrt(diag(vdelta))
    htab <- cbind(ml=betaml,ind.lik=betaind,diff=delta,s.e.=se,z=delta/se)
    rownames(htab) <- names(betaml)
    vdelta <- mrpinv(vdelta)
    hall <- sum(delta*(vdelta%*%delta))
    hall <- c(stat=hall,df=attr(vdelta,"df"),
              p.value=pchisq(hall,attr(vdelta,"df"),lower.tail=FALSE))
    list(overall=hall,parameters=htab,
         options=c(method=method,nboot=nboot,nrep=nrep))    
}

dh1 <- function(m, mind, sim, nboot) {
    mind$estimate <- coef(m)[m$ibeta]
    nb <- Hind <- 0
    sml <- matrix(0,nboot,length(coef(m))) 
    sind <- matrix(0,nboot,length(coef(mind))) 
    while ( nb < nboot ) {
        nb <- nb+1
        m$y[m$not.na,] <- mind$y[mind$not.na,] <- sim()
        sml[nb,] <- colSums(mrjhess(m, only.jac=TRUE)$jac)
        mind <- mrjhess(mind, only.jac=FALSE)
        sind[nb,] <- colSums(mind$jac) 
        Hind <- Hind + (mind$hessian-Hind)/nb
    }
    (sml%*%mrpinv(crossprod(sml/sqrt(nboot))))[,m$ibeta]-sind%*%mrpinv(Hind)
}

dhf <- function(m, mind, sim, nboot) {
    theta <- m$estimate
    m$options$opt <- mind$options$opt <- function(start,loglik) {
        ans <- nlminb(start,function(x) -sum(loglik(x)),
                      control=list(rel.tol=1e-4,iter.max=3*length(theta)))
        list(estimate=ans$par,maximum=ans$objective)
    }
    d <- matrix(0,nboot,length(coef(mind)))
    nb <- 0
    while (nb < nboot) {
        nb <- nb+1
        m$y[m$not.na,] <- mind$y[mind$not.na,] <- sim()
        m$estimate <- theta
        m <- mrtruefit(m)
        mind$estimate <- mind$marginal$start(mind$y,mind$x)
        mind <- mrtruefit(mind)
        d[nb,] <- coef(m)[m$ibeta]-coef(mind)
    }
    d
}
