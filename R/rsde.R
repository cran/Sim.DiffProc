## Fri Mar 07 18:39:01 2014
## Original file Copyright © 2014 A.C. Guidoum, K. Boukhetala
## This file is part of the R package Sim.DiffProc
## Department of Probabilities & Statistics
## Faculty of Mathematics
## University of Science and Technology Houari Boumediene
## BP 32 El-Alia, U.S.T.H.B, Algeris
## Algeria

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## A copy of the GNU General Public License is available at
## http://www.r-project.org/Licenses/
## Unlimited use and distribution (see LICENCE).
###################################################################################################



#####
##### rsde1d

rsde1d <- function(N, ...)  UseMethod("rsde1d")

rsde1d.default <- function(N =100,M=10,x0=0,t0=0,T=1,Dt,tau=0.5,drift,diffusion,alpha=0.5,mu=0.5,
                     type=c("ito","str"), method=c("euler","milstein","predcorr",
                     "smilstein","taylor","heun","rk1","rk2","rk3"),...)
                     {
    if (!is.numeric(x0)) stop("'x0' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T))) stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0)) stop(" 'M' must be a positive integer ")
    if (any(!is.expression(drift) || !is.expression(diffusion) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions")
    if (missing(type)) type <- "ito"
    method <- match.arg(method)
    if (method =="predcorr"){
    if(any(alpha > 1 || alpha < 0)) stop("please use '0 <= alpha <= 1' ")
    if(any(mu > 1 || mu < 0))       stop("please use '0 <= mu <= 1' ")
                            }
    if (any(t0 < 0 || T < 0 || T <= t0) ) stop(" please use positive times! (0 <= t0 < T) ")
    if (missing(Dt)) {
        t <- seq(t0, T, length = N + 1)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
        T <- t[N + 1]
    }
    if (any(T < tau || t0 > tau) )   stop( " please use 't0 <= tau <= T'")
    Dt <- (T - t0)/N 
    Y <- snssde1d(N,M,x0,t0,T,Dt,drift,diffusion,alpha,mu,type, method,...)
    X <- as.vector(Y$X[which(t==tau),])
    if (length(X) == 0){
    F   <- lapply(1:M,function(i) approxfun(time(Y),Y$X[,i]))
    X   <- sapply(1:length(F),function(i) F[[i]](tau)) 
    }
    plot(Y,plot.type="single")
    points(rep(tau,M),X,col = 2, pch = 16,cex=0.8)
    Axis(at = tau, side=1,col=2, labels = bquote(tau))
    lines(c(tau,tau),c(0,max(X)),lty=2,col=2)
    lines(c(tau,tau),c(0,min(X)),lty=2,col=2)
    structure(list(t=tau,x=X),class="rsde1d")
}

###

bconfint.rsde1d <- function(x,level=0.95,...)
                 {
     class(x) <-"rsde1d"
     bconfint(x$x,level=level,...)
}

skewness.rsde1d <- function(x,...)
                    {
    class(x) <- "rsde1d"
    skewness(x$x) 
}

kurtosis.rsde1d <- function(x,...)
                    {
    class(x) <- "rsde1d"
    kurtosis(x$x) 
}

median.rsde1d <- function(x,...)
                    {
    class(x) <- "rsde1d"
    median(x$x) 
}

mean.rsde1d <- function(x,...)
                    {
    class(x) <- "rsde1d"
    mean(x$x) 
}

quantile.rsde1d <- function(x,...)
             {
    class(x) <- "rsde1d"
    quantile(x$x,...)
}

moment.rsde1d <- function(x,order = 2,...)
             {
    class(x) <- "rsde1d"
    Mom <- data.frame(do.call("cbind",lapply(1:length(order), function(j) moment(x$x,order=order[j]))))
    row.names(Mom) <- paste("X(t=",x$t,")",sep="")
    names(Mom) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    return(Mom)
}

summary.rsde1d <- function(object, ...)
           {
	x <- object	   
    class(x) <- "rsde1d"
    cat("\n\tMonte-Carlo Statistics for X(t) at t = ",x$t,"\n",
        sep="")
    res <- data.frame(matrix(c(mean(x$x),var(x$x),median(x$x),quantile(x$x,0.25),quantile(x$x,0.75),
                               skewness(x$x),kurtosis(x$x),moment(x$x,order=2),moment(x$x,order=3),
                               moment(x$x,order=4),moment(x$x,order=5),bconfint(x$x)[1],bconfint(x$x)[2]),
                               ncol=1))
    row.names(res) <- paste(c("Mean","Variance","Median","First quartile","Third quartile",
                              "Skewness","Kurtosis","Moment of order 2","Moment of order 3",
                              "Moment of order 4","Moment of order 5","Bound conf Inf (95%)","Bound conf Sup (95%)"),sep="")
    names(res) <- paste(c(""),sep="")
    print(res, quote = FALSE, right = TRUE,...)
    invisible(x)
}



#####
##### rsde2d

rsde2d <- function(N, ...)  UseMethod("rsde2d")

rsde2d.default <- function(N =100,M=10,x0=0,y0=0,t0=0,T=1,Dt,tau=0.5,driftx,diffx,drifty,diffy,alpha=0.5,mu=0.5,
                     type=c("ito","str"), method=c("euler","milstein","predcorr",
                     "smilstein","taylor","heun","rk1","rk2","rk3"),...)
                    { 
    if (any(!is.numeric(x0) || !is.numeric(y0))) stop("'x0' and 'y0' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T)))  stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0)) stop(" 'M' must be a positive integer ")
    if (any(!is.expression(driftx) || !is.expression(diffx) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions")
    if (any(!is.expression(drifty) || !is.expression(diffy) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions")
    if (missing(type)) type <- "ito"
    method <- match.arg(method)
    if (method =="predcorr"){
    if(any(alpha > 1 || alpha < 0)) stop("please use '0 <= alpha <= 1' ")
    if(any(mu > 1 || mu < 0))       stop("please use '0 <= mu <= 1' ")
                            }
    if (any(t0 < 0 || T < 0 || T <= t0) ) stop(" please use positive times! (0 <= t0 < T) ")
    if (missing(Dt)) {
        t <- seq(t0, T, length = N + 1)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
        T <- t[N + 1]
    }
    if (any(T < tau || t0 > tau) )   stop( " please use 't0 <= tau <= T'")
    Dt <- (T - t0)/N 
    Sim <- lapply(1:M,function(i) snssde2d(N,x0,y0,t0,T,Dt,driftx,diffx,drifty,diffy,alpha,mu,type, method,...)$XY)
    xy  <- do.call("rbind",lapply(1:length(Sim),function(i) Sim[[i]][which(t==tau),]))
    if (length(xy[,1])== 0){
        F1   <- lapply(1:M,function(i) approxfun(t,Sim[[i]][,1]))
        F2   <- lapply(1:M,function(i) approxfun(t,Sim[[i]][,2]))
        x    <- do.call("rbind",lapply(1:length(F1),function(i) F1[[i]](tau)))
        y    <- do.call("rbind",lapply(1:length(F2),function(i) F2[[i]](tau)))
        xy   <- matrix(c(x,y),nrow=M)
         }
    structure(list(t=tau,x=xy[,1],y=xy[,2]),class="rsde2d")
}

##
##

bconfint.rsde2d <- function(x,level=0.95,...)
             {
    class(x) <- "rsde2d"
    Bcon <- t(data.frame(do.call("cbind",lapply(2:3,function(i) bconfint(x[i][[1]],level=level)))))
    row.names(Bcon) <- paste(c("X(t=","Y(t="),c(rep(x$t,2)),c(rep(")",2)),sep="")
    return(Bcon)
}

skewness.rsde2d <- function(x,...)
             {
    class(x) <- "rsde2d"
    Skew <- data.frame(do.call("cbind",lapply(2:3,function(i) skewness(x[i][[1]]))),row.names = "")
    names(Skew) <- paste(c("X(t=","Y(t="),c(rep(x$t,2)),c(rep(")",2)),sep="")
    return(Skew)
}

kurtosis.rsde2d <- function(x,...)
             {
    class(x) <- "rsde2d"
    Kurt <- data.frame(do.call("cbind",lapply(2:3,function(i) kurtosis(x[i][[1]]))),row.names = "")
    names(Kurt) <- paste(c("X(t=","Y(t="),c(rep(x$t,2)),c(rep(")",2)),sep="")
    return(Kurt)
}

median.rsde2d <- function(x,...)
             {
    class(x) <- "rsde2d"
    Med <- data.frame(do.call("cbind",lapply(2:3,function(i) median(x[i][[1]]))),row.names = "")
    names(Med) <- paste(c("X(t=","Y(t="),c(rep(x$t,2)),c(rep(")",2)),sep="")
    return(Med)
}

mean.rsde2d <- function(x,...)
             {
    class(x) <- "rsde2d"
    Mean <- data.frame(do.call("cbind",lapply(2:3,function(i) mean(x[i][[1]]))),row.names = "")
    names(Mean) <- paste(c("X(t=","Y(t="),c(rep(x$t,2)),c(rep(")",2)),sep="")
    return(Mean)
}

quantile.rsde2d <- function(x,...)
             {
    class(x) <- "rsde2d"
    Qun <- t(data.frame(do.call("cbind",lapply(2:3,function(i) quantile(x[i][[1]],...)))))
    row.names(Qun) <- paste(c("X(t=","Y(t="),c(rep(x$t,2)),c(rep(")",2)),sep="")
    return(Qun)
}

moment.rsde2d <- function(x,order = 2,...)
             {
    class(x) <- "rsde2d"
    Mom <- data.frame(do.call("cbind",lapply(1:length(order), function(j) sapply(2:3,function(i) moment(x[i][[1]],order=order[j])))))
    row.names(Mom) <- paste(c("X(t=","Y(t="),c(rep(x$t,2)),c(rep(")",2)),sep="")
    names(Mom) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    return(Mom)
}

summary.rsde2d <- function(object, ...)
           {
	x <- object	   
    class(x) <- "rsde2d"
    cat("\n\tMonte-Carlo Statistics for X(t) and Y(t) at t = ",x$t,"\n\n",
        sep="")
    res <- data.frame(matrix(c(mean(x$x),var(x$x),median(x$x),quantile(x$x,0.25),quantile(x$x,0.75),
                               skewness(x$x),kurtosis(x$x),moment(x$x,order=2),moment(x$x,order=3),
                               moment(x$x,order=4),moment(x$x,order=5),bconfint(x$x)[1],bconfint(x$x)[2],
                               mean(x$y),var(x$y),median(x$y),quantile(x$y,0.25),quantile(x$y,0.75),
                               skewness(x$y),kurtosis(x$y),moment(x$y,order=2),moment(x$y,order=3),
                               moment(x$y,order=4),moment(x$y,order=5),bconfint(x$y)[1],bconfint(x$y)[2]),
                               ncol=2))
    row.names(res) <- paste(c("Mean","Variance","Median","First quartile","Third quartile",
                              "Skewness","Kurtosis","Moment of order 2","Moment of order 3",
                              "Moment of order 4","Moment of order 5","Bound conf Inf (95%)","Bound conf Sup (95%)"),sep="")
    names(res) <- paste(c("X","Y"),sep="")
    print(res, quote = FALSE, right = TRUE,...)
    invisible(x)
}


#####
##### rsde3d

rsde3d <- function(N, ...)  UseMethod("rsde3d")

rsde3d.default <- function(N =100,M=10,x0=0,y0=0,z0=0,t0=0,T=1,Dt,tau=0.5,driftx,diffx,drifty,diffy,driftz,diffz,
                             alpha=0.5,mu=0.5,type=c("ito","str"), method=c("euler","milstein",
                             "predcorr","smilstein","taylor","heun","rk1","rk2","rk3"),...)
                    { 
    if (any(!is.numeric(x0) || !is.numeric(y0) || !is.numeric(z0))) stop("'x0','y0' and 'z0' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T))) stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0)) stop(" 'M' must be a positive integer ")
    if (any(!is.expression(driftx) || !is.expression(diffx) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions")
    if (any(!is.expression(drifty) || !is.expression(diffy) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions")
    if (any(!is.expression(driftz) || !is.expression(diffz) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions")
    if (missing(type)) type <- "ito"
    method <- match.arg(method)
    if (method =="predcorr"){
    if(any(alpha > 1 || alpha < 0)) stop("please use '0 <= alpha <= 1' ")
    if(any(mu > 1 || mu < 0))       stop("please use '0 <= mu <= 1' ")
                            }
    if (any(t0 < 0 || T < 0 || T <= t0) ) stop(" please use positive times! (0 <= t0 < T) ")
    if (missing(Dt)) {
        t <- seq(t0, T, length = N + 1)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
        T <- t[N + 1]
    }
    if (any(T < tau || t0 > tau) )   stop( " please use 't0 <= tau <= T'")
    Dt <- (T - t0)/N 
    Sim <- lapply(1:M,function(i) snssde3d(N,x0,y0,z0,t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,alpha,mu,type, method,...)$XYZ)
    xyz <- do.call("rbind",lapply(1:length(Sim),function(i) Sim[[i]][which(t==tau),]))
    if (length(xyz[,1])== 0){
         F1   <- lapply(1:M,function(i) approxfun(time(Sim[[1]]),Sim[[i]][,1]))
         F2   <- lapply(1:M,function(i) approxfun(time(Sim[[1]]),Sim[[i]][,2]))
         F3   <- lapply(1:M,function(i) approxfun(time(Sim[[1]]),Sim[[i]][,3]))
         x    <- do.call("rbind",lapply(1:length(F1),function(i) F1[[i]](tau)))
         y    <- do.call("rbind",lapply(1:length(F2),function(i) F2[[i]](tau)))
         z    <- do.call("rbind",lapply(1:length(F3),function(i) F3[[i]](tau)))
         xyz  <- matrix(c(x,y,z),nrow=M)
       }
    structure(list(t=tau,x=xyz[,1],y=xyz[,2],z=xyz[,3]),class="rsde3d")
}

###

bconfint.rsde3d <- function(x,level=0.95,...)
             {
    class(x) <- "rsde3d"
    Bcon <- t(data.frame(do.call("cbind",lapply(2:4,function(i) bconfint(x[i][[1]],level=level)))))
    row.names(Bcon) <- paste(c("X(t=","Y(t=","Z(t="),c(rep(x$t,3)),c(rep(")",3)),sep="")
    return(Bcon)
}

skewness.rsde3d <- function(x,...)
             {
    class(x) <- "rsde3d"
    Skew <- data.frame(do.call("cbind",lapply(2:4,function(i) skewness(x[i][[1]]))),row.names = "")
    names(Skew) <- paste(c("X(t=","Y(t="),c(rep(x$t,3)),c(rep(")",3)),sep="")
    return(Skew)
}

kurtosis.rsde3d <- function(x,...)
             {
    class(x) <- "rsde3d"
    Kurt <- data.frame(do.call("cbind",lapply(2:4,function(i) kurtosis(x[i][[1]]))),row.names = "")
    names(Kurt) <- paste(c("X(t=","Y(t=","Z(t="),c(rep(x$t,3)),c(rep(")",3)),sep="")
    return(Kurt)
}

median.rsde3d <- function(x,...)
             {
    class(x) <- "rsde3d"
    Med <- data.frame(do.call("cbind",lapply(2:4,function(i) median(x[i][[1]]))),row.names = "")
    names(Med) <- paste(c("X(t=","Y(t="),c(rep(x$t,3)),c(rep(")",3)),sep="")
    return(Med)
}

mean.rsde3d <- function(x,...)
             {
    class(x) <- "rsde3d"
    Mean <- data.frame(do.call("cbind",lapply(2:4,function(i) mean(x[i][[1]]))),row.names = "")
    names(Mean) <- paste(c("X(t=","Y(t=","Z(t="),c(rep(x$t,3)),c(rep(")",3)),sep="")
    return(Mean)
}

quantile.rsde3d <- function(x,...)
             {
    class(x) <- "rsde3d"
    Qun <- t(data.frame(do.call("cbind",lapply(2:4,function(i) quantile(x[i][[1]],...)))))
    row.names(Qun) <- paste(c("X(t=","Y(t="),c(rep(x$t,3)),c(rep(")",3)),sep="")
    return(Qun)
}

moment.rsde3d <- function(x,order = 2,...)
             {
    class(x) <- "rsde3d"
    Mom <- data.frame(do.call("cbind",lapply(1:length(order), function(j) sapply(2:4,function(i) moment(x[i][[1]],order=order[j])))))
    row.names(Mom) <- paste(c("X(t=","Y(t=","Z(t="),c(rep(x$t,3)),c(rep(")",3)),sep="")
    names(Mom) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    return(Mom)
}

summary.rsde3d <- function(object, ...)
           {
	x <- object	   
    class(x) <- "rsde3d"
    cat("\n\tMonte-Carlo Statistics for X(t), Y(t) and Z(t) at t = ",x$t,"\n\n",
        sep="")
    res <- data.frame(matrix(c(mean(x$x),var(x$x),median(x$x),quantile(x$x,0.25),quantile(x$x,0.75),
                               skewness(x$x),kurtosis(x$x),moment(x$x,order=2),moment(x$x,order=3),
                               moment(x$x,order=4),moment(x$x,order=5),bconfint(x$x)[1],bconfint(x$x)[2],
                               mean(x$y),var(x$y),median(x$y),quantile(x$y,0.25),quantile(x$y,0.75),
                               skewness(x$y),kurtosis(x$y),moment(x$y,order=2),moment(x$y,order=3),
                               moment(x$y,order=4),moment(x$y,order=5),bconfint(x$y)[1],bconfint(x$y)[2],
                               mean(x$z),var(x$z),median(x$z),quantile(x$z,0.25),quantile(x$z,0.75),
                               skewness(x$z),kurtosis(x$z),moment(x$z,order=2),moment(x$z,order=3),
                               moment(x$z,order=4),moment(x$z,order=5),bconfint(x$z)[1],bconfint(x$z)[2]),
                               ncol=3))
    row.names(res) <- paste(c("Mean","Variance","Median","First quartile","Third quartile",
                              "Skewness","Kurtosis","Moment of order 2","Moment of order 3",
                              "Moment of order 4","Moment of order 5","Bound conf Inf (95%)","Bound conf Sup (95%)"),sep="")
    names(res) <- paste(c("X","Y","Z"),sep="")
    print(res, quote = FALSE, right = TRUE,...)
    invisible(x)
}
