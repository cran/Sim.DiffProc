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
##### snssde1D

snssde1d <- function(N, ...)  UseMethod("snssde1d")

snssde1d.default <- function(N =100,M=1,x0=0,t0=0,T=1,Dt,drift,diffusion,alpha=0.5,mu=0.5,
                     type=c("ito","str"), method=c("euler","milstein","predcorr",
                     "smilstein","taylor","heun","rk1","rk2","rk3"),...)
        {
    if (!is.numeric(x0)) stop("'x0' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T))) stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0))  stop(" 'M' must be a positive integer ")
    if (any(!is.expression(drift) || !is.expression(diffusion) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions")
    if (missing(type)) type <- "ito"
    method <- match.arg(method)
    if (method =="predcorr"){
    if (any(alpha > 1 || alpha < 0)) stop("please use '0 <= alpha <= 1' ")
    if (any(mu > 1 || mu < 0))       stop("please use '0 <= mu <= 1' ")
                            }
    if (any(t0 < 0 || T < 0 || T <= t0) ) stop(" please use positive times! (0 <= t0 < T) ")
    if (missing(Dt)) {
        t <- seq(t0, T, length = N + 1)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
        T <- t[N + 1]
    }
    Dt <- (T - t0)/N    
    if (method=="euler")         {res <- .Euler1D(N,M,x0,t0=t0,T=T,Dt=Dt,drift,diffusion,type)}
    else if (method=="predcorr") {res <- .PredCorr1D(N,M,x0,t0,T,Dt,alpha,mu,drift,diffusion,type)}
    else if (method=="milstein") {res <- .Milstein1D(N,M,x0,t0,T,Dt,drift,diffusion,type)}
    else if (method=="smilstein"){res <- .SMilstein1D(N,M,x0,t0,T,Dt,drift,diffusion,type)}
    else if (method=="taylor")   {res <- .STS1D(N,M,x0,t0,T,Dt,drift,diffusion,type)}
    else if (method=="heun")     {res <- .Heun1D(N,M,x0,t0,T,Dt,drift,diffusion,type)}
    else if (method=="rk1")      {res <- .RK1D(N,M,x0,t0,T,Dt,drift,diffusion,type,order=1)}
    else if (method=="rk2")      {res <- .RK1D(N,M,x0,t0,T,Dt,drift,diffusion,type,order=2)}
    else if (method=="rk3")      {res <- .RK1D(N,M,x0,t0,T,Dt,drift,diffusion,type,order=3)}
    structure(list(X=res$X,drift=drift[[1]], diffusion=diffusion[[1]],type=type,method=method, 
                   x0=x0, N=N, M=M,Dt=Dt,t0=t0,T=T),class="snssde1d")
}

###

print.snssde1d <- function(x, digits=NULL, ...)
           {
    class(x) <- "snssde1d"
    if (x$method=="euler")         {sch <- "Euler scheme of order 0.5"}
    else if (x$method=="milstein") {sch <- "Milstein scheme of order 1"}
    else if (x$method=="predcorr") {sch <- "Predictor-corrector method of order 1"}
    else if (x$method=="smilstein"){sch <- "Second Milstein scheme"}
    else if (x$method=="taylor")   {sch <- "Ito-Taylor scheme of order 1.5"}
    else if (x$method=="heun")     {sch <- "Heun scheme of order 2"}
    else if (x$method=="rk1")      {sch <- "Runge-Kutta method of order 1"}
    else if (x$method=="rk2")      {sch <- "Runge-Kutta method of order 2"}
    else if (x$method=="rk3")      {sch <- "Runge-Kutta method of order 3"}
    if(x$type=="ito"){
    cat("Ito Sde 1D:","\n",        
        "\t| dx = ", deparse(x$drift)," * dt + ", deparse(x$diffusion)," * dw","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        "\t| Number of simulation","\t| M  = ",format(x$M,digits=digits),".","\n",
        "\t| Initial value","\t\t| x0 = ",format(x$x0,digits=digits),".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",       
        sep="")}else{
    cat("Stratonovich Sde 1D:","\n",
        "\t| dx = ", deparse(x$drift)," * dt + ", deparse(x$diffusion)," o dw","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        "\t| Number of simulation","\t| M  = ",format(x$M,digits=digits),".","\n",
        "\t| Initial value","\t\t| x0 = ",format(x$x0,digits=digits),".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}
    invisible(x)
}

##
## Plot

plot.snssde1d <- function(x,...)
                 {
    class(x) <- "snssde1d"
    X <- x$X
    plot(X,...)
}

lines.snssde1d <- function(x,...)
                 {
    class(x) <- "snssde1d"
    X <- x$X
    for (i in 1:dim(X)[2]){
    lines(time(x),X[,i],...)}
}

points.snssde1d <- function(x,...)
                 {
    class(x) <- "snssde1d"
    X <- x$X
    for (i in 1:dim(X)[2]){
    points(time(x),X[,i],...)}
}

add.mean <- function(x,lty=NULL,lwd=NULL,col=NULL,cex=NULL,...)
                 {
    class(x) <- "snssde1d"
    X <- x$X
    if (is.null(lty)) {lty = 1}
    if (is.null(lwd)) {lwd = 1}
    if (is.null(col)) {col = 2}
    if (is.null(cex)) {cex = 0.8}
    lines(as.vector(time(X)),rowMeans(X),lwd=lwd,lty=lty,col=col,...)
    legend("topright",c("mean path"),inset = .01,lty=lty,col=col,lwd=lwd,cex=cex,...)
}


add.bconfint.snssde1d <- function(x,level=0.95,lty=NULL,lwd=NULL,col=NULL,cex=NULL,...)
                 {
    class(x) <- "snssde1d"
    if (is.null(lty)) {lty = 1}
    if (is.null(lwd)) {lwd = 1}
    if (is.null(col)) {col = 4}
    if (is.null(cex)) {cex = 0.8}
    lines(time(x),bconfint(x,level)[,1],lwd=lwd,lty=lty,col=col,...)
    lines(time(x),bconfint(x,level)[,2],lwd=lwd,lty=lty,col=col,...)
    legend("topleft",c(paste("bound of",level*100,"% confidence")),inset = .01,lty=lty,col=col,lwd=lwd,cex=cex,...)
}

##
## summary

summary.snssde1d <- function(object, ...)
           {
	x <- object	   
    class(x) <- "snssde1d"
    cat("\n\tMonte-Carlo Statistics for X(t) at final time T = ",x$T,"\n\n",        
        "| Process mean","\t\t\t = ",mean(x$X[which(time(x)==x$T),]),"\n",
        "| Process variance","\t\t = ",var(x$X[which(time(x)==x$T),]),"\n",
        "| Process median","\t\t = ",median(x$X[which(time(x)==x$T),]),"\n",
        "| Process first quartile","\t = ",quantile(x$X[which(time(x)==x$T),],0.25),"\n",
        "| Process third quartile","\t = ",quantile(x$X[which(time(x)==x$T),],0.75),"\n",
        "| Process skewness","\t\t = ",skewness(x$X[which(time(x)==x$T),]),"\n",
        "| Process kurtosis","\t\t = ",kurtosis(x$X[which(time(x)==x$T),]),"\n",
        "| Process moment of order 2","\t = ",moment(x$X[which(time(x)==x$T),],order=2),"\n",
        "| Process moment of order 3","\t = ",moment(x$X[which(time(x)==x$T),],order=3),"\n",
        "| Process moment of order 4","\t = ",moment(x$X[which(time(x)==x$T),],order=4),"\n",
        "| Process moment of order 5","\t = ",moment(x$X[which(time(x)==x$T),],order=5),"\n",
        "| Bound of confidence (95%)","\t = [",bconfint(x$X[which(time(x)==x$T),])[1],",",bconfint(x$X[which(time(x)==x$T),])[2],"]","\n",
        "  for the trajectories","\n",
        sep="")
    invisible(x)
}

mean.snssde1d <- function(x,...)
                    {
    class(x) <- "snssde1d"
    rowMeans(x$X,...)
}

skewness.snssde1d <- function(x,...)
                    {
    class(x) <- "snssde1d"
    Skew <- data.frame(sapply(1:(x$N+1),function(i) skewness(x$X[i,]) ))
    row.names(Skew) <- paste("X(t=",time(x),")",sep="")
    names(Skew) <- paste(c("Skewness"),sep="")
    return(Skew)
}

kurtosis.snssde1d <- function(x,...)
                    {
    class(x) <- "snssde1d"
    kurt <- data.frame(sapply(1:(x$N+1),function(i) kurtosis(x$X[i,]) ))
    row.names(kurt) <- paste("X(t=",time(x),")",sep="")
    names(kurt) <- paste(c("Kurtosis"),sep="")
    return(kurt)
}

median.snssde1d <- function(x,...)
                    {
    class(x) <- "snssde1d"
    Med <- data.frame(sapply(1:(x$N+1),function(i) median(x$X[i,]) ))
    row.names(Med) <- paste("X(t=",time(x),")",sep="")
    names(Med) <- paste(c("Median"),sep="")
    return(Med)
}

quantile.snssde1d <- function(x,...)
                    {
    class(x) <- "snssde1d"
    Qun <- t(data.frame(do.call("cbind",lapply(1:(x$N+1),function(i) quantile(x$X[i,],...) ))))
    row.names(Qun) <- paste("X(t=",time(x),")",sep="")
    return(Qun)
}

moment.snssde1d <- function(x,order = 2,...)
                    {
    if (any(!is.numeric(order)  || (order - floor(order) > 0) || order < 1)) stop(" 'order' must be a positive integer")
    class(x) <- "snssde1d"
    Mom <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
               sapply(1:length(order), function(j) moment(x$X[i,],order=order[j],...)))))
    row.names(Mom) <- paste("X(t=",time(x),")",sep="")
    names(Mom) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    return(Mom)
}

bconfint.snssde1d <- function(x,level = 0.95,...)
                    {
    class(x) <- "snssde1d"
    conf <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
                      quantile(x$X[i,], c(0.5*(1-level), 1-0.5*(1-level)),type=8) ) ) )
    row.names(conf) <- paste("X(t=",time(x),")",sep="")
    names(conf) <- paste(c(0.5*(1-level)*100,(1-(1-level)/2)*100),c(" %"," %"),sep="")
    return(conf)
}

#var.x <- function(x,...)
#                    {
#    class(x) <- "snssde1d" 
#    apply(x$X,1,var)
#}

time.snssde1d <- function(x,...)
                    {
    class(x) <- "snssde1d"
    as.vector(time(x$X))
}


################################################################################
################################################################################
##### snssde2d

snssde2d <- function(N, ...)  UseMethod("snssde2d")

snssde2d.default <- function(N =100,x0=0,y0=0,t0=0,T=1,Dt,driftx,diffx,drifty,diffy,alpha=0.5,mu=0.5,
                     type=c("ito","str"), method=c("euler","milstein","predcorr",
                     "smilstein","taylor","heun","rk1","rk2","rk3"),...)
        {
    if (any(!is.numeric(x0) || !is.numeric(y0))) stop("'x0' and 'y0' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T))) stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) stop(" 'N' must be a positive integer ")
    if (any(!is.expression(driftx) || !is.expression(diffx) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions")
    if (any(!is.expression(drifty) || !is.expression(diffy) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions")
    if (missing(type)) type <- "ito"
    method <- match.arg(method)
    if (method =="predcorr"){
    if (any(alpha > 1 || alpha < 0)) stop("please use '0 <= alpha <= 1' ")
    if (any(mu > 1 || mu < 0))       stop("please use '0 <= mu <= 1' ")
                            }
    if (any(t0 < 0 || T < 0 || T <= t0) ) stop(" please use positive times! (0 <= t0 < T) ")
    if (missing(Dt)) {
        t <- seq(t0, T, length = N + 1)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
        T <- t[N + 1]
    }
    Dt <- (T - t0)/N    
    if (method=="euler")         {res <- .Euler2D(N,x0,y0,t0,T,Dt,driftx,diffx,drifty,diffy,type)}
    else if (method=="predcorr") {res <- .PredCorr2D(N,x0,y0,t0,T,Dt,alpha,mu,driftx,diffx,drifty,diffy,type)}
    else if (method=="milstein") {res <- .Milstein2D(N,x0,y0,t0,T,Dt,driftx,diffx,drifty,diffy,type)}
    else if (method=="smilstein"){res <- .SMilstein2D(N,x0,y0,t0,T,Dt,driftx,diffx,drifty,diffy,type)}
    else if (method=="taylor")   {res <- .STS2D(N,x0,y0,t0,T,Dt,driftx,diffx,drifty,diffy,type)}
    else if (method=="heun")     {res <- .Heun2D(N,x0,y0,t0,T,Dt,driftx,diffx,drifty,diffy,type)}
    else if (method=="rk1")      {res <- .RK2D(N,x0,y0,t0,T,Dt,driftx,diffx,drifty,diffy,type,order=1)}
    else if (method=="rk2")      {res <- .RK2D(N,x0,y0,t0,T,Dt,driftx,diffx,drifty,diffy,type,order=2)}
    else if (method=="rk3")      {res <- .RK2D(N,x0,y0,t0,T,Dt,driftx,diffx,drifty,diffy,type,order=3)}
    structure(list(XY=res$X, driftx=driftx[[1]], diffx=diffx[[1]],drifty=drifty[[1]], diffy=diffy[[1]],type=type,method=method, 
                   x0=x0,y0=y0, N=N,Dt=Dt,t0=t0,T=T),class="snssde2d")
}

###

print.snssde2d <- function(x, digits=NULL, ...)
           {
    class(x) <- "snssde2d"
    if (x$method=="euler")         {sch <- "Euler scheme of order 0.5"}
    else if (x$method=="milstein") {sch <- "Milstein scheme of order 1"}
    else if (x$method=="predcorr") {sch <- "Predictor-corrector method of order 1"}
    else if (x$method=="smilstein"){sch <- "Second Milstein scheme of order 1.5"}
    else if (x$method=="taylor")   {sch <- "Ito-Taylor scheme of order 2"}
    else if (x$method=="heun")     {sch <- "Heun scheme of order 2"}
    else if (x$method=="rk1")      {sch <- "Runge-Kutta method of order 1"}
    else if (x$method=="rk2")      {sch <- "Runge-Kutta method of order 2"}
    else if (x$method=="rk3")      {sch <- "Runge-Kutta method of order 3"}
    if(x$type=="ito"){
    cat("Ito Sde 2D:","\n",
        "\t| dx = ", deparse(x$driftx)," * dt + ", deparse(x$diffx)," * dw1","\n", 
        "\t| dy = ", deparse(x$drifty)," * dt + ", deparse(x$diffy)," * dw2","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        "\t| Initial values","\t| (x0,y0) = ","(",format(x$x0,digits=digits),",",format(x$y0,digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}else{
    cat("Stratonovich Sde 2D:","\n",
        "\t| dx = ", deparse(x$driftx)," * dt + ", deparse(x$diffx)," o dw1","\n", 
        "\t| dy = ", deparse(x$drifty)," * dt + ", deparse(x$diffy)," o dw2","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        "\t| Initial values","\t| (x0,y0) = ","(",format(x$x0,digits=digits),",",format(x$y0,digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}
    invisible(x)
}

##
## Plot

plot.ssnssde2d <- function(x,plot.type = NULL,col = NULL,lty = NULL,lwd=NULL,main=NULL,...)
                 {
    class(x) <- "snssde2d"
    X <- x$XY
    if (is.null(lty)) {lty = 1}
    if (is.null(lwd)) {lwd = 1}
    if (is.null(col)) {col =1:2}
    if (is.null(main)){main=""}
    if (!is.null(plot.type)){
    plot(X,plot.type="single",col=col,lty=lty,lwd=lwd,main=main,...)
    legend("topright",c(expression(X[t]),expression(Y[t])),inset = .01,col=col,lty=lty,lwd=lwd)
    }else{
    plot(X,col=col,lty=lty,lwd=lwd,main=main,...)
     }
}

plot.snssde2d <- function(x,...) plot.ssnssde2d(x,...)

lines.snssde2d <- function(x,...)
                 {
    class(x) <- "snssde2d"
    X <- x$XY
    for (i in 1:dim(X)[2]){
    lines(time(x),X[,i],...)}
}

points.snssde2d <- function(x,...)
                 {
    class(x) <- "snssde2d"
    X <- x$XY
    for (i in 1:dim(X)[2]){
    points(time(x),X[,i],...)}
}


plot2d.snssde2d <- function(x,...)
                 {
    class(x) <- "snssde2d"
    X <- x$XY[,1]
    Y <- x$XY[,2]
    plot2d(X,Y,...)
    for(i in 3:4) axis(i)
}

lines2d.snssde2d <- function(x,...)
                 {
    class(x) <- "snssde2d"
    X <- x$XY[,1]
    Y <- x$XY[,2]
    lines2d(X,Y,...)
}

points2d.snssde2d <- function(x,...)
                 {
    class(x) <- "snssde2d"
    X <- x$XY[,1]
    Y <- x$XY[,2]
    points2d(X,Y,...)
}

##
## summary

summary.snssde2d <- function(object,...)
                    {
	x <- object				
    class(x) <- "snssde2d"
    summary(x$XY,...)
}

time.snssde2d <- function(x,...)
                    {
    class(x) <- "snssde2d"
    as.vector(time(x$XY))
}




################################################################################
################################################################################
##### snssde3d


snssde3d <- function(N, ...)  UseMethod("snssde3d")

snssde3d.default <- function(N =100,x0=0,y0=0,z0=0,t0=0,T=1,Dt,driftx,diffx,drifty,diffy,driftz,diffz,
                             alpha=0.5,mu=0.5,type=c("ito","str"), method=c("euler","milstein",
                             "predcorr","smilstein","taylor","heun","rk1","rk2","rk3"),...)
        {
    if (any(!is.numeric(x0) || !is.numeric(y0) || !is.numeric(z0))) stop("'x0','y0' and 'z0' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T))) stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) stop(" 'N' must be a positive integer ")
    if (any(!is.expression(driftx) || !is.expression(diffx) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions")
    if (any(!is.expression(drifty) || !is.expression(diffy) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions")
    if (any(!is.expression(driftz) || !is.expression(diffz) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions")
    if (missing(type)) type <- "ito"
    method <- match.arg(method)
    if (method =="predcorr"){
    if (any(alpha > 1 || alpha < 0)) stop("please use '0 <= alpha <= 1' ")
    if (any(mu > 1 || mu < 0))       stop("please use '0 <= mu <= 1' ")
                            }
    if (any(t0 < 0 || T < 0 || T <= t0) ) stop(" please use positive times! (0 <= t0 < T) ")
    if (missing(Dt)) {
        t <- seq(t0, T, length = N + 1)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
        T <- t[N + 1]
    }
    Dt <- (T - t0)/N    
    if (method=="euler")         {res <- .Euler3D(N,x0,y0,z0,t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,type)}
    else if (method=="predcorr") {res <- .PredCorr3D(N,x0,y0,z0,t0,T,Dt,alpha,mu,driftx,diffx,drifty,diffy,driftz,diffz,type)}
    else if (method=="milstein") {res <- .Milstein3D(N,x0,y0,z0,t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,type)}
    else if (method=="smilstein"){res <- .SMilstein3D(N,x0,y0,z0,t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,type)}
    else if (method=="taylor")   {res <- .STS3D(N,x0,y0,z0,t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,type)}
    else if (method=="heun")     {res <- .Heun3D(N,x0,y0,z0,t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,type)}
    else if (method=="rk1")      {res <- .RK3D(N,x0,y0,z0,t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,type,order=1)}
    else if (method=="rk2")      {res <- .RK3D(N,x0,y0,z0,t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,type,order=2)}
    else if (method=="rk3")      {res <- .RK3D(N,x0,y0,z0,t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,type,order=3)}
    structure(list(XYZ=res$X, driftx=driftx[[1]], diffx=diffx[[1]],drifty=drifty[[1]], diffy=diffy[[1]],driftz=driftz[[1]], 
                   diffz=diffz[[1]],type=type,method=method,x0=x0,y0=y0,z0=z0, N=N,Dt=Dt,t0=t0,T=T),class="snssde3d")
}


###

print.snssde3d <- function(x, digits=NULL, ...)
           {
    class(x) <- "snssde3d"
    if (x$method=="euler")         {sch <- "Euler scheme of order 0.5"}
    else if (x$method=="milstein") {sch <- "Milstein scheme of order 1"}
    else if (x$method=="predcorr") {sch <- "Predictor-corrector method of order 1"}
    else if (x$method=="smilstein"){sch <- "Second Milstein scheme of order 1.5"}
    else if (x$method=="taylor")   {sch <- "Ito-Taylor scheme of order 1.5"}
    else if (x$method=="heun")     {sch <- "Heun scheme of order 2"}
    else if (x$method=="rk1")      {sch <- "Runge-Kutta method of order 1"}
    else if (x$method=="rk2")      {sch <- "Runge-Kutta method of order 2"}
    else if (x$method=="rk3")      {sch <- "Runge-Kutta method of order 3"}
    if(x$type=="ito"){
    cat("Ito Sde 3D:","\n",
        "\t| dx = ", deparse(x$driftx)," * dt + ", deparse(x$diffx)," * dw1","\n", 
        "\t| dy = ", deparse(x$drifty)," * dt + ", deparse(x$diffy)," * dw2","\n",
        "\t| dz = ", deparse(x$driftz)," * dt + ", deparse(x$diffz)," * dw3","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        "\t| Initial values","\t| (x0,y0,z0) = ","(",format(x$x0,digits=digits),",",format(x$y0,digits=digits),",",format(x$z0,digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}else{
    cat("Stratonovich Sde 3D:","\n",
        "\t| dx = ", deparse(x$driftx)," * dt + ", deparse(x$diffx)," o dw1","\n", 
        "\t| dy = ", deparse(x$drifty)," * dt + ", deparse(x$diffy)," o dw2","\n",
        "\t| dz = ", deparse(x$driftz)," * dt + ", deparse(x$diffz)," o dw3","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        "\t| Initial values","\t| (x0,y0,z0) = ","(",format(x$x0,digits=digits),",",format(x$y0,digits=digits),",",format(x$z0,digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}
    invisible(x)
}

##
## Plot

plot.ssnssde3d <- function(x,plot.type = NULL,col = NULL,lty = NULL,lwd=NULL,main=NULL,...)
                 {
    class(x) <- "snssde3d"
    X <- x$XYZ
    if (is.null(lty)) {lty = 1}
    if (is.null(lwd)) {lwd = 1}
    if (is.null(main)){main=""}
    if (is.null(col)) {col =1:3}
    if (!is.null(plot.type)){
    plot(X,plot.type="single",col=col,lty=lty,lwd=lwd,main=main,...)
    legend("topright",c(expression(X[t]),expression(Y[t]),expression(Z[t])),inset = .01,col=col,lty=lty,lwd=lwd)
    }else{
    plot(X,col=col,lty=lty,lwd=lwd,main=main,...)
     }
}

plot.snssde3d <- function(x,...) plot.ssnssde3d(x,...)

lines.snssde3d <- function(x,...)
                 {
    class(x) <- "snssde3d"
    X <- x$XYZ
    for (i in 1:dim(X)[2]){
    lines(time(x),X[,i],...)}
}

points.snssde3d <- function(x,...)
                 {
    class(x) <- "snssde3d"
    X <- x$XYZ
    for (i in 1:dim(X)[2]){
    points(time(x),X[,i],...)}
}

plot3D.snssde3d <- function(x,display = c("persp", "rgl"),...)
                {
    class(x) <- "snssde3d"
    X <- x$XYZ[,1]
    Y <- x$XYZ[,2]
    Z <- x$XYZ[,3]
    plot3D(X,Y,Z,display,...)
}



##
## summary

summary.snssde3d <- function(object,...)
                    {
	x <- object				
    class(x) <- "snssde3d"
    summary(x$XYZ,...)
}

time.snssde3d <- function(x,...)
                    {
    class(x) <- "snssde3d"
    as.vector(time(x$XYZ))
}

