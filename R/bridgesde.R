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
##### bridgesde1d

bridgesde1d <- function(N, ...)  UseMethod("bridgesde1d")

bridgesde1d.default <- function(N =1000,x0=0,y=0,t0=0,T=1,Dt,drift,diffusion,
                              alpha=0.5,mu=0.5,type=c("ito","str"), method=c(
                              "euler","milstein","predcorr","smilstein","taylor",
                              "heun","rk1","rk2","rk3"),...)
                     {
    if (!is.numeric(x0) || !is.numeric(y)) stop("'x0' and 'y' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T))) stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) stop(" 'N' must be a positive integer ")
    ## if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0))  stop(" 'M' must be a positive integer ")
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
    ## bridgesde <- function()
    ##         {
    done <- FALSE
    while (!done){
    X1 <- snssde1d(N,M=1,x0,t0,T,Dt,drift,diffusion,alpha,mu,type, method,...)$X
    X2 <- apply(snssde1d(N,M=1,x0=y,t0,T,Dt,drift,diffusion,alpha,mu,type, method,...)$X,2,rev)
    G <- Inf
    if (X1[1] >= X2[1]){
           if (!all(X1 > X2))
           G <- min(which((X1-X2) <= 0)) - 1
    }else{ if (!all(X1 < X2)) 
           G <- min(which((X1-X2) >= 0)) - 1}
    if (G == 0 || G == length(X1) || G == Inf){
           stop( "A crossing has been no realized,trying again (Repeat)..." )
           done <- FALSE
    }else{ done <- TRUE }
            } 
    X <- c(X1[1:G],X2[-(1:G)])
    ## }   
    ## res <- data.frame(sapply(1:1,function(i) bridgesde()))
    ## names(res) <- paste("X",1:1,sep="")
    x <- ts(X, start = t0,end = T, deltat = Dt)
    structure(list(X=x,drift=drift[[1]], diffusion=diffusion[[1]],type=type,method=method, 
                   x0=x0,y=y, N=N,Dt=Dt,t0=t0,T=T),class="bridgesde1d")
}

###

print.bridgesde1d <- function(x, digits=NULL, ...)
           {
    class(x) <- "bridgesde1d"
    if (x$method=="euler")         {sch <- "Euler scheme of order 0.5"}
    else if (x$method=="milstein") {sch <- "Milstein scheme of order 1"}
    else if (x$method=="predcorr") {sch <- "Predictor-corrector method of order 1"}
    else if (x$method=="smilstein"){sch <- "Second Milstein scheme of order 2"}
    else if (x$method=="taylor")   {sch <- "Ito-Taylor scheme of order 1.5"}
    else if (x$method=="heun")     {sch <- "Heun scheme of order 2"}
    else if (x$method=="rk1")      {sch <- "Runge-Kutta method of order 1"}
    else if (x$method=="rk2")      {sch <- "Runge-Kutta method of order 2"}
    else if (x$method=="rk3")      {sch <- "Runge-Kutta method of order 3"}
    if(x$type=="ito"){
    cat("Ito Bridges Sde 1D:","\n",        
        "\t| dx = ", deparse(x$drift)," * dt + ", deparse(x$diffusion)," * dw","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        #"\t| Number of simulation","\t| M  = ",format(x$M,digits=digits),".","\n",
        "\t| Initial value","\t\t| x0 = ",format(x$x0,digits=digits),".","\n",
        "\t| Final value","\t\t| y = ",format(x$y,digits=digits),".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",       
        sep="")}else{
    cat("Stratonovich Bridges Sde 1D:","\n",
        "\t| dx = ", deparse(x$drift)," * dt + ", deparse(x$diffusion)," o dw","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        #"\t| Number of simulation","\t| M  = ",format(x$M,digits=digits),".","\n",
        "\t| Initial value","\t\t| x0 = ",format(x$x0,digits=digits),".","\n",
        "\t| Final value","\t\t| y = ",format(x$y,digits=digits),".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}
    invisible(x)
}

time.bridgesde1d <- function(x,...)
                    {
    class(x) <- "bridgesde1d"
    as.vector(time(x$X))
}

summary.bridgesde1d <- function(object,...)
                    {
	x <- object				
    class(x) <- "bridgesde1d"
    if (is.null(digits)) digits <- 5
    X <- x$X
    summary(X,...)
}

##
## Plot

plot.bridgesde1d <- function(x,...)
                 {
    class(x) <- "bridgesde1d"
    X <- x$X
    plot(X,...)
}

lines.bridgesde1d <- function(x,...)
                 {
    class(x) <- "bridgesde1d"
    X <- x$X
    #for (i in 1:dim(X)[2]){
    lines(time(x),X,...)
    #}
}

points.bridgesde1d <- function(x,...)
                 {
    class(x) <- "bridgesde1d"
    X <- x$X
    #for (i in 1:dim(X)[2]){
    points(time(x),X,...)
    #}
}


#####
##### bridgesde2d

bridgesde2d <- function(N, ...)  UseMethod("bridgesde2d")

bridgesde2d.default <- function(N =1000,x0=c(0,0),y=c(1,1),t0=0,T=1,Dt,driftx,diffx,drifty,
                              diffy,alpha=0.5,mu=0.5,type=c("ito","str"), method=c(
                              "euler","milstein","predcorr","smilstein","taylor",
                              "heun","rk1","rk2","rk3"),...)
                     {
    if (!is.numeric(x0) || length(x0) != 2) stop("'x0' must be numeric, and length(x0) = 2 ")
    if (!is.numeric(y)  || length(y) != 2) stop("'y' must be numeric, and length(y) = 2 ")
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
    done <- FALSE
    while (!done){
    X1 <- snssde2d(N,x0=x0[1],y0=x0[2],t0,T,Dt,driftx,diffx,drifty,diffy,alpha,mu,type, method,...)$XY
    X2 <- apply(snssde2d(N,x0=y[1],y0=y[2],t0,T,Dt,driftx,diffx,drifty,diffy,alpha,mu,type, method,...)$XY,2,rev)
    G <- rep(Inf,2)
    for (i in 1:2){
        if (X1[1,i] >= X2[1,i]){
             if (!all(X1[,i] > X2[,i]))
                 G[i] <- min(which((X1[,i]-X2[,i]) <= 0)) - 1
            }else{ if (!all(X1[,i] < X2[,i])) 
                 G[i] <- min(which((X1[,i]-X2[,i]) >= 0)) - 1
           }
                 }
         if (G[1] == 0 || G[1] == length(X1[,1]) || G[1] == Inf){
                 stop( "A crossing has been no realized,trying again (Repeat)..." )
                 done <- FALSE
         }else if (G[2] == 0 || G[2] == length(X1[,2]) || G[2] == Inf){
                 stop( "A crossing has been no realized,trying again (Repeat)..." )
                 done <- FALSE
         }else{ done <- TRUE }        
            } 
    NX <- c(X1[,1][1:G[1]],X2[,1][-(1:G[1])])
    NY <- c(X1[,2][1:G[2]],X2[,2][-(1:G[2])])
    res <- data.frame(NX,NY)
    names(res) <- paste(c("X","Y"),sep="")
    X <- ts(res, start = t0, deltat = Dt)
    structure(list(XY=X, driftx=driftx[[1]], diffx=diffx[[1]],drifty=drifty[[1]], diffy=diffy[[1]],type=type,method=method, 
                   x0=x0,y=y, N=N,Dt=Dt,t0=t0,T=T),class="bridgesde2d")
}

###

print.bridgesde2d <- function(x, digits=NULL, ...)
           {
    class(x) <- "bridgesde2d"
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
    cat("Ito Bridges Sde 2D:","\n",
        "\t| dx = ", deparse(x$driftx)," * dt + ", deparse(x$diffx)," * dw1","\n", 
        "\t| dy = ", deparse(x$drifty)," * dt + ", deparse(x$diffy)," * dw2","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        "\t| Initial values","\t| x0 = c","(",format(x$x0[1],digits=digits),",",format(x$x0[2],digits=digits),")",".","\n",
        "\t| Final values","\t\t| y  = c","(",format(x$y[1],digits=digits),",",format(x$y[2],digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}else{
    cat("Stratonovich Bridges Sde 2D:","\n",
        "\t| dx = ", deparse(x$driftx)," * dt + ", deparse(x$diffx)," o dw1","\n", 
        "\t| dy = ", deparse(x$drifty)," * dt + ", deparse(x$diffy)," o dw2","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        "\t| Initial values","\t| x0 = c","(",format(x$x0[1],digits=digits),",",format(x$x0[2],digits=digits),")",".","\n",
        "\t| Final values","\t\t| y  = c","(",format(x$y[1],digits=digits),",",format(x$y[2],digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}
    invisible(x)
}
##
## Plot

plot.bbridgesde2d <- function(x,plot.type = NULL,col = NULL,lty = NULL,lwd=NULL,main=NULL,...)
                 {
    class(x) <- "bridgesde2d"
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

plot.bridgesde2d <- function(x,...) plot.bbridgesde2d(x,...)

lines.bridgesde2d <- function(x,...)
                 {
    class(x) <- "bridgesde2d"
    X <- x$XY
    for (i in 1:dim(X)[2]){
    lines(time(x),X[,i],...)}
}

points.bridgesde2d <- function(x,...)
                 {
    class(x) <- "bridgesde2d"
    X <- x$XY
    for (i in 1:dim(X)[2]){
    points(time(x),X[,i],...)}
}

plot2d.bridgesde2d <- function(x,...)
                 {
    class(x) <- "bridgesde2d"
    X <- x$XY[,1]
    Y <- x$XY[,2]
    plot2d(X,Y,...)
    for(i in 3:4) axis(i)
}

lines2d.bridgesde2d <- function(x,...)
                 {
    class(x) <- "bridgesde2d"
    X <- x$XY[,1]
    Y <- x$XY[,2]
    lines2d(X,Y,...)
}

points2d.bridgesde2d <- function(x,...)
                 {
    class(x) <- "bridgesde2d"
    X <- x$XY[,1]
    Y <- x$XY[,2]
    points2d(X,Y,...)
}

##
## summary

summary.bridgesde2d <- function(object,...)
                    {
	x <- object				
    class(x) <- "bridgesde2d"
    summary(x$XY,...)
}

time.bridgesde2d <- function(x,...)
                    {
    class(x) <- "bridgesde2d"
    as.vector(time(x$XY))
}

#####
##### bridgesde3d

bridgesde3d <- function(N, ...)  UseMethod("bridgesde3d")

bridgesde3d.default <- function(N =1000,x0=c(0,0,0),y=c(1,-1,2),t0=0,T=1,Dt,driftx,diffx,drifty,
                              diffy,driftz,diffz,alpha=0.5,mu=0.5,type=c("ito","str"), method=c(
                              "euler","milstein","predcorr","smilstein","taylor",
                              "heun","rk1","rk2","rk3"),...)
                     {
    if (!is.numeric(x0) || length(x0) != 3) stop("'x0' must be numeric, and length(x0) = 3 ")
    if (!is.numeric(y)  || length(y) != 3) stop("'y' must be numeric, and length(y) = 3 ")
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
    done <- FALSE
    while (!done){
    X1 <- snssde3d(N,x0=x0[1],y0=x0[2],z0=x0[3],t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,alpha,mu,type, method,...)$XYZ
    X2 <- apply(snssde3d(N,x0=y[1],y0=y[2],z0=y[3],t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,alpha,mu,type, method,...)$XYZ,2,rev)
    G <- rep(Inf,3)
    for (i in 1:3){
        if (X1[1,i] >= X2[1,i]){
             if (!all(X1[,i] > X2[,i]))
                 G[i] <- min(which((X1[,i]-X2[,i]) <= 0)) - 1
            }else{ if (!all(X1[,i] < X2[,i])) 
                 G[i] <- min(which((X1[,i]-X2[,i]) >= 0)) - 1
           }
                 }
         if (G[1] == 0 || G[1] == length(X1[,1]) || G[1] == Inf ||
             G[2] == 0 || G[2] == length(X1[,2]) || G[2] == Inf ||
             G[3] == 0 || G[3] == length(X1[,3]) || G[3] == Inf){
                 stop( "A crossing has been no realized,trying again (Repeat)..." )
                 done <- FALSE
         #}else if (G[2] == 0 || G[2] == length(X1[,2]) || G[2] == Inf){
          #       stop( "A crossing has been no realized,trying again (Repeat)..." )
          #       done <- FALSE
         #}else if (G[3] == 0 || G[3] == length(X1[,3]) || G[3] == Inf){
          #       stop( "A crossing has been no realized,trying again (Repeat)..." )
          #       done <- FALSE
         }else{ done <- TRUE }        
            } 
    NX <- c(X1[,1][1:G[1]],X2[,1][-(1:G[1])])
    NY <- c(X1[,2][1:G[2]],X2[,2][-(1:G[2])])
    NZ <- c(X1[,2][1:G[3]],X2[,3][-(1:G[3])])
    res <- data.frame(NX,NY,NZ)
    names(res) <- paste(c("X","Y","Z"),sep="")
    X <- ts(res, start = t0, deltat = Dt)
    structure(list(XYZ=X, driftx=driftx[[1]], diffx=diffx[[1]],drifty=drifty[[1]], diffy=diffy[[1]],driftz=driftz[[1]],diffz=diffz[[1]],
                   type=type,method=method, x0=x0,y=y, N=N,Dt=Dt,t0=t0,T=T),class="bridgesde3d")
}

###

print.bridgesde3d <- function(x, digits=NULL, ...)
           {
    class(x) <- "bridgesde3d"
    if (x$method=="euler")         {sch <- "Euler scheme of order 0.5"}
    else if (x$method=="milstein") {sch <- "Milstein scheme of order 1"}
    else if (x$method=="predcorr") {sch <- "Predictor-corrector method of order 1"}
    else if (x$method=="smilstein"){sch <- "Second Milstein scheme of order 2"}
    else if (x$method=="taylor")   {sch <- "Ito-Taylor scheme of order 1.5"}
    else if (x$method=="heun")     {sch <- "Heun scheme of order 2"}
    else if (x$method=="rk1")      {sch <- "Runge-Kutta method of order 1"}
    else if (x$method=="rk2")      {sch <- "Runge-Kutta method of order 2"}
    else if (x$method=="rk3")      {sch <- "Runge-Kutta method of order 3"}
    if(x$type=="ito"){
    cat("Ito Bridges Sde 3D:","\n",
        "\t| dx = ", deparse(x$driftx)," * dt + ", deparse(x$diffx)," * dw1","\n", 
        "\t| dy = ", deparse(x$drifty)," * dt + ", deparse(x$diffy)," * dw2","\n",
        "\t| dz = ", deparse(x$driftz)," * dt + ", deparse(x$diffz)," * dw3","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        "\t| Initial values","\t| x0 = c","(",format(x$x0[1],digits=digits),",",format(x$x0[2],digits=digits),",",format(x$x0[3],digits=digits),")",".","\n",
        "\t| Final values","\t\t| y  = c","(",format(x$y[1],digits=digits),",",format(x$y[2],digits=digits),",",format(x$y[3],digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}else{
    cat("Stratonovich Bridges Sde 3D:","\n",
        "\t| dx = ", deparse(x$driftx)," * dt + ", deparse(x$diffx)," o dw1","\n", 
        "\t| dy = ", deparse(x$drifty)," * dt + ", deparse(x$diffy)," o dw2","\n",
        "\t| dz = ", deparse(x$driftz)," * dt + ", deparse(x$diffz)," o dw3","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        "\t| Initial values","\t| x0 = c","(",format(x$x0[1],digits=digits),",",format(x$x0[2],digits=digits),",",format(x$x0[3],digits=digits),")",".","\n",
        "\t| Final values","\t\t| y  = c","(",format(x$y[1],digits=digits),",",format(x$y[2],digits=digits),",",format(x$y[3],digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}
    invisible(x)
}

##
## summary

summary.bridgesde3d <- function(object,...)
                    {
	x <- object				
    class(x) <- "bridgesde3d"
    summary(x$XYZ,...)
}

time.bridgesde3d <- function(x,...)
                    {
    class(x) <- "bridgesde3d"
    as.vector(time(x$XYZ))
}

##
## Plot

plot.bbridgesde3d <- function(x,plot.type = NULL,col = NULL,lty = NULL,lwd=NULL,main=NULL,...)
                 {
    class(x) <- "bridgesde3d"
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

plot.bridgesde3d <- function(x,...) plot.bbridgesde3d(x,...)

lines.bridgesde3d  <- function(x,...)
                 {
    class(x) <- "bridgesde3d"
    X <- x$XYZ
    for (i in 1:dim(X)[2]){
    lines(time(x),X[,i],...)}
}

points.bridgesde3d <- function(x,...)
                 {
    class(x) <- "bridgesde3d"
    X <- x$XYZ
    for (i in 1:dim(X)[2]){
    points(time(x),X[,i],...)}
}

plot3D.bridgesde3d <- function(x,display = c("persp", "rgl"),...)
                {
    class(x) <- "bridgesde3d"
    X <- x$XYZ[,1]
    Y <- x$XYZ[,2]
    Z <- x$XYZ[,3]
    plot3D(X,Y,Z,display,...)
}


