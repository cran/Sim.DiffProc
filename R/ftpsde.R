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
##### fptsde1d

fptsde1d <- function(N, ...)  UseMethod("fptsde1d")

fptsde1d.default <- function(N =100,M=10,x0=0,t0=0,T=1,Dt,c=0.5,drift,diffusion,alpha=0.5,mu=0.5,
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
    Dt <- (T - t0)/N 
    res <- snssde1d(N,M,x0,t0,T,Dt,drift,diffusion,alpha,mu,type, method,...)
    R <- data.frame(res$X)
    if (x0 > c){
          X1  <- sapply(1:M,function(i) ifelse(length(which(R[, i] <= c)) == 0, NA,R[,i][min(which(R[,i] < c))-1] ))
          X2  <- sapply(1:M,function(i) ifelse(length(which(R[, i] <= c)) == 0, NA,R[,i][min(which(R[,i] < c))] ))
          fpt <- sapply(1:M,function(i) ifelse(length(which(R[, i] <= c)) == 0, NA,(c-X1[i])*Dt /(X2[i]-X1[i]) + time(res)[min(which(R[,i] < c))-1] ))                                    
    }else{
          X1  <- sapply(1:M,function(i) ifelse(length(which(R[, i] >= c)) == 0, NA,R[,i][min(which(R[,i] > c))-1] ))
          X2  <- sapply(1:M,function(i) ifelse(length(which(R[, i] >= c)) == 0, NA,R[,i][min(which(R[,i] > c))] ))
          fpt <- sapply(1:M,function(i) ifelse(length(which(R[, i] >= c)) == 0, NA,(c-X1[i])*Dt /(X2[i]-X1[i]) + time(res)[min(which(R[,i] > c))-1] ))                                                                       
         }
    plot(res,plot.type="single")
    if (length(which(is.na(fpt)==TRUE)) < M){
    points(fpt,rep(c,M),col = 2, pch = 16,cex=0.8)
    Axis(at = c, side=2,col = 2, labels = bquote(c))
    lines(c(0,max(fpt,na.rm=TRUE)),c(c,c),lty=2,col=2)}
    structure(list(x0=x0,c=c,tau=fpt),class="fptsde1d")
}

###

bconfint.fptsde1d <- function(x,level=0.95,...)
                 {
    class(x) <-"fptsde1d"
    x <- x$tau[!is.na(x$tau)]
    bconfint(x,level=level,...)
}

skewness.fptsde1d <- function(x,...)
                    {
    class(x) <- "fptsde1d"
    x <- x$tau[!is.na(x$tau)]
    skewness(x) 
}

kurtosis.fptsde1d <- function(x,...)
                    {
    class(x) <- "fptsde1d"
    x <- x$tau[!is.na(x$tau)]
    kurtosis(x) 
}

median.fptsde1d <- function(x,...)
                    {
    class(x) <- "fptsde1d"
    x <- x$tau[!is.na(x$tau)]
    median(x) 
}

mean.fptsde1d <- function(x,...)
                    {
    class(x) <- "fptsde1d"
    x <- x$tau[!is.na(x$tau)]
    mean(x) 
}

quantile.fptsde1d <- function(x,...)
             {
    class(x) <- "fptsde1d"
    x <- x$tau[!is.na(x$tau)]
    quantile(x,...)
}

moment.fptsde1d <- function(x,order = 2,...)
             {
    class(x) <- "fptsde1d"
    x <- x$tau[!is.na(x$tau)]
    Mom <- data.frame(do.call("cbind",lapply(1:length(order), function(j) moment(x,order=order[j]))))
    row.names(Mom) <- paste("x",sep="")
    names(Mom) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    return(Mom)
}

summary.fptsde1d <- function(object, ...)
           {
	x <- object	   
    class(x) <- "fptsde1d"
    if (x$x0 > x$c){
    cat("\n   Monte-Carlo Statistics for tau(c) = inf{t >= 0 : X(t) <= ", x$c," }","\n",
        sep="")
    }else{
    cat("\n   Monte-Carlo Statistics for tau(c) = inf{t >= 0 : X(t) >= ", x$c," }","\n",
        sep="")
    }
    x <- x$tau[!is.na(x$tau)]
    res <- data.frame(matrix(c(mean(x),var(x),median(x),quantile(x,0.25),quantile(x,0.75),
                               skewness(x),kurtosis(x),moment(x,order=2),moment(x,order=3),
                               moment(x,order=4),moment(x,order=5),bconfint(x)[1],bconfint(x)[2]),
                               ncol=1))
    row.names(res) <- paste(c("Mean","Variance","Median","First quartile","Third quartile",
                              "Skewness","Kurtosis","Moment of order 2","Moment of order 3",
                              "Moment of order 4","Moment of order 5","Bound conf Inf (95%)","Bound conf Sup (95%)"),sep="")
    names(res) <- paste(c(""),sep="")
    print(res, quote = FALSE, right = TRUE,...)
    invisible(x)
}


################################################################################
################################################################################
#####
##### fptsde2d

fptsde2d <- function(N, ...)  UseMethod("fptsde2d")


fptsde2d.default <- function(N =100,M=10,x0=0,y0=0,t0=0,T=1,Dt,c=0.5,driftx,diffx,
                     drifty,diffy,alpha=0.5,mu=0.5,type=c("ito","str"), method=c("euler",
                     "milstein","predcorr","smilstein","taylor","heun","rk1","rk2","rk3"),...)
                     {
    if (any(!is.numeric(x0) || !is.numeric(y0))) stop("'x0' and 'y0' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T))) stop(" 't0' and 'T' must be numeric")
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
    Dt <- (T - t0)/N 
    R <- lapply(1:M,function(i) data.frame(snssde2d(N,x0,y0,t0,T,Dt,driftx,diffx,drifty,diffy,alpha,mu,type, method,...)$XY))
       if (x0 > c){
             X1   <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 1] <= c)) == 0, NA,R[[i]][,1][min(which(R[[i]][,1] < c))-1] ))
             X2   <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 1] <= c)) == 0, NA,R[[i]][,1][min(which(R[[i]][,1] < c))] ))
             fptx <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 1] <= c)) == 0, NA,(c-X1[i])*Dt /(X2[i]-X1[i]) + t[min(which(R[[i]][,1] < c))-1] ))    
       }else{
             X1   <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 1] >= c)) == 0, NA,R[[i]][,1][min(which(R[[i]][,1] > c))-1] ))
             X2   <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 1] >= c)) == 0, NA,R[[i]][,1][min(which(R[[i]][,1] > c))] ))
             fptx <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 1] >= c)) == 0, NA,(c-X1[i])*Dt /(X2[i]-X1[i]) + t[min(which(R[[i]][,1] > c))-1] ))
               }
       if (y0 > c){
             Y1   <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 2] <= c)) == 0, NA,R[[i]][,2][min(which(R[[i]][,2] < c))-1] ))
             Y2   <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 2] <= c)) == 0, NA,R[[i]][,2][min(which(R[[i]][,2] < c))] ))
             fpty <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 2] <= c)) == 0, NA,(c-Y1[i])*Dt /(Y2[i]-Y1[i]) + t[min(which(R[[i]][,2] < c))-1] ))    
       }else{
             Y1   <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 2] >= c)) == 0, NA,R[[i]][,2][min(which(R[[i]][,2] > c))-1] ))
             Y2   <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 2] >= c)) == 0, NA,R[[i]][,2][min(which(R[[i]][,2] > c))] ))
             fpty <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 2] >= c)) == 0, NA,(c-Y1[i])*Dt /(Y2[i]-Y1[i]) + t[min(which(R[[i]][,2] > c))-1] ))
       }  
    structure(list(x0=x0,y0=y0,c=c,tau_x=fptx,tau_y=fpty),class="fptsde2d")
}

###

bconfint.fptsde2d <- function(x,level=0.95,...)
             {
    class(x) <- "fptsde2d"
    Bcon <- t(data.frame(do.call("cbind",lapply(4:5,function(i) bconfint(x[[i]][!is.na(x[[i]])],level=level)))))
    row.names(Bcon) <- paste(c("tau(x)","tau(y)"),sep="")
    return(Bcon)
}

skewness.fptsde2d <- function(x,...)
             {
    class(x) <- "fptsde2d"
    Skew <- data.frame(do.call("cbind",lapply(4:5,function(i) skewness(x[[i]][!is.na(x[[i]])]))),row.names = "")
    names(Skew) <-  paste(c("tau(x)","tau(y)"),sep="")
    return(Skew)
}

kurtosis.fptsde2d <- function(x,...)
             {
    class(x) <- "fptsde2d"
    Kurt <- data.frame(do.call("cbind",lapply(4:5,function(i) kurtosis(x[[i]][!is.na(x[[i]])]))),row.names = "")
    names(Kurt) <- paste(c("tau(x)","tau(y)"),sep="")
    return(Kurt)
}

median.fptsde2d <- function(x,...)
             {
    class(x) <- "fptsde2d"
    Med <- data.frame(do.call("cbind",lapply(4:5,function(i) median(x[[i]][!is.na(x[[i]])]))),row.names = "")
    names(Med) <- paste(c("tau(x)","tau(y)"),sep="")
    return(Med)
}

mean.fptsde2d <- function(x,...)
             {
    class(x) <- "fptsde2d"
    Mean <- data.frame(do.call("cbind",lapply(4:5,function(i) mean(x[[i]][!is.na(x[[i]])]))),row.names = "")
    names(Mean) <- paste(c("tau(x)","tau(y)"),sep="")
    return(Mean)
}

quantile.fptsde2d <- function(x,...)
             {
    class(x) <- "fptsde2d"
    Qun <- t(data.frame(do.call("cbind",lapply(4:5,function(i) quantile(x[[i]][!is.na(x[[i]])],...)))))
    row.names(Qun) <- paste(c("tau(x)","tau(y)"),sep="")
    return(Qun)
}

moment.fptsde2d <- function(x,order = 2,...)
             {
    class(x) <- "fptsde2d"
    Mom <- data.frame(do.call("cbind",lapply(1:length(order), function(j) sapply(4:5,function(i) moment(x[[i]][!is.na(x[[i]])],order=order[j])))))
    row.names(Mom) <- paste(c("tau(x)","tau(y)"),sep="")
    names(Mom) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    return(Mom)
}


summary.fptsde2d <- function(object, ...)
           {
	X <- object	   
    class(X) <- "fptsde2d"
    taux <- ifelse(X$x0 > X$c,paste("tau(x) = inf{t >= 0 : X(t) <= ", X$c," }"),
            paste("tau(x) = inf{t >= 0 : X(t) >= ", X$c," }"))
    tauy <- ifelse(X$y0 > X$c,paste("tau(y) = inf{t >= 0 : Y(t) <= ", X$c," }"),
            paste("tau(y) = inf{t >= 0 : Y(t) >= ", X$c," }"))
    cat("\nMonte-Carlo Statistics for :","\n",
        "\t| ",taux,"\n",
        "\t| ",tauy,"\n\n",
       sep="")
    x <- X$tau_x[!is.na(X$tau_x)]
    y <- X$tau_y[!is.na(X$tau_y)]
    res <- data.frame(matrix(c(mean(x),var(x),median(x),quantile(x,0.25),quantile(x,0.75),
                               skewness(x),kurtosis(x),moment(x,order=2),moment(x,order=3),
                               moment(x,order=4),moment(x,order=5),bconfint(x)[1],bconfint(x)[2],
                               mean(y),var(y),median(y),quantile(y,0.25),quantile(y,0.75),
                               skewness(y),kurtosis(y),moment(y,order=2),moment(y,order=3),
                               moment(y,order=4),moment(y,order=5),bconfint(y)[1],bconfint(y)[2]),
                               ncol=2))
    row.names(res) <- paste(c("Mean","Variance","Median","First quartile","Third quartile",
                              "Skewness","Kurtosis","Moment of order 2","Moment of order 3",
                              "Moment of order 4","Moment of order 5","Bound conf Inf (95%)","Bound conf Sup (95%)"),sep="")
    names(res) <- paste(c("tau(x)","tau(y)"),sep="")
    print(res, quote = FALSE, right = TRUE,...)
    invisible(X)
}


################################################################################
################################################################################
#####
##### fptsde3d

fptsde3d <- function(N, ...)  UseMethod("fptsde3d")

fptsde3d.default <- function(N =100,M=10,x0=0,y0=0,z0=0,t0=0,T=1,Dt,c=0.5,driftx,diffx,drifty,diffy,
                     driftz,diffz,alpha=0.5,mu=0.5,type=c("ito","str"), method=c("euler",
                     "milstein","predcorr","smilstein","taylor","heun","rk1","rk2","rk3"),...)
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
    Dt <- (T - t0)/N 
    R <- lapply(1:M,function(i) data.frame(snssde3d(N,x0,y0,z0,t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,alpha,mu,type, method,...)$XYZ))
    ##res <- lapply(1:M,function(i) snssde3d(N,x0,y0,z0,t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,alpha,mu,type, method,...))
    ##R <- lapply(1:M,function(i) data.frame(res[[i]]$XYZ))
       if (x0 > c){
             X1   <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 1] <= c)) == 0, NA,R[[i]][,1][min(which(R[[i]][,1] < c))-1] ))
             X2   <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 1] <= c)) == 0, NA,R[[i]][,1][min(which(R[[i]][,1] < c))] ))
             fptx <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 1] <= c)) == 0, NA,(c-X1[i])*Dt /(X2[i]-X1[i]) + t[min(which(R[[i]][,1] < c))-1] ))    
       }else{
             X1   <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 1] >= c)) == 0, NA,R[[i]][,1][min(which(R[[i]][,1] > c))-1] ))
             X2   <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 1] >= c)) == 0, NA,R[[i]][,1][min(which(R[[i]][,1] > c))] ))
             fptx <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 1] >= c)) == 0, NA,(c-X1[i])*Dt /(X2[i]-X1[i]) + t[min(which(R[[i]][,1] > c))-1] ))
               }
       if (y0 > c){
             Y1   <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 2] <= c)) == 0, NA,R[[i]][,2][min(which(R[[i]][,2] < c))-1] ))
             Y2   <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 2] <= c)) == 0, NA,R[[i]][,2][min(which(R[[i]][,2] < c))] ))
             fpty <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 2] <= c)) == 0, NA,(c-Y1[i])*Dt /(Y2[i]-Y1[i]) + t[min(which(R[[i]][,2] < c))-1] ))    
       }else{
             Y1   <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 2] >= c)) == 0, NA,R[[i]][,2][min(which(R[[i]][,2] > c))-1] ))
             Y2   <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 2] >= c)) == 0, NA,R[[i]][,2][min(which(R[[i]][,2] > c))] ))
             fpty <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 2] >= c)) == 0, NA,(c-Y1[i])*Dt /(Y2[i]-Y1[i]) + t[min(which(R[[i]][,2] > c))-1] ))
                  }  
       if (z0 > c){
             Z1   <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 3] <= c)) == 0, NA,R[[i]][,3][min(which(R[[i]][,3] < c))-1] ))
             Z2   <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 3] <= c)) == 0, NA,R[[i]][,3][min(which(R[[i]][,3] < c))] ))
             fptz <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 3] <= c)) == 0, NA,(c-Z1[i])*Dt /(Z2[i]-Z1[i]) + t[min(which(R[[i]][,3] < c))-1] ))    
       }else{
             Z1   <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 3] >= c)) == 0, NA,R[[i]][,3][min(which(R[[i]][,3] > c))-1] ))
             Z2   <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 3] >= c)) == 0, NA,R[[i]][,3][min(which(R[[i]][,3] > c))] ))
             fptz <- sapply(1:M,function(i) ifelse(length(which(R[[i]][, 3] >= c)) == 0, NA,(c-Z1[i])*Dt /(Z2[i]-Z1[i]) + t[min(which(R[[i]][,3] > c))-1] ))
                  }
    structure(list(x0=x0,y0=y0,z0=z0,c=c,tau_x=fptx,tau_y=fpty,tau_z=fptz),class="fptsde3d")
}

###

bconfint.fptsde3d <- function(x,level=0.95,...)
             {
    class(x) <- "fptsde3d"
    Bcon <- t(data.frame(do.call("cbind",lapply(5:7,function(i) bconfint(x[[i]][!is.na(x[[i]])],level=level)))))
    row.names(Bcon) <- paste(c("tau(x)","tau(y)","tau(z)"),sep="")
    return(Bcon)
}

skewness.fptsde3d <- function(x,...)
             {
    class(x) <- "fptsde3d"
    Skew <- data.frame(do.call("cbind",lapply(5:7,function(i) skewness(x[[i]][!is.na(x[[i]])]))),row.names = "")
    names(Skew) <- paste(c("tau(x)","tau(y)","tau(z)"),sep="")
    return(Skew)
}

kurtosis.fptsde3d <- function(x,...)
             {
    class(x) <- "fptsde3d"
    Kurt <- data.frame(do.call("cbind",lapply(5:7,function(i) kurtosis(x[[i]][!is.na(x[[i]])]))),row.names = "")
    names(Kurt) <-paste(c("tau(x)","tau(y)","tau(z)"),sep="")
    return(Kurt)
}

median.fptsde3d <- function(x,...)
             {
    class(x) <- "fptsde3d"
    Med <- data.frame(do.call("cbind",lapply(5:7,function(i) median(x[[i]][!is.na(x[[i]])]))),row.names = "")
    names(Med) <-paste(c("tau(x)","tau(y)","tau(z)"),sep="")
    return(Med)
}

mean.fptsde3d <- function(x,...)
             {
    class(x) <- "fptsde3d"
    Mean <- data.frame(do.call("cbind",lapply(5:7,function(i) mean(x[[i]][!is.na(x[[i]])]))),row.names = "")
    names(Mean) <-paste(c("tau(x)","tau(y)","tau(z)"),sep="")
    return(Mean)
}

quantile.fptsde3d <- function(x,...)
             {
    class(x) <- "fptsde3d"
    Qun <- t(data.frame(do.call("cbind",lapply(5:7,function(i) quantile(x[[i]][!is.na(x[[i]])],...)))))
    row.names(Qun) <-paste(c("tau(x)","tau(y)","tau(z)"),sep="")
    return(Qun)
}

moment.fptsde3d <- function(x,order = 2,...)
             {
    class(x) <- "fptsde3d"
    Mom <- data.frame(do.call("cbind",lapply(1:length(order), function(j) sapply(5:7,function(i) moment(x[[i]][!is.na(x[[i]])],order=order[j])))))
    row.names(Mom) <-paste(c("tau(x)","tau(y)","tau(z)"),sep="")
    names(Mom) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    return(Mom)
}

summary.fptsde3d <- function(object, ...)
           {
	X <- object	   
    class(X) <- "fptsde3d"
    taux <- ifelse(X$x0 > X$c,paste("tau(x) = inf{t >= 0 : X(t) <= ", X$c," }"),
            paste("tau(x) = inf{t >= 0 : X(t) >= ", X$c," }"))
    tauy <- ifelse(X$y0 > X$c,paste("tau(y) = inf{t >= 0 : Y(t) <= ", X$c," }"),
            paste("tau(y) = inf{t >= 0 : Y(t) >= ", X$c," }"))
    tauz <- ifelse(X$z0 > X$c,paste("tau(z) = inf{t >= 0 : Z(t) <= ", X$c," }"),
            paste("tau(z) = inf{t >= 0 : Z(t) >= ", X$c," }"))
    cat("\nMonte-Carlo Statistics for :","\n",
        "\t| ",taux,"\n",
        "\t| ",tauy,"\n",
        "\t| ",tauz,"\n\n",
       sep="")
    x <- X$tau_x[!is.na(X$tau_x)]
    y <- X$tau_y[!is.na(X$tau_y)]
    z <- X$tau_z[!is.na(X$tau_z)]
    res <- data.frame(matrix(c(mean(x),var(x),median(x),quantile(x,0.25),quantile(x,0.75),
                               skewness(x),kurtosis(x),moment(x,order=2),moment(x,order=3),
                               moment(x,order=4),moment(x,order=5),bconfint(x)[1],bconfint(x)[2],
                               mean(y),var(y),median(y),quantile(y,0.25),quantile(y,0.75),
                               skewness(y),kurtosis(y),moment(y,order=2),moment(y,order=3),
                               moment(y,order=4),moment(y,order=5),bconfint(y)[1],bconfint(y)[2],
                               mean(z),var(z),median(z),quantile(z,0.25),quantile(z,0.75),
                               skewness(z),kurtosis(z),moment(z,order=2),moment(z,order=3),
                               moment(z,order=4),moment(z,order=5),bconfint(z)[1],bconfint(z)[2]),
                               ncol=3))
    row.names(res) <- paste(c("Mean","Variance","Median","First quartile","Third quartile",
                              "Skewness","Kurtosis","Moment of order 2","Moment of order 3",
                              "Moment of order 4","Moment of order 5","Bound conf Inf (95%)","Bound conf Sup (95%)"),sep="")
    names(res) <- paste(c("tau(x)","tau(y)","tau(z)"),sep="")
    print(res, quote = FALSE, right = TRUE,...)
    invisible(X)
}

