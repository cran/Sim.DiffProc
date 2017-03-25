## Fri Mar 24 23:26:35 2017
## Original file Copyright © 2017 A.C. Guidoum, K. Boukhetala
## This file is part of the R package Sim.DiffProc
## Department of Probabilities & Statistics
## Faculty of Mathematics
## University of Science and Technology Houari Boumediene
## BP 32 El-Alia, U.S.T.H.B, Algiers
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
##### dsde1d


dsde1d <- function(object, ...)  UseMethod("dsde1d")

dsde1d.default <- function(object,at,...)
                     {
    class(object) <- "snssde1d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)}else{X = object$X}
    x   <- as.vector(X[which(time(object)==at),])
    if (length(x) == 0){
	if (as.numeric(object$M==1)){ F   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
               F   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    structure(list(ech=x,res=density.default(x,na.rm = TRUE,...),at=at),class="dsde1d")
}

print.dsde1d <- function(x, digits=NULL, ...)
           {
    class(x) <- "dsde1d"
    cat("\n Density of X(t)|X(0)=x0 at time t = ",as.numeric(x$at),"\n",
        sep="")
    cat(
	"\nData: ", x$res$data.name, " (", x$res$n, " obs.);",
	"\tBandwidth 'bw' = ", formatC(x$res$bw, digits = digits), "\n\n", sep = "")
    out <- as.data.frame(x$res[c("x","y")])
    names(out) <- c("x","f(x)")
    print(summary(out, digits = digits), ...)
    invisible(x)
}

##

plot.dsde1d <- function(x,...) .plot.dsde1d(x,...)


#####
##### dsde2d

dsde2d <- function(object, ...)  UseMethod("dsde2d")

dsde2d.default <- function(object,pdf=c("Joint","Marginal"),at,...)
                     {
    class(object) <- "snssde2d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	                     Y = matrix(object$Y,nrow=length(object$Y),ncol=1)}else{
				  X = object$X
				  Y = object$Y}
    x   <- as.vector(X[which(as.vector(time(object))==at),])
    y   <- as.vector(Y[which(time(object)==at),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
                      Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
                      x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
                      Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
                      y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    pdf <- match.arg(pdf)
    if (pdf=="Marginal"){
    structure(list(ech=data.frame(x,y),resx=density.default(x,na.rm = TRUE,...),resy=density.default(y,na.rm = TRUE,...),at=at,pdf=pdf,SDE=object),class="dsde2d")
    }else{
    structure(list(ech=data.frame(x,y),res=MASS::kde2d(x, y,n = 100, ...),at=at,pdf=pdf,SDE=object),class="dsde2d")
    }
}

print.dsde2d <- function(x, digits=NULL, ...)
           {
    class(x) <- "dsde2d"
    if (x$pdf=="Joint") {
    cat("\n Joint density of (X(t),Y(t)|X(0)=x0,Y(0)=y0) at time t = ",as.numeric(x$at),"\n",
        sep="")
    cat(
	"\nData: (x,y)", " (2 x ", dim(x$ech)[1], " obs.);", "\n\n", sep = "")
	##"\tBandwidth 'bw' = ", formatC(MASS::bandwidth.nrd(x$res$x), digits = digits),"~~", formatC(MASS::bandwidth.nrd(x$res$y), digits = digits), "\n\n", sep = "")
    out3 <- list(x=x$res$x,y=x$res$y,z=as.vector(x$res$z))
    names(out3) <- c("x","y","f(x,y)")
    print(summary.data.frame(out3, digits = digits), ...)}else{
    cat("\n Marginal density of X(t)|X(0)=x0 at time t = ",as.numeric(x$at),"\n",
        sep="")
    cat(
	"\nData: ", x$resx$data.name, " (", x$resx$n, " obs.);",
	"\tBandwidth 'bw' = ", formatC(x$resx$bw, digits = digits), "\n\n", sep = "")
    out1 <- as.data.frame(x$resx[c("x","y")])
    names(out1) <- c("x","f(x)")
    print(summary(out1, digits = digits), ...)

    cat("\n Marginal density of Y(t)|Y(0)=y0 at time t = ",as.numeric(x$at),"\n",
        sep="")
    cat(
	"\nData: ", x$resy$data.name, " (", x$resy$n, " obs.);",
	"\tBandwidth 'bw' = ", formatC(x$resy$bw, digits = digits), "\n\n", sep = "")
    out2 <- as.data.frame(x$resy[c("x","y")])
    names(out2) <- c("y","f(y)")
    print(summary(out2, digits = digits), ...)
   }
    invisible(x)
}

plot.dsde2d <- function(x,display=c("persp","rgl","image","contour"),hist=FALSE,...) .plot.dsde2d(x,display,hist,...)

#####
##### dsde3d

dsde3d <- function(object, ...)  UseMethod("dsde3d")

dsde3d.default <- function(object,pdf=c("Joint","Marginal"),at,...)
                     {
    class(object) <- "snssde3d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  stop( " please use 't0 <= at <= T'")		
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)
	              Y = matrix(object$Y,nrow=length(object$Y),ncol=1)
				  Z = matrix(object$Z,nrow=length(object$Z),ncol=1)}else{
				  X = object$X
				  Y = object$Y
				  Z = object$Z}
    x   <- as.vector(X[which(time(object)==at),])
    y   <- as.vector(Y[which(time(object)==at),])
    z   <- as.vector(Z[which(time(object)==at),])
    if (length(x) == 0){
	if (as.numeric(object$M)==1){ Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
               Fx   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
               x   <- sapply(1:length(Fx),function(i) Fx[[i]](at)) 
    }
    if (length(y) == 0){
	if (as.numeric(object$M)==1){ Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y))}else{
               Fy   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Y[,i]))}
               y   <- sapply(1:length(Fy),function(i) Fy[[i]](at)) 
    }
    if (length(z) == 0){
	if (as.numeric(object$M)==1){ Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z))}else{
               Fz   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$Z[,i]))}
               z   <- sapply(1:length(Fz),function(i) Fz[[i]](at)) 
    }    
    pdf <- match.arg(pdf)
    if (pdf=="Marginal"){
    structure(list(ech=data.frame(x,y,z),resx=density.default(x,na.rm = TRUE,...),resy=density.default(y,na.rm = TRUE,...),resz=density.default(z,na.rm = TRUE,...),pdf=pdf,at=at),class="dsde3d")
	}else{
    structure(list(ech=data.frame(x,y,z),res=ks::kde(data.frame(x,y,z), ...),at=at,pdf=pdf,SDE=object),class="dsde3d")
    }
}

print.dsde3d <- function(x, digits=NULL, ...)
           {
    class(x) <- "dsde3d"
    if (x$pdf=="Joint") {
    cat("\n Joint density for (X(t),Y(t),Z(t)|X(0)=x0,Y(0)=y0,Z(0)=z0) at time t = ",as.numeric(x$at),"\n",
        sep="")
    cat(
	"\nData: (x,y,z)", " (3 x ", dim(x$ech)[1], " obs.);", 
	"\tBandwidth 'bw' = c(", formatC(x$res$H[1,1], digits = 3),",",formatC(x$res$H[2,2], digits = 3),",",formatC(x$res$H[3,3], digits = 3), ")\n\n", sep = "")
    out3 <- list(x=x$res$eval.points[[1]],y=x$res$eval.points[[2]],z=x$res$eval.points[[3]],d=as.vector(x$res$estimate))
    names(out3) <- c("x","y","z","f(x,y,z)")
    print(summary.data.frame(out3, digits = digits), ...)}else{
    cat("\n Marginal density for X(t)|X(0)=x0 at time t = ",as.numeric(x$at),"\n",
        sep="")
    cat(
	"\nData: ", x$resx$data.name, " (", x$resx$n, " obs.);",
	"\tBandwidth 'bw' = ", formatC(x$resx$bw, digits = digits), "\n\n", sep = "")
    out1 <- as.data.frame(x$resx[c("x","y")])
    names(out1) <- c("x","f(x)")
    print(summary(out1, digits = digits), ...)

    cat("\n Marginal density for Y(t)|Y(0)=y0 at time t = ",as.numeric(x$at),"\n",
        sep="")
    cat(
	"\nData: ", x$resy$data.name, " (", x$resy$n, " obs.);",
	"\tBandwidth 'bw' = ", formatC(x$resy$bw, digits = digits), "\n\n", sep = "")
    out2 <- as.data.frame(x$resy[c("x","y")])
    names(out2) <- c("y","f(y)")
    print(summary(out2, digits = digits), ...)

    cat("\n Marginal density for Z(t)|Z(0)=z0 at time t = ",as.numeric(x$at),"\n",
        sep="")
    cat(
	"\nData: ", x$resz$data.name, " (", x$resz$n, " obs.);",
	"\tBandwidth 'bw' = ", formatC(x$resz$bw, digits = digits), "\n\n", sep = "")
    out3 <- as.data.frame(x$resz[c("x","y")])
    names(out3) <- c("z","f(z)")
    print(summary(out3, digits = digits), ...)
	}
    invisible(x)
}

plot.dsde3d <- function(x,display=c("rgl","persp"),hist=FALSE,...) .plot.dsde3d(x,display,hist,...)