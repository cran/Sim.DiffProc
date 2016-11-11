## Mon Oct 03 16:04:10 2016
## Original file Copyright © 2016 A.C. Guidoum, K. Boukhetala
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
##### rsde1d



rsde1d <- function(object, ...)  UseMethod("rsde1d")

rsde1d.default <- function(object,at,...)
                     {
    class(object) <- "snssde1d"
    if (missing(at)) {at = as.numeric(object$T)}
    if (any(as.numeric(object$T) < at || as.numeric(object$t0) > at) )  stop( " please use 't0 <= at <= T'")
    if (as.numeric(object$M) == 1){  X = matrix(object$X,nrow=length(object$X),ncol=1)}else{X = object$X}
    xx   <- as.vector(X[which(time(object)==at),])
    if (length(xx) == 0){
	if (as.numeric(object$M==1)){ F   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X))}else{
               F   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),object$X[,i]))}
               xx   <- sapply(1:length(F),function(i) F[[i]](at)) 
    }
    return(xx)
}

###


#####
##### rsde2d

rsde2d <- function(object, ...)  UseMethod("rsde2d")

rsde2d.default <- function(object,at,...)
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
    return(data.frame(x,y))
}
##
##


#####
##### rsde3d

rsde3d <- function(object, ...)  UseMethod("rsde3d")

rsde3d.default <- function(object,at,...)
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
    return(data.frame(x,y,z))
}

###

