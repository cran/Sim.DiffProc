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
##### fptsde1d

fptsde1d <- function(object, ...)  UseMethod("fptsde1d")

fptsde1d.default <- function(object,boundary,...)
                     {
    class(object) <- "snssde1d"
    if (any(!is.expression(boundary))) stop(" must be expressions of a constant or time-dependent boundary ")
    #if (length(all.vars(boundary)) > 1 )  stop("boundary must depend on 't' or constant")
    #if (length(all.vars(boundary)) == 1 && all.vars(boundary) != "t" ) stop("boundary must depend on 't' or constant")
    R   <- data.frame(object$X)
    Bn  <- function(t)  eval(boundary)
    F   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),as.vector(R[,i])) )
      if (as.numeric(object$x0) > Bn(as.numeric(object$t0))){
           v1  <- sapply(1:length(F),function(i) ifelse(length(which(R[,i] <= Bn(time(object))))==0, NA,min(which(R[,i] <= Bn(time(object))))))
           fpt <- sapply(1:length(F),function(i) ifelse(is.na(v1[i]),NA,uniroot(f = function(x) (F[[i]](x)- Bn(x)),lower=time(object)[v1[i]-1],upper=time(object)[v1[i]],tol= .Machine$double.eps)$root))
      }else if (as.numeric(object$x0) < Bn(as.numeric(object$t0))){
           v2  <- sapply(1:length(F),function(i) ifelse(length(which(R[,i] <= Bn(time(object))))==as.numeric(object$N)+1, NA,min(which(R[,i] >= Bn(time(object))))))
           fpt <- sapply(1:length(F),function(i) ifelse(is.na(v2[i]),NA,uniroot(f = function(x) (F[[i]](x)- Bn(x)),lower=time(object)[v2[i]-1],upper=time(object)[v2[i]],tol= .Machine$double.eps)$root))
      }else{
           fpt <- rep(0,as.numeric(object$M))
                     }
    structure(list(SDE=object,boundary=boundary[[1]],fpt=fpt),class="fptsde1d")
}

###

summary.fptsde1d <- function(object, digits=5,...)
           {
    x <- object	   
    class(x) <- "fptsde1d"
    S <- function(t) eval(x$boundary)
    if (as.numeric(x$SDE$x0) > S(as.numeric(x$SDE$t0))){
    cat("\n\t Monte-Carlo Statistics for the F.P.T of X(t)","\n", "| T(S,X) = inf{t >= ",x$SDE$t0 ," : X(t) <= ", deparse(x$boundary),"}","\n",
        sep="")
    }else{
    cat("\n\t Monte-Carlo Statistics for the F.P.T of X(t)","\n", "| T(S,X) = inf{t >= ",x$SDE$t0 ," : X(t) >= ", deparse(x$boundary),"}","\n",
        sep="")
    }
    res <- as.data.frame(matrix(c(format(length(which(is.na(x$fpt))),digits=1),round(digits=digits,mean(x$fpt,na.rm = TRUE)),round(digits=digits,var(x$fpt,na.rm = TRUE)),
                               round(digits=digits,median(x$fpt,na.rm = TRUE)),round(digits=digits,quantile(x$fpt,0.25,na.rm = TRUE)),
							   round(digits=digits,quantile(x$fpt,0.75,na.rm = TRUE)),round(digits=digits,skewness(x$fpt)),round(digits=digits,kurtosis(x$fpt)),
							   round(digits=digits,moment(x$fpt,order=3)),round(digits=digits,moment(x$fpt,order=4)),
							   round(digits=digits,moment(x$fpt,order=5)),round(digits=digits,bconfint(x$fpt)[1]),round(digits=digits,bconfint(x$fpt)[2])),
                               ncol=1))
	dimnames(res) <- list(c("NA's","Mean","Variance","Median","First quartile","Third quartile",
                              "Skewness","Kurtosis","Moment of order 3",
                              "Moment of order 4","Moment of order 5","Int.conf Inf (95%)","Int.conf Sup (95%)"),c("T(S,X)"))
    print(res, quote = FALSE, right = TRUE,...)
    invisible(x)
}

plot.fptsde1d <- function(x,...) .plot.fptsde1d(x,...)

################################################################################
################################################################################
#####
##### fptsde2d

fptsde2d <- function(object, ...)  UseMethod("fptsde2d")


fptsde2d.default <- function(object,boundary,...)
                     {
    class(object) <- "snssde2d"
    if (any(!is.expression(boundary))) stop(" must be expression of a constant or time-dependent boundary ")
    R1   <- data.frame(object$X)
    R2   <- data.frame(object$Y)
    Bn  <- function(t)  eval(boundary) 
    F1   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),as.vector(R1[,i])) )
    F2   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),as.vector(R2[,i])) )
    if (as.numeric(object$x0) > Bn(as.numeric(object$t0))){
           v11  <- sapply(1:length(F1),function(i) ifelse(length(which(R1[,i] <= Bn(time(object))))==0, NA,min(which(R1[,i] <= Bn(time(object))))))
           fptx <- sapply(1:length(F1),function(i) ifelse(is.na(v11[i]),NA,uniroot(f = function(x) (F1[[i]](x)- Bn(x)),lower=time(object)[v11[i]-1],upper=time(object)[v11[i]],tol= .Machine$double.eps)$root))
      }else if (as.numeric(object$x0) < Bn(as.numeric(object$t0))){
           v12  <- sapply(1:length(F1),function(i) ifelse(length(which(R1[,i] <= Bn(time(object))))==as.numeric(object$N)+1, NA,min(which(R1[,i] >= Bn(time(object))))))
           fptx <- sapply(1:length(F1),function(i) ifelse(is.na(v12[i]),NA,uniroot(f = function(x) (F1[[i]](x)- Bn(x)),lower=time(object)[v12[i]-1],upper=time(object)[v12[i]],tol= .Machine$double.eps)$root))
      }else{
           fptx <- rep(0,as.numeric(object$M))
                     }
    if (as.numeric(object$y0) > Bn(as.numeric(object$t0))){
           v21  <- sapply(1:length(F2),function(i) ifelse(length(which(R2[,i] <= Bn(time(object))))==0, NA,min(which(R2[,i] <= Bn(time(object))))))
           fpty <- sapply(1:length(F2),function(i) ifelse(is.na(v21[i]),NA,uniroot(f = function(x) (F2[[i]](x)- Bn(x)),lower=time(object)[v21[i]-1],upper=time(object)[v21[i]],tol= .Machine$double.eps)$root))
      }else if (as.numeric(object$y0) < Bn(as.numeric(object$t0))){
           v22  <- sapply(1:length(F2),function(i) ifelse(length(which(R2[,i] <= Bn(time(object))))==as.numeric(object$N)+1, NA,min(which(R2[,i] >= Bn(time(object))))))
           fpty <- sapply(1:length(F2),function(i) ifelse(is.na(v22[i]),NA,uniroot(f = function(x) (F2[[i]](x)- Bn(x)),lower=time(object)[v22[i]-1],upper=time(object)[v22[i]],tol= .Machine$double.eps)$root))
      }else{
           fpty <- rep(0,as.numeric(object$M))
                     }
    structure(list(SDE=object,boundary=boundary[[1]],fptx=fptx,fpty=fpty),class="fptsde2d")
}
###

summary.fptsde2d <- function(object,digits=5, ...)
           {  
    class(object) <- "fptsde2d"
    S <- function(t) eval(object$boundary)
    if (as.numeric(object$SDE$x0) > S(as.numeric(object$SDE$t0))){
    fpt_x <- paste("T(S,X) = inf{t >= ",as.numeric(object$SDE$t0) ," : X(t) <= ", deparse(object$boundary),"}")
    }else{
    fpt_x <- paste("T(S,X) = inf{t >= ",as.numeric(object$SDE$t0) ," : X(t) >= ", deparse(object$boundary),"}")
    }
    if (as.numeric(object$SDE$y0) > S(as.numeric(object$SDE$t0))){
    fpt_y <- paste("T(S,Y) = inf{t >= ",as.numeric(object$SDE$t0) ," : Y(t) <= ", deparse(object$boundary),"}")
    }else{
    fpt_y <- paste("T(S,Y) = inf{t >= ",as.numeric(object$SDE$t0) ," : Y(t) >= ", deparse(object$boundary),"}")
    }
    cat("\n\t Monte-Carlo Statistics for the F.P.T of (X(t),Y(t))","\n",
        "| ",fpt_x,"\n",
        "| ",fpt_y,"\n",
       sep="")
    x <- object$fptx##[!is.na(X$fptx)]
    y <- object$fpty##[!is.na(X$fpty)]
    res <- as.data.frame(matrix(c(format(length(which(is.na(x))),digits=1),round(digits=digits,mean(x,na.rm = TRUE)),round(digits=digits,var(x,na.rm = TRUE)),round(digits=digits,median(x,na.rm = TRUE)),
                               round(digits=digits,quantile(x,0.25,na.rm = TRUE)),round(digits=digits,quantile(x,0.75,na.rm = TRUE)),
                               round(digits=digits,skewness(x)),round(digits=digits,kurtosis(x)),round(digits=digits,moment(x,order=3)),
                               round(digits=digits,moment(x,order=4)),round(digits=digits,moment(x,order=5)),round(digits=digits,bconfint(x)[1]),round(digits=digits,bconfint(x)[2]),
                               format(length(which(is.na(y))),digits=1),round(digits=digits,mean(y,na.rm = TRUE)),round(digits=digits,var(y,na.rm = TRUE)),round(digits=digits,median(y,na.rm = TRUE)),
                               round(digits=digits,quantile(y,0.25,na.rm = TRUE)),round(digits=digits,quantile(y,0.75,na.rm = TRUE)),
                               round(digits=digits,skewness(y)),round(digits=digits,kurtosis(y)),round(digits=digits,moment(y,order=3)),
                               round(digits=digits,moment(y,order=4)),round(digits=digits,moment(y,order=5)),round(digits=digits,bconfint(y)[1]),round(digits=digits,bconfint(y)[2])),
                               ncol=2))
	dimnames(res) <- list(c("NA's","Mean","Variance","Median","First quartile","Third quartile",
                              "Skewness","Kurtosis","Moment of order 3",
                              "Moment of order 4","Moment of order 5","Int.conf Inf (95%)","Int.conf Sup (95%)"),c("T(S,X)","T(S,Y)"))
    print(res, quote = FALSE, right = TRUE,...)
    invisible(object)
}


plot.fptsde2d <- function(x,...) .plot.fptsde2d(x,...)
	

################################################################################
################################################################################
#####
##### fptsde3d

fptsde3d <- function(object, ...)  UseMethod("fptsde3d")

fptsde3d.default <- function(object,boundary,...)
                     {
    class(object) <- "snssde3d"
    if (any(!is.expression(boundary))) stop(" must be expression of a constant or time-dependent boundary ")
    R1   <- data.frame(object$X)
    R2   <- data.frame(object$Y)
    R3   <- data.frame(object$Z)	
    Bn  <- function(t)  eval(boundary)
    F1   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),as.vector(R1[,i])) )
    F2   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),as.vector(R2[,i])) )
    F3   <- lapply(1:as.numeric(object$M),function(i) approxfun(time(object),as.vector(R3[,i])) )
    if (as.numeric(object$x0) > Bn(as.numeric(object$t0))){
           v11  <- sapply(1:length(F1),function(i) ifelse(length(which(R1[,i] <= Bn(time(object))))==0, NA,min(which(R1[,i] <= Bn(time(object))))))
           fptx <- sapply(1:length(F1),function(i) ifelse(is.na(v11[i]),NA,uniroot(f = function(x) (F1[[i]](x)- Bn(x)),lower=time(object)[v11[i]-1],upper=time(object)[v11[i]],tol= .Machine$double.eps)$root))
      }else if (as.numeric(object$x0) < Bn(as.numeric(object$t0))){
           v12  <- sapply(1:length(F1),function(i) ifelse(length(which(R1[,i] <= Bn(time(object))))==as.numeric(object$N)+1, NA,min(which(R1[,i] >= Bn(time(object))))))
           fptx <- sapply(1:length(F1),function(i) ifelse(is.na(v12[i]),NA,uniroot(f = function(x) (F1[[i]](x)- Bn(x)),lower=time(object)[v12[i]-1],upper=time(object)[v12[i]],tol= .Machine$double.eps)$root))
      }else{
           fptx <- rep(0,as.numeric(object$M))
                     }
    if (as.numeric(object$y0) > Bn(as.numeric(object$t0))){
           v21  <- sapply(1:length(F2),function(i) ifelse(length(which(R2[,i] <= Bn(time(object))))==0, NA,min(which(R2[,i] <= Bn(time(object))))))
           fpty <- sapply(1:length(F2),function(i) ifelse(is.na(v21[i]),NA,uniroot(f = function(x) (F2[[i]](x)- Bn(x)),lower=time(object)[v21[i]-1],upper=time(object)[v21[i]],tol= .Machine$double.eps)$root))
      }else if (as.numeric(object$y0) < Bn(as.numeric(object$t0))){
           v22  <- sapply(1:length(F2),function(i) ifelse(length(which(R2[,i] <= Bn(time(object))))==as.numeric(object$N)+1, NA,min(which(R2[,i] >= Bn(time(object))))))
           fpty <- sapply(1:length(F2),function(i) ifelse(is.na(v22[i]),NA,uniroot(f = function(x) (F2[[i]](x)- Bn(x)),lower=time(object)[v22[i]-1],upper=time(object)[v22[i]],tol= .Machine$double.eps)$root))
      }else{
           fpty <- rep(0,as.numeric(object$M))
                     }
    if (as.numeric(object$z0) > Bn(as.numeric(object$t0))){
           v31  <- sapply(1:length(F3),function(i) ifelse(length(which(R3[,i] <= Bn(time(object))))==0, NA,min(which(R3[,i] <= Bn(time(object))))))
           fptz <- sapply(1:length(F3),function(i) ifelse(is.na(v31[i]),NA,uniroot(f = function(x) (F3[[i]](x)- Bn(x)),lower=time(object)[v31[i]-1],upper=time(object)[v31[i]],tol= .Machine$double.eps)$root))
      }else if (object$z0 < Bn(as.numeric(object$t0))){
           v32  <- sapply(1:length(F3),function(i) ifelse(length(which(R3[,i] <= Bn(time(object))))==as.numeric(object$N)+1, NA,min(which(R3[,i] >= Bn(time(object))))))
           fptz <- sapply(1:length(F3),function(i) ifelse(is.na(v32[i]),NA,uniroot(f = function(x) (F3[[i]](x)- Bn(x)),lower=time(object)[v32[i]-1],upper=time(object)[v32[i]],tol= .Machine$double.eps)$root))
      }else{
           fptz <- rep(0,as.numeric(object$M))
                     }
    structure(list(SDE=object,boundary=boundary[[1]],fptx=fptx,fpty=fpty,fptz=fptz),class="fptsde3d")
}

###

summary.fptsde3d <- function(object,digits=5, ...)
           {  
    class(object) <- "fptsde3d"
    S <- function(t) eval(object$boundary)
    if (object$SDE$x0 > S(as.numeric(object$SDE$t0))){
    fpt_x <- paste("T(S,X) = inf{t >= ",as.numeric(object$SDE$t0) ," : X(t) <= ", deparse(object$boundary),"}")
    }else{
    fpt_x <- paste("T(S,X) = inf{t >= ",as.numeric(object$SDE$t0) ," : X(t) >= ", deparse(object$boundary),"}")
    }
    if (object$SDE$y0 > S(as.numeric(object$SDE$t0))){
    fpt_y <- paste("T(S,Y) = inf{t >= ",as.numeric(object$SDE$t0) ," : Y(t) <= ", deparse(object$boundary),"}")
    }else{
    fpt_y <- paste("T(S,Y) = inf{t >= ",as.numeric(object$SDE$t0) ," : Y(t) >= ", deparse(object$boundary),"}")
    }
    if (object$SDE$z0 > S(as.numeric(object$SDE$t0))){
    fpt_z <- paste("T(S,Z) = inf{t >= ",as.numeric(object$SDE$t0) ," : Z(t) <= ", deparse(object$boundary),"}")
    }else{
    fpt_z <- paste("T(S,Z) = inf{t >= ",as.numeric(object$SDE$t0) ," : Z(t) >= ", deparse(object$boundary),"}")
    }
    cat("\n\tMonte-Carlo Statistics for the F.P.T of (X(t),Y(t),Z(t))","\n",
        "| ",fpt_x,"\n",
        "| ",fpt_y,"\n",
        "| ",fpt_z,"\n",
       sep="")
    x <- object$fptx##[!is.na(X$fptx)]
    y <- object$fpty##[!is.na(X$fpty)]
    z <- object$fptz##[!is.na(X$fpty)]
    res <- as.data.frame(matrix(c(format(length(which(is.na(x))),digits=1),round(digits=digits,mean(x,na.rm = TRUE)),round(digits=digits,var(x,na.rm = TRUE)),round(digits=digits,median(x,na.rm = TRUE)),
                               round(digits=digits,quantile(x,0.25,na.rm = TRUE)),round(digits=digits,quantile(x,0.75,na.rm = TRUE)),
                               round(digits=digits,skewness(x)),round(digits=digits,kurtosis(x)),round(digits=digits,moment(x,order=3)),
                               round(digits=digits,moment(x,order=4)),round(digits=digits,moment(x,order=5)),round(digits=digits,bconfint(x)[1]),round(digits=digits,bconfint(x)[2]),
                               format(length(which(is.na(y))),digits=1),round(digits=digits,mean(y,na.rm = TRUE)),round(digits=digits,var(y,na.rm = TRUE)),round(digits=digits,median(y,na.rm = TRUE)),
                               round(digits=digits,quantile(y,0.25,na.rm = TRUE)),round(digits=digits,quantile(y,0.75,na.rm = TRUE)),
                               round(digits=digits,skewness(y)),round(digits=digits,kurtosis(y)),round(digits=digits,moment(y,order=3)),
                               round(digits=digits,moment(y,order=4)),round(digits=digits,moment(y,order=5)),round(digits=digits,bconfint(y)[1]),round(digits=digits,bconfint(y)[2]),
							   format(length(which(is.na(z))),digits=1),round(digits=digits,mean(z,na.rm = TRUE)),round(digits=digits,var(z,na.rm = TRUE)),round(digits=digits,median(z,na.rm = TRUE)),
                               round(digits=digits,quantile(z,0.25,na.rm = TRUE)),round(digits=digits,quantile(z,0.75,na.rm = TRUE)),
                               round(digits=digits,skewness(z)),round(digits=digits,kurtosis(z)),round(digits=digits,moment(z,order=3)),
                               round(digits=digits,moment(z,order=4)),round(digits=digits,moment(z,order=5)),round(digits=digits,bconfint(z)[1]),round(digits=digits,bconfint(z)[2])),
                               ncol=3))
	dimnames(res) <- list(c("NA's","Mean","Variance","Median","First quartile","Third quartile",
                              "Skewness","Kurtosis","Moment of order 3",
                              "Moment of order 4","Moment of order 5","Int.conf Inf (95%)","Int.conf Sup (95%)"),c("T(S,X)","T(S,Y)","T(S,Z)"))
    print(res, quote = FALSE, right = TRUE,...)
    invisible(object)
}

###

plot.fptsde3d <- function(x,...) .plot.fptsde3d(x,...)

