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

rfptsde1d <- function(object, ...)  UseMethod("rfptsde1d")

rfptsde1d.default <- function(object,boundary,...)
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
	#if (length(fpt != object$M)) {cat("missing output are removed\n")}
	fpt <- fpt[which(!is.na(fpt))]
    #structure(list(SDE=object,boundary=boundary[[1]],fpt=fpt),class="rfptsde1d")
	return(fpt)
}

###

# print.rfptsde1d <- function(x, digits=NULL, ...)
           # {
    # class(x) <- "rfptsde1d"
    # print(x$fpt,...)
# }

# summary.rfptsde1d <- function(object, digits=5,...)
           # {
    # x <- object	   
    # class(x) <- "rfptsde1d"
    # S <- function(t) eval(x$boundary)
    # if (as.numeric(x$SDE$x0) > S(as.numeric(x$SDE$t0))){
    # cat("\n Monte-Carlo Statistics for the F.P.T of X(t)","\n", " T(S,X) = inf{t >= ",x$SDE$t0 ," : X(t) <= ", deparse(x$boundary),"}","\n",
        # sep="")
    # }else{
    # cat("\n Monte-Carlo Statistics for the F.P.T of X(t)","\n", " T(S,X) = inf{t >= ",x$SDE$t0 ," : X(t) >= ", deparse(x$boundary),"}","\n",
        # sep="")
    # }
    # res <- as.data.frame(matrix(c(format(length(which(is.na(x$fpt))),digits=1),round(digits=digits,mean(x$fpt,na.rm = TRUE)),round(digits=digits,var(x$fpt,na.rm = TRUE)),
                               # round(digits=digits,median(x$fpt,na.rm = TRUE)),round(digits=digits,quantile(x$fpt,0.25,na.rm = TRUE)),
							   # round(digits=digits,quantile(x$fpt,0.75,na.rm = TRUE)),round(digits=digits,skewness(x$fpt)),round(digits=digits,kurtosis(x$fpt)),
							   # round(digits=digits,moment(x$fpt,order=3)),round(digits=digits,moment(x$fpt,order=4)),
							   # round(digits=digits,moment(x$fpt,order=5)),round(digits=digits,bconfint(x$fpt)[1]),round(digits=digits,bconfint(x$fpt)[2])),
                               # ncol=1))
	# dimnames(res) <- list(c("NA's","Mean","Variance","Median","First quartile","Third quartile",
                              # "Skewness","Kurtosis","Moment of order 3",
                              # "Moment of order 4","Moment of order 5","Int.conf Inf (95%)","Int.conf Sup (95%)"),c("T(S,X)"))
    # print(res, quote = FALSE, right = TRUE,...)
    # invisible(x)
# }

#plot.rfptsde1d <- function(x,...) .plot.rfptsde1d(x,...)



## density fpt

dfptsde1d <- function(object, ...)  UseMethod("dfptsde1d")

dfptsde1d.default <- function(object,boundary,...)
                     {
    class(object) <- "snssde1d"
    if (any(!is.expression(boundary))) stop(" must be expressions of a constant or time-dependent boundary ")
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
      structure(list(SDE=object,ech=fpt,res=density.default(fpt,na.rm = TRUE,...),boundary=boundary[[1]]),class="dfptsde1d")
}

print.dfptsde1d <- function(x, digits=NULL, ...)
           {
    class(x) <- "dfptsde1d"
    S <- function(t) eval(x$boundary)
   if (as.numeric(x$SDE$x0) > S(as.numeric(x$SDE$t0))){
    cat("\n Kernel density for the F.P.T of X(t)","\n", " T(S,X) = inf{t >= ",x$SDE$t0 ," : X(t) <= ", deparse(x$boundary),"}","\n",
        sep="")
    }else{
    cat("\n Kernel density for the F.P.T of X(t)","\n", " T(S,X) = inf{t >= ",x$SDE$t0 ," : X(t) >= ", deparse(x$boundary),"}","\n",
        sep="")}
    cat(
	"\nData: ", x$res$data.name, " (", x$res$n, " obs.);",
	"\tBandwidth 'bw' = ", formatC(x$res$bw, digits = digits), "\n\n", sep = "")
    out <- as.data.frame(x$res[c("x","y")])
    names(out) <- c("x","f(x)")
    print(summary(out, digits = digits), ...)
    invisible(x)
}

plot.dfptsde1d <- function(x,...) .plot.dfptsde1d(x,...)


################################################################################
################################################################################
#####
##### fptsde2d

rfptsde2d <- function(object, ...)  UseMethod("rfptsde2d")


rfptsde2d.default <- function(object,boundary,...)
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
    #structure(list(SDE=object,boundary=boundary[[1]],fpt=na.omit(data.frame(fptx,fpty))),class="rfptsde2d")	
	if (length(which(is.na(fptx))) != length(which(is.na(fpty)))) cat("missing output are removed\n")
	fpt <- na.omit(data.frame(fptx,fpty))
	out <- data.frame(fpt$fptx,fpt$fpty)
	names(out) <- c("x","y")
	return(out)
}


# print.rfptsde2d <- function(x, digits=NULL, ...)
           # {
    # class(x) <- "rfptsde2d"
    # out <- data.frame(x$fptx,x$fpty)
    # names(out) <- c("fptx","fpty")
    # print(out, digits = digits, ...)
    # invisible(x)
# }

###

# summary.rfptsde2d <- function(object,digits=5, ...)
           # {  
    # class(object) <- "rfptsde2d"
    # S <- function(t) eval(object$boundary)
    # if (as.numeric(object$SDE$x0) > S(as.numeric(object$SDE$t0))){
    # fpt_x <- paste("T(S,X) = inf{t >= ",as.numeric(object$SDE$t0) ," : X(t) <= ", deparse(object$boundary),"}")
    # }else{
    # fpt_x <- paste("T(S,X) = inf{t >= ",as.numeric(object$SDE$t0) ," : X(t) >= ", deparse(object$boundary),"}")
    # }
    # if (as.numeric(object$SDE$y0) > S(as.numeric(object$SDE$t0))){
    # fpt_y <- paste("T(S,Y) = inf{t >= ",as.numeric(object$SDE$t0) ," : Y(t) <= ", deparse(object$boundary),"}")
    # }else{
    # fpt_y <- paste("T(S,Y) = inf{t >= ",as.numeric(object$SDE$t0) ," : Y(t) >= ", deparse(object$boundary),"}")
    # }
    # cat("\n Monte-Carlo Statistics for the F.P.T of (X(t),Y(t))","\n",
        # " ",fpt_x,"\n",
        # " ",fpt_y,"\n",
       # sep="")
    # x <- object$fpt[,1]##[!is.na(X$fptx)]
    # y <- object$fpt[,2]##[!is.na(X$fpty)]
    # res <- as.data.frame(matrix(c(format(length(which(is.na(x))),digits=1),round(digits=digits,mean(x,na.rm = TRUE)),round(digits=digits,var(x,na.rm = TRUE)),round(digits=digits,median(x,na.rm = TRUE)),
                               # round(digits=digits,quantile(x,0.25,na.rm = TRUE)),round(digits=digits,quantile(x,0.75,na.rm = TRUE)),
                               # round(digits=digits,skewness(x)),round(digits=digits,kurtosis(x)),round(digits=digits,moment(x,order=3)),
                               # round(digits=digits,moment(x,order=4)),round(digits=digits,moment(x,order=5)),round(digits=digits,bconfint(x)[1]),round(digits=digits,bconfint(x)[2]),
                               # format(length(which(is.na(y))),digits=1),round(digits=digits,mean(y,na.rm = TRUE)),round(digits=digits,var(y,na.rm = TRUE)),round(digits=digits,median(y,na.rm = TRUE)),
                               # round(digits=digits,quantile(y,0.25,na.rm = TRUE)),round(digits=digits,quantile(y,0.75,na.rm = TRUE)),
                               # round(digits=digits,skewness(y)),round(digits=digits,kurtosis(y)),round(digits=digits,moment(y,order=3)),
                               # round(digits=digits,moment(y,order=4)),round(digits=digits,moment(y,order=5)),round(digits=digits,bconfint(y)[1]),round(digits=digits,bconfint(y)[2])),
                               # ncol=2))
	# dimnames(res) <- list(c("NA's","Mean","Variance","Median","First quartile","Third quartile",
                              # "Skewness","Kurtosis","Moment of order 3",
                              # "Moment of order 4","Moment of order 5","Int.conf Inf (95%)","Int.conf Sup (95%)"),c("T(S,X)","T(S,Y)"))
    # print(res, quote = FALSE, right = TRUE,...)
    # invisible(object)
# }


#plot.rfptsde2d <- function(x,...) .plot.rfptsde2d(x,...)
	
## density fpt2d

dfptsde2d <- function(object, ...)  UseMethod("dfptsde2d")

dfptsde2d.default <- function(object,boundary,pdf=c("Joint","Marginal"),...)
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
    if (length(which(is.na(fptx))) != length(which(is.na(fpty)))) cat("missing output are removed\n")
	fpt <- na.omit(data.frame(fptx,fpty))
	out <- data.frame(fpt$fptx,fpt$fpty)
	names(out) <- c("x","y")
    pdf <- match.arg(pdf)
    if (pdf=="Marginal"){
    structure(list(fpt=out,resx=density.default(out[,"x"],na.rm = TRUE,...),resy=density.default(out[,"y"],na.rm = TRUE,...),boundary=boundary[[1]],pdf=pdf,SDE=object),class="dfptsde2d")
    }else{
    structure(list(fpt=out,res=MASS::kde2d(out$x, out$y,n = 100, ...),boundary=boundary[[1]],pdf=pdf,SDE=object),class="dfptsde2d")
    }
}

print.dfptsde2d <- function(x, digits=NULL, ...)
           {
    class(x) <- "dfptsde2d"
    S <- function(t) eval(x$boundary)
    if (x$pdf=="Joint") {
    if (as.numeric(x$SDE$x0) > S(as.numeric(x$SDE$t0))){
     fpt_x <- paste(" : X(t) <= ", deparse(x$boundary)," ")
     }else{
     fpt_x <- paste(" : X(t) >= ", deparse(x$boundary)," ")
     }
     if (as.numeric(x$SDE$y0) > S(as.numeric(x$SDE$t0))){
     fpt_y <- paste(" Y(t) <= ", deparse(x$boundary))
     }else{
     fpt_y <- paste(" Y(t) >= ", deparse(x$boundary))
     }
     cat("\n Joint density for the F.P.T of (X(t),Y(t))","\n",
         " T(S,X,Y) = inf{t >= ",x$SDE$t0 ,fpt_x,"and",fpt_y,"}","\n",
        sep="")
    cat(
	"\nData: (x,y)", " (2 x ", dim(x$fpt)[1], " obs.);", "\n\n", sep = "")
	##"\tBandwidth 'bw' = ", formatC(MASS::bandwidth.nrd(x$res$x), digits = digits),"~~", formatC(MASS::bandwidth.nrd(x$res$y), digits = digits), "\n\n", sep = "")
    out3 <- list(x=x$res$x,y=x$res$y,z=as.vector(x$res$z))
    names(out3) <- c("x","y","f(x,y)")
    print(summary.data.frame(out3, digits = digits), ...)}else{
    if (as.numeric(x$SDE$x0) > S(as.numeric(x$SDE$t0))){
    cat("\n Marginal density for the F.P.T of X(t)","\n", " T(S,X) = inf{t >= ",x$SDE$t0 ," : X(t) <= ", deparse(x$boundary),"}","\n",
        sep="")
    }else{
    cat("\n Marginal density for the F.P.T of X(t)","\n", " T(S,X) = inf{t >= ",x$SDE$t0 ," : X(t) >= ", deparse(x$boundary),"}","\n",
        sep="")}
    cat(
	"\nData: ", x$resx$data.name, " (", x$resx$n, " obs.);",
	"\tBandwidth 'bw' = ", formatC(x$resx$bw, digits = digits), "\n\n", sep = "")
    out1 <- as.data.frame(x$resx[c("x","y")])
    names(out1) <- c("x","f(x)")
    print(summary(out1, digits = digits), ...)

    if (as.numeric(x$SDE$y0) > S(as.numeric(x$SDE$t0))){
    cat("\n Marginal density for the F.P.T of Y(t)","\n", " T(S,Y) = inf{t >= ",x$SDE$t0 ," : Y(t) <= ", deparse(x$boundary),"}","\n",
        sep="")
    }else{
    cat("\n Marginal density for the F.P.T of Y(t)","\n", " T(S,Y) = inf{t >= ",x$SDE$t0 ," : Y(t) >= ", deparse(x$boundary),"}","\n",
        sep="")}
    cat(
	"\nData: ", x$resy$data.name, " (", x$resy$n, " obs.);",
	"\tBandwidth 'bw' = ", formatC(x$resy$bw, digits = digits), "\n\n", sep = "")
    out2 <- as.data.frame(x$resy[c("x","y")])
    names(out2) <- c("y","f(y)")
    print(summary(out2, digits = digits), ...)
   }
    invisible(x)
}


plot.dfptsde2d <- function(x,display=c("persp","rgl","image","contour"),hist=FALSE,...) .plot.dfptsde2d(x,display,hist,...)

################################################################################
################################################################################
#####
##### fptsde3d

rfptsde3d <- function(object, ...)  UseMethod("rfptsde3d")

rfptsde3d.default <- function(object,boundary,...)
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
    #structure(list(SDE=object,boundary=boundary[[1]],fptx=fptx,fpty=fpty,fptz=fptz),class="fptsde3d")
	if (length(which(is.na(fptx))) != length(which(is.na(fpty))) | length(which(is.na(fptx))) != length(which(is.na(fptz))) | length(which(is.na(fpty))) != length(which(is.na(fptz)))) cat("missing output are removed\n")
	fpt <- na.omit(data.frame(fptx,fpty,fptz))
	out <- data.frame(fpt$fptx,fpt$fpty,fpt$fptz)
	names(out) <- c("x","y","z")
	return(out)
}

###

# summary.fptsde3d <- function(object,digits=5, ...)
           # {  
    # class(object) <- "fptsde3d"
    # S <- function(t) eval(object$boundary)
    # if (object$SDE$x0 > S(as.numeric(object$SDE$t0))){
    # fpt_x <- paste("T(S,X) = inf{t >= ",as.numeric(object$SDE$t0) ," : X(t) <= ", deparse(object$boundary),"}")
    # }else{
    # fpt_x <- paste("T(S,X) = inf{t >= ",as.numeric(object$SDE$t0) ," : X(t) >= ", deparse(object$boundary),"}")
    # }
    # if (object$SDE$y0 > S(as.numeric(object$SDE$t0))){
    # fpt_y <- paste("T(S,Y) = inf{t >= ",as.numeric(object$SDE$t0) ," : Y(t) <= ", deparse(object$boundary),"}")
    # }else{
    # fpt_y <- paste("T(S,Y) = inf{t >= ",as.numeric(object$SDE$t0) ," : Y(t) >= ", deparse(object$boundary),"}")
    # }
    # if (object$SDE$z0 > S(as.numeric(object$SDE$t0))){
    # fpt_z <- paste("T(S,Z) = inf{t >= ",as.numeric(object$SDE$t0) ," : Z(t) <= ", deparse(object$boundary),"}")
    # }else{
    # fpt_z <- paste("T(S,Z) = inf{t >= ",as.numeric(object$SDE$t0) ," : Z(t) >= ", deparse(object$boundary),"}")
    # }
    # cat("\n\tMonte-Carlo Statistics for the F.P.T of (X(t),Y(t),Z(t))","\n",
        # "| ",fpt_x,"\n",
        # "| ",fpt_y,"\n",
        # "| ",fpt_z,"\n",
       # sep="")
    # x <- object$fptx##[!is.na(X$fptx)]
    # y <- object$fpty##[!is.na(X$fpty)]
    # z <- object$fptz##[!is.na(X$fpty)]
    # res <- as.data.frame(matrix(c(format(length(which(is.na(x))),digits=1),round(digits=digits,mean(x,na.rm = TRUE)),round(digits=digits,var(x,na.rm = TRUE)),round(digits=digits,median(x,na.rm = TRUE)),
                               # round(digits=digits,quantile(x,0.25,na.rm = TRUE)),round(digits=digits,quantile(x,0.75,na.rm = TRUE)),
                               # round(digits=digits,skewness(x)),round(digits=digits,kurtosis(x)),round(digits=digits,moment(x,order=3)),
                               # round(digits=digits,moment(x,order=4)),round(digits=digits,moment(x,order=5)),round(digits=digits,bconfint(x)[1]),round(digits=digits,bconfint(x)[2]),
                               # format(length(which(is.na(y))),digits=1),round(digits=digits,mean(y,na.rm = TRUE)),round(digits=digits,var(y,na.rm = TRUE)),round(digits=digits,median(y,na.rm = TRUE)),
                               # round(digits=digits,quantile(y,0.25,na.rm = TRUE)),round(digits=digits,quantile(y,0.75,na.rm = TRUE)),
                               # round(digits=digits,skewness(y)),round(digits=digits,kurtosis(y)),round(digits=digits,moment(y,order=3)),
                               # round(digits=digits,moment(y,order=4)),round(digits=digits,moment(y,order=5)),round(digits=digits,bconfint(y)[1]),round(digits=digits,bconfint(y)[2]),
							   # format(length(which(is.na(z))),digits=1),round(digits=digits,mean(z,na.rm = TRUE)),round(digits=digits,var(z,na.rm = TRUE)),round(digits=digits,median(z,na.rm = TRUE)),
                               # round(digits=digits,quantile(z,0.25,na.rm = TRUE)),round(digits=digits,quantile(z,0.75,na.rm = TRUE)),
                               # round(digits=digits,skewness(z)),round(digits=digits,kurtosis(z)),round(digits=digits,moment(z,order=3)),
                               # round(digits=digits,moment(z,order=4)),round(digits=digits,moment(z,order=5)),round(digits=digits,bconfint(z)[1]),round(digits=digits,bconfint(z)[2])),
                               # ncol=3))
	# dimnames(res) <- list(c("NA's","Mean","Variance","Median","First quartile","Third quartile",
                              # "Skewness","Kurtosis","Moment of order 3",
                              # "Moment of order 4","Moment of order 5","Int.conf Inf (95%)","Int.conf Sup (95%)"),c("T(S,X)","T(S,Y)","T(S,Z)"))
    # print(res, quote = FALSE, right = TRUE,...)
    # invisible(object)
# }

###

# plot.fptsde3d <- function(x,...) .plot.fptsde3d(x,...)

## density fpt3d

dfptsde3d <- function(object, ...)  UseMethod("dfptsde3d")

dfptsde3d.default <- function(object,boundary,...)
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
    #structure(list(SDE=object,boundary=boundary[[1]],fptx=fptx,fpty=fpty,fptz=fptz),class="fptsde3d")
	if (length(which(is.na(fptx))) != length(which(is.na(fpty))) | length(which(is.na(fptx))) != length(which(is.na(fptz))) | length(which(is.na(fpty))) != length(which(is.na(fptz)))) cat("missing output are removed\n")
	fpt <- na.omit(data.frame(fptx,fpty,fptz))
	out <- data.frame(fpt$fptx,fpt$fpty,fpt$fptz)
	names(out) <- c("x","y","z") 
    structure(list(fpt=out,resx=density.default(out[,"x"],na.rm = TRUE,...),resy=density.default(out[,"y"],na.rm = TRUE,...),resz=density.default(out[,"z"],na.rm = TRUE,...),boundary=boundary[[1]],SDE=object),class="dfptsde3d")
}

print.dfptsde3d <- function(x, digits=NULL, ...)
           {
    class(x) <- "dfptsde3d"
    S <- function(t) eval(x$boundary)
    if (as.numeric(x$SDE$x0) > S(as.numeric(x$SDE$t0))){
    cat("\n Marginal density for the F.P.T of X(t)","\n", " T(S,X) = inf{t >= ",x$SDE$t0 ," : X(t) <= ", deparse(x$boundary),"}","\n",
        sep="")
    }else{
    cat("\n Marginal density for the F.P.T of X(t)","\n", " T(S,X) = inf{t >= ",x$SDE$t0 ," : X(t) >= ", deparse(x$boundary),"}","\n",
        sep="")}
    cat(
	"\nData: ", x$resx$data.name, " (", x$resx$n, " obs.);",
	"\tBandwidth 'bw' = ", formatC(x$resx$bw, digits = digits), "\n\n", sep = "")
    out1 <- as.data.frame(x$resx[c("x","y")])
    names(out1) <- c("x","f(x)")
    print(summary(out1, digits = digits), ...)

    if (as.numeric(x$SDE$y0) > S(as.numeric(x$SDE$t0))){
    cat("\n Marginal density for the F.P.T of Y(t)","\n", " T(S,Y) = inf{t >= ",x$SDE$t0 ," : Y(t) <= ", deparse(x$boundary),"}","\n",
        sep="")
    }else{
    cat("\n Marginal density for the F.P.T of Y(t)","\n", " T(S,Y) = inf{t >= ",x$SDE$t0 ," : Y(t) >= ", deparse(x$boundary),"}","\n",
        sep="")}
    cat(
	"\nData: ", x$resy$data.name, " (", x$resy$n, " obs.);",
	"\tBandwidth 'bw' = ", formatC(x$resy$bw, digits = digits), "\n\n", sep = "")
    out2 <- as.data.frame(x$resy[c("x","y")])
    names(out2) <- c("y","f(y)")
    print(summary(out2, digits = digits), ...)

    if (as.numeric(x$SDE$z0) > S(as.numeric(x$SDE$t0))){
    cat("\n Marginal density for the F.P.T of Z(t)","\n", " T(S,Z) = inf{t >= ",x$SDE$t0 ," : Z(t) <= ", deparse(x$boundary),"}","\n",
        sep="")
    }else{
    cat("\n Marginal density for the F.P.T of Z(t)","\n", " T(S,Z) = inf{t >= ",x$SDE$t0 ," : Z(t) >= ", deparse(x$boundary),"}","\n",
        sep="")}
    cat(
	"\nData: ", x$resz$data.name, " (", x$resz$n, " obs.);",
	"\tBandwidth 'bw' = ", formatC(x$resz$bw, digits = digits), "\n\n", sep = "")
    out3 <- as.data.frame(x$resz[c("x","y")])
    names(out3) <- c("z","f(z)")
    print(summary(out3, digits = digits), ...)
    invisible(x)
}

plot.dfptsde3d <- function(x,hist=FALSE,...) .plot.dfptsde3d(x,hist,...)
