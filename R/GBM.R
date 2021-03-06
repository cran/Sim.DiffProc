## Fri Apr 03 08:54:53 2020
## Original file Copyright © 2020 A.C. Guidoum, K. Boukhetala
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


GBM <- function(N, ...)  UseMethod("GBM")

GBM.default <- function(N =1000,M=1,x0=1,t0=0,T=1,Dt=NULL,theta=1,sigma=1,...)
             {
    if (any(!is.numeric(x0) || x0 <= 0 )) 
	    stop("'x0' must be numeric > 0")
    if (any(!is.numeric(t0) || !is.numeric(T))) 
	    stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) 
	    stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0)) 
	    stop(" 'M' must be a positive integer ")
    if (any(!is.numeric(sigma) || sigma <= 0) ) 
	    stop(" 'sigma' must be a numeric > 0 ")
    if (any(t0 < 0 || T < 0 || T <= t0) ) 
	        stop(" please use positive times! (0 <= t0 < T) ")
    if (is.null(Dt)) {
        Dt <- (T - t0)/N
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    }
    res <- data.frame(sapply(1:M,function(i) x0*exp((theta - 0.5*sigma^2)*t + sigma*c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt)))))))
    names(res) <- paste("X",1:M,sep="")
    X <- ts(res, start = t0, deltat = Dt)
    return(X)
}
