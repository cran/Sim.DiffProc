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
##### PredCorr1D

.PredCorr1D <- function(N =100,M=1,x0=2,t0=0,T=1,Dt,alpha=0.5,mu=0.5,
                               drift,diffusion,type=c("ito","str"),...)
                       {
    if (missing(Dt)) {
        t <- seq(t0, T, length = N + 1)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
        T <- t[N + 1]
    } 
    if (missing(type)) type <- "ito"
    DSx  <- D(diffusion,"x")
    if (type=="ito"){A    <- function(t,x)  eval(drift)}else{
    A <- function(t,x) eval(drift) - 0.5 * eval(diffusion) * eval(DSx)}
    S    <- function(t,x)  eval(diffusion)
    Sx   <- function(t,x)  eval(DSx)
    SS   <- function(t,x)  A(t,x) - mu * S(t,x) * Sx(t,x)
    predcorr1d <- function()
           {
    w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    dw   <- diff(w)
    X    <- Y <- numeric()
    X[1] <- Y[1] <- x0 
    for (i in 2:(N+1)){
       Y[i] = Y[i-1] + A(t[i-1],Y[i-1])*Dt + S(t[i-1],Y[i-1])*dw[i-1]
       X[i] = X[i-1] +(alpha*SS(t[i],Y[i])+(1-alpha)*SS(t[i-1],X[i-1]))*Dt+(mu*S(t[i],Y[i])+(1-mu)*S(t[i-1],X[i-1]))*dw[i-1]
             }
    X
         }
    res <- data.frame(replicate(M,predcorr1d()))
    names(res) <- paste("X",1:M,sep="")
    X <- ts(res, start = t0, deltat = Dt)
    structure(list(X=X, drift=drift[[1]], diffusion=diffusion[[1]],type=type, x0=x0, N=N, M = M, 
                   Dt=Dt,t0=t0,T=T),class="snssde1d")
}


#####
##### PredCorr2D

.PredCorr2D <- function(N =100,x0=2,y0=1,t0=0,T=1,Dt,alpha=0.5,mu=0.5,driftx,diffx,
                        drifty,diffy,type=c("ito","str"),...)
                       {
    if (missing(Dt)) {
        t <- seq(t0, T, length = N + 1)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
        T <- t[N + 1]
    } 
    if (missing(type)) type <- "ito" 
    DSx  <- D(diffx,"x")
    DSy  <- D(diffy,"y")
    if (type=="ito"){
    Ax    <- function(t,x,y)  eval(driftx)
    Ay    <- function(t,x,y)  eval(drifty)
    }else{
    Ax <- function(t,x,y) eval(driftx) - 0.5 * eval(diffx) * eval(DSx)
    Ay <- function(t,x,y) eval(drifty) - 0.5 * eval(diffy) * eval(DSy)
         }
    Sx    <- function(t,x,y)  eval(diffx)
    dSx   <- function(t,x,y)  eval(DSx)
    SSx   <- function(t,x,y)  Ax(t,x,y) - mu * Sx(t,x,y) * dSx(t,x,y)
    Sy    <- function(t,x,y)  eval(diffy)
    dSy   <- function(t,x,y)  eval(DSy)
    SSy   <- function(t,x,y)  Ay(t,x,y) - mu * Sy(t,x,y) * dSy(t,x,y)
    wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    dwx   <- diff(wx)
    dwy   <- diff(wy)
    X = Y = XX = YY <- numeric()
    X[1] = XX[1] <- x0 
    Y[1] = YY[1] <- y0
    for (i in 2:(N+1)){
       XX[i] = XX[i-1] + Ax(t[i-1],XX[i-1],YY[i-1])*Dt + Sx(t[i-1],XX[i-1],YY[i-1])*dwx[i-1]
       YY[i] = YY[i-1] + Ay(t[i-1],XX[i-1],YY[i-1])*Dt + Sy(t[i-1],XX[i-1],YY[i-1])*dwy[i-1]
       X[i] = X[i-1] +(alpha*SSx(t[i],XX[i],YY[i])+(1-alpha)*SSx(t[i-1],X[i-1],Y[i-1]))*Dt+(mu*Sx(t[i],XX[i],YY[i])+(1-mu)*Sx(t[i-1],X[i-1],Y[i-1]))*dwx[i-1]
       Y[i] = Y[i-1] +(alpha*SSy(t[i],XX[i],YY[i])+(1-alpha)*SSy(t[i-1],X[i-1],Y[i-1]))*Dt+(mu*Sy(t[i],XX[i],YY[i])+(1-mu)*Sy(t[i-1],X[i-1],Y[i-1]))*dwy[i-1]
      }
    res <- data.frame(X,Y)
    names(res) <- paste(c("X","Y"),sep="")
    X <- ts(res, start = t0, deltat = Dt)
    structure(list(X=X, driftx=driftx[[1]], diffx=diffx[[1]],drifty=drifty[[1]], diffy=diffy[[1]],
                  type=type, x0=x0,y0=y0, N=N, Dt=Dt,t0=t0,T=T),class="snssde2d")
} 


#####
##### PredCorr3D

.PredCorr3D <- function(N =100,x0=2,y0=1,z0=1,t0=0,T=1,Dt,alpha=0.5,mu=0.5,driftx,diffx,drifty,diffy,
                     driftz,diffz,type=c("ito","str"),...)
                       {
    if (missing(Dt)) {
        t <- seq(t0, T, length = N + 1)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
        T <- t[N + 1]
    } 
    if (missing(type)) type <- "ito" 
    DSx  <- D(diffx,"x")
    DSy  <- D(diffy,"y")
    DSz  <- D(diffz,"z")
    if (type=="ito"){
    Ax    <- function(t,x,y,z)  eval(driftx)
    Ay    <- function(t,x,y,z)  eval(drifty)
    Az    <- function(t,x,y,z)  eval(driftz)
    }else{
    Ax <- function(t,x,y,z) eval(driftx) - 0.5 * eval(diffx) * eval(DSx)
    Ay <- function(t,x,y,z) eval(drifty) - 0.5 * eval(diffy) * eval(DSy)
    Az <- function(t,x,y,z) eval(driftz) - 0.5 * eval(diffz) * eval(DSz)
         }
    Sx    <- function(t,x,y,z)  eval(diffx)
    dSx   <- function(t,x,y,z)  eval(DSx)
    SSx   <- function(t,x,y,z)  Ax(t,x,y,z) - mu * Sx(t,x,y,z) * dSx(t,x,y,z)
    Sy    <- function(t,x,y,z)  eval(diffy)
    dSy   <- function(t,x,y,z)  eval(DSy)
    SSy   <- function(t,x,y,z)  Ay(t,x,y,z) - mu * Sy(t,x,y,z) * dSy(t,x,y,z)
    Sz    <- function(t,x,y,z)  eval(diffz)
    dSz   <- function(t,x,y,z)  eval(DSz)
    SSz   <- function(t,x,y,z)  Az(t,x,y,z) - mu * Sz(t,x,y,z) * dSz(t,x,y,z)
    wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    wz = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    dwx   <- diff(wx)
    dwy   <- diff(wy)
    dwz   <- diff(wz)
    X = Y = Z = XX = YY = ZZ <- numeric()
    X[1] = XX[1] <- x0 
    Y[1] = YY[1] <- y0
    Z[1] = ZZ[1] <- z0
    for (i in 2:(N+1)){
       XX[i] = XX[i-1] + Ax(t[i-1],XX[i-1],YY[i-1],ZZ[i-1])*Dt + Sx(t[i-1],XX[i-1],YY[i-1],ZZ[i-1])*dwx[i-1]
       YY[i] = YY[i-1] + Ay(t[i-1],XX[i-1],YY[i-1],ZZ[i-1])*Dt + Sy(t[i-1],XX[i-1],YY[i-1],ZZ[i-1])*dwy[i-1]
       ZZ[i] = ZZ[i-1] + Az(t[i-1],XX[i-1],YY[i-1],ZZ[i-1])*Dt + Sz(t[i-1],XX[i-1],YY[i-1],ZZ[i-1])*dwz[i-1]
       X[i] = X[i-1] +(alpha*SSx(t[i],XX[i],YY[i],ZZ[i])+(1-alpha)*SSx(t[i-1],X[i-1],Y[i-1],Z[i-1]))*Dt+(mu*Sx(t[i],XX[i],YY[i],ZZ[i])+(1-mu)*Sx(t[i-1],X[i-1],Y[i-1],Z[i-1]))*dwx[i-1]
       Y[i] = Y[i-1] +(alpha*SSy(t[i],XX[i],YY[i],ZZ[i])+(1-alpha)*SSy(t[i-1],X[i-1],Y[i-1],Z[i-1]))*Dt+(mu*Sy(t[i],XX[i],YY[i],ZZ[i])+(1-mu)*Sy(t[i-1],X[i-1],Y[i-1],Z[i-1]))*dwy[i-1]
       Z[i] = Z[i-1] +(alpha*SSz(t[i],XX[i],YY[i],ZZ[i])+(1-alpha)*SSz(t[i-1],X[i-1],Y[i-1],Z[i-1]))*Dt+(mu*Sz(t[i],XX[i],YY[i],ZZ[i])+(1-mu)*Sz(t[i-1],X[i-1],Y[i-1],Z[i-1]))*dwz[i-1]
      }
    res <- data.frame(X,Y,Z)
    names(res) <- paste(c("X","Y","Z"),sep="")
    X <- ts(res, start = t0, deltat = Dt)
    structure(list(X=X, driftx=driftx[[1]], diffx=diffx[[1]],drifty=drifty[[1]], diffy=diffy[[1]],
                  driftz=driftz[[1]], diffz=diffz[[1]],type=type, x0=x0,y0=y0,z0=z0, N=N, Dt=Dt,t0=t0,T=T),class="snssde3d")
}


