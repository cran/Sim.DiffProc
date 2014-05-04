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
##### Milstein1D

.Milstein1D <- function(N =100,M=1,x0=2,t0=0,T=1,Dt,drift,diffusion,
                          type=c("ito","str"),...)
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
    A  <- function(t,x)  eval(drift) - 0.5 * eval(diffusion) * eval(DSx)}
    S  <- function(t,x)  eval(diffusion)
    Sx <- function(t,x)  eval(DSx)
    milstein1d <- function()
           {
    w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    dw   <- diff(w)
    X    <- numeric()
    X[1] <- x0 
    for (i in 2:(N+1)){
        X[i] = X[i-1]+ A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*dw[i-1]+0.5 *S(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])*((dw[i-1])^2 -Dt)
       }
    X
         }
    res <- data.frame(replicate(M,milstein1d()))
    names(res) <- paste("X",1:M,sep="")
    X <- ts(res, start = t0, deltat = Dt)
    structure(list(X=X, drift=drift[[1]], diffusion=diffusion[[1]],type=type, x0=x0, N=N, M = M, 
                   Dt=Dt,t0=t0,T=T),class="snssde1d")
}           

#####
##### Milstein2D

.Milstein2D <- function(N =100,x0=2,y0=1,t0=0,T=1,Dt,driftx,diffx,drifty,diffy,
                          type=c("ito","str"),...)
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
    if (type=="ito"){Ax    <- function(t,x,y) eval(driftx)
                     Ay    <- function(t,x,y) eval(drifty)}else{
    Ax  <- function(t,x,y) eval(driftx) - 0.5 * eval(diffx) * eval(DSx)
    Ay  <- function(t,x,y) eval(drifty) - 0.5 * eval(diffy) * eval(DSy)}
    Sx  <- function(t,x,y) eval(diffx)
    Sy  <- function(t,x,y) eval(diffy) 
    dSx <- function(t,x,y) eval(DSx)
    dSy <- function(t,x,y) eval(DSy)
    wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    dwx   <- diff(wx)
    dwy   <- diff(wy)
    X = Y <- numeric()
    X[1] <- x0 
    Y[1] <- y0
    for (i in 2:(N+1)){
        X[i] = X[i-1]+ Ax(t[i-1],X[i-1],Y[i-1])*Dt + Sx(t[i-1],X[i-1],Y[i-1])*dwx[i-1]+0.5 *Sx(t[i-1],X[i-1],Y[i-1])*dSx(t[i-1],X[i-1],Y[i-1])*((dwx[i-1])^2 -Dt)
        Y[i] = Y[i-1]+ Ay(t[i-1],X[i-1],Y[i-1])*Dt + Sy(t[i-1],X[i-1],Y[i-1])*dwy[i-1]+0.5 *Sy(t[i-1],X[i-1],Y[i-1])*dSy(t[i-1],X[i-1],Y[i-1])*((dwy[i-1])^2 -Dt)
       }
    res <- data.frame(X,Y)
    names(res) <- paste(c("X","Y"),sep="")
    X <- ts(res, start = t0, deltat = Dt)
    structure(list(X=X, driftx=driftx[[1]], diffx=diffx[[1]],drifty=drifty[[1]], diffy=diffy[[1]],
                  type=type, x0=x0,y0=y0, N=N, Dt=Dt,t0=t0,T=T),class="snssde2d")
} 

#####
##### Milstein3D

.Milstein3D <- function(N =100,x0=2,y0=1,z0=1,t0=0,T=1,Dt,driftx,diffx,drifty,diffy,
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
    if (type=="ito"){Ax    <- function(t,x,y,z) eval(driftx)
        Ay    <- function(t,x,y,z) eval(drifty)
        Az    <- function(t,x,y,z) eval(driftz)}else{
    Ax <- function(t,x,y,z) eval(driftx) - 0.5 * eval(diffx) * eval(D(diffx,"x"))
    Ay <- function(t,x,y,z) eval(drifty) - 0.5 * eval(diffy) * eval(D(diffy,"y"))
    Az <- function(t,x,y,z) eval(driftz) - 0.5 * eval(diffz) * eval(D(diffz,"z"))}
    Sx  <- function(t,x,y,z) eval(diffx)
    Sy  <- function(t,x,y,z) eval(diffy)
    Sz  <- function(t,x,y,z) eval(diffz) 
    dSx <- function(t,x,y,z) eval(DSx)
    dSy <- function(t,x,y,z) eval(DSy)
    dSz <- function(t,x,y,z) eval(DSz)
    wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    wz = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    dwx   <- diff(wx)
    dwy   <- diff(wy)
    dwz   <- diff(wz)
    X = Y = Z <- numeric()
    X[1] <- x0 
    Y[1] <- y0
    Z[1] <- z0
    for (i in 2:(N+1)){
        X[i] = X[i-1]+ Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt + Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])*dwx[i-1]+0.5 *Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])*dSx(t[i-1],X[i-1],Y[i-1],Z[i-1])*((dwx[i-1])^2 -Dt)
        Y[i] = Y[i-1]+ Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt + Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])*dwy[i-1]+0.5 *Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])*dSy(t[i-1],X[i-1],Y[i-1],Z[i-1])*((dwy[i-1])^2 -Dt)
        Z[i] = Z[i-1]+ Az(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt + Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])*dwz[i-1]+0.5 *Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])*dSz(t[i-1],X[i-1],Y[i-1],Z[i-1])*((dwz[i-1])^2 -Dt)
       }
    res <- data.frame(X,Y,Z)
    names(res) <- paste(c("X","Y","Z"),sep="")
    X <- ts(res, start = t0, deltat = Dt)
    structure(list(X=X, driftx=driftx[[1]], diffx=diffx[[1]],drifty=drifty[[1]], diffy=diffy[[1]],
                  driftz=driftz[[1]], diffz=diffz[[1]],type=type, x0=x0,y0=y0,z0=z0, N=N, Dt=Dt,t0=t0,T=T),class="snssde3d")
} 


