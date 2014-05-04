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
##### STS1D

.STS1D <- function(N =100,M=1,x0=2,t0=0,T=1,Dt,drift,diffusion,
                          type=c("ito","str"),...)
                       {
    if (missing(Dt)) {
        t <- seq(t0, T, length = N + 1)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
        T <- t[N + 1]
    } 					   
    if (missing(type)) type <- "ito"  					   
    if (type=="ito"){
    A    <- function(t,x)  eval(drift)
    Ax   <- function(t,x)  eval(D(drift,"x"))
    Axx  <- function(t,x)  eval(D(D(drift,"x"),"x"))
    }else{
    A    <- function(t,x)  eval(drift) - 0.5 * eval(diffusion) * eval(D(diffusion,"x"))
    Ax   <- function(t,x)  eval(D(drift,"x")) - 0.5 * (eval(D(diffusion,"x")) * eval(D(diffusion,"x"))+ eval(diffusion) * eval(D(D(diffusion,"x"),"x")))
    Axx  <- function(t,x)  eval(D(D(drift,"x"),"x")) - 0.5 * ( eval(D(D(diffusion,"x"),"x")) * eval(D(diffusion,"x"))+ eval(D(diffusion,"x")) * eval(D(D(diffusion,"x"),"x"))+
                           eval(D(diffusion,"x")) * eval(D(D(diffusion,"x"),"x")) + eval(diffusion) * eval(D(D(D(diffusion,"x"),"x"),"x")) )
                  }
    DSx  <- D(diffusion,"x")  
    DSxx <- D(DSx,"x")
    S    <- function(t,x)  eval(diffusion)
    Sx   <- function(t,x)  eval(DSx)
    Sxx  <- function(t,x)  eval(DSxx)
    sts1d <- function()
           {
    w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    dw   <- diff(w)
    dz   <- rnorm(N,0,sqrt((1/3)*Dt^3))
    X    <- numeric()
    X[1] <- x0 
    for (i in 2:(N+1)){
       X[i]=X[i-1]+A(t[i-1],X[i-1])*Dt+S(t[i-1],X[i-1])*dw[i-1]+ 0.5*S(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])*((dw[i-1]^2)-Dt)+
            Ax(t[i-1],X[i-1])*S(t[i-1],X[i-1])*dz[i-1]+0.5*(A(t[i-1],X[i-1])*Ax(t[i-1],X[i-1])+0.5*(S(t[i-1],X[i-1])^2)*
            Axx(t[i-1],X[i-1]))*(Dt^2)+(A(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])+0.5*(S(t[i-1],X[i-1])^2)*Sxx(t[i-1],X[i-1]))*
            (dw[i-1]*Dt-dz[i-1])+0.5*S(t[i-1],X[i-1])*(S(t[i-1],X[i-1])*Sxx(t[i-1],X[i-1])+(Sx(t[i-1],X[i-1])^2))*((1/3)*(dw[i-1]^2)-Dt)*dw[i-1]
                      }
    X
         }
    res <- data.frame(replicate(M,sts1d()))
    names(res) <- paste("X",1:M,sep="")
    X <- ts(res, start = t0, deltat = Dt)
    structure(list(X=X, drift=drift[[1]], diffusion=diffusion[[1]],type=type, x0=x0, N=N, M = M, 
                   Dt=Dt,t0=t0,T=T),class="snssde1d")
}           

#####
##### STS2D

.STS2D <- function(N =100,x0=2,y0=1,t0=0,T=1,Dt,driftx,diffx,drifty,diffy,
                          type=c("ito","str"),...)
                       {
    if (missing(Dt)) {
        t <- seq(t0, T, length = N + 1)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
        T <- t[N + 1]
    } 					   
    if (missing(type)) type <- "ito"  					   
    if (type=="ito"){
    Ax    <- function(t,x,y)  eval(driftx)
    dAx   <- function(t,x,y)  eval(D(driftx,"x"))
    dAxx  <- function(t,x,y)  eval(D(D(driftx,"x"),"x"))
    Ay    <- function(t,x,y)  eval(drifty)
    dAy   <- function(t,x,y)  eval(D(drifty,"y"))
    dAyy  <- function(t,x,y)  eval(D(D(drifty,"y"),"y"))
    }else{
    Ax    <- function(t,x,y)  eval(driftx) - 0.5 * eval(diffx) * eval(D(diffx,"x"))
    dAx   <- function(t,x,y)  eval(D(driftx,"x")) - 0.5 * (eval(D(diffx,"x")) * eval(D(diffx,"x"))+ eval(diffx) * eval(D(D(diffx,"x"),"x")))
    dAxx  <- function(t,x,y)  eval(D(D(driftx,"x"),"x")) - 0.5 * ( eval(D(D(diffx,"x"),"x")) * eval(D(diffx,"x"))+ eval(D(diffx,"x")) * eval(D(D(diffx,"x"),"x"))+
                              eval(D(diffx,"x")) * eval(D(D(diffx,"x"),"x")) + eval(diffx) * eval(D(D(D(diffx,"x"),"x"),"x")) )
    Ay    <- function(t,x,y)  eval(drifty) - 0.5 * eval(diffy) * eval(D(diffy,"y"))
    dAy   <- function(t,x,y)  eval(D(drifty,"y")) - 0.5 * (eval(D(diffy,"y")) * eval(D(diffy,"y"))+ eval(diffy) * eval(D(D(diffy,"y"),"y")))
    dAyy  <- function(t,x,y)  eval(D(D(drifty,"y"),"y")) - 0.5 * ( eval(D(D(diffy,"y"),"y")) * eval(D(diffy,"y"))+ eval(D(diffy,"y")) * eval(D(D(diffy,"y"),"y"))+
                              eval(D(diffy,"y")) * eval(D(D(diffy,"y"),"y")) + eval(diffy) * eval(D(D(D(diffy,"y"),"y"),"y")) )
                  }
    DSx  <- D(diffx,"x")  
    DSxx <- D(DSx,"x")
    Sx    <- function(t,x,y)  eval(diffx)
    dSx   <- function(t,x,y)  eval(DSx)
    dSxx  <- function(t,x,y)  eval(DSxx)
    DSy  <- D(diffy,"y")  
    DSyy <- D(DSy,"y")
    Sy    <- function(t,x,y)  eval(diffy)
    dSy   <- function(t,x,y)  eval(DSy)
    dSyy  <- function(t,x,y)  eval(DSyy)
    wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    dwx   <- diff(wx)
    dwy   <- diff(wy)
    dzx= rnorm(N,0,sqrt((1/3)*Dt^3))
    dzy= rnorm(N,0,sqrt((1/3)*Dt^3))
    X = Y <- numeric()
    X[1] <- x0 
    Y[1] <- y0
    for (i in 2:(N+1)){
       X[i]=X[i-1]+Ax(t[i-1],X[i-1],Y[i-1])*Dt+Sx(t[i-1],X[i-1],Y[i-1])*dwx[i-1]+ 0.5*Sx(t[i-1],X[i-1],Y[i-1])*dSx(t[i-1],X[i-1],Y[i-1])*((dwx[i-1]^2)-Dt)+
            dAx(t[i-1],X[i-1],Y[i-1])*Sx(t[i-1],X[i-1],Y[i-1])*dzx[i-1]+0.5*(Ax(t[i-1],X[i-1],Y[i-1])*dAx(t[i-1],X[i-1],Y[i-1])+0.5*(Sx(t[i-1],X[i-1],Y[i-1])^2)*
            dAxx(t[i-1],X[i-1],Y[i-1]))*(Dt^2)+(Ax(t[i-1],X[i-1],Y[i-1])*dSx(t[i-1],X[i-1],Y[i-1])+0.5*(Sx(t[i-1],X[i-1],Y[i-1])^2)*dSxx(t[i-1],X[i-1],Y[i-1]))*
            (dwx[i-1]*Dt-dzx[i-1])+0.5*Sx(t[i-1],X[i-1],Y[i-1])*(Sx(t[i-1],X[i-1],Y[i-1])*dSxx(t[i-1],X[i-1],Y[i-1])+(dSx(t[i-1],X[i-1],Y[i-1])^2))*((1/3)*(dwx[i-1]^2)-Dt)*dwx[i-1]
       Y[i]=Y[i-1]+Ay(t[i-1],X[i-1],Y[i-1])*Dt+Sy(t[i-1],X[i-1],Y[i-1])*dwy[i-1]+ 0.5*Sy(t[i-1],X[i-1],Y[i-1])*dSy(t[i-1],X[i-1],Y[i-1])*((dwy[i-1]^2)-Dt)+
            dAy(t[i-1],X[i-1],Y[i-1])*Sy(t[i-1],X[i-1],Y[i-1])*dzy[i-1]+0.5*(Ay(t[i-1],X[i-1],Y[i-1])*dAy(t[i-1],X[i-1],Y[i-1])+0.5*(Sy(t[i-1],X[i-1],Y[i-1])^2)*
            dAyy(t[i-1],X[i-1],Y[i-1]))*(Dt^2)+(Ay(t[i-1],X[i-1],Y[i-1])*dSy(t[i-1],X[i-1],Y[i-1])+0.5*(Sy(t[i-1],X[i-1],Y[i-1])^2)*dSyy(t[i-1],X[i-1],Y[i-1]))*
            (dwy[i-1]*Dt-dzy[i-1])+0.5*Sy(t[i-1],X[i-1],Y[i-1])*(Sy(t[i-1],X[i-1],Y[i-1])*dSyy(t[i-1],X[i-1],Y[i-1])+(dSy(t[i-1],X[i-1],Y[i-1])^2))*((1/3)*(dwy[i-1]^2)-Dt)*dwy[i-1]
                      }
    res <- data.frame(X,Y)
    names(res) <- paste(c("X","Y"),sep="")
    X <- ts(res, start = t0, deltat = Dt)
    structure(list(X=X, driftx=driftx[[1]], diffx=diffx[[1]],drifty=drifty[[1]], diffy=diffy[[1]],
                  type=type, x0=x0,y0=y0, N=N, Dt=Dt,t0=t0,T=T),class="snssde2d")
} 


#####
##### STS3D

.STS3D <- function(N =100,x0=2,y0=1,z0=1,t0=0,T=1,Dt,driftx,diffx,drifty,diffy,
                     driftz,diffz,type=c("ito","str"),...)
                       {
    if (missing(Dt)) {
        t <- seq(t0, T, length = N + 1)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
        T <- t[N + 1]
    } 					   
    if (missing(type)) type <- "ito"  					   
    if (type=="ito"){
    Ax    <- function(t,x,y,z)  eval(driftx)
    dAx   <- function(t,x,y,z)  eval(D(driftx,"x"))
    dAxx  <- function(t,x,y,z)  eval(D(D(driftx,"x"),"x"))
    Ay    <- function(t,x,y,z)  eval(drifty)
    dAy   <- function(t,x,y,z)  eval(D(drifty,"y"))
    dAyy  <- function(t,x,y,z)  eval(D(D(drifty,"y"),"y"))
    Az    <- function(t,x,y,z)  eval(driftz)
    dAz   <- function(t,x,y,z)  eval(D(driftz,"z"))
    dAzz  <- function(t,x,y,z)  eval(D(D(driftz,"z"),"z"))}else{
    Ax    <- function(t,x,y,z)  eval(driftx) - 0.5 * eval(diffx) * eval(D(diffx,"x"))
    dAx   <- function(t,x,y,z)  eval(D(driftx,"x")) - 0.5 * (eval(D(diffx,"x")) * eval(D(diffx,"x"))+ eval(diffx) * eval(D(D(diffx,"x"),"x")))
    dAxx  <- function(t,x,y,z)  eval(D(D(driftx,"x"),"x")) - 0.5 * ( eval(D(D(diffx,"x"),"x")) * eval(D(diffx,"x"))+ eval(D(diffx,"x")) * eval(D(D(diffx,"x"),"x"))+
                                eval(D(diffx,"x")) * eval(D(D(diffx,"x"),"x")) + eval(diffx) * eval(D(D(D(diffx,"x"),"x"),"x")) )
    Ay    <- function(t,x,y,z)  eval(drifty) - 0.5 * eval(diffy) * eval(D(diffy,"y"))
    dAy   <- function(t,x,y,z)  eval(D(drifty,"y")) - 0.5 * (eval(D(diffy,"y")) * eval(D(diffy,"y"))+ eval(diffy) * eval(D(D(diffy,"y"),"y")))
    dAyy  <- function(t,x,y,z)  eval(D(D(drifty,"y"),"y")) - 0.5 * ( eval(D(D(diffy,"y"),"y")) * eval(D(diffy,"y"))+ eval(D(diffy,"y")) * eval(D(D(diffy,"y"),"y"))+
                                eval(D(diffy,"y")) * eval(D(D(diffy,"y"),"y")) + eval(diffy) * eval(D(D(D(diffy,"y"),"y"),"y")) )
    Az    <- function(t,x,y,z)  eval(driftz) - 0.5 * eval(diffz) * eval(D(diffz,"z"))
    dAz   <- function(t,x,y,z)  eval(D(driftz,"z")) - 0.5 * (eval(D(diffz,"z")) * eval(D(diffz,"z"))+ eval(diffz) * eval(D(D(diffz,"z"),"z")))
    dAzz  <- function(t,x,y,z)  eval(D(D(driftz,"z"),"z")) - 0.5 * ( eval(D(D(diffz,"z"),"z")) * eval(D(diffz,"z"))+ eval(D(diffz,"z")) * eval(D(D(diffz,"z"),"z"))+
                                eval(D(diffz,"z")) * eval(D(D(diffz,"z"),"z")) + eval(diffz) * eval(D(D(D(diffz,"z"),"z"),"z")) )
                  }
    DSx  <- D(diffx,"x")  
    DSxx <- D(DSx,"x")
    Sx    <- function(t,x,y,z)  eval(diffx)
    dSx   <- function(t,x,y,z)  eval(DSx)
    dSxx  <- function(t,x,y,z)  eval(DSxx)
    DSy  <- D(diffy,"y")  
    DSyy <- D(DSy,"y")
    Sy    <- function(t,x,y,z)  eval(diffy)
    dSy   <- function(t,x,y,z)  eval(DSy)
    dSyy  <- function(t,x,y,z)  eval(DSyy)
    DSz  <- D(diffz,"z")  
    DSzz <- D(DSz,"z")
    Sz    <- function(t,x,y,z)  eval(diffz)
    dSz   <- function(t,x,y,z)  eval(DSz)
    dSzz  <- function(t,x,y,z)  eval(DSzz)
    wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    wz = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    dwx   <- diff(wx)
    dwy   <- diff(wy)
    dwz   <- diff(wz)
    dzx= rnorm(N,0,sqrt((1/3)*Dt^3))
    dzy= rnorm(N,0,sqrt((1/3)*Dt^3))
    dzz= rnorm(N,0,sqrt((1/3)*Dt^3))
    X = Y = Z <- numeric()
    X[1] <- x0 
    Y[1] <- y0
    Z[1] <- z0
    for (i in 2:(N+1)){
       X[i]=X[i-1]+Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])*dwx[i-1]+ 0.5*Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])*dSx(t[i-1],X[i-1],Y[i-1],Z[i-1])*((dwx[i-1]^2)-Dt)+
            dAx(t[i-1],X[i-1],Y[i-1],Z[i-1])*Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])*dzx[i-1]+0.5*(Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])*dAx(t[i-1],X[i-1],Y[i-1],Z[i-1])+0.5*(Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])^2)*
            dAxx(t[i-1],X[i-1],Y[i-1],Z[i-1]))*(Dt^2)+(Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])*dSx(t[i-1],X[i-1],Y[i-1],Z[i-1])+0.5*(Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])^2)*dSxx(t[i-1],X[i-1],Y[i-1],Z[i-1]))*
            (dwx[i-1]*Dt-dzx[i-1])+0.5*Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])*(Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])*dSxx(t[i-1],X[i-1],Y[i-1],Z[i-1])+(dSx(t[i-1],X[i-1],Y[i-1],Z[i-1])^2))*((1/3)*(dwx[i-1]^2)-Dt)*dwx[i-1]
       Y[i]=Y[i-1]+Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])*dwy[i-1]+ 0.5*Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])*dSy(t[i-1],X[i-1],Y[i-1],Z[i-1])*((dwy[i-1]^2)-Dt)+
            dAy(t[i-1],X[i-1],Y[i-1],Z[i-1])*Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])*dzy[i-1]+0.5*(Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])*dAy(t[i-1],X[i-1],Y[i-1],Z[i-1])+0.5*(Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])^2)*
            dAyy(t[i-1],X[i-1],Y[i-1],Z[i-1]))*(Dt^2)+(Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])*dSy(t[i-1],X[i-1],Y[i-1],Z[i-1])+0.5*(Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])^2)*dSyy(t[i-1],X[i-1],Y[i-1],Z[i-1]))*
            (dwy[i-1]*Dt-dzy[i-1])+0.5*Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])*(Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])*dSyy(t[i-1],X[i-1],Y[i-1],Z[i-1])+(dSy(t[i-1],X[i-1],Y[i-1],Z[i-1])^2))*((1/3)*(dwy[i-1]^2)-Dt)*dwy[i-1]
       Z[i]=Z[i-1]+Az(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])*dwz[i-1]+ 0.5*Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])*dSz(t[i-1],X[i-1],Y[i-1],Z[i-1])*((dwz[i-1]^2)-Dt)+
            dAz(t[i-1],X[i-1],Y[i-1],Z[i-1])*Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])*dzz[i-1]+0.5*(Az(t[i-1],X[i-1],Y[i-1],Z[i-1])*dAz(t[i-1],X[i-1],Y[i-1],Z[i-1])+0.5*(Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])^2)*
            dAzz(t[i-1],X[i-1],Y[i-1],Z[i-1]))*(Dt^2)+(Az(t[i-1],X[i-1],Y[i-1],Z[i-1])*dSz(t[i-1],X[i-1],Y[i-1],Z[i-1])+0.5*(Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])^2)*dSzz(t[i-1],X[i-1],Y[i-1],Z[i-1]))*
            (dwz[i-1]*Dt-dzz[i-1])+0.5*Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])*(Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])*dSzz(t[i-1],X[i-1],Y[i-1],Z[i-1])+(dSz(t[i-1],X[i-1],Y[i-1],Z[i-1])^2))*((1/3)*(dwz[i-1]^2)-Dt)*dwz[i-1]
                      }
    res <- data.frame(X,Y,Z)
    names(res) <- paste(c("X","Y","Z"),sep="")
    X <- ts(res, start = t0, deltat = Dt)
    structure(list(X=X, driftx=driftx[[1]], diffx=diffx[[1]],drifty=drifty[[1]], diffy=diffy[[1]],
                  driftz=driftz[[1]], diffz=diffz[[1]],type=type, x0=x0,y0=y0,z0=z0, N=N, Dt=Dt,t0=t0,T=T),class="snssde3d")
} 




