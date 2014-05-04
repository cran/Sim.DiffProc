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
##### RK1D

.RK1D <- function(N =100,M=1,x0=2,t0=0,T=1,Dt,drift,diffusion,
                          type=c("ito","str"),order=c(1,2,3),...)
                       { 
    if (missing(Dt)) {
        t <- seq(t0, T, length = N + 1)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
        T <- t[N + 1]
    } 					   
    if (missing(type)) type <- "ito" 					   
    if (type=="ito"){A    <- function(t,x)  eval(drift)}else{
    A  <- function(t,x) eval(drift) - 0.5 * eval(diffusion) * eval(D(diffusion,"x"))}
    S  <- function(t,x) eval(diffusion) 
    rk1d <- function()
           {
    w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    dw   <- diff(w)
    X    <- Y <- Z <- numeric()
    X[1] <- x0 
    for (i in 2:(N+1)){
    if (order==1){
      X[i] = X[i-1]+ A(t[i-1],X[i-1])*Dt+S(t[i-1],X[i-1])*dw[i-1]+(0.5/sqrt(Dt)) * (S(t[i-1]+Dt,X[i-1]+S(t[i-1],X[i-1])*sqrt(Dt))-S(t[i-1],X[i-1])) * (dw[i-1]^2 - Dt)}
    else if (order==2){
      X[i] = X[i-1]+ 0.5 * (A(t[i-1],X[i-1])+A(t[i-1]+Dt,X[i-1]+A(t[i-1],X[i-1])*Dt+S(t[i-1],X[i-1])*dw[i-1])) *Dt + 0.25* (2*S(t[i-1],X[i-1])+S(t[i-1]+Dt,X[i-1]+A(t[i-1],X[i-1])*Dt+
             sqrt(Dt)*S(t[i-1],X[i-1]))+S(t[i-1]+Dt,X[i-1]+A(t[i-1],X[i-1])*Dt-sqrt(Dt)*S(t[i-1],X[i-1])) ) *dw[i-1]+0.25 * (S(t[i-1]+Dt,X[i-1]+A(t[i-1],X[i-1])*Dt-sqrt(Dt)*S(t[i-1],X[i-1])) - 
             S(t[i-1]+Dt,X[i-1]+A(t[i-1],X[i-1])*Dt+sqrt(Dt)*S(t[i-1],X[i-1]))) *(sqrt(Dt) - (dw[i-1])^2 / sqrt(Dt))}
    else if (order==3){
              Y[i-1]=X[i-1]+0.5*Dt*A(t[i-1],X[i-1])+S(t[i-1],X[i-1])*dw[i-1]
              Z[i-1]=X[i-1]-A(t[i-1],X[i-1])*Dt+2*Dt*A(t[i-1]+0.5*Dt,Y[i-1])+(2*S(t[i-1]+0.5*Dt,Y[i-1])-S(t[i-1],X[i-1]))*dw[i-1]
              X[i] = X[i-1]+(Dt/6)*(A(t[i-1],X[i-1])+4*A(t[i-1]+0.5*Dt,Y[i-1])+A(t[i-1]+Dt,Z[i-1]))+(1/6)*(S(t[i-1],X[i-1])+4*S(t[i-1]+0.5*Dt,Y[i-1])+S(t[i-1]+Dt,Z[i-1]))*dw[i-1]}
       }
    X
         }
    res <- data.frame(replicate(M,rk1d()))
    names(res) <- paste("X",1:M,sep="")
    X <- ts(res, start = t0, deltat = Dt)
    structure(list(X=X, drift=drift[[1]], diffusion=diffusion[[1]],type=type, x0=x0, N=N, M = M, 
                   Dt=Dt,t0=t0,T=T),class="snssde1d")
}           

#####
##### RK2D

.RK2D <- function(N =100,x0=2,y0=1,t0=0,T=1,Dt,driftx,diffx,drifty,diffy,
                          type=c("ito","str"),order=c(1,2,3),...)
                       {
    if (missing(Dt)) {
        t <- seq(t0, T, length = N + 1)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
        T <- t[N + 1]
    } 					   
    if (missing(type)) type <- "ito" 					   
    if (type=="ito"){
    Ax <- function(t,x,y)  eval(driftx)
    Ay <- function(t,x,y)  eval(drifty) }else{
    Ax <- function(t,x,y) eval(driftx) - 0.5 * eval(diffx) * eval(D(diffx,"x"))
    Ay <- function(t,x,y) eval(drifty) - 0.5 * eval(diffy) * eval(D(diffy,"y"))
                         }
    Sx <- function(t,x,y) eval(diffx)
    Sy <- function(t,x,y) eval(diffy) 
    wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    dwx   <- diff(wx)
    dwy   <- diff(wy)
    X = Y <- numeric()
    XX = YY = XXX = YYY   <- numeric()
    X[1] <- x0 
    Y[1] <- y0
    for (i in 2:(N+1)){
    if (order==1){
      X[i] = X[i-1]+ Ax(t[i-1],X[i-1],Y[i-1])*Dt+Sx(t[i-1],X[i-1],Y[i-1])*dwx[i-1]+(0.5/sqrt(Dt)) * (Sx(t[i-1]+Dt,X[i-1]+Sx(t[i-1],X[i-1],Y[i-1])*sqrt(Dt),Y[i-1])-Sx(t[i-1],X[i-1],Y[i-1])) * (dwx[i-1]^2 - Dt)
      Y[i] = Y[i-1]+ Ay(t[i-1],X[i-1],Y[i-1])*Dt+Sy(t[i-1],X[i-1],Y[i-1])*dwy[i-1]+(0.5/sqrt(Dt)) * (Sy(t[i-1]+Dt,X[i-1],Y[i-1]+Sy(t[i-1],X[i-1],Y[i-1])*sqrt(Dt))-Sy(t[i-1],X[i-1],Y[i-1])) * (dwy[i-1]^2 - Dt)
	  }
    else if (order==2){
      X[i] = X[i-1]+ 0.5 * (Ax(t[i-1],X[i-1],Y[i-1])+Ax(t[i-1]+Dt,X[i-1]+Ax(t[i-1],X[i-1],Y[i-1])*Dt+Sx(t[i-1],X[i-1],Y[i-1])*dwx[i-1],Y[i-1])) *Dt + 
			  0.25* (2*Sx(t[i-1],X[i-1],Y[i-1])+Sx(t[i-1]+Dt,X[i-1]+Ax(t[i-1],X[i-1],Y[i-1])*Dt+sqrt(Dt)*Sx(t[i-1],X[i-1],Y[i-1]),Y[i-1])+Sx(t[i-1]+Dt,X[i-1]+Ax(t[i-1],X[i-1],Y[i-1])*Dt-sqrt(Dt)*Sx(t[i-1],X[i-1],Y[i-1]),Y[i-1]) ) *dwx[i-1]+
              0.25 * (Sx(t[i-1]+Dt,X[i-1]+Ax(t[i-1],X[i-1],Y[i-1])*Dt-sqrt(Dt)*Sx(t[i-1],X[i-1],Y[i-1]),Y[i-1]) - Sx(t[i-1]+Dt,X[i-1]+Ax(t[i-1],X[i-1],Y[i-1])*Dt+sqrt(Dt)*Sx(t[i-1],X[i-1],Y[i-1]),Y[i-1])) *(sqrt(Dt) - (dwx[i-1])^2 / sqrt(Dt))
      Y[i] = Y[i-1]+ 0.5 * (Ay(t[i-1],X[i-1],Y[i-1])+Ay(t[i-1]+Dt,X[i-1],Y[i-1]+Ay(t[i-1],X[i-1],Y[i-1])*Dt+Sy(t[i-1],X[i-1],Y[i-1])*dwy[i-1])) *Dt + 
			  0.25* (2*Sy(t[i-1],X[i-1],Y[i-1])+Sy(t[i-1]+Dt,X[i-1],Y[i-1]+Ay(t[i-1],X[i-1],Y[i-1])*Dt+sqrt(Dt)*Sy(t[i-1],X[i-1],Y[i-1]))+Sy(t[i-1]+Dt,X[i-1],Y[i-1]+Ay(t[i-1],X[i-1],Y[i-1])*Dt-sqrt(Dt)*Sy(t[i-1],X[i-1],Y[i-1])) ) *dwy[i-1]+
              0.25 * (Sy(t[i-1]+Dt,X[i-1],Y[i-1]+Ay(t[i-1],X[i-1],Y[i-1])*Dt-sqrt(Dt)*Sy(t[i-1],X[i-1],Y[i-1]))-Sy(t[i-1]+Dt,X[i-1],Y[i-1]+Ay(t[i-1],X[i-1],Y[i-1])*Dt+sqrt(Dt)*Sy(t[i-1],X[i-1],Y[i-1]))) *(sqrt(Dt) - (dwy[i-1])^2 / sqrt(Dt))
			  }
    else if (order==3){
      XX[i-1] =X[i-1]+0.5*Dt*Ax(t[i-1],X[i-1],Y[i-1])+Sx(t[i-1],X[i-1],Y[i-1])*dwx[i-1]
      YY[i-1] =Y[i-1]+0.5*Dt*Ay(t[i-1],X[i-1],Y[i-1])+Sy(t[i-1],X[i-1],Y[i-1])*dwy[i-1]
      XXX[i-1]=X[i-1]-Ax(t[i-1],X[i-1],Y[i-1])*Dt+2*Dt*Ax(t[i-1]+0.5*Dt,XX[i-1],Y[i-1])+(2*Sx(t[i-1]+0.5*Dt,XX[i-1],Y[i-1])-Sx(t[i-1],X[i-1],Y[i-1]))*dwx[i-1]
      YYY[i-1]=Y[i-1]-Ay(t[i-1],X[i-1],Y[i-1])*Dt+2*Dt*Ay(t[i-1]+0.5*Dt,X[i-1],YY[i-1])+(2*Sy(t[i-1]+0.5*Dt,X[i-1],YY[i-1])-Sy(t[i-1],X[i-1],Y[i-1]))*dwy[i-1]
      X[i] = X[i-1]+(Dt/6)*(Ax(t[i-1],X[i-1],Y[i-1])+4*Ax(t[i-1]+0.5*Dt,XX[i-1],Y[i-1])+Ax(t[i-1]+Dt,XXX[i-1],Y[i-1]))+(1/6)*(Sx(t[i-1],X[i-1],Y[i-1])+4*Sx(t[i-1]+0.5*Dt,XX[i-1],Y[i-1])+Sx(t[i-1]+Dt,XXX[i-1],Y[i-1]))*dwx[i-1]
      Y[i] = Y[i-1]+(Dt/6)*(Ay(t[i-1],X[i-1],Y[i-1])+4*Ay(t[i-1]+0.5*Dt,X[i-1],YY[i-1])+Ay(t[i-1]+Dt,X[i-1],YYY[i-1]))+(1/6)*(Sy(t[i-1],X[i-1],Y[i-1])+4*Sy(t[i-1]+0.5*Dt,X[i-1],YY[i-1])+Sy(t[i-1]+Dt,X[i-1],YYY[i-1]))*dwy[i-1]
			  }
       }
    res <- data.frame(X,Y)
    names(res) <- paste(c("X","Y"),sep="")
    X <- ts(res, start = t0, deltat = Dt)
    structure(list(X=X, driftx=driftx[[1]], diffx=diffx[[1]],drifty=drifty[[1]], diffy=diffy[[1]],
                  type=type, x0=x0,y0=y0, N=N, Dt=Dt,t0=t0,T=T),class="snssde2d")
} 


#####
##### RK3D

.RK3D <- function(N =100,x0=2,y0=1,z0=1,t0=0,T=1,Dt,driftx,diffx,drifty,diffy,
                     driftz,diffz,type=c("ito","str"),order=c(1,2,3),...)
                       { 
    if (missing(Dt)) {
        t <- seq(t0, T, length = N + 1)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
        T <- t[N + 1]
    } 					 
    if (missing(type)) type <- "ito" 					   
    if (type=="ito"){
    Ax <- function(t,x,y,z)  eval(driftx)
    Ay <- function(t,x,y,z)  eval(drifty) 
    Az <- function(t,x,y,z)  eval(driftz)}else{
    Ax <- function(t,x,y,z)  eval(driftx) - 0.5 * eval(diffx) * eval(D(diffx,"x"))
    Ay <- function(t,x,y,z)  eval(drifty) - 0.5 * eval(diffy) * eval(D(diffy,"y"))
    Az <- function(t,x,y,z)  eval(driftz) - 0.5 * eval(diffz) * eval(D(diffz,"z"))
                         }
    Sx <- function(t,x,y,z) eval(diffx)
    Sy <- function(t,x,y,z) eval(diffy) 
    Sz <- function(t,x,y,z) eval(diffz)
    wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    wz = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
    dwx   <- diff(wx)
    dwy   <- diff(wy)
    dwz   <- diff(wz)
    X = Y = Z <- numeric()
    XX = YY = ZZ = XXX = YYY = ZZZ <- numeric()
    X[1] <- x0 
    Y[1] <- y0
    Z[1] <- z0
    for (i in 2:(N+1)){
    if (order==1){
      X[i] = X[i-1]+ Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])*dwx[i-1]+(0.5/sqrt(Dt)) * (Sx(t[i-1]+Dt,X[i-1]+Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])*sqrt(Dt),Y[i-1],Z[i-1])-Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])) * (dwx[i-1]^2 - Dt)
	  Y[i] = Y[i-1]+ Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])*dwy[i-1]+(0.5/sqrt(Dt)) * (Sy(t[i-1]+Dt,X[i-1],Y[i-1]+Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])*sqrt(Dt),Z[i-1])-Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])) * (dwy[i-1]^2 - Dt)			  
	  Z[i] = Z[i-1]+ Az(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])*dwz[i-1]+(0.5/sqrt(Dt)) * (Sz(t[i-1]+Dt,X[i-1],Y[i-1],Z[i-1]+Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])*sqrt(Dt))-Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])) * (dwz[i-1]^2 - Dt)
	  }
    else if (order==2){
      X[i] = X[i-1]+ 0.5 * (Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])+Ax(t[i-1]+Dt,X[i-1]+Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])*dwx[i-1],Y[i-1],Z[i-1])) *Dt + 
			 0.25* (2*Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])+Sx(t[i-1]+Dt,X[i-1]+Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+sqrt(Dt)*Sx(t[i-1],X[i-1],Y[i-1],Z[i-1]),Y[i-1],Z[i-1])+ Sx(t[i-1]+Dt,X[i-1]+Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt-sqrt(Dt)*Sx(t[i-1],X[i-1],Y[i-1],Z[i-1]),Y[i-1],Z[i-1]) ) *dwx[i-1]+
             0.25 * (Sx(t[i-1]+Dt,X[i-1]+Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt-sqrt(Dt)*Sx(t[i-1],X[i-1],Y[i-1],Z[i-1]),Y[i-1],Z[i-1]) - Sx(t[i-1]+Dt,X[i-1]+Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+sqrt(Dt)*Sx(t[i-1],X[i-1],Y[i-1],Z[i-1]),Y[i-1],Z[i-1])) *(sqrt(Dt) - (dwx[i-1])^2 / sqrt(Dt))              
	  Y[i] = Y[i-1]+ 0.5 * (Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])+Ay(t[i-1]+Dt,X[i-1],Y[i-1]+Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])*dwy[i-1],Z[i-1])) *Dt + 
			 0.25* (2*Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])+Sy(t[i-1]+Dt,X[i-1],Y[i-1]+Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+sqrt(Dt)*Sy(t[i-1],X[i-1],Y[i-1],Z[i-1]),Z[i-1])+Sy(t[i-1]+Dt,X[i-1],Y[i-1]+Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt-sqrt(Dt)*Sy(t[i-1],X[i-1],Y[i-1],Z[i-1]),Z[i-1]) ) *dwy[i-1]+
             0.25 * (Sy(t[i-1]+Dt,X[i-1],Y[i-1]+Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt-sqrt(Dt)*Sy(t[i-1],X[i-1],Y[i-1],Z[i-1]),Z[i-1])-Sy(t[i-1]+Dt,X[i-1],Y[i-1]+Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+sqrt(Dt)*Sy(t[i-1],X[i-1],Y[i-1],Z[i-1]),Z[i-1])) *(sqrt(Dt) - (dwy[i-1])^2 / sqrt(Dt))              
	  Z[i] = Z[i-1]+ 0.5 * (Az(t[i-1],X[i-1],Y[i-1],Z[i-1])+Az(t[i-1]+Dt,X[i-1],Y[i-1],Z[i-1]+Az(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])*dwz[i-1])) *Dt + 
			 0.25* (2*Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])+Sz(t[i-1]+Dt,X[i-1],Y[i-1],Z[i-1]+Az(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+sqrt(Dt)*Sz(t[i-1],X[i-1],Y[i-1],Z[i-1]))+Sz(t[i-1]+Dt,X[i-1],Y[i-1],Z[i-1]+Az(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt-sqrt(Dt)*Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])) ) *dwz[i-1]+
             0.25 * (Sz(t[i-1]+Dt,X[i-1],Y[i-1],Z[i-1]+Az(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt-sqrt(Dt)*Sz(t[i-1],X[i-1],Y[i-1],Z[i-1]))-Sz(t[i-1]+Dt,X[i-1],Y[i-1],Z[i-1]+Az(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+sqrt(Dt)*Sz(t[i-1],X[i-1],Y[i-1],Z[i-1]))) *(sqrt(Dt) - (dwz[i-1])^2 / sqrt(Dt))  
			  }
    else if (order==3){
      XX[i-1]=X[i-1]+0.5*Dt*Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])+Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])*dwx[i-1]
      YY[i-1]=Y[i-1]+0.5*Dt*Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])+Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])*dwy[i-1]
      ZZ[i-1]=Z[i-1]+0.5*Dt*Az(t[i-1],X[i-1],Y[i-1],Z[i-1])+Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])*dwz[i-1]
      XXX[i-1]=X[i-1]-Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+2*Dt*Ax(t[i-1]+0.5*Dt,XX[i-1],Y[i-1],Z[i-1])+(2*Sx(t[i-1]+0.5*Dt,XX[i-1],Y[i-1],Z[i-1])-Sx(t[i-1],X[i-1],Y[i-1],Z[i-1]))*dwx[i-1]
      YYY[i-1]=Y[i-1]-Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+2*Dt*Ay(t[i-1]+0.5*Dt,X[i-1],YY[i-1],Z[i-1])+(2*Sy(t[i-1]+0.5*Dt,X[i-1],YY[i-1],Z[i-1])-Sy(t[i-1],X[i-1],Y[i-1],Z[i-1]))*dwy[i-1]
      ZZZ[i-1]=Z[i-1]-Az(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+2*Dt*Az(t[i-1]+0.5*Dt,X[i-1],Y[i-1],ZZ[i-1])+(2*Sz(t[i-1]+0.5*Dt,X[i-1],Y[i-1],ZZ[i-1])-Sz(t[i-1],X[i-1],Y[i-1],Z[i-1]))*dwz[i-1]
      X[i] = X[i-1]+(Dt/6)*(Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])+4*Ax(t[i-1]+0.5*Dt,XX[i-1],Y[i-1],Z[i-1])+Ax(t[i-1]+Dt,XXX[i-1],Y[i-1],Z[i-1]))+
             (1/6)*(Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])+4*Sx(t[i-1]+0.5*Dt,XX[i-1],Y[i-1],Z[i-1])+Sx(t[i-1]+Dt,XXX[i-1],Y[i-1],Z[i-1]))*dwx[i-1]
      Y[i] = Y[i-1]+(Dt/6)*(Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])+4*Ay(t[i-1]+0.5*Dt,X[i-1],YY[i-1],Z[i-1])+Ay(t[i-1]+Dt,X[i-1],YYY[i-1],Z[i-1]))+
             (1/6)*(Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])+4*Sy(t[i-1]+0.5*Dt,X[i-1],YY[i-1],Z[i-1])+Sy(t[i-1]+Dt,X[i-1],YYY[i-1],Z[i-1]))*dwy[i-1]
      Z[i] = Z[i-1]+(Dt/6)*(Az(t[i-1],X[i-1],Y[i-1],Z[i-1])+4*Az(t[i-1]+0.5*Dt,X[i-1],Y[i-1],ZZ[i-1])+Az(t[i-1]+Dt,X[i-1],Y[i-1],ZZZ[i-1]))+
             (1/6)*(Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])+4*Sz(t[i-1]+0.5*Dt,X[i-1],Y[i-1],ZZ[i-1])+Sz(t[i-1]+Dt,X[i-1],Y[i-1],ZZZ[i-1]))*dwz[i-1]
			  }
       }
    res <- data.frame(X,Y,Z)
    names(res) <- paste(c("X","Y","Z"),sep="")
    X <- ts(res, start = t0, deltat = Dt)
    structure(list(X=X, driftx=driftx[[1]], diffx=diffx[[1]],drifty=drifty[[1]], diffy=diffy[[1]],
                  driftz=driftz[[1]], diffz=diffz[[1]],type=type, x0=x0,y0=y0,z0=z0, N=N, Dt=Dt,t0=t0,T=T),class="snssde3d")
} 



