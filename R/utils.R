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




###
### Computes the bound of the confidence

bconfint <- function(x, ...)  UseMethod("bconfint")

bconfint.default <- function(x,level = 0.95,...)
            {
   return(quantile(x,c(0.5*(1-level), 1-0.5*(1-level)),type=8))
         }
###

add.bconfint <- function(x,...) UseMethod("add.bconfint")

add.bconfint.default <- function(x,level = 0.95,...)
            {
    lines(time(x),bconfint(x,level)[,1],...)
    lines(time(x),bconfint(x,level)[,2],...)
         }

###
### Computes the moment

moment <- function(x, ...)  UseMethod("moment")

moment.default <- function(x, order = 2,...)
            {
    if (any(!is.numeric(order)  || (order - floor(order) > 0) || order < 1)) stop(" 'order' must be a positive integer ")
    ifelse(order == 1,return(0),return(mean((x - mean(x))^order)))
         }

###
### Computes the sample skewness

skewness <- function(x,...) UseMethod("skewness")

skewness.default <- function(x,...)
                   {
     return(mean((x-mean(x))^3)/sd(x)^3)
}

###
### Computes the sample kurtosis

kurtosis <- function(x,...) UseMethod("kurtosis")

kurtosis.default <- function(x,...)
                   {
     return(mean((x-mean(x))^4)/sd(x)^4)
}

###
### Plot

plot2d   <- function(x,...) UseMethod("plot2d")
lines2d  <- function(x,...) UseMethod("lines2d")
points2d <- function(x,...) UseMethod("points2d")
plot3D   <- function(x, ...)  UseMethod("plot3D")

plot2d.default <- function(x,...)
        {
    class(x) <- "plot2d"
    plot(x,...)
}

lines2d.default <- function(x,...)
        {
    class(x) <- "lines2d"
    lines(x,...)
}

points2d.default <- function(x,...)
        {
    class(x) <- "points2d"
    points(x,...)
}


plot3DD.default <- function(X,Y,Z,display = c("persp","rgl"),col=NULL,lwd=NULL,pch=NULL,
                            type = NULL,cex=NULL,main=NULL,sub=NULL,xlab=NULL,ylab=NULL,
                            zlab=NULL,grid=NULL,angle=NULL,...)
                 {
    display <- match.arg(display)
    if ((display == "rgl") && !(require(rgl)) ) 
                 stop("The 'rgl' package is not available.")        
    else if ((display == "persp") && !(require(scatterplot3d)) ) 
                 stop("The 'scatterplot3d' package is not available.")

    if (is.null(lwd))  {lwd = 1}
    if (is.null(col))  {col = 2}
    if (is.null(pch))  {pch = 16}
    if (is.null(cex))  {cex = 0.6}
    if (is.null(type)) {type = "l"}
    if (is.null(main)) {main = ""}
    if (is.null(sub))  {sub = ""}
    if (is.null(xlab)) {xlab = expression(X[t])}
    if (is.null(ylab)) {ylab = expression(Y[t])}
    if (is.null(zlab)) {zlab = expression(Z[t])} 
    if (is.null(grid)) {grid = TRUE}   
    if (is.null(angle)){angle = 140} 
    if (display=="persp"){
         scatterplot3d(X,Y,Z,angle =angle,color=col,lwd=lwd,type=type,pch=pch,
             main = main, sub = sub,xlab = xlab, ylab = ylab, zlab = zlab,
             grid = grid,cex.symbols=cex,...)
    }else{
         plot3d(X,Y,Z,col=col,lwd=lwd,type=type,pch=pch,main = main, sub = sub,
             xlab = xlab, ylab = ylab, zlab = zlab,size=cex,...)
         }
}

plot3D.default <- function(x,display = c("persp", "rgl"),...) plot3DD.default(x,display,...)

####
#### Char2expression

.Char2exp <- function(expr)
         {
y <- gsub(pattern = 'theta[1]', replacement = 'theta1', x = expr, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[2]', replacement = 'theta2', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[3]', replacement = 'theta3', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[4]', replacement = 'theta4', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[5]', replacement = 'theta5', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[6]', replacement = 'theta6', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[7]', replacement = 'theta7', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[8]', replacement = 'theta8', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[9]', replacement = 'theta9', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[10]', replacement = 'theta10', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[11]', replacement = 'theta11', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[12]', replacement = 'theta12', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[13]', replacement = 'theta13', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[14]', replacement = 'theta14', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[15]', replacement = 'theta15', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[16]', replacement = 'theta16', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[17]', replacement = 'theta17', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[18]', replacement = 'theta18', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[19]', replacement = 'theta19', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[20]', replacement = 'theta20', x = y, ignore.case = F,fixed = T)
Y <- parse(text=y,srcfile=y)
return(Y)
}


.Exp2char <- function(expr)
         {
y <- gsub(pattern = 'theta1', replacement = 'theta[1]', x = expr, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta2', replacement = 'theta[2]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta3', replacement = 'theta[3]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta4', replacement = 'theta[4]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta5', replacement = 'theta[5]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta6', replacement = 'theta[6]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta7', replacement = 'theta[7]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta8', replacement = 'theta[8]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta9', replacement = 'theta[9]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta10', replacement = 'theta[10]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta11', replacement = 'theta[11]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta12', replacement = 'theta[12]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta13', replacement = 'theta[13]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta14', replacement = 'theta[14]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta15', replacement = 'theta[15]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta16', replacement = 'theta[16]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta17', replacement = 'theta[17]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta18', replacement = 'theta[18]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta19', replacement = 'theta[19]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta20', replacement = 'theta[20]', x = y, ignore.case = F,fixed = T)
Y <- parse(text=y,srcfile=y)
return(Y)
}

#########
#########
#####
##### dc.sde1d

.dcEuler <- function(x,t,x0,t0,theta,drift,diffusion,log=FALSE)
              {
	 A   <- function(t,x,theta)  eval(drift)
     S   <- function(t,x,theta)  eval(diffusion)
     dt  <- t-t0
     lik <- dnorm(x,mean=x0+A(t0,x0,theta)*dt,sd=sqrt(dt)*S(t0,x0,theta),log=log)
     lik[is.infinite(lik)] <- NA
     lik
       }
	   
.dcOzaki <- function(x, t, x0, t0, theta, drift,diffusion, log=FALSE)
             {  
     A   <- function(t,x,theta)  eval(drift)
     S   <- function(t,x,theta)  eval(diffusion)
     Ax  <- function(t,x,theta)  eval(.Exp2char(as.expression(D(.Char2exp(drift),"x"))))
     dt  <- t-t0
     K   <- log(1+(A(t0,x0,theta)/(x0*Ax(t0,x0,theta)))*(exp(Ax(t0,x0,theta)*dt)-1))/dt
     E   <- x0 + (A(t0,x0,theta)/Ax(t0,x0,theta))*(exp(Ax(t0,x0,theta)*dt)-1)
     V   <- S(t0,x0,theta)^2 * ((exp(2*K*dt) -1)/(2*K))
     lik <- dnorm(x, mean=E, sd=sqrt(V),log=log)
     lik[is.infinite(lik)] <- NA
     lik      
}

.dcShoji <- function(x, t, x0, t0, theta, drift,diffusion, log=FALSE)
            {
     A   <- function(t,x,theta)  eval(drift)
     S   <- function(t,x,theta)  eval(diffusion)
	 At  <- function(t,x,theta)  eval(.Exp2char(as.expression(D(.Char2exp(drift),"t"))))
	 Ax  <- function(t,x,theta)  eval(.Exp2char(as.expression(D(.Char2exp(drift),"x"))))
	 Axx <- function(t,x,theta)  eval(.Exp2char(as.expression(D(D(.Char2exp(drift),"x"),"x"))))
     dt  <- t-t0
     E   <- x0 + A(t0,x0,theta)*(exp(Ax(t0,x0,theta)*dt)-1)/Ax(t0,x0,theta) + (0.5 * S(t0,x0,theta)^2 * Axx(t0,x0,theta) + At(t0,x0,theta))*(exp(Ax(t0,x0,theta)*dt)-1-Ax(t0,x0,theta)*dt)/Ax(t0,x0,theta)^2 
     V   <- S(t0,x0,theta)^2 * (exp(2*Ax(t0,x0,theta)*dt)-1)/(2*Ax(t0,x0,theta))
     lik <- dnorm(x, mean=E, sd=sqrt(V),log=log) 
	 lik[is.infinite(lik)] <- NA
     lik   
}

# .dcElerian <- function(x, t, x0, t0, theta, drift,diffusion, log=FALSE)
            # {
     # test <- .Exp2char(as.expression(D(.Char2exp(diffusion),"x")))[[1]]
     # if (test == 0) stop("The approximation is not valid, because 'deriv(diffusion,x) = 0'")
     # A   <- function(t,x,theta)  eval(drift)
     # S   <- function(t,x,theta)  eval(diffusion)
	 # Sx  <- function(t,x,theta)  eval(.Exp2char(as.expression(D(.Char2exp(diffusion),"x"))))	 
     # dt <- t-t0
     # E <- 0.5*S(t0, x0, theta)*Sx(t0, x0, theta)*dt
     # B <- -0.5* (S(t0, x0,theta)/Sx(t0, x0,theta)) + x0 + A(t0, x0, theta)*dt - E
     # C <- 1/((S(t0, x0,theta)^2) * dt)
     # z <- (x-B)/E
     # z[z<0] <- NA
     # lik <- ( exp(-0.5*(C+z)) * cosh(sqrt(C*z)) )/( sqrt(z) *abs(E)* sqrt(2*pi) )	 
     # if(log) lik <- log(lik)
     # lik[is.infinite(lik)] <- NA
     # lik
# }

.dcKessler <- function(x, t, x0, t0, theta, drift,diffusion, log=FALSE)
           {
     A   <- function(t,x,theta)  eval(drift)
	 Ax  <- function(t,x,theta)  eval(.Exp2char(as.expression(D(.Char2exp(drift),"x"))))
	 Axx <- function(t,x,theta)  eval(.Exp2char(as.expression(D(D(.Char2exp(drift),"x"),"x"))))	
     S   <- function(t,x,theta)  eval(diffusion)
	 Sx  <- function(t,x,theta)  eval(.Exp2char(as.expression(D(.Char2exp(diffusion),"x"))))
	 Sxx <- function(t,x,theta)  eval(.Exp2char(as.expression(D(D(.Char2exp(diffusion),"x"),"x"))))		   
     dt <- t-t0 
     E   <- x0 + A(t0,x0,theta)*dt + 0.5*(dt)^2 * (A(t0,x0,theta)*Ax(t0,x0,theta)+0.5*(S(t0,x0,theta))^2 * Axx(t0,x0,theta))
     V   <- x0^2 + (2*A(t0,x0,theta)*x0+(S(t0,x0,theta))^2)*dt + 0.5*(dt)^2 *(2*A(t0,x0,theta)*(Ax(t0,x0,theta)*x0+A(t0,x0,theta)+S(t0,x0,theta)*Sx(t0,x0,theta))+(S(t0,x0,theta))^2 *(Axx(t0,x0,theta)*x0+2*Ax(t0,x0,theta)+(Sx(t0,x0,theta))^2 + S(t0,x0,theta)*Sxx(t0,x0,theta))) - E^2
	 V[V < 0] <- NA
     lik <- dnorm(x, mean=E, sd=sqrt(V),log=log) 
	 lik[is.infinite(lik)] <- NA
     lik 
}


.onLoad <- function(libname, pkgname)
   packageStartupMessage("Package 'Sim.DiffProc' version 2.6 loaded. help(Sim.DiffProc) for summary information.\nA.C. Guidoum and K. Boukhetala.\n")

