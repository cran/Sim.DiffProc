\name{MEM.sde}
\alias{MEM.sde}
\alias{MEM.sde.default}
\alias{print.MEM.sde}
\alias{summary.MEM.sde}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Moment Equations Methods for SDE's
}
\description{
Calculate and numerical approximation of moment equations (Symbolic ODE's of means and variances-covariance) at any time for SDE's (1,2 and 3 dim) for the two cases Ito and Stratonovich interpretations.
}
\usage{
MEM.sde(drift, diffusion, \dots)

\method{MEM.sde}{default}(drift, diffusion, corr = NULL, type = c("ito", "str"), solve = FALSE, 
      parms = NULL, init = NULL, time = NULL, ...)
					
\method{summary}{MEM.sde}(object, at , \dots)					
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{drift}{
drift coefficient: an \code{\link[base]{expression}} 1-dim(t,x), 2-dim(t,x,y) or 3-dim(t,x,y,z).
}
  \item{diffusion}{
diffusion coefficient: an \code{\link[base]{expression}} 1-dim(t,x), 2-dim(t,x,y) or 3-dim(t,x,y,z).
}
\item{corr}{
the correlation coefficient '|corr|<=1' of W1(t) and W2(t) (2d) must be an expression length equal 1. And for 3d (W1(t),W2(t),W3(t)) an expressions length equal 3. See examples.
}
\item{type}{
type of process \code{"ito"} or \code{"Stratonovich"}; the default \code{type="ito"}.
}
  \item{solve}{
if \code{solve=TRUE} solves a system of ordinary differential equations.
}
  \item{parms}{
parameters passed to \code{drift} and \code{diffusion}.
}
  \item{init}{
the initial (state) values for the ODE system. for 1-dim (m=x0,S=0), 2-dim (m1=x0,m2=y0,S1=0,S2=0,C12=0) and 
for 3-dim (m1=x0,m2=y0,m3=z0,S1=0,S2=0,S3=0,C12=0,C13=0,C23=0), see examples.
}
  \item{time}{	
time sequence (vector) for which output is wanted; the first value of time must be the initial time.
}
  \item{object, at}{
  an object inheriting from class \code{"MEM.sde"} and summaries at any time \code{at}.
}
  \item{\dots}{
potentially arguments to be passed to methods, such as \code{\link[deSolve]{ode}} for solver for ODE's.
}
}
\details{
The stochastic transition is approximated by the moment equations, and the numerical treatment is required to solve these equations from above with given initial conditions.

An overview of this package, see \code{\link[utils]{browseVignettes}('Sim.DiffProc')} for more informations.
}
\value{
Symbolic ODE's of means and variances-covariance. If \code{solve=TRUE} approximate the moment of SDE's at any time. 
}
\references{

   Guidoum AC, Boukhetala K (2020). 
   "Performing Parallel Monte Carlo and Moment Equations Methods for Itô and Stratonovich Stochastic Differential Systems: R Package Sim.DiffProc". 
   \emph{Journal of Statistical Software}, \bold{96}(2), 1--82. doi:10.18637/jss.v096.i02
   
  Rodriguez R, Tuckwell H (2000). 
  \emph{A dynamical system for the approximate moments of nonlinear stochastic models of spiking neurons and networks.} 
  Mathematical and Computer Modelling, 31(4), 175--180.
  
  Alibrandi U, Ricciardi G (2012). 
  \emph{Stochastic Methods in Nonlinear Structural Dynamics},
  3--60. Springer Vienna, Vienna. ISBN 978-3-7091-1306-6.
}
\author{
A.C. Guidoum, K. Boukhetala.
}

\seealso{
\code{\link{MCM.sde}} Monte-Carlo methods for SDE's.
}
\examples{
library(deSolve)
## Example 1: 1-dim
## dX(t) = mu * X(t) * dt + sigma * X(t) * dW(t)
## Symbolic ODE's of mean and variance
f <- expression(mu*x)
g <- expression(sigma*x)
res1 <- MEM.sde(drift=f,diffusion=g,type="ito")
res2 <- MEM.sde(drift=f,diffusion=g,type="str")
res1
res2
## numerical approximation of mean and variance
para <- c(mu=2,sigma=0.5)
t    <- seq(0,1,by=0.001)
init <- c(m=1,S=0)
res1 <- MEM.sde(drift=f,diffusion=g,solve=TRUE,init=init,parms=para,time=t)
res1
matplot.0D(res1$sol.ode,main="Mean and Variance of X(t), type Ito")
plot(res1$sol.ode,select=c("m","S"))
## approximation at time = 0.75
summary(res1,at=0.75)

##
res2 <- MEM.sde(drift=f,diffusion=g,solve=TRUE,init=init,parms=para,time=t,type="str")
res2
matplot.0D(res2$sol.ode,main="Mean and Variance of X(t), type Stratonovich")
plot(res2$sol.ode,select=c("m","S"))
## approximation at time = 0.75
summary(res2,at=0.75)

## Comparison:

plot(res1$sol.ode, res2$sol.ode,ylab = c("m(t)"),select="m",xlab = "Time",
     col = c("red", "blue"))
plot(res1$sol.ode, res2$sol.ode,ylab = c("S(t)"),select="S",xlab = "Time",
     col = c("red", "blue"))
	 
## Example2: 2-dim
## dX(t) = 1/mu*(theta-X(t)) dt + sqrt(sigma) * dW1(t),
## dY(t) = X(t) dt + 0 * dW2(t)	 
\dontrun{
para=c(mu=0.75,sigma=0.1,theta=2)
init=c(m1=0,m2=0,S1=0,S2=0,C12=0)
t <- seq(0,10,by=0.001)
f <- expression(1/mu*(theta-x), x)  
g <- expression(sqrt(sigma),0)
res2d <- MEM.sde(drift=f,diffusion=g,solve=TRUE,init=init,parms=para,time=t)
res2d

## Exact moment

mu=0.75;sigma=0.1;theta=2;x0=0;y0=0
E_x <- function(t) theta+(x0-theta)*exp(-t/mu)
V_x <- function(t) 0.5*sigma*mu *(1-exp(-2*(t/mu)))
E_y <- function(t) y0+theta*t+(x0-theta)*mu*(1-exp(-t/mu))
V_y <- function(t) sigma*mu^3*((t/mu)-2*(1-exp(-t/mu))+0.5*(1-exp(-2*(t/mu))))
cov_xy <- function(t) 0.5*sigma*mu^2 *(1-2*exp(-t/mu)+exp(-2*(t/mu)))

## 
summary(res2d,at=5)
E_x(5);E_y(5);V_x(5);V_y(5);cov_xy(5)

matplot.0D(res2d$sol.ode,select=c("m1"))
curve(E_x,add=TRUE,col="red") 

## plot

plot(res2d$sol.ode)
matplot.0D(res2d$sol.ode,select=c("S1","S2","C12"))   
plot(res2d$sol.ode[,"m1"], res2d$sol.ode[,"m2"], xlab = "m1(t)",
  ylab = "m2(t)", type = "l",lwd = 2)
hist(res2d$sol.ode,select=c("m1","m2"), col = c("darkblue", "red", "orange", "black"))

## Example3: 2-dim with correlation
## Heston model
## dX(t) = mu*X(t) dt + sqrt(Y(t))*X(t) * dW1(t),
## dY(t) = lambda*(theta-Y(t)) dt + sigma*sqrt(Y(t)) * dW2(t)	
## with E(dw1dw2)=rho


f <- expression( mu*x, lambda*(theta-y) )
g <- expression( sqrt(y)*x, sigma*sqrt(y) )
RHO  <- expression(rho)
res2d <- MEM.sde(drift=f,diffusion=g,corr=RHO)
res2d

## Numerical approximation
RHO <- expression(0.5)
para=c(mu=1,lambda=3,theta=0.5,sigma=0.1)
ini=c(m1=10,m2=2,S1=0,S2=0,C12=0)
res2d = MEM.sde(drift=f,diffusion=g,solve=TRUE,parms=para,init=ini,time=seq(0,1,by=0.01))
res2d

matplot.0D(res2d$sol.ode,select=c("m1","m2")) 
matplot.0D(res2d$sol.ode,select=c("S1","S2","C12")) 
	 
## Example4: 3-dim
## dX(t) = sigma*(Y(t)-X(t)) dt + 0.1 * dW1(t)
## dY(t) = (rho*X(t)-Y(t)-X(t)*Z(t)) dt + 0.1 * dW2(t)
## dZ(t) = (X(t)*Y(t)-bet*Z(t)) dt + 0.1 * dW3(t)
## with E(dw1dw2)=rho1, E(dw1dw3)=rho2 and E(dw2dw3)=rho3


f <- expression(sigma*(y-x),rho*x-y-x*z,x*y-bet*z)
g <- expression(0.1,0.1,0.1)
RHO <- expression(rho1,rho2,rho3)
## Symbolic moments equations
res3d = MEM.sde(drift=f,diffusion=g,corr=RHO)
res3d

## Numerical approximation
RHO <- expression(0.5,0.2,-0.7)
para=c(sigma=10,rho=28,bet=8/3)
ini=c(m1=1,m2=1,m3=1,S1=0,S2=0,S3=0,C12=0,C13=0,C23=0)
res3d = MEM.sde(drift=f,diffusion=g,solve=T,parms=para,init=ini,time=seq(0,1,by=0.01))
res3d

summary(res3d,at=0.25)
summary(res3d,at=0.50)
summary(res3d,at=0.75)

plot(res3d$sol.ode)
matplot.0D(res3d$sol.ode,select=c("m1","m2","m3")) 
matplot.0D(res3d$sol.ode,select=c("S1","S2","S3")) 
matplot.0D(res3d$sol.ode,select=c("C12","C13","C23")) 

plot3D(res3d$sol.ode[,"m1"], res3d$sol.ode[,"m2"],res3d$sol.ode[,"m3"], xlab = "m1(t)",
  ylab = "m2(t)",zlab="m3(t)", type = "l",lwd = 2,box=F)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{moment equations}
\keyword{sde}
