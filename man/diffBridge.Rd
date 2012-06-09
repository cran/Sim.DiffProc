\name{diffBridge}
\alias{diffBridge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Diffusion Bridge Models (by Euler Scheme)
}
\description{
Simulation of diffusion bridge models by euler scheme.
}
\usage{
diffBridge(N, t0, T, x, y, drift, diffusion, Output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of process.
}
  \item{t0}{
initial time.
}
  \item{T}{
final time.
}
  \item{x}{
initial value of the process at time \code{t0}.
}
  \item{y}{
terminal value of the process at time \code{T}.
}
  \item{drift}{
drift coefficient: an expression of two variables \code{t} and \code{x}.
}
  \item{diffusion}{
diffusion coefficient: an expression of two variables \code{t} and \code{x}.
}
  \item{Output}{
if \code{Output = TRUE} write a \code{Output} to an Excel (.csv).
}
}
\details{
The function \code{diffBridge} returns a trajectory of the diffusion bridge starting at \code{x} at time \code{t0} and ending at \code{y} at time \code{T}, the discretization \code{dt = (T-t0)/N}. 
}
\value{
data.frame(time,x) and plot of process.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{CEV}} Constant Elasticity of Variance Models, \code{\link{CKLS}} Chan-Karolyi-Longstaff-Sanders Models, \code{\link{CIR}} Cox-Ingersoll-Ross Models, \code{\link{CIRhy}} modified CIR and hyperbolic Process, \code{\link{DWP}} Double-Well Potential Model, \code{\link{GBM}} Model of Black-Scholes, \code{\link{HWV}} Hull-White/Vasicek Models, \code{\link{INFSR}} Inverse of Feller s Square Root models, \code{\link{JDP}} Jacobi Diffusion Process, \code{\link{PDP}} Pearson Diffusions Process, \code{\link{ROU}} Radial Ornstein-Uhlenbeck Process, \code{\link{snssde}} Simulation Numerical Solution of SDE.
}
\examples{
\donttest{
## example 1 : Ornstein-Uhlenbeck Bridge Model (x0=1,t0=0,y=3,T=1)
 drift     <- expression( (3*(2-x)) )
 diffusion <- expression( (2) )
 diffBridge(N=1000,t0=0,T=1,x=1,y=1,drift,diffusion) 


## example 2 : Brownian Bridge Model (x0=0,t0=0,y=1,T=1)
 drift     <- expression( 0)
 diffusion <- expression( 1 )
 diffBridge(N=1000,t0=0,T=1,x=0,y=0,drift,diffusion) 

## example 3 : Geometric Brownian Bridge Model (x0=1,t0=1,y=3,T=3)
 drift     <- expression( (3*x) )
 diffusion <- expression( (2*x) )
 diffBridge(N=1000,t0=0,T=10,x=1,y=1,drift,diffusion)

## example 4 : sde\\ dX(t)=(0.03*t*X(t)-X(t)^3)*dt+0.1*dW(t) (x0=0,t0=0,y=2,T=100)
 drift     <- expression( (0.03*t*x-x^3) )
 diffusion <- expression( (0.1) )
 diffBridge(N=1000,t0=0,T=100,x=1,y=1,drift,diffusion)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
\keyword{ Stochastic Differential Equation }
\keyword{Numerical Solution of Stochastic Differential Equation}
