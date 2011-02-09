\name{snssde}
\alias{snssde}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Numerical Solution of One-Dimensional SDE
}
\description{
Different methods of simulation of solutions to stochastic differential equations one-dimensional. 
}
\usage{
snssde(N, M, T = 1, t0, x0, Dt, drift, diffusion, Output = FALSE,                
       Methods = c("SchEuler", "SchMilstein", "SchMilsteinS",
       "SchTaylor", "SchHeun", "SchRK3"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of process.
}
  \item{M}{
number of trajectories.
}
  \item{T}{
final time.
}
  \item{t0}{
initial time.
}
  \item{x0}{
initial value of the process at time \code{t0}.
}
  \item{Dt}{
time step of the simulation \code{(discretization)}.
}
  \item{drift}{
drift coefficient: an expression of two variables \code{t} and \code{x}.
}
  \item{diffusion}{
diffusion coefficient: an expression of two variables \code{t} and \code{x}.
}
  \item{Output}{
if \code{Output = TRUE} write a \code{Output} to an Excel 2007.
}
  \item{Methods}{
method of simulation ,\code{see details}.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
The function \code{snssde} returns a trajectory of the process; i.e., \code{x0} and the new \code{N} simulated values if \code{M = 1}. For \code{M > 1}, an mts (multidimensional trajectories) is returned, which means that \code{M} independent trajectories are simulated. \code{Dt} the best discretization \code{Dt = (T-t0)/N}. 

Simulation methods are usually based on discrete approximations of the continuous solution to a stochastic differential equation. The methods of approximation
are classified according to their different properties. Mainly two criteria of optimality are used in the literature: the strong and the weak (orders of) convergence. The \code{methods} of simulation can be one among: \code{Euler Order 0.5} , \code{Milstein Order 1} , \code{Milstein Second-Order} , \code{Ito-Taylor Order 1.5} , \code{Heun Order 2} , \code{Runge-Kutta Order 3}. 
}
\value{
data.frame(time,x) and plot of process.
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\note{
\itemize{\item If \code{methods} is not specified, it is assumed to be the \code{Euler Scheme}. \item If \code{T} and \code{t0} specified, the best discretization \code{Dt = (T-t0)/N}.}
}
\seealso{
\code{\link{diffBridge}} Creating Diffusion Bridge Models.\code{\link{PredCorr}} Predictor-Corrector Method for one-dimensional SDE. \code{\link{snssde2D}} numerical solution of two-dimensional SDE. \code{\link{PredCorr2D}} predictor-corrector method for two-dimensional SDE.

}
\examples{
## example 1
## Hull-White/Vasicek Model
## T = 1 , t0 = 0 and N = 1000 ===> Dt = 0.001 
drift     <- expression( (3*(2-x)) )
diffusion <- expression( (2) )
snssde(N=1000,M=1,T=1,t0=0,x0=10,Dt=0.001,
drift,diffusion,Output=FALSE)
## Multiple trajectories of the OU process by Euler Scheme
snssde(N=1000,M=5,T=1,t0=0,x0=10,Dt=0.001,
drift,diffusion,Output=FALSE)

## example 2
## Black-Scholes models
## T = 1 , t0 = 0 and N = 1000 ===> Dt = 0.001 
drift     <- expression( (3*x) )
diffusion <- expression( (2*x) )
snssde(N=1000,M=1,T=1,t0=0,x0=10,Dt=0.001,drift,
diffusion,Output=FALSE,Methods="SchMilstein")

## example 3
## Constant Elasticity of Variance (CEV) Models
## T = 1 , t0 = 0 and N = 1000 ===> Dt = 0.001 
drift     <- expression( (0.3*x) )
diffusion <- expression( (0.2*x^0.75) )
snssde(N=1000,M=1,T=1,t0=0,x0=1,Dt=0.001,drift,
diffusion,Output=FALSE,Methods="SchMilsteinS")

## example 4
## sde\\ dX(t)=(0.03*t*X(t)-X(t)^3)*dt+0.1*dW(t)
## T = 100 , t0 = 0 and N = 1000 ===> Dt = 0.1 
drift     <- expression( (0.03*t*x-x^3) )
diffusion <- expression( (0.1) )
snssde(N=1000,M=1,T=100,t0=0,x0=0,Dt=0.1,drift,
diffusion,Output=FALSE,Methods="SchTaylor")

## example 5
## sde\\ dX(t)=cos(t*x)*dt+sin(t*x)*dW(t) by Heun Scheme
drift     <- expression( (cos(t*x)) )
diffusion <- expression( (sin(t*x)) )
snssde(N=1000,M=1,T=100,t0=0,x0=0,Dt=0.1,drift,
diffusion,Output=FALSE,Methods="SchHeun")

## example 6
## sde\\ dX(t)=exp(t)*dt+tan(t)*dW(t) by Runge-Kutta Scheme
drift     <- expression( (exp(t)) )
diffusion <- expression( (tan(t)) )
snssde(N=1000,M=1,T=1,t0=0,x0=1,Dt=0.001,drift,
diffusion,Output=FALSE,Methods="SchRK3")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
\keyword{ Stochastic Differential Equation }
\keyword{Numerical Solution of Stochastic Differential Equation}
