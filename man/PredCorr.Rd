\name{PredCorr}
\alias{PredCorr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predictor-Corrector Method
}
\description{
Predictor-Corrector method of simulation numerical solution of stochastic differential equation.
}
\usage{
PredCorr(N, M, T = 1, t0, x0, Dt, alpha = 0.5,
         mu = 0.5, drift, diffusion, output = FALSE)
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
time step of the simulation (\code{discretization}).
}
  \item{alpha}{
weight \code{alpha} of the predictor-corrector scheme.
}
  \item{mu}{
weight \code{mu} of the predictor-corrector scheme.
}
  \item{drift}{
drift coefficient: an expression of two variables \code{t} and \code{x}.
}
  \item{diffusion}{
diffusion coefficient: an expression of two variables \code{t} and \code{x}.
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel 2007.
}
}
\details{
The function returns a trajectory of the process; i.e., \code{x0} and the new \code{N} simulated values if \code{ M = 1}. For \code{M > 1}, an mts (multidimensional trajectories) is returned, which means that \code{M} independent trajectories are simulated. If \code{Dt} is not specified, then \code{Dt = (T-t0)/N}. If \code{Dt} is specified, then \code{N} values of the solution of the sde are generated and the time horizon \code{T} is adjusted to be \code{T = N * Dt}.
  
The method we present here just tries to approximate
the states of the process first. This method is of weak convergence order \code{1}. 

The predictor-corrector algorithm is as follows. First consider the simple approximation (the predictor), Then choose two weighting coefficients \code{alpha} and \code{mu} in \code{[0,1]} and calculate the corrector.
}
\value{
data.frame(time,x) and plot of process.
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\note{
\itemize{\item Note that the predictor-corrector method falls back to the standard Euler method for \code{alpha = mu = 0}.\item The function by default implements the predictor corrector
method with \code{alpha = mu = 0.5}.}
}
\seealso{
\code{\link{diffBridge}} Creating Diffusion Bridge Models.\code{\link{snssde}} Simulation numerical solution of stochastic differential equation.
}
\examples{
## example 1
## Hull-White/Vasicek Model
## T = 1 , t0 = 0 and N = 1000 ===> Dt = 0.001 
drift     <- expression( (3*(2-x)) )
diffusion <- expression( (2) )
PredCorr(N=1000, M=1, T = 1, t0=0, x0=10, Dt=0.001, alpha = 0.5,
         mu = 0.5, drift, diffusion, output = FALSE)
## Multiple trajectories of the OU process by Euler Scheme
PredCorr(N=1000,M=5,T=1,t0=0,x0=10,Dt=0.001,alpha = 0.5,
         mu = 0.5,drift,diffusion,output=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
\keyword{ Stochastic Differential Equation }
\keyword{Numerical Solution of Stochastic Differential Equation}
