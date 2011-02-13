\name{PredCorr2D}
\alias{PredCorr2D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predictor-Corrector Method For Two-Dimensional SDE
}
\description{
Predictor-Corrector method of simulation numerical solution of Two dimensional stochastic differential equation.
}
\usage{
PredCorr2D(N, T = 1, t0, x0, y0, Dt, alpha = 0.5, mu = 0.5, driftx, 
           drifty, diffx, diffy, Step = FALSE, Output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of process.
}
  \item{T}{
final time.
}
  \item{t0}{
initial time.
}
  \item{x0}{
initial value of the process \code{X(t)} at time \code{t0}.
}
  \item{y0}{
initial value of the process \code{Y(t)} at time \code{t0}.
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
  \item{driftx}{
drift coefficient of process \code{X(t)}: an expression of three variables \code{t} , \code{x} and \code{y}.
}
  \item{drifty}{
drift coefficient of process \code{Y(t)}: an expression of three variables \code{t} , \code{x} and \code{y}.
}
  \item{diffx}{
diffusion coefficient of process \code{X(t)}: an expression of three variables \code{t} , \code{x} and \code{y}.
}
  \item{diffy}{
diffusion coefficient of process \code{Y(t)}: an expression of three variables \code{t} , \code{x} and \code{y}.
}
  \item{Step}{
if \code{Step = TRUE} ploting step by step.
}
  \item{Output}{
if \code{output = TRUE} write a \code{output} to an Excel 2007.
}
}
\details{
the system for stochastic differential equation Two dimensional is :\emph{\deqn{dX(t) = ax(t,X(t),Y(t))* dt + bx(t,X(t),Y(t))* dWx(t)}}
\emph{\deqn{dY(t) = ay(t,X(t),Y(t))* dt + by(t,X(t),Y(t))* dWy(t)}}
with \code{driftx=ax(t,X(t),Y(t))}, \code{drifty=ay(t,X(t),Y(t))} and 
\code{diffx=bx(t,X(t),Y(t))}, \code{diffy=by(t,X(t),Y(t))}.

The method we present here just tries to approximate
the states of the process first. This method is of weak convergence order \code{1}. \code{dW1(t)} and \code{dW2(t)} are brownian motions independent. 

The predictor-corrector algorithm is as follows. First consider the simple approximation (the predictor), Then choose two weighting coefficients \code{alpha} and \code{mu} in \code{[0,1]} and calculate the corrector.
}
\value{
data.frame(time,X(t),Y(t)) and plot of process 2-D.
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\note{
\itemize{\item Note that the predictor-corrector method falls back to the standard Euler method for \code{alpha = mu = 0}.\item The function by default implements the predictor corrector
method with \code{alpha = mu = 0.5}.
}
}
\seealso{
\code{\link{diffBridge}} Creating Diffusion Bridge Models. \code{\link{snssde}} numerical solution of one-dimensional SDE. \code{\link{snssde2D}} numerical solution of Two-dimensional SDE. \code{\link{PredCorr}} predictor-corrector method for one-dimensional SDE.
}
\examples{
## Example 1
driftx <- expression(cos(t*x*y))
drifty <- expression(cos(t))
diffx <- expression(0.1)
diffy <- expression(0.1)
PredCorr2D(N=5000, T = 1, t0=0, x0=0, y0=0, Dt=0.001, alpha = 0.5, 
           mu = 0.5, driftx, drifty, diffx, diffy, Step = FALSE, 
           Output = FALSE)
## ploting Step by Step
##PredCorr2D(N=5000, T = 1, t0=0, x0=0, y0=0, Dt=0.001, alpha = 0.5, 
##           mu = 0.5, driftx, drifty, diffx, diffy, Step = TRUE, 
##           Output = FALSE)

## Example 2
## BM 2-D
driftx <- expression(0)
drifty <- expression(0)
diffx <- expression(1)
diffy <- expression(1)
PredCorr2D(N=5000, T = 1, t0=0, x0=0, y0=0, Dt=0.001, alpha = 0.5, 
           mu = 0.5, driftx, drifty, diffx, diffy, Step = FALSE, 
           Output = FALSE)
## ploting Step by Step
##PredCorr2D(N=5000, T = 1, t0=0, x0=0, y0=0, Dt=0.001, alpha = 0.5, 
##           mu = 0.5, driftx, drifty, diffx, diffy, Step = TRUE, 
##           Output = FALSE)


## Example 3
driftx <- expression(0.03*t*x-x^3)
drifty <- expression(0.03*t*y-y^3)
diffx <- expression(0.1)
diffy <- expression(0.1)
PredCorr2D(N=5000, T = 1, t0=0, x0=0, y0=0, Dt=0.001, alpha = 0.5, 
           mu = 0.5, driftx, drifty, diffx, diffy, Step = FALSE, 
           Output = FALSE)
## ploting Step by Step
## PredCorr2D(N=5000, T = 1, t0=0, x0=0, y0=0, Dt=0.001, alpha = 0.5, 
##           mu = 0.5, driftx, drifty, diffx, diffy, Step = FALSE, 
##           Output = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process Multidimensional}
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Stochastic Differential Equation Multidimensional}
\keyword{ Numerical Solution of Stochastic Differential Equation Multidimensional}
