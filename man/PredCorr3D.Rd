\name{PredCorr3D}
\alias{PredCorr3D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predictor-Corrector Method For Three-Dimensional SDE
}
\description{
Predictor-Corrector method of simulation numerical solution of Three dimensional stochastic differential equations.
}
\usage{
PredCorr3D(N, T = 1, t0, x0, y0, z0, Dt, alpha = 0.5, mu = 0.5, driftx, drifty, driftz, diffx, diffy, diffz, Step = FALSE, Output = FALSE)
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
  \item{z0}{
initial value of the process \code{Z(t)} at time \code{t0}.
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
drift coefficient of process \code{X(t)}: an expression of three variables \code{t} , \code{x} and \code{y}, \code{z}.
}
  \item{drifty}{
drift coefficient of process \code{Y(t)}: an expression of three variables \code{t} , \code{x} and \code{y}, \code{z}.
}
  \item{driftz}{
drift coefficient of process \code{Z(t)}: an expression of three variables \code{t} , \code{x} and \code{y}, \code{z}.
}
  \item{diffx}{
diffusion coefficient of process \code{X(t)}: an expression of three variables \code{t} , \code{x} and \code{y}, \code{z}.
}
  \item{diffy}{
diffusion coefficient of process \code{Y(t)}: an expression of three variables \code{t} , \code{x} and \code{y}, \code{z}.
}
  \item{diffz}{
diffusion coefficient of process \code{Z(t)}: an expression of three variables \code{t} , \code{x} and \code{y}, \code{z}.
}
  \item{Step}{
if \code{Step = TRUE} ploting step by step.
}
  \item{Output}{
if \code{output = TRUE} write a \code{output} to an Excel 2007.
}
}
\details{
the system for stochastic differential equation Three dimensional is :\emph{\deqn{dX(t) = ax(t,X(t),Y(t),Z(t))* dt + bx(t,X(t),Y(t),Z(t))* dWx(t)}}
\emph{\deqn{dY(t) = ay(t,X(t),Y(t),Z(t))* dt + by(t,X(t),Y(t),Z(t))* dWy(t)}}
\emph{\deqn{dZ(t) = az(t,X(t),Y(t),Z(t))* dt + bz(t,X(t),Y(t),Z(t))* dWz(t)}}
with \code{driftx=ax(t,X(t),Y(t),Z(t))}, \code{drifty=ay(t,X(t),Y(t),Z(t))}, \code{driftz=az(t,X(t),Y(t),Z(t))} and 
\code{diffx=bx(t,X(t),Y(t),Z(t))}, \code{diffy=by(t,X(t),Y(t),Z(t))},\code{diffz=bz(t,X(t),Y(t),Z(t))}

The method we present here just tries to approximate
the states of the process first. This method is of weak convergence order \code{1}. \code{dW1(t)} and \code{dW2(t)}, \code{dW3(t)} are brownian motions independent. 

The predictor-corrector algorithm is as follows. First consider the simple approximation (the predictor), Then choose two weighting coefficients \code{alpha} and \code{mu} in \code{[0,1]} and calculate the corrector.
}
\value{
data.frame(time,X(t),Y(t),Z(t)) and plot of process 3-D.
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\note{
\itemize{\item Note that the predictor-corrector method falls back to the standard Euler method for \code{alpha = mu = 0}.\item The function by default implements the predictor corrector
method with \code{alpha = mu = 0.5}.}
}
\seealso{
\code{\link{snssde}} numerical solution of one-dimensional SDE. \code{\link{snssde2D}} numerical solution of Two-dimensional SDE. \code{\link{snssde3D}} numerical solution of Three-dimensional SDE. \code{\link{PredCorr}} predictor-corrector method for one-dimensional SDE. \code{\link{PredCorr2D}} predictor-corrector method for Two-dimensional SDE.
}
\examples{
driftx <- expression(0)
drifty <- expression(0)
driftz <- expression(0)
diffx <- expression(1)
diffy <- expression(1)
diffz <- expression(1)
PredCorr3D(N=1000, T = 1, t0=0, x0=0, y0=0, z0=0, Dt=0.001, alpha = 0.5, mu = 0.5, 
    driftx, drifty, driftz, diffx, diffy, diffz)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process Multidimensional}
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Stochastic Differential Equation Multidimensional}
\keyword{ Numerical Solution of Stochastic Differential Equation Multidimensional}
