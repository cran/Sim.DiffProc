\name{snssde3D}
\alias{snssde3D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Numerical Solution of Three-Dimensional SDE
}
\description{
Different methods of simulation of solutions to stochastic differential equations Three-dimensional.
}
\usage{
snssde3D(N, T = 1, t0, x0, y0, z0, Dt, driftx, drifty, driftz, diffx, 
       diffy, diffz, Step = FALSE, Output = FALSE, Methods = c        
       ("SchEuler", "SchMilstein", "SchMilsteinS", 
        "SchTaylor", "SchHeun", "SchRK3"), ...)
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
  \item{driftx}{
drift coefficient of process \code{X(t)}: an expression of variables \code{t} , \code{x} and \code{y}, \code{z}.
}
  \item{drifty}{
drift coefficient of process \code{Y(t)}: an expression of variables \code{t} , \code{x} and \code{y}, \code{z}.
}
  \item{driftz}{
drift coefficient of process \code{Z(t)}: an expression of variables \code{t} , \code{x} and \code{y}, \code{z}.
}
  \item{diffx}{
diffusion coefficient of process \code{X(t)}: an expression of variables \code{t} , \code{x} and \code{y}, \code{z}.
}
  \item{diffy}{
diffusion coefficient of process \code{Y(t)}: an expression of variables \code{t} , \code{x} and \code{y}, \code{z}.
}
  \item{diffz}{
diffusion coefficient of process \code{Z(t)}: an expression of variables \code{t} , \code{x} and \code{y}, \code{z}.
}
  \item{Step}{
if \code{Step = TRUE} ploting step by step.
}
  \item{Output}{
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
  \item{Methods}{
method of simulation ,\code{see details}.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
the system for stochastic differential equation Two dimensional is :\emph{\deqn{dX(t) = ax(t,X(t),Y(t),Z(t))* dt + bx(t,X(t),Y(t),Z(t))* dW1(t)}}
\emph{\deqn{dY(t) = ay(t,X(t),Y(t),Z(t))* dt + by(t,X(t),Y(t),Z(t))* dW2(t)}}
\emph{\deqn{dZ(t) = az(t,X(t),Y(t),Z(t))* dt + bz(t,X(t),Y(t),Z(t))* dW3(t)}}
 

Simulation methods are usually based on discrete approximations of the continuous solution to a stochastic differential equation. The methods of approximation
are classified according to their different properties. Mainly two criteria of optimality are used in the literature: the strong and the weak (orders of) convergence. The \code{methods} of simulation can be one among: \code{Euler Order 0.5} , \code{Milstein Order 1} , \code{Milstein Second-Order} , \code{Ito-Taylor Order 1.5} , \code{Heun Order 2} , \code{Runge-Kutta Order 3}.
}
\value{
data.frame(time,X(t),Y(t),Z(t)) and plot of process 3-D.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\note{
\itemize{\item If \code{methods} is not specified, it is assumed to be the \code{Euler Scheme}. \item If \code{T} and \code{t0} specified, the best discretization \code{Dt = (T-t0)/N}.}
}
\seealso{
\code{\link{snssde}} numerical solution of one-dimensional SDE. \code{\link{snssde2D}} numerical solution of Two-dimensional SDE. \code{\link{PredCorr}} predictor-corrector method for one-dimensional SDE.  \code{\link{PredCorr2D}} predictor-corrector method for Two-dimensional SDE. \code{\link{PredCorr3D}} predictor-corrector method for Three-dimensional SDE.
}
\examples{
\donttest{
 driftx <- expression(0)
 drifty <- expression(0)
 driftz <- expression(0)
 diffx <- expression(1)
 diffy <- expression(1)
 diffz <- expression(1)
 snssde3D(N=1000, T = 1, t0=0, x0=0, y0=0, z0=0, Dt=0.001,  
          driftx, drifty, driftz, diffx,diffy, diffz)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process Multidimensional}
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Stochastic Differential Equation Multidimensional}
\keyword{ Numerical Solution of Stochastic Differential Equation Multidimensional}
