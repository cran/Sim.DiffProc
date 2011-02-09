\name{snssde2D}
\alias{snssde2D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Numerical Solution of Two-Dimensional SDE
}
\description{
Different methods of simulation of solutions to stochastic differential equations Two-dimensional. 
}
\usage{
snssde2D(N, T = 1, t0, x0, y0, Dt, driftx, drifty, diffx, diffy, 
        Step = FALSE, Output = FALSE, Methods = c("SchEuler", 
        "SchMilstein", "SchMilsteinS", "SchTaylor", "SchHeun", 
        "SchRK3"), ...)
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
  \item{Methods}{
method of simulation ,\code{see details}.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
the system for stochastic differential equation Two dimensional is :\emph{\deqn{dX(t) = ax(t,X(t),Y(t))* dt + bx(t,X(t),Y(t))* dW1(t)}}
\emph{\deqn{dY(t) = ay(t,X(t),Y(t))* dt + by(t,X(t),Y(t))* dW2(t)}}
with \code{driftx=ax(t,X(t),Y(t))}, \code{drifty=ay(t,X(t),Y(t))} and 
\code{diffx=bx(t,X(t),Y(t))}, \code{diffy=by(t,X(t),Y(t))}.
\code{dW1(t)} and \code{dW2(t)} are brownian motions independent. 

Simulation methods are usually based on discrete approximations of the continuous solution to a stochastic differential equation. The methods of approximation
are classified according to their different properties. Mainly two criteria of optimality are used in the literature: the strong and the weak (orders of) convergence. The \code{methods} of simulation can be one among: \code{Euler Order 0.5} , \code{Milstein Order 1} , \code{Milstein Second-Order} , \code{Ito-Taylor Order 1.5} , \code{Heun Order 2} , \code{Runge-Kutta Order 3}.
}
\value{
data.frame(time,X(t),Y(t)) and plot of process 2-D.
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\note{
\itemize{\item If \code{methods} is not specified, it is assumed to be the \code{Euler Scheme}. \item If \code{T} and \code{t0} specified, the best discretization \code{Dt = (T-t0)/N}.}
}
\seealso{
\code{\link{diffBridge}} Creating Diffusion Bridge Models. \code{\link{snssde}} numerical solution of one-dimensional SDE. \code{\link{snssde}} numerical solution of one-dimensional SDE. \code{\link{PredCorr}} predictor-corrector method for one-dimensional SDE.  \code{\link{PredCorr2D}} predictor-corrector method for Two-dimensional SDE.
}
\examples{

## Example 1 
driftx <- expression(cos(t*x))
drifty <- expression(cos(t*y))
diffx <- expression(sin(t*x))
diffy <- expression(sin(t*y))
snssde2D(N=1000, T = 1, t0=0, x0=0, y0=0, Dt=0.001, driftx, 
         drifty, diffx, diffy, Step = FALSE, Output = FALSE,
         Methods="SchTaylor")

## Example 2
driftx <- expression(cos(t*x*y))
drifty <- expression(sin(t*y*y))
diffx <- expression(atan2(y, x))
diffy <- expression(atan2(y, x))
snssde2D(N=5000, T = 1, t0=0, x0=1, y0=1, Dt=0.001, driftx, 
         drifty, diffx, diffy, Step = FALSE, Output = FALSE,
         Methods="SchHeun")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process Multidimensional}
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Stochastic Differential Equation Multidimensional}
\keyword{ Numerical Solution of Stochastic Differential Equation Multidimensional}
