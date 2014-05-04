\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{bridgesde2d}
\alias{bridgesde2d}
\alias{bridgesde2d.default}
\alias{print.bridgesde2d}
\alias{time.bridgesde2d}
\alias{summary.bridgesde2d}
\alias{plot.bridgesde2d}
\alias{points.bridgesde2d}
\alias{lines.bridgesde2d}
\alias{plot2d.bridgesde2d}
\alias{points2d.bridgesde2d}
\alias{lines2d.bridgesde2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of 2-Dim Diffusion Bridge}
\description{The (S3) generic function \code{bridgesde2d} for simulation of 2-dim diffusion bridge.}
\usage{
bridgesde2d(N, \dots)
\method{bridgesde2d}{default}(N = 1000, x0 = c(0, 0), y = c(1, 1), t0 = 0, T = 1, Dt, 
   driftx, diffx, drifty, diffy, alpha = 0.5, mu = 0.5, type = c("ito", "str"),
   method = c("euler", "milstein", "predcorr", "smilstein", "taylor", 
   "heun", "rk1", "rk2", "rk3"), \dots)
							  

\method{summary}{bridgesde2d}(object, \dots)							  
\method{time}{bridgesde2d}(x, \dots)
\method{plot}{bridgesde2d}(x, \dots)
\method{lines}{bridgesde2d}(x, \dots)
\method{points}{bridgesde2d}(x, \dots)	
\method{plot2d}{bridgesde2d}(x, \dots)
\method{lines2d}{bridgesde2d}(x, \dots)
\method{points2d}{bridgesde2d}(x, \dots)								  
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{N}{number of simulation steps.}
  \item{x0}{initial value (numeric vector of length 2) of the process \eqn{X_t}{X(t)} and \eqn{Y_t}{Y(t)} at time \eqn{t_0}{t0}.}
  \item{y}{terminal value (numeric vector of length 2) of the process \eqn{X_t}{X(t)} and \eqn{Y_t}{Y(t)} at time \eqn{T}{T}.}
  \item{t0}{initial time.}
  \item{T}{final time.} 
  \item{Dt}{time step of the simulation (discretization). If it is \code{\link[base]{missing}} a default \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}.}
  \item{driftx}{drift coefficient: an \code{\link[base]{expression}} of three variables \code{t}, \code{x} and \code{y} for process \eqn{X_t}{X(t)}.}
  \item{diffx}{diffusion coefficient: an \code{\link[base]{expression}} of three variables \code{t}, \code{x} and \code{y} for process \eqn{X_t}{X(t)}.}
  \item{drifty}{drift coefficient: an \code{\link[base]{expression}} of three variables \code{t}, \code{x} and \code{y} for process \eqn{Y_t}{Y(t)}.}
  \item{diffy}{diffusion coefficient: an \code{\link[base]{expression}} of three variables \code{t}, \code{x} and \code{y} for process \eqn{Y_t}{Y(t)}.}
  \item{alpha}{weight \code{alpha} of the predictor-corrector scheme; the default \code{alpha = 0.5}.}
  \item{mu}{weight \code{mu} of the predictor-corrector scheme; the default \code{mu = 0.5}.}
  \item{type}{if \code{type="ito"} simulation diffusion bridge of Ito type, else \code{type="str"} simulation diffusion bridge of Stratonovich type; the default \code{type="ito"}.}
  \item{method}{numerical methods of simulation, the default \code{method = "euler"}; see \code{\link{snssde2d}}.}
  \item{x, object}{an object inheriting from class \code{"bridgesde2d"}.}
  \item{\dots}{further arguments for (non-default) methods.}
}
\details{
The function \code{bridgesde2d} returns a \code{mts} of the diffusion bridge starting at \code{x} at time \code{t0} and ending at \code{y} at time \code{T}.

The methods of approximation are classified according to their different properties. Mainly two criteria of optimality are used in the literature: the strong 
and the weak (orders of) convergence. The \code{method} of simulation can be one among: \code{Euler-Maruyama Order 0.5}, \code{Milstein Order 1}, \code{Milstein Second-Order}, 
\code{Predictor-Corrector method}, \code{Ito-Taylor Order 1.5}, \code{Heun Order 2} and \code{Runge-Kutta Order 1, 2 and 3}.

For more details see \code{vignette("SDEs")}.
}
\value{
\code{bridgesde2d} returns an object inheriting from \code{\link[base]{class}} \code{"bridgesde2d"}.
\item{XY }{an invisible \code{\link[stats]{ts}} object \eqn{(X_{t},Y_{t})}{(X(t),Y(t))}.}
}
\references{
   Bladt, M. and Sorensen, M. (2007).
   Simple simulation of diffusion bridges with application to likelihood inference for diffusions.
   \emph{Working Paper, University of Copenhagen}. 
   Available at \url{http://www.math.ku.dk/~michael/diffusionbridgepreprint.pdf}

   Iacus, S.M. (2008). 
   \emph{Simulation and inference for stochastic differential equations: with R examples}. 
   Springer-Verlag, New York
}
\author{
A.C. Guidoum, K. Boukhetala.
}

\seealso{
\code{\link{bridgesde1d}} for simulation of 1-dim diffusion bridge. \code{\link[sde]{DBridge}} in package \CRANpkg{sde}.
}
\examples{
## dX(t) = 4*(-2-X(t)) * dt + Y(t) * dW1(t)          
## dY(t) = 4*(-2-Y(t)) * dt + X(t) * dW2(t)
## x01 = 1 , y01 = 1
## x02 = 1, y02 = 1 
## W1(t) and W2(t) two independent Brownian motion

fx <- expression(4*(-2-x))
gx <- expression(y)
fy <- expression(4*(-2-y))
gy <- expression(x)

XY <- bridgesde2d(x0=c(1,1),y=c(1,1),driftx=fx,diffx=gx,drifty=fy,diffy=gy)
XY
summary(XY)
plot(XY)
dev.new()
plot(XY,plot.type="single")
dev.new()
plot2d(XY,type="n")
points2d(XY,col=rgb(0,100,0,50,maxColorValue=255), pch=16)
}
\keyword{sde}
\keyword{ts}
\keyword{mts}