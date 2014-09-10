\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{bridgesde2d}
\alias{bridgesde2d}
\alias{bridgesde2d.default}
\alias{print.bridgesde2d}
\alias{time.bridgesde2d}
\alias{mean.bridgesde2d}
\alias{median.bridgesde2d}
\alias{quantile.bridgesde2d}
\alias{kurtosis.bridgesde2d}
\alias{skewness.bridgesde2d}
\alias{moment.bridgesde2d}
\alias{bconfint.bridgesde2d}
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
\method{bridgesde2d}{default}(N = 1000, M = 1, x0 = c(0, 0), y = c(1, 1), t0 = 0, T = 1, Dt, 
   driftx, diffx, drifty, diffy, alpha = 0.5, mu = 0.5, type = c("ito", "str"),
   method = c("euler", "milstein", "predcorr", "smilstein", "taylor", 
   "heun", "rk1", "rk2", "rk3"), \dots)
							  
							  
\method{time}{bridgesde2d}(x, \dots)
\method{mean}{bridgesde2d}(x, \dots)
\method{median}{bridgesde2d}(x, \dots)
\method{quantile}{bridgesde2d}(x, \dots)
\method{kurtosis}{bridgesde2d}(x, \dots)
\method{skewness}{bridgesde2d}(x, \dots)
\method{moment}{bridgesde2d}(x, order = 2, \dots)
\method{bconfint}{bridgesde2d}(x, level=0.95, \dots)
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
  \item{M}{number of trajectories.}
  \item{x0}{initial value (numeric vector of length 2) of the process \eqn{X_t}{X(t)} and \eqn{Y_t}{Y(t)} at time \eqn{t_0}{t0}.}
  \item{y}{terminal value (numeric vector of length 2) of the process \eqn{X_t}{X(t)} and \eqn{Y_t}{Y(t)} at time \eqn{T}{T}.}
  \item{t0}{initial time.}
  \item{T}{final time.} 
  \item{Dt}{time step of the simulation (discretization). If it is \code{\link[base]{missing}} a default \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}.}
  \item{driftx, drifty}{drift coefficient: an \code{\link[base]{expression}} of three variables \code{t}, \code{x} and \code{y} for process \eqn{X_t}{X(t)} and \eqn{Y_t}{Y(t)}.}
  \item{diffx, diffy}{diffusion coefficient: an \code{\link[base]{expression}} of three variables \code{t}, \code{x} and \code{y} for process \eqn{X_t}{X(t)} and \eqn{Y_t}{Y(t)}.}
  \item{alpha, mu}{weight of the predictor-corrector scheme; the default \code{alpha = 0.5} and \code{mu = 0.5}.}
  \item{type}{if \code{type="ito"} simulation diffusion bridge of Ito type, else \code{type="str"} simulation diffusion bridge of Stratonovich type; the default \code{type="ito"}.}
  \item{method}{numerical methods of simulation, the default \code{method = "euler"}; see \code{\link{snssde2d}}.}
  \item{x}{an object inheriting from class \code{"bridgesde2d"}.}
  \item{order}{order of moment.}
  \item{level}{the confidence level required.} 
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
\item{X, Y }{an invisible \code{mts} (2-dim) object (X(t),Y(t)).}
\item{driftx, drifty }{drift coefficient of X(t) and Y(t).}
\item{diffx, diffy }{diffusion coefficient of X(t) and Y(t).}
\item{Cx, Cy}{numbre of crossing realized of X(t) (Y(t)).}
\item{type }{type of sde.}
\item{method }{the numerical method used.}
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
\code{\link{bridgesde1d}} for simulation of 1-dim diffusion bridge. 

\code{\link[sde]{DBridge}} in package \CRANpkg{sde}.
}
\examples{
## dX(t) = 4*(-1-X(t))*Y(t) dt + 0.2 dW1(t)
## dY(t) = 4*(1-Y(t))*X(t) dt + 0.2 dW2(t)
## x01 = 0 , y01 = 1
## x02 = -1, y02 = 0 
## W1(t) and W2(t) two independent Brownian motion

fx <- expression(4*(-1-x)*y)
gx <- expression(0.2)
fy <- expression(4*(1-y)*x)
gy <- expression(0.2)

res <- bridgesde2d(x0=c(0,-1),y=c(1,0),driftx=fx,diffx=gx,drifty=fy,diffy=gy,M=50)
res
plot(res)
dev.new()
plot2d(res,type="n")
points2d(res,col=rgb(0,100,0,50,maxColorValue=255), pch=16)
}
\keyword{sde}
\keyword{ts}
\keyword{mts}