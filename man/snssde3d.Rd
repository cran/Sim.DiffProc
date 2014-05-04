\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{snssde3d}
\alias{snssde3d}
\alias{snssde3d.default}
\alias{summary.snssde3d}
\alias{print.snssde3d}
\alias{time.snssde3d}
\alias{plot.snssde3d}
\alias{points.snssde3d}
\alias{lines.snssde3d}
\alias{plot3D.snssde3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of 3-Dim Stochastic Differential Equation}
\description{The (S3) generic function \code{snssde3d} of simulation of solutions to 3-dim stochastic differential equations of Ito or Stratonovich type, with different methods.
}
\usage{
snssde3d(N, \dots)
\method{snssde3d}{default}(N = 100, x0 = 0, y0 = 0, z0 = 0, t0 = 0, T = 1, Dt, 
   driftx, diffx, drifty, diffy, driftz, diffz, alpha = 0.5, mu = 0.5, 
   type = c("ito", "str"), method = c("euler", "milstein","predcorr", 
   "smilstein", "taylor", "heun", "rk1", "rk2", "rk3"), \dots)
						   
						   
\method{summary}{snssde3d}(object, \dots)
\method{time}{snssde3d}(x, \dots)
\method{plot}{snssde3d}(x, \dots)
\method{lines}{snssde3d}(x, \dots)
\method{points}{snssde3d}(x, \dots)
\method{plot3D}{snssde3d}(x, display = c("persp","rgl"), \dots)
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{N}{number of simulation steps.}
  \item{x0}{initial value of the process \eqn{X_{t}}{X(t)} at time \code{t0}.}
  \item{y0}{initial value of the process \eqn{Y_{t}}{Y(t)} at time \code{t0}.}
  \item{z0}{initial value of the process \eqn{Z_{t}}{Z(t)} at time \code{t0}.}
  \item{t0}{initial time.}
  \item{T}{final time.} 
  \item{Dt}{time step of the simulation (discretization). If it is \code{\link[base]{missing}} a default \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}.}
  \item{driftx}{drift coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{X_t}{X(t)}.}
  \item{diffx}{diffusion coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{X_t}{X(t)}.}
  \item{drifty}{drift coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{Y_t}{Y(t)}.}
  \item{diffy}{diffusion coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{Y_t}{Y(t)}.}
  \item{driftz}{drift coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{Z_t}{Z(t)}.}
  \item{diffz}{diffusion coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{Z_t}{Z(t)}.}
  \item{alpha}{weight \code{alpha} of the predictor-corrector scheme; the default \code{alpha = 0.5}.}
  \item{mu}{weight \code{mu} of the predictor-corrector scheme; the default \code{mu = 0.5}.}
  \item{type}{if \code{type="ito"} simulation sde of Ito type, else \code{type="str"} simulation sde of Stratonovich type; the default \code{type="ito"}.}
  \item{method}{numerical methods of simulation, the default \code{method = "euler"}.}
  \item{x, object}{an object inheriting from class \code{"snssde3d"}.}
  \item{display}{\code{"persp"} perspective or \code{"rgl"} plots.} 
  \item{\dots}{further arguments for (non-default) methods.}
}
\details{
The function \code{snssde3d} returns a \code{mts} x of length \code{N+1}; i.e. solution of the 3-dim sde \eqn{(X_{t},Y_{t},Z_{t})}{(X(t),Y(t),Z(t))}
of Ito or Stratonovich types; If \code{Dt} is not specified, then the best discretization \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}. 

The 3-dim Ito stochastic differential equation is:
\deqn{dX(t) = a(t,X(t),Y(t),Z(t)) dt + b(t,X(t),Y(t),Z(t)) dW_{1}(t)}{dX(t) = a(t,X(t),Y(t),Z(t))*dt + b(t,X(t),Y(t),Z(t))*dW1(t)}
\deqn{dY(t) = a(t,X(t),Y(t),Z(t)) dt + b(t,X(t),Y(t),Z(t)) dW_{2}(t)}{dY(t) = a(t,X(t),Y(t),Z(t))*dt + b(t,X(t),Y(t),Z(t))*dW2(t)}
\deqn{dZ(t) = a(t,X(t),Y(t),Z(t)) dt + b(t,X(t),Y(t),Z(t)) dW_{3}(t)}{dZ(t) = a(t,X(t),Y(t),Z(t))*dt + b(t,X(t),Y(t),Z(t))*dW3(t)}
3-dim Stratonovich sde :
\deqn{dX(t) = a(t,X(t),Y(t),Z(t)) dt + b(t,X(t),Y(t),Z(t)) \circ dW_{1}(t)}{dX(t) = a(t,X(t),Y(t),Z(t))*dt + b(t,X(t),Y(t),Z(t)) o dW1(t)}
\deqn{dY(t) = a(t,X(t),Y(t),Z(t)) dt + b(t,X(t),Y(t),Z(t)) \circ dW_{2}(t)}{dY(t) = a(t,X(t),Y(t),Z(t))*dt + b(t,X(t),Y(t),Z(t)) o dW2(t)}
\deqn{dZ(t) = a(t,X(t),Y(t),Z(t)) dt + b(t,X(t),Y(t),Z(t)) \circ dW_{3}(t)}{dZ(t) = a(t,X(t),Y(t),Z(t))*dt + b(t,X(t),Y(t),Z(t)) o dW3(t)}
\eqn{W_{1}(t), W_{2}(t), W_{3}(t)}{W1(t), W2(t) and W3(t)} three standard Brownian motion independent.

The methods of approximation are classified according to their different properties. Mainly two criteria of optimality are used in the literature: the strong 
and the weak (orders of) convergence. The \code{method} of simulation can be one among: \code{Euler-Maruyama Order 0.5}, \code{Milstein Order 1}, \code{Milstein Second-Order}, 
\code{Predictor-Corrector method}, \code{Ito-Taylor Order 1.5}, \code{Heun Order 2} and \code{Runge-Kutta Order 1, 2 and 3}. 

For more details see \code{vignette("SDEs")}.
}
\value{
\code{snssde3d} returns an object inheriting from \code{\link[base]{class}} \code{"snssde3d"}.
\item{XYZ }{an invisible \code{mts} (3-dim) object (X(t),Y(t),Z(t)).}
\item{driftx, drifty, driftz }{drift coefficient of X(t), Y(t) and Z(t).}
\item{diffx, diffy, diffz }{diffusion coefficient of X(t), Y(t) and Z(t).}
\item{type }{type of sde.}
\item{method }{the numerical method used.}
}
\references{
  Friedman, A. (1975). 
  \emph{Stochastic differential equations and applications}.
  Volume 1, ACADEMIC PRESS.
  
  Henderson, D. and Plaschko,P. (2006). 
  \emph{Stochastic differential equations in science and engineering}.
  World Scientific.
  
  Allen, E. (2007).
  \emph{Modeling with Ito stochastic differential equations}. 
  Springer-Verlag.

  Jedrzejewski, F. (2009). 
  \emph{Modeles aleatoires et physique probabiliste}. 
  Springer-Verlag.
  
  Iacus, S.M. (2008). 
  \emph{Simulation and inference for stochastic differential equations: with R examples}. 
  Springer-Verlag, New York.
  
  Kloeden, P.E, and Platen, E. (1989). 
  A survey of numerical methods for stochastic differential equations. 
  \emph{Stochastic Hydrology and Hydraulics}, \bold{3}, 155--178.

  Kloeden, P.E, and Platen, E. (1991a).
  Relations between multiple ito and stratonovich integrals.
  \emph{Stochastic Analysis and Applications}, \bold{9}(3), 311--321.

  Kloeden, P.E, and Platen, E. (1991b).
  Stratonovich and ito stochastic taylor expansions. 
  \emph{Mathematische Nachrichten}, \bold{151}, 33--50.

  Kloeden, P.E, and Platen, E. (1995).
  \emph{Numerical Solution of Stochastic Differential Equations}. 
  Springer-Verlag, New York.

  Oksendal, B. (2000).
  \emph{Stochastic Differential Equations: An Introduction with Applications}.
  5th edn. Springer-Verlag, Berlin.

  Platen, E. (1980).
  Weak convergence of approximations of ito integral equations. 
  \emph{Z Angew Math Mech}. \bold{60}, 609--614.

  Platen, E. and Bruti-Liberati, N. (2010).
  \emph{Numerical Solution of Stochastic Differential Equations with Jumps in Finance.} 
  Springer-Verlag, New York

  Saito, Y, and Mitsui, T. (1993).
  Simulation of Stochastic Differential Equations. 
  \emph{The Annals of the Institute of Statistical Mathematics}, \bold{3}, 419--432.
}
\author{
A.C. Guidoum, K. Boukhetala.
}
\note{
If you use higher order schemes as Heun or Runge-Kutta, you may have to wait long time.
}

\seealso{
\code{\link[sde]{sde.sim}} in package \CRANpkg{sde}. \code{\link[yuima]{simulate}} in package \CRANpkg{yuima}. \code{\link{snssde1d}} for 1-dim sde.
}
\examples{

## Example 1: Ito sde
## dX(t) = W2(t) dt + W3(t) * dW1(t) 
## dY(t) = dW2(t) 
## dZ(t) = dW3(t)        
## W1(t), W2(t) and W3(t) three independent Brownian motion

fx <- expression(y)
gx <- expression(z)
fy <- expression(0)
gy <- expression(1)
fz <- expression(0)
gz <- expression(1)

res <- snssde3d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,driftz=fz,diffz=gz,N=1000)
res
## Sim <- res$XYZ
dev.new()
plot(res)
dev.new()
plot(res,plot.type="single")
dev.new()
plot3D(res,display="persp") ## in space (O,X,Y,Z)

## Example 2: Stratonovich sde
## dX(t) = 2*(3-X(t))* dt + (Y(t)+Z(t)) o dW1(t)          
## dY(t) = 2*(3-Y(t))* dt + (X(t)+Z(t)) o dW2(t)
## dZ(t) = 2*(3-Z(t))* dt + (X(t)+Y(t)) o dW3(t)

fx <- expression(2*(3-x))
gx <- expression(y+z)
fy <- expression(2*(3-y))
gy <- expression(x+z)
fz <- expression(2*(3-z))
gz <- expression(x+y)

res <- snssde3d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,driftz=fz,diffz=gz,
                ,N=1000,type="str")
res
dev.new()
plot(res,plot.type="single")
dev.new()
plot(res)
dev.new()
plot3D(res,display="rgl") ## in space (O,X,Y,Z)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sde}
\keyword{ts}
\keyword{mts}
