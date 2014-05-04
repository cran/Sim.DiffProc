\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{fptsde2d}
\alias{fptsde2d}
\alias{fptsde2d.default}
\alias{summary.fptsde2d}
\alias{mean.fptsde2d}
\alias{median.fptsde2d}
\alias{quantile.fptsde2d}
\alias{kurtosis.fptsde2d}
\alias{skewness.fptsde2d}
\alias{moment.fptsde2d}
\alias{bconfint.fptsde2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{First Passage Time in 2-Dim SDE}
\description{The (S3) generic function \code{fptsde2d} for simulate first-passage-time (f.p.t) in 2-dim stochastic differential equations.}
\usage{
fptsde2d(N, \dots)
\method{fptsde2d}{default}(N = 100, M = 10, x0 = 0, y0 = 0, t0 = 0, T = 1, Dt, c = 0.5, 
   driftx, diffx, drifty, diffy, alpha = 0.5, mu = 0.5, type = c("ito", "str"), 
   method = c("euler", "milstein", "predcorr", "smilstein", "taylor", 
   "heun", "rk1", "rk2", "rk3"), \dots)		

   
\method{summary}{fptsde2d}(object, \dots)	
\method{mean}{fptsde2d}(x, \dots)
\method{median}{fptsde2d}(x, \dots)
\method{quantile}{fptsde2d}(x, \dots)
\method{kurtosis}{fptsde2d}(x, \dots)
\method{skewness}{fptsde2d}(x, \dots)
\method{moment}{fptsde2d}(x, order = 2, \dots)
\method{bconfint}{fptsde2d}(x, level=0.95, \dots)						  
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{N}{size of sde.}
  \item{M}{size of fpt.}
  \item{x0}{initial value of the process \eqn{X_{t}}{X(t)} at time \code{t0}.}
  \item{y0}{initial value of the process \eqn{Y_{t}}{Y(t)} at time \code{t0}.}
  \item{t0}{initial time.}
  \item{T}{final time.} 
  \item{Dt}{time step of the simulation (discretization). If it is \code{\link[base]{missing}} a default \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}.}
  \item{c}{boundary or barrier (threshold).}
  \item{driftx}{drift coefficient: an \code{\link[base]{expression}} of three variables \code{t}, \code{x} and \code{y} for process \eqn{X_t}{X(t)}.}
  \item{diffx}{diffusion coefficient: an \code{\link[base]{expression}} of three variables \code{t}, \code{x} and \code{y} for process \eqn{X_t}{X(t)}.}
  \item{drifty}{drift coefficient: an \code{\link[base]{expression}} of three variables \code{t}, \code{x} and \code{y} for process \eqn{Y_t}{Y(t)}.}
  \item{diffy}{diffusion coefficient: an \code{\link[base]{expression}} of three variables \code{t}, \code{x} and \code{y} for process \eqn{Y_t}{Y(t)}.}
  \item{alpha}{weight \code{alpha} of the predictor-corrector scheme; the default \code{alpha = 0.5}.}
  \item{mu}{weight \code{mu} of the predictor-corrector scheme; the default \code{mu = 0.5}.}
  \item{type}{sde of the type Ito or Stratonovich.}
  \item{method}{numerical methods of simulation, the default \code{method = "euler"}; see \code{\link{snssde2d}}.}
  \item{x, object}{an object inheriting from class \code{"fptsde2d"}.}
  \item{order}{order of moment.}
  \item{level}{the confidence level required.}    
  \item{\dots}{further arguments for (non-default) methods.}
}
\details{
The function \code{fptsde2d} returns a random variables \eqn{\tau_{c}(x)}{tau_x(c)} and \eqn{\tau_{c}(y)}{tau_y(c)} "first passage time" for 
\eqn{(X(t),Y(t))}{(X(t),Y(t))}, defined by : 
\deqn{\tau_{c}(x) = \{ t \geq 0 ; X_{t} \geq c \},\quad if \quad (x_{0} < c)}{tau_x(c)={t>=0; X(t) >= c}, if (x0 < c)}
\deqn{\tau_{c}(y) = \{ t \geq 0 ; Y_{t} \geq c \},\quad if \quad (y_{0} < c)}{tau_y(c)={t>=0; Y(t) >= c}, if (y0 < c)}
and 
\deqn{\tau_{c}(x) = \{ t \geq 0 ; X_{t} \leq c \},\quad if \quad (x_{0} > c)}{tau_x(c)={t>=0; X(t) <= c}, if (x0 > c)}
\deqn{\tau_{c}(y) = \{ t \geq 0 ; Y_{t} \leq c \},\quad if \quad (y_{0} > c)}{tau_y(c)={t>=0; Y(t) <= c}, if (y0 > c)}
with \eqn{c}{c} is a fixed boundary or barrier.
}
\value{
\code{fptsde2d} returns an object inheriting from \code{\link[base]{class}} \code{"fptsde2d"}.
\item{tau_x, tau_y }{a vector of couple 'fpt' \eqn{(\tau_{c}(x),\tau_{c}(y))}{tau_x(c),tau_y(c)}.}
}
\references{
  Argyrakisa, P. and G.H. Weiss (2006).
  A first-passage time problem for many random walkers.
  \emph{Physica A}. \bold{363}, 343--347.
  
  Aytug H., G. J. Koehler (2000). 
  New stopping criterion for genetic algorithms. 
  \emph{European Journal of Operational Research}, \bold{126}, 662--674.
  
  Boukhetala, K. (1996)
  Modelling and simulation of a dispersion pollutant with attractive centre.
  ed by Computational Mechanics Publications, Southampton ,U.K and Computational Mechanics Inc, Boston, USA, 245--252.
  
  Boukhetala, K. (1998a). 
  Estimation of the first passage time distribution for a simulated diffusion process.
  \emph{Maghreb Math.Rev}, \bold{7}(1), 1--25.
  
  Boukhetala, K. (1998b). 
  Kernel density of the exit time in a simulated diffusion. 
  \emph{les Annales Maghrebines De L ingenieur}, \bold{12}, 587--589.
  
  Ding, M. and G. Rangarajan. (2004). 
  First Passage Time Problem: A Fokker-Planck Approach. 
  \emph{New Directions in Statistical Physics}. ed by L. T. Wille. Springer. 31--46.
  
   Roman, R.P., Serrano, J. J., Torres, F. (2008).
   First-passage-time location function: Application to determine first-passage-time densities in diffusion processes.
   \emph{Computational Statistics and Data Analysis}. \bold{52}, 4132--4146.
  
   Roman, R.P., Serrano, J. J., Torres, F. (2012). 
   An R package for an efficient approximation of first-passage-time densities for diffusion processes based on the FPTL function. 
   \emph{Applied Mathematics and Computation}, \bold{218}, 8408--8428.
   
   Roman, R.P., Serrano, J. J., Torres, F. (2013). 
   \CRANpkg{fptdApprox}: Approximation of first-passage-time densities for diffusion processes. 
   \emph{\R package version 2.0}. 
  
  Gardiner, C. W. (1997). 
  \emph{Handbook of Stochastic Methods}.
  Springer-Verlag, New York.
}
\author{
A.C. Guidoum, K. Boukhetala.
}

\seealso{
\code{\link[Sim.DiffProc]{fptsde1d}} for simulation fpt in sde 1-dim. \code{\link[fptdApprox]{FPTL}} in package \CRANpkg{fptdApprox} for computes values of the first passage time location. 
}
\examples{

## Example 1: Ito sde 
## dX(t) = 5*(-1-Y(t))*X(t) * dt + 0.5 * dW1(t)          
## dY(t) = 5*(-1-X(t))*Y(t) * dt + 0.5 * dW2(t)
## x0 = 2, y0 = 2, and barrier c = 0.
## W1(t) and W2(t) two independent Brownian motion

fx <- expression(5*(-1-y)*x)
gx <- expression(0.5)
fy <- expression(5*(-1-x)*y)
gy <- expression(0.5)

res1 <- fptsde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,x0=2,y0=2,c=0,M=30,N=1000)
res1
## taux <- res1$tau_x
## tauy <- res1$tau_y
summary(res1)
bconfint(res1,level=0.95)
moment(res1,order=c(2,3,4,5))
X1 <- cbind(res1$tau_x,res1$tau_y)
## library(sm)
## sm.density(X1,display="persp")

## Example 2: Stratonovich sde 
## dX(t) = 5*(-1-Y(t))*X(t) * dt + 0.5 o dW1(t)          
## dY(t) = 5*(-1-X(t))*Y(t) * dt + 0.5 o dW2(t)
## x0 = 2, y0 = 2, and barrier c = 0.
## W1(t) and W2(t) two independent Brownian motion

fx <- expression(5*(-1-y)*x)
gx <- expression(0.5)
fy <- expression(5*(-1-x)*y)
gy <- expression(0.5)

res2 <- fptsde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,x0=2,y0=2,c=0,
                 M=30,N=1000,type="str")
res2
summary(res2)
moment(res2,order=c(2,3,4,5))
bconfint(res2,level=0.95)
X2 <- cbind(res2$tau_x,res2$tau_y)
## sm.density(X2,display="persp")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{fpt}
\keyword{sde}
\keyword{ts}
\keyword{mts}
