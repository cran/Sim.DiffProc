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
\alias{plot.fptsde2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{First Passage Time in 2-Dim SDE}
\description{The (S3) generic function \code{fptsde2d} for simulate first-passage-time (f.p.t) in 2-dim stochastic differential equations.}
\usage{
fptsde2d(N, \dots)
\method{fptsde2d}{default}(N = 1000, M = 100, x0 = 0, y0 = 0, t0 = 0, T = 1, Dt, 
   boundary, driftx, diffx, drifty, diffy, alpha = 0.5, mu = 0.5, type = 
   c("ito", "str"), method = c("euler", "milstein", "predcorr", "smilstein", 
   "taylor", "heun", "rk1", "rk2", "rk3"), \dots)		

   
\method{summary}{fptsde2d}(object, \dots)	
\method{mean}{fptsde2d}(x, \dots)
\method{median}{fptsde2d}(x, \dots)
\method{quantile}{fptsde2d}(x, \dots)
\method{kurtosis}{fptsde2d}(x, \dots)
\method{skewness}{fptsde2d}(x, \dots)
\method{moment}{fptsde2d}(x, order = 2, \dots)
\method{bconfint}{fptsde2d}(x, level=0.95, \dots)
\method{plot}{fptsde2d}(x, \dots)						  
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{N}{size of sde.}
  \item{M}{size of fpt.}
  \item{x0, y0}{initial value of the process \eqn{X_{t}}{X(t)} and \eqn{Y_{t}}{Y(t)} at time \code{t0}.}
  \item{t0}{initial time.}
  \item{T}{final time.} 
  \item{Dt}{time step of the simulation (discretization). If it is \code{\link[base]{missing}} a default \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}.}
  \item{boundary}{an \code{\link[base]{expression}} of a constant or time-dependent boundary.}
  \item{driftx, drifty}{drift coefficient: an \code{\link[base]{expression}} of three variables \code{t}, \code{x} and \code{y} for process \eqn{X_t}{X(t)} and \eqn{Y_t}{Y(t)}.}
  \item{diffx, diffy}{diffusion coefficient: an \code{\link[base]{expression}} of three variables \code{t}, \code{x} and \code{y} for process \eqn{X_t}{X(t)} and \eqn{Y_t}{Y(t)}.}
  \item{alpha, mu}{weight of the predictor-corrector scheme; the default \code{alpha = 0.5} and \code{mu = 0.5}.}
  \item{type}{sde of the type Ito or Stratonovich.}
  \item{method}{numerical methods of simulation, the default \code{method = "euler"}; see \code{\link{snssde2d}}.}
  \item{x, object}{an object inheriting from class \code{"fptsde2d"}.}
  \item{order}{order of moment.}
  \item{level}{the confidence level required.}    
  \item{\dots}{further arguments for (non-default) methods.}
}
\details{
The function \code{fptsde1d} returns a random variable \eqn{(\tau_{(X(t),S(t))},\tau_{(Y(t),S(t))})}{(tau(X(t),S(t)),tau(Y(t),S(t)))} "first passage time", is defined as : 
\deqn{\tau_{(X(t),S(t))} = \{ t \geq 0 ; X_{t} \geq S(t) \},\quad if \quad X(t_{0}) < S(t_{0})}{tau(X(t),S(t))={t>=0; X(t) >= S(t)}, if X(t0) < S(t0)}
\deqn{\tau_{(Y(t),S(t))} = \{ t \geq 0 ; Y_{t} \geq S(t) \},\quad if \quad Y(t_{0}) < S(t_{0})}{tau(Y(t),S(t))={t>=0; Y(t) >= S(t)}, if Y(t0) < S(t0)}
and:
\deqn{\tau_{(X(t),S(t))} = \{ t \geq 0 ; X_{t} \leq S(t) \},\quad if \quad X(t_{0}) > S(t_{0})}{tau(X(t),S(t))={t>=0; X(t) <= S(t)}, if X(t0) > S(t0)}
\deqn{\tau_{(Y(t),S(t))} = \{ t \geq 0 ; Y_{t} \leq S(t) \},\quad if \quad Y(t_{0}) > S(t_{0})}{tau(Y(t),S(t))={t>=0; Y(t) <= S(t)}, if Y(t0) > S(t0)}
with \eqn{S(t)}{S(t)} is through a continuous boundary (barrier).
}
\value{
\code{fptsde2d} returns an object inheriting from \code{\link[base]{class}} \code{"fptsde2d"}.
\item{fptx, fpty }{a vector of couple 'fpt' \eqn{(\tau_{(X(t),S(t))},\tau_{(Y(t),S(t))})}{(tau(X(t),S(t)),tau(Y(t),S(t)))}.}
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
    
  Gardiner, C. W. (1997). 
  \emph{Handbook of Stochastic Methods}.
  Springer-Verlag, New York.
}
\author{
A.C. Guidoum, K. Boukhetala.
}

\seealso{
\code{\link[Sim.DiffProc]{fptsde1d}} for simulation fpt in sde 1-dim. 

\code{\link[fptdApprox]{FPTL}} for computes values of the first passage time location (FPTL) function, and \code{\link[fptdApprox]{Approx.fpt.density}} 
for approximate first-passage-time (f.p.t.) density in package \CRANpkg{fptdApprox}. 
}
\examples{

## dX(t) = 5*(-1-Y(t))*X(t) * dt + 0.5 * dW1(t)          
## dY(t) = 5*(-1-X(t))*Y(t) * dt + 0.5 * dW2(t)
## x0 = 2, y0 = -2, and barrier -3+5*t.
## W1(t) and W2(t) two independent Brownian motion
set.seed(1234)

fx <- expression(5*(-1-y)*x)
gx <- expression(0.5)
fy <- expression(5*(-1-x)*y)
gy <- expression(0.5)

St <- expression(-3+5*t)

res <- fptsde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,boundary=St,
                x0=2,y0=-2)
res
summary(res)
plot(res)
##

fptx <- res$fptx
fpty <- res$fpty
X1 <- cbind(fptx,fpty)
## library(sm)
## sm.density(X1,display="persp")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{fpt}
\keyword{sde}
\keyword{ts}
\keyword{mts}
