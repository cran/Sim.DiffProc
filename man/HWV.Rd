\name{HWV}
\alias{HWV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Hull-White/Vasicek (HWV) Gaussian Diffusion Models
}
\description{
Simulation the Hull-White/Vasicek or gaussian diffusion models.
}
\usage{
HWV(N, t0, T, x0, theta, r, sigma, output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of process.
}
  \item{t0}{
initial time.
}
  \item{T}{
final time.
}
  \item{x0}{
initial value of the process at time \code{t0}.
}
  \item{theta}{
constant (\code{theta is the long-run equilibrium value
of the process} and \code{r*(theta -X(t)) :drift coefficient}).
}
  \item{r}{
constant positive (\code{r is speed of reversion} and  \code{r*(theta -X(t)):drift coefficient}). 
}
  \item{sigma}{
constant positive (\code{sigma (volatility) :diffusion coefficient}).
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\details{
The Hull-White/Vasicek \code{(HWV)} short rate class derives directly from SDE with mean-reverting drift: \emph{\deqn{dX(t) = r * (theta -  X(t)) * dt + sigma * dW(t)}}
With \code{r *(theta- X(t)) :drift coefficient} and \code{sigma  : diffusion coefficient}, \code{W(t)} is Wiener process, the discretization \code{dt = (T-t0)/N}.

The process is also ergodic, and its invariant law is the \code{Gaussian density}.
}
\value{
data.frame(time,x) and plot of process.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{HWVF}} Flow of Gaussian Diffusion Models, \code{\link{PEOUG}} Parametric Estimation of Hull-White/Vasicek Models, \code{\link{snssde}} Simulation Numerical Solution of SDE.
}
\examples{
\donttest{
## Hull-White/Vasicek Models
## dX(t) = 4 * (2.5 - X(t)) * dt + 1 *dW(t)
 HWV(N=1000,t0=0,T=1,x0=10,theta=2.5,r=4,sigma=1)
## if theta = 0 than "OU" = "HWV"
## dX(t) = 4 * ( 0 - X(t)) * dt + 1 *dW(t)
 system.time(OU(N=10^4,t0=0,T=1,x0=10,r=4,sigma=1))
 system.time(HWV(N=10^4,t0=0,T=1,x0=10,theta=0,r=4,sigma=1))}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
\keyword{ Stochastic Differential Equation }
