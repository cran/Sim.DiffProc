\name{HWVF}
\alias{HWVF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Flow of Hull-White/Vasicek (HWV) Gaussian Diffusion Models
}
\description{
Simulation flow of the Hull-White/Vasicek or gaussian diffusion models.
}
\usage{
HWVF(N, M, t0, T, x0, theta, r, sigma, output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of process.
}
  \item{M}{
number of trajectories.
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
\code{\link{HWV}} Hull-White/Vasicek Models, \code{\link{PEOUG}} Parametric Estimation of Hull-White/Vasicek Models, \code{\link{snssde}} Simulation Numerical Solution of SDE.
}
\examples{
\donttest{
## flow of Hull-White/Vasicek Models
## dX(t) = 4 * (2.5 - X(t)) * dt + 1 *dW(t)
 HWVF(N=1000,M=10,t0=0,T=1,x0=10,theta=2.5,r=4,sigma=1)
## if theta = 0 than "OUF" = "HWVF"
## dX(t) = 4 * ( 0 - X(t)) * dt + 1 *dW(t)
 system.time(HWVF(N=1000,M=10,t0=0,T=1,x0=10,theta=0,r=4,sigma=1))
 system.time(OUF(N=1000,M=5,t0=0,T=1,x0=10,r=4,sigma=1))}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
\keyword{ Stochastic Differential Equation }
