\name{GBMF}
\alias{GBMF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Flow of Geometric Brownian Motion Models
}
\description{
Simulation flow of geometric brownian motion or Black-Scholes models.
}
\usage{
GBMF(N, M, t0, T, x0, theta, sigma, output = FALSE)
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
initial value of the process at time \code{t0} (\code{x0 > 0}).
}
  \item{theta}{
constant (\code{theta is the constant interest rateand} and \code{theta * X(t) :drift coefficient}).
}
  \item{sigma}{
constant positive (\code{sigma is volatility of risky activities} and \code{sigma * X(t):diffusion coefficient}).
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\details{
This process is sometimes called the Black-Scholes-Merton model after its introduction in the finance context to model asset prices. 

The process is the solution to the stochastic differential equation : \emph{\deqn{dX(t) = theta * X(t)* dt + sigma * X(t)* dW(t)}}
With \code{theta * X(t) :drift coefficient} and \code{sigma * X(t)  : diffusion coefficient}, \code{W(t)} is Wiener process, the discretization \code{dt = (T-t0)/N}.

\code{sigma > 0}, the parameter \code{theta} is interpreted as the constant interest rate and \code{sigma} as the volatility of risky activities. 

The explicit solution is : \emph{\deqn{X(t) = x0 * exp((theta - 0.5 * sigma^2) * t + sigma * W(t))}} 

The conditional density function is \code{log-normal}.
}
\value{
data.frame(time,x) and plot of process.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{GBM}} Geometric Brownian Motion, \code{\link{PEBS}} Parametric Estimation of Model Black-Scholes, \code{\link{snssde}} Simulation Numerical Solution of SDE.
}
\examples{
\donttest{
## Flow of Black-Scholes Models
## dX(t) = 4 * X(t) * dt + 2 * X(t) *dW(t)
 GBMF(N=1000,M=5,T=1,t0=0,x0=1,theta=4,sigma=2)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
\keyword{ Stochastic Differential Equation }
