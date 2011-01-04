\name{OU}
\alias{OU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Ornstein-Uhlenbeck Process
}
\description{
Simulation the ornstein-uhlenbeck or Hull-White/Vasicek model.
}
\usage{
OU(N, t0, T, x0, r, sigma, output = FALSE)
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
  \item{r}{
constant positive (\code{r is speed of reversion} and  \code{-r * X(t) :drift coefficient}). 
}
  \item{sigma}{
constant positive (\code{sigma (volatility) :diffusion coefficient}).
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel 2007.
}
}
\details{
The Ornstein-Uhlenbeck or Vasicek process is the unique solution to the following stochastic differential equation : \emph{\deqn{dX(t) = - r * X(t) * dt + sigma * dW(t)}}  With \code{-r * X(t) :drift coefficient} and \code{sigma  : diffusion coefficient}, \code{W(t)} is Wiener process, the discretization \code{dt = (T-t0)/N}.

Please note that the process is stationary only if \code{r > 0}. 
}
\value{
data.frame(time,x) and plot of process.
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\seealso{
\code{\link{OUF}} Flow of Ornstein-Uhlenbeck Process, \code{\link{PEOU}} Parametric Estimation of Ornstein-Uhlenbeck Model, \code{\link{PEOUexp}} Explicit Estimators of Ornstein-Uhlenbeck Model, \code{\link{snssde}} Simulation Numerical Solution of SDE.
}
\examples{
## Ornstein-Uhlenbeck Process
## dX(t) = -2 * X(t) * dt + 1 *dW(t)
OU(N=1000,t0=0,T=10,x0=10,r=2,sigma=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
\keyword{ Stochastic Differential Equation }
