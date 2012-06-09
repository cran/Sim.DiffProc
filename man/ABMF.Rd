\name{ABMF}
\alias{ABMF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Flow of The Arithmetic Brownian Motion Model
}
\description{
Simulation flow of the arithmetic brownian motion model.
}
\usage{
ABMF(N, M, t0, T, x0, theta, sigma, output = FALSE)
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
constant (\code{Coefficient of drift}).
}
  \item{sigma}{
constant positive (\code{Coefficient of diffusion}).
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\details{
The function \code{ABMF} returns a flow of the Arithmetic Brownian motion starting at \code{x0} at time \code{t0}, than the discretization \code{dt = (T-t0)/N}.

The stochastic differential equation of the Arithmetic Brownian motion is : \emph{\deqn{dX(t) = theta * dt + sigma * dW(t)}} 
With \code{theta :drift coefficient} and \code{sigma :diffusion coefficient}, \code{W(t)} is Wiener process.
}
\value{
data.frame(time,x) and plot of process. 
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{ABM}} creating the arithmetic brownian motion model.
}
\examples{
## Flow of Arithmetic Brownian Motion Model
## dX(t) = 3 * dt + 2 * dW(t) ; x0 = 0 and t0 = 0
\donttest{ ABMF(N=1000,M=5,t0=0,T=1,x0=0,theta=3,sigma=2)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
\keyword{ Diffusion Process }
\keyword{ Actuarial Modeling }

