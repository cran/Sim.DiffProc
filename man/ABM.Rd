\name{ABM}
\alias{ABM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Arithmetic Brownian Motion Model
}
\description{
Simulation of the arithmetic brownian motion model.
}
\usage{
ABM(N, t0, T, x0, theta, sigma, output = FALSE)
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
The function \code{ABM} returns a trajectory of the Arithmetic Brownian motion starting at \code{x0} at time \code{t0}, than the Discretization \code{dt = (T-t0)/N}.

The stochastic differential equation of the Arithmetic Brownian motion is : \emph{\deqn{dX(t) = theta * dt + sigma * dW(t)}} 
with \code{theta :drift coefficient} and \code{sigma :diffusion coefficient},\code{W(t)} is Wiener process.
}
\value{
data.frame(time,x) and plot of process. 
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{ABMF}} creating flow of the arithmetic brownian motion model.
}
\examples{
## Arithmetic Brownian Motion Model
## dX(t) = 3 * dt + 2 * dW(t) ; x0 = 0 and t0 = 0
\donttest{ ABM(N=1000,t0=0,T=1,x0=0,theta=3,sigma=2)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Actuarial Modeling }
\keyword{ Financial models }
\keyword{ Environment R }
\keyword{ Simulation }


