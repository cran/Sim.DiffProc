\name{FBD}
\alias{FBD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Feller Branching Diffusion
}
\description{
Simulation the Feller Branching diffusion.
}
\usage{
FBD(N, M, t0, T, x0, mu, sigma, output = FALSE)
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
  \item{mu}{
constant (\code{mu * X(t) :drift coefficient}). 
}
  \item{sigma}{
constant positive (\code{sigma * sqrt(X(t)) :diffusion coefficient}).
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\details{
The Feller Branching diffusion model also derives directly from the linear drift class, the discretization \code{dt = (T-t0)/N}. 

A simple branching process is a model in which individuals reproduce independently
of each other and of the history of the process. The continuous
approximation to branching process is the branching diffusion. It is given by
the stochastic differential equation for the population size X(t), 0 < X(t) < +Inf : \emph{\deqn{dX(t) = mu * X(t)* dt + sigma * sqrt(X(t)) *dW(t)}} with \code{mu * X(t) :drift coefficient} and \code{sigma * sqrt(X(t)) :diffusion coefficient}, \code{W(t)} is Wiener process.
}
\value{
data.frame(time,x) and plot of process. 
}
\references{
Fima C Klebaner. Introduction to stochastic calculus with application (Second Edition), Imperial College Press (ICP), 2005.
}
\author{
Guidoum Arsalane.
}
\seealso{
\code{\link{SLVM}} Stochastic Lotka-Volterra, \code{\link{WFD}} Wright-Fisher Diffusion.
}
\examples{
\donttest{
 FBD(N=1000,M=1,t0=0,T=1,x0=1,mu=2,sigma=0.5,output=FALSE)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Biology modeling }