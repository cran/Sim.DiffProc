\name{Hyproc}
\alias{Hyproc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating The Hyperbolic Process (by Milstein Scheme)
}
\description{
Simulation hyperbolic process by milstein scheme.
}
\usage{
Hyproc(N, M, t0, T, x0, theta, output = FALSE)
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
constant positive.
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\details{
A process \code{X} satisfying : \emph{\deqn{dX(t) = (-theta*X(t)/sqrt(1+X(t)^2)) *dt + dW(t)}}
With \code{(-theta*X(t)/sqrt(1+X(t)^2)) :drift coefficient} and \code{1 :diffusion coefficient}, \code{W(t)} is Wiener process, discretization \code{dt = (T-t0)/N}.
 
Constraints: \code{theta > 0}. 
}
\value{
data.frame(time,x) and plot of process.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{Hyprocg}} General Hyperbolic Diffusion, \code{\link{CIRhy}} modified CIR and hyperbolic Process, \code{\link{snssde}} Simulation Numerical Solution of SDE.
}
\examples{
\donttest{
## Hyperbolic Process 
## dX(t) = (-2*X(t)/sqrt(1+X(t)^2)) *dt + dW(t)
## One trajectorie
 Hyproc(N=1000,M=1,T=100,t0=0,x0=3,theta=2)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
\keyword{ Stochastic Differential Equation }
\keyword{Numerical Solution of Stochastic Differential Equation}

