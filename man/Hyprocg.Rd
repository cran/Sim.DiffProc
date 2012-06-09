\name{Hyprocg}
\alias{Hyprocg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating The General Hyperbolic Diffusion (by Milstein Scheme)
}
\description{
Simulation the general hyperbolic diffusion by milstein scheme.
}
\usage{
Hyprocg(N, M, t0, T, x0, beta, gamma, theta, mu, sigma, output = FALSE)
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
  \item{beta}{
constant (\code{0.5*sigma^2*(beta-(gamma*X(t))/sqrt(theta^2+(X(t)-mu)^2)):drift coefficient}). 
}
  \item{gamma}{
constant positive (\code{0.5*sigma^2*(beta-(gamma*X(t))/sqrt(theta^2+(X(t)-mu)^2)):drift coefficient}). 
}
  \item{theta}{
constant positive (\code{0.5*sigma^2*(beta-(gamma*X(t))/sqrt(theta^2+(X(t)-mu)^2)):drift coefficient}). 
}
  \item{mu}{
constant (\code{0.5*sigma^2*(beta-(gamma*X(t))/sqrt(theta^2+(X(t)-mu)^2)):drift coefficient}). 
}
  \item{sigma}{
constant positive ( \code{sigma :diffusion coefficient}).
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\details{
A process \code{X} satisfying : \emph{\deqn{dX(t) = (0.5*sigma^2*(beta-(gamma*X(t))/sqrt(theta^2+(X(t)-mu)^2))*dt + dW(t)}}
With \code{(0.5*sigma^2*(beta-(gamma*X(t))/sqrt(theta^2+(X(t)-mu)^2)):drift coefficient} and \code{sigma :diffusion coefficient}, \code{W(t)} is Wiener process, discretization \code{dt = (T-t0)/N}.
 
The parameters \code{gamma > 0} and \code{0 <=  abs(beta)  < gamma} determine the shape of the distribution, and \code{theta >= 0}, and \code{mu} are, respectively, the scale and location parameters of the distribution. 

Constraints: \code{gamma > 0} , \code{0 <=  abs(beta)  < gamma} , \code{theta >= 0} , \code{sigma > 0}.   
}
\value{
data.frame(time,x) and plot of process.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{Hyproc}} Hyperbolic Process, \code{\link{CIRhy}} modified CIR and hyperbolic Process, \code{\link{snssde}} Simulation Numerical Solution of SDE.
}
\examples{
\donttest{
## Hyperbolic Process 
## dX(t) = 0.5 * (2)^2*(0.25-(0.5*X(t))/sqrt(2^2+(X(t)-1)^2)) *dt + 2* dW(t)
## One trajectorie
 Hyprocg(N=1000,M=1,T=100,t0=0,x0=-10,beta=0.25,gamma=0.5,theta=2,mu=1,sigma=2)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
\keyword{ Stochastic Differential Equation }
\keyword{Numerical Solution of Stochastic Differential Equation}

