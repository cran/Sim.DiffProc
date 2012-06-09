\name{Srayle}
\alias{Srayle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stochastic Rayleigh oscillator
}
\description{
The stochastic Rayleigh oscillator is much like the stochastic van Der Pol oscillator save one key difference: as voltage increases, the van Der Pol oscillator increases in frequency while the Rayleigh oscillator increases in amplitude.
You can see from this simulation the stochastic Rayleigh oscillator, the temporal graph and the phase portrait.
}
\usage{
Srayle(N, T, x0, v0, a, omega, sigma, Step = FALSE, Output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of process.
}
  \item{T}{
final time.
}
  \item{x0}{
Initial conditions, position (mm).
}
  \item{v0}{
Initial conditions, speed (mm/s).
}
  \item{a}{
Amortization (1/s).
}
  \item{omega}{
Angular frequency (rad/s).
}
  \item{sigma}{
Dark random excitation (>= 0).
}
  \item{Step}{
if \code{Step = TRUE} ploting step by step.
}
  \item{Output}{
If \code{Output = yes} write a output to an Excel (.csv).
}
}
\details{
Here is the second order differential equation for the stochastic Rayleigh oscillator :
\emph{\deqn{x'' - a * ( 1 - x'^2 ) * x' + omega^2 * x = e(t) }}
where \code{a,omega > 0} and \code{sigma >= 0}.

Like the stochastic van Der Pol oscillator \code{\link{Svandp}}, \code{omega} controls how much voltage is injected into the system. \code{a} controls the way in which voltage flows through the system. 
}
\value{
data.frame(time,X(t)), plot of process X(t) in the phase portrait (2D) and temporal evolution of stochastic Rayleigh equation.
}
\references{
Fima C Klebaner. Introduction to stochastic calculus with application (Second Edition), Imperial College Press (ICP), 2005.
}
\author{
Guidoum Arsalane.
}
\note{
\itemize{\item If \code{sigma = 0} is a determinist system.
         \item Time step of the simulation \code{T/N}.}
}
\seealso{
\code{\link{Spendu}} stochastic pendulum, \code{\link{Sharosc}} stochastic harmonic oscillator, \code{\link{Svandp}} stochastic Van der Pol oscillator,
\code{\link{SSCPP}} stochastic system with a cylindric phase plane, \code{\link{Sosadd}} stochastic oscillator with additive noise.
}
\examples{
\donttest{
## a= 4 , omega= 1, sigma =0.1
Srayle(N=5000, T=50, x0=3, v0=0, a=4, omega=1, sigma=0.1)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Stochastics Oscillators  }
