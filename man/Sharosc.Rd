\name{Sharosc}
\alias{Sharosc}
\title{
Stochastic harmonic oscillator
}
\description{
The simulation shows the oscillations of a mass suspended from a spring. The graphs show the time evolution and the phase portrait.
}
\usage{
Sharosc(N, T, x0, v0, lambda, omega, sigma, Step = FALSE, Output = FALSE)
}

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
  \item{lambda}{
Amortization (1/s).
}
  \item{omega}{
Angular frequency (rad/s).
}
  \item{sigma}{
Dark random excitation.
}
  \item{Step}{
if \code{Step = TRUE} ploting step by step.
}
  \item{Output}{
If \code{Output = yes} write a output to an Excel (.csv).
}
}
\details{
Cursors used to vary the parameters of the oscillator (the damping \code{lambda} and natural frequency \code{omega}) and the initial conditions (position and velocity). 
To vary \code{lambda} and \code{omega}, \code{sigma} and observe the different regimes of the oscillator (pseudo-periodic, critical, supercritical).

Stochastic perturbations of the harmonic oscillator equation, and random excitations force of such systems by White noise \code{e(t)}, with delta-type correlation functions
\code{E(e(t)e(t+h))=sigma*deltat(h)}: \emph{\deqn{x'' + 2*lambda*x' +omega^2 *x = e(t) }}
where \code{lambda,sigma >= 0} and \code{omega > 0}. 
}
\value{
data.frame(time,X(t)), plot of process X(t) in the phase portrait (2D) and temporal evolution of stochastic harmonic oscillator.
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
\code{\link{Spendu}} stochastic pendulum, \code{\link{Svandp}} stochastic Van der Pol oscillator, \code{\link{Srayle}} stochastic Rayleigh oscillator,
\code{\link{SSCPP}} stochastic system with a cylindric phase plane, \code{\link{Sosadd}} stochastic oscillator with additive noise.
}
\examples{
\donttest{
## lambda = 0.1, omega = 1.5, sigma = 2.
 Sharosc(N=5000, T=50, x0=100, v0=0, lambda=0.1, omega=1.5, sigma=2)}
}
\keyword{ Stochastics Oscillators }
