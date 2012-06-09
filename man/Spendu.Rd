\name{Spendu}
\alias{Spendu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stochastic pendulum
}
\description{
You can see from this simulation the stochastic pendulum, the temporal graph and the phase portrait.
}
\usage{
Spendu(N, T, theta0, theta1, lambda, omega, sigma, Step = FALSE, Output = FALSE)
}
\arguments{
  \item{N}{
size of process.
}
  \item{T}{
final time.
}
  \item{theta0}{
Initial conditions, position (rad), \code{-pi < theta0 < pi}.
}
  \item{theta1}{
Initial conditions, speed (rad/s).
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
On this simulation the movement of the stochastic pendulum as well as the temporal graph and the portrait of phase. 
Cursors make it possible to modify the parameters of the oscillator (\code{lambda}: damping and \code{omega}: own pulsation), as well as the initial conditions.

Stochastic perturbations of the pendulum equation, and random excitations force of such systems by White noise \code{e(t)}, with delta-type correlation functions
\code{E(e(t)e(t+h))=sigma*deltat(h)}: \emph{\deqn{x'' + 2*lambda*x' +omega^2 * sin(x) = e(t) }}
where \code{lambda,sigma >= 0} and \code{omega > 0}. 

To observe the evolution of the portrait of phase when the initial conditions are modified:
\itemize{\item When the amplitude is large, one notices the difference in behavior with the harmonic oscillator: lengthening of the period, deformation of the graphs.
\item When the pendulum does one (or several) turn, the portrait of phase opens: is it closed again?}
}
\value{
data.frame(time,X(t)), plot of process X(t) in the phase portrait (2D) and temporal evolution of stochastic pendulum.
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
\code{\link{Sharosc}} stochastic harmonic oscillator, \code{\link{Svandp}} stochastic Van der Pol oscillator, \code{\link{Srayle}} stochastic Rayleigh oscillator,
\code{\link{SSCPP}} stochastic system with a cylindric phase plane, \code{\link{Sosadd}} stochastic oscillator with additive noise.

}
\examples{
\donttest{
## theta0= 3, theta1 = 0, lambda=0.1, omega=2, sigma=0.1
 Spendu(N=5000, T=50, theta0=3, theta1=0, lambda=0.1, omega=2, sigma=0.1)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Stochastics Oscillators  }
