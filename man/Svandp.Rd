\name{Svandp}
\alias{Svandp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stochastic Van der Pol oscillator
}
\description{
The stochastic Van Der Pol equation is used to model oscillators maintained. It is not linear and has no explicit solution, you can see from this simulation the stochastic Van der Pol equation, the temporal graph and the phase portrait.
}
\usage{
Svandp(N, T, x0, v0, a, b, omega, sigma, Step = FALSE, Output = FALSE)
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
reaction parameter (>= 0).
}
  \item{b}{
control parameter (> 0).
}
  \item{omega}{
Angular frequency (>= 0).
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
The stochastic equation of Van Der pol, is used to model maintained oscillators. It is not linear and does not have an explicit solution.
In this simulation makes it possible to vary these parameters (cursors), as well as the initial conditions \code{x0} and \code{v0} (ringed points).

Stochastic perturbations of the  Van Der pol equation, and random excitations force of such systems by White noise \code{e(t)}, with delta-type correlation functions
\code{E(e(t)e(t+h))=sigma*deltat(h)}: \emph{\deqn{x'' + a * x' * ( x^2 /b  - 1 ) + omega^2 * x = e(t) }}
where \code{a,omega,sigma >= 0} and \code{b > 0 }.

\itemize{\item Influence initial conditions: the oscillations occur even with low values, then are stabilized. The portrait of phase shows a limiting cycle, which does not depend on the initial conditions.
\item Influence reaction: when \code{a=0} one obtains the stochastic harmonic oscillator \code{\link{Sharosc}}; the amplitude of the oscillations depends on the initial conditions. By increasing \code{a} one notes an increasingly important deformation of the oscillations and portrait of phase. 
\item Influence control: the coefficient \code{b} determines the amplitude of the oscillations: when \code{|x|<b}, the reaction is positive and the amplitude increases. When \code{|x|>b} it is the reverse which occurs. The amplitude is stabilized around \code{2b}.}
}
\value{
data.frame(time,X(t)), plot of process X(t) in the phase portrait (2D) and temporal evolution of stochastic van der Pol equation.
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
\code{\link{Spendu}} stochastic pendulum, \code{\link{Sharosc}} stochastic harmonic oscillator, \code{\link{Srayle}} stochastic Rayleigh oscillator,
\code{\link{SSCPP}} stochastic system with a cylindric phase plane, \code{\link{Sosadd}} stochastic oscillator with additive noise.
}
\examples{
\donttest{
## a = 0, b = 0.3, omega= 2.5, sigma=0.1
 Svandp(N=10000, T=100, x0=1, v0=0, a=0, b=0.3, omega=2.5, sigma=0.1)
## a = 3
 Svandp(N=10000, T=100, x0=1, v0=0, a=3, b=0.3, omega=2.5, sigma=0.1)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Stochastics Oscillators  }
