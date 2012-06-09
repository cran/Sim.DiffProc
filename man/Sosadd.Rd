\name{Sosadd}
\alias{Sosadd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stochastic oscillator with additive noise
}
\description{
You can see from this simulation the stochastic oscillator with additive noise and the temporal graph and the phase portrait, and 3D plot for Fokker-Planck equation.
}
\usage{
Sosadd(N, T, x0, v0, a, omega, sigma, K0 = 1, Step = FALSE, Output = FALSE)
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
Initial conditions, position.
}
  \item{v0}{
Initial conditions, speed.
}
  \item{a}{
Constant (>= 0).
}
  \item{omega}{
Angular frequency (>= 0).
}
\item{sigma}{
Dark random excitation (>= 0).
}
  \item{K0}{
Constant for Fokker-Planck equation (\code{K0 > 0}).
}
  \item{Step}{
if \code{Step = TRUE} ploting step by step.
}
  \item{Output}{
If \code{Output = yes} write a output to an Excel (.csv).
}
}
\details{
Stochastic perturbations of oscillator with additive noise, and random excitations force of such systems by White noise \code{e(t)}, with delta-type correlation functions
\code{E(e(t)e(t+h))=sigma*deltat(h)}: \emph{\deqn{x'' - a * (1 - x^2 - x'^2) * x' + omega^2 * x = e(t) }}
where \code{a,omega,sigma >= 0}.

The Fokker-Planck equation of this system:
\emph{\deqn{P(s,x,t,y) = exp( - a * (x^2 + y^2 )^2 / (2*pi*K0) ) }}

}
\value{
data.frame(time,X(t)), plot of process X(t) in the phase portrait (2D) and temporal evolution of stochastic oscillator with additive noise. 3D plot for Fokker-Planck equation.
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
\code{\link{Srayle}} stochastic Rayleigh oscillator, \code{\link{SSCPP}} Stochastic system with a cylindric phase plane.
}
\examples{
\donttest{
## a = 0.1, omega= 1, sigma = 0.2 ,K0 = 0.5.
 Sosadd(N=5000, T=50, x0=3, v0=0, a=0.1, omega=1, sigma=0.2, K0 = 0.5)
## a = 3
 Sosadd(N=5000, T=50, x0=3, v0=0, a=3, omega=1, sigma=0.2, K0 = 0.5)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Stochastics Oscillators }
\keyword{Fokker-Planck equation}
