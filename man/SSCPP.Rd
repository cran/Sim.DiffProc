\name{SSCPP}
\alias{SSCPP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stochastic system with a cylindric phase plane
}
\description{
You can see from this simulation the stochastic system with a cylindric phase plane and the temporal graph and the phase portrait, and 3D plot for Fokker-Planck equation.
}
\usage{
SSCPP(N, T, theta0, theta1, a, b, omega, sigma, K0 = 1, Prd = 6, Step = FALSE, Output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{a}{
Amortization (>= 0).
}
  \item{b}{
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
  \item{Prd}{
Period for plot 3D  (\code{Prd > 0}).
}
  \item{Step}{
if \code{Step = TRUE} ploting step by step.
}
  \item{Output}{
If \code{Output = yes} write a output to an Excel (.csv).
}
}
\details{
Stochastic perturbations of the system with a cylindric phase plane equation, and random excitations force of such systems by White noise \code{e(t)}, with delta-type correlation functions
\code{E(e(t)e(t+h))=sigma*deltat(h)}: \emph{\deqn{x'' + a * x' + b +omega^2 * sin(x) = e(t) }}
where \code{a,b,omega,sigma >= 0}.

The Fokker-Planck equation of this system:
\emph{\deqn{P(s,x,t,y) = exp( -a * (y^2 + 2 * b *x - 2 * omega^2 * cos(x)) / (2*pi*K0) ) }}
}
\value{
data.frame(time,X(t)), plot of process X(t) in the phase portrait (2D) and temporal evolution of stochastic Rayleigh equation. 3D plot for Fokker-Planck equation.
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
\code{\link{Srayle}} stochastic Rayleigh oscillator, \code{\link{Sosadd}} stochastic oscillator with additive noise.
}
\examples{
\donttest{
## a = 0.1, b = 0.15, omega= 2, sigma = 0.2, K0 = 3, Prd = 6
SSCPP(N=5000, T=50, theta0=3, theta1=0, a=0.1, b=0.15, omega=2, sigma=0.2, K0 = 3)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Stochastics Oscillators }
\keyword{Fokker-Planck equation}
