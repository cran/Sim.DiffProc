\name{SLVM}
\alias{SLVM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stochastic Lotka-Volterra Model
}
\description{
Simulation the stochastic Lotka-Volterra model.
}
\usage{
SLVM(N, t0, T, x0, y0, a, b, c, d, sigma, Step = FALSE, Output = FALSE)
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
initial value of the process at time \code{t0} (\code{x0 > 0}).
}
  \item{y0}{
initial value of the process at time \code{t0} (\code{y0 > 0}).
}
  \item{a}{
positive parameter.
}
  \item{b}{
positive parameter.
}
  \item{c}{
positive parameter.
}
  \item{d}{
positive parameter.
}
  \item{sigma}{
positive parameter.
}
  \item{Step}{
if \code{Step = TRUE} ploting step by step.
}
  \item{Output}{
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\details{
The Lotka-Volterra system of stochastics differential equations, (Lotka (1925),Volterra (1926)):
\emph{\deqn{dX(t) = (a * X(t) - b * X(t) * Y(t)) dt + sigma* dW1(t)}}
\emph{\deqn{dY(t) = (c * X(t) * Y(t) - d * Y(t)) dt + sigma* dW2(t)}}

with positive \code{x0, y0} and positive parameters \code{a, b, c, d} describes a behaviour
of a prey-predator system in terms of the prey and predator (intensities) \code{X(t)} and \code{Y(t)}.

Here, \code{a} is the rate of increase of prey in the absence of predator, \code{d} is
a rate of decrease of predator in the absence of prey while the rate of decrease
in prey is proportional to the number of predators \code{b* Y(t)}, and similarly the rate
of increase in predator is proportional to the number of prey \code{c* X(t)}.

The system possesses the first integral which is a closed orbit in
the first quadrant of phase plane \code{x, y}. It is given by :
\emph{\deqn{r(x,y) = c * x - d * log(x) + b * y - a * log(y) + r0}}
}
\value{
data.frame(time,x,y), plot 1D and 2D of the process. 
}
\references{
Fima C Klebaner. Introduction to stochastic calculus with application (Second Edition), Imperial College Press (ICP), 2005.
}
\author{
Guidoum Arsalane.
}
\seealso{
\code{\link{WFD}} Feller Branching Diffusion, \code{\link{FBD}} Feller Branching Diffusion.
}
\examples{
\donttest{
 SLVM(N=5000,t0=0,T=100,x0=1,y0=1,a=1,b=2,c=0.5,d=0.25,sigma=0.01)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Biology modeling }
