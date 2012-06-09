\name{WFD}
\alias{WFD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wright-Fisher Diffusion
}
\description{
Simulation the Wright-Fisher diffusion.
}
\usage{
WFD(N, M, t0, T, x0, gamma1, gamma2, sigma, output = FALSE)
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
initial value of the process at time \code{t0} (\code{0 < x0 < 1}).
}
  \item{gamma1}{
constant (\code{-gamma1 * X(t) + gamma2 * (1 - X(t)) :drift coefficient}), (\code{gamma1 >= 0}). 
}
  \item{gamma2}{
constant (\code{-gamma1 * X(t) + gamma2 * (1 - X(t)) :drift coefficient}). (\code{gamma2 >= 0})  
}
  \item{sigma}{
constant positive (\code{sigma * sqrt(X(t)*(1-X(t))) :diffusion coefficient}).
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\details{
The Wright-Fisher diffusion also derives directly from the linear drift class, the discretization \code{dt = (T-t0)/N}.

In population dynamics frequencies of genes or alleles are studied. It is assumed
for simplicity that the population size \code{N} is fixed and individuals are of
two types: \code{A} and \code{a}. If individuals of type \code{A} mutate to type \code{a} with the rate \code{gamma1/N} and individuals of type \code{a} mutate to type
\code{A} with the rate \code{gamma2/N}, then it is possible to approximate the frequency of type \code{A} individuals \code{X(t)} by the
Wright-Fisher diffusion, given by the stochastic equation : 
\emph{\deqn{dX(t) = (-gamma1 * X(t) + gamma2 * ( 1 - X(t)) ) * dt + sigma * sqrt(X(t)*(1-X(t))) *dW(t)}} with \code{(-gamma1 * X(t) + gamma2 * ( 1 - X(t)) )  :drift coefficient} and \code{sigma * sqrt(X(t)*(1-X(t))) :diffusion coefficient}, \code{W(t)} is Wiener process. 
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
\code{\link{SLVM}} Stochastic Lotka-Volterra, \code{\link{FBD}} Feller Branching Diffusion.
}
\examples{
\donttest{
WFD(N=1000,M=1,t0=0,T=1,x0=0.5,gamma1=0,gamma2=0.5,sigma=0.2)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Biology modeling }
