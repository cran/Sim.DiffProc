\name{RadialP3D_2}
\alias{RadialP3D_2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Three-Dimensional Attractive Model Model(S >= 2,Sigma)
}
\description{
Simulation 3-dimensional attractive model (S >= 2).
}
\usage{
RadialP3D_2(N, t0, Dt, T = 1, X0, Y0, Z0, v, K, s, Sigma, 
            Output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of process.
}
  \item{t0}{
initial time.
}
  \item{Dt}{
time step of the simulation (\code{discretization}).
}
  \item{T}{
final time.
}
  \item{X0}{
initial value of the process \code{X(t)} at time \code{t0}.
}
  \item{Y0}{
initial value of the process \code{Y(t)} at time \code{t0}.
}
  \item{Z0}{
initial value of the process \code{Z(t)} at time \code{t0}.
}
  \item{v}{
threshold. \code{0 < v < sqrt(X0^2 + Y0 ^2 + Z0^2)}
}
  \item{K}{
constant \code{K > 0}.
}
  \item{s}{
constant \code{s >= 2}.
}
  \item{Sigma}{
constant \code{Sigma > 0}.
}
  \item{Output}{
if \code{Output = TRUE} write a \code{Output} to an Excel 2007.
}
}
\details{
The attractive models is defined by the system for stochastic differential equation three-dimensional : 
\emph{\deqn{dX(t) = (-K * X(t)/(sqrt(X(t)^2 + Y(t)^2 + Z(t)^2))^(S+1) )* dt + Sigma* dW1(t)}}
\emph{\deqn{dY(t) = (-K * Y(t)/(sqrt(X(t)^2 + Y(t)^2 + Z(t)^2))^(S+1) )* dt + Sigma* dW2(t)}}
\emph{\deqn{dZ(t) = (-K * Z(t)/(sqrt(X(t)^2 + Y(t)^2 + Z(t)^2))^(S+1) )* dt + Sigma* dW3(t)}}
\code{dW1(t)}, \code{dW2(t)} and \code{dW3(t)} are brownian motions independent. 

For more detail consulted \code{References}.
}
\value{
data.frame(time,X(t),Y(t),Z(t)) and plot of process 3-D.
}
\references{
\enumerate{\item K.Boukhetala, Estimation of the first passage time distribution for a simulated diffusion process, Maghreb Math.Rev, Vol.7, No 1, Jun 1998, pp. 1-25. \item K.Boukhetala, Simulation study of a dispersion about an attractive centre. In proceedings of 11th Symposium Computational Statistics, edited by R.Dutter and W.Grossman, Wien , Austria, 1994, pp. 128-130. \item K.Boukhetala,Modelling and simulation of a dispersion pollutant with attractive centre, Edited by Computational Mechanics Publications, Southampton ,U.K and Computational Mechanics Inc, Boston, USA, pp. 245-252.\item K.Boukhetala, Kernel density of the exit time in a simulated diffusion, les Annales Maghrebines De L ingenieur, Vol , 12, N Hors Serie. Novembre 1998, Tome II, pp 587-589. 
}
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\note{
\itemize{\item \code{2*K > Sigma^2}.
}
}
\seealso{
\code{\link{RadialP3D_1}}.
}
\examples{
RadialP3D_2(N=1000, t0=0, Dt=0.001, T = 1, X0=1, Y0=0.5, Z0=0.5,
            v=0.2,K=3,s=2,Sigma=0.2, Output = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Attractive Model}
\keyword{ Diffusion Process Multidimensional}
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Stochastic Differential Equation Multidimensional}
