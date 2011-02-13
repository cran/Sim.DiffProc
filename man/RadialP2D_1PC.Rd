\name{RadialP2D_1PC}
\alias{RadialP2D_1PC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Two-Dimensional Attractive Model in Polar Coordinates Model(S = 1,Sigma)
}
\description{
Simulation 2-dimensional attractive model (S = 1) in polar coordinates.
}
\usage{
RadialP2D_1PC(N, R0, t0, T, ThetaMax, K, sigma, output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of process.
}
  \item{R0}{
initial value\code{R0 > 0} at time \code{t0}.
}
  \item{t0}{
initial time.
}
  \item{T}{
final time.
}
  \item{ThetaMax}{
polar coordinates, example \code{ThetaMax = 2*pi}.
}
  \item{K}{
constant \code{K > 0}.
}
  \item{sigma}{
constant \code{sigma > 0}.
}
  \item{output}{
if \code{Output = TRUE} write a \code{Output} to an Excel 2007.
}
}
\details{
The attractive models is defined by the system for stochastic differential equation Two-dimensional : 
\emph{\deqn{dX(t) = (-K * X(t)/(sqrt(X(t)^2 + Y(t)^2))^(S+1) )* dt + Sigma* dW1(t)}}
\emph{\deqn{dY(t) = (-K * Y(t)/(sqrt(X(t)^2 + Y(t)^2))^(S+1) )* dt + Sigma* dW2(t)}}
\code{dW1(t)} and \code{dW2(t)} are brownian motions independent. 

Using Ito transform, it is shown that the \code{Radial Process R(t)} with \code{R(t)=||(X(t),Y(t))||} is a markovian diffusion, solution of the stochastic differential equation one-dimensional:
\emph{\deqn{dR(t) = ((0.5 * Sigma^2 * R(t)^(S-1) - K)/ R(t)^S )* dt + Sigma* dW(t)}}

If \code{S = 1} (ie \code{M(S=1,Sigma)}) the \code{R(t)} is :
\emph{\deqn{dR(t) = ((0.5*Sigma^2 -K )/R(t) ) * dt + Sigma* dW(t)}}
Where ||.|| is the Euclidean norm and \code{dW(t)} is a determined brownian motions.

\code{R(t)=sqrt(X(t)^2 + Y(t)^2)} it is distance between \code{X(t)} and \code{Y(t)}, then \code{X(t)=R(t)*cos(theta(t))}
and \code{Y(t)=R(t)*sin(theta(t))}, 

For more detail consulted \code{References}.
}
\value{
data.frame(time,R(t),theta(t)) and plot of process 2-D in polar coordinates.
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
\code{\link{snssde2D}}, \code{\link{PredCorr2D}}, \code{\link{RadialP2D_2PC}}, \code{\link{RadialP3D_1}}, \code{\link{tho_M1}}, \code{\link{fctgeneral}}, \code{\link{hist_general}}, \code{\link{Kern_meth}}.
}
\examples{
RadialP2D_1PC(N=1000, R0=3, t0=0, T=1, ThetaMax=4*pi, K=2, sigma=1, 
              output = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Attractive Model}
\keyword{ Diffusion Process Multidimensional}
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Stochastic Differential Equation Multidimensional}
