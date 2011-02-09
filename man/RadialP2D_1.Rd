\name{RadialP2D_1}
\alias{RadialP2D_1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Two-Dimensional Attractive Model Model(S = 1,Sigma)  
}
\description{
Simulation 2-dimensional attractive model (S = 1).
}
\usage{
RadialP2D_1(N, t0, Dt, T = 1, X0, Y0, v, K, Sigma, Output = FALSE)
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
  \item{v}{
threshold. \code{0 < v < sqrt(X0^2 + Y0 ^2)}
}
  \item{K}{
constant \code{K > 0}.
}
  \item{Sigma}{
constant \code{Sigma > 0}.
}
  \item{Output}{
if \code{Output = TRUE} write a \code{Output} to an Excel 2007.
}
}
\details{
The attractive models is defined by the system for stochastic differential equation Two-dimensional : 
\emph{\deqn{dX(t) = (-K * X(t)/(sqrt(X(t)^2 + Y(t)^2))^(S+1) )* dt + Sigma* dW1(t)}}
\emph{\deqn{dY(t) = (-K * Y(t)/(sqrt(X(t)^2 + Y(t)^2))^(S+1) )* dt + Sigma* dW2(t)}}
\code{dW1(t)} and \code{dW2(t)} are brownian motions independent. 

If \code{S = 1} (ie \code{M(S=1,Sigma)}) the system SDE is :
\emph{\deqn{dX(t) = (-K * X(t)/(X(t)^2 + Y(t)^2) )* dt + Sigma* dW1(t)}}
\emph{\deqn{dY(t) = (-K * Y(t)/(X(t)^2 + Y(t)^2) )* dt + Sigma* dW2(t)}}

For more detail consulted \code{References}.
}
\value{
data.frame(time,X(t),Y(t)) and plot of process 2-D.
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
\code{\link{snssde2D}}, \code{\link{PredCorr2D}}, \code{\link{RadialP2D_1PC}}, \code{\link{RadialP3D_1}}, \code{\link{tho_M1}}, \code{\link{fctgeneral}}, \code{\link{hist_general}}, \code{\link{Kern_meth}}.
}
\examples{
RadialP2D_1(N=1000, t0=0, Dt=0.001, T = 1, X0=2, Y0=1, v=0.3, 
            K=3, Sigma=0.2, Output = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Attractive Model}
\keyword{ Diffusion Process Multidimensional}
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Stochastic Differential Equation Multidimensional}
