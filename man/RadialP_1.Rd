\name{RadialP_1}
\alias{RadialP_1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Radial Process Model(S = 1,Sigma) Or Attractive Model 
}
\description{
Simulation the radial process one-dimensional (S = 1).
}
\usage{
RadialP_1(N, t0, Dt, T = 1, R0, K, Sigma, Output = FALSE, 
          Methods = c("Euler", "Milstein", "MilsteinS", 
          "Ito-Taylor", "Heun", "RK3"), ...)
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
  \item{R0}{
initial value of the process at time \code{t0} ,(R0 > 0).
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
  \item{Methods}{
method of simulation ,see details \code{\link{snssde}}.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
The attractive models is defined by the system for stochastic differential equation two-dimensional : 
\emph{\deqn{dX(t) = (-K * X(t)/(sqrt(X(t)^2 + Y(t)^2))^(S+1) )* dt + Sigma* dW1(t)}}
\emph{\deqn{dY(t) = (-K * Y(t)/(sqrt(X(t)^2 + Y(t)^2))^(S+1) )* dt + Sigma* dW2(t)}}
\code{dW1(t)} and \code{dW2(t)} are brownian motions independent. 

Using Ito transform, it is shown that the \code{Radial Process R(t)} with \code{R(t)=||(X(t),Y(t))||} is a markovian diffusion, solution of the stochastic differential equation one-dimensional:
\emph{\deqn{dR(t) = ((0.5 * Sigma^2 * R(t)^(S-1) - K)/ R(t)^S )* dt + Sigma* dW(t)}}

If \code{S = 1} (ie \code{M(S=1,Sigma)}) the \code{R(t)} is :
\emph{\deqn{dR(t) = ((0.5*Sigma^2 -K )/R(t) ) * dt + Sigma* dW(t)}}
Where ||.|| is the Euclidean norm and \code{dW(t)} is a determined brownian motions.

For more detail consulted \code{References}.
}
\value{
data.frame(time,R(t)) and plot of process R(t).
}
\references{
\enumerate{\item K.Boukhetala, Estimation of the first passage time distribution for a simulated diffusion process, Maghreb Math.Rev, Vol.7, No 1, Jun 1998, pp. 1-25. \item K.Boukhetala, Simulation study of a dispersion about an attractive centre. In proceedings of 11th Symposium Computational Statistics, edited by R.Dutter and W.Grossman, Wien , Austria, 1994, pp. 128-130. \item K.Boukhetala,Modelling and simulation of a dispersion pollutant with attractive centre, Edited by Computational Mechanics Publications, Southampton ,U.K and Computational Mechanics Inc, Boston, USA, pp. 245-252.\item K.Boukhetala, Kernel density of the exit time in a simulated diffusion, les Annales Maghrebines De L ingenieur, Vol , 12, N Hors Serie. Novembre 1998, Tome II, pp 587-589. 
}
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\note{
\itemize{\item If \code{methods} is not specified, it is assumed to be the \code{Euler Scheme}. \item If \code{T} and \code{t0} specified, the best discretization \code{Dt = (T-t0)/N}. \item \code{2*K > Sigma^2}.
}
}
\seealso{
\code{\link{RadialP2D_1}}, \code{\link{RadialP2D_1PC}}, \code{\link{RadialP3D_1}}, \code{\link{tho_M1}}, \code{\link{fctgeneral}}, \code{\link{hist_general}}, \code{\link{Kern_meth}}.
}
\examples{
## Example 1
RadialP_1(N=1000, t0=0, Dt=0.001, T = 1, R0=2, K=2, 
          Sigma=0.5, Output = FALSE)
## Example 2
RadialP_1(N=1000, t0=0, Dt=0.001, T = 1, R0=3, K=3.5, 
          Sigma=0.5, Output = FALSE,Methods="Heun")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Attractive Model}
\keyword{ Diffusion Process Multidimensional}
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Stochastic Differential Equation Multidimensional}

