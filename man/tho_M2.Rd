\name{tho_M2}
\alias{tho_M2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation The First Passage Time FPT For Attractive Model(S >= 2,Sigma)
}
\description{
simulation M-sample for the first passage time "FPT" for attractive model(S >= 2,Sigma).
}
\usage{
tho_M2(N, M, t0, T, R0, v, K, s, Sigma, Output = FALSE,
       Methods = c("Euler", "Milstein", "MilsteinS", 
                   "Ito-Taylor", "Heun", "RK3"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of the diffusion process.
}
  \item{M}{
size of the FPT.
}
  \item{t0}{
initial time.
}
  \item{T}{
final time.
}
  \item{R0}{
initial value of the process at time \code{t0} ,(R0 > 0).
}
  \item{v}{
threshold. \code{see detail}.
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
if \code{Output = TRUE} write a \code{Output} to an Excel (.csv).
}
  \item{Methods}{
method of simulation ,see details \code{\link{snssde}}.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
Using Ito transform, it is shown that the \code{Radial Process R(t)} with \code{R(t)=||(X(t),Y(t))||} is a markovian diffusion, solution of the stochastic differential equation one-dimensional:
\emph{\deqn{dR(t) = ((0.5 * Sigma^2 * R(t)^(S-1) - K)/ R(t)^S )* dt + Sigma* dW(t)}}

We take interest in the random variable \code{FPT} "first passage time", is defined by : \emph{\deqn{FPT = inf(t>=0 \ R(t) <= v )}} with \code{v is the threshold}.

For more detail consulted \code{References}.
}
\value{
M-sample for FPTT.
}
\references{
\enumerate{\item K.Boukhetala, Estimation of the first passage time distribution for a simulated diffusion process, Maghreb Math.Rev, Vol.7, No 1, Jun 1998, pp. 1-25. \item K.Boukhetala, Simulation study of a dispersion about an attractive centre. In proceedings of 11th Symposium Computational Statistics, edited by R.Dutter and W.Grossman, Wien , Austria, 1994, pp. 128-130. \item K.Boukhetala,Modelling and simulation of a dispersion pollutant with attractive centre, Edited by Computational Mechanics Publications, Southampton ,U.K and Computational Mechanics Inc, Boston, USA, pp. 245-252.\item K.Boukhetala, Kernel density of the exit time in a simulated diffusion, les Annales Maghrebines De L ingenieur, Vol , 12, N Hors Serie. Novembre 1998, Tome II, pp 587-589. 
}
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\note{
\itemize{\item \code{2*K > Sigma^2}.
}

o system.time(tho_M2(N=1000, M=100, t0=0, T=1, R0=2, v=0.05, K=3,
                   s=2,Sigma=0.3,Output = FALSE,Methods="Euler"))

utilisateur     systeme      ecoule 

       9.58        0.14        9.74

o system.time(tho_M2(N=1000, M=100, t0=0, T=1, R0=2, v=0.05, K=3,
                  s=2,Sigma=0.3,Output = FALSE,Methods="RK3"))

utilisateur     systeme      ecoule 

      51.29        0.36       52.79
}
\seealso{
\code{\link{AnaSimFPT}} Simulation The First Passage Time FPT For A Simulated Diffusion Process.
}
\examples{
\donttest{
 tho_M2(N=1000, M=50, t0=0, T=1, R0=2, v=0.05, K=3,s=2,
        Sigma=0.3,Output = FALSE,Methods="Euler")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Attractive models }
\keyword{ Stochastic Differential Equation }
\keyword{ Statistical Analysis }
