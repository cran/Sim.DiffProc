\name{BB}
\alias{BB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Brownian Bridge Model
}
\description{
Simulation of brownian bridge model.
}
\usage{
BB(N, t0, T, x0, y, output = FALSE)
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
initial value of the process at time \code{t0}.
}
  \item{y}{
terminal value of the process at time \code{T}.
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\details{
The function returns a trajectory of the brownian bridge starting at \code{x0} at time \code{t0} and ending at \code{y} at time \code{T}.

It is defined as : \emph{\deqn{Xt(t0,x0,T,y) = x0 + W(t-t0) - (t-t0/T-t0)*(W(T-t0)-y+x0)}}
This process is easily simulated using the simulated trajectory of the Wiener process \code{W(t)}.
}
\value{
data.frame(time,x) and plot of process. 
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{BBF}} simulation flow of brownian bridge Model, \code{\link{diffBridge}} Diffusion Bridge Models, \code{\link{BMN}} simulation brownian motion by the Normal Distribution , \code{\link{BMRW}} simulation brownian motion by a Random Walk, \code{\link{GBM}} simulation geometric brownian motion, \code{\link{ABM}} simulation arithmetic brownian motion, \code{\link{snssde}} Simulation Numerical Solution of SDE. 
}
\examples{
\donttest{
##brownian bridge model
##starting at x0 =0 at time t0=0 and ending at y =3 at time T =1.
BB(N=1000,t0=0,T=1,x0=0,y=3)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
\keyword{ Stochastic Differential Equation }


