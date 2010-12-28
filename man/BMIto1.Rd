\name{BMIto1}
\alias{BMIto1}
\alias{ItovsStra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Properties of the stochastic integral and Ito Process [1]
}
\description{
Simulation of the Ito \code{integral(W(s)dW(s),0,t)}. 
}
\usage{
BMIto1(N, T, output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of process.
}
  \item{T}{
final time.
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel 2007.
}
}
\details{
However the Ito integral also has the peculiar property, amongst others, that : \emph{\deqn{integral(W(s)dW(s),0,t) = 0.5 * (W(t)^2 - t )}}
from classical calculus for Ito integral with \code{w(0) = 0}.

The  follows from the algebraic rearrangement : \emph{\deqn{integral(W(s)dW(s),0,t) = sum ( W(t)*(W(t+1)-W(t)),0,t)}}
}
\value{
data frame(time,Ito,sum.Ito) and plot of the Ito integral.
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\seealso{
\code{\link{BMIto2}} simulation of the Ito integral[2], \code{\link{BMItoC}} properties of the stochastic integral and Ito processes[3], \code{\link{BMItoP}} properties of the stochastic integral and Ito processes[4], \code{\link{BMItoT}} properties of the stochastic integral and Ito processes[5].
}
\examples{
##
BMIto1(N=1000,T=1)
## comparison with BMIto2
system.time(BMIto1(N=10^4,T=1))
system.time(BMIto2(N=10^4,T=1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Stochastic integral }
