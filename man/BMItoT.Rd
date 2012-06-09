\name{BMItoT}
\alias{BMItoT}
\alias{ItovsStraT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Properties of the stochastic integral and Ito Process [5]
}
\description{
Simulation of the Ito \code{integral(s*dW(s),0,t)}. 
}
\usage{
BMItoT(N, T, output = FALSE)
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
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\details{
However the Ito integral also has the peculiar property, amongst others, that : \emph{\deqn{integral(s*dW(s),0,t) =  t * W(t) - integral(W(s)*ds,0,t)}}
from classical calculus for Ito integral with \code{w(0) = 0}.

The  follows from the algebraic rearrangement : \emph{\deqn{integral(s*dW(s),0,t) = sum ( t *(W(t+1)-W(t)),0,t)}}
}
\value{
data frame(time,Ito,sum.Ito) and plot of the Ito integral.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{BMIto1}} simulation of the Ito integral[1], \code{\link{BMIto2}} simulation of the Ito integral[2], \code{\link{BMItoC}} properties of the stochastic integral and Ito processes[3], \code{\link{BMItoP}} properties of the stochastic integral and Ito processes[4].}
\examples{
\donttest{
 BMItoT(N=1000,T=1)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Stochastic integral }
