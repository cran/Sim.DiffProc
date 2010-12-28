\name{BMItoC}
\alias{BMItoC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Properties of the stochastic integral and Ito Process [3]
}
\description{
Simulation of the Ito \code{integral(alpha*dW(s),0,t)}. 
}
\usage{
BMItoC(N, T, alpha, output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of process.
}
  \item{T}{
final time.
}
  \item{alpha}{
constant.
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel 2007.
}
}
\details{
However the Ito integral also has the peculiar property, amongst others, that : \emph{\deqn{integral(alpha*dW(s),0,t) = alpha * W(t)}}
from classical calculus for Ito integral with \code{w(0) = 0}.

The  follows from the algebraic rearrangement : \emph{\deqn{integral(alpha * dW(s),0,t) = sum ( alpha*(W(t+1)-W(t)),0,t)}} 
}
\value{
data frame(time,Ito,sum.Ito) and plot of the Ito integral.
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\seealso{
\code{\link{BMIto1}} simulation of the Ito integral[1], \code{\link{BMIto2}} simulation of the Ito integral[2], \code{\link{BMItoP}} properties of the stochastic integral and Ito processes[4], \code{\link{BMItoT}} properties of the stochastic integral and Ito processes[5].
}
\examples{
##
BMItoC(N=1000,T=1,alpha=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Stochastic integral }
