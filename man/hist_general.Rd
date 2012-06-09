\name{hist_general}
\alias{hist_general}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adjustment the Density of Random Variable X by Histograms Methods  
}
\description{
Adjusted your density of random variable X by histograms methods with Different number of cells.
}
\usage{
hist_general(Data, Breaks, Law = c("exp", "GAmma", "chisq", "Beta", 
             "fisher", "student", "weibull", "Normlog", "Norm"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
a numeric vector of the observed values.
}
  \item{Breaks}{
one of:
 o a vector giving the breakpoints between histogram cells.
 o a single number giving the number of cells for the histogram.
 o a function to compute the number of cells.
 o Breaks = c('scott','Sturges','FD') or manual.
}
  \item{Law}{
distribution function with Adjusted. see details \code{Distributions} (R >= 2.12.1)
}
}
\details{
Ajusted the density for random variable X by histograms methods with Different number of cells see details \code{nclass.scott}, 
ajusted with the Distribution c("dexp","dgamma",
"dchisq","dbeta","df","dt","dweibull",
"dlnorm","dnorm").
}
\value{
plot.histogram with Adjustment and Estimation.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\note{
Choose your best distribution with minimum AIC.
}
\seealso{
\code{\link{fctgeneral}} empirical distribution, \code{\link{Kern_general}} Kernel Methods.
}
\examples{
 
 X <- rexp(1000,2)
 par(mfrow=c(2,2))
 hist_general(Data=X, Breaks='FD', Law="exp")
 hist_general(Data=X, Breaks='scott', Law="exp")
 hist_general(Data=X, Breaks='Sturges', Law="exp")
 hist_general(Data=X, Breaks=60, Law="exp")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Parametric Estimation }
\keyword{ Statistical Analysis }