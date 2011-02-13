\name{fctgeneral}
\alias{fctgeneral}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adjustment the Empirical Distribution of Random Variable X  
}
\description{
Adjusted your empirical distribution of Random Variable X.
}
\usage{
fctgeneral(Data, Law = c("exp", "GAmma", "chisq", "Beta", "fisher",
                      "student", "weibull", "Normlog", "Norm"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
a numeric vector of the observed values.
}
  \item{Law}{
distribution function with Adjusted. see details \code{Distributions} (R >= 2.12.1)
}
}
\details{
calculating the empirical distribution \code{F[i]= (1/n)*Sum(V[i])} with \code{V[i]= 1 if x[i] <= X else V[i] = 0}. 

And ajusted with the Distribution c("pexp","pgamma","pchisq",
"pbeta","pf","pt","pweibull","plnorm","pnorm") 
}
\value{
Plot the empirical distribution with Adjustment and Estimation.
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\note{
Choose your best distribution with minimum AIC.
}
\seealso{
\code{\link{hist_general}} Histograms Methods,\code{\link{Kern_general}} Kernel Methods.
}
\examples{
## Example 
## X <- rgamma(100,1,4)
## par(mfrow=c(2,2))
## fctgeneral(Data=X,Law=("exp"))
## fctgeneral(Data=X,Law=("GAmma"))
## fctgeneral(Data=X,Law=("weibull"))
## fctgeneral(Data=X,Law=("Normlog"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Parametric Estimation }
\keyword{ Statistical Analysis }
