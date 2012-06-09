\name{hist_meth}
\alias{hist_meth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Histograms of Random Variable X  
}
\description{
The generic function hist_meth computes a histogram of the given data values.
}
\usage{
hist_meth(X, Breaks, Prob = c("TRUE", "FALSE"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a numeric vector of the observed values.
}
  \item{Breaks}{
one of:
 o a vector giving the breakpoints between histogram cells.
 o a single number giving the number of cells for the histogram.
 o a function to compute the number of cells.
 o Breaks = c('scott','Sturges','FD') or manual.

}
  \item{Prob}{
logical; if TRUE, the histogram graphic is a representation of frequencies, the counts component of the result; if FALSE, probability densities, component density, are plotted (so that the histogram has a total area of one). Defaults to TRUE if and only if breaks are equidistant (and probability is not specified).
}
}
\details{
The definition of histogram differs by source (with country-specific biases). R's default with equi-spaced breaks (also the default) is to plot the counts in the cells defined by breaks. Thus the height of a rectangle is proportional to the number of points falling into the cell, as is the area provided the breaks are equally-spaced. 
}
\value{
plot.histogram for the random variable X.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{Kern_meth}} Kernel Density,\code{\link{fctrep_Meth}} Empirical Distribution.
}
\examples{

 X <- rexp(1000,2)
 X11()
 hist_meth(X, Breaks='scott', Prob ="TRUE")
 curve(dexp(x, 2), col = 2, lwd = 2, add = TRUE)
 X11()
 hist_meth(X, Breaks='FD', Prob ="TRUE")
 curve(dgamma(x,1, 2), col = 2, lwd = 2, add = TRUE)
 X11()
 hist_meth(X, Breaks=100, Prob ="TRUE")
 curve(dweibull(x,1, 0.5),col=2, lwd = 2, add = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Parametric Estimation }
\keyword{ Statistical Analysis }
