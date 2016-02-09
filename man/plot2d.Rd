\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{plot2d}
\alias{plot2d}
\alias{lines2d}
\alias{points2d}
\alias{plot3D}
\alias{plot2d.default}
\alias{lines2d.default}
\alias{points2d.default}
\alias{plot3D.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting for Class SDE}
\description{Generic function for plotting.}
\usage{
\method{plot2d}{default}(x, \dots)
\method{lines2d}{default}(x, \dots)
\method{points2d}{default}(x, \dots)
\method{plot3D}{default}(x, display = c("persp","rgl"), \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object inheriting from class \code{\link{snssde2d}}, \code{\link{snssde3d}}, \code{\link{bridgesde2d}} and \code{\link{bridgesde3d}}.}
  \item{display}{\code{"persp"} perspective or \code{"rgl"} plots.}   
  \item{\dots}{other graphics parameters, see \code{\link[graphics]{par}} in package "graphics".}
}
\details{
The 2 and 3-dim plot of class sde.
}
\author{
A.C. Guidoum, K. Boukhetala.
}
\examples{
## Example 1: 
set.seed(1234)

fx <- expression(0)
gx <- expression(1)
fy <- expression(0)
gy <- expression(1)

res <- snssde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,N=10000)
plot2d(res,type="l")

## Example 2: 
set.seed(1234)

fx <- expression(0)
gx <- expression(1)
fy <- expression(0)
gy <- expression(1)
fz <- expression(0)
gz <- expression(1)

res <- snssde3d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,driftz=fz,diffz=gz,N=10000)
plot3D(res,display="persp")
}
