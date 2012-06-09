\name{TwoDiffAtra3D}
\alias{TwoDiffAtra3D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Three-Dimensional Attractive Model for Two-Diffusion Processes V(1) and V(2)
}
\description{
simulation 3-dimensional attractive model for 2-diffusion processes V(1)=(X1(t),X2(t),X3(t)) and V(2)=c(Y1(t),Y2(t),Y3(t)).
}
\usage{
TwoDiffAtra3D(N, t0, Dt, T = 1, X1_0, X2_0, X3_0, Y1_0,
              Y2_0, Y3_0, v, K, m, Sigma, Output = FALSE)
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
  \item{X1_0}{
initial value of the process \code{X1(t)} at time \code{t0}.
}
  \item{X2_0}{
initial value of the process \code{X2(t)} at time \code{t0}.
}
  \item{X3_0}{
initial value of the process \code{X3(t)} at time \code{t0}.
}
  \item{Y1_0}{
initial value of the process \code{Y1(t)} at time \code{t0}.
}
  \item{Y2_0}{
initial value of the process \code{Y2(t)} at time \code{t0}.
}
  \item{Y3_0}{
initial value of the process \code{Y3(t)} at time \code{t0}.
}
  \item{v}{
threshold. \code{see detail}
}
  \item{K}{
constant \code{K > 0}.
}
  \item{m}{
constant \code{m > 0}.
}
  \item{Sigma}{
constant \code{Sigma > 0}.
}
  \item{Output}{
if \code{Output = TRUE} write a \code{Output} to an Excel (.csv).
}
}
\details{
The 3-dimensional attractive models for 2-diffusion processes V(1)=(X1(t),X2(t),X3(t)) and V(2)=c(Y1(t),Y2(t),Y3(t)) is defined by the Two (02) system for stochastic differential equation three-dimensional :
 
\emph{\deqn{dV1(t) = dV2(t)+Mu(m+1)(||D(t)||)*D(t)*dt +Sigma I(3*3)*dW1(t)}}
\emph{\deqn{dV2(t) = Sigma*I(3*3)*dW2(t)}}

with:

\emph{\deqn{ D(t) = V1(t) - V2(t)}}
\emph{\deqn{ Mu(m)(||d||) = -K/||d||^m}}

Where ||.|| is the Euclidean norm and I(3*3) is identity matrix, \code{dW1(t)} and \code{dW2(t)} are brownian motions independent. 

\code{D(t)=sqrt((X1(t)^2 - Y1(t)^2)+(X2(t)^2 - Y2(t)^2)+(X3(t)^2-Y3(t)^2))} it is distance between \code{V1(t)} and \code{V2(t)} .

And the random variable \code{tau} "first passage time", is defined by : \emph{\deqn{tau(V1(t),V2(t)) = inf(t>=0 \ ||D(t)|| <= v )}} with \code{v is the threshold}.
}
\value{
data.frame(time,X1(t),X2(t),X3(t),Y1(t),Y2(t),Y3(t),D(t)) and plot of process 3-D.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\note{
\itemize{\item \code{2*K > Sigma^2}.
}
}
\seealso{
\code{\link{TwoDiffAtra2D}}, \code{\link{tho_02diff}}.
}
\examples{
\donttest{
 TwoDiffAtra3D(N=500, t0=0, Dt=0.001, T = 1, X1_0=0.5, X2_0=0.25,
               X3_0=0.1,Y1_0=-0.5,Y2_0=-1, Y3_0=0.25, v=0.01, K=5,
               m=0.2, Sigma=0.1, Output = FALSE)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Attractive Model}
\keyword{ Diffusion Process Multidimensional}
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Stochastic Differential Equation Multidimensional}
