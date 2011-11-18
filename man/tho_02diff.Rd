\name{tho_02diff}
\alias{tho_02diff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation The First Passage Time FPT For Attractive Model for Two-Diffusion Processes V(1) and V(2)
}
\description{
simulation M-sample for the first passage time "FPT" for attractive for 2-diffusion processes V(1)=c(X1(t),X2(t)) and V(2)=c(Y1(t),Y2(t)) or V(1)=c(X1(t),X2(t),X3(t)) and V(2)=c(Y1(t),Y2(t),Y3(t)). 
}
\usage{
tho_02diff(N, M, t0, Dt, T = 1, X1_0, X2_0, Y1_0, Y2_0,
           v, K, m, Sigma,Output=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of the diffusion process V1(t) and V2(t).
}
  \item{M}{
size of the FPT.
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
  \item{Y1_0}{
initial value of the process \code{Y1(t)} at time \code{t0}.
}
  \item{Y2_0}{
initial value of the process \code{Y2(t)} at time \code{t0}.
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
if \code{Output = TRUE} write a \code{Output} to an Excel 2007.
}
}
\details{
The 2-dimensional attractive models for 2-diffusion processes V(1)=(X1(t),X2(t)) and V(2)=c(Y1(t),Y2(t)) is defined by the Two (02) system for stochastic differential equation Two-dimensional :
 
\emph{\deqn{dV1(t) = dV2(t)+Mu(m+1)(||D(t)||)*D(t)*dt +Sigma I(2*2)*dW1(t)}}
\emph{\deqn{dV2(t) = Sigma*I(2*2)*dW2(t)}}

with:

\emph{\deqn{ D(t) = V1(t) - V2(t)}}
\emph{\deqn{ Mu(m)(||d||) = -K/||d||^m}}

Where ||.|| is the Euclidean norm and I(2*2) is identity matrix, \code{dW1(t)} and \code{dW2(t)} are brownian motions independent. 

\code{D(t)=sqrt((X1(t)^2 - Y1(t)^2)+(X2(t)^2 - Y2(t)^2))} it is distance between \code{V1(t)} and \code{V2(t)} .

And the random variable \code{tau} "first passage time FPT", is defined by : \emph{\deqn{tau(V1(t),V2(t)) = inf(t>=0 \ ||D(t)|| <= v )}} with \code{v is the threshold}.
}
\value{
Random variable tau "FPT".
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\note{
\itemize{\item \code{2*K > Sigma^2}.
}
}
\seealso{
\code{\link{TwoDiffAtra3D}}, \code{\link{TwoDiffAtra2D}}, \code{\link{fctgeneral}}, \code{\link{hist_general}}, \code{\link{Kern_meth}},\code{\link{AnaSimFPT}} Simulation The First Passage Time FPT For A Simulated Diffusion Process.
}
\examples{
tho_02diff(N=500, M=50, t0=0, Dt=0.001, T = 1,
           X1_0=1, X2_0=1,Y1_0=0.5, Y2_0=0.5,v=0.05, 
           K=4, m=0.2, Sigma=0.2)
summary(FPT)
hist(FPT,prob=TRUE)
plot(density(FPT,kernel="gaussian"),col="red")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Attractive models }
\keyword{ Stochastic Differential Equation }
\keyword{ Statistical Analysis }
