\name{AnaSimX}
\alias{AnaSimX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation M-Samples of Random Variable X(v[t]) For A Simulated Diffusion Process
}
\description{
Simulation M-samples of the random variable X(v(t)) by a simulated diffusion process with a fixed the time v , \code{ v = k * Dt} with \code{k integer, 1 <= k <= N }.
}
\usage{
AnaSimX(N, M, t0, Dt, T = 1, X0, v, drift, diff, Output = FALSE, 
        Methods = c("Euler", "Milstein", "MilsteinS", "Ito-Taylor", 
                    "Heun", "RK3"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of the diffusion process.
}
  \item{M}{
size of the random variable.
}
  \item{t0}{
initial time.
}
  \item{Dt}{
time step of the simulation \code{(discretization)}.
}
  \item{T}{
final time.
}
  \item{X0}{
initial value of the process at time \code{t0}.
}
  \item{v}{
moment (time) between \code{t0} and \code{T} ,\code{ v = k * Dt} with \code{k integer, 1 <= k <= N }.
}
  \item{drift}{
drift coefficient: an expression of two variables \code{t} and \code{x}.
}
  \item{diff}{
diffusion coefficient: an expression of two variables \code{t} and \code{x}.
}
  \item{Output}{
if \code{Output = TRUE} write a \code{Output} to an Excel (.csv).
}
  \item{Methods}{
method of simulation ,see details \code{\link{snssde}}.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
The stochastic differential equation of  is  : \emph{\deqn{dX(t) = a(t,X(t)) *dt + b(t,X(t)) *dW(t)}} with \code{a(t,X(t)) :drift coefficient} and \code{b(t,X(t)) :diffusion coefficient}, \code{W(t)} is Wiener process.

We take interest in the random variable \code{X(v)}, is defined by : \emph{\deqn{X =(t>=0 \ X = X(v))}} with \code{v} is the time between \code{t0} and \code{T} ,\code{ v = k * Dt} with \code{k integer, 1 <= k <= N }.
}
\value{
Random variable "X(v(t))".
}
\references{
\enumerate{\item K.Boukhetala, Estimation of the first passage time distribution for a simulated diffusion process, Maghreb Math.Rev, Vol.7, No 1, Jun 1998, pp. 1-25. \item K.Boukhetala, Simulation study of a dispersion about an attractive centre. In proceedings of 11th Symposium Computational Statistics, edited by R.Dutter and W.Grossman, Wien , Austria, 1994, pp. 128-130. \item K.Boukhetala,Modelling and simulation of a dispersion pollutant with attractive centre, Edited by Computational Mechanics Publications, Southampton ,U.K and Computational Mechanics Inc, Boston, USA, pp. 245-252.\item K.Boukhetala, Kernel density of the exit time in a simulated diffusion, les Annales Maghrebines De L ingenieur, Vol , 12, N Hors Serie. Novembre 1998, Tome II, pp 587-589. 
}
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\note{
Time of Calculating

The Ornstein-Uhlenbeck Process (example)
drift <- expression(-5*x)
diff <- expression(1)

system.time(AnaSimX(N=1000,M=30,t0=0,Dt=0.001,T=1,X0=0,
            v=0.5,drift,diff,Output=FALSE))

utilisateur     systeme      ecoule 

       1.88        0.56        2.59

system.time(AnaSimX(N=1000,M=30,t0=0,Dt=0.001,T=1,X0=0,
            v=0.5,drift,diff,Output=FALSE,Methods="RK3"))

utilisateur     systeme      ecoule 

       8.64        0.72        9.24
}
\seealso{
\code{\link{AnaSimFPT}} Simulation The First Passage Time FPT For A Simulated Diffusion Process, \code{\link{tho_M1}} Simulation The FPT For Attractive Model(S = 1,Sigma), \code{\link{tho_M1}} Simulation The FPT For Attractive Model(S >= 2,Sigma), \code{\link{tho_02diff}} Simulation FPT For Attractive Model for 2-Diffusion Processes.
}
\examples{
\donttest{
## Example 1:  BM
## v = k * Dt with k integer , 1 <= k <= N .
## k = 500 nombre for discretization
## Dt = 0.001 ===> v = 500 * 0.001 = 0.5

 drift <- expression(0)
 diff <- expression(1)
 AnaSimX(N=1000,M=30,t0=0,Dt=0.001,T=1,X0=0,v=0.5,drift,diff,Output=FALSE,Methods="Euler")
 summary(X)
 hist(X)
 v=0.5
 plot(density(X,kernel ="gaussian"),col="red")
 x <- seq(min(X),max(X),length=1000)
 curve(dnorm(x,0,v), col = 3, lwd = 2, add = TRUE,
      panel.first=grid(col="gray"))


## Example 2: BMG or BS
## v = k * Dt with k integer , 1 <= k <= N .
## k = 800 nombre for discretization
## Dt = 0.001 ===> v = 800 * 0.001 = 0.8

 drift <- expression(2*x)
 diff <- expression(x)
 AnaSimX(N=1000,M=30,t0=0,Dt=0.001,T=1,X0=1,v=0.8,drift,diff,Output=FALSE,Methods="Euler")
 summary(X)
 hist(X)
 plot(density(X,kernel ="gaussian"),col="red")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Stochastic Differential Equation }
\keyword{ Statistical Analysis }
