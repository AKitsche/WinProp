\name{Nrel}
\alias{Nrel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate sample size per group required to observe a relative effect of at least phi
}
\description{
The function calculates the sample size per group required to observe a relative effect of at least phi with probability 1-beta, if the approximate sample size formula for the t-test is used. The sample size is caluculated by iteratively searching for n according to Equation (2) in Kieser et al. (2013). The relative effect phi defines a population measure on the original scale but also includes a success criterion for the individual patients.
}
\usage{
Nrel(nrel_start, phi, sigma, delta, beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nrel_start}{
starting sample size for iterative search
}
  \item{phi}{
pre-specified threshold value
}
  \item{sigma}{
standard deviation
}
  \item{delta}{
difference in expectation
}
  \item{beta}{
type-II error rate
}
  \item{alpha}{
  type-I error rate for sample size  calculation to achieve significance
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

%%  ~Describe the value returned
%%  If it is a LIST, use
\item{nrel}{per group sample size to claim clinical relevance}
\item{nsig}{per group sample size to claim significant effect}
\item{phi}{prespecified minimum relative effect phi}
\item{sigma}{standard deviation}
\item{delta}{difference in expectation}
\item{theta}{estimated relative effect, theta= Phi((delta/sigma)/sqrt(2)), where Phi is the standard normal distribution}
\item{beta}{type-II error rate}
\item{power}{1-beta}
}
\references{
Kieser, M.,Friede, T., Gondan, M. Assessment of statistical significance and clinical relevance. Statistics in Medicine (2013) 32(10): 1707-1719
}
\author{
Andreas Kitsche
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
Nrel(nrel_start=100, phi=0.7, sigma=10, delta=-10, beta=0.1, sigma=0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
