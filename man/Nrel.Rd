\name{Nrel}
\alias{Nrel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate sample size per group required to observe a relative effect of at least phi
}
\description{
The function calculates the sample size per group required to observe  relative effect of at least phi with probability 1-beta, if te approximate sample size formula for the t-test is used. The sam,ple size is caluculated by iteratively searching for nrel according to Equation (2) in Kieser et al. (2013).
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

%%  ~Describe the value returned
%%  If it is a LIST, use
\item{nrel}{per group sample size}
\item{phi}{relative effect phi}
\item{sigma}{standard error}
\item{delta}{difference in expectation}
\item{beta}{type II error rate}
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
