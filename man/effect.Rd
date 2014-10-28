\name{effect}
\alias{effect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to convert a p-value and associate sample size into alternative efffect size measurements
}
\description{
This function provides more meaningful measurement of the effect size associated with the p-value from the two-sample t-test. The focus is on the calculation of the relative effect defined by the probability that an outcome under the experimental treatment is better than an outcome under the control P(X>Y). For each effect size measure the corresponding 95\% confidence intervals are given

- Cohens effect size (mean(X)-mean(Y))/s

- relative effect P(X>Y) (aka stochastic superiority, probabilistic index, Win probability)

- odds P(X>Y)/(1-P(X>Y))
}
\usage{
effect(pvalue, n1, n2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pvalue}{
p-value from the two-sample t-test
}
  \item{n1}{
sample size of the first group
}
  \item{n2}{
sample size of the second group
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{pvalue}{p-value associated with the two-sample t-test}
\item{n1}{sample size of the first group}
\item{n1}{sample size of the second group}
\item{delta}{observed Cohens effect size (mean(x)-mean(y))/s and its corresponding confidence interval}
\item{lambda}{observed non-centrality parameter and its corresponding confidence interval}
\item{W}{relative effect P(X>Y) and its corresponding confidence intervals}
\item{Psi}{Odds of X being greater than Y W/(1-W) and its corresponding confidence interval}
%% ...
}
\references{
A. J. Hayter (2013): Inferences on the difference between future observations for comparing two treatments, 
Journal of Applied Statistics, 40:4, 887-900

Richard H. Browne (2010) The t-Test p Value and Its Relationship to the Effect Size and P(X>Y), 
The American Statistician, 64:1, 30-33
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
effect(pvalue=0.01, n1=10,n2=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
