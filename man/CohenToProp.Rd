\name{CohenToProp}
\alias{CohenToProp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to convert Cohens effect size into Probability that X>Y and the corresponding odds ratio
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CohenToProp(delta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{
Cohens effect size (mean(x)-mean(y))/s
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
A. J. Hayter (2013): Inferences on the difference between future observations for comparing two treatments, 
Journal of Applied Statistics, 40:4, 887-900
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
CohenToProp(4.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
