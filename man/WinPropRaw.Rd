\name{WinPropRaw}
\alias{WinPropRaw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confidence Interval for Win Probability using Raw data
}
\description{
This function calculates several effect size measurements for the comparison of two treatments:

- t-Test Statistic and the corresponding p-value

- Confidence interval for the difference of the two treatments

- Confidence interval for the Win Probability 

- Confidence interval for the threshold c, that specifies a threshold to which a potential future   
   observation 

- Odds W/(1-W) and corresponding confidence interval
}
\usage{
WinPropRaw(x, y, alpha=0.05, beta=0.95, c=0, var.equal=TRUE, 
alternative=c("two.sided", "greater","less"), m1=1, m2=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a numeric vector of data values from the first group
}
  \item{Y}{
a numeric vector of data values from the second group
}
\item{alpha}{
a number indicating the type-I error rate
}
  \item{beta}{
a number indicating the chance that the first group is different from the second group
}
  \item{c}{
a number indicating the value c of interest, WinProp(c) directly assesses the probability that the difference between potential future single observations will be at least that value
}
  \item{var.equal}{
a logical variable indicating whether to treat the two variances as being equal. 
}
\item{alternative}{
a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
}
  \item{m1}{
a number indicating the number of potential future observations from treatment 1
}
  \item{m2}{
a number indicating the number of potential future observations from treatment 2
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with class "winprop" containing the following components:
\item{pvalue}{p-value associated with the two-sample t-test}
\item{n1}{sample size of the first group}
\item{n1}{sample size of the second group}
\item{delta}{observed Cohens effect size (mean(x)-mean(y))/s and its corresponding confidence interval}
\item{lambda}{observed non-centrality parameter and its corresponding confidence interval}
\item{W}{relative effect P(X>Y) and its corresponding confidence intervals}
\item{Psi}{Odds of X being greater than Y W/(1-W) and its corresponding confidence interval}
\item{Diff}{Difference in means}
\item{TStat}{the value of the t-statistic}
\item{pvalue}{the p-value for the t-test}
\item{df}{the degrees of freedom for the t-statistic}
\item{CI}{a confidence interval for the mean appropriate to the specified alternative hypothesis}
\item{PI}{a predicition interval for the difference in means}
\item{W}{Win Probability and corresponding confidence interval}
\item{Cbeta}{confidence interval for the threshold c}
\item{Phi}{odds W/(1-W) and corresponding confidence interval}
}
\references{
A. J. Hayter (2013): Inferences on the difference between future observations for comparing two treatments, Journal of Applied Statistics, 40:4, 887-900
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
#The following example is adopted from Hayter (2013).

#########################################################################################
#############                Example 5           ########################################
#########################################################################################
X <- c(76, 57, 71, 57, 65, 64, 65, 64, 70, 59)
Y <- c(52, 53, 40, 58, 46)


WinPropRaw(X=X, Y=Y, alpha=0.05, c=0, var.equal=FALSE, alternative="greater", m1=1, m2=1)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
