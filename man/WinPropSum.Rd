\name{WinPropSum}
\alias{WinPropSum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confidence Interval for Win Probability using Summary data
}
\description{
This function calculates the Win Probability and the corresponding confidence interval from some summary statistics (mean of the groups, standard deviation and sample size)
}
\usage{
WinPropSum(X, Y, s, n1, n2, alpha = 0.05, beta = 0.95, 
c = 0, var.equal = TRUE, s2 = 1, 
alternative = c("two.sided", "greater", "less"), m1 = 1, m2 = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
mean of the first group
}
  \item{Y}{
mean of the second group
}
  \item{s}{
pooled standard deviation, if var.equal=TRUE

standard deviation of the first group, if var.equal=FALSE
}
  \item{n1}{
sample size of the first group
}
  \item{n2}{
sample size of the second group
}
  \item{alpha}{
type-I error rate
}
  \item{beta}{
chance that the first group is different from the second group
}
  \item{c}{
value c of interest, WinProp(c) directly assesses the probability that the difference between potential future single observations will be at least that value
}
  \item{var.equal}{
specifies if homogeneous variances between the two groups is assumed
}
  \item{s2}{
standard deviation of the second group, if var.equal=FALSE
}
  \item{alternative}{
specifies the direction of the hypothesis
}
  \item{m1}{
number of potential future observations from treatment 1
}
  \item{m2}{
number of potential future observations from treatment 2
}
}
\details{
This function is used internally in the function \code{WinPropRaw} after the corresponding summary statistics are claculated.
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
\code{\link{WinPropSum}}
}
\examples{
#The following examples are adopted from Hayter (2013).


#########################################################################################
#############                Example 1           ########################################
#########################################################################################
#n1=n2=8
#S1=S2=10
#X-Y=8.5
WinProp(X=18.5, Y=10, s=10, n1=8, n2=8, alpha=0.05, beta=0.95, c=0, alternative="two.sided", m1=1, m2=1)

#########################################################################################
#############                Example 2           ########################################
#########################################################################################
#n1=n2=32
#X-Y=1.2
#s=10
WinProp(X=11.2, Y=10, s=10, n1=32, n2=32, alpha=0.05,c=-15, alternative="two.sided", m1=1, m2=1)

#########################################################################################
#############                Example 3           ########################################
#########################################################################################
#n1=n2=8
#X-Y=15
#s=10
WinProp(X=25, Y=10, s=10, n1=8, n2=8, alpha=0.05, c=-15, alternative="two.sided", m1=1, m2=1)


#########################################################################################
#############                Example 4           ########################################
#########################################################################################
#n1=n2=32
#X-Y=15
#s=10
WinProp(X=25, Y=10, s=10, n1=32, n2=32, alpha=0.05, c=-15, alternative="two.sided", m1=1, m2=1)


#########################################################################################
#############                Example 6           ########################################
#########################################################################################
#n1=n2=9
#X=3.9
#Y=-0.4
#s=5.26
WinProp(X=3.9, Y=-0.4, s=5.26, n1=9, n2=9, alpha=0.05, c=0, var.equal=TRUE, alternative="two.sided", m1=1, m2=1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
