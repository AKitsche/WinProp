\name{WinProp}
\alias{WinProp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confidence interval for Win Probability
}
\description{
This function calculates the confidence interval for a Win Probability and the confidence interval for the the threshold c 
}
\usage{
WinProp(X, Y, s, n1, n2, alpha = 0.05, beta = 0.95, 
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
hance that the first group is different from the second group
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
pecifies the direction of the hypothesis
}
  \item{m1}{
number of potential future observations from treatment 1
}
  \item{m2}{
number of potential future observations from treatment 2
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
List of
\item{Wl}{lower confidence limit of the Win Probability}
\item{Wu}{upper confidence limit of the Win Probability}
\item{Cbetal}{lower confidence limit of the threshold c}
\item{Cbetau}{upper confidence limit of the threshold c}
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
#############                Example 5           ########################################
#########################################################################################
X <- c(76, 57, 71, 57, 65, 64, 65, 64, 70, 59)
Y <- c(52, 53, 40, 58, 46)

t.test(X,Y, alternative="greater")
WinProp(X=mean(X), Y=mean(Y), s=sd(X), n1=length(X), n2=length(Y), alpha=0.05, c=0, s2=sd(Y), var.equal=FALSE, alternative="greater")


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
