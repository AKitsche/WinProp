#Function to calculate effect sizes and the corresponding confidence intervals 
#for two sample comparisons

#This is a wrapper function that uses the function WinPropSum() after calculating summary statistics

#Function arguments:
#x - vector of observations from the first group
#y - vector of observations from the second group
#alpha - type-I error rate
#beta - chance that the first group is different from the second group
#c - value c of interest, WinProp(c) directly assesses the probability that the difference between potential future single observations will be at least that value
#var.equal - specifies if homogeneous variances between the two groups is assumed
#alternative - specifies the direction of the hypothesis
#m1 - number of potential future observations from treatment 1
#m2 - number of potential future observations from treatment 2

WinPropRaw <- function(x, y, alpha=0.05, beta=0.95, c=0, var.equal=TRUE, alternative=c("two.sided", "greater","less"), m1=1, m2=1){
  options(warn=-1)
  #checks
  if(length(x) == 1 | length(y) == 1 | !is.numeric(x) | !is.numeric(y)){
    stop("x and y must be numeric vectors with length greater than 1")
  }
    if(alternative !=  "two.sided" & alternative !=  "greater" & alternative != "less") {
    stop("alternative must be either two.sided, greater or less")
  }
  if(length(c) != 1 | !is.numeric(c)) {
    stop("c must be a single numeric value")
  }
  if(alpha < 0 || alpha > 1) {
    stop("alpha must be a numeric value between 0 and 1")
  }
  if(beta < 0 || beta > 1) {
    stop("alpha must be a numeric value between 0 and 1")
  }
  #function to check for whole numbers
  is.whole <- function(x){
    is.numeric(x) && floor(x)==x 
  }
  if(length(m1) != 1 | length(m2) != 1 | !is.whole(m1) | !is.whole(m2)){
    stop("m1 and m2 must be positive whole numbers")
  }
  
  #calcuating the sample sizes
  n1 <- length(x)
  n2 <- length(y)
  #calculating the expected values
  X <- mean(x)
  Y <- mean(y)
  #calcuating the standard deviation
  if(var.equal==TRUE){
    df <- n1+n2-2#degrees of freedom
    s <- sqrt(((var(x)*(n1-1))+ (var(y)*(n2-1)))/(n1+n2-2))
    s2 <- s
  }else{
    s <- sd(x)
    s2 <- sd(y)
    df <- ((s^2/n1)+(s2^2/n2))^2/(s^4/(n1^2*(n1-1))+s2^4/(n2^2*(n2-1)))
  }

  
  #calculating Win probability and corresponding confidence interval and confidence interval for c_[beta]
  outWinPropSum <- WinPropSum(X=X, Y=Y, s=s, n1=n1, n2=n2, alpha=alpha, beta=beta, c=c, var.equal=var.equal, s2=s2, alternative=alternative, m1=m1, m2=m2)

  #output
  out <- list(Diff=outWinPropSum$Diff,
              TStat=outWinPropSum$TStat,
              pvalue=outWinPropSum$pvalue,
              df=outWinPropSum$df,
              CI = outWinPropSum[[5]],
              PI = outWinPropSum[[6]],
              W =  outWinPropSum[[7]],
              Cbeta = outWinPropSum[[8]],
               Phi = outWinPropSum[[9]])
  class(out) <- "winprop"
  out
}