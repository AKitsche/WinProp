#print functions for objects WinPropRaw and 
print.SampleSize <- function(x, ...){
  cat("required sample size per group for demonstrating statistical significance with alpha =",x[[6]] ,"\n")
  cat("   n     = ",x[[2]], "\n")
  cat("   power       = ",1-x[[7]], "\n \n")
  cat("required sample size per group to observe a relative effect of at least phi P(X > Y)=",x[[3]] ,"\n")
  cat("   n     = ",x[[1]], "\n")
  cat("   power       = ",1-x[[7]], "\n \n")
}