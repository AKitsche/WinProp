#print functions for objects WinPropRaw and 
print.winprop <- function(x, ...){
  cat("Results from the two-sample t-Test: \n")
  cat("   x-y     = ",x[[1]], "\n")
  cat("   t       = ",x[[2]], "\n")
  cat("   df      =  ",x[[4]], "\n")
  cat("   p-value =  ",x[[3]], "\n  \n")
  cat("95 percent confidence interval: \n")
  cat("Conf.Int.: ",x[[5]][[1]], x[[5]][[2]], "\n \n")
  cat("Win Probability (W=P[X>Y]) and confidence interval: \n")
  cat("   W     = ",x[[7]][[1]], "\n")
  cat("Conf.Int.: ",x[[7]][[2]],x[[7]][[3]], "\n \n")
  cat("Odds of X being greater than Y and confidence interval: \n")
  cat("W/(1-W)  = ",x[[9]][[1]], "\n")
  cat("Conf.Int.: ",x[[9]][[2]],x[[9]][[3]], "\n \n")
  cat("Confidence interval for c_[beta]: \n")
  cat("Conf.Int.: ",x[[8]][[1]], x[[8]][[2]], "\n \n")
}
