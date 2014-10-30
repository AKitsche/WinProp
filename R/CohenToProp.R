#Function to convert Cohens effect size into Probability that X>Y and the corresponding odds ratio

#arguments
#delta - Cohens effect size

CohenToProp <- function(delta){
  #checks
  if(!is.numeric(delta)) {
    stop("delta must be a numeric value")
  }
  
  #Calculating W=P(X>Y)
  W <- pnorm(q=delta/sqrt(2))
  #calculating the odds ratio
  Phi <- W/(1-W)
  
  #output
  return(list(W=W,
              Phi=Phi))
}