#Function to calculate sample size per group required to observe a relative effect of at least phi
#According to Eq. 2 in Kiser et al. (2013) Statistics in Medicine 32 1707–1719

#arguments:
#nrel_start - starting sample size for iterative search
#phi - pre-specified threshold value
#sigma - standard deviation
#delta - difference in expectation
#beta - type-II error
Nrel <- function(nrel_start, phi, sigma, delta, beta){
  #checks
  if(phi <= 0.5) {
    stop("phi must be a single numeric value greater than 0.5")
  }
  if(sigma < 0) {
    stop("sigma must be a positive numeric value")
  }
  if(beta < 0 || beta > 1) {
    stop("beta must be a numeric value between 0 and 1")
  }
  if(nrel_start < 0 ) {
    stop("nrel_start must be a positive numeric value")
  }
  
  #iteratively searching the sample size until condition nrel_start=nrel is fulfilled
  #formula corresponds to Eq. 2 in Kieser et al. (2012)
  d <- 1
  nrel=nrel_start
  while(abs(d) > 0.001){
    nrel_start <- nrel
    alphaRel <- 1-pnorm(sqrt(nrel_start)*qnorm(phi))
    nrel <- 2 * (qnorm(1-alphaRel) + qnorm(1-beta))^2 * (sigma/delta)^2
    nrel
    #thetaHat <- pnorm(sqrt(1/n)*qt(1-Test$p.value,df=2*n-2))
    #  nrel2 <- (qnorm(1-alphaRel) + qnorm(1-beta))^2 * (1/(qnorm(thetaHat)^2))
    #  nrel
    d <- nrel - nrel_start
  }
  return(list(nrel=as.integer(round(nrel,0)),
              phi=phi, 
              sigma=sigma, 
              delta=delta, 
              beta=beta,
              power=1-beta))
}
