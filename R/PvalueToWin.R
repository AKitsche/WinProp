#Function to convert a p-value and associated sample size into estimate of Cohens effect size

#Refernces
#A. J. Hayter (2013): Inferences on the difference between future observations for comparing two treatments, 
#Journal of Applied Statistics, 40:4, 887-900

#Richard H. Browne (2010) The t-Test p Value and Its Relationship to the Effect Size and P(X>Y), 
#The American Statistician, 64:1, 30-33

#Function arguments:
#pvalue - p-value from the two-sample t-test
#n1 - sample size of the first group
#n2 - sample size of the second group

PvalueToWin <- function(pvalue, n1, n2, alternative=c("two.sided", "greater","less"), alpha=0.05){
  #checks
  if(alternative !=  "two.sided" & alternative !=  "greater" & alternative != "less") {
    stop("alternative must be either two.sided, greater or less")
  }
  if(pvalue < 0 || pvalue > 1) {
    stop("pvalue must be a numeric value between 0 and 1")
  }
  if(alpha < 0 || alpha > 1) {
    stop("alpha must be a numeric value between 0 and 1")
  }
  #function to check for whole numbers
  is.whole <- function(x){
    is.numeric(x) && floor(x)==x 
  }
  #check for whole numbers of sample sizes  
  if(!is.whole(n1) | !is.whole(n2) | n1<0 | n2<0){
    stop("n1 and n2 must be positive whole numbers")
  }
  
  #calculating observed Cohens effect size 
  df <- n1+n2-2#degrees of freedom
  p <- pvalue#observed p-value
  
  alternative <- match.arg(alternative)
  switch(alternative, 
         two.sided = {
           tobs <- qt(p=1-(p/2), df=df, ncp=0)},#observed test statistic
         greater = {
           tobs <- qt(p=1-p, df=df, ncp=0)},#observed test statistic
         less = {
           tobs <- qt(p=1-p, df=df, ncp=0)})#observed test statistic
  
  delta <- tobs*(sqrt((1/n1)+(1/n2)))#observed Cohens effect size (mean(x)-mean(y))/s
  lambda <- delta/sqrt((1/n1)+(1/n2))
  #calculating exact 95% confidence intervals for lambda (non-centrality parameter)
  BoundLambda <- function(lambda=0, cl, df, tobs){#cl - confidence level
    sum(pt(q=tobs,ncp=lambda, df=df), -cl)
  }
  switch(alternative, 
         two.sided = {
           #lambdau <- uniroot(BoundLambda, c(-100,100), df=df, cl=alpha/2, tobs=tobs)$root#upper 95% confidence limit of lambda
           #lambdal <- uniroot(BoundLambda, c(-100,100), df=df, cl=1-alpha/2, tobs=tobs)$root#lower 95% confidence limit of lambda
           lambdau <- qt(p=1-alpha/2, df=df, ncp=tobs)#upper 95% confidence limit of lambda
           lambdal <- qt(p=alpha/2, df=df, ncp=tobs)#lower 95% confidence limit of lambda
           #calculating the confidence limits on delta
           deltal <- lambdal*sqrt((1/n1)+(1/n2))
           deltau <- lambdau*sqrt((1/n1)+(1/n2))
           #calculating the probability X>Y (W)
           #corresponds to the Win Probability P(X*>Y*) where X* and Y* denote future observations
           W <- pnorm(q=delta/sqrt(2))
           Wu <- pnorm(deltau/sqrt(2))
           Wl <- pnorm(deltal/sqrt(2))
           #calculating the odds of X being greater than Y
           Psi <- W/(1-W)
           Psil <- Wl/(1-Wl)
           Psiu <- Wu/(1-Wu)
         },
         greater = {
           #lambdau <- NA#uniroot(BoundLambda, c(-10,10), df=df, cl=1-alpha/2, tobs=tobs)$root#upper 95% confidence limit of delta
           #lambdal <- uniroot(BoundLambda, c(-100,100), df=df, cl=1-alpha, tobs=tobs)$root#lower 95% confidence limit of delta
           lambdau <- 1#qt(p=1-alpha/2, df=df, ncp=tobs)#upper 95% confidence limit of lambda
           lambdal <- qt(p=alpha, df=df, ncp=tobs)#lower 95% confidence limit of lambda
           #calculating the confidence limits on delta
           deltal <- lambdal*sqrt((1/n1)+(1/n2))
           deltau <- NA#lambdau*sqrt((1/n1)+(1/n2))
           #calculating the probability X>Y (W)
           #corresponds to the Win Probability P(X*>Y*) where X* and Y* denote future observations
           W <- pnorm(q=delta/sqrt(2))
           Wu <- NA#pnorm(deltau/sqrt(2))
           Wl <- pnorm(deltal/sqrt(2))
           #calculating the odds of X being greater than Y
           Psi <- W/(1-W)
           Psil <- Wl/(1-Wl)
           Psiu <- NA#Wu/(1-Wu)
         },
         less = {
           #lambdau <- uniroot(BoundLambda, c(-100,100), df=df, cl=alpha, tobs=tobs)$root#upper 95% confidence limit of delta
           #lambdal <- NA#uniroot(BoundLambda, c(-10,10), df=df, cl=alpha, tobs=tobs)$root#lower 95% confidence limit of delta
           lambdau <- qt(p=1-alpha, df=df, ncp=tobs)#upper 95% confidence limit of lambda
           lambdal <- 0#qt(p=alpha/2, df=df, ncp=tobs)#lower 95% confidence limit of lambda
           #calculating the confidence limits on delta
           deltal <- NA#lambdal*sqrt((1/n1)+(1/n2))
           deltau <- lambdau*sqrt((1/n1)+(1/n2))
           #calculating the probability X>Y (W)
           #corresponds to the Win Probability P(X*>Y*) where X* and Y* denote future observations
           W <- pnorm(q=delta/sqrt(2))
           Wu <- pnorm(deltau/sqrt(2))
           Wl <- NA#pnorm(deltal/sqrt(2))
           #calculating the odds of X being greater than Y
           Psi <- W/(1-W)
           Psil <- NA#Wl/(1-Wl)
           Psiu <- Wu/(1-Wu)
         })

  #output
  out <- list(Diff = NA,
              TStat = tobs,
              pvalue=pvalue,
              df = df,
              CI = list(CIl=NA,
                        CIu=NA),
              PI = list(PIl=NA,
                        PIu=NA),
              Cohen=list(Cohenl=deltal,
                         Cohen=delta,
                         Cohenu=deltau),
#              lambda=list(lambdal=lambdal,
#                          lambda=lambda,
#                          lambdau=lambdau),
              W = list(Wl=Wl,
                       W=W,
                       Wu=Wu),
              Cbeta = list(Cbetal=NA,
                           Cbetau=NA),
              Psi = list(Psil=Psil,
                         Psi=Psi,
                         Psiu=Psiu))
  class(out) <- "winprop"
  out
}

#effect(pvalue=0.05, n1=10,n2=10, alternative="two.sided", alpha=0.05)$W
