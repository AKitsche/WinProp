#Calculating the confidence interval for the Win Probability P(X-Y > c) according to Hayter 
#Calculating the confidence interval for margin c with probability beta 

#References
#A. J. Hayter (2013): Inferences on the difference between future observations for comparing two treatments, 
#Journal of Applied Statistics, 40:4, 887-900

#Function arguments:
#X - mean of the first group
#Y - mean of the second group
#s - pooled standard deviation, if var.equal=TRUE
#  - standard deviation of the first group, if var.equal=FALSE
#n1 - sample size of the first group
#n2 - sample size of the second group
#alpha - type-I error rate
#beta - chance that the first group is different from the second group
#c - value c of interest, WinProp(c) directly assesses the probability that the difference between potential future single observations will be at least that value
#var.equal - specifies if homogeneous variances between the two groups is assumed
#s1 - standard deviation of the second group, if var.equal=FALSE
#alternative - specifies the direction of the hypothesis
#m1 - number of potential future observations from treatment 1
#m2 - number of potential future observations from treatment 2

WinPropSum <- function(X, Y, s, n1, n2, alpha=0.05, beta=0.95, c=0, var.equal=TRUE, s2=1, alternative=c("two.sided", "greater","less"), m1=1, m2=1){  
  #checks
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
  #check for whole numbers of sample sizes  
  if(!is.whole(n1) | !is.whole(n2) | n1<0 | n2<0){
    stop("n1 and n2 must be positive whole numbers")
  }
  if(length(s) != 1 | !is.numeric(s) | s < 0) {
    stop("s must be a single numeric value")}
  
  if(length(m1) != 1 | length(m2) != 1 | !is.whole(m1) | !is.whole(m2)){
    stop("m1 and m2 must be be positive whole numbers")
  }
  
  alternative <- match.arg(alternative)
  #calculating the t-Teststatistic
  Diff <- X-Y#difference in means
  CohenSize <- Diff/s
  if(var.equal==TRUE){
    df <- n1+n2-2#degrees of freedom
    TStat <- (Diff-c)/(s*sqrt((1/n1)+(1/n2)))
  }else{
    df <- ((s^2/n1)+(s2^2/n2))^2/(s^4/(n1^2*(n1-1))+s2^4/(n2^2*(n2-1)))
    TStat <- (Diff-c)/sqrt((s^2/n1)+(s2^2/n2))
  }
  #calculating the p-value
  alternative <- match.arg(alternative)
  switch(alternative, 
         two.sided = {
           pvalue <- 2*(1-pt(q=abs(TStat), df=df))
         },
         greater = {
           pvalue <- 1-pt(q=TStat, df=df)
         },
         less = {
           pvalue <- pt(q=TStat, df=df)
         })
  #calculating the confidence interval and prediction intervals for the diference of means, and Cohense effect size
  switch(alternative, 
         two.sided = {
           tq <- qt(p=1-(alpha/2), df=df)#quantile of the t-distribution
           #s <- sqrt(((var(x)*(n1-1))+ (var(y)*(n2-1)))/(n1+n2-2))
           CIl <- Diff-(tq*sqrt((s^2/n1)+(s2^2/n2)))#lower confidence limit
           CIu <- Diff+(tq*sqrt((s^2/n1)+(s2^2/n2)))#upper confidence limit
           PIl <- Diff-(tq*(sqrt(((s^2*(n1))+ (s2^2*(n2)))/(n1+n2))*sqrt(2+(1/n1)+(1/n2))))#lower limit
           PIu <- Diff+(tq*(sqrt(((s^2*(n1))+ (s2^2*(n2)))/(n1+n2))*sqrt(2+(1/n1)+(1/n2))))#upper limit
           lambdau <- qt(p=1-alpha/2, df=df, ncp=TStat)#upper 95% confidence limit of lambda
           lambdal <- qt(p=alpha/2, df=df, ncp=TStat)#lower 95% confidence limit of lambda
           #calculating the confidence limits on delta
           deltal <- lambdal*sqrt((1/n1)+(1/n2))
           deltau <- lambdau*sqrt((1/n1)+(1/n2))
         },
         greater = {
           tq <- qt(p=1-(alpha), df=df)#quantile of the t-distribution
           #s <- sqrt(((var(x)*(n1-1))+ (var(y)*(n2-1)))/(n1+n2-2))
           CIl <- Diff-(tq*sqrt((s^2/n1)+(s2^2/n2)))#lower confidence limit
           CIu <- Inf#Diff+(tq*(s*sqrt((1/n1)+(1/n2))))#upper confidence limit
           PIl <- Diff-(tq*(sqrt(((s^2*(n1))+ (s2^2*(n2)))/(n1+n2))*sqrt(2+(1/n1)+(1/n2))))#lower limit
           PIu <- Inf#Diff+(tq*(sqrt(((var(x)*(n1-1))+ (var(y)*(n2-1)))/(n1+n2-2))*sqrt(2+(1/n1)+(1/n2))))#upper limit
           lambdau <- NA#qt(p=1-alpha, df=df, ncp=TStat)#upper 95% confidence limit of lambda
           lambdal <- qt(p=alpha, df=df, ncp=TStat)#lower 95% confidence limit of lambda
           #calculating the confidence limits on delta
           deltal <- lambdal*sqrt((1/n1)+(1/n2))
           deltau <- -Inf#lambdau*sqrt((1/n1)+(1/n2))
         },
         less = {
           tq <- qt(p=1-alpha, df=df)#quantile of the t-distribution
           #s <- sqrt(((var(x)*(n1-1))+ (var(y)*(n2-1)))/(n1+n2-2))
           CIl <- -Inf#Diff-(tq*(s*sqrt((1/n1)+(1/n2))))#lower confidence limit
           CIu <- Diff+(tq*sqrt((s^2/n1)+(s2^2/n2)))#upper confidence limit
           PIl <- -Inf#Diff-(tq*(sqrt(((var(x)*(n1-1))+ (var(y)*(n2-1)))/(n1+n2-2))*sqrt(2+(1/n1)+(1/n2))))#lower limit
           PIu <- Diff+(tq*(sqrt(((s^2*(n1))+ (s2^2*(n2)))/(n1+n2))*sqrt(2+(1/n1)+(1/n2))))#upper limit
           lambdau <- qt(p=1-alpha, df=df, ncp=TStat)#upper 95% confidence limit of lambda
           lambdal <- NA#qt(p=alpha, df=df, ncp=TStat)#lower 95% confidence limit of lambda
           #calculating the confidence limits on delta
           deltal <- Inf#lambdal*sqrt((1/n1)+(1/n2))
           deltau <- lambdau*sqrt((1/n1)+(1/n2))
         })
  
  switch(alternative, 
         two.sided = {
           if(var.equal==TRUE){
             df <- n1+n2-2#degrees of freedom
             Stat <- (X-Y-c)/(s*sqrt((1/n1)+(1/n2)))
             deltal2 <- qt(p=alpha/2, df=df, ncp=Stat) 
             deltau2 <- qt(p=1-(alpha/2), df=df, ncp=Stat) 
             W <- pnorm((X-Y-c)/(sqrt(1/m1 + 1/m2)*s))
             Wl <- pnorm((sqrt((1/n1)+(1/n2)))/(sqrt((1/m1)+(1/m2)))*deltal2)
             Wu <- pnorm((sqrt((1/n1)+(1/n2)))/(sqrt((1/m1)+(1/m2)))*deltau2)
             Cbetal <- (X-Y)- (qt(p=1-alpha/2, df=df)*(s*sqrt((1/n1)+(1/n2))))-(qnorm(p=beta)*sqrt((s^2/m1) + (s^2/m2)))
             Cbetau <- (X-Y)+ (qt(p=1-alpha/2, df=df)*(s*sqrt((1/n1)+(1/n2))))-(qnorm(p=beta)*sqrt((s^2/m1) + (s^2/m2)))
           }else{
             se <- sqrt((s^2/n1)+(s2^2/n2))/sqrt((s^2/m1)+(s2^2/m2))
             Stat <- (X-Y-c)/sqrt((s^2/m1) + (s2^2/m2))
             df <- ((s^2/n1)+(s2^2/n2))^2/(s^4/(n1^2*(n1-1))+s2^4/(n2^2*(n2-1)))
             deltal2 <- qt(p=alpha/2, df=df, ncp=Stat) 
             deltau2 <- qt(p=1-(alpha/2), df=df, ncp=Stat) 
             W <- pnorm((X-Y-c)/sqrt((s^2/m1) + (s2^2/m2)))
             Wl <- pnorm(Stat-qt(p=1-alpha/2, df=df)*se)
             Wu <- pnorm(Stat+qt(p=1-alpha/2, df=df)*se)
             Cbetal <- (X-Y)-(qt(p=1-alpha/2, df=df)*sqrt((s^2/n1)+(s2^2/n2)))-(qnorm(p=beta)*sqrt((s^2/m1) + (s2^2/m2)))
             Cbetau <- (X-Y)+(qt(p=1-alpha/2, df=df)*sqrt((s^2/n1)+(s2^2/n2)))-(qnorm(p=beta)*sqrt((s^2/m1) + (s2^2/m2)))
           }
         },
         greater = {
           if(var.equal==TRUE){
             df <- n1+n2-2#degrees of freedom
             Stat <- (X-Y-c)/(s*sqrt((1/n1)+(1/n2)))
             deltal2 <- qt(p=alpha, df=df, ncp=Stat) 
             #deltau <- Inf#qt(p=1-(alpha/2), df=df, ncp=Stat) 
             W <- pnorm((X-Y-c)/(sqrt(1/m1 + 1/m2)*s))
             Wl <- pnorm((sqrt((1/n1)+(1/n2)))/(sqrt((1/m1)+(1/m2)))*deltal2)
             Wu <- 1#pnorm((sqrt((1/n1)+(1/n2)))/(sqrt((1/m1)+(1/m2)))*deltau)
             Cbetal <- (X-Y)- (qt(p=1-alpha, df=df)*(s*sqrt((1/n1)+(1/n2))))-(qnorm(p=beta)*sqrt((s^2/m1) + (s^2/m2)))
             Cbetau <- Inf#(X-Y)+ (qt(p=1-alpha/2, df=df)*(s*sqrt((1/n1)+(1/n2))))-(qnorm(p=beta)*sqrt((s^2/m1) + (s^2/m2)))
           }else{
             se <- sqrt((s^2/n1)+(s2^2/n2))/sqrt((s^2/m1)+(s2^2/m2))
             Stat <- (X-Y-c)/sqrt((s^2/m1) + (s2^2/m2))
             df <- ((s^2/n1)+(s2^2/n2))^2/(s^4/(n1^2*(n1-1))+s2^4/(n2^2*(n2-1)))
             W <- pnorm((X-Y-c)/sqrt((s^2/m1) + (s2^2/m2)))
             Wl <- pnorm(Stat-qt(p=1-alpha, df=df)*se)
             Wu <- 1#dnorm(Stat+qt(p=alpha/2, df=df)*se)
             Cbetal <- (X-Y)-(qt(p=1-alpha, df=df)*sqrt((s^2/n1)+(s2^2/n2)))-(qnorm(p=beta)*sqrt((s^2/m1) + (s2^2/m2)))
             Cbetau <- Inf#(X-Y)+(qt(p=1-alpha/2, df=df)*sqrt((s^2/n1)+(s2^2/n2)))-(qnorm(p=beta)*sqrt((s^2/m1) + (s2^2/m2)))
           }
         },
         less = {
           if(var.equal==TRUE){
             df <- n1+n2-2#degrees of freedom
             Stat <- (X-Y-c)/(s*sqrt((1/n1)+(1/n2)))
             #deltal <- -Inf#qt(p=alpha/2, df=df, ncp=Stat) 
             deltau2 <- qt(p=1-alpha, df=df, ncp=Stat) 
             W <- pnorm((X-Y-c)/(sqrt(1/m1 + 1/m2)*s))
             Wl <- 0#pnorm((sqrt((1/n1)+(1/n2)))/(sqrt((1/m1)+(1/m2)))*deltal)
             Wu <- pnorm((sqrt((1/n1)+(1/n2)))/(sqrt((1/m1)+(1/m2)))*deltau2)
             Cbetal <- -Inf#(X-Y)- (qt(p=1-alpha/2, df=df)*(s*sqrt((1/n1)+(1/n2))))-(qnorm(p=beta)*sqrt((s^2/m1) + (s^2/m2)))
             Cbetau <- (X-Y)+ (qt(p=1-alpha, df=df)*(s*sqrt((1/n1)+(1/n2))))-(qnorm(p=beta)*sqrt((s^2/m1) + (s^2/m2)))
           }else{
             se <- sqrt((s^2/n1)+(s2^2/n2))/sqrt((s^2/m1)+(s2^2/m2))
             Stat <- (X-Y-c)/sqrt((s^2/m1) + (s2^2/m2))
             df <- ((s^2/n1)+(s2^2/n2))^2/(s^4/(n1^2*(n1-1))+s2^4/(n2^2*(n2-1)))
             #deltal <- -Inf#qt(p=alpha/2, df=df, ncp=Stat) 
             #deltau <- qt(p=1-alpha, df=df, ncp=Stat) 
             W <- pnorm((X-Y-c)/sqrt((s^2/m1) + (s2^2/m2)))
             Wl <- 0#dnorm(Stat-qt(p=alpha/2, df=df)*se)
             Wu <- pnorm(Stat+qt(p=1-alpha, df=df)*se)
             Cbetal <- -Inf#(X-Y)-(qt(p=1-alpha/2, df=df)*sqrt((s^2/n1)+(s2^2/n2)))-(qnorm(p=beta)*sqrt((s^2/m1) + (s2^2/m2)))
             Cbetau <- (X-Y)+(qt(p=1-alpha, df=df)*sqrt((s^2/n1)+(s2^2/n2)))-(qnorm(p=beta)*sqrt((s^2/m1) + (s2^2/m2)))
           }
         })
  #calculating the odds W/(1-W) and corresponding confidence interval
  Phi <- W/(1-W)
  Phil <- Wl/(1-Wl)
  Phiu <- Wu/(1-Wu)
  out <- list(Diff=Diff,
              TStat=TStat,
              pvalue=pvalue,
              df=df,
              CI = list(CIl=CIl,
                        CIu=CIu),
              PI = list(PIl=PIl,
                        PIu=PIu),
              Cohen=list(Cohenl=deltal,
                         Cohen=CohenSize,
                         Cohenu=deltau),
              W=list(W=W,
                     Wl=Wl,
                     Wu=Wu),
              Cbeta=list(Cbetal=Cbetal,
                         Cbetau=Cbetau),
              Phi = list(Phi=Phi,
                         Phil=Phil,
                         Phiu=Phiu))
  class(out) <- "winprop"
  out
}
