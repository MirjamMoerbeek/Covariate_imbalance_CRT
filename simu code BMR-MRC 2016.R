require(nlme)

f.imbal=function(beta1, beta2, n1, ICC, imbalance)
	{
	#################################################################################################################################
	### input parameters of function
	#################################################################################################################################
	#beta1                 # unstandardized effect size. Note: coding is 0 control and 1 intervention
	#beta2                 # unstandardized effect covariate. Note: covariate follows standard normal distribution
	#n1                    # cluster size
	#ICC                   # intraclass correlation coefficient
	#imbalance             # choose quantiles from values 0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975

	#################################################################################################################################
	### parameters of simulation study
	#################################################################################################################################
	nr.iter=5000          # number of iterations
	results.matrix=matrix(0,nr.iter,40) # matrix to store results
	alpha=0.05            # type I error rate
	beta=0.2              # type II error rate

	#################################################################################################################################
	### calculate variance components
	#################################################################################################################################
	var.e=1-ICC            # variance person level in a model with both predictors: condition and covariate
	var.u=ICC              # variance cluster level in a model with both predictors: condition and covariate

	#################################################################################################################################
	### calculate number of clusters to achieve desired power level in a test with a two-sided alternative hypothesis
	#################################################################################################################################
	n2=(4*(var.e+n1*var.u)/n1) * ((qnorm(1-alpha/2)+qnorm(1-beta))/(beta1))^2
	n2=4*ceiling(n2/4)      # n2 should be even number
	if(n2<4)
	n2=4
	var=4*(var.e+n1*var.u)/(n1*n2)
	lambda=beta1/sqrt(var)
	power=1-pf(q=qf(p=0.95,df1=1,df2=n2-3),df1=1,df2=n2-3,ncp=lambda^2)

	while(power<0.8)
    		{
		n2=n2+4
    		var=4*(var.e+n1*var.u)/(n1*n2)
    		lambda=beta1/sqrt(var)
		power=1-pf(q=qf(p=0.95,df1=1,df2=n2-3),df1=1,df2=n2-3,ncp=lambda^2)
		}
	N=n1*n2               # total sample size

	#################################################################################################################################
	### generate data that does not vary across iterations
	#################################################################################################################################
	cluster=rep(seq(1,n2),each=n1)
      condition=rep(c(0,1),each=N/2)
      number=qhyper(imbalance,0.5*n2,0.5*n2,0.5*n2)
	if(number==0)                 # constraints to avoid multicollinearity
      number=1
      if(number==0.5*n2)
      number=0.5*n2-1

      covariate=rep(c(0,1,0,1),times=c(number,0.5*n2-number,0.5*n2-number,number))
      covariate=rep(covariate,each=n1)

	#################################################################################################################################
	### repeat nr.iter times
	#################################################################################################################################
	begin=date()
	for (ii in 1:nr.iter)
	      {
		#################################################################################################################################
		### generate data that varies across iterations
		#################################################################################################################################
	      uu=rnorm(n2,0,sqrt(var.u))
	      uu=rep(uu,each=n1)
	      ee=rnorm(N,0,sqrt(var.e))
	      response=beta1*condition+beta2*covariate+uu+ee
		
		#################################################################################################################################
		### analyse data based on ANCOVA model
		#################################################################################################################################
		analysis <- lme(response~1 + condition + covariate, random = ~1|cluster)

	      beta.vector <- fixef(analysis)
	      varcov <- vcov(analysis, useScale = FALSE)
	      se.vector <- sqrt(diag(varcov))
	      test.stat0 <- beta.vector[1] / sqrt(varcov[1,1])
	      p.value0 <- 2 * min(pt(test.stat0,N-n2,lower.tail=FALSE),pt(test.stat0,N-n2,lower.tail=TRUE),1)
	      if (p.value0>0.05) {reject0<-0} else{reject0<-1}
	      test.stat1 <- beta.vector[2] / sqrt(varcov[2,2])
	      p.value1 <- 2 * min(pt(test.stat1,n2-3,lower.tail=FALSE),pt(test.stat1,n2-3,lower.tail=TRUE),1)
	      if (p.value1>0.05) {reject1<-0} else{reject1<-1}
	      test.stat2 <- beta.vector[3] / sqrt(varcov[3,3])
	      p.value2 <- 2 * min(pt(test.stat2,n2-3,lower.tail=FALSE),pt(test.stat2,n2-3,lower.tail=TRUE),1)
	      if (p.value2>0.05) {reject2<-0} else{reject2<-1}

	      var.e.est <- as.numeric(VarCorr(analysis)[2])
	      var.u.est <- as.numeric(VarCorr(analysis)[1])
	      ICC.est <- var.u.est/(var.e.est + var.u.est)

	      results.vector.ANCOVA =  c(beta.vector[1], se.vector[1], test.stat0, p.value0, reject0, beta.vector[2], se.vector[2], test.stat1, p.value1, reject1, beta.vector[3], se.vector[3], test.stat2, p.value2, reject2, analysis$logLik, var.e.est, var.u.est, ICC.est)

		#################################################################################################################################
		### analyse data based on ANOVA model
		#################################################################################################################################
	      analysis <- lme(response~1 + condition, random = ~1|cluster)

	      beta.vector <- fixef(analysis)
	      varcov <- vcov(analysis, useScale = FALSE)
	      se.vector <- sqrt(diag(varcov))
	      test.stat0 <- beta.vector[1] / sqrt(varcov[1,1])
	      p.value0 <- 2 * min(pt(test.stat0,N-n2,lower.tail=FALSE),pt(test.stat0,N-n2,lower.tail=TRUE),1)
	      if (p.value0>0.05) {reject0<-0} else{reject0<-1}
	      test.stat1 <- beta.vector[2] / sqrt(varcov[2,2])
	      p.value1 <- 2 * min(pt(test.stat1,n2-2,lower.tail=FALSE),pt(test.stat1,n2-2,lower.tail=TRUE),1)
	      if (p.value1>0.05) {reject1<-0} else{reject1<-1}

	      var.e.est <- as.numeric(VarCorr(analysis)[2])
	      var.u.est <- as.numeric(VarCorr(analysis)[1])
	      ICC.est <- var.u.est/(var.e.est + var.u.est)

	      results.vector.ANOVA =  c(beta.vector[1], se.vector[1], test.stat0, p.value0, reject0, beta.vector[2], se.vector[2], test.stat1, p.value1, reject1, analysis$logLik, var.e.est, var.u.est, ICC.est)

		#################################################################################################################################
		### store results in row ii of results.matrix
		#################################################################################################################################
	    	results.matrix[ii,]=c(beta1, beta2, n1, ICC, imbalance,n2,power,results.vector.ANCOVA,results.vector.ANOVA)
	     	}
		return(results.matrix)
	}

