
InvLogit = function(beta) {
	exp(beta)
}


TransDrawMean = function(drawval, currMonth, currPop, currmcmc) {

	u1x = median(currmcmc[,which(dimnames(currmcmc)[[2]] == paste("u[1,",currMonth,"]",sep=''))])
	u2x = median(currmcmc[,which(dimnames(currmcmc)[[2]] == paste("u[2,",currMonth,"]",sep=''))])
	u3x = median(currmcmc[,which(dimnames(currmcmc)[[2]] == paste("u[3,",currMonth,"]",sep=''))])

	beta1x = median(currmcmc[,which(dimnames(currmcmc)[[2]] == paste("beta[1,",currPop,"]",sep=''))])
	beta2x = median(currmcmc[,which(dimnames(currmcmc)[[2]] == paste("beta[2,",currPop,"]",sep=''))])
	
	t1 = exp(beta1x + u1x )
	t2 = exp(beta2x + u2x )
	t3 = 1 
	tt = t1 + t2 + t3
	
	return(c(t1/tt, t2/tt, t3/tt))
	
}

TransDraw = function(drawval, currMonth, currPop, currmcmc) {

	u1x = currmcmc[drawval,which(dimnames(currmcmc)[[2]] == paste("u[1,",currMonth,"]",sep=''))]
	u2x = currmcmc[drawval,which(dimnames(currmcmc)[[2]] == paste("u[2,",currMonth,"]",sep=''))]
	u3x = currmcmc[drawval,which(dimnames(currmcmc)[[2]] == paste("u[3,",currMonth,"]",sep=''))]

	beta1x = currmcmc[drawval,which(dimnames(currmcmc)[[2]] == paste("beta[1,",currPop,"]",sep=''))]
	beta2x = currmcmc[drawval,which(dimnames(currmcmc)[[2]] == paste("beta[2,",currPop,"]",sep=''))]
	
	t1 = exp(beta1x + u1x )
	t2 = exp(beta2x + u2x )
	t3 = 1 
	tt = t1 + t2 + t3
	
	return(c(t1/tt, t2/tt, t3/tt))
	
}


##constructs the fecundities
fecundDraw = function(numRepro, currMonth, currPop, fdraw) {
  if(numRepro==0) {return(0)}
  
  currRow 	= which(dimnames(fdraw)[[2]] == paste('mu[', currPop, ',', currMonth, ']', sep=''))
  currDraw 	= sample(1:dim(fdraw)[1], 1, replace=TRUE)
  
  return(fdraw[currDraw,currRow]/numRepro)
}

fecundDrawMean = function(numRepro, currMonth, currPop, fdraw, fVal) {
  
  currRow 	= which(dimnames(fdraw)[[2]] == paste('mu[', currPop, ',', currMonth, ']', sep=''))
#   print(currRow)
  return(median(fdraw[,currRow])/(fVal))
}


multinomDraw <- function(popState, tmat, index) {
	
	probs = c(tmat[,index], 1-sum(tmat[,index]))
	temp = rmultinom(n=1, size=popState[index], prob=probs)
	return(temp[1:4,1])

}


constructData = function(state1, state2, master) {

	covMat = NULL

	for(tstep in 6:(dim(master)[2]-1)) {
		temp1 = which(master[,tstep] == state1)	
		temp2 = which(master[,tstep] == state1 & master[,tstep+1] == state2)

		state = vector('numeric', length(temp1))
		for(i in 1:length(temp1)) {
			state[i] = any(temp2==temp1[i])
		}
		monthvec = tstep-5 #(rep(tstep-5,length(temp1))+8)%%12+1
		seasonvec = rep(as.numeric(any(monthvec[1] == c(5:11))), length(monthvec))
		
		newvec = data.frame(Pop=as.factor(master[temp1,1]), Plot=as.factor(master[temp1,2]), Tree=as.factor(master[temp1,3]), Individual=as.factor(master[temp1,5]), Month=as.factor(monthvec), Season=as.factor(seasonvec), State=state)
		covMat = rbind(covMat, newvec)
	}

	return(covMat)

}


MatrixBuild = function(state1, state2, master) {
	
	tableSS = NULL
	temp = constructData(state1, state1, master)
	temp = cbind(Transition=rep(1, dim(temp)[1]), temp)
	rm.vec = which(temp$State == 0)
	temp = temp[-rm.vec,]
	tableSS = rbind(tableSS, temp)

	temp = constructData(state1, state2, master)
	temp = cbind(Transition=rep(2, dim(temp)[1]), temp)
	rm.vec = which(temp$State == 0)
	temp = temp[-rm.vec,]
	tableSS = rbind(tableSS, temp)

	temp = constructData(state1, "Dead", master)
	temp = cbind(Transition=rep(3, dim(temp)[1]), temp)
	rm.vec = which(temp$State == 0)
	temp = temp[-rm.vec,]
	tableSS = rbind(tableSS, temp)

	return(tableSS)

}


constructData = function(state1, state2, master) {

	covMat = NULL

	for(tstep in 6:(dim(master)[2]-1)) {
		temp1 = which(master[,tstep] == state1)	
		temp2 = which(master[,tstep] == state1 & master[,tstep+1] == state2)

		state = vector('numeric', length(temp1))
		for(i in 1:length(temp1)) {
			state[i] = any(temp2==temp1[i])
		}
		monthvec = (rep(tstep-5,length(temp1))) 
		seasonvec = rep(as.numeric(any(monthvec[1] == c(5:11))), length(monthvec))
		
		newvec = data.frame(Pop=as.factor(master[temp1,1]), Plot=as.factor(master[temp1,2]), Tree=as.factor(master[temp1,3]), Individual=as.factor(master[temp1,5]), Month=as.factor(monthvec), Season=as.factor(seasonvec), State=state)
		covMat = rbind(covMat, newvec)
	}

	return(covMat)

}


constructFecund = function(master) {

	covMat = NULL
	dataMat = matrix(NA, 6, 24)
	dimnames(dataMat)[[2]] = 1:24
	dimnames(dataMat)[[1]]= 1:6

	for(tstep in 6:(dim(master)[2]-1)) {

		temp1 = which(master[,tstep] == "Reproductive")	
		temp2 = which(master[,tstep] == "Seedling")

		monthval = (tstep-5+8)%%12+1
		for(i in 1:6) {
			dimnames(dataMat)[[2]][tstep-5] = monthval
			dataMat[i,tstep-5] = length(which(master[temp2,1]==i))/length(which(master[temp1,1]==i))
		}
	}

	if(any(!is.finite(dataMat))) {dataMat[!is.finite(dataMat)]=0}
	return(dataMat)

}


MatrixBuild = function(state1, state2, master) {
	
	tableSS = NULL
	temp = constructData(state1, state1, master)
	temp = cbind(Transition=rep(1, dim(temp)[1]), temp)
	rm.vec = which(temp$State == 0)
	temp = temp[-rm.vec,]
	tableSS = rbind(tableSS, temp)

	temp = constructData(state1, state2, master)
	temp = cbind(Transition=rep(2, dim(temp)[1]), temp)
	rm.vec = which(temp$State == 0)
	temp = temp[-rm.vec,]
	tableSS = rbind(tableSS, temp)

	temp = constructData(state1, "Dead", master)
	temp = cbind(Transition=rep(3, dim(temp)[1]), temp)
	rm.vec = which(temp$State == 0)
	temp = temp[-rm.vec,]
	tableSS = rbind(tableSS, temp)

	return(tableSS)

}


JAGSruns = function(dataTable, samplesize=1e2, thin=5, name) {
	
	y = matrix(0, dim(transTable)[1], 3)
	for(i in 1:dim(transTable)[1]) {
		y[i,transTable$Transition[i]] = 1
	}

	source("logisticMultinomial.bug")
	
	r4 = run.jags(model=logisticREMonth, monitor=c('tau', 'u', 'beta', 'sigma'), n.chains=3, burnin=1e3, sample=samplesize, thin=thin, data = list('Trans'=as.numeric(dataTable$Transition), 'Month'=dataTable$Month, 'Pop'=dataTable$Pop, 'y'=y, 'N' =length(dataTable$Pop)), silent=TRUE, method="parallel")
	dic4 = extract(r4, "DIC")
	cat('***logisticReMonth***\n')
	print(dic4)
	cat('\n')

	save(r4, dic4, file=paste(name,'.Rdata',sep=''))

	return()

	r1 = run.jags(model=logistic, monitor=c('beta'), n.chains=3, burnin=1e3, sample=samplesize, thin=thin, data = list('Trans'=as.numeric(dataTable$Transition), 'Pop'=dataTable$Pop, 'y'=y, 'N' =length(dataTable$Pop)), silent.jags=TRUE, method="parallel")
	dic1 = extract(r1, "DIC")
	cat('***logistic***\n')
	print(dic1)
	cat('\n')
	r2 = run.jags(model=logisticSeason, monitor=c('gamma','beta'), n.chains=3, burnin=1e3, sample=samplesize, thin=thin, data = list('Trans'=as.numeric(dataTable$Transition), 'Month'=dataTable$Month, 'Season'=as.numeric(dataTable$Season),'Pop'=dataTable$Pop, 'y'=y, 'N' =length(dataTable$Pop)), silent=TRUE, method="parallel")
	dic2 = extract(r2, "DIC")
	cat('***logisticSeason***\n')
	print(dic2)
	cat('\n')
	r3 = run.jags(model=logisticMonth, monitor=c('alpha','beta'), n.chains=3, burnin=1e3, sample=samplesize, thin=thin, data = list('Trans'=as.numeric(dataTable$Transition), 'Month'=dataTable$Month, 'Season'=as.numeric(dataTable$Season),'Pop'=dataTable$Pop, 'y'=y, 'N' =length(dataTable$Pop)), silent=TRUE, method="parallel")
	dic3 = extract(r3, "DIC")
	cat('***logisticMonth***\n')
	print(dic3)
	cat('\n')
	r4 = run.jags(model=logisticREMonth, monitor=c('tau', 'u', 'beta'), n.chains=3, burnin=1e3, sample=samplesize, thin=thin, data = list('Trans'=as.numeric(dataTable$Transition), 'Month'=dataTable$Month, 'Pop'=dataTable$Pop, 'y'=y, 'N' =length(dataTable$Pop)), silent=TRUE, method="parallel")
	dic4 = extract(r4, "DIC")
	cat('***logisticReMonth***\n')	
	print(dic4)
	cat('\n')
	r5 = run.jags(model=logisticSeasonREMonth, monitor=c('gamma', 'u', 'tau', 'beta'), n.chains=3, burnin=1e3, sample=samplesize, thin=thin, data = list('Trans'=as.numeric(dataTable$Transition), 'Month'=dataTable$Month, 'Season'=as.numeric(dataTable$Season), 'Pop'=dataTable$Pop, 'y'=y, 'N' =length(dataTable$Pop)), silent=TRUE, method="parallel")
	dic5 = extract(r5, "DIC")
	cat('***logisticSeasonReMonth***\n')	
	print(dic5)
	cat('\n')
	r6 = run.jags(model=logisticRESeason, monitor=c('gamma', 'tau', 'beta'), n.chains=3, burnin=1e3, sample=samplesize, thin=thin, data = list('Trans'=as.numeric(dataTable$Transition), 'Month'=dataTable$Month, 'Season'=as.numeric(dataTable$Season), 'Pop'=dataTable$Pop, 'y'=y, 'N' =length(dataTable$Pop)), silent=TRUE, method="parallel")
	dic6 = extract(r6, "DIC")
	cat('***logisticReSeason***\n')
	print(dic6)
	cat('\n')
	r7 = run.jags(model=logisticRESeasonREMonth, monitor=c('beta','gamma','u', 'tau', 'tau2'), n.chains=3, burnin=1e3, sample=samplesize, thin=thin, data = list('Trans'=as.numeric(dataTable$Transition), 'Month'=dataTable$Month, 'Season'=as.numeric(dataTable$Season), 'Pop'=dataTable$Pop, 'y'=y, 'N' =length(dataTable$Pop)), silent=TRUE, method="parallel")
	dic7 = extract(r7, "DIC")
	cat('***logisticReSeasonReMonth***\n')	
	print(dic7)
	cat('\n')
	
	print('$$$$$$$$$$$$$$$$$$$$$$$$')
	print(c(sum(dic1$dev+dic1$penalty), sum(dic2$dev+dic2$penalty), sum(dic3$dev+dic3$penalty), sum(dic4$dev+dic4$penalty), sum(dic5$dev+dic5$penalty), sum(dic6$dev+dic6$penalty), sum(dic7$dev+dic7$penalty)))
	print('$$$$$$$$$$$$$$$$$$$$$$$$')
	save(r1, dic1, r2, dic2, r3, dic3, r4, dic4, r5, dic5, r6, dic6, r7, dic7, file=paste(name,'.Rdata',sep=''))

}

JAGSfecund = function(fDat, rDat, samplesize=1e3, thin=5, name) {

	source("fecundity.bug")

	r0 = run.jags(model=Poisson, monitor=c('mu'), n.chains=3, burnin=1e3, sample=samplesize, thin=thin, data = list('fDat'=fDat, 'rDat'=rDat, "x"= matrix(as.numeric(fDat>0), 6, 24),"P"=6, "M"=24), silent.jags=TRUE, method="parallel")
	dic0 = extract(r0, "DIC")
	print(dic0)
	
	r1 = run.jags(model=PoissonPop, monitor=c('mu'), n.chains=3, burnin=1e3, sample=samplesize, thin=thin, data = list('fDat'=fDat, 'rDat'=rDat, "x"= matrix(as.numeric(fDat>0), 6, 24),"P"=6, "M"=24), silent.jags=TRUE, method="parallel")
	dic1 = extract(r1, "DIC")
	print(dic1)
	
	r2 = run.jags(model=PoissonPopRemonth, monitor=c('mu'), n.chains=3, burnin=1e3, sample=samplesize, thin=thin, data = list('fDat'=fDat, 'rDat'=rDat, "x"= matrix(as.numeric(fDat>0), 6, 24),"P"=6, "M"=24), silent.jags=TRUE, method="parallel")
	dic2 = extract(r2, "DIC")
	print(dic2)
	
	r3 = run.jags(model=PoissonPopMonth, monitor=c('mu'), n.chains=3, burnin=1e3, sample=samplesize, thin=thin, data = list('fDat'=fDat, 'rDat'=rDat, "x"= matrix(as.numeric(fDat>0), 6, 24),"P"=6, "M"=24), silent.jags=TRUE, method="parallel")
	dic3 = extract(r3, "DIC")
	print(dic3)

	r4 = run.jags(model=zipPoisson, monitor=c('mu'), n.chains=3, burnin=1e3, sample=samplesize, thin=thin, data = list('fDat'=fDat, 'rDat'=rDat, "x"= matrix(as.numeric(fDat>0), 6, 24),"P"=6, "M"=24), silent.jags=TRUE, method="parallel")
	dic4 = extract(r4, "DIC")
	print(dic4)
	
	r5 = run.jags(model=zipPoissonPop, monitor=c('mu'), n.chains=3, burnin=1e3, sample=samplesize, thin=thin, data = list('fDat'=fDat, 'rDat'=rDat, "x"= matrix(as.numeric(fDat>0), 6, 24),"P"=6, "M"=24), silent.jags=TRUE, method="parallel")
	dic5 = extract(r5, "DIC")
	print(dic5)
	
	r6 = run.jags(model=zipPoissonPopRemonth, monitor=c('mu'), n.chains=3, burnin=1e3, sample=samplesize, thin=thin, data = list('fDat'=fDat, 'rDat'=rDat, "x"= matrix(as.numeric(fDat>0), 6, 24),"P"=6, "M"=24), silent.jags=TRUE, method="parallel")
	dic6 = extract(r6, "DIC")
	print(dic6)
	
	r7 = run.jags(model=zipPoissonPopMonth, monitor=c('mu'), n.chains=3, burnin=1e3, sample=samplesize, thin=thin, data = list('fDat'=fDat, 'rDat'=rDat, "x"= matrix(as.numeric(fDat>0), 6, 24),"P"=6, "M"=24), silent.jags=TRUE, method="parallel")
	dic7 = extract(r7, "DIC")
	print(dic7)

	save(r0, dic0, r1, dic1, r2, dic2, r3, dic3, r4, dic4, r5, dic5, r6, dic6, r7, dic7, file=paste(name,'.Rdata',sep=''))

}
