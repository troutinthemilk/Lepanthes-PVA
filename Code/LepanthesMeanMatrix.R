##constructs the tranasition matrix
source('MatrixFuncs.R')

load("Seedling.Rdata")
seed.mcmc = rbind(r4$mcmc[[1]], r4$mcmc[[2]], r4$mcmc[[3]])
load("Juvenile.Rdata")
juv.mcmc = rbind(r4$mcmc[[1]], r4$mcmc[[2]], r4$mcmc[[3]])
load("Reproductive.Rdata")
reproductive.mcmc = rbind(r4$mcmc[[1]], r4$mcmc[[2]], r4$mcmc[[3]])
load("Adult.Rdata")
adult.mcmc = rbind(r4$mcmc[[1]], r4$mcmc[[2]], r4$mcmc[[3]])
load("Fecundity.Rdata")
fecundity.mcmc = rbind(r5$mcmc[[1]], r5$mcmc[[2]], r5$mcmc[[3]])

#run to 24 months
whichPop 	= 6
master 		= read.csv(file="../Data/LepanthesMasterNew.csv", na.strings="0")

pop=master$POPULATION
year.sims = 1

popState = array(NA, c(4, 24*year.sims))

dimnames(popState)[[1]] = c('Seedling', 'Juvenile', 'Reproductive', 'Adult')

	  popState[,i] = currState
	  
	  tmat = matrix(0, 4, 4)
	  dimnames(tmat)[[1]] = c("Seedling", "Juvenile", "Reproductive", "Adult")
	  dimnames(tmat)[[2]] = c("Seedling", "Juvenile", "Reproductive", "Adult")

	  sT = TransDrawMean(drawval=NA, currMonth=i, currPop=whichPop, currmcmc=seed.mcmc)
	  jT = TransDrawMean(drawval=NA, currMonth=i, currPop=whichPop, currmcmc=juv.mcmc)
	  fT = TransDrawMean(drawval=NA, currMonth=i, currPop=whichPop, currmcmc=reproductive.mcmc)
	  aT = TransDrawMean(drawval=NA, currMonth=i, currPop=whichPop, currmcmc=adult.mcmc)
	  
	  fecundT = fecundDrawMean(numRepro=rDat[whichPop,i], currMonth=i, currPop=whichPop, fdraw=fecundity.mcmc, fVal=fDat[whichPop,i])
	  
	  tmat[1,1] = sT[1]
	  tmat[2,1] = sT[2]
	  tmat[2,2] = jT[1]
	  tmat[3,2] = jT[2]
	  tmat[3,3] = fT[1]
	  tmat[4,3] = fT[2]
	  tmat[4,4] = aT[1]
	  tmat[3,4] = aT[2]
	  tmat[1,3] = fecundT
	  
	  meanMat[[i]] = tmat
	  names(meanMat)[[i]] = paste("T", i, sep='')
}

save(meanMat, popState, file=paste('Pop', whichPop, 'Matrix.Rdata', sep=''))

