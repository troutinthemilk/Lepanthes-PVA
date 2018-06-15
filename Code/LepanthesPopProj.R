library('popbio')
library('popdemo')
library('modeest')
graphics.off()

args=(commandArgs(TRUE))
if(length(args)==0) {
    whichPop = 1
} else {
    for(i in 1:length(args)) {
      eval(parse(text=args[[i]]))
    }
}

load(file=paste('Pop', whichPop, 'Matrix.Rdata', sep=''))
print(meanMat)
meanT 	= mean(meanMat)
print(meanT)
sdT 	= sqrt(var2(meanMat))
print(sdT)
lambda 	= lambda(meanT)
cat('lambda', lambda, '\n')
ss 		= stable.stage(meanT)
cat('stable stage', ss, '\n')
rv 		= reproductive.value(meanT)
cat('reproductive value', rv, '\n')
elas 	= elasticity(meanT)
cat('elasticity', '\n')
print(elas)

cat('transfer function analysis',  '\n')
ele.mat <- matrix(NA, 4,4)
ele.mat[1,3] <- "F"
ele.mat[1,1] <- "P"
ele.mat[2,1] <- "P"
ele.mat[2,2] <- "P"
ele.mat[3,2] <- "P"
ele.mat[3,3] <- "P"
ele.mat[4,3] <- "P"
ele.mat[3,4] <- "P"
ele.mat[4,4] <- "P"

#tfml <- tfam_lambda(meanT, elementtype=ele.mat, Flim=c(-0.96241630, 0.9624163), Plim=c(-1, 1))
#plot(tfml)
par(mfrow=c(4,4))
par(mar=c(4.1, 4.1, 1.1, 1.1))
zeros <- c(0, 0, 0, 0)
for(i in 1:4) {
	for(j in 1:4) {
		if(meanT[i,j] <= 1e-10) {
			plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
			next()
		}
		d1 	  <- e1 <- zeros
		d1[i] <- 1
		e1[j] <- 1
		delta <- 1*seq(-meanT[i,j], meanT[i,j], 0.01)
		tf1 <- tfa_lambda(meanT, d = d1, e = e1, prange = delta)
		plot(tf1)#, xlab="Perturbation", ylab=expression(lambda))

		s76 <- sens(meanT)[i, j]
		abline(eigs(meanT, "lambda"), s76, lty = 2)
	}
}

stop()

popSGR 			= stoch.growth.rate(meanMat, verbose=F)
lambdaCI 		= list()
lambdaCI$approx	= exp(popSGR$approx)
lambdaCI$sim 	= exp(popSGR$sim)
lambdaCI$sim.CI	= exp(popSGR$sim.CI)
cat('lambda\n')
print(lambdaCI)

popProject	= stoch.projection(matrices=meanMat, n0=popState[,24], tmax=240, nmax=1000, nreps=1e4, prob=rep(1/24, length=24), verbose=FALSE)
totPopSize 	= apply(popProject, 1, sum)

pdf(file=paste("../Figures/FinalPop_",whichPop,".pdf", sep=""))
hist(totPopSize, breaks=100, col="black", border="white", xlab="Final population size at t=240 (24 years)", main=paste('Sub-population', whichPop))
graphics.off()

mlvPop = mlv(totPopSize, method = "naive")
extin=c()
popQuant = matrix(NA, 3, 241)
popQuant[,1] = rep(sum(popState[,24]), 3)
popMean	= vector('numeric', 241)
popMean[1] = sum(popState[,24])
popSD	= vector('numeric', 241)
popSD[1] = NA

for(i in 1:240){
matriz=c()
interactionX=100000
quasi=2
popExtinctProject<-stoch.projection(matrices=meanMat, n0=popState[,24], tmax=i, nmax=1000, nreps = interactionX, prob = rep(1/24, length=24), verbose=FALSE)
for(ii in 1:interactionX){
		a=sum(popExtinctProject[ii,])
		matriz<-rbind(matriz,c(a))
		}
	vv=matriz[matriz<quasi]
	s=length(vv)
	extin=c(extin,s)
	
	popQuant[,i+1] <- quantile(apply(popExtinctProject, 1, sum), probs=c(0.025, 0.5, 0.975))
	
	popMean[i+1] 	= mean(apply(popExtinctProject, 1, sum))
	popSD[i+1] 		= sd(apply(popExtinctProject, 1, sum))

#	popQuant[,i+1] = mean(apply(popExtinctProject, 1, sum)) - sd(apply(popExtinctProject, 1, sum))
#	popQuant[,i+1] = mean(apply(popExtinctProject, 1, sum)) + sd(apply(popExtinctProject, 1, sum))

}

pdf(file=paste("../Figures/PopQuant_",whichPop,".pdf", sep=""))
plot(popMean, xlab="Time (months)", ylab="Number of individuals (n)", type='l', main=paste("Sub-pop", whichPop), lwd=2, ylim=c(0,max(popQuant, na.rm=T)))
lines(popQuant[1,])
lines(popQuant[3,])
graphics.off() 

pdf(file=paste("../Figures/Extinct_",whichPop,".pdf", sep=""))
plot(extin/interactionX, type='l', lwd=2, xlab="Months into future", main=paste("Extinction probability for sub-population", whichPop))
graphics.off()

save.image(file=paste('Pop', whichPop, 'Output.Rdata', sep=''))
q('no')

