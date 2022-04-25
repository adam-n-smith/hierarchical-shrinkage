wchown = as.vector(diag(p)==1)
end = Mcmc$R/Mcmc$keep
burn = 0.5*end

par(mfrow=c(2,1))

beta.sparse.ridge = apply(out.sparse.ridge$betadraws[burn:end,],2,mean)
hist(beta.sparse.ridge[wchown], main="Own Elasticities",xlim=c(-8,4))
summary(beta.sparse.ridge[wchown])
# hist(beta.sparse.ridge[-wchown], main="Cross Elasticities",xlim=c(-2,2),breaks=20)
# summary(beta.sparse.ridge[-wchown])

beta.ridge.ridge = apply(out.ridge.ridge$betadraws[burn:end,],2,mean)
hist(beta.ridge.ridge[wchown], main="Own Elasticities",xlim=c(-8,4))
summary(beta.ridge.ridge[wchown])
# hist(beta.ridge.ridge[-wchown], main="Cross Elasticities",xlim=c(-2,2),breaks=20)
# summary(beta.ridge.ridge[-wchown])

plot(beta.sparse.ridge,beta.ridge.ridge,col=wchown+1,pch=16)
abline(0,1)

# top
wchtop = (cumsum(npar_own[-1])[1]+1):cumsum(npar_own[-1])[2]
matplot(out.ridge.ridge$thetaowndraws[,wchtop],type="l")
data.frame(name = unique(tree_names$LARGE_CATEGORY),
           theta = apply(as.matrix(out.ridge.ridge$thetaowndraws[burn:end,wchtop]),2,mean))

# middle
wchmid = 1:cumsum(npar_own[-1])[1]
matplot(out.ridge.ridge$thetaowndraws[,wchmid],type="l")
data.frame(name = unique(tree_names$SMALL_CATEGORY),
           theta = apply(out.ridge.ridge$thetaowndraws[burn:end,wchmid],2,mean))

# UPC
matplot(out.sparse.ridge$betadraws[,wchown],type="l")
summary(apply(out.sparse.ridge$betadraws[,wchown],2,mean))
matplot(out.ridge.ridge$betadraws[,wchown],type="l")
summary(apply(out.ridge.ridge$betadraws[,wchown],2,mean))

# tausq (own)
matplot(log(out.sparse.ridge$tausqowndraws),type="l")
matplot(log(out.ridge.ridge$tausqowndraws),type="l")

# tausq (cross)
matplot(log(out.sparse.ridge$tausqdraws),type="l")
matplot(log(out.ridge.ridge$tausqdraws),type="l")
