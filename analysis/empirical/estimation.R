library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(here)

source(here("src","shrinkage-functions.R"))
sourceCpp(here("src","shrinkage-mcmc.cpp"))
source(here("analysis","empirical","estimation-build.R"))

# --------------------------------------------------------- #
# hierarchical models
# --------------------------------------------------------- #

# data
Data = list(
  Y = Y,
  X = X,
  Clist = Clist,
  tree = tree,
  parindextree = parindextree,
  parindextree_own = parindextree_own,
  npar = npar,
  npar_own = npar_own
)

# priors
Prior = list(
  thetabar_cross = 0,
  thetabar_own = 0,
  Aphi = .01*diag(nphi),
  phibar = double(nphi),
  a = 5,
  b = 5
)

# mcmc
Mcmc = list(
  R = 500,
  keep = 1
)

out.sparse.ridge = rSURshrinkage(Data,Prior,Mcmc,Shrinkage="ridge",print=TRUE)
out.sparse.lasso = rSURshrinkage(Data,Prior,Mcmc,Shrinkage="lasso",print=TRUE)
out.sparse.horseshoe = rSURshrinkage(Data,Prior,Mcmc,Shrinkage="horseshoe",print=TRUE)

out.ridge.ridge = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="ridge",group="ridge"),print=TRUE)
out.ridge.lasso = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="lasso",group="ridge"),print=TRUE)
out.ridge.horseshoe = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="horseshoe",group="ridge"),print=TRUE)

out.horseshoe.ridge = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="ridge",group="horseshoe"),print=TRUE)
out.horseshoe.lasso = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="lasso",group="horseshoe"),print=TRUE)
out.horseshoe.horseshoe = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="horseshoe",group="horseshoe"),print=TRUE)

# --------------------------------------------------------- #
# visualization
# --------------------------------------------------------- #

wchown = as.vector(diag(p)==1)
end = Mcmc$R/Mcmc$keep
burn = 0.5*end

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

# phi (promos)
phi = apply(out.sparse.ridge$phidraws[burn:end,cumsum(nphivec)[nphivec==max(nphivec)]],2,mean)
hist(phi)
phi = apply(out.ridge.ridge$phidraws[burn:end,cumsum(nphivec)[nphivec==max(nphivec)]],2,mean)
hist(phi)