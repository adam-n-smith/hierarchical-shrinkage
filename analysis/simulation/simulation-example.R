library(tidyverse)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(here)

source(here("src","simulation-functions.R"))
source(here("src","shrinkage-functions.R"))
sourceCpp(here("src","shrinkage-mcmc.cpp"))

# ---------------------------------------------------------------------------- #
# generate data
# ---------------------------------------------------------------------------- #

# dimensions
n = 50
p = 100
d = 1

# tree
tree = matrix(c(1:p,
                rep(1:10,each=p/10),
                rep(1:5,each=p/5)),nrow=p,ncol=3)
objects = create_treeobjects(tree)
parindextree = objects$parindextree
parindextree_own = objects$parindextree_own
npar = objects$npar
npar_own = objects$npar_own

# beta generated from hierarchical prior
settings = list(
  type = "hierarchical",
  tree = tree,
  parindextree = parindextree,
  npar = npar,
  prop_sparse = c(0,0,0),
  transform = FALSE
)

set.seed(1)
data = simdata(n,p,d,settings)

# ---------------------------------------------------------------------------- #
# estimation
# ---------------------------------------------------------------------------- #

# data
Data = list(
  Y = data$Y,
  X = data$X, 
  Clist = data$Clist, 
  tree = tree,
  parindextree = parindextree,
  parindextree_own = parindextree_own,
  npar = npar,
  npar_own = npar_own
)

# priors
Prior = list(
  thetabar_own = 0, 
  thetabar_cross = 0, 
  Aphi = .01*diag(p), 
  phibar = double(p),
  a = 5, 
  b = 5
)

# mcmc
Mcmc = list(
  R = 1000,
  keep = 1
)

out = rSURshrinkage(Data, Prior, Mcmc, Shrinkage="ridge",print=TRUE)
out = rSURshrinkage(Data, Prior, Mcmc, Shrinkage="lasso",print=TRUE)
out = rSURshrinkage(Data, Prior, Mcmc, Shrinkage="horseshoe",print=TRUE)

out = rSURhiershrinkage(Data, Prior, Mcmc, Shrinkage=list(product="ridge",group="ridge"),print=TRUE)
out = rSURhiershrinkage(Data, Prior, Mcmc, Shrinkage=list(product="lasso",group="ridge"),print=TRUE)
out = rSURhiershrinkage(Data, Prior, Mcmc, Shrinkage=list(product="horseshoe",group="ridge"),print=TRUE)

out = rSURhiershrinkage(Data, Prior, Mcmc, Shrinkage=list(product="ridge",group="horseshoe"),print=TRUE)
out = rSURhiershrinkage(Data, Prior, Mcmc, Shrinkage=list(product="lasso",group="horseshoe"),print=TRUE)
out = rSURhiershrinkage(Data, Prior, Mcmc, Shrinkage=list(product="horseshoe",group="horseshoe"),print=TRUE)

# ---------------------------------------------------------------------------- #
# model summary
# ---------------------------------------------------------------------------- #

end = Mcmc$R/Mcmc$keep
burn = 0.5*end
cumnpar = cumsum(npar[-1])

# theta - top
matplot(out$thetadraws[,(cumnpar[1]+1):cumnpar[2]],type="l",col=rainbow(npar[3]))
abline(h=data$thetalist[[3]],col=rainbow(npar[3]))
plot(apply(out$thetadraws[,(cumnpar[1]+1):cumnpar[2]],2,mean),data$thetalist[[3]])
abline(0,1)

# theta - middle
matplot(out$thetadraws[,1:cumnpar[1]],type="l",col=rainbow(npar[2]))
abline(h=data$thetalist[[2]],col=rainbow(npar[2]))
plot(apply(out$thetadraws[,1:cumnpar[1]],2,mean),data$thetalist[[2]])
abline(0,1)

# beta
matplot(out$betadraws,type="l",col=rainbow(p^2))
abline(h=data$B,col=rainbow(p^2))
plot(data$B,apply(out$betadraws,2,mean),xlab="true",ylab="estimated")
abline(0,1)

# beta: check share of correct signs
mean(sign(data$B) == sign(apply(out$betadraws[burn:end,],2,mean)))

# tausq
matplot(log(out$tausqdraws),type="l",col=1:ncol(tree))

# phi
matplot(out$phidraws,type="l",col=1:length(data$phivec))
abline(h=data$phivec,col=1:length(data$phivec))
plot(data$phivec,apply(out$phidraws[burn:end,],2,mean),xlab="true",ylab="estimated")
abline(0,1)

# plot shrinkage
kappa = apply(1/(1+out$tausqdraws[burn:end,1]*out$Psidraws[burn:end,1:p^2]),2,mean)
hist(kappa)