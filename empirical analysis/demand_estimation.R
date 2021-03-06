library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(here)

source(here("functions","shrinkage_functions.R"))
sourceCpp(here("functions","shrinkage_mcmc.cpp"))

# load data
load(here("build","output","store_panel.RData"))

# build tree objects
childrencounts = countchildren_cpp(tree)
treeindex = createindex(tree)
index = treeindex$index
list = treeindex$list
list_own = treeindex$list_own
npar = treeindex$npar
npar_own = treeindex$npar_own
L = ncol(tree)
p = nrow(tree)

# --------------------------------------------------------- #
# split data into training/test
# --------------------------------------------------------- #

nweeks = nrow(demand)
inweeks = demand %>%
  select(IRI_KEY,WEEK) %>%
  rowid_to_column() %>%
  filter(WEEK<=max(WEEK)-26) %>%
  pull(rowid)

Y = demand %>%
  slice(inweeks) %>%
  select(starts_with("UNITS")) %>%
  as.matrix()
Ytest = demand %>%
  slice(-inweeks) %>%
  select(starts_with("UNITS")) %>%
  as.matrix()
X = prices %>%
  slice(inweeks) %>%
  select(starts_with("PRICE")) %>%
  as.matrix()
Xtest = prices %>%
  slice(-inweeks) %>%
  select(starts_with("PRICE")) %>%
  as.matrix()

n = nrow(Y)
ntest = nrow(Ytest)

Clist = NULL
Clisttest = NULL
nphi = 0
nphivec = double(p)
for(i in 1:p){
  
  feati = feature[inweeks,2+i]
  dispi = display[inweeks,2+i]
  promo = feati + dispi
  
  # training
  C = cbind(matrix(1,nrow=n),as.matrix(seasonal[inweeks,c("SUMMER","HOLIDAY")]))
  if(any(promo!=0) & n_distinct(promo)>1){
    C = cbind(C,as.matrix(promo))
  }
  colnames(C) = NULL
  nphi = nphi + ncol(C)
  nphivec[i] = ncol(C)
  Clist[[i]] = C
  
  # test
  C = matrix(1,nrow=ntest)
  C = cbind(matrix(1,nrow=ntest),as.matrix(seasonal[-inweeks,c("SUMMER","HOLIDAY")]))
  if(any(promo!=0) & n_distinct(promo)>1){
    C =  cbind(C,as.matrix(feature[-inweeks,2+i] + display[-inweeks,2+i]))
  }
  colnames(C) = NULL
  Clisttest[[i]] = C
  
}
cumnphi = c(0,cumsum(nphivec))

# --------------------------------------------------------- #
# hierarchical models
# --------------------------------------------------------- #

# data
Data = list(
  Y = Y,
  X = X,
  Clist = Clist,
  tree = tree,
  childrencounts = childrencounts,
  list = list,
  list_own = list_own,
  npar = npar,
  npar_own = npar_own
)

# priors
Prior = list(
  thetabar_cross = 0,
  thetabar_own = 0,
  Aphi = .1*diag(nphi),
  phibar = double(nphi),
  a = 5,
  b = 5
)

# mcmc
Mcmc = list(
  R = 20000,
  initial_run = 0,
  keep = 10
)

out.sparse.ridge = rSURshrinkage(Data,Prior,Mcmc,Shrinkage="ridge",print=TRUE)
out.sparse.lasso = rSURshrinkage(Data,Prior,Mcmc,Shrinkage="lasso",print=TRUE)
out.sparse.horse = rSURshrinkage(Data,Prior,Mcmc,Shrinkage="horseshoe",print=TRUE)

out.ridge.ridge = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="ridge",group="ridge"),print=TRUE)
out.ridge.lasso = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="lasso",group="ridge"),print=TRUE)
out.ridge.horseshoe = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="horseshoe",group="ridge"),print=TRUE)

out.horseshoe.ridge = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="ridge",group="horseshoe"),print=TRUE)
out.horseshoe.lasso = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="lasso",group="horseshoe"),print=TRUE)
out.horseshoe.horseshoe = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="horseshoe",group="horseshoe"),print=TRUE)

# --------------------------------------------------------- #
# visualization
# --------------------------------------------------------- #

end = Mcmc$R/Mcmc$keep
burn = 0.5*end

# out = out.sparse_ridge
# out = out.sparse.lasso
# out = out.sparse.horseshoe
# out = out.ridge.ridge
# out = out.ridge.lasso
# out = out.ridge.horseshoe
# out = out.horseshoe.ridge
# out = out.horseshoe.lasso
# out = out.horseshoe.horseshoe

# theta (own effects)
matplot(out$thetaowndraws[,1:npar_own[1]],type="l",col=rainbow(npar_own[1]))
matplot(out$thetaowndraws[,(npar_own[1]+1):cumsum(npar_own)[2]],type="l",col=rainbow(npar_own[2]))
matplot(out$thetaowndraws[,(cumsum(npar_own)[2]+1):cumsum(npar_own)[3]],type="l",col=rainbow(p))

# thetas (cross effects)
matplot(out$thetadraws[,1:npar[1]],type="l",col=1:npar[1])
for(i in 1:npar[1]){
  plot(out$thetadraws[,i],type="l")
  abline(h=0)
}
matplot(out$thetadraws[,(cumsum(npar)[1]+1):cumsum(npar)[2]],type="l")
matplot(out$thetadraws[,284:300],type="l")

# beta
wchown = as.vector(diag(p)==1)
hist(apply(out$betadraws[,wchown],2,mean), main="Own Elasticities")
hist(apply(out$betadraws[,!wchown],2,mean), main="Cross Elasticities")
data.frame(out$betadraws[,wchown]) %>%
  summarise(across(everything(),
                   list(mean=mean, low=~quantile(.,0.025), high=~quantile(.,0.975)))) %>%
  pivot_longer(everything(),
               names_to = c("set",".value"),
               names_pattern = "(.*)_(.*)") %>%
  mutate(significant = ifelse(low<0 & high>0,"no","yes")) %>%
  ggplot(aes(x=mean, fill=significant)) +
  geom_histogram()

# sigmasq
matplot(sqrt(out$sigmasqdraws),type="l")
hist(apply(out$sigmasqdraws,2,mean))

# tau
matplot(sqrt(out$taudraws),type="l")
matplot(sqrt(out$tauowndraws),type="l")
matplot(log(out$taudraws),type="l")
matplot(log(out$tauowndraws),type="l")

# plot shrinkage
# lambdahat = apply(out$lambdadraws[burn:end,1:npar[1]],2,mean)
# tauhat = mean(out$taudraws[burn:end,1])
# lambdahat = apply(out$lambdadraws[burn:end,(cumsum(npar)[1]+1):cumsum(npar)[2]],2,mean)
# tauhat = mean(out$taudraws[burn:end,2])
lambdahat = apply(out$lambdadraws[burn:end,(cumsum(npar)[2]+1):cumsum(npar)[3]],2,mean)
tauhat = mean(out$taudraws[burn:end,3])
kappa = 1/(1+lambdahat*tauhat)
hist(kappa)

# phi (promos)
phi = apply(out$phidraws[burn:end,cumsum(nphivec)[nphivec==max(nphivec)]],2,mean)
hist(phi)