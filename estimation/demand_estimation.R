library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(here)

source(here("functions","shrinkage_functions.R"))
sourceCpp(here("functions","shrinkage_mcmc.cpp"))

load(here("build","output","store_panel.RData"))

childrencounts = countchildren_cpp(tree)
treeindex = createindex(tree)
list = treeindex$list
index = treeindex$index
npar = unlist(lapply(list,length))
L = ncol(tree)
p = nrow(tree)

# --------------------------------------------------------- #
# split data into training/test
# --------------------------------------------------------- #

set.seed(1)
nweeks = nrow(demand)
inweeks = sample(1:nweeks,80)
# inweeks = which(demand$WEEK<unique(demand$WEEK)[41])

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
  
  # training
  # C = matrix(1,nrow=n)
  C = cbind(matrix(1,nrow=n),as.matrix(seasonal[inweeks,]))
  if(any(feati!=0) & n_distinct(feati)>1){
    C = cbind(C,as.matrix(feati))
  }
  if(any(dispi!=0) & n_distinct(dispi)>1){
    C = cbind(C,as.matrix(dispi))
  }
  colnames(C) = NULL
  nphi = nphi + ncol(C)
  nphivec[i] = ncol(C)
  Clist[[i]] = C
  
  # test
  # C = matrix(1,nrow=ntest)
  C = cbind(matrix(1,nrow=ntest),as.matrix(seasonal[-inweeks,]))
  if(any(feati!=0) & n_distinct(feati)>1){
    C = cbind(C,as.matrix(feature[-inweeks,2+i]))
  }
  if(any(dispi!=0) & n_distinct(dispi)>1){
    C = cbind(C,as.matrix(display[-inweeks,2+i]))
  }
  colnames(C) = NULL
  Clisttest[[i]] = C
  
}
cumnphi = c(0,cumsum(nphivec))

# --------------------------------------------------------- #
# sparse models
# --------------------------------------------------------- #

# data
Data = list(
  Y = Y,
  X = X,
  Clist = Clist
)

# priors
Prior = list(
  thetabar = 0,
  betabarii = -3,
  taubarii = 1,
  Aphi = .01*diag(nphi),
  phibar = double(nphi),
  a = 5,
  b = 5
)

# mcmc
Mcmc = list(
  R = 1000,
  keep = 1
)

out_ridge = rSURshrinkage(Data,Prior,Mcmc,Shrinkage="ridge",print=TRUE)
# out_lasso = rSURshrinkage(Data,Prior,Mcmc,Shrinkage="lasso",print=TRUE)
out_horse = rSURshrinkage(Data,Prior,Mcmc,Shrinkage="horseshoe",print=TRUE)

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
  npar = npar
)

# priors
Prior = list(
  thetabar = 0,
  betabarii = -3,
  taubarii = 1,
  Aphi = .01*diag(nphi),
  phibar = double(nphi),
  a = 5,
  b = 5
)

# mcmc
Mcmc = list(
  R = 1000,
  initial_run = 100,
  keep = 1,
  burn_pct = 0.5
)

out_hierridge = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="ridge",group="ridge"),print=TRUE)
# out_hierlasso = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="lasso",group="ridge"),print=TRUE)
out_hierhorse = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="horseshoe",group="ridge"),print=TRUE)

out_hierhorseridge = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="ridge",group="horseshoe"),print=TRUE)
out_hierhorselasso = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="lasso",group="horseshoe"),print=TRUE)
out_hierhorsehorse = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="horseshoe",group="horseshoe"),print=TRUE)

# thetas
matplot(out_hierhorse$thetadraws[,1:npar[1]],type="l",col=1:npar[1])
for(i in 1:npar[1]){
  plot(out_hierhorse$thetadraws[,i],type="l")
  abline(h=0)
}
matplot(out_hierhorse$thetadraws[,(cumsum(npar)[1]+1):cumsum(npar)[2]],type="l")

# beta
wchown = as.vector(diag(p)==1)
hist(apply(out_hierhorse$betadraws,2,mean)[wchown])
hist(apply(out_hierhorse$betadraws,2,mean)[!wchown])

# sigmasq
hist(apply(out_hierhorse$sigmasqdraws,2,mean))
matplot(sqrt(out_hierhorse$sigmasqdraws),type="l")

# tau
matplot(sqrt(out_hierhorse$taudraws),type="l")
