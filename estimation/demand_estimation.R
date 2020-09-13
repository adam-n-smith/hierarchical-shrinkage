library(gtools)
library(prodlim)
library(Rcpp)
library(RcppArmadillo)
library(here)

sourceCpp(here("analysis","code","horse_mcmc.cpp"))
load(here("build","output","store_panel.RData"))

createindex = function(tree){
  
  # number of levels
  L = ncol(tree)
  
  # index positions of each element from current parameter vector
  K = max(tree[,1])
  indexpairlast = permutations(K,2,1:K,2,repeats.allowed=TRUE)
  
  # fill in parameter index
  out = NULL
  out[[1]] = rep(1,nrow(indexpairlast))
  
  # repeat for the remaining levels
  for(i in 2:L){
    subtree = unique(tree[,1:i])
    K = max(subtree)
    
    # group-level parameters assumed to be assumtric, use "combinations" to assume symmetry
    if(i<L){
      
      # index positions of each element from current parameter vector
      indexpair = permutations(K,2,1:K,repeats.allowed=TRUE)
      
      # index positions of each element from last level's parameter vector
      indexpairlast = permutations(max(tree[,i-1]),2,1:max(tree[,i-1]),repeats.allowed=TRUE)
      
      # match current positions with last positions
      index = matrix(subtree[indexpair,i-1],ncol=2)
      match = apply(index,1,function(x)row.match(x,indexpairlast))
      out[[i]] = match
    }
    
    # SKU-level parameters
    else{
      
      # index positions of each element from current parameter vector
      indexpair = permutations(K,2,1:K,repeats.allowed=TRUE)
      
      # index positions of each element from last level's parameter vector
      indexpairlast = permutations(max(tree[,i-1]),2,1:max(tree[,i-1]),repeats.allowed=TRUE)
      
      # match current positions with last positions
      index = matrix(subtree[indexpair,i-1],ncol=2)
      match = apply(index,1,function(x)row.match(x,indexpairlast))
      
      # ignore own elasticities
      own = apply(indexpair,1,function(x)x[1]==x[2])
      out[[i]] = match[!own]
    }
  }
  
  return(out)
  
}

childrencounts = countchildren_cpp(tree)
list = createindex(tree)
npar = unlist(lapply(list,length))
p = nrow(tree)

# training/test data
set.seed(999)
nweeks = nrow(demand)
inweeks = sample(1:nweeks,80)

Y = demand %>%
  slice(inweeks) %>%
  select(starts_with("UNITS")) %>%
  as.matrix()
Ytest = demand %>%
  slice(-inweeks) %>%
  select(starts_with("UNITS")) %>%
  as.matrix()
X = prices %>%
  # slice(inweeks) %>%
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
  taubar = 1,
  betabarii = -3,
  taubarii = 1,
  Aphi = .01*diag(nphi),
  phibar = double(nphi),
  a = 5,
  b = 5
)

# Mcmc
Mcmc = list(
  R = 1000,
  keep = 1
)

# sourceCpp("code/horse_mcmc.cpp")
out.hierridge = rSURhiershrinkage(Data,Prior,Mcmc,shrinkage="ridge",print=TRUE)
# out.hierlasso = rSURhiershrinkage(Data,Prior,Mcmc,shrinkage="lasso",print=TRUE)
out.hierhorse = rSURhiershrinkage(Data,Prior,Mcmc,shrinkage="horseshoe",print=TRUE)
# out.hierhorse = rSURhorse(Data,Prior,Mcmc,hierarchical=TRUE,propagate=FALSE,print=TRUE)

# thetas
matplot(out.hierhorse$thetadraws[,1:npar[1]],type="l",col=1:npar[1])
for(i in 1:npar[1]){
  plot(out.hierhorse$thetadraws[,i],type="l")
  abline(h=0)
}
matplot(out.hierhorse$thetadraws[,(cumsum(npar)[1]+2):cumsum(npar)[2]],type="l")

# beta
wchown = as.vector(diag(p)==1)
hist(apply(out.hierhorse$betadraws[burn:end,],2,mean)[wchown])
hist(apply(out.hierhorse$betadraws[burn:end,],2,mean)[!wchown])

# sigmasq
hist(apply(out.hierhorse$sigmasqdraws,2,mean))
matplot(sqrt(out.hierhorse$sigmasqdraws),type="l")

# tau
matplot(sqrt(out.hierhorse$taudraws),type="l")

# --------------------------------------------------------- #
# benchmark models
# --------------------------------------------------------- #

# data
Data = list(
  Y = Y,
  X = X,
  Clist = Clist
)

# priors
Prior = list(
  betabarii = -3,
  taubarii = 1,
  Aphi = .01*diag(nphi),
  phibar = double(nphi),
  a = 5,
  b = 5
)

# Mcmc
Mcmc = list(
  R = 100,
  keep = 1
)

# sourceCpp("code/horse_mcmc.cpp")
out.ridge = rSURshrinkage(Data,Prior,Mcmc,shrinkage="ridge",print=TRUE)
out.lasso = rSURshrinkage(Data,Prior,Mcmc,shrinkage="lasso",print=TRUE)
out.horse = rSURshrinkage(Data,Prior,Mcmc,shrinkage="horseshoe",print=TRUE)

# beta
wchown = as.vector(diag(p)==1)
hist(apply(out.ridge$betadraws,2,mean)[wchown])
hist(apply(out.ridge$betadraws,2,mean)[!wchown])

# sigmasq
hist(apply(out.horse$sigmasqdraws,2,mean))

# tau
matplot(sqrt(out.horse$taudraws),type="l")

# sigmasq
matplot(sqrt(out.horse$lambdadraws),type="l")
