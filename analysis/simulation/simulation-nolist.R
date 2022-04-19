library(Rcpp)
library(RcppArmadillo)
library(here)
library(gtools)
library(prodlim)

sourceCpp(here("src","shrinkage-mcmc-nolist.cpp"))

create_treeobjects = function(tree){
  
  # number of levels
  L = ncol(tree)
  p = nrow(tree)

  # number of parameters per level
  npar = apply(tree,2,max)^2
  npar[1] = npar[1] - p
  npar = c(npar,1)
  
  # number of parameters per level
  npar_own = apply(tree,2,max)
  npar_own = c(npar_own,1)

  # reverse
  tree = tree[,L:1]
 
  # index positions of each element from current parameter vector
  K = max(tree[,1])
  indexpairlast = permutations(K,2,1:K,2,repeats.allowed=TRUE)
  
  # fill in parameter index
  out = NULL
  out[[1]] = rep(1,nrow(indexpairlast))
  
  # repeat for the remaining levels
  if(L>1){
    
    # group-level parameters assumed to be asymmetric, use "combinations" to assume symmetry
    for(i in 2:L){
      
      subtree = unique(tree[,1:i])
      K = max(subtree)
      
      # index positions of each element from current parameter vector
      indexpair = permutations(K,2,1:K,repeats.allowed=TRUE)
      
      # index positions of each element from last level's parameter vector
      indexpairlast = permutations(max(tree[,i-1]),2,1:max(tree[,i-1]),repeats.allowed=TRUE)
      
      # match current positions with last positions
      index = matrix(subtree[indexpair,i-1],ncol=2)
      match = apply(index,1,function(x)row.match(x,data.frame(indexpairlast)))
      
      if(i<L){
        out[[i]] = match
      }
      else{
        own = apply(indexpair,1,function(x)x[1]==x[2])
        out[[i]] = match[!own]
      }
      
    }
    
  }
  
  # list for own parameters
  out_own = list(L)
  out_own[[1]] = rep(1,max(tree[,1]))
  for(l in 1:(L-1)){
    out_own[[l+1]] = unique(tree[,1:(l+1)])[,l]
  }
  
  return(list(npar=npar,npar_own=npar_own,parindextree=rev(out),parindextree_own=rev(out_own)))
  
}

simdata = function(n,p,d,settings){
  
  # observational error
  sigmasq = rep(1,p)
  Sigma = diag(sigmasq)
  
  thetalist = list()
  lambdalist = list()
  
  if(settings$type=="hierarchical"){
    
    tree = settings$tree
    counts = settings$counts
    parindextree = settings$parindextree
    npar = settings$npar
    transform = settings$transform
    L = ncol(tree)
    
    # global variances and initial mean
    tau = rep(1,L)
    mean = 0
    
    # top level
    theta = runif(npar[L],1,3)*sample(c(1,-1),npar[L],replace=TRUE)
    thetalist[[L]] = theta
    lambdalist[[L]] = rep(1,npar[L])
    Psi = lambdalist[[L]] 
    
    # middle levels
    for(ell in (L-1):1){
      
      # dummy variables indicating sparse elements
      sparse_ind = rbinom(npar[ell],1,1-settings$prop_sparse[ell])
      
      # product of previous lambda parameters
      # Psi = create_Psi_ellmone_cpp(lambdalist,counts,list,npar,ell)
      Psi = replace_cpp(parindextree[[ell]],Psi)
      
      # mean as theta from last level
      mean = theta[parindextree[[ell]]]
      
      # draw local variances for current level
      lambda = rep(1,npar[ell])
      
      # construct theta for current level
      if(settings$transform){
        theta = mean + sqrt(lambda*Psi*tau[ell])*rnorm(npar[ell])*sparse_ind
      }
      else{
        theta = (mean + sqrt(lambda*Psi*tau[ell])*rnorm(npar[ell]))*sparse_ind
      }
      
      # save
      thetalist[[ell]] = theta
      lambdalist[[ell]] = lambda*sparse_ind
    }
    
  }
  
  if(settings$type=="sparse"){
    
    mean = 0
    tau = 1
    theta = runif(p^2-p,1,3)*sample(c(1,-1),p^2-p,replace=TRUE)*rbinom(p^2-p,1,1-settings$prop_sparse)

  }
  
  # elasticity matrix
  B = matrix(NA,p,p)
  wchown = which(diag(p)==1)
  betaown = runif(p,-5,-1)
  B[wchown] = betaown
  B[-wchown] = theta
  
  # intercepts and other controls
  Clist = NULL
  phivec = NULL
  Cphi = NULL
  np = double(p)
  for(i in 1:p){
    C = cbind(rep(1,n),matrix(runif(n*(d-1),-1,1),nrow=n))
    np[i] = ncol(C)
    phi = double(d)
    Clist[[i]] = C
    Cphi = c(Cphi,C%*%phi) 
    phivec = c(phivec,phi)
  }
  nphi = sum(np)
  Cphi = matrix(Cphi,n,p)
  
  # data
  X = matrix(rnorm(n*p),n,p)
  error = matrix(rnorm(n*p),n,p)%*%chol(Sigma)
  Y = X%*%B + Cphi + error
  
  # errors at product level
  E = matrix(NA,p,p)
  E[!(diag(p)==1)] = theta - mean
  Bbar = matrix(0,p,p)
  Bbar[!(diag(p)==1)] = mean
  
  return(list(Y=Y,X=X,B=B,Clist=Clist,
              thetalist=thetalist,lambdalist=lambdalist,tau=tau,sigmasq=sigmasq,
              phivec=phivec,nphi=nphi,E=E,Bbar=Bbar))
  
}

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
  R = 100,
  keep = 1
)

sourceCpp(here("src","shrinkage-mcmc-nolist.cpp"))
out = rSURhiershrinkage(Data, Prior, Mcmc, Shrinkage=list(product="horseshoe",group="horseshoe"),print=TRUE)
out = rSURhiershrinkage(Data, Prior, Mcmc, Shrinkage=list(product="ridge",group="horseshoe"),print=TRUE)
out = rSURhiershrinkage(Data, Prior, Mcmc, Shrinkage=list(product="horseshoe",group="ridge"),print=TRUE)
out = rSURhiershrinkage(Data, Prior, Mcmc, Shrinkage=list(product="ridge",group="ridge"),print=TRUE)

# tausq
matplot(sqrt(out$tausqdraws),type="l")
matplot(log(out$tausqdraws),type="l")

cumnpar = cumsum(npar)
matplot(out$lambdasqdraws[,(cumnpar[1]+1):cumnpar[2]],type="l",col=rainbow(npar[2]))
matplot(out$lambdasqdraws[,(cumnpar[2]+1):cumnpar[3]],type="l",col=rainbow(npar[3]))

matplot(log(out$lambdasqdraws[,(cumnpar[1]+1):cumnpar[2]]),type="l",col=rainbow(npar[2]))
matplot(log(out$lambdasqdraws[,(cumnpar[2]+1):cumnpar[3]]),type="l",col=rainbow(npar[3]))
matplot(log(out$Psidraws[,(cumnpar[1]+1):cumnpar[2]]),type="l",col=rainbow(npar[3]))
matplot(log(out$Psidraws[,(cumnpar[2]+1):cumnpar[3]]),type="l",col=rainbow(npar[3]))

cumnpar = cumsum(npar)
wch = c(1,(cumsum(npar)+1)[1:2])

out$lambdasqdraws[10,wch]
out$Psidraws[10,wch]
rev(cumprod(rev(out$lambdasqdraws[10,wch])))

cumnpar = cumsum(npar[-1])

# top
matplot(out$thetadraws[,(cumnpar[1]+1):cumnpar[2]],type="l",col=rainbow(npar[3]))
abline(h=data$thetalist[[3]],col=rainbow(npar[3]))
plot(apply(out$thetadraws[,(cumnpar[1]+1):cumnpar[2]],2,mean),data$thetalist[[3]])
abline(0,1)

# middle
matplot(out$thetadraws[,1:cumnpar[1]],type="l",col=rainbow(npar[2]))
abline(h=data$thetalist[[2]],col=rainbow(npar[2]))
plot(apply(out$thetadraws[,1:cumnpar[1]],2,mean),data$thetalist[[2]])
abline(0,1)

# own
matplot(out$thetaowndraws,type="l")

# matplot(out$betadraws,type="l")
# abline(h=data$B,col=rainbow(p^2))
plot(data$B,apply(out$betadraws,2,mean))
abline(0,1)

# tausq
matplot(log(out$tausqdraws),type="l",col=1:ncol(tree))
matplot(log(out$tausqowndraws),type="l")

