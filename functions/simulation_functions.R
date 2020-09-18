library(gtools)
library(prodlim)
library(extraDistr)

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
  if(L>1){
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
    
  }
  
  # save parameter index matrix
  l = out[[L]]
  index = l
  if(L>1){
    for(i in (L-1):1){
      l = out[[i]][l]
      index = rbind(l,index)
    }
  }
  rownames(index) = NULL
  
  # number of parameters per level
  npar = apply(index,1,max)
  if(L>1){
    npar = c(npar[-1],ncol(index))
  }
  
  return(list(list=out,index=index,npar=npar))
  
}
simdata_dense = function(n,p,d,tree,childrencounts,npar,tausq0){
  
  L = ncol(tree)
  
  # observational error
  sigmasq = runif(p,0,1)
  Sigma = diag(sigmasq)
  
  # global variances and initial mean
  tau = rep(1,L)
  tau[L] = tausq0
  mean = 0
  
  # level 1
  lambda = rep(1,npar[1])
  theta = mean + sqrt(tau[1])*rnorm(npar[1])
  thetalist = list()
  lambdalist = list()
  thetalist[[1]] = theta
  lambdalist[[1]] = lambda
  
  # levels 2-L
  for(ell in 2:L){
    
    # product of previous lambda parmaeters
    Psi = create_Psi_ellmone_cpp(lambdalist,childrencounts,list,npar,ell)
    
    # mean as theta from last level
    mean = theta[list[[ell]]]
    
    # draw local variances for current level
    # lambda = runif(npar[ell])
    lambda = 1/rgamma(npar[ell],2,2)
    # lambda = rep(1,npar[ell])
    
    # construct theta for current level
    thetasd = sqrt(lambda)*sqrt(Psi)*sqrt(tau[ell])
    if(ell==L){
      thetasd = rep(sqrt(sigmasq),each=p-1) * thetasd
    }
    theta = mean + thetasd*rnorm(npar[ell])
    
    # save
    thetalist[[ell]] = theta
    lambdalist[[ell]] = lambda
  }
  
  # elasticity matrix
  B = matrix(NA,p,p)
  wchown = which(diag(p)==1)
  betaown = -exp(rnorm(p))
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
  
  # training data
  X = matrix(rnorm(n*p),n,p)
  error = matrix(rnorm(n*p),n,p)%*%chol(Sigma)
  Y = X%*%B + Cphi + error
  
  return(list(Y=Y,X=X,B=B,Clist=Clist,thetalist=thetalist,lambdalist=lambdalist,tau=tau))
  
}
simdata_sparse = function(n,p,d){
  
  npar = 
    
    # observational error
    sigmasq = runif(p,0,1)
  Sigma = diag(sigmasq)
  
  # global variances and initial mean
  tau = 1/rgamma(1,2,2)
  mean = 0
  
  # elasticity matrix
  B = matrix(NA,p,p)
  wchown = which(diag(p)==1)
  betaown = -exp(rnorm(p))
  B[wchown] = betaown
  theta = rnorm(p^2-p,0,sqrt(tau))*rbinom(p^2-p,1,.05)
  B[-wchown] = theta
  
  # intercepts and other controls
  Clist = NULL
  phivec = NULL
  Cphi = NULL
  d = 1
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
  
  # training data
  X = matrix(rnorm(n*p),n,p)
  error = matrix(rnorm(n*p),n,p)%*%chol(Sigma)
  Y = X%*%B + Cphi + error
  
  return(list(Y=Y,X=X,B=B,Clist=Clist,npar=npar))
  
}
