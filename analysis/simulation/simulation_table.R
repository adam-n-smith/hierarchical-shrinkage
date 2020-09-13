library(gtools)
library(prodlim)
library(extraDistr)
library(ggplot2)
library(reshape2)
library(forcats)
library(Rcpp)
library(RcppArmadillo)

sourceCpp("code/horse_mcmc.cpp")

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
simdata = function(n,p,d,tree,sparse=FALSE){
  
  # observational error
  sigmasq = runif(p,0,1)
  Sigma = diag(sigmasq)
  
  # global variances and initial mean
  # tau = rep(1,L)
  tau = 1/rgamma(L,2,2)
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
  if(sparse){
    theta = rnorm(p^2-p,0,sqrt(tau))*rbinom(p^2-p,1,.05)
    B[-wchown] = theta
  }
  else{
    B[-wchown] = theta
  }
  
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
  
  # test data
  ntest = n*0.5
  Cphitest = NULL
  for(i in 1:p){
    C = cbind(rep(1,ntest),matrix(runif(ntest*(d-1),-1,1),nrow=ntest))
    Cphitest = c(Cphitest,C%*%phivec[(i*d-d+1):(i*d)]) 
  }
  Cphitest = matrix(Cphitest,ntest,p)
  Xtest = matrix(rnorm(ntest*p),ntest,p)
  error = matrix(rnorm(ntest*p),ntest,p)%*%chol(Sigma)
  Ytest = Xtest%*%B + Cphitest + error
  
  return(list(Y=Y,X=X,B=B,Clist=Clist,Ytest=Ytest,Xtest=Xtest))
  
}

nlist = c(50)
plist = c(50,100)

nrep = 25
R = 500
keep = 1
end = R/keep
burn = 0.75*end

yrmse.mean = array(0,dim=c(6,length(nlist),length(plist)))
yrmse.sd = array(0,dim=c(6,length(nlist),length(plist)))
Brmse.mean = array(0,dim=c(6,length(nlist),length(plist)))
Brmse.sd = array(0,dim=c(6,length(nlist),length(plist)))

set.seed(999)
nindex = 1
for(n in nlist){
  pindex = 1
  for(p in plist){
    cat("n =",n," / p =",p," / rep =")
    
    yrmse.ridge = matrix(0,end-burn,nrep)
    yrmse.lasso = matrix(0,end-burn,nrep)
    yrmse.horse = matrix(0,end-burn,nrep)
    yrmse.hierridge = matrix(0,end-burn,nrep)
    yrmse.hierlasso = matrix(0,end-burn,nrep)
    yrmse.hierhorse = matrix(0,end-burn,nrep)
    Brmse.ridge = matrix(0,end-burn,nrep)
    Brmse.lasso = matrix(0,end-burn,nrep)
    Brmse.horse = matrix(0,end-burn,nrep)
    Brmse.hierridge = matrix(0,end-burn,nrep)
    Brmse.hierlasso = matrix(0,end-burn,nrep)
    Brmse.hierhorse = matrix(0,end-burn,nrep)
    
    # tree
    # tree = matrix(c(rep(1:2,each=p/2),
    #                 rep(1:4,each=p/4),
    #                 1:p),nrow=p,ncol=3)
    tree = matrix(c(rep(1:5,each=p/5),
                    rep(1:10,each=p/10),
                    1:p),nrow=p,ncol=3)
    
    childrencounts = countchildren_cpp(tree)
    treeindex = createindex(tree)
    list = treeindex$list
    index = treeindex$index
    npar = unlist(lapply(list,length))
    L = ncol(tree)
    
    for(b in 1:nrep){
      
      # simulate data
      d = 1
      data = simdata(n,p,d,tree)
      # data = simdata(n,p,d,tree,sparse=TRUE)
      Y = data$Y
      X = data$X
      Clist = data$Clist
      Ytest = data$Ytest
      Xtest = data$Xtest
      nphi = p*d
      ntest = n*0.5
      
      # priors
      Prior = list(
        thetabar = 0,
        taubar = 10,
        betabarii = 0,
        taubarii = 10,
        Aphi = .01*diag(nphi),
        phibar = double(nphi),
        a = 5,
        b = 5
      )
      
      # Mcmc
      Mcmc = list(
        R = R,
        keep = keep
      )
      
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
      
      # hierarchical horseshoe
      out.hierridge = rSURhiershrinkage(Data,Prior,Mcmc,shrinkage="ridge",print=FALSE)
      # out.hierlasso = rSURhiershrinkage(Data,Prior,Mcmc,shrinkage="lasso",print=FALSE)
      # out.hierhorse = rSURhorse(Data,Prior,Mcmc,hierarchical=TRUE,propagate=TRUE,print=TRUE)
      out.hierhorse = rSURhiershrinkage(Data,Prior,Mcmc,shrinkage="horseshoe",print=FALSE)
      
      # data
      Data = list(
        Y = Y,
        X = X,
        Clist = Clist
      )
      out.horse = rSURshrinkage(Data,Prior,Mcmc,shrinkage="horseshoe",print=FALSE)
      out.lasso = rSURshrinkage(Data,Prior,Mcmc,shrinkage="lasso",print=FALSE)
      out.ridge = rSURshrinkage(Data,Prior,Mcmc,shrinkage="ridge",print=FALSE)
      
      for(r in 1:(end-burn)){
        
        wchcross = diag(p)!=1
        
        # ridge
        B = matrix(out.ridge$betadraws[burn+r,],p,p)
        Cphi = matrix(out.ridge$phidraws[burn+r,],ntest,p,byrow=TRUE)
        yrmse.ridge[r,b] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))
        Brmse.ridge[r,b] = sqrt(mean((B[wchcross] - data$B[wchcross])^2))

        # lasso
        B = matrix(out.lasso$betadraws[burn+r,],p,p)
        Cphi = matrix(out.lasso$phidraws[burn+r,],ntest,p,byrow=TRUE)
        yrmse.lasso[r,b] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))
        Brmse.lasso[r,b] = sqrt(mean((B[wchcross] - data$B[wchcross])^2))

        # horseshoe
        B = matrix(out.horse$betadraws[burn+r,],p,p)
        Cphi = matrix(out.horse$phidraws[burn+r,],ntest,p,byrow=TRUE)
        yrmse.horse[r,b] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))
        Brmse.horse[r,b] = sqrt(mean((B[wchcross] - data$B[wchcross])^2))
        
        # hierarchical ridge
        B = matrix(out.hierridge$betadraws[burn+r,],p,p)
        Cphi = matrix(out.hierridge$phidraws[burn+r,],ntest,p,byrow=TRUE)
        yrmse.hierridge[r,b] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))
        Brmse.hierridge[r,b] = sqrt(mean((B[wchcross] - data$B[wchcross])^2))

        # hierarchical lasso
        # B = matrix(out.hierlasso$betadraws[burn+r,],p,p)
        # Cphi = matrix(out.hierlasso$phidraws[burn+r,],ntest,p,byrow=TRUE)
        # yrmse.hierlasso[r,b] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))
        # Brmse.hierlasso[r,b] = sqrt(mean((B[wchcross] - data$B[wchcross])^2))

        # hierarchical horseshoe
        B = matrix(out.hierhorse$betadraws[burn+r,],p,p)
        Cphi = matrix(out.hierhorse$phidraws[burn+r,],ntest,p,byrow=TRUE)
        yrmse.hierhorse[r,b] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))
        Brmse.hierhorse[r,b] = sqrt(mean((B[wchcross] - data$B[wchcross])^2))
      }
      
      if(b==nrep){
        cat("",b,"\n")
      }
      else{
        cat("",b)
      }
    }
    
    # save rmse
    yrmse.mean[1,nindex,pindex] = mean(yrmse.ridge)
    yrmse.mean[2,nindex,pindex] = mean(yrmse.lasso)
    yrmse.mean[3,nindex,pindex] = mean(yrmse.horse)
    yrmse.mean[4,nindex,pindex] = mean(yrmse.hierridge)
    yrmse.mean[5,nindex,pindex] = mean(yrmse.hierlasso)
    yrmse.mean[6,nindex,pindex] = mean(yrmse.hierhorse)
    yrmse.sd[1,nindex,pindex] = sd(yrmse.ridge)
    yrmse.sd[2,nindex,pindex] = sd(yrmse.lasso)
    yrmse.sd[3,nindex,pindex] = sd(yrmse.horse)
    yrmse.sd[4,nindex,pindex] = sd(yrmse.hierridge)
    yrmse.sd[5,nindex,pindex] = sd(yrmse.hierlasso)
    yrmse.sd[6,nindex,pindex] = sd(yrmse.hierhorse)
    Brmse.mean[1,nindex,pindex] = mean(Brmse.ridge)
    Brmse.mean[2,nindex,pindex] = mean(Brmse.lasso)
    Brmse.mean[3,nindex,pindex] = mean(Brmse.horse)
    Brmse.mean[4,nindex,pindex] = mean(Brmse.hierridge)
    Brmse.mean[5,nindex,pindex] = mean(Brmse.hierlasso)
    Brmse.mean[6,nindex,pindex] = mean(Brmse.hierhorse)
    Brmse.sd[1,nindex,pindex] = sd(Brmse.ridge)
    Brmse.sd[2,nindex,pindex] = sd(Brmse.lasso)
    Brmse.sd[3,nindex,pindex] = sd(Brmse.horse)
    Brmse.sd[4,nindex,pindex] = sd(Brmse.hierridge)
    Brmse.sd[5,nindex,pindex] = sd(Brmse.hierlasso)
    Brmse.sd[6,nindex,pindex] = sd(Brmse.hierhorse)
    
    # next
    pindex = pindex+1
    
  }
  
  # next
  nindex = nindex+1
  
}

# mean
apply(yrmse.mean,3,c)
apply(Brmse.mean,3,c)

# sd
apply(yrmse.sd,3,c)
apply(Brmse.sd,3,c)

