library(extraDistr)

# simulate data from hierarchical prior
simdata_tree = function(n,p,d,tree,counts,list,npar,prop_sparse){
  
  L = ncol(tree)
  
  # observational error
  sigmasq = runif(p,0,1)
  Sigma = diag(sigmasq)
  
  # global variances and initial mean
  tau = rep(1,L)
  tau[L] = 1/rgamma(1,2,2)
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
    
    # product of previous lambda parameters
    Psi = create_Psi_ellmone_cpp(lambdalist,counts,list,npar,ell)
    
    # mean as theta from last level
    mean = theta[list[[ell]]]
    
    # draw local variances for current level
    lambda = rep(1,npar[ell])
    
    # construct theta for current level
    thetasd = sqrt(lambda*Psi*tau[ell])
    theta = mean + rbinom(npar[ell],1,1-prop_sparse[ell])*thetasd*rnorm(npar[ell])
    
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
  Y = X%*%B + error
  
  # errors at product level
  E = matrix(NA,p,p)
  E[!(diag(p)==1)] = theta - mean
  
  return(list(Y=Y,X=X,B=B,Clist=Clist,thetalist=thetalist,
              lambdalist=lambdalist,tau=tau,sigmasq=sigmasq,
              E=E))
  
}

# simulate data from sparse prior
simdata_sparse = function(n,p,d,prop_sparse){
  
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
  B[-wchown] = rnorm(p^2-p,0,sqrt(tau))*rbinom(p^2-p,1,1-prop_sparse)
  
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
  
  return(list(Y=Y,X=X,B=B,Clist=Clist,tau=tau))
  
}

# function to permute a tree
permute_tree = function(tree){
  
  L = ncol(tree)
  for(i in 1:(L-1)){
    groups = unique(tree[,i])
    tree[,i] = rep(groups,times=nrow(tree)/length(groups))
  }
  
  return(tree)
}

# simulate a batch of data sets to estimate in parallel
simdata_batch = function(pvec,nrep,prior,prop_sparse){
  
  d = 1
  
  # generate data
  simdata = NULL
  j = 1
  for(p in p_vec) {
    
    # tree for DGP
    tree_dgp = matrix(c(rep(1:5,each=p/5),
                        rep(1:10,each=p/10),
                        1:p),nrow=p,ncol=3)
    childrencounts_dgp = countchildren_cpp(tree_dgp)
    treeindex_dgp = createindex(tree_dgp)
    list_dgp = treeindex_dgp$list
    index_dgp = treeindex_dgp$index
    npar_dgp = unlist(lapply(list_dgp,length))
    L = ncol(tree_dgp)
    
    if(prior=="tree"){
      
      # wrong tree
      tree_not = permute_tree(tree_dgp)
      childrencounts_not = countchildren_cpp(tree_not)
      treeindex_not = createindex(tree_not)
      list_not = treeindex_not$list
      index_not = treeindex_not$index
      npar_not = unlist(lapply(list_not,length))
      
      simdata[j]  = list(nrep)
      tmp = NULL
      for(b in 1:nrep){
        
        # simulate data
        tmp[[b]] = simdata_tree(n,p,d,tree_dgp,childrencounts_dgp,list_dgp,npar_dgp,prop_sparse)
        tmp[[b]]$npar = npar_dgp
        tmp[[b]]$tree = tree_dgp
        tmp[[b]]$list = list_dgp
        tmp[[b]]$childrencounts = childrencounts_dgp
        tmp[[b]]$npar_not = npar_not
        tmp[[b]]$tree_not = tree_not
        tmp[[b]]$list_not = list_not
        tmp[[b]]$childrencounts_not = childrencounts_not
        
      }
      
      simdata[[j]] = tmp
      j = j+1 
    }
    
    if(prior=="sparse"){
      
      simdata[j]  = list(nrep)
      tmp = NULL
      for(b in 1:nrep){
        
        # simulate data
        tmp[[b]] = simdata_sparse(n,p,d,prop_sparse)
        tmp[[b]]$npar = npar_dgp
        tmp[[b]]$tree = tree_dgp
        tmp[[b]]$list = list_dgp
        tmp[[b]]$childrencounts = childrencounts_dgp
        
      }
      
      simdata[[j]] = tmp
      j = j+1
      
    }
  }
  
  return(simdata)
  
}

# fit models on batch data in parallel
fit_parallel = function(simdata,Mcmc,p_vec,models,dgp){
  
  
  nrep = length(simdata[[1]])
  end = Mcmc$R/Mcmc$keep
  burn = Mcmc$burn_pct*end
  
  itime = proc.time()[3]
  out = foreach(j=1:length(p_vec), .combine='cbind') %:% 
    foreach(b=1:nrep, .combine='rbind') %:% 
    foreach(m=1:nrow(models), .combine='rbind', 
            .packages=c("Rcpp","RcppArmadillo","here"),
            .options.snow=opts) %dopar% { 
              
              sourceCpp(here("functions","shrinkage_mcmc.cpp"))
              
              data = simdata[[j]][[b]]
              p = p_vec[j]
              
              # priors
              Prior = list(
                thetabar=0, 
                taubar=1, 
                betabarii=0, 
                taubarii=10,
                Aphi=.01*diag(p), 
                phibar=double(p), 
                a=5, 
                b=5
              )
              
              # fit model
              if(models[m,2]=="sparse"){
                Data = list(
                  Y=data$Y, 
                  X=data$X, 
                  Clist=data$Clist, 
                  tree=data$tree
                )
                fit = rSURshrinkage(Data,Prior,Mcmc,Shrinkage=models[m,1],print=FALSE)
              }
              else{
                Data = list(
                  Y=data$Y, 
                  X=data$X, 
                  Clist=data$Clist, 
                  tree=data$tree,
                  childrencounts=data$childrencounts, 
                  list=data$list, 
                  npar=data$npar
                )
                fit = rSURhiershrinkage(Data,Prior,Mcmc,
                                        Shrinkage=list(product=models[m,1],
                                                       group=models[m,2]),
                                        print=FALSE)
              }

              # compute rmse
              Brmse = double(end-burn)
              for(r in 1:(end-burn)){
                B = matrix(fit$betadraws[burn+r,],p,p)
                Brmse[r] = sqrt(mean((data$B[diag(p)!=1] - B[diag(p)!=1])^2))
              }
              
              # save 
              mean(Brmse)
              
            }
  
  colnames(out) = paste0("p_",p_vec,"_",dgp)
  
  ctime = proc.time()[3]
  cat("\n Total Time:",round((ctime - itime)/60,2),"minutes",fill = TRUE)
  
  return(out)
  
}


