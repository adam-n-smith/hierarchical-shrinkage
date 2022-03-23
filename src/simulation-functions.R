library(extraDistr)

# simulate data
simdata = function(n,p,d,settings){

  # observational error
  # sigmasq = runif(p,0,1)
  sigmasq = runif(p,1,1)
  Sigma = diag(sigmasq)
  
  thetalist = list()
  lambdalist = list()
  
  if(settings$type=="hierarchical"){
    
    tree = settings$tree
    counts = settings$counts
    list = settings$list
    npar = settings$npar
    transform = settings$transform
    L = ncol(tree)

    # global variances and initial mean
    tau = rep(1,L)
    mean = 0
    
    # level 1
    theta = runif(npar[1],1,3)*sample(c(1,-1),npar[1],replace=TRUE)
    thetalist[[1]] = theta
    lambdalist[[1]] = rep(1,npar[1])
    
    # levels 2-L
    for(ell in 2:L){
      
      # dummy variables indicating sparse elements
      sparse_ind = rbinom(npar[ell],1,1-settings$prop_sparse[ell])
      
      # product of previous lambda parameters
      Psi = create_Psi_ellmone_cpp(lambdalist,counts,list,npar,ell)
      
      # mean as theta from last level
      mean = theta[list[[ell]]]
      
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
    # theta = runif(p^2-p,-3,3)*rbinom(p^2-p,1,1-settings$prop_sparse)
    
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

# simulate a batch of data sets to estimate in parallel
simdata_batch = function(n,p,rep,type,prop_sparse,transform=NULL){
  
  d = 1
  
  # generate data
  data = NULL

  # tree for DGP
  if(type=="hierarchical"){
    
    tree = matrix(c(rep(1:5,each=p/5),
                    rep(1:10,each=p/10),
                    1:p),nrow=p,ncol=3)
    counts = countchildren_cpp(tree)
    treeindex = createindex(tree)
    list = treeindex$list
    list_own = treeindex$list_own
    npar = treeindex$npar
    npar_own = treeindex$npar_own
    
    settings = list(
      type = type,
      tree = tree,
      counts = counts,
      list = list,
      list_own = list_own,
      npar = npar,
      npar_own = npar_own,
      prop_sparse = prop_sparse,
      transform = transform
    )
    
  }
  
  if(type=="sparse"){
    settings = list(
      type = type,
      prop_sparse = prop_sparse
    )
  }
  
  data = list(rep)
  tmp = NULL
  for(b in 1:rep){
    
    # simulate data
    data[[b]] = simdata(n,p,d,settings)
    
    if(type=="hierarchical"){
      
      data[[b]]$npar = npar
      data[[b]]$npar_own = npar_own
      data[[b]]$tree = tree
      data[[b]]$list = list
      data[[b]]$list_own = list_own
      data[[b]]$childrencounts = counts
      
    }
    
  }
  
  return(data)
  
}


# fit models on batch data in parallel
fit_parallel = function(DataList,Prior,Mcmc,p,models,dgp,tree=NULL){
  
  rep = length(DataList)
  end = Mcmc$R/Mcmc$keep
  burn = Mcmc$burn_pct*end
  
  # true tree
  tree_dgp = data[[which(str_detect(names(data),"dense"))[1]]][[1]]$tree
  childrencounts_dgp = countchildren_cpp(tree_dgp)
  treeindex = createindex(tree_dgp)
  index_dgp = treeindex$index
  list_dgp = treeindex$list
  list_own_dgp = treeindex$list_own
  npar_dgp = treeindex$npar
  npar_own_dgp = treeindex$npar_own
  L_dgp = ncol(tree_dgp)
  
  # misspecified tree
  tree_mis = tree_dgp
  for(k in 1:(ncol(tree_mis)-1)){
    K = max(tree_mis[,k])
    tree_mis[,k] = rep(1:K,times=p/K)
  }
  childrencounts_mis = countchildren_cpp(tree_mis)
  treeindex = createindex(tree_mis)
  index_mis = treeindex$index
  list_mis = treeindex$list
  list_own_mis = treeindex$list_own
  npar_mis = treeindex$npar
  npar_own_mis = treeindex$npar_own
  L_mis = ncol(tree_mis)
  
  # start estimation loop
  pb = txtProgressBar(max=rep*nrow(models), style=3)
  progress = function(itr) setTxtProgressBar(pb, itr)
  opts = list(progress=progress)
  itime = proc.time()[3]
  out = foreach(b=1:rep, .combine='rbind') %:% 
    foreach(m=1:nrow(models), .combine='rbind', 
            .packages=c("Rcpp","RcppArmadillo","here"),
            .options.snow=opts) %dopar% { 
              
              sourceCpp(here("src","shrinkage-mcmc.cpp"))
              
              data = DataList[[b]]

              # fit model
              if(models[m,2]=="sparse"){
                Data = list(
                  Y=data$Y,
                  X=data$X,
                  Clist=data$Clist
                )
                fit = rSURshrinkage(Data,Prior,Mcmc,Shrinkage=models[m,1],print=FALSE)
              }
              else{
                
                if(is.null(tree)){
                  Data = list(
                    Y=data$Y,
                    X=data$X,
                    Clist=data$Clist,
                    tree=tree_dgp,
                    childrencounts=childrencounts_dgp,
                    list=list_dgp,
                    list_own=list_own_dgp,
                    npar=npar_dgp,
                    npar_own=npar_own_dgp
                  )
                }
                else if(tree=="misspecified"){
                  Data = list(
                    Y=data$Y,
                    X=data$X,
                    Clist=data$Clist,
                    tree=tree_mis,
                    childrencounts=childrencounts_mis,
                    list=list_mis,
                    list_own=list_own_mis,
                    npar=npar_mis,
                    npar_own=npar_own_mis
                  )
                }
                
                fit = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product=models[m,1],group=models[m,2]),
                                        print=FALSE)
              }
              
              # compute rmse
              Brmse = double(end-burn)
              for(r in 1:(end-burn)){
                B = matrix(fit$betadraws[burn+r,],p,p)
                Brmse[r] = sqrt(mean((data$B-B)^2))
              }
              rmse = mean(Brmse)
              
              # compute share of correct signs
              B = apply(fit$betadraws[burn:end,],2,mean)
              sign = mean(sign(data$B)==sign(B))
              
              # save
              c(rmse,sign)
            }
  
  colnames(out) = c(paste0("rmse_",p,"_",dgp),paste0("signs_",p,"_",dgp))
  
  ctime = proc.time()[3]
  cat("\n Total Time:",round((ctime - itime)/60,2),"minutes",fill = TRUE)
  
  return(out)
  
}


