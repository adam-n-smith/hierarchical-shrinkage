library(extraDistr)

# simulate data
simdata = function(n,p,d,settings){
  
  # observational error
  sigmasq = rep(1,p)
  Sigma = diag(sigmasq)
  
  thetalist = list()
  lambdalist = list()
  
  if(settings$type=="hierarchical"){
    
    tree = settings$tree
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
      
      # draw local variances for current level
      lambda = rep(1,npar[ell])

      # product of previous lambda parameters
      # Psi = create_Psi_ellmone_cpp(lambdalist,counts,list,npar,ell)
      Psi = replace_cpp(parindextree[[ell]],Psi)*lambda*sparse_ind
      
      # mean as theta from last level
      mean = theta[parindextree[[ell]]]
      
      # construct theta for current level
      if(settings$transform){
        theta = mean + sqrt(Psi*tau[ell])*rnorm(npar[ell])*sparse_ind
      }
      else{
        theta = (mean + sqrt(Psi*tau[ell])*rnorm(npar[ell]))*sparse_ind
      }
      
      # save
      thetalist[[ell]] = theta
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
  Bres = matrix(NA,p,p)
  Bres[!(diag(p)==1)] = theta - mean
  Bbar = matrix(0,p,p)
  Bbar[!(diag(p)==1)] = mean
  
  return(list(Y=Y,X=X,B=B,Clist=Clist,
              thetalist=thetalist,tau=tau,sigmasq=sigmasq,
              phivec=phivec,nphi=nphi,Bres=Bres,Bbar=Bbar))
  
}

# simulate a batch of data sets to estimate in parallel
simdata_batch = function(n,p,rep,type,prop_sparse,transform=NULL){
  
  d = 1
  
  # generate data
  data = NULL
  
  # tree for DGP
  if(type=="hierarchical"){
    
    tree = matrix(c(1:p,
                    rep(1:10,each=p/10),
                    rep(1:5,each=p/5)),nrow=p,ncol=3)
    objects = create_treeobjects(tree)
    parindextree = objects$parindextree
    parindextree_own = objects$parindextree_own
    npar = objects$npar
    npar_own = objects$npar_own
    
    settings = list(
      type = type,
      tree = tree,
      npar = npar,
      npar_own = npar_own,
      parindextree = parindextree,
      parindextree_own = parindextree_own,
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
      
      data[[b]]$tree = tree
      data[[b]]$npar = npar
      data[[b]]$npar_own = npar_own
      data[[b]]$parindextree = parindextree
      data[[b]]$parindextree_own = parindextree_own

    }
    
  }
  
  return(data)
  
}

# fit models on batch data in parallel
fit_parallel = function(DataList,Prior,Mcmc,models,dgp,tree=NULL){
  
  rep = length(DataList)
  end = Mcmc$R/Mcmc$keep
  burn = Mcmc$burn_pct*end
  
  # tree
  if(is.null(tree)){
    tree = DataList[[1]]$tree
  }
  objects = create_treeobjects(tree)
  parindextree = objects$parindextree
  parindextree_own = objects$parindextree_own
  npar = objects$npar
  npar_own = objects$npar_own
  L = ncol(tree)
  p = nrow(tree)
  
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
              sign = mean(sign(data$B[data$B!=0])==sign(B[as.vector(data$B!=0)]))

              # save
              c(rmse,sign)
            }
  
  colnames(out) = c(paste0("rmse_",p,"_",dgp),paste0("signs_",p,"_",dgp))
  
  ctime = proc.time()[3]
  cat("\n Total Time:",round((ctime - itime)/60,2),"minutes",fill = TRUE)
  
  return(out)
  
}