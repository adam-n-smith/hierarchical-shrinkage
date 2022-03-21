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
simdata_batch = function(pvec,nrep,type,prop_sparse,transform=NULL){
  
  d = 1
  
  # generate data
  data = NULL
  j = 1
  for(p in p_vec) {
    
    tree = matrix(c(rep(1:5,each=p/5),
                    rep(1:10,each=p/10),
                    1:p),nrow=p,ncol=3)
    counts = countchildren_cpp(tree)
    treeindex = createindex(tree)
    list = treeindex$list
    list_own = treeindex$list_own
    npar = treeindex$npar
    npar_own = treeindex$npar_own
    
    # tree for DGP
    if(type=="hierarchical"){
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
    
    data[j] = list(nrep)
    tmp = NULL
    for(b in 1:nrep){
      
      # simulate data
      tmp[[b]] = simdata(n,p,d,settings)
      tmp[[b]]$npar = npar
      tmp[[b]]$npar_own = npar_own
      tmp[[b]]$tree = tree
      tmp[[b]]$list = list
      tmp[[b]]$list_own = list_own
      tmp[[b]]$childrencounts = counts
      
    }
    
    data[[j]] = tmp
    j = j+1 
    
  }
  
  return(data)
  
}

# fit models on batch data in parallel
fit_parallel = function(DataList,Mcmc,p_vec,models,dgp){
  
  
  nrep = length(DataList[[1]])
  end = Mcmc$R/Mcmc$keep
  burn = Mcmc$burn_pct*end
  
  itime = proc.time()[3]
  out = foreach(j=1:length(p_vec), .combine='cbind') %:% 
    foreach(b=1:nrep, .combine='rbind') %:% 
    foreach(m=1:nrow(models), .combine='rbind', 
            .packages=c("Rcpp","RcppArmadillo","here"),
            .options.snow=opts) %dopar% { 
              
              sourceCpp(here("functions","shrinkage_mcmc.cpp"))

              data = DataList[[j]][[b]]
              p = p_vec[j]

              # Prior
              Prior = list(
                thetabar_own=0,
                thetabar_cross=0,
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
                  Clist=data$Clist
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
                  list_own=data$list_own,
                  npar=data$npar,
                  npar_own=data$npar_own
                )
                fit = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product=models[m,1],group=models[m,2]),
                                        print=FALSE)
              }

              # compute rmse
              Brmse = double(end-burn)
              for(r in 1:(end-burn)){
                B = matrix(fit$betadraws[burn+r,],p,p)
                Brmse[r] = sqrt(mean((data$B - B)^2))
              }
              rmse = mean(Brmse)
              
              # compute share of correct signs
              B = apply(fit$betadraws[burn:end,],2,mean)
              signs = mean(sign(B)==sign(data$B))

              # save
              c(rmse,signs)
            }
  
  colnames(out) = c(paste0("rmse_",p_vec,"_",dgp),paste0("signs_",p_vec,"_",dgp))
  
  ctime = proc.time()[3]
  cat("\n Total Time:",round((ctime - itime)/60,2),"minutes",fill = TRUE)
  
  return(out)
  
}


