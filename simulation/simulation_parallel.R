library(Rcpp)
library(RcppArmadillo)
library(foreach)
library(doParallel)
library(abind)

source(here("functions","simulation_functions.R"))
sourceCpp(here("functions","horse_mcmc.cpp"))

n = 50
p_vec = c(50,100)
nrep = 1

models = c("ridge","horseshoe")
R = 500
keep = 1
end = R/keep
burn = 0.75*end

# ------------------------------------------------------------- #
# DGP: dense + structured, tree: correctly specified
# ------------------------------------------------------------- #

tausq0 = 0.05

# generate data
simdata = NULL
j = 1
for(p in p_vec) {
  
  tree = matrix(c(rep(1:5,each=p/5),
                  rep(1:10,each=p/10),
                  1:p),nrow=p,ncol=3)
  childrencounts = countchildren_cpp(tree)
  treeindex = createindex(tree)
  list = treeindex$list
  index = treeindex$index
  npar = unlist(lapply(list,length))
  L = ncol(tree)
  
  simdata[j]  = list(nrep)
  tmp = NULL
  for(b in 1:nrep){
    
      # simulate data
      d = 1
      tmp[[b]] = simdata_dense(n,p,d,tree,childrencounts,npar,tausq0)
      tmp[[b]]$npar = npar
      tmp[[b]]$tree = tree
      tmp[[b]]$list = list
      tmp[[b]]$childrencounts = childrencounts

  }

  simdata[[j]] = tmp
  j = j+1
  
}

# fit models
numCores = detectCores()
registerDoParallel(numCores)
out = foreach(j=1:length(p_vec), .combine='abind') %:% 
  foreach(m=models, .combine='cbind') %:% 
  foreach(b=1:nrep, .combine='rbind') %dopar% {
    
    data = simdata[[j]][[b]]
    Brmse = double(end-burn)
    p = p_vec[j]
    
    # priors
    Prior = list(thetabar=0, taubar=10, betabarii=0, taubarii=10,
                 Aphi=.01*diag(p), phibar=double(p), a=5, b=5)
    
    # mcmc
    Mcmc = list(R = R,keep = keep)
    
    # data
    Data = list(Y=data$Y, X=data$X, Clist=data$Clist, tree=data$tree,
                childrencounts=data$childrencounts, list=data$list, npar=data$npar)
    
    fit = rSURhiershrinkage(Data,Prior,Mcmc,shrinkage=m,print=FALSE)
    
    # compute rmse
    for(r in 1:(end-burn)){
      B = matrix(fit$betadraws[burn+r,],p,p)
      Brmse[r] = sqrt(mean((data$B[diag(p)!=1] - B[diag(p)!=1])^2))
    }
    
    # save 
    mean(Brmse)

  }

fit_struct_tree = apply(out,2,mean)
dim(fit_struct_tree) = c(length(models),length(p_vec))
fit_struct_tree

# ------------------------------------------------------------- #
# DGP: dense + structured, tree: misspecified
# ------------------------------------------------------------- #

tausq0 = 0.05

# generate data
simdata = NULL
j = 1
for(p in p_vec) {
  
  # tree for DGP
  tree_dgp = matrix(c(rep(1:5,each=p/5),
                      rep(1:10,each=p/10),
                      1:p),nrow=p,ncol=3)
  childrencounts_dgp = countchildren_cpp(tree)
  treeindex_dgp = createindex(tree)
  list_dgp = treeindex$list
  index_dgp = treeindex$index
  npar_dgp = unlist(lapply(list,length))
  L = ncol(tree_dgp)
  
  # tree for estimation
  tree = tree_dgp
  tree[which(tree_dgp[,2]>6),2] = c(rep(7:8,times=p/10),
                                    rep(9:10,times=p/10))
  childrencounts = countchildren_cpp(tree)
  treeindex = createindex(tree)
  list = treeindex$list
  index = treeindex$index
  npar = unlist(lapply(list,length))
  
  simdata[j]  = list(nrep)
  tmp = NULL
  for(b in 1:nrep){
    
    # simulate data
    d = 1
    tmp[[b]] = simdata_dense(n,p,d,tree_dgp,childrencounts_dgp,npar_dgp,tausq0)
    tmp[[b]]$npar = npar
    tmp[[b]]$tree = tree
    tmp[[b]]$list = list
    tmp[[b]]$childrencounts = childrencounts
    
  }
  
  simdata[[j]] = tmp
  j = j+1
  
}

# fit models
models = c("ridge","horseshoe")
R = 100
keep = 1
end = R/keep
burn = 0.75*end
numCores = detectCores()
registerDoParallel(numCores)
out = foreach(j=1:length(p_vec), .combine='abind') %:% 
  foreach(m=models, .combine='cbind') %:% 
  foreach(b=1:nrep, .combine='rbind') %dopar% {
    
    data = simdata[[j]][[b]]
    yrmse = double(end-burn)
    Brmse = double(end-burn)
    p = p_vec[j]
    
    # priors
    Prior = list(thetabar=0, taubar=10, betabarii=0, taubarii=10,
                 Aphi=.01*diag(p), phibar=double(p), a=5, b=5)
    
    # mcmc
    Mcmc = list(R = R,keep = keep)
    
    # data
    Data = list(Y=data$Y, X=data$X, Clist=data$Clist, tree=data$tree,
                childrencounts=data$childrencounts, list=data$list, npar=data$npar)
    
    fit = rSURhiershrinkage(Data,Prior,Mcmc,shrinkage=m,print=FALSE)
    
    # compute rmse
    for(r in 1:(end-burn)){
      B = matrix(fit$betadraws[burn+r,],p,p)
      Brmse[r] = sqrt(mean((data$B[diag(p)!=1] - B[diag(p)!=1])^2))
    }
    
    # save 
    mean(Brmse)
    
  }

fit_table = apply(out,2,mean)
dim(fit_table) = c(length(models),length(p_vec))
fit_table

# ------------------------------------------------------------- #
# DGP: sparse
# ------------------------------------------------------------- #

simdata = NULL

# generate data
j = 1
for(p in p_vec) {

  tree = matrix(c(rep(1:5,each=p/5),
                  rep(1:10,each=p/10),
                  1:p),nrow=p,ncol=3)
  childrencounts = countchildren_cpp(tree)
  treeindex = createindex(tree)
  list = treeindex$list
  index = treeindex$index
  npar = unlist(lapply(list,length))
  L = ncol(tree)
  
  simdata[j]  = list(nrep)
  tmp = NULL
  for(b in 1:nrep){
    
    # simulate data
    d = 1
    tmp[[b]] = simdata_sparse(n,p,d)
    tmp[[b]]$npar = npar
    tmp[[b]]$tree = tree
    tmp[[b]]$list = list
    tmp[[b]]$childrencounts = childrencounts
    
  }
  
  simdata[[j]] = tmp
  j = j+1
  
}


# fit models
models = c("ridge","horseshoe")
R = 100
keep = 1
end = R/keep
burn = 0.75*end
numCores = detectCores()
registerDoParallel(numCores)
out = foreach(j=1:length(p_vec), .combine='abind') %:% 
  foreach(m=models, .combine='cbind') %:% 
  foreach(b=1:nrep, .combine='rbind') %dopar% {
    
    data = simdata[[j]][[b]]
    yrmse = double(end-burn)
    Brmse = double(end-burn)
    p = p_vec[j]
    
    # priors
    Prior = list(thetabar=0, taubar=10, betabarii=0, taubarii=10,
                 Aphi=.01*diag(p), phibar=double(p), a=5, b=5)
    
    # mcmc
    Mcmc = list(R = R,keep = keep)
    
    # data
    Data = list(Y=data$Y, X=data$X, Clist=data$Clist, tree=data$tree,
                childrencounts=data$childrencounts, list=data$list, npar=data$npar)
    
    fit = rSURhiershrinkage(Data,Prior,Mcmc,shrinkage=m,print=FALSE)
    
    # compute rmse
    for(r in 1:(end-burn)){
      B = matrix(fit$betadraws[burn+r,],p,p)
      Brmse[r] = sqrt(mean((data$B[diag(p)!=1] - B[diag(p)!=1])^2))
    }
    
    # save 
    mean(Brmse)
    
  }

fit_table = apply(out,2,mean)
dim(fit_table) = c(length(models),length(p_vec))
fit_table

stopCluster(cl)
