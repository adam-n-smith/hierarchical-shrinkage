library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(foreach)
library(doParallel)
library(doSNOW)
library(here)

sourceCpp(here("src","shrinkage-mcmc.cpp"))
source(here("analysis","empirical","estimation-build.R"))

# where to save output
folder = "output"

# --------------------------------------------------------- #
# functions
# --------------------------------------------------------- #

# print start time
print_start = function(){
  cat("Start time -",format(Sys.time(), "%I:%M%p"), "\n")
}

# wrapper
fit_parallel = function(Data,Prior,Mcmc,beta,theta="sparse",path){
  
  # start time
  itime = proc.time()[3]
  
  # run
  if(theta=="sparse"){
    out = rSURshrinkage(Data,Prior,Mcmc,Shrinkage=beta,print=FALSE)
  }
  else{
    out = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product=beta,group=theta),print=FALSE)
  }
  
  # save
  outname = paste("out",theta,beta,sep=".")
  assign(outname, out)
  save(list=outname, file=paste0(path,"/out-",theta,"-",beta,".RData"))
  
  # print
  return(paste("Total time -",round((proc.time()[3]-itime)/60,2),"minutes"))
  
}

# --------------------------------------------------------- #
# run
# --------------------------------------------------------- #

# data
Data = list(
  Y = Y,
  X = X,
  Clist = Clist,
  tree = tree,
  parindextree = parindextree,
  parindextree_own = parindextree_own,
  npar = npar,
  npar_own = npar_own
)

# priors
Prior = list(
  thetabar_cross = 0,
  thetabar_own = 0,
  Aphi = .01*diag(nphi),
  phibar = double(nphi),
  a = 5,
  b = 5
)

# mcmc
Mcmc = list(
  R = 100000,
  keep = 100
)

# models
Models = matrix(c("ridge","sparse",
                  "lasso","sparse",
                  "horseshoe","sparse",
                  "ridge","ridge",
                  "lasso","ridge",
                  "horseshoe","ridge",
                  "ridge","horseshoe",
                  "lasso","horseshoe",
                  "horseshoe","horseshoe"), ncol=2, byrow=TRUE)

# run
cl = makeSOCKcluster(9)
registerDoSNOW(cl)
print_start()
foreach(i=1:nrow(Models), .packages=c("Rcpp","RcppArmadillo","here"), .noexport=c("rSURshrinkage","rSURhiershrinkage")) %dopar% {
  sourceCpp(here("src","shrinkage-mcmc.cpp"))
  fit_parallel(Data,Prior,Mcmc,beta=Models[i,1],theta=Models[i,2],path=here("analysis",folder))
}
stopCluster(cl)
