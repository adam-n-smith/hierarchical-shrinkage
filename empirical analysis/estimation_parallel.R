library(foreach)
library(doParallel)

# --------------------------------------------------------- #
# functions
# --------------------------------------------------------- #

# print start time
print_start = function(){
  cat("Start time -",format(Sys.time(), "%I:%M%p"), "\n")
}

# wrapper
fit_parallel = function(Data,Prior,Mcmc,beta,theta="sparse",folder){
  
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
  save(list=outname, file=here("empirical analysis",folder,paste0(theta,"_",beta,".RData")))
  
  # print
  return(paste("Total time -",round((proc.time()[3]-itime)/60,2),"minutes"))

}

# --------------------------------------------------------- #
# run
# --------------------------------------------------------- #

# models
beta_priors = c("ridge","lasso","horseshoe")

# mcmc
Mcmc = list(
  R = 20000,
  initial_run = 100,
  keep = 10
)
# Mcmc = list(
#   R = 10000,
#   initial_run = 100,
#   keep = 5
# )

# where to save output
folder = "tmp"

# --------------------------------------------------------- #
# sparse models
# --------------------------------------------------------- #

# start
registerDoParallel(3)

# data
Data = list(
  Y = Y,
  X = X,
  Clist = Clist
)

# priors
Prior = list(
  Aphi = .1*diag(nphi),
  phibar = double(nphi),
  a = 5,
  b = 5
)

# sparse
print_start()
foreach (i=1:length(beta_priors)) %dopar% {
  fit_parallel(Data,Prior,Mcmc,beta=beta_priors[i],folder=folder)
}

# stop
stopImplicitCluster()

# --------------------------------------------------------- #
# hierarchical models
# --------------------------------------------------------- #

# start
registerDoParallel(3)
folder = "tmp"

# data
Data = list(
  Y = Y,
  X = X,
  Clist = Clist,
  tree = tree,
  childrencounts = childrencounts,
  list = list,
  list_own = list_own,
  npar = npar,
  npar_own = npar_own
)

# priors
Prior = list(
  thetabar_cross = 0,
  thetabar_own = 0,
  Aphi = .1*diag(nphi),
  phibar = double(nphi),
  a = 5,
  b = 5
)

# hierarchical (ridge)
print_start()
foreach (i=1:length(beta_priors)) %dopar% {
  fit_parallel(Data,Prior,Mcmc,beta=beta_priors[i],theta="ridge",folder=folder)
}

# hierarchical (horseshoe)
print_start()
foreach (i=1:length(beta_priors)) %dopar% {
  fit_parallel(Data,Prior,Mcmc,beta=beta_priors[i],theta="horseshoe",folder)
}

# stop
stopImplicitCluster()
