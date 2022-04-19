drawtree = function(tree){
  out = tree
  n = nrow(tree)
  L = ncol(tree)
  out[,1] = 1:n
  flip = sample(0:1,1,)
  mean = flip*tree[,2] + (1-flip)*(1:n)
  out[,2] = rlsp(1,mean,runif(1,0,1))
  for(ell in 3:L){
    K = max(out[,ell-1])
    draw = rlsp(1,1:K,10)
    out[,ell] = plyr::mapvalues(out[,ell-1],1:K,draw)
  }
  return(out)
}


set.seed(1)
n = 50
p = 50
rep = 1
data = list(
  dense.dense = simdata_batch(n,p,rep,"hierarchical",prop_sparse=c(0,0,0),transform=FALSE),
  dense.transform = simdata_batch(n,p,rep,"hierarchical",prop_sparse=c(0.75,0,0),transform=TRUE),
  sparse.sparse = simdata_batch(n,p,rep,"hierarchical",prop_sparse=c(0,0.75,0),transform=FALSE),
  sparse = simdata_batch(n,p,rep,"sparse",prop_sparse=0.95)
)


# Prior
Prior = list(
  thetabar_own = 0, 
  thetabar_cross = 0, 
  Aphi = .01*diag(p), 
  phibar = double(p),
  a = 5, 
  b = 5
)

# Mcmc
Mcmc = list(
  R = 1000,
  keep = 1,
  burn_pct = 0.5
)

# models to fit
models = matrix(c("ridge","ridge",
                  "lasso","ridge",
                  "horseshoe","ridge",
                  "ridge","horseshoe",
                  "lasso","horseshoe",
                  "horseshoe","horseshoe"),
                ncol=2,byrow=TRUE)

TreeList = list(10)
tree_true = data$dense.dense[[1]]$tree
for(i in 1:10){
  TreeList[[i]] = drawtree(tree_true)
}

# fit models across all dgps
fitmis1 = fit_parallel_misspecified(data$dense.dense,TreeList,Prior,Mcmc,models,dgp="dense.dense")
fitmis2 = fit_parallel_misspecified(data$dense.transform,TreeList,Prior,Mcmc,models,dgp="dense.transform")
fitmis3 = fit_parallel_misspecified(data$sparse.sparse,TreeList,Prior,Mcmc,models,dgp="sparse.sparse")
fitmis = cbind(fitmis1,fitmis2,fitmis3)

data.frame(fitmis1) %>%
  # add model labels
  mutate(shrinkage = rep(paste0(to_any_case(models[,2],case="upper_camel"),"/",
                                to_any_case(models[,1],case="upper_camel")),rep),
         shrinkage = factor(shrinkage,levels=unique(shrinkage)),
         rep = rep(1:(nrow(fitmis1)/nrow(models)),each=nrow(models))) %>%
  # compute means across data replicates
  # group_by(shrinkage) %>%
  # summarise(across(everything(),mean)) %>%
  # restructure
  pivot_longer(-c(shrinkage,rep),names_pattern = "(.*)_(.*)_(.*)",names_to=c("variable","p","dgp")) %>%
  pivot_wider(names_from=c(variable,p,dgp),values_from=value) %>%
  ggplot(aes(x=sim_50_dense.dense,y=rmse_50_dense.dense,color=shrinkage)) +
  geom_point() +
  geom_smooth(method="lm")
