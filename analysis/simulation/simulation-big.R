library(tidyverse)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(doSNOW)
library(here)
library(xtable)
library(snakecase)

source(here("src","simulation-functions.R"))
source(here("src","shrinkage-functions.R"))
sourceCpp(here("src","shrinkage-mcmc.cpp"))

# generate data
set.seed(1)
n = 100
p = 300
rep = 1
data_big = list(
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
  R = 40000,
  keep = 40,
  burn_pct = 0.75
)

# models to fit
models = matrix(c("ridge","sparse",
                  "lasso","sparse",
                  "horseshoe","sparse",
                  "ridge","ridge",
                  "lasso","ridge",
                  "horseshoe","ridge",
                  "ridge","horseshoe",
                  "lasso","horseshoe",
                  "horseshoe","horseshoe"),
                ncol=2,byrow=TRUE)

# initialize clusters
ncores = min(parallel::detectCores(),15)
cl = makeSOCKcluster(ncores)
registerDoSNOW(cl)

# fit models across all dgps
fitbig1 = fit_parallel(data_big$dense.dense,Prior,Mcmc,models,dgp="dense.dense")
fitbig2 = fit_parallel(data_big$dense.transform,Prior,Mcmc,models,dgp="dense.transform")
fitbig3 = fit_parallel(data_big$sparse.sparse,Prior,Mcmc,models,dgp="sparse.sparse")
fitbig4 = fit_parallel(data_big$sparse,Prior,Mcmc,models,dgp="sparse",tree=data_big$dense.dense[[1]]$tree)
fitbig = cbind(fitbig1,fitbig2,fitbig3,fitbig4)

# print results
out_big = data.frame(fitbig) %>%
  # add model labels
  mutate(shrinkage = rep(paste0(to_any_case(models[,2],case="upper_camel"),"/",
                                to_any_case(models[,1],case="upper_camel")),rep),
         shrinkage = factor(shrinkage,levels=unique(shrinkage))) %>%
  # compute means across data replicates
  group_by(shrinkage) %>%
  summarise(across(everything(),mean)) %>%
  # restructure
  pivot_longer(-shrinkage,names_pattern = "(.*)_(.*)_(.*)",names_to=c("variable","p","dgp")) %>%
  pivot_wider(names_from=c(variable,p,dgp),values_from=value) %>%
  # reorder columns
  select(shrinkage,starts_with("rmse"),starts_with("signs")) %>%
  # add bold to 
  mutate(across(starts_with("rmse_"),function(x)ifelse(x==min(x),paste("\\bf",sprintf("%.3f",x)),sprintf("%.3f",x))))

out_big
print(xtable(out_big),include.rownames=FALSE, sanitize.text.function=identity)
save(n,p,rep,data_big,models,Prior,Mcmc,fitbig,fitbig1,fitbig2,fitbig3,fitbig4,out_big,
     file="analysis/simulation/output/sim-big.RData")

# stop
stopCluster(cl)

