library(tidyverse)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(doSNOW)
library(here)
library(xtable)
library(snakecase)
library(mclust)

source(here("src","simulation-functions.R"))
source(here("src","shrinkage-functions.R"))
sourceCpp(here("src","shrinkage-mcmc.cpp"))

# ------------------------------------------------------------- #
# generate data sets (n=50, p=100)
# ------------------------------------------------------------- #

set.seed(1)
n = 50
p = 100
rep = 25
data = list(
  dense.dense = simdata_batch(n,p,rep,"hierarchical",prop_sparse=c(0,0,0),transform=FALSE),
  dense.transform = simdata_batch(n,p,rep,"hierarchical",prop_sparse=c(0.75,0,0),transform=TRUE),
  sparse.sparse = simdata_batch(n,p,rep,"hierarchical",prop_sparse=c(0,0.75,0),transform=FALSE)
)
tree_true = data$dense.dense[[1]]$tree

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
  R = 5000,
  keep = 5,
  burn_pct = 0.75
)

# initialize clusters
ncores = min(parallel::detectCores(),15)
cl = makeSOCKcluster(ncores)
registerDoSNOW(cl)

# ------------------------------------------------------------- #
# misspecification A
# ------------------------------------------------------------- #

# misspecified tree (no common structure)
tree_mis = tree_true
for(k in 2:ncol(tree_mis)){
  K = max(tree_mis[,k])
  tree_mis[,k] = rep(1:K,times=p/K)
}
adjustedRandIndex(tree_true[,-1],tree_mis[,-1])

# fit models across all dgps
fitmis1a = fit_parallel(data$dense.dense,Prior,Mcmc,models,dgp="dense.dense",tree=tree_mis)
fitmis2a = fit_parallel(data$dense.transform,Prior,Mcmc,models,dgp="dense.transform",tree=tree_mis)
fitmis3a = fit_parallel(data$sparse.sparse,Prior,Mcmc,models,dgp="sparse.sparse",tree=tree_mis)
fitmisa = cbind(fitmis1a,fitmis2a,fitmis3a)

# print results
out_misa = data.frame(fitmisa) %>%
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

out_misa
print(xtable(out_misa),include.rownames=FALSE, sanitize.text.function=identity)
save(n,p,rep,data,models,Prior,Mcmc,fitmisa,fitmis1a,fitmis2a,fitmis3a,out_misa,tree_mis,tree_true,
     file="analysis/simulation/output/sim-misspecified-a.RData")

# stop
# stopCluster(cl)

# ------------------------------------------------------------- #
# misspecification B
# ------------------------------------------------------------- #

# misspecified tree (half has no common structure)
tree_mis = tree_true
wchrow = max(which(tree_mis[,2]==4))
tree_mis[1:wchrow,2] = rep(1:4,times=p/10)
tree_mis[1:wchrow,3] = rep(1:2,times=p/5)
adjustedRandIndex(tree_true[,-1],tree_mis[,-1])

# fit models across all dgps
fitmis1b = fit_parallel(data$dense.dense,Prior,Mcmc,models,dgp="dense.dense",tree=tree_mis)
fitmis2b = fit_parallel(data$dense.transform,Prior,Mcmc,models,dgp="dense.transform",tree=tree_mis)
fitmis3b = fit_parallel(data$sparse.sparse,Prior,Mcmc,models,dgp="sparse.sparse",tree=tree_mis)
fitmisb = cbind(fitmis1b,fitmis2b,fitmis3b)

# print results
out_misb = data.frame(fitmisb) %>%
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

out_misb
print(xtable(out_misb),include.rownames=FALSE, sanitize.text.function=identity)
save(n,p,rep,data,models,Prior,Mcmc,fitmisb,fitmis1b,fitmis2b,fitmis3b,out_misb,tree_mis,tree_true,
     file="analysis/simulation/output/sim-misspecified-b.RData")

# stop
# stopCluster(cl)

# ------------------------------------------------------------- #
# misspecification C
# ------------------------------------------------------------- #

# misspecified tree (right structure but more groups than necessary)
tree_mis = tree_true
wchrow = max(which(tree_mis[,2]==4))
tree_mis[1:wchrow,2] = rep(1:8,each=5)
tree_mis[(wchrow+1):p,2] = tree_mis[(wchrow+1):p,2] + 4
adjustedRandIndex(tree_true[,-1],tree_mis[,-1])

# fit models across all dgps
fitmis1c = fit_parallel(data$dense.dense,Prior,Mcmc,models,dgp="dense.dense",tree=tree_mis)
fitmis2c = fit_parallel(data$dense.transform,Prior,Mcmc,models,dgp="dense.transform",tree=tree_mis)
fitmis3c = fit_parallel(data$sparse.sparse,Prior,Mcmc,models,dgp="sparse.sparse",tree=tree_mis)
fitmisc = cbind(fitmis1c,fitmis2c,fitmis3c)

# print results
out_misc = data.frame(fitmisc) %>%
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

out_misc
print(xtable(out_misc),include.rownames=FALSE, sanitize.text.function=identity)
save(n,p,rep,data,models,Prior,Mcmc,fitmisc,fitmis1c,fitmis2c,fitmis3c,out_misc,tree_mis,tree_true,
     file="analysis/simulation/output/sim-misspecified-c.RData")

# stop
stopCluster(cl)
