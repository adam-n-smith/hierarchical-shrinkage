library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(here)
library(doSNOW)

source(here("functions","simulation_functions.R"))
sourceCpp(here("functions","shrinkage_mcmc.cpp"))

n = 50
p_vec = c(50,100)
nrep = 25
# p_vec = 100
# nrep = 25

set.seed(1)
simdata_sparsetree = simdata_batch(p_vec,nrep,"tree",prop_sparse=0.95)

# shrinkage = c("ridge","lasso","horseshoe")
shrinkage = c("ridge","horseshoe")
Mcmc = list(
  R = 2000,
  initial_run=100,
  keep = 1,
  burn_pct = 0.5
)

# initialize clusters
cl <- makeSOCKcluster(3)
registerDoSNOW(cl)
pb <- txtProgressBar(max=nrep, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

out = fit_parallel(simdata_sparsetree,Mcmc,p_vec,shrinkage,model="tree_true",dgp="sparsetree")

stopCluster(cl)

data.frame(out) %>%
  mutate(shrinkage=rep(shrinkage,nrep)) %>%
  group_by(shrinkage) %>%
  summarise(across(everything(),mean))
