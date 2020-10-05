library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(doSNOW)
library(here)
library(xtable)
library(snakecase)

source(here("functions","simulation_functions.R"))
source(here("functions","shrinkage_functions.R"))
sourceCpp(here("functions","shrinkage_mcmc.cpp"))

n = 50
p_vec = c(50,100)
# p_vec = 100
nrep = 25

Mcmc = list(
  R = 2000,
  initial_run = 100,
  keep = 1,
  burn_pct = 0.5
)


# ------------------------------------------------------------- #
# generate data sets
# ------------------------------------------------------------- #

set.seed(1)
data_dense_dense = simdata_batch(p_vec,nrep,"tree",prop_sparse=c(0,0,0))
data_sparse_dense = simdata_batch(p_vec,nrep,"tree",prop_sparse=c(0,0.95,0))
data_dense_sparse = simdata_batch(p_vec,nrep,"tree",prop_sparse=c(0,0,0.95))
data_sparse_sparse = simdata_batch(p_vec,nrep,"tree",prop_sparse=c(0,0.95,0.95))
data_sparse = simdata_batch(p_vec,nrep,"sparse",prop_sparse=0.95)

# bind together data
data = list(dense_dense = data_dense_dense,
            sparse_dense = data_sparse_dense,
            dense_sparse = data_dense_sparse,
            sparse_sparse = data_sparse_sparse,
            sparse = data_sparse)

# data generating processes
dgps = c("dense_dense",
         "sparse_dense",
         "dense_sparse",
         "sparse_sparse",
         "sparse")

# ------------------------------------------------------------- #
# fit models
# ------------------------------------------------------------- #

# initialize clusters
cl = makeSOCKcluster(4)
registerDoSNOW(cl)
pb = txtProgressBar(max=nrep, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

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

# fit models across all dgps
fit = NULL
for(i in 1:length(dgps)){
  fit = cbind(fit,fit_parallel(data[[i]],Mcmc,p_vec,models,dgp=dgps[i]),999)
  print(dgps[i])
}

# stop
stopCluster(cl)

# ------------------------------------------------------------- #
# print results
# ------------------------------------------------------------- #

out = fit %>%
  data.frame() %>%
  mutate(shrinkage = rep(paste0(models[,1],"/",models[,2]),nrep),
         shrinkage = factor(shrinkage,levels=unique(shrinkage))) %>%
  group_by(shrinkage) %>%
  summarise(across(everything(),mean)) %>%
  mutate(across(where(is.numeric),~ifelse(.x==999,"",.x)),
         V0 = "") %>%
  relocate(V0,.after=shrinkage)

out

print(xtable(out),include.rownames=FALSE, sanitize.text.function=identity)
