library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(doSNOW)
library(here)
library(xtable)
library(snakecase)

source(here("simulation","simulation_functions.R"))
source(here("functions","shrinkage_functions.R"))
sourceCpp(here("functions","shrinkage_mcmc.cpp"))

n = 50
# p_vec = c(50,100,150)
# nrep = 25
# p_vec = c(50,100)
# nrep = 10
p_vec = 50
nrep = 10

Mcmc = list(
  R = 1000,
  initial_run = 100,
  keep = 1,
  burn_pct = 0.5
)


# ------------------------------------------------------------- #
# generate data sets
# ------------------------------------------------------------- #

set.seed(1)
data_dense_dense = simdata_batch(p_vec,nrep,"hierarchical",prop_sparse=c(0,0,0),transform=FALSE)
data_sparse_dense = simdata_batch(p_vec,nrep,"hierarchical",prop_sparse=c(0,0,0.75),transform=TRUE)
data_sparse_sparse = simdata_batch(p_vec,nrep,"hierarchical",prop_sparse=c(0,0,0.75),transform=FALSE)
data_sparse = simdata_batch(p_vec,nrep,"sparse",prop_sparse=0.75)

# bind together data
data = list(dense_dense = data_dense_dense,
            sparse_dense = data_sparse_dense,
            sparse_sparse = data_sparse_sparse,
            sparse = data_sparse)

# data generating processes (beta_theta)
# dgps = c("dense_dense",
#          "sparse_dense",
#          "sparse_sparse",
#          "sparse")
dgps = "sparse_sparse"

# ------------------------------------------------------------- #
# fit models
# ------------------------------------------------------------- #

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
# models = matrix(c("ridge","sparse",
#                   "horseshoe","sparse",
#                   "horseshoe","ridge",
#                   "ridge","horseshoe",
#                   "horseshoe","horseshoe"),
#                 ncol=2,byrow=TRUE)
# 
# models = matrix(c("ridge","sparse",
#                   "horseshoe","sparse"),
#                 ncol=2,byrow=TRUE)
# 
# models = matrix(c("horseshoe","sparse",
#                   "horseshoe","ridge",
#                   "horseshoe","horseshoe"),
#                 ncol=2,byrow=TRUE)

# initialize clusters
cl = makeSOCKcluster(4)
registerDoSNOW(cl)
pb = txtProgressBar(max=nrep, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

# fit models across all dgps
fit = NULL
for(i in 1:length(dgps)){
  print(dgps[i])
  tmp = fit_parallel(data[[which(names(data)==dgps[i])]],Mcmc,p_vec,models,dgp=dgps[i])
  fit = cbind(fit,tmp,999)
}

# stop
stopCluster(cl)

# ------------------------------------------------------------- #
# print results
# ------------------------------------------------------------- #

out = data.frame(fit[,-ncol(fit)]) %>%
  mutate(shrinkage = rep(paste0(to_any_case(models[,2],case="upper_camel"),"/",
                                to_any_case(models[,1],case="upper_camel")),nrep),
         shrinkage = factor(shrinkage,levels=unique(shrinkage))) %>%
  group_by(shrinkage) %>%
  summarise(across(everything(),mean)) %>%
  mutate(across(where(is.numeric),~ifelse(.x==999,"",.x)),
         V0 = "") %>%
  relocate(V0,.after=shrinkage) %>%
  mutate(shrinkage = str_replace(shrinkage,"/Sparse",""),
         shrinkage = str_replace(shrinkage,"Horseshoe","HS"),
         shrinkage = str_replace(shrinkage,"/Horseshoe","/HS"))

out

print(xtable(out, digits=3),include.rownames=FALSE, sanitize.text.function=identity)
