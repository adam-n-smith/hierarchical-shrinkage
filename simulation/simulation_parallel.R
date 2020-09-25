library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(doSNOW)
library(here)
library(xtable)
library(snakecase)

# 30 minutes for p=100, nrep=2, R=1000
# 50 minutes for p=100, nrep=4, R=1000

source(here("functions","simulation_functions.R"))
sourceCpp(here("functions","shrinkage_mcmc.cpp"))

n = 50
p_vec = c(50,100)
nrep = 25
# p_vec = 100
# nrep = 10

# shrinkage = c("ridge","lasso","horseshoe")
shrinkage = c("ridge","horseshoe")
Mcmc = list(
  R = 2000,
  initial_run=100,
  keep = 1,
  burn_pct = 0.5
)


# ------------------------------------------------------------- #
# generate data sets
# ------------------------------------------------------------- #

set.seed(1)
simdata_densetree = simdata_batch(p_vec,nrep,"tree",prop_sparse=0)
simdata_sparsetree = simdata_batch(p_vec,nrep,"tree",prop_sparse=0.95)
simdata_sparse = simdata_batch(p_vec,nrep,"sparse",prop_sparse=0.95)

# ------------------------------------------------------------- #
# fit models
# ------------------------------------------------------------- #

# storage objects
out = vector(mode="list",length=3)
names(out) = c("sparse","tree_true","tree_not")

# initialize clusters
cl = makeSOCKcluster(3)
registerDoSNOW(cl)
pb = txtProgressBar(max=nrep, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)


# DGP: dense tree
models = c("sparse","tree_true","tree_not")
for(i in 1:length(models)){
  print(models[i])
  out[[i]] = fit_parallel(simdata_densetree,Mcmc,p_vec,shrinkage,model=models[i],dgp="densetree")
}

# DGP: sparse tree
models = c("sparse","tree_true","tree_not")
for(i in 1:length(models)){
  print(models[i])
  tmp = fit_parallel(simdata_sparsetree,Mcmc,p_vec,shrinkage,model=models[i],dgp="sparsetree")
  out[[i]] = cbind(out[[i]],99,tmp)
}

# DGP: sparse
models = c("sparse","tree_true")
for(i in 1:length(models)){
  print(models[i])
  tmp = fit_parallel(simdata_sparse,Mcmc,p_vec,shrinkage,model=models[i],dgp="sparse")
  out[[i]] = cbind(out[[i]],99,tmp)
}

# stop
stopCluster(cl)

# ------------------------------------------------------------- #
# merge results
# ------------------------------------------------------------- #

df_sparse = out %>%
  pluck("sparse") %>%
  data.frame() %>%
  mutate(shrinkage = rep(shrinkage,nrep)) %>%
  group_by(shrinkage) %>%
  summarise(across(everything(),mean)) %>%
  mutate(across(everything(),~ifelse(.x==99,"",.x)),
         space = "") %>%
  relocate(space,.after=shrinkage) %>%
  arrange(desc(shrinkage)) %>%
  mutate(shrinkage = paste("\\qquad",to_any_case(shrinkage, case="upper_camel")))
df_treetrue = out %>%
  pluck("tree_true") %>%
  data.frame() %>%
  mutate(shrinkage = rep(shrinkage,nrep)) %>%
  group_by(shrinkage) %>%
  summarise(across(everything(),mean)) %>%
  mutate(across(everything(),~ifelse(.x==99,"",.x)),
         space = "") %>%
  relocate(space,.after=shrinkage) %>%
  arrange(desc(shrinkage))%>%
  mutate(shrinkage = paste("\\qquad",to_any_case(shrinkage, case="upper_camel")))
df_treenot = out %>%
  pluck("tree_not") %>%
  data.frame() %>%
  mutate(shrinkage = rep(shrinkage,nrep)) %>%
  group_by(shrinkage) %>%
  summarise(across(everything(),mean)) %>%
  mutate(across(everything(),~ifelse(.x==99,"",.x)),
         space = "") %>%
  relocate(space,.after=shrinkage) %>%
  arrange(desc(shrinkage)) %>%
  mutate(shrinkage = paste("\\qquad",to_any_case(shrinkage, case="upper_camel")))

print(xtable(df_sparse),include.rownames=FALSE, sanitize.text.function=identity)
print(xtable(df_treetrue),include.rownames=FALSE, sanitize.text.function=identity)
print(xtable(df_treenot),include.rownames=FALSE, sanitize.text.function=identity)
