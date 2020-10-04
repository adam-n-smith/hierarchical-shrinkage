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
# p_vec = c(50,100)
# nrep = 25
p_vec = c(50,60)
nrep = 2

Mcmc = list(
  R = 200,
  initial_run=0,
  keep = 1,
  burn_pct = 0.5
)


# ------------------------------------------------------------- #
# generate data sets
# ------------------------------------------------------------- #

set.seed(1)
# simdata_densetree = simdata_batch(p_vec,nrep,"tree",prop_sparse=0)
# simdata_sparsetree = simdata_batch(p_vec,nrep,"tree",prop_sparse=0.95)
# simdata_sparse = simdata_batch(p_vec,nrep,"sparse",prop_sparse=0.95)
data_dense_dense = simdata_batch(p_vec,nrep,"tree",prop_sparse=c(0,0,0))
data_sparse_dense = simdata_batch(p_vec,nrep,"tree",prop_sparse=c(0,0.95,0))
data_dense_sparse = simdata_batch(p_vec,nrep,"tree",prop_sparse=c(0,0,0.95))
data_sparse_sparse = simdata_batch(p_vec,nrep,"tree",prop_sparse=c(0,0.95,0.95))
data_sparse = simdata_batch(p_vec,nrep,"sparse",prop_sparse=0.95)


# ------------------------------------------------------------- #
# fit models
# ------------------------------------------------------------- #

# storage objects
out = vector(mode="list",length=3)
names(out) = c("sparse","tree_true","tree_not")

# initialize clusters
cl = makeSOCKcluster(2)
registerDoSNOW(cl)
pb = txtProgressBar(max=nrep, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

# data generating processes
dgps = c("dense_dense",
         "sparse_dense",
         "dense_sparse",
         "sparse_sparse",
         "sparse")

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
  fit = cbind(fit,fit_parallel(data_dense_dense,Mcmc,p_vec,models,dgp=dgps[i]))
  print(dgps[i])
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
