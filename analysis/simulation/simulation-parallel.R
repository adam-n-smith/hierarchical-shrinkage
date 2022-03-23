library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(doSNOW)
library(here)
library(xtable)
library(snakecase)

source(here("src","simulation-functions.R"))
source(here("src","shrinkage-functions.R"))
sourceCpp(here("src","shrinkage-mcmc.cpp"))

# ------------------------------------------------------------- #
# generate data sets
# ------------------------------------------------------------- #

set.seed(1)
n = 50
p = 100
rep = 25
data = list(
  dense.dense = simdata_batch(n,p,rep,"hierarchical",prop_sparse=c(0,0,0),transform=FALSE),
  dense.transform = simdata_batch(n,p,rep,"hierarchical",prop_sparse=c(0,0,0.75),transform=TRUE),
  sparse.sparse = simdata_batch(n,p,rep,"hierarchical",prop_sparse=c(0,0.75,0),transform=FALSE),
  sparse = simdata_batch(n,p,rep,"sparse",prop_sparse=0.75)
)

# ------------------------------------------------------------- #
# Simulation 1: n=50, p=100, rep=25
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

# initialize clusters
cl = makeSOCKcluster(6)
registerDoSNOW(cl)

# Prior
Prior = list(
  thetabar_own = 0, 
  thetabar_cross = 0, 
  taubarii = 10,
  Aphi = .01*diag(p), 
  phibar = double(p),
  a = 5, 
  b = 5
)

# Mcmc
Mcmc = list(
  R = 1000,
  initial_run = 100,
  keep = 1,
  burn_pct = 0.75
)

# fit models across all dgps
dgps = c("dense.dense","dense.transform","sparse.sparse","sparse")
fit = NULL
for(dgp in dgps){
  print(dgp)
  tmp = fit_parallel(data[[which(names(data)==dgp)]],Prior,Mcmc,p,models,dgp=dgp)
  fit = cbind(fit,tmp)
}

# stop
stopCluster(cl)

# print results
out = data.frame(fit) %>%
  # add model labels
  mutate(shrinkage = rep(paste0(to_any_case(models[,2],case="upper_camel"),"/",
                                to_any_case(models[,1],case="upper_camel")),rep),
         shrinkage = factor(shrinkage,levels=unique(shrinkage))) %>%
  # compute means across data replicates
  group_by(shrinkage) %>%
  summarise(across(everything(),mean)) %>%
  # restructure
  pivot_longer(-shrinkage,names_pattern = "(.*)_(.*)_(.*)",names_to=c("variable","p","dgp")) %>%
  filter(!(variable=="signs" & str_detect(dgp,"sparse"))) %>%
  pivot_wider(names_from=c(variable,p,dgp),values_from=value) %>%
  # reorder columns
  select(shrinkage,starts_with("rmse"),starts_with("signs")) %>%
  # add bold to 
  mutate(across(starts_with("rmse_"),function(x)ifelse(x==min(x),paste("\\bf",sprintf("%.3f",x)),sprintf("%.3f",x))))

out

print(xtable(out),include.rownames=FALSE, sanitize.text.function=identity)

# ------------------------------------------------------------- #
# Simulation 2: n=50, p=100, rep=25, misspecified tree
# ------------------------------------------------------------- #

# fit models across all dgps
fit_mis = NULL
dgps = c("dense.dense","dense.transform","sparse.sparse")
for(dgp in dgps){
  print(dgp)
  tmp = fit_parallel(data[[which(names(data)==dgp)]],Prior,Mcmc,p,models,dgp=dgp,tree="misspecified")
  fit_mis = cbind(fit_mis,tmp)
}

# print results
out_mis = data.frame(fit_mis) %>%
  # add model labels
  mutate(shrinkage = rep(paste0(to_any_case(models[,2],case="upper_camel"),"/",
                                to_any_case(models[,1],case="upper_camel")),rep),
         shrinkage = factor(shrinkage,levels=unique(shrinkage))) %>%
  # compute means across data replicates
  group_by(shrinkage) %>%
  summarise(across(everything(),mean)) %>%
  # restructure
  pivot_longer(-shrinkage,names_pattern = "(.*)_(.*)_(.*)",names_to=c("variable","p","dgp")) %>%
  filter(!(variable=="signs" & str_detect(dgp,"sparse"))) %>%
  pivot_wider(names_from=c(variable,p,dgp),values_from=value) %>%
  # reorder columns
  select(shrinkage,starts_with("rmse"),starts_with("signs")) %>%
  # add bold to 
  mutate(across(starts_with("rmse_"),function(x)ifelse(x==min(x),paste("\\bf",sprintf("%.3f",x)),sprintf("%.3f",x))))

out_mis

# ------------------------------------------------------------- #
# Simulation 3: n=100, p=300, rep=1
# ------------------------------------------------------------- #

# generate data
set.seed(1)
n = 100
p = 300
rep = 1
data_big = list(
  dense.dense = simdata_batch(n,p,rep,"hierarchical",prop_sparse=c(0,0,0),transform=FALSE),
  dense.transform = simdata_batch(n,p,rep,"hierarchical",prop_sparse=c(0,0,0.75),transform=TRUE),
  sparse.sparse = simdata_batch(n,p,rep,"hierarchical",prop_sparse=c(0,0.75,0),transform=FALSE),
  sparse = simdata_batch(n,p,rep,"sparse",prop_sparse=0.75)
)

# Prior
Prior = list(
  thetabar_own = 0, 
  thetabar_cross = 0, 
  taubarii = 10,
  Aphi = .01*diag(p), 
  phibar = double(p),
  a = 5, 
  b = 5
)

# Mcmc
Mcmc = list(
  R = 5000,
  initial_run = 100,
  keep = 1,
  burn_pct = 0.75
)

# fit models across all dgps
fit_big = NULL
dgps = c("dense.dense","dense.transform","sparse.sparse")
for(dgp in dgps){
  print(dgp)
  tmp = fit_parallel(data[[which(names(data)==dgp)]],Prior,Mcmc,p,models,dgp=dgp,tree="misspecified")
  fit_big = cbind(fit_big,tmp)
}

# print results
out_big = data.frame(fit_big) %>%
  # add model labels
  mutate(shrinkage = rep(paste0(to_any_case(models[,2],case="upper_camel"),"/",
                                to_any_case(models[,1],case="upper_camel")),rep),
         shrinkage = factor(shrinkage,levels=unique(shrinkage))) %>%
  # compute means across data replicates
  group_by(shrinkage) %>%
  summarise(across(everything(),mean)) %>%
  # restructure
  pivot_longer(-shrinkage,names_pattern = "(.*)_(.*)_(.*)",names_to=c("variable","p","dgp")) %>%
  filter(!(variable=="signs" & str_detect(dgp,"sparse"))) %>%
  pivot_wider(names_from=c(variable,p,dgp),values_from=value) %>%
  # reorder columns
  select(shrinkage,starts_with("rmse"),starts_with("signs")) %>%
  # add bold to 
  mutate(across(starts_with("rmse_"),function(x)ifelse(x==min(x),paste("\\bf",sprintf("%.3f",x)),sprintf("%.3f",x))))

out_big

# ------------------------------------------------------------- #
# plot elasticity matrices
# ------------------------------------------------------------- #

path = "/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/"
for(i in 1:length(data)){
  tmp = melt(data[[i]][[1]][[rep]]$B) %>%
    filter(Var1!=Var2) %>%
    mutate(Var1=as.factor(Var1),
           Var2=as.factor(Var2)) %>%
    ggplot(aes(x=Var1,y=fct_rev(Var2),fill=value)) +
    geom_tile(color="white") +
    labs(x="",y="") +
    scale_x_discrete(position = "top") +
    scale_fill_gradient2(low="#F8766D",mid="white",high="#00BFC4") +
    theme_minimal() +
    theme(legend.position="none",
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA))
    ggsave(filename=paste0(path,"simulated_elasticities_",i,".png"),
           tmp,
           height=5,width=5,
           bg = "transparent")
}
