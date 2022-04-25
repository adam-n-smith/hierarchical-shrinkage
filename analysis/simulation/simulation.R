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
  sparse.sparse = simdata_batch(n,p,rep,"hierarchical",prop_sparse=c(0,0.75,0),transform=FALSE),
  sparse = simdata_batch(n,p,rep,"sparse",prop_sparse=0.95)
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
# main simulations table
# ------------------------------------------------------------- #

# fit models across all dgps
fit1 = fit_parallel(data$dense.dense,Prior,Mcmc,models,dgp="dense.dense",tree=tree_true)
fit2 = fit_parallel(data$dense.transform,Prior,Mcmc,models,dgp="dense.transform",tree=tree_true)
fit3 = fit_parallel(data$sparse.sparse,Prior,Mcmc,models,dgp="sparse.sparse",tree=tree_true)
fit4 = fit_parallel(data$sparse,Prior,Mcmc,models,dgp="sparse.sparse",tree=tree_true)
fit = cbind(fit1,fit2,fit3,fit4)

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
  pivot_wider(names_from=c(variable,p,dgp),values_from=value) %>%
  # reorder columns
  select(shrinkage,starts_with("rmse"),starts_with("signs")) %>%
  # add bold to 
  mutate(across(starts_with("rmse_"),function(x)ifelse(x==min(x),paste("\\bf",sprintf("%.3f",x)),sprintf("%.3f",x))))

out
print(xtable(out),include.rownames=FALSE, sanitize.text.function=identity)
save(n,p,rep,data,models,Prior,Mcmc,fit,fit1,fit2,fit3,fit4,out,tree_true,
     file="analysis/simulation/output/sim-main.RData")

# stop
stopCluster(cl)

# ------------------------------------------------------------- #
# plot elasticity matrices
# ------------------------------------------------------------- #

path = "/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/"
for(i in 1:length(data)){
  tmp = melt(data[[i]][[rep]]$B) %>%
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