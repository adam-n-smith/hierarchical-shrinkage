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
p_vec = 100
nrep = 25

# ------------------------------------------------------------- #
# generate data sets
# ------------------------------------------------------------- #

set.seed(1)
data_dense_dense = simdata_batch(p_vec,nrep,"hierarchical",prop_sparse=c(0,0,0),transform=FALSE)
data_dense_transform = simdata_batch(p_vec,nrep,"hierarchical",prop_sparse=c(0,0,0.75),transform=TRUE)
data_sparse_sparse = simdata_batch(p_vec,nrep,"hierarchical",prop_sparse=c(0,0.75,0),transform=FALSE)
data_sparse = simdata_batch(p_vec,nrep,"sparse",prop_sparse=0.75)

# bind together data
data = list(dense.dense = data_dense_dense,
            dense.transform = data_dense_transform,
            sparse.sparse = data_sparse_sparse,
            sparse = data_sparse)

# data generating processes (beta_theta)
dgps = c("dense.dense",
         "dense.transform",
         "sparse.sparse",
         "sparse")

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

# models to fit
models = matrix(c("ridge","ridge",
                  "lasso","ridge",
                  "horseshoe","ridge",
                  "ridge","horseshoe",
                  "lasso","horseshoe",
                  "horseshoe","horseshoe"),
                ncol=2,byrow=TRUE)

# initialize clusters
cl = makeSOCKcluster(8)
registerDoSNOW(cl)
pb = txtProgressBar(max=nrep, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

# Mcmc
Mcmc = list(
  R = 1000,
  initial_run = 100,
  keep = 1,
  burn_pct = 0.75
)

# fit models across all dgps
fit = NULL
for(i in 1:length(dgps)){
  print(dgps[i])
  tmp = fit_parallel(data[[which(names(data)==dgps[i])]],Mcmc,p_vec,models,dgp=dgps[i])
  fit = cbind(fit,tmp)
}

# stop
stopCluster(cl)

# ------------------------------------------------------------- #
# print results
# ------------------------------------------------------------- #

out = data.frame(fit) %>%
  # add model labels
  mutate(shrinkage = rep(paste0(to_any_case(models[,2],case="upper_camel"),"/",
                                to_any_case(models[,1],case="upper_camel")),nrep),
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
# plot elasticity matrices
# ------------------------------------------------------------- #

path = "/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/"
for(i in 1:length(data)){
  tmp = melt(data[[i]][[1]][[nrep]]$B) %>%
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
