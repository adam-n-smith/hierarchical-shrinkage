library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(here)
library(reshape2)

# source(here("functions","simulation_functions.R"))
source(here("simulation","simulation_functions.R"))
source(here("functions","shrinkage_functions.R"))
sourceCpp(here("functions","shrinkage_mcmc.cpp"))

# dimensions
n = 50
p = 100
d = 1

# tree
tree = matrix(c(rep(1:5,each=p/5),
                rep(1:10,each=p/10),
                1:p),nrow=p,ncol=3)
# tree = matrix(c(rep(1:2,each=p/2),
#                 rep(1:4,each=p/4),
#                 1:p),nrow=p,ncol=3)
childrencounts = countchildren_cpp(tree)
treeindex = createindex(tree)
list = treeindex$list
index = treeindex$index
npar = unlist(lapply(list,length))
L = ncol(tree)

# simulate data
# debug(simdata_tree)
set.seed(1)
prior = list(
  type = "hierarchical",
  tree = tree,
  counts = countchildren_cpp(tree),
  list = treeindex$list,
  npar = unlist(lapply(list,length)),
  prop_sparse = c(0,0,0.75),
  transform = FALSE
)
prior = list(
  type = "sparse",
  prop_sparse = 0.75
)
data = simdata(n,p,d,prior)
# data = simdata_tree(n,p,d,tree,childrencounts,list,npar,c(0,0,.95),transform=TRUE)
# data = data_sparse_sparse[[1]][[1]]
data$npar = npar
data$tree = tree
data$list = list
data$childrencounts = childrencounts

mean(abs(data$E==0),na.rm=TRUE)
hist(data$B[diag(p)==0])

# plot elasticities
melt(data$B) %>%
  mutate(Var1=as.factor(Var1),
         Var2=as.factor(Var2)) %>%
  filter(Var1!=Var2) %>%
  ggplot(aes(x=Var1,y=fct_rev(Var2),fill=value)) + 
  geom_tile(color="white") + 
  labs(x="",y="") +
  scale_x_discrete(position = "top") +
  scale_fill_gradient2(low="#F8766D",mid="white",high="#00BFC4") +
  theme_minimal() +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.text.y=element_blank())

# plot centered elasticities
melt(data$E) %>%
  mutate(Var1=as.factor(Var1),
         Var2=as.factor(Var2)) %>%
  filter(Var1!=Var2) %>%
  ggplot(aes(x=Var1,y=fct_rev(Var2),fill=value)) + 
  geom_tile(color="white") + 
  labs(x="",y="") +
  scale_x_discrete(position = "top") +
  scale_fill_gradient2(low="#F8766D",mid="white",high="#00BFC4") +
  theme_minimal() +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.text.y=element_blank())

# data
Data = list(
  Y=data$Y, 
  X=data$X, 
  Clist=data$Clist, 
  tree=data$tree,
  childrencounts=data$childrencounts, 
  list=data$list, 
  npar=data$npar
)

# priors
Prior = list(
  thetabar=0, 
  betabarii=0, 
  taubarii=10,
  Aphi=.01*diag(p), 
  phibar=double(p),
  a=5, 
  b=5
)

# mcmc
Mcmc = list(
  R = 1000,
  initial_run=100,
  keep = 1,
  burn_pct = 0.5
)

sourceCpp(here("functions","shrinkage_mcmc.cpp"))

out = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="ridge",group="ridge"),print=TRUE)
out = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="lasso",group="ridge"),print=TRUE)
out = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="horseshoe",group="ridge"),print=TRUE)

out = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="ridge",group="horseshoe"),print=TRUE)
out = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="lasso",group="horseshoe"),print=TRUE)
out = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="horseshoe",group="horseshoe"),print=TRUE)

end = Mcmc$R/Mcmc$keep
burn = Mcmc$burn_pct*end
# burn = 1

# compute rmse
Brmse = double(end)
for(r in 1:end){
  B = matrix(out$betadraws[r,],p,p)
  Brmse[r] = sqrt(mean((data$B[diag(p)!=1] - B[diag(p)!=1])^2))
}
plot(Brmse)
mean(Brmse[burn:end])

# compute rmse
# rmse = double(end)
# for(r in 1:end){
#   B = matrix(out$betadraws[r,],p,p)
#   rmse[r] = sqrt(mean((data_out$X%*%B - data_out$X%*%data_out$B)^2))
# }
# plot(rmse)
# mean(rmse[burn:end])






# theta (second level)
matplot(out$thetadraws[,(npar[1]+1):cumsum(npar)[2]],type="l",col=rainbow(cumsum(npar)[2]))
abline(h=data$thetalist[[2]],col=rainbow(cumsum(npar)[2]))
plot(data$thetalist[[2]],apply(out$thetadraws[burn:end,(npar[1]+1):cumsum(npar)[2]],2,mean),xlab="true",ylab="estimated")
abline(0,1)

# theta (first level)
matplot(out$thetadraws[,1:npar[1]],type="l",col=rainbow(npar[1]))
abline(h=data$thetalist[[1]],col=rainbow(npar[1]))
plot(data$thetalist[[1]],apply(out$thetadraws[burn:end,1:npar[1]],2,mean),
     xlab="true",ylab="estimated",xlim=c(-3,3),ylim=c(-3,3))
abline(0,1)

# tau
matplot(sqrt(out$taudraws),type="l",col=1:L)
abline(h=sqrt(data$tau),col=1:L)

# # plot shrinkage
# lambdahat = apply(out$lambdadraws[burn:end,(cumsum(npar)[2]+1):cumsum(npar)[3]],2,mean)
# tauhat = mean(out$taudraws[burn:end,3])
# kappa = 1/(1+lambdahat*tauhat)
# hist(kappa)
