library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(here)
library(reshape2)

source(here("simulation","simulation_functions.R"))
source(here("functions","shrinkage_functions.R"))
sourceCpp(here("functions","shrinkage_mcmc.cpp"))

# dimensions
n = 50
p = 50
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
index = treeindex$index
list = treeindex$list
list_own = treeindex$list_own
npar = treeindex$npar
npar_own = treeindex$npar_own
L = ncol(tree)

# simulate data
settings = list(
  type = "hierarchical",
  tree = tree,
  counts = countchildren_cpp(tree),
  list = treeindex$list,
  npar = unlist(lapply(list,length)),
  # prop_sparse = c(0,0,0.9),
  # transform = TRUE
  prop_sparse = c(0,0,0),
  transform = FALSE
)
# settings = list(
#   type = "sparse",
#   prop_sparse = 0.75
# )
set.seed(1)
data = simdata(n,p,d,settings)
# data$npar = npar
# data$npar_own = npar_own
# data$tree = tree
# data$childrencounts = childrencounts

# plot elasticities
melt(data$B) %>%
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
        axis.text.y=element_blank())

# plot centered elasticities
melt(data$E) %>%
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
        axis.text.y=element_blank())

Pilist = list(L)
for(i in 1:L){
  Pi = matrix(0,p,max(tree[,i]))
  Pi[cbind(1:p,tree[,i])] = 1
  Pilist[[i]] = Pi
}

# data
Data = list(
  Y = data$Y,
  thetalist = data$thetalist,
  X = data$X, 
  Clist = data$Clist, 
  tree = tree,
  childrencounts = childrencounts, 
  list = treeindex$list, 
  list_own = treeindex$list_own,
  npar = treeindex$npar,
  npar_own = treeindex$npar_own,
  Pilist = Pilist
)

# priors
Prior = list(
  thetabar_own=0, 
  thetabar_cross=0, 
  taubarii=10,
  Aphi=.01*diag(p), 
  phibar=double(p),
  a=5, 
  b=5
)

# mcmc
Mcmc = list(
  R = 10000,
  initial_run=0,
  keep = 10,
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

matplot(out$thetadraws[,(npar[1]+1):cumsum(npar)[2]],type="l",col=rainbow(npar[2]))
abline(h=data$thetalist[[2]],col=rainbow(npar[2]))
plot(data$thetalist[[2]],apply(out$thetadraws[burn:end,(npar[1]+1):cumsum(npar)[2]],2,mean),xlab="true",ylab="estimated",col=diag(10)+1,pch=16)
abline(0,1)
boxplot(out$thetadraws[burn:end,(npar[1]+1):cumsum(npar)[2]])
points(data$thetalist[[2]],col=2,pch=16)

# compute rmse
Brmse = double(end)
for(r in 1:end){
  B = matrix(out$betadraws[r,],p,p)
  Brmse[r] = sqrt(mean((data$B - B)^2))
}
plot(Brmse,ylim=c(0,.7))
mean(Brmse[burn:end])

# theta (own effects)
matplot(out$thetaowndraws[,1:npar_own[1]],type="l",col=rainbow(npar_own[1]))
matplot(out$thetaowndraws[,(npar_own[1]+1):cumsum(npar_own)[2]],type="l",col=rainbow(npar_own[2]))
matplot(out$thetaowndraws[,(cumsum(npar_own)[2]+1):cumsum(npar_own)[3]],type="l",col=rainbow(p))
abline(h=diag(data$B),col=rainbow(p))

matplot(out$tauowndraws,type="l")

# theta (highest level)
matplot(out$thetadraws[,1:npar[1]],type="l",col=rainbow(npar[1]))
abline(h=data$thetalist[[1]],col=rainbow(npar[1]))
plot(data$thetalist[[1]],apply(out$thetadraws[burn:end,1:npar[1]],2,mean),
     xlab="true",ylab="estimated",xlim=c(-3,3),ylim=c(-3,3))
abline(0,1)

# theta (middle level)
matplot(out$thetadraws[,(npar[1]+1):cumsum(npar)[2]],type="l",col=rainbow(npar[2]))
abline(h=data$thetalist[[2]],col=rainbow(npar[2]))
plot(data$thetalist[[2]],apply(out$thetadraws[burn:end,(npar[1]+1):cumsum(npar)[2]],2,mean),xlab="true",ylab="estimated")
abline(0,1)

# theta (lowest level)
matplot(out$betadraws,type="l",col=rainbow(p^2))
abline(h=data$B,col=rainbow(p^2))
plot(data$thetalist[[3]],apply(out$thetadraws[burn:end,(cumsum(npar)[2]+1):cumsum(npar)[3]],2,mean),xlab="true",ylab="estimated")
abline(0,1)

# tau
matplot(sqrt(out$taudraws),type="l",col=1:L)
abline(h=sqrt(data$tau),col=1:L)
apply(out$taudraws[burn:end,],2,mean)

# plot shrinkage
lambdahat = apply(out$lambdadraws[burn:end,(cumsum(npar)[2]+1):cumsum(npar)[3]],2,mean)
tauhat = mean(out$taudraws[burn:end,3])
kappa = 1/(1+lambdahat*tauhat)
hist(kappa)

Kappa = matrix(NA,p,p)
Kappa[diag(p)!=1] = kappa
tmp = melt(Kappa)

# phi
matplot(out$phidraws,type="l",col=1:length(data$phivec))
abline(h=data$phivec,col=1:length(data$phivec))