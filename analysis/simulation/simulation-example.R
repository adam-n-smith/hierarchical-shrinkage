library(tidyverse)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(here)

source(here("src","simulation-functions.R"))
source(here("src","shrinkage-functions.R"))
sourceCpp(here("src","shrinkage-mcmc.cpp"))

# ---------------------------------------------------------------------------- #
# generate data
# ---------------------------------------------------------------------------- #

# dimensions
n = 50
p = 20
d = 1

# tree
tree = matrix(c(rep(1:5,each=p/5),
                rep(1:10,each=p/10),
                1:p),nrow=p,ncol=3)
childrencounts = countchildren_cpp(tree)
treeindex = createindex(tree)
index = treeindex$index
list = treeindex$list
list_own = treeindex$list_own
npar = treeindex$npar
npar_own = treeindex$npar_own
L = ncol(tree)

# beta generated from hierarchical prior
settings = list(
  type = "hierarchical",
  tree = tree,
  counts = childrencounts,
  list = list,
  npar = npar,
  # prop_sparse = c(0,0,0.9),
  # transform = TRUE
  prop_sparse = c(0,0,0),
  transform = FALSE
)

# sparse beta
settings = list(
  type = "sparse",
  prop_sparse = 0.75
)

# simulate data
set.seed(1)
data = simdata(n,p,d,settings)

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

#
# Pilist = list(L)
# for(i in 1:L){
#   Pi = matrix(0,p,max(tree[,i]))
#   Pi[cbind(1:p,tree[,i])] = 1
#   Pilist[[i]] = Pi
# }


# ---------------------------------------------------------------------------- #
# estimation
# ---------------------------------------------------------------------------- #

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
  npar_own = treeindex$npar_own
)

# priors
Prior = list(
  thetabar_own = 0, 
  thetabar_cross = 0, 
  taubarii = 10,
  Aphi = .01*diag(p), 
  phibar = double(p),
  a = 5, 
  b = 5
)

# mcmc
Mcmc = list(
  R = 10000,
  initial_run=0,
  keep = 10,
  burn_pct = 0.5
)

# sparse
out = rSURshrinkage(Data,Prior,Mcmc,Shrinkage="ridge",print=TRUE)
out = rSURshrinkage(Data,Prior,Mcmc,Shrinkage="lasso",print=TRUE)
out = rSURshrinkage(Data,Prior,Mcmc,Shrinkage="horseshoe",print=TRUE)

# hierarchical (theta-ridge)
out = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="ridge",group="ridge"),print=TRUE)
out = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="lasso",group="ridge"),print=TRUE)
out = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="horseshoe",group="ridge"),print=TRUE)

# hierarchical (theta-horseshoe)
out = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="ridge",group="horseshoe"),print=TRUE)
out = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="lasso",group="horseshoe"),print=TRUE)
out = rSURhiershrinkage(Data,Prior,Mcmc,Shrinkage=list(product="horseshoe",group="horseshoe"),print=TRUE)

# ---------------------------------------------------------------------------- #
# summary
# ---------------------------------------------------------------------------- #

end = Mcmc$R/Mcmc$keep
burn = Mcmc$burn_pct*end

# estimation RMSE
Brmse = double(end-burn)
for(r in 1:(end-burn)){
  B = matrix(out$betadraws[burn+r,],p,p)
  Brmse[r] = sqrt(mean((data$B - B)^2))
}
mean(Brmse)

# theta-own
matplot(out$thetaowndraws[,1:npar_own[1]],type="l",col=rainbow(npar_own[1]))
matplot(out$thetaowndraws[,(npar_own[1]+1):cumsum(npar_own)[2]],type="l",col=rainbow(npar_own[2]))

# theta-cross (highest level)
matplot(out$thetadraws[,1:npar[1]],type="l",col=rainbow(npar[1]))
abline(h=data$thetalist[[1]],col=rainbow(npar[1]))
plot(data$thetalist[[1]],apply(out$thetadraws[burn:end,1:npar[1]],2,mean),
     xlab="true",ylab="estimated",xlim=c(-3,3),ylim=c(-3,3))
abline(0,1)

# theta-cross (middle level)
matplot(out$thetadraws[,(npar[1]+1):cumsum(npar)[2]],type="l",col=rainbow(npar[2]))
abline(h=data$thetalist[[2]],col=rainbow(npar[2]))
plot(data$thetalist[[2]],apply(out$thetadraws[burn:end,(npar[1]+1):cumsum(npar)[2]],2,mean),xlab="true",ylab="estimated")
abline(0,1)

# beta
matplot(out$betadraws,type="l",col=rainbow(p^2))
abline(h=data$B,col=rainbow(p^2))
plot(data$thetalist[[3]],apply(out$thetadraws[burn:end,(cumsum(npar)[2]+1):cumsum(npar)[3]],2,mean),xlab="true",ylab="estimated")
abline(0,1)

# check share of correct signs for beta
mean(sign(data$B) == sign(apply(out$betadraws[burn:end,],2,mean)))

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
plot(data$phivec,apply(out$phidraws[burn:end,],2,mean),xlab="true",ylab="estimated")
abline(0,1)
