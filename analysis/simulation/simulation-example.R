library(tidyverse)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(here)

source(here("src","simulation-functions.R"))
source(here("src","shrinkage-functions.R"))
source(here("src","summary-functions.R"))
sourceCpp(here("src","shrinkage-mcmc.cpp"))

# ---------------------------------------------------------------------------- #
# generate data
# ---------------------------------------------------------------------------- #

# dimensions
n = 50
p = 100
d = 1

# tree
tree = matrix(c(1:p,
                rep(1:10,each=p/10),
                rep(1:5,each=p/5)),nrow=p,ncol=3)
objects = create_treeobjects(tree)
parindextree = objects$parindextree
parindextree_own = objects$parindextree_own
npar = objects$npar
npar_own = objects$npar_own

# beta generated from hierarchical prior
settings = list(
  type = "hierarchical",
  tree = tree,
  parindextree = parindextree,
  npar = npar,
  prop_sparse = c(0,0,0),
  transform = FALSE
)

set.seed(1)
data = simdata(n,p,d,settings)

# ---------------------------------------------------------------------------- #
# estimation
# ---------------------------------------------------------------------------- #

# data
Data = list(
  Y = data$Y,
  X = data$X, 
  Clist = data$Clist, 
  tree = tree,
  parindextree = parindextree,
  parindextree_own = parindextree_own,
  npar = npar,
  npar_own = npar_own
)

# priors
Prior = list(
  thetabar_own = 0, 
  thetabar_cross = 0, 
  Aphi = .01*diag(p), 
  phibar = double(p),
  a = 5, 
  b = 5
)

# mcmc
Mcmc = list(
  R = 10000,
  keep = 10
)
end = Mcmc$R/Mcmc$keep
burn = 0.5*end
cumnpar = cumsum(npar[-1])

out.sparse.ridge = rSURshrinkage(Data, Prior, Mcmc, Shrinkage="ridge",print=TRUE)
out.sparse.lasso = rSURshrinkage(Data, Prior, Mcmc, Shrinkage="lasso",print=TRUE)
out.sparse.horseshoe = rSURshrinkage(Data, Prior, Mcmc, Shrinkage="horseshoe",print=TRUE)

out.ridge.ridge = rSURhiershrinkage(Data, Prior, Mcmc, Shrinkage=list(product="ridge",group="ridge"),print=TRUE)
out.ridge.lasso = rSURhiershrinkage(Data, Prior, Mcmc, Shrinkage=list(product="lasso",group="ridge"),print=TRUE)
out.ridge.horseshoe = rSURhiershrinkage(Data, Prior, Mcmc, Shrinkage=list(product="horseshoe",group="ridge"),print=TRUE)

out.horseshoe.ridge = rSURhiershrinkage(Data, Prior, Mcmc, Shrinkage=list(product="ridge",group="horseshoe"),print=TRUE)
out.horseshoe.laso = rSURhiershrinkage(Data, Prior, Mcmc, Shrinkage=list(product="lasso",group="horseshoe"),print=TRUE)
out.horseshoe.horseshoe = rSURhiershrinkage(Data, Prior, Mcmc, Shrinkage=list(product="horseshoe",group="horseshoe"),print=TRUE)

# ---------------------------------------------------------------------------- #
# model summary
# ---------------------------------------------------------------------------- #

# choose model
out = out.ridge.ridge

# rmse
Brmse = double(end-burn)
for(r in 1:(end-burn)){
  B = matrix(out$betadraws[burn+r,],p,p)
  Brmse[r] = sqrt(mean((data$B-B)^2))
}
mean(Brmse)

# theta - top
matplot(out$thetadraws[,(cumnpar[1]+1):cumnpar[2]],type="l",col=rainbow(npar[3]))
abline(h=data$thetalist[[3]],col=rainbow(npar[3]))
plot(apply(out$thetadraws[,(cumnpar[1]+1):cumnpar[2]],2,mean),data$thetalist[[3]])
abline(0,1)

# theta - middle
matplot(out$thetadraws[,1:cumnpar[1]],type="l",col=rainbow(npar[2]))
abline(h=data$thetalist[[2]],col=rainbow(npar[2]))
plot(apply(out$thetadraws[,1:cumnpar[1]],2,mean),data$thetalist[[2]])
abline(0,1)

# beta
matplot(out$betadraws,type="l",col=rainbow(p^2))
abline(h=data$B,col=rainbow(p^2))
plot(data$B,apply(out$betadraws,2,mean),xlab="true",ylab="estimated")
abline(0,1)

# beta: check share of correct signs
mean(sign(data$B) == sign(apply(out$betadraws[burn:end,],2,mean)))

# tausq
matplot(log(out$tausqdraws),type="l",col=1:ncol(tree))

# phi
matplot(out$phidraws,type="l",col=1:length(data$phivec))
abline(h=data$phivec,col=1:length(data$phivec))
plot(data$phivec,apply(out$phidraws[burn:end,],2,mean),xlab="true",ylab="estimated")
abline(0,1)

# plot shrinkage
kappa = apply(1/(1+out$tausqdraws[burn:end,1]*out$Psidraws[burn:end,1:p^2]),2,mean)
hist(kappa)

# ---------------------------------------------------------------------------- #
# plot means and confidence intervals (hierarchical vs. standard)
# ---------------------------------------------------------------------------- #

model_hierarchical = list(out.ridge.ridge = out.ridge.ridge)
model_standard = list(out.sparse.ridge = out.sparse.ridge)

# compute means and credible intervals
par_top = summarize_elasticities_higher(model_hierarchical,burn:end,p,(cumnpar[1]+1):(cumnpar[2])) 
par_mid = summarize_elasticities_higher(model_hierarchical,burn:end,p,1:cumnpar[1]) 
par_bot = summarize_elasticities(model_hierarchical,burn:end,p) 
par_std = summarize_elasticities(model_standard,burn:end,p) 
par_true = c(data$thetalist[[3]],data$thetalist[[2]],data$B)
level_names = c("Hierarchical Shrinkage\ntheta (top level)",
                "Hierarchical Shrinkage\ntheta (middle level)",
                "Hierarchical Shrinkage\nbeta")

# tables
df.ridge.ridge = par_top %>%
  bind_rows(par_mid) %>%
  bind_rows(par_bot) %>%
  mutate(true = rep(par_true,each=length(models)),
         level = rep(level_names, length(models)*c(npar[3],npar[2],p^2))) %>%
  formatout(.,"plot","wrap")
df.sparse.ridge = par_std %>%
  mutate(true = as.vector(data$B),
         level = "Standard Shrinkage\nbeta") %>%
  formatout(.,"plot","wrap")
df = df.ridge.ridge %>%
  bind_rows(df.sparse.ridge) %>%
  mutate(level = factor(level,levels=unique(level)))

# plot
ggplot(df,aes(x=true,y=mean)) + 
  geom_segment(aes(x=true,xend=true,y=lower,yend=upper),color=4,lwd=1,alpha=0.25) +
  geom_point(size=0.25) +
  geom_abline(aes(intercept=0,slope=1)) +
  labs(x="true value", y="posterior mean") +
  facet_wrap(.~level, scales="free",nrow=1, labeller = ) +
  theme_shrink() + 
  theme(strip.text.x = element_text(size=11))

