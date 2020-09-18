library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(here)
library(reshape2)

source(here("functions","simulation_functions.R"))
sourceCpp(here("functions","horse_mcmc.cpp"))

# dimensions
n = 50
p = 100
d = 1

# tree
tree = matrix(c(rep(1:2,each=p/2),
                rep(1:4,each=p/4),
                1:p),nrow=p,ncol=3)
childrencounts = countchildren_cpp(tree)
treeindex = createindex(tree)
list = treeindex$list
index = treeindex$index
npar = unlist(lapply(list,length))
L = ncol(tree)

# simulate data
tausq0 = .01
set.seed(1)
data = simdata_dense(n,p,d,tree,childrencounts,npar,tausq0)
data$npar = npar
data$tree = tree
data$list = list
data$childrencounts = childrencounts

# plot elasticities
B = data$B
B[diag(p)==1] = NA
melt(B) %>%
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

# priors
Prior = list(
  thetabar=0, 
  taubar=1, 
  betabarii=0, 
  taubarii=10,
  Aphi=.01*diag(p), 
  phibar=double(p),
  a=5, 
  b=5
)

# mcmc
Mcmc = list(
  R = 500,
  keep = 1
)

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

sourceCpp(here("functions","shrinkage_mcmc.cpp"))
out = rSURhiershrinkage(Data,Prior,Mcmc,shrinkage="ridge",print=TRUE)
# out = rSURhiershrinkage(Data,Prior,Mcmc,shrinkage="lasso",print=TRUE)
out = rSURhiershrinkage(Data,Prior,Mcmc,shrinkage="horseshoe",print=TRUE)

end = Mcmc$R/Mcmc$keep
burn = 0.75*end

# compute rmse
Brmse = double(end-burn)
for(r in 1:(end-burn)){
  B = matrix(out$betadraws[burn+r,],p,p)
  Brmse[r] = sqrt(mean((data$B[diag(p)!=1] - B[diag(p)!=1])^2))
}
mean(Brmse)

# theta (first level)
# plot(data$thetalist[[1]],apply(out$thetadraws[burn:end,1:npar[1]],2,mean),
#      xlab="true",ylab="estimated",xlim=c(-3,3),ylim=c(-3,3))
# abline(0,1)
matplot(out$thetadraws[,1:npar[1]],type="l",col=rainbow(npar[1]))
abline(h=data$thetalist[[1]],col=rainbow(npar[1]))

# theta (second level)
# plot(data$thetalist[[2]],apply(out$thetadraws[burn:end,(npar[1]+1):cumsum(npar)[2]],2,mean),xlab="true",ylab="estimated")
# abline(0,1)
matplot(out$thetadraws[,(npar[1]+1):cumsum(npar)[2]],type="l",col=rainbow(cumsum(npar)[2]))
abline(h=data$thetalist[[2]],col=rainbow(cumsum(npar)[2]))

# tau
matplot(sqrt(out$taudraws),type="l",col=rainbow(L))
abline(h=sqrt(data$tau),col=rainbow(L))

# # credible intervals
# df = data.frame(true=c(true.thetalist[[1]],true.thetalist[[2]]),
#                 mean = c(apply(out$thetadraws[burn:end,1:npar[1]],2,mean),
#                          apply(out$thetadraws[burn:end,(npar[1]+1):cumsum(npar)[2]],2,mean)),
#                 lb = c(apply(out$thetadraws[burn:end,1:npar[1]],2,quantile,.025),
#                        apply(out$thetadraws[burn:end,(npar[1]+1):cumsum(npar)[2]],2,quantile,.025)),
#                 ub = c(apply(out$thetadraws[burn:end,1:npar[1]],2,quantile,.975),
#                        apply(out$thetadraws[burn:end,(npar[1]+1):cumsum(npar)[2]],2,quantile,.975)),
#                 level = rep(c("Level 3 (Top)","Level 2 (Middle)"),npar[1:2]))
# ggplot(df,aes(true,mean)) + 
#   geom_point(shape=1) + 
#   labs(x="\n true values", y="estimated values\n") +
#   facet_wrap(.~ level) +
#   geom_segment(aes(x=true,xend=true,y=lb,yend=ub)) + 
#   geom_segment(aes(x=true-0.1,xend=true+0.1,y=lb,yend=lb)) + 
#   geom_segment(aes(x=true-0.1,xend=true+0.1,y=ub,yend=ub)) + 
#   geom_abline(aes(intercept=0,slope=1)) + 
#   theme_minimal() +
#   theme(strip.text = element_text(size=12))
# dev.copy2pdf(file="figures/simulation_intervals.pdf",width=6,height=3)
