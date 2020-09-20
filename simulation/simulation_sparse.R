library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(here)
library(reshape2)

source(here("functions","simulation_functions.R"))
sourceCpp(here("functions","shrinkage_mcmc.cpp"))

# dimensions
n = 50
p = 100
d = 1

# simulate data
set.seed(1)
data = simdata_sparse(n,p,d,0.95)

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

# data
Data = list(
  Y=data$Y, 
  X=data$X, 
  Clist=data$Clist
)

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
  R = 200,
  keep = 1
)

sourceCpp(here("functions","shrinkage_mcmc.cpp"))
out = rSURshrinkage(Data,Prior,Mcmc,shrinkage="ridge",print=TRUE)
out = rSURshrinkage(Data,Prior,Mcmc,shrinkage="horseshoe",print=TRUE)

end = Mcmc$R/Mcmc$keep
burn = 0.5*end

# compute rmse
Brmse = double(end-burn)
for(r in 1:(end-burn)){
  B = matrix(out$betadraws[burn+r,],p,p)
  Brmse[r] = sqrt(mean((data$B[diag(p)!=1] - B[diag(p)!=1])^2))
}
mean(Brmse)
plot(Brmse)

# plot shrinkage
lambdahat = apply(out$lambdadraws[burn:end,],2,mean)
tauhat = mean(out$taudraws[burn:end])
kappa = 1/(1+tauhat*lambdahat)
plot(kappa,data$B[diag(p)!=1])
