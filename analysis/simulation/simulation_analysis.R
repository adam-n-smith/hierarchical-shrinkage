library(Rcpp)
library(RcppArmadillo)
library(ggplot2)

# data
Data = list(
  Y = Y,
  X = X,
  Clist = Clist,
  tree = tree,
  childrencounts = childrencounts,
  list = list,
  npar = npar
)

# priors
Prior = list(
  thetabar = 0,
  taubar = 10,
  betabarii = 0,
  taubarii = 10,
  Aphi = .01*diag(nphi),
  phibar = double(nphi),
  a = 5,
  b = 5
)

# Mcmc
Mcmc = list(
  R = 200,
  keep = 1
)

sourceCpp("code/horse_mcmc.cpp")
out = rSURhiershrinkage(Data,Prior,Mcmc,shrinkage="horseshoe",print=TRUE)
out = rSURhiershrinkage(Data,Prior,Mcmc,shrinkage="lasso",print=TRUE)
out = rSURhiershrinkage(Data,Prior,Mcmc,shrinkage="ridge",print=TRUE)
# out = rSURhorse(Data,Prior,Mcmc,hierarchical=TRUE,propagate=TRUE,print=TRUE)
# out = rSURhorse(Data,Prior,Mcmc,hierarchical=FALSE,propagate=FALSE,print=TRUE)

end = Mcmc$R/Mcmc$keep
burn = 0.5*end

# theta (first level)
plot(true.thetalist[[1]],apply(out$thetadraws[burn:end,1:npar[1]],2,mean),
     xlab="true",ylab="estimated",xlim=c(-3,3),ylim=c(-3,3))
abline(0,1)
matplot(out$thetadraws[,1:npar[1]],type="l",col=1:npar[1])
abline(h=true.thetalist[[1]],col=1:npar[1])

# theta (second level)
plot(true.thetalist[[2]],apply(out$thetadraws[burn:end,(npar[1]+1):cumsum(npar)[2]],2,mean),xlab="true",ylab="estimated")
abline(0,1)
matplot(out$thetadraws[,(npar[1]+1):cumsum(npar)[2]],type="l",col=1:cumsum(npar)[2])
abline(h=true.thetalist[[2]],col=1:npar[2])

# beta
plot(as.vector(B),apply(out$betadraws[burn:end,],2,mean),xlab="true",ylab="estimated")
abline(0,1)
matplot(out$betadraws[,1:10],type="l",col=1:p^2)
abline(h=as.vector(B)[1:10],col=1:npar[2])

# tau
matplot(sqrt(out$taudraws),type="l",col=1:L)
abline(h=sqrt(true.tau),col=1:L)

# phi
plot(phivec,apply(out$phidraws[burn:end,],2,mean),xlab="true",ylab="estimated")
abline(0,1)
matplot(out$phidraws,type="l",col=1:nphi)
abline(h=phivec,col=1:nphi)

# sigmasq 
plot(sqrt(sigmasq),apply(sqrt(out$sigmasqdraws[burn:end,]),2,mean),xlab="true",ylab="estimated")
abline(0,1)
matplot(sqrt(out$sigmasqdraws),type="l",col=1:p)
abline(h=sqrt(sigmasq),col=1:p)

# credible intervals
df = data.frame(true=c(true.thetalist[[1]],true.thetalist[[2]]),
                mean = c(apply(out$thetadraws[burn:end,1:npar[1]],2,mean),
                         apply(out$thetadraws[burn:end,(npar[1]+1):cumsum(npar)[2]],2,mean)),
                lb = c(apply(out$thetadraws[burn:end,1:npar[1]],2,quantile,.025),
                       apply(out$thetadraws[burn:end,(npar[1]+1):cumsum(npar)[2]],2,quantile,.025)),
                ub = c(apply(out$thetadraws[burn:end,1:npar[1]],2,quantile,.975),
                       apply(out$thetadraws[burn:end,(npar[1]+1):cumsum(npar)[2]],2,quantile,.975)),
                level = rep(c("Level 3 (Top)","Level 2 (Middle)"),npar[1:2]))
ggplot(df,aes(true,mean)) + 
  geom_point(shape=1) + 
  labs(x="\n true values", y="estimated values\n") +
  facet_wrap(.~ level) +
  geom_segment(aes(x=true,xend=true,y=lb,yend=ub)) + 
  geom_segment(aes(x=true-0.1,xend=true+0.1,y=lb,yend=lb)) + 
  geom_segment(aes(x=true-0.1,xend=true+0.1,y=ub,yend=ub)) + 
  geom_abline(aes(intercept=0,slope=1)) + 
  theme_minimal() +
  theme(strip.text = element_text(size=12))
dev.copy2pdf(file="figures/simulation_intervals.pdf",width=6,height=3)


# --------------------------------------------------------- #
# model fit
# --------------------------------------------------------- #

# data
Data = list(
  Y = Y,
  X = X,
  Clist = Clist,
  tree = tree,
  childrencounts = childrencounts,
  list = list,
  npar = npar
)

# priors
Prior = list(
  thetabar = 0,
  taubar = 10,
  betabarii = 0,
  taubarii = 10,
  Aphi = .01*diag(nphi),
  phibar = double(nphi),
  a = 5,
  b = 5
)

# Mcmc
Mcmc = list(
  R = 500,
  keep = 1
)

sourceCpp("code/horse_mcmc.cpp")

out.hierhorse = rSURhorse(Data,Prior,Mcmc,hierarchical=TRUE,propagate=TRUE,print=TRUE)
out.horse = rSURshrinkage(Data,Prior,Mcmc,shrinkage="horseshoe",print=TRUE)
out.lasso = rSURshrinkage(Data,Prior,Mcmc,shrinkage="lasso",print=TRUE)
out.ridge = rSURshrinkage(Data,Prior,Mcmc,shrinkage="ridge",print=TRUE)

end = Mcmc$R/Mcmc$keep
burn = 0.5*end

rmse.ridge = double(end-burn)
rmse.lasso = double(end-burn)
rmse.horse = double(end-burn)
rmse.hierhorse = double(end-burn)
for(i in 1:(end-burn)){
  # ridge
  B = matrix(out.ridge$betadraws[burn+i,],p,p)
  Cphi = matrix(out.ridge$phidraws[burn+i,],ntest,p,byrow=TRUE)
  rmse.ridge[i] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))
  
  # lasso
  B = matrix(out.lasso$betadraws[burn+i,],p,p)
  Cphi = matrix(out.lasso$phidraws[burn+i,],ntest,p,byrow=TRUE)
  rmse.lasso[i] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))

  # horseshoe
  B = matrix(out.horse$betadraws[burn+i,],p,p)
  Cphi = matrix(out.horse$phidraws[burn+i,],ntest,p,byrow=TRUE)
  rmse.horse[i] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))
  
  # hierarchical horseshoe
  B = matrix(out.hierhorse$betadraws[burn+i,],p,p)
  Cphi = matrix(out.hierhorse$phidraws[burn+i,],ntest,p,byrow=TRUE)
  rmse.hierhorse[i] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))
}

mean(rmse.ridge)
mean(rmse.lasso)
mean(rmse.horse)
mean(rmse.hierhorse)
