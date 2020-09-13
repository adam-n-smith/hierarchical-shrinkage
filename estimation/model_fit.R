# --------------------------------------------------------- #
# fit statistics
# --------------------------------------------------------- #

end = Mcmc$R/Mcmc$keep
burn = end*0.8

rmse.ridge = double(end-burn)
rmse.lasso = double(end-burn)
rmse.horse = double(end-burn)
rmse.hierridge = double(end-burn)
rmse.hierlasso = double(end-burn)
rmse.hierhorse = double(end-burn)
for(i in 1:(end-burn)){
  # ridge
  B = matrix(out.ridge$betadraws[burn+i,],p,p)
  # Cphi = matrix(out.ridge$phidraws[burn+i,],ntest,p,byrow=TRUE)
  Cphi = matrix(0,ntest,p)
  for(j in 1:p){
    Cphi[,j] = Clisttest[[j]]%*%out.ridge$phidraws[burn+i,(cumnphi[j]+1):cumnphi[j+1]]
  }
  rmse.ridge[i] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))
  
  # lasso
  B = matrix(out.lasso$betadraws[burn+i,],p,p)
  # Cphi = matrix(out.lasso$phidraws[burn+i,],ntest,p,byrow=TRUE)
  Cphi = matrix(0,ntest,p)
  for(j in 1:p){
    Cphi[,j] = Clisttest[[j]]%*%out.lasso$phidraws[burn+i,(cumnphi[j]+1):cumnphi[j+1]]
  }
  rmse.lasso[i] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))
  
  # horseshoe
  B = matrix(out.horse$betadraws[burn+i,],p,p)
  # Cphi = matrix(out.horse$phidraws[burn+i,],ntest,p,byrow=TRUE)
  Cphi = matrix(0,ntest,p)
  for(j in 1:p){
    Cphi[,j] = Clisttest[[j]]%*%out.horse$phidraws[burn+i,(cumnphi[j]+1):cumnphi[j+1]]
  }
  rmse.horse[i] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))
  
  # hierarchical ridge
  B = matrix(out.hierridge$betadraws[burn+i,],p,p)
  # Cphi = matrix(out.hierridge$phidraws[burn+i,],ntest,p,byrow=TRUE)
  Cphi = matrix(0,ntest,p)
  for(j in 1:p){
    Cphi[,j] = Clisttest[[j]]%*%out.hierridge$phidraws[burn+i,(cumnphi[j]+1):cumnphi[j+1]]
  }
  rmse.hierridge[i] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))
  
  # # hierarchical lasso
  # # B = matrix(out.hierlasso$betadraws[burn+i,],p,p)
  # # # Cphi = matrix(out.hierlasso$phidraws[burn+i,],ntest,p,byrow=TRUE)
  # # Cphi = matrix(0,ntest,p)
  # # for(j in 1:p){
  # #   Cphi[,j] = Clisttest[[j]]%*%out.hierlasso$phidraws[burn+i,(cumnphi[j]+1):cumnphi[j+1]]
  # # }
  # # rmse.hierlasso[i] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))
  # 
  # hierarchical horseshoe
  B = matrix(out.hierhorse$betadraws[burn+i,],p,p)
  # Cphi = matrix(out.hierhorse$phidraws[burn+i,],ntest,p,byrow=TRUE)
  Cphi = matrix(0,ntest,p)
  for(j in 1:p){
    Cphi[,j] = Clisttest[[j]]%*%out.hierhorse$phidraws[burn+i,(cumnphi[j]+1):cumnphi[j+1]]
  }
  rmse.hierhorse[i] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))
}

mean(rmse.ridge)
mean(rmse.lasso)
mean(rmse.horse)
mean(rmse.hierridge)
mean(rmse.hierhorse)
sd(rmse.ridge)
sd(rmse.lasso)
sd(rmse.horse)
sd(rmse.hierridge)
sd(rmse.hierhorse)

matplot(cbind(rmse.ridge,rmse.lasso,rmse.horse,
              rmse.hierridge,rmse.hierlasso,rmse.hierhorse),type="l",ylab="")
