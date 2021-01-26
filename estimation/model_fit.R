# --------------------------------------------------------- #
# fit statistics
# --------------------------------------------------------- #

end = Mcmc$R/Mcmc$keep
burn = end*0.5

rmse.ridge = double(end-burn)
rmse.lasso = double(end-burn)
rmse.horse = double(end-burn)
rmse.hierridge = double(end-burn)
rmse.hierlasso = double(end-burn)
rmse.hierhorse = double(end-burn)
rmse.hierhorseridge = double(end-burn)
rmse.hierhorselasso = double(end-burn)
rmse.hierhorsehorse = double(end-burn)
for(i in 1:(end-burn)){
  # ridge
  B = matrix(out_ridge$betadraws[burn+i,],p,p)
  Cphi = matrix(0,ntest,p)
  for(j in 1:p){
    Cphi[,j] = Clisttest[[j]]%*%out_ridge$phidraws[burn+i,(cumnphi[j]+1):cumnphi[j+1]]
  }
  rmse.ridge[i] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))
  
  # lasso
  # B = matrix(out_lasso$betadraws[burn+i,],p,p)
  # Cphi = matrix(0,ntest,p)
  # for(j in 1:p){
  #   Cphi[,j] = Clisttest[[j]]%*%out_lasso$phidraws[burn+i,(cumnphi[j]+1):cumnphi[j+1]]
  # }
  # rmse.lasso[i] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))
  
  # horseshoe
  B = matrix(out_horse$betadraws[burn+i,],p,p)
  Cphi = matrix(0,ntest,p)
  for(j in 1:p){
    Cphi[,j] = Clisttest[[j]]%*%out_horse$phidraws[burn+i,(cumnphi[j]+1):cumnphi[j+1]]
  }
  rmse.horse[i] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))
  
  # hierarchical ridge
  B = matrix(out_hierridge$betadraws[burn+i,],p,p)
  Cphi = matrix(0,ntest,p)
  for(j in 1:p){
    Cphi[,j] = Clisttest[[j]]%*%out_hierridge$phidraws[burn+i,(cumnphi[j]+1):cumnphi[j+1]]
  }
  rmse.hierridge[i] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))

  # hierarchical lasso
  # B = matrix(out_hierlasso$betadraws[burn+i,],p,p)
  # Cphi = matrix(0,ntest,p)
  # for(j in 1:p){
  #   Cphi[,j] = Clisttest[[j]]%*%out_hierlasso$phidraws[burn+i,(cumnphi[j]+1):cumnphi[j+1]]
  # }
  # rmse.hierlasso[i] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))

  # hierarchical horseshoe
  B = matrix(out_hierhorse$betadraws[burn+i,],p,p)
  Cphi = matrix(0,ntest,p)
  for(j in 1:p){
    Cphi[,j] = Clisttest[[j]]%*%out_hierhorse$phidraws[burn+i,(cumnphi[j]+1):cumnphi[j+1]]
  }
  rmse.hierhorse[i] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))
  
  # hierarchical horseshoe-ridge
  B = matrix(out_hierhorseridge$betadraws[burn+i,],p,p)
  Cphi = matrix(0,ntest,p)
  for(j in 1:p){
    Cphi[,j] = Clisttest[[j]]%*%out_hierhorseridge$phidraws[burn+i,(cumnphi[j]+1):cumnphi[j+1]]
  }
  rmse.hierhorseridge[i] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))
  
}

mean(rmse.ridge)
# mean(rmse.lasso)
mean(rmse.horse)
mean(rmse.hierridge)
# mean(rmse.hierlasso)
mean(rmse.hierhorse)
mean(rmse.hierhorseridge)
