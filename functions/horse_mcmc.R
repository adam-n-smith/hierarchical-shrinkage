R = 500
keep = 1

# storage
betadraws = matrix(0,R/keep,p*p)
thetadraws = matrix(0,R/keep,sum(npar[1:L]))
lambdadraws = matrix(0,R/keep,sum(npar[1:L]))
taudraws = matrix(0,R/keep,L)
sigmasqdraws = matrix(0,R/keep,q)

# initial
wchown = which(diag(p)[,1]==1)
betabar = double(p*q)
tau = true.tau
vbeta = rep(tau[L],p*q)
xitau = tau
thetalist = NULL
lambdalist = NULL
xilambdalist = NULL
for(ell in 1:L){
  thetalist[[ell]] = double(npar[ell])
  lambdalist[[ell]] = rep(1,npar[ell])
  xilambdalist[[ell]] = rep(1,npar[ell])
}
sigmasq = Sigma[1,1]
sigmasqvec = rep(sigmasq,each=p*q-q)
beta = B[,1]
XpX = crossprod(X)

for(rep in 1:R){
  
  # vbeta = rep(tau[L],p*q)
  # 
  # for(i in 1:q){
  #   XpXlami = XpX + diag(1/vbeta[(i*p-p+1):(i*p)])
  #   iroot = backsolve(chol(XpXlami),diag(p))
  #   X_XpXlami_Xp = X%*%(iroot%*%t(iroot))%*%t(X)
  #   rate = 0.5 * t(Y[,i]) %*% (diag(n) - X_XpXlami_Xp) %*% Y[,i]
  #   sigmasq[i] = 1/rgamma(1,0.5*n,rate)
  #   sigmasqvec[(i*p-p+1):(i*p-1)] = sigmasq[i]
  # }

  # ------------------------------------------------ #
  # SKU-LEVEL
  # ------------------------------------------------ #
  
  # product of previous lambda parmaeters
  Psi_Lmone = create_Psi_ellmone_cpp(lambdalist,childrencounts,list,npar,L)
  # 
  # # beta ------------------------------------------------
  # 
  # # mean and variance of beta
  # betabar[wchown] = betabarii
  # betabar[-wchown] = thetalist[[L-1]][list[[L]]]
  # vbeta[wchown] = tauii
  # vbeta[-wchown] = lambdalist[[L]]*Psi_Lmone*tau[L]
  # 
  # # draw
  # u = betabar + sqrt(vbeta)*rnorm(p*q)
  # delta = rnorm(n*q)
  # v = Xt%*%u + delta
  # iroot = matrix(0,n*q,n*q)
  # cptiroot = matrix(0,n*q,n*q)
  # for(i in 1:q){
  #   vari = vbeta[(p*i-p+1):(p*i)]
  #   r = chol((X/sqrt(sigmasq[i]))%*%diag(vari)%*%t(X/sqrt(sigmasq[i]))+diag(n))
  #   ir = backsolve(r,diag(n))
  #   iroot[(n*i-n+1):(n*i),(n*i-n+1):(n*i)] = ir
  #   cptiroot[(n*i-n+1):(n*i),(n*i-n+1):(n*i)] = crossprod(t(ir))
  # }
  # w = cptiroot%*%(yt-v)
  # beta = u + diag(vbeta)%*%crossprod(Xt,w)
  thetalist[[L]] = beta[-wchown]
  # 
  # # lambda ------------------------------------------------
  # 
  # # rate: compute sums of scaled differences in theta
  # tmp = (thetalist[[L]]-thetalist[[L-1]][list[[L]]])^2/(Psi_Lmone*tau[L])
  # rate = 1/xilambdalist[[L]] + 0.5*tmp
  # lambda = 1/rgamma(npar[L],shape=1,rate=rate)
  # xilambda = 1/rgamma(npar[L],shape=1,rate=1+1/lambda)
  # 
  # # store draws
  # lambdalist[[L]] = lambda
  # xilambdalist[[L]] = xilambda
  # 
  # # tau ------------------------------------------------
  levelsum = sum((thetalist[[L]]-thetalist[[L-1]][list[[L]]])^2/(sigmasqvec * lambdalist[[L]] * Psi_Lmone))
  # levelsum = sum((thetalist[[L]]-thetalist[[L-1]][list[[L]]])^2/(sigmasqvec * lambdalist[[L]]))
  shape = (npar[L]+1)/2
  rate = 1/xitau[L] + 0.5*levelsum
  tau[L] = 1/rgamma(1,shape=shape,rate=rate)
  xitau[L] = 1/rgamma(1,shape=1,rate=1+1/tau[L])
  
  # save
  if(rep%%keep==0){
    mkeep = rep/keep
    # betadraws[mkeep,] = as.vector(beta)
    # thetadraws[mkeep,] = unlist(thetalist)
    # lambdadraws[mkeep,] = unlist(lambdalist)
    taudraws[mkeep,] = tau
    sigmasqdraws[mkeep,] = sigmasq
  }
}  


matplot(sqrt(sigmasqdraws),type="l")
abline(h=sqrt(Sigma[1,1]))

matplot(sqrt(taudraws),type="l")
abline(h=sqrt(true.tau[L]))
