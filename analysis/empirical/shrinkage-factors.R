# estimates of tau
median(out.sparse.ridge$taudraws[burn:end])
median(out.sparse.lasso$taudraws[burn:end])
median(out.sparse.horseshoe$taudraws[burn:end])
apply(out.ridge.ridge$taudraws[burn:end,],2,median)
apply(out.ridge.lasso$taudraws[burn:end,],2,median)
apply(out.ridge.horseshoe$taudraws[burn:end,],2,median)
apply(out.horseshoe.ridge$taudraws[burn:end,],2,median)
apply(out.horseshoe.lasso$taudraws[burn:end,],2,median)
apply(out.horseshoe.horseshoe$taudraws[burn:end,],2,median)

# sparse lasso
lamtau = out.sparse.lasso$lambdadraws[burn:end,] *
  matrix(out.sparse.lasso$taudraws[burn:end],length(burn:end),npar[3])
kappa = 1/(1+lamtau)
kappa = apply(kappa,2,median)
min(kappa)
quantile(kappa,0.01)
quantile(kappa,0.99)

# sparse horseshoe
lamtau = out.sparse.horseshoe$lambdadraws[burn:end,] *
  matrix(out.sparse.horseshoe$taudraws[burn:end],length(burn:end),npar[3])
kappa = 1/(1+lamtau)
kappa = apply(kappa,2,median)
min(kappa)
quantile(kappa,0.01)
quantile(kappa,0.99)

# hierarchical ridge-lasso
lamtau = out.ridge.lasso$lambdadraws[burn:end,(cumsum(npar)[2]+1):cumsum(npar)[3]] *
  matrix(out.ridge.lasso$taudraws[burn:end,3],length(burn:end),npar[3])
kappa = 1/(1+lamtau)
kappa = apply(kappa,2,median)
min(kappa)
quantile(kappa,0.01)
quantile(kappa,0.99)

# hierarchical ridge-horseshoe
lamtau = out.ridge.horseshoe$lambdadraws[burn:end,(cumsum(npar)[2]+1):cumsum(npar)[3]] *
  matrix(out.ridge.horseshoe$taudraws[burn:end,3],length(burn:end),npar[3])
kappa = 1/(1+lamtau)
kappa = apply(kappa,2,median)
min(kappa)
quantile(kappa,0.01)
quantile(kappa,0.99)

# hierarchical horseshoe-ridge
lamtau = matrix(0,length(burn:end),npar[3])
for(i in 1:length(burn:end)){
  lamtau[i,] = replace_cpp(list[[3]],out.horseshoe.ridge$lambdadraws[burn-1+i,(cumsum(npar)[1]+1):cumsum(npar)[2]]) *
      out.horseshoe.ridge$taudraws[burn-1+i,3]
}
kappa = 1/(1+lamtau)
kappa = apply(kappa,2,median)
min(kappa)
quantile(kappa,0.01)
quantile(kappa,0.99)

# hierarchical horseshoe-lasso
lamtau = matrix(0,length(burn:end),npar[3])
for(i in 1:length(burn:end)){
  lamtau[i,] = replace_cpp(list[[3]],out.horseshoe.lasso$lambdadraws[burn-1+i,(cumsum(npar)[1]+1):cumsum(npar)[2]]) *
    out.horseshoe.lasso$lambdadraws[burn-1+i,(cumsum(npar)[2]+1):cumsum(npar)[3]] *
    out.horseshoe.lasso$taudraws[burn-1+i,3]
}
kappa = 1/(1+lamtau)
kappa = apply(kappa,2,median)
min(kappa)
quantile(kappa,0.01)
quantile(kappa,0.99)

# hierarchical horseshoe-horseshoe
lamtau = matrix(0,length(burn:end),npar[3])
for(i in 1:length(burn:end)){
  lamtau[i,] = replace_cpp(list[[3]],out.horseshoe.horseshoe$lambdadraws[burn-1+i,(cumsum(npar)[1]+1):cumsum(npar)[2]]) *
    out.horseshoe.horseshoe$lambdadraws[burn-1+i,(cumsum(npar)[2]+1):cumsum(npar)[3]] *
    out.horseshoe.horseshoe$taudraws[burn-1+i,3]
}
kappa = 1/(1+lamtau)
kappa = apply(kappa,2,median)
min(kappa)
quantile(kappa,0.01)
quantile(kappa,0.99)

# hierarchical horseshoe-ridge (subcategory)
lamtau = out.horseshoe.ridge$lambdadraws[burn:end,(cumsum(npar)[1]+1):cumsum(npar)[2]] *
  matrix(out.horseshoe.ridge$taudraws[burn:end,2],length(burn:end),npar[2])
kappa = 1/(1+lamtau)
kappa = apply(kappa,2,median)
min(kappa)
quantile(kappa,0.01)
quantile(kappa,0.99)

# hierarchical horseshoe-lasso (subcategory)
lamtau = out.horseshoe.lasso$lambdadraws[burn:end,(cumsum(npar)[1]+1):cumsum(npar)[2]] *
  matrix(out.horseshoe.lasso$taudraws[burn:end,2],length(burn:end),npar[2])
kappa = 1/(1+lamtau)
kappa = apply(kappa,2,median)
min(kappa)
quantile(kappa,0.01)
quantile(kappa,0.99)

# hierarchical horseshoe-horseshoe (subcategory)
lamtau = out.horseshoe.horseshoe$lambdadraws[burn:end,(cumsum(npar)[1]+1):cumsum(npar)[2]] *
  matrix(out.horseshoe.horseshoe$taudraws[burn:end,2],length(burn:end),npar[2])
kappa = 1/(1+lamtau)
kappa = apply(kappa,2,median)
min(kappa)
quantile(kappa,0.01)
quantile(kappa,0.99)
