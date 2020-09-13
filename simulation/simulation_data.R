library(gtools)
library(prodlim)
library(extraDistr)
library(ggplot2)
library(reshape2)
library(forcats)
library(Rcpp)
library(RcppArmadillo)

sourceCpp("code/horse_mcmc.cpp")

createindex = function(tree){
  
  # number of levels
  L = ncol(tree)
  
  # index positions of each element from current parameter vector
  K = max(tree[,1])
  indexpairlast = permutations(K,2,1:K,2,repeats.allowed=TRUE)
  
  # fill in parameter index
  out = NULL
  out[[1]] = rep(1,nrow(indexpairlast))
  
  # repeat for the remaining levels
  for(i in 2:L){
    subtree = unique(tree[,1:i])
    K = max(subtree)
    
    # group-level parameters assumed to be assumtric, use "combinations" to assume symmetry
    if(i<L){
      
      # index positions of each element from current parameter vector
      indexpair = permutations(K,2,1:K,repeats.allowed=TRUE)
      
      # index positions of each element from last level's parameter vector
      indexpairlast = permutations(max(tree[,i-1]),2,1:max(tree[,i-1]),repeats.allowed=TRUE)
      
      # match current positions with last positions
      index = matrix(subtree[indexpair,i-1],ncol=2)
      match = apply(index,1,function(x)row.match(x,indexpairlast))
      out[[i]] = match
    }
    
    # SKU-level parameters
    else{
      
      # index positions of each element from current parameter vector
      indexpair = permutations(K,2,1:K,repeats.allowed=TRUE)
      
      # index positions of each element from last level's parameter vector
      indexpairlast = permutations(max(tree[,i-1]),2,1:max(tree[,i-1]),repeats.allowed=TRUE)
      
      # match current positions with last positions
      index = matrix(subtree[indexpair,i-1],ncol=2)
      match = apply(index,1,function(x)row.match(x,indexpairlast))
      
      # ignore own elasticities
      own = apply(indexpair,1,function(x)x[1]==x[2])
      out[[i]] = match[!own]
    }
  }
  
  # save parameter index matrix
  l = out[[L]]
  index = l
  for(i in (L-1):1){
    l = out[[i]][l]
    index = rbind(l,index)
  }
  rownames(index) = NULL
  
  # number of parameters per level
  npar = apply(index,1,max)
  npar = c(npar[-1],ncol(index))
  
  return(list(list=out,index=index,npar=npar))
  
}

n = 50
p = 100
# tree = matrix(c(rep(1:4,each=p/4),
#                 1:p),nrow=p,ncol=2)
# tree = matrix(c(rep(1:2,each=p/2),
#                 rep(1:4,each=p/4),
#                 1:p),nrow=p,ncol=3)
tree = matrix(c(rep(1:5,each=p/5),
                rep(1:10,each=p/10),
                1:p),nrow=p,ncol=3)

# debug(createindex)
childrencounts = countchildren_cpp(tree)
treeindex = createindex(tree)
list = treeindex$list
index = treeindex$index
# npar = treeindex$npar
npar = unlist(lapply(list,length))
L = ncol(tree)

set.seed(1)

# sigmasq = runif(p,1,1)
sigmasq = runif(p,0,1)
Sigma = diag(sigmasq)

# global variances and initial mean
# tau = rep(4,L)
tau = 1/rgamma(L,2,2)
# tau = c(7,5,3)
mean = 0

# level 1
lambda = rep(1,npar[1])
theta = mean + sqrt(tau[1])*rnorm(npar[1])
# theta = c(3,-1,1,-3)
thetalist = list()
lambdalist = list()
thetalist[[1]] = theta
lambdalist[[1]] = lambda

# levels 2-L
for(ell in 2:L){
  
  # product of previous lambda parmaeters
  Psi = create_Psi_ellmone_cpp(lambdalist,childrencounts,list,npar,ell)
  
  # mean as theta from last level
  mean = theta[list[[ell]]]

  # draw local variances for current level
  # lambda = rhcauchy(npar[ell],sigma=1)
  # lambda = runif(npar[ell])
  lambda = 1/rgamma(npar[ell],10,10)
  # lambda = rep(1,npar[ell])

  # construct theta for current level
  thetasd = sqrt(lambda)*sqrt(Psi)*sqrt(tau[ell])
  if(ell==L){
    thetasd = rep(sqrt(sigmasq),each=p-1) * thetasd
  }
  theta = mean + thetasd*rnorm(npar[ell])

  # save
  thetalist[[ell]] = theta
  lambdalist[[ell]] = lambda
}

true.thetalist = thetalist
true.lambdalist = lambdalist
true.tau = tau

# elasticity matrix
wchown = which(diag(p)==1)
B = matrix(NA,p,p)
B[-wchown] = theta

# plot elasticities
df = melt(B)
df$Var1 = as.factor(df$Var1)
df$Var2 = as.factor(df$Var2)
ggplot(df, aes(x=Var1,y=fct_rev(Var2),fill=value)) + 
  geom_tile(color="white") + 
  labs(x="",y="") +
  scale_x_discrete(position = "top") +
  scale_fill_gradient2(low="#F8766D",mid="white",high="#00BFC4") +
  theme_minimal() +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.text.y=element_blank())

# intercepts and other controls
Clist = NULL
phivec = NULL
Cphi = NULL
d = 1
np = double(p)
for(i in 1:p){
  C = cbind(rep(1,n),matrix(runif(n*(d-1),-1,1),nrow=n))
  np[i] = ncol(C)
  # phi = runif(d,-5,5)
  phi = double(d)
  Clist[[i]] = C
  Cphi = c(Cphi,C%*%phi) 
  phivec = c(phivec,phi)
}
nphi = sum(np)
Cphi = matrix(Cphi,n,p)

# training data
betaown = -exp(rnorm(p))
B[which(diag(p)==1)] = betaown
X = matrix(rnorm(n*p),n,p)
# rho = runif(p)
# X = matrix(0,n,p)
# X[1,] = rnorm(p)
# for(t in 2:n){
#   X[t,] = X[t-1,] + rho*rnorm(p)
# }
# B = B*matrix(rbinom(p^2,1,.05),p,p)
error = matrix(rnorm(n*p),n,p)%*%chol(Sigma)
Y = X%*%B + Cphi + error

# test data
ntest = 25
Cphitest = NULL
for(i in 1:p){
  C = cbind(rep(1,ntest),matrix(runif(ntest*(d-1),-1,1),nrow=ntest))
  Cphitest = c(Cphitest,C%*%phivec[(i*d-d+1):(i*d)]) 
}
Cphitest = matrix(Cphitest,ntest,p)
Xtest = matrix(rnorm(ntest*p),ntest,p)
error = matrix(rnorm(ntest*p),ntest,p)%*%chol(Sigma)
Ytest = Xtest%*%B + Cphitest + error

