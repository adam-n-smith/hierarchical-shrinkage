library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(tidyverse)
library(here)

sourceCpp(here("functions","shrinkage_mcmc.cpp"))

simdata_scale = function(n,p){
  Clist = NULL
  phivec = NULL
  Cphi = NULL
  d = 1
  np = double(p)
  for(i in 1:p){
    C = cbind(rep(1,n),matrix(rnorm(n*(d-1)),nrow=n))
    np[i] = ncol(C)
    phi = double(d)
    Clist[[i]] = C
    Cphi = c(Cphi,C%*%phi) 
    phivec = c(phivec,phi)
  }
  nphi = sum(np)
  Cphi = matrix(Cphi,n,p)
  B = matrix(p*p,p,p)
  X = matrix(rnorm(n*p),n,p)
  Sigma = diag(p)
  error = matrix(rnorm(n*p),n,p)%*%chol(Sigma)
  Y = X%*%B + Cphi + error
  
  return(list(Y=Y,X=X))
}

n = 100
ntimes = 1
plist = seq(100,1000,by=100)
runtime = double(length(plist))
runtime.fast = double(length(plist))
j = 1
for(p in plist){
  
  set.seed(1)
  data = simdata_scale(n,p)
  Y = data$Y
  X = data$X
  out = microbenchmark(drawbeta(Y,X,FALSE),times=ntimes)
  out.fast = microbenchmark(drawbeta(Y,X,TRUE),times=ntimes)
  
  runtime[j] = mean(out$time)*1e-9
  runtime.fast[j] = mean(out.fast$time)*1e-9
  print(p)
  j = j+1
}

df = data.frame(p=rep(plist,times=2),
                runtime=c(runtime,runtime.fast),
                method=rep(c("standard","fast"),each=length(plist)))
ggplot(df,aes(x=p,y=log(runtime),group=method)) +
  geom_line(aes(color=method, linetype=method)) +
  scale_x_continuous(breaks=plist) +
  geom_point(aes(fill=method, shape=method, color=method), size=2, stroke=1) + 
  labs(x="\n number of products (p)",y="time per iteration\n (in log seconds)\n") +
  theme_minimal()
dev.copy2pdf(file="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/scalability.pdf",width=7,height=4)
