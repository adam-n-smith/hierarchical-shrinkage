library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(tidyverse)
library(here)

sourceCpp(here("src","shrinkage-mcmc.cpp"))

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
pvec = seq(100,300,by=100)
runtime = double(length(pvec))
runtime.fast = double(length(pvec))
j = 1
for(p in pvec){
  
  set.seed(1)
  data = simdata_scale(n,p)
  Y = data$Y
  X = data$X
  out = microbenchmark(drawbeta(Y,X,FALSE), times=ntimes, setup=set.seed(1))
  out.fast = microbenchmark(drawbeta(Y,X,TRUE), times=ntimes, setup=set.seed(1))
  
  runtime[j] = mean(out$time)*1e-9
  runtime.fast[j] = mean(out.fast$time)*1e-9
  print(p)
  j = j+1
}

beta1 = drawbeta(Y,X,FALSE)
beta2 = drawbeta(Y,X,TRUE)
abs(beta1-beta2)<1e-5

df = data.frame(p=rep(pvec,times=2),
                runtime=c(runtime,runtime.fast),
                method=rep(c("standard","fast"),each=length(pvec)))
ggplot(df,aes(x=p,y=log(runtime),group=method)) +
  geom_line(aes(color=method, linetype=method)) +
  scale_x_continuous(breaks=pvec) +
  geom_point(aes(fill=method, shape=method, color=method), size=2, stroke=1) + 
  labs(x="\n number of products (p)",y="time per iteration\n (in log seconds)\n") +
  ylim(-4,4) +
  theme_minimal()
# dev.copy2pdf(file="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/scalability.pdf",width=7,height=4)
