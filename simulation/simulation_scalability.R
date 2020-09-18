library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(foreach)
library(doParallel)
library(ggplot2)

sourceCpp(here("functions","horse_mcmc.cpp"))

simdata = function(n,p){
  Clist = NULL
  phivec = NULL
  Cphi = NULL
  d = 1
  np = double(p)
  for(i in 1:p){
    C = cbind(rep(1,n),matrix(runif(n*(d-1),-1,1),nrow=n))
    np[i] = ncol(C)
    phi = runif(d,-5,5)
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

n = 50
ntimes = 10
plist = c(50,100,250,500,1000)
numCores = detectCores()
registerDoParallel(numCores)
out = foreach(p=plist, .combine=rbind) %dopar% {
  
  data = simdata(n,p)
  Y = data$Y
  X = data$X
  runtime = mean(microbenchmark(drawbeta(Y,X,FALSE),times=ntimes)$time*1e-9)
  runtime_fast = mean(microbenchmark(drawbeta(Y,X,TRUE),times=ntimes)$time*1e-9)
  c(runtime,runtime_fast)
}

df = data.frame(p=rep(plist,times=2),
                runtime=c(as.vector(out)),
                method=rep(c("standard","fast"),each=length(plist)))
ggplot(df,aes(x=p,y=log(runtime),group=method)) +
  geom_line(aes(color=method,linetype=method)) +
  geom_point(aes(fill=method,color=method),shape=21, size=2, stroke=1) +
  labs(x="\n number of products (p)",y="time per iteration\n (in log seconds)\n") +
  theme_minimal()
dev.copy2pdf(file=here("figures","scalability.pdf"),width=6,height=3)
