library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(tidyverse)
library(here)

source(here("src","summary-functions.R"))
sourceCpp(here("src","shrinkage-mcmc.cpp"))


simdata_scale = function(n,p){
  B = matrix(rnorm(p*p),p,p)
  X = matrix(rnorm(n*p),n,p)
  Sigma = diag(p)
  error = matrix(rnorm(n*p),n,p)%*%chol(Sigma)
  Y = X%*%B + error
  return(list(Y=Y,X=X))
}

set.seed(1)
n = 100
ntimes = 10
pvec = seq(100,1000,by=100)
runtime = matrix(0,length(pvec),3)
runtime.fast = matrix(0,length(pvec),3)
j = 1
for(p in pvec){
  
  data = simdata_scale(n,p)
  Y = data$Y
  X = data$X
  out = summary(microbenchmark(draw_beta(Y,X,fast=FALSE), 
                               draw_beta(Y,X,fast=TRUE), 
                               times=ntimes, unit="s", setup=set.seed(1)))
  runtime[j,1] = out$lq[1]
  runtime[j,2] = out$mean[1]
  runtime[j,3] = out$uq[1]
  runtime.fast[j,1] = out$lq[2]
  runtime.fast[j,2] = out$mean[2]
  runtime.fast[j,3] = out$uq[2]
  
  print(p)
  j = j+1
  
}

# gains
cbind(pvec,runtime[,2]/runtime.fast[,2])

# table
timetable = rbind(runtime,runtime.fast)
colnames(timetable) = c("lq","mean","uq")
df = data.frame(p=rep(pvec,times=2),
                timetable,
                method=rep(c("standard (inverts X'X)","fast (inverts XX')"),each=length(pvec)))

# plot
ggplot(df,aes(x=p,y=log(mean),group=method)) +
  geom_line(aes(color=method, linetype=method)) +
  scale_x_continuous(breaks=pvec) +
  geom_point(aes(fill=method, shape=method, color=method), size=2, stroke=1) +
  labs(x="number of products",y="time per iteration\n (in log seconds)") +
  theme_shrink()
ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/scalability.png",
       height=3.5,width=7)
# dev.copy2pdf(file="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/scalability.pdf",
#              width=6,height=3)
