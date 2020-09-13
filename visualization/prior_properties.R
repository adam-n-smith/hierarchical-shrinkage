library(extraDistr)
library(reshape2)
library(tidyverse)

hc = function(x,R,log=FALSE){
  
  lambda = rhcauchy(R)
  den = matrix(0,length(x),R)
  for(i in 1:length(x)){
    den[i,] = dnorm(x[i],0,sd=lambda,log=FALSE)
  }
  den = apply(den,1,mean)
  out = ifelse(log,log(den),den)
  if(log){
    out = log(den)
  }
  else{
    out = den
  }
  return(out)
}

sphc = function(x,tau,R,log=FALSE){
  
  q = length(tau)
  lambda = matrix(rhcauchy(R*q),q,R)
  var = matrix(apply(lambda^2,2,cumprod),ncol=R)
  sd = apply(sqrt(matrix(tau,q,R)*var),2,sum)
  den = matrix(0,length(x),R)
  for(i in 1:length(x)){
    den[i,] = dnorm(x[i],0,sd,log=FALSE)
  }
  den = apply(den,1,mean)
  
  if(log){
    out = log(den)
  }
  else{
    out = den
  }
  return(out)
}

R = 10000

# origin
x = seq(-2,2, length=1001)
out.main = data.frame(
  x = x,
  Normal = dnorm(x),
  Laplace = dlaplace(x),
  HC = hc(x,R),
  SPHC2 = sphc(x,c(1,1),R),
  SPHC3 = sphc(x,c(1,1,1),R),
  SPHC4 = sphc(x,c(1,1,1,1),R)
)
df.main = pivot_longer(out.main,cols=-x,names_to="prior",values_to="density")

# tails
x = exp(seq(log(2), log(8), length=1001))
out.tail = data.frame(
  x = x,
  Normal = dnorm(x),
  Laplace = dlaplace(x),
  HC = hc(x,R),
  SPHC2 = sphc(x,c(1,1),R),
  SPHC3 = sphc(x,c(1,1,1),R),
  SPHC4 = sphc(x,c(1,1,1,1),R)
)
df.tail = pivot_longer(out.tail,cols=-x,names_to="prior",values_to="density")

# plot
df = df.main %>%
  bind_rows(df.tail) %>%
  mutate(plot=rep(c("Near the Origin","Tail Region"),c(nrow(df.main),nrow(df.tail))),
         density=ifelse(density>1.05,NA,density))
df %>%
  mutate(prior = factor(prior,levels=unique(prior))) %>%
  mutate(prior = recode(prior,SPHC2="SPHC (q=2)",SPHC3="SPHC (q=3)",SPHC4="SPHC (q=4)")) %>%
  ggplot(.,aes(x,density,color=prior)) + 
  geom_line(aes(linetype=prior),size=0.5) + 
  facet_wrap(.~plot,scales="free")+
  labs(x="\n delta",y="p(delta) \n") +
  theme_minimal() + 
  theme(legend.title=element_blank(),
        # legend.position="bottom",
        strip.text.x = element_text(size = 12),
        plot.title=element_text(hjust=0.5))
# theme(legend.title=element_blank(),
#         # legend.position="bottom",
#         plot.title=element_text(hjust=0.5),
#         panel.grid=element_blank(),
#         panel.border = element_rect(fill = NA),
#         axis.ticks = element_line(size = 0.25),
#         axis.ticks.length.x = unit(.2, "cm"),
#         axis.ticks.length.y = unit(.2, "cm"))
dev.copy2pdf(file="figures/prior_density.pdf",width=7,height=3)
