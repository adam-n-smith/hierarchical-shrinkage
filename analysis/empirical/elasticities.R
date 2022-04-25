library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(GGally)

sourceCpp(here("src","shrinkage-mcmc.cpp"))
source(here("analysis","empirical","estimation-build.R"))

# function to compute summaries of price elasticities
summarize_elasticities = function(out_list,index,p){
  
  mean = matrix(0,p^2,length(out_list))
  lower = matrix(0,p^2,length(out_list))
  upper = matrix(0,p^2,length(out_list))
  for(m in 1:length(out_list)){
    mean[,m] = apply(out_list[[m]]$betadraws[index,],2,mean)
    lower[,m] = apply(out_list[[m]]$betadraws[index,],2,function(x)quantile(x,0.025))
    upper[,m] = apply(out_list[[m]]$betadraws[index,],2,function(x)quantile(x,0.975))
  }
  cols = str_remove(names(out_list),"out.")
  colnames(mean) = paste("mean",cols,sep="_")
  colnames(lower) = paste("lower",cols,sep="_")
  colnames(upper) = paste("upper",cols,sep="_")
  
  out = data.frame(cbind(mean,lower,upper)) %>%
    mutate(own = as.vector(diag(p)==1),
           id = 1:p^2) %>%
    pivot_longer(-c(own,id), names_to = c("variable","theta","beta"), names_pattern="(.*)_(.*)[.](.*)") %>%
    mutate(model=paste(theta,beta,sep="_")) %>%
    pivot_wider(names_from=variable,values_from=value) %>%
    mutate(beta = case_when(beta == "ridge" ~ paste("beta","-ridge"),
                            beta == "lasso" ~ paste("beta","-lasso"),
                            beta == "horseshoe" ~ paste("beta","-horseshoe")),
           beta = factor(beta, levels=unique(beta)),
           theta = case_when(theta == "ridge" ~ paste("theta","-ridge"),
                             theta == "horseshoe" ~ paste("theta","-horseshoe"),
                             theta == "sparse" ~ paste("sparse")),
           theta = factor(theta, levels=unique(theta)))
  
  return(out)
  
}

# function to compute summaries of promotion effects
summarize_promos = function(out_list,index,p){
  
  mean = matrix(nrow=p,ncol=length(out_list))
  wchcols = cumnphi[-1][which(nphivec==max(nphivec))]
  for(m in 1:length(out_list)){
    mean[which(nphivec==max(nphivec)),m] = apply(out_list[[m]]$phidraws[index,wchcols],2,mean)
  }
  cols = str_remove(names(out_list),"out.")
  colnames(mean) = paste("mean",cols,sep="_")
  
  out = data.frame(mean) %>%
    mutate(id = 1:p) %>%
    pivot_longer(-id,names_to = c("variable","theta","beta"), names_pattern="(.*)_(.*)[.](.*)") %>%
    mutate(model=paste(theta,beta,sep="_")) %>%
    pivot_wider(names_from=variable,values_from=value) %>%
    mutate(beta = case_when(beta == "ridge" ~ paste("beta","-ridge"),
                            beta == "lasso" ~ paste("beta","-lasso"),
                            beta == "horseshoe" ~ paste("beta","-horseshoe")),
           beta = factor(beta, levels=unique(beta)),
           theta = case_when(theta == "ridge" ~ paste("theta","-ridge"),
                             theta == "horseshoe" ~ paste("theta","-horseshoe"),
                             theta == "sparse" ~ paste("sparse")),
           theta = factor(theta, levels=unique(theta)))
  
  return(out)
  
}

# function to compute summaries of markups
summarize_markups = function(out_list,index,p){
  
  mean = matrix(0,p,length(out_list))
  mkp = matrix(0,p,length(out_list))
  wchown = as.vector(diag(matrix(1:p^2,p,p)))
  for(m in 1:length(out_list)){
    b = out_list[[m]]$betadraws[index,wchown]
    mean[,m] = apply(b,2,mean)
    mkp[,m] = apply(b/(1+b),2,mean)
  }
  cols = str_remove(names(out_list),"out.")
  colnames(mean) = paste("mean",cols,sep="_")
  colnames(mkp) = paste("mkp",cols,sep="_")

  out = data.frame(cbind(mean,mkp)) %>%
    mutate(id = 1:p) %>%
    pivot_longer(-id,names_to = c("variable","theta","beta"), names_pattern="(.*)_(.*)[.](.*)") %>%
    mutate(model=paste(theta,beta,sep="_")) %>%
    pivot_wider(names_from=variable,values_from=value) %>%
    mutate(beta = case_when(beta == "ridge" ~ paste("beta","-ridge"),
                            beta == "lasso" ~ paste("beta","-lasso"),
                            beta == "horseshoe" ~ paste("beta","-horseshoe")),
           beta = factor(beta, levels=unique(beta)),
           theta = case_when(theta == "ridge" ~ paste("theta","-ridge"),
                             theta == "horseshoe" ~ paste("theta","-horseshoe"),
                             theta == "sparse" ~ paste("sparse")),
           theta = factor(theta, levels=unique(theta)))
  
  return(out)
  
}

# load data
all_files = dir(here("analysis","output"))
out_files = paste0("analysis/output/",all_files[!str_detect(all_files,"sim-")])
lapply(out_files,load,.GlobalEnv)

# choose models
models = mget(rev(ls()[str_detect(ls(),"out[.]")]))

# compute
end = nrow(models[[1]]$betadraws)
burn = 0.5*end
elasticities = summarize_elasticities(models,burn:end,p)
promotions = summarize_promos(models,burn:end,p)
markups = summarize_markups(models,burn:end,p)

# summary stats (own)
elasticities %>%
  filter(own) %>%
  mutate(sig = !(lower<0 & upper>0)) %>%
  group_by(model) %>%
  summarise(mean_all = mean(mean), 
            mean_sig = mean(mean[sig]),
            n_sig = sum(sig),
            share_neg = mean(mean<0))

# summary stats (cross)
elasticities %>%
  filter(!own) %>%
  mutate(sig = !(lower<0 & upper>0)) %>%
  group_by(model) %>%
  summarise(mean_all = mean(mean), 
            mean_sig = mean(mean[sig]),
            n_sig = sum(sig),
            share_neg = mean(mean<0))

# own
elasticities %>%
  filter(own) %>%
  ggplot(aes(x=mean)) +
  geom_histogram(alpha=0.75) +
  xlim(-5,3) +
  facet_grid(theta~beta, labeller = label_parsed, scales="free") +
  labs(x="") +
  theme_minimal() +
  theme(plot.title = element_text(face="bold",hjust=0.5))
# ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/elasticities-beta-own.png",
#        height=4,width=7)

# cross
elasticities %>%
  filter(!own) %>%
  # mutate(significant = ifelse(lower<0 & upper>0,"no","yes")) %>%
  # ggplot(aes(x=mean, fill=significant)) + 
  ggplot(aes(x=mean)) +
  geom_histogram(alpha=0.75) +
  xlim(-0.25,0.25) +
  labs(x="") +
  theme_minimal() +
  facet_grid(theta~beta, labeller = label_parsed, scales="free_y") +
  theme(plot.title = element_text(face="bold",hjust=0.5))
# ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/elasticities-beta-cross.png",
#        height=4,width=7)

# own promos
promotions %>%
  ggplot(aes(x=mean)) +
  geom_histogram(alpha=0.75) +
  labs(x="") +
  theme_minimal() +
  facet_grid(theta~beta, labeller = label_parsed, scales="free_y") +
  theme(plot.title = element_text(face="bold",hjust=0.5))
# ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/promotions.png",
#        height=4,width=7)


# pairwise own
# elasticities %>%
#   filter(own) %>%
#   mutate(model = paste0(theta,"\n",beta)) %>%
#   select(id,model,mean) %>%
#   pivot_wider(names_from=model,values_from=mean) %>% 
#   select(-id) %>%
#   ggpairs(aes(alpha = 0.4)) +
#   theme_minimal()
# ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/elasticities-pairwise.png",
#        height=10,width=10)

# table
library(xtable)
elasticities %>%
  group_by(theta,beta,own) %>%
  summarise(neg = mean(mean<0),
            low = quantile(mean,0.25),
            mid = quantile(mean,0.5),
            high = quantile(mean,0.75)) %>%
  pivot_wider(names_from=own,values_from=c(neg,low,mid,high)) %>%
  select(theta,beta,contains("TRUE"),contains("FALSE")&!contains("neg")) %>%
  xtable(digits=3)


# effects of hierarchical shrinkage
qtlow = quantile(apply(prices[inweeks,-(1:2)],2,sd),0.25)
qthigh = quantile(apply(prices[inweeks,-(1:2)],2,sd),0.75)
sds = data.frame(id=1:p^2,
                 Var1=rep(1:p,each=p),
                 Var2=rep(1:p,times=p),
                 sd=rep(apply(prices[inweeks,-(1:2)],2,sd),times=p))
elasticities %>%
  left_join(sds) %>%
  mutate(sd = ifelse(sd<qtlow,"Q1",ifelse(sd>qthigh,"Q3","(Q1,Q3)")),
         sd = relevel(factor(sd),"Q1")) %>% 
  filter(own) %>%
  select(id,beta,theta,mean,sd) %>%
  pivot_wider(names_from=theta,values_from=mean) %>%
  pivot_longer(-c(id,beta,sd,sparse),names_to="theta",values_to="mean") %>%
  arrange(desc(sd)) %>%
  ggplot(aes(x=sparse,y=mean)) +
  geom_point(aes(color=sd,shape=sd)) +
  geom_abline(slope=1,intercept=0) +
  scale_color_manual(values=c("#F8766D","grey90","#00BFC4")) +
  scale_shape_manual(values=c(17,16,15)) +
  xlim(-3,1) +
  ylim(-3,1) +
  labs(x="beta (sparse prior)",y="beta (hierarchical prior)",shape="SD(price)",color="SD(price)") +
  facet_grid(theta~beta, labeller=label_parsed) +
  theme_minimal() 
ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/elasticities-beta-scatter.png",
       height=5,width=8)

# markups
markups %>%
  mutate(mkp = ifelse(mkp>3,3,mkp)) %>%
  ggplot(aes(x=mean,y=mkp)) +
  geom_point() +
  geom_vline(xintercept=-1) +
  labs(x="") +
  ylim(0,3) +
  theme_minimal() +
  facet_grid(theta~beta, labeller = label_parsed, scales="free_y") +
  theme(plot.title = element_text(face="bold",hjust=0.5))
