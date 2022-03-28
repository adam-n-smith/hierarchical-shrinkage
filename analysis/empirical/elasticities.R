library(GGally)

# function to compute elasticities
elast = function(out_list,index,p){
  
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
    pivot_wider(names_from=variable,values_from=value)
  
  return(out)
  
}

# load data
out_files = paste0("analysis/output/",dir(here("analysis","output")))
lapply(out_files,load,.GlobalEnv)

# choose models
models = mget(rev(ls()[str_detect(ls(),"out[.]")]))

# compute elasticities
end = nrow(models[[1]]$betadraws)
burn = 0.5*end
elasticities = elast(models,burn:end,p)

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

# cross
elasticities %>%
  filter(!own) %>%
  mutate(beta = case_when(beta == "ridge" ~ paste("beta","-ridge"),
                          beta == "lasso" ~ paste("beta","-lasso"),
                          beta == "horseshoe" ~ paste("beta","-horseshoe")),
         theta = case_when(theta == "ridge" ~ paste("theta","-ridge"),
                           theta == "horseshoe" ~ paste("theta","-horseshoe"),
                           theta == "sparse" ~ paste("theta","-sparse"))) %>%
  mutate(beta = factor(beta, levels=unique(beta)),
         theta = factor(theta, levels=unique(theta))) %>%
  # mutate(significant = ifelse(lower<0 & upper>0,"no","yes")) %>%
  # ggplot(aes(x=mean, fill=significant)) + 
  ggplot(aes(x=mean)) +
  geom_histogram(alpha=0.75) +
  xlim(-0.25,0.25) +
  labs(x="",y="") +
  theme_minimal() +
  facet_grid(theta~beta, labeller = label_parsed, scales="free_y") +
  theme(plot.title = element_text(face="bold",hjust=0.5))
ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/elasticities-cross.png",
       height=5,width=7)

# own
elasticities %>%
  filter(own) %>%
  mutate(beta = case_when(beta == "ridge" ~ paste("beta","-ridge"),
                          beta == "lasso" ~ paste("beta","-lasso"),
                          beta == "horseshoe" ~ paste("beta","-horseshoe")),
         theta = case_when(theta == "ridge" ~ paste("theta","-ridge"),
                           theta == "horseshoe" ~ paste("theta","-horseshoe"),
                           theta == "sparse" ~ paste("theta","-sparse"))) %>%
  mutate(beta = factor(beta, levels=unique(beta)),
         theta = factor(theta, levels=unique(theta))) %>%
  ggplot(aes(x=mean)) +
  geom_histogram(alpha=0.75) +
  xlim(-4,2) +
  facet_grid(theta~beta, labeller = label_parsed, scales="free") +
  labs(x="",y="") +
  theme_minimal() +
  theme(plot.title = element_text(face="bold",hjust=0.5))
ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/elasticities-own.png",
       height=5,width=7)

# pairwise own
elasticities %>%
  filter(own) %>%
  mutate(model = paste0(theta,"\n",beta)) %>%
  select(id,model,mean) %>%
  pivot_wider(names_from=model,values_from=mean) %>% 
  select(-id) %>%
  ggpairs(aes(alpha = 0.4)) +
  theme_minimal()
ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/elasticities-pairwise.png",
       height=10,width=10)













# effects of hierarchical shrinkage

qtlow = quantile(apply(prices[inweeks,-(1:2)],2,sd),0.2)
qthigh = quantile(apply(prices[inweeks,-(1:2)],2,sd),0.8)
sds = data.frame(id=1:p^2,
                 Var1=rep(1:p,each=p),
                 Var2=rep(1:p,times=p),
                 sd=rep(apply(prices[inweeks,-(1:2)],2,sd),times=p))
elasticities %>%
  left_join(sds) %>%
  mutate(sd = ifelse(sd<qtlow,"10",ifelse(sd>qthigh,"90","10-90"))) %>% 
  filter(own) %>%
  select(id,beta,theta,mean,sd) %>%
  pivot_wider(names_from=theta,values_from=mean) %>%
  pivot_longer(ridge:horseshoe,names_to="theta",values_to="mean") %>%
  ggplot(aes(x=sparse,y=mean)) +
  geom_point(data = . %>% filter(sd=="10-90"),aes(color=sd,shape=sd),size=1) +
  geom_point(data = . %>% filter(sd=="90"),aes(color=sd,shape=sd),size=1) +
  geom_point(data = . %>% filter(sd=="10"),aes(color=sd,shape=sd),size=1) +
  geom_abline(slope=1,intercept=0) +
  scale_color_manual(values=c("#F8766D","grey90","#00BFC4")) +
  scale_shape_manual(values=c(17,16,15)) +
  xlim(-3,1) +
  labs(x="beta from sparse model",y="beta from hierarchical model") +
  facet_grid(theta~beta) +
  theme_minimal()
ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/elasticities-beta-own.png",
       height=5,width=8)

# qts = quantile(sds %>% filter(Var1==Var2) %>% pull(sd),seq(0,1,by=0.1))
# qts = quantile(sds %>% filter(Var1>Var2) %>% pull(sd),seq(0,1,by=0.1))
qtlow = quantile(apply(prices[inweeks,-(1:2)],2,sd),0.01)
qthigh = quantile(apply(prices[inweeks,-(1:2)],2,sd),0.99)
elasticities %>%
  left_join(sds) %>%
  mutate(sd = ifelse(sd<qtlow,"10",ifelse(sd>qthigh,"90","10-90"))) %>%
  # mutate(sd = cut(sd,breaks=qts,include.lowest=TRUE,labels=names(qts)[-1])) %>%
  filter(!own) %>%
  select(id,beta,theta,mean,sd) %>%
  pivot_wider(names_from=theta,values_from=mean) %>%
  pivot_longer(ridge:horseshoe,names_to="theta",values_to="mean") %>%
  ggplot(aes(x=sparse,y=mean)) +
  geom_point(data = . %>% filter(sd=="10-90"),aes(color=sd,shape=sd),size=1) +
  geom_point(data = . %>% filter(sd=="90"),aes(color=sd,shape=sd),size=1) +
  geom_point(data = . %>% filter(sd=="10"),aes(color=sd,shape=sd),size=1) +
  geom_abline(slope=1,intercept=0) +
  # scale_colour_hue(h = c(270, 360)) +
  scale_color_manual(values=c("#F8766D","grey90","#00BFC4")) +
  scale_shape_manual(values=c(17,16,15)) +
  labs(x="beta from sparse model",y="beta from hierarchical model") +
  xlim(-0.2,0.2) +
  ylim(-0.2,0.2) +
  facet_grid(theta~beta) +
  theme_minimal()
ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/elasticities-beta-cross.png",
       height=5,width=8)
