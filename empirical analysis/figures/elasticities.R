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

# choose models
out_list = mget(rev(ls()[str_detect(ls(),"out[.]")]))

# compute elasticities
end = Mcmc$R/Mcmc$keep
burn = 0.5*end
elasticities = elast(out_list,burn:end,p)

# summary stats (own)
elasticities %>%
  filter(own) %>%
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
