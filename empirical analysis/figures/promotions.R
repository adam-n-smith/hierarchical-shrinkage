# function to compute elasticities
promos = function(out_list,index,p){
  
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
    pivot_wider(names_from=variable,values_from=value)
  
  return(out)
  
}

promotions = promos(out_list,burn:end,p)


# cross
promotions %>%
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
  # xlim(-0.25,0.25) +
  labs(x="",y="") +
  theme_minimal() +
  facet_grid(theta~beta, labeller = label_parsed, scales="free_y") +
  theme(plot.title = element_text(face="bold",hjust=0.5))
ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/promotions.png",
       height=5,width=7)
