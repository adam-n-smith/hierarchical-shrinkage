# --------------------------------------------------------- #
# predictive fit statistics
# --------------------------------------------------------- #

# function to compute rmse
compute_rmse = function(out_list,index,n,p,Clisttest,cumnphi){

  rmse = matrix(0,length(index),length(out_list))
  i = 1
  for(r in index){
    for(m in 1:length(out_list)){
      B = matrix(out_list[[m]]$betadraws[r,],p,p)
      Cphi = matrix(0,n,p)
      for(j in 1:p){
        Cphi[,j] = Clisttest[[j]]%*%out_list[[m]]$phidraws[r,(cumnphi[j]+1):cumnphi[j+1]]
      }
      # rmse[i,m] = mean((Ytest - Xtest%*%B - Cphi)^2)
      rmse[i,m] = sqrt(mean((Ytest - Xtest%*%B - Cphi)^2))
    }
    i = i + 1
  }
  colnames(rmse) = str_remove(names(out_list),"out.")
  return(data.frame(rmse))

}

# load data
out_files = paste0("empirical analysis/output/",dir(here("empirical analysis","output")))
out_files = out_files[!str_detect(out_files,"old")]
lapply(out_files,load,.GlobalEnv)
# out_list = lapply(out_files,function(x)unlist(mget(load(x)),recursive=FALSE))

# choose models
out_list = mget(rev(ls()[str_detect(ls(),"out[.]")]))

# compute rmse
end = Mcmc$R/Mcmc$keep
burn = 0.5*end
rmse = compute_rmse(out_list,burn:end,ntest,p,Clisttest,cumnphi)

# plot
rmse %>%
  pivot_longer(everything()) %>%
  mutate(name = factor(name,levels=unique(name))) %>%
  group_by(name) %>%
  mutate(position = round(quantile(value,0.8),3),
         label = round(mean(value),3)) %>%
  ggplot(aes(x=name,y=value))  +
  geom_boxplot() + 
  geom_text(aes(x=name,y=position,label=label))+
  labs(x="",y="") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90,hjust=1))

# means
mean = apply(rmse,2,function(x)sprintf("%.3f",mean(x)))
sd = apply(rmse,2,function(x)sprintf("%.4f",sd(x)))
cbind(mean,sd)

# matplot
matplot(rmse,type="l")
