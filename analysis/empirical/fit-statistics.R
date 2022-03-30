sourceCpp(here("src","shrinkage-mcmc.cpp"))
source(here("analysis","empirical","estimation-build.R"))

# --------------------------------------------------------- #
# predictive fit statistics
# --------------------------------------------------------- #

# function to compute rmse
compute_rmse = function(out_list,index,n,p,Clisttest,cumnphi,wchprod=NULL){

  if(is.null(wchprod)){
    wchprod = 1:p
  }
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
      Yhat = Xtest%*%B + Cphi
      rmse[i,m] = sqrt(mean((Ytest[,wchprod] - Yhat[,wchprod])^2))
    }
    i = i + 1
  }
  colnames(rmse) = str_remove(names(out_list),"out.")
  return(data.frame(rmse))

}

# load data
out_files = paste0("analysis/output/",dir(here("analysis","output")))
out_files = out_files[!str_detect(out_files,"march")]
lapply(out_files,load,.GlobalEnv)

# choose models
models = mget(rev(ls()[str_detect(ls(),"out[.]")]))

# compute rmse
end = nrow(models[[1]]$betadraws)
burn = 0.5*end
rmse = compute_rmse(models,burn:end,ntest,p,Clisttest,cumnphi)

# products with limited price variation
sds = apply(prices[inweeks,-(1:2)],2,sd)
rmse_psd = compute_rmse(models,burn:end,ntest,p,Clisttest,cumnphi,wchprod=which(sds<quantile(sds,0.25)))

# products with new prices in test data
newlow = apply(prices[-inweeks,-(1:2)],2,min) < apply(prices[inweeks,-(1:2)],2,min)
newhigh = apply(prices[-inweeks,-(1:2)],2,max) > apply(prices[inweeks,-(1:2)],2,max)
newprice = newlow|newhigh
rmse_new = compute_rmse(models,burn:end,ntest,p,Clisttest,cumnphi,wchprod=which(newprice))

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
mean = apply(rmse_new,2,function(x)sprintf("%.3f",mean(x)))
sd = apply(rmse_new,2,function(x)sprintf("%.4f",sd(x)))
cbind(mean,sd)

# matplot
matplot(rmse,type="l")

rmse %>%
  mutate(iteration=1:nrow(rmse)) %>%
  pivot_longer(-iteration,names_to="model") %>%
  ggplot(aes(x=iteration,y=value)) +
  geom_line(aes(color=model)) +
  facet_wrap(vars(model),scales="free")



