library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(here)

source(here("analysis","empirical","estimation-build.R"))
source(here("src","summary-functions.R"))

# load models
out_files = paste0("analysis/output/",dir(here("analysis","output")))
out_files = out_files[str_detect(out_files,"out-")]
lapply(out_files,load,.GlobalEnv)
models = mget(rev(ls()[str_detect(ls(),"out[.]")]))

# --------------------------------------------------------- #
# predictive fit statistics
# --------------------------------------------------------- #

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

# bind together
rmsetable = cbind(var=rep(c("all","sd","new"),each=nrow(rmse)),
                  rbind(rmse,rmse_psd,rmse_new))

# plot
rmsetable %>%
  pivot_longer(-var) %>%
  mutate(name = factor(name,levels=unique(name))) %>%
  group_by(name,var) %>%
  mutate(position = round(quantile(value,0.8),3),
         label = round(mean(value),3)) %>%
  ggplot(aes(x=name,y=value))  +
  geom_boxplot() + 
  geom_text(aes(x=name,y=position,label=label),color="red")+
  facet_grid(var~.,scales="free") +
  labs(x="",y="") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90,hjust=1))

tbl = rmsetable %>%
  pivot_longer(-var) %>%
  group_by(name,var) %>%
  summarise(mean = mean(value),
            sd = sd(value), .groups="drop") %>%
  pivot_wider(names_from=var, values_from=c(mean,sd)) %>%
  mutate(theta=word(name,1,sep=fixed(".")),
         beta=word(name,2,sep=fixed("."))) %>%
  formatout(.,"table") %>%
  select(theta,beta,mean_all,sd_all,mean_new,sd_new,mean_sd,sd_sd) %>%
  arrange(rev(theta),rev(beta)) %>%
  mutate(beta = as.character(beta),
         beta = ifelse(beta=="$\\quad\\beta$-Ridge",paste0(theta,beta),beta)) %>%
  select(-theta) %>%
  mutate(across(starts_with("sd_"), ~paste0("(",sprintf("%.3f",.x),")")),
         across(starts_with("mean_"),function(x)ifelse(x==min(x),paste("\\bf",sprintf("%.3f",x)),sprintf("%.3f",x))))

print(xtable(tbl, digits=3), 
      include.rownames=FALSE, 
      sanitize.text.function=function(x){x})

