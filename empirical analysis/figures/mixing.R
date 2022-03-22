library(tidyverse)
library(reshape2)
library(coda)

getacf = function(models,lags,quantiles,burn,wchpar){
  
  # number of models
  m = length(models)
  
  # format model names
  model_names = names(models) %>% 
    str_remove(., "out.") %>%
    str_replace(.,"\\.","\n")
  
  # storage
  out = array(dim=c(m,length(quantiles),length(lags)), 
              dimnames=list(model=model_names,quantile=quantiles,lag=lags))
  
  # start loop over models
  pb = txtProgressBar(min=1, max=m, initial=1,style=3) 
  for(i in 1:m){
    
    #  format draws (check if parameter is beta or indices for theta)
    if(all(wchpar=="beta")){
      draws = mcmc(models[[i]]$betadraws[-(1:burn),])
    }
    else{
      draws = mcmc(models[[i]]$thetadraws[-(1:burn),wchpar])
    }
    
    # compute autocorrelation function 
    acf = autocorr.diag(draws,lags=lags)
    
    # compute quantiles
    out[i,,] = apply(acf,1,quantile,quantiles)
    setTxtProgressBar(pb,i)
  }
  
  close(pb)
  
  # reformat array as tbl
  out = melt(out, value.name="acf") %>%
    mutate(quantile = factor(quantile,levels=rev(quantiles)))
  
  return(out)
  
}


models = mget(ls()[str_detect(ls(),"out.")])

end = Mcmc$R/Mcmc$keep
burn = Mcmc$burn_pct*end
quantiles = c(0.1,0.5,0.9)
lags = seq(0,50,by=1)
colors = RColorBrewer::brewer.pal(8,"Blues")[c(4,6,8)]

# thetas (upper level)
acf.theta.top = getacf(models,lags,quantiles,burn,wchpar=1:npar[1])
acf.theta.top %>%
  ggplot(aes(x=lag,y=acf,group=quantile)) +
  geom_line(aes(color=quantile)) +
  ylim(-0.2,1) +
  facet_wrap(vars(model)) +
  scale_color_manual(values=colors) +
  theme_minimal()

# thetas (middle level)
acf.theta.mid = getacf(models,lags,quantiles,burn,wchpar=(npar[1]+1):cumsum(npar)[2])
acf.theta.mid %>%  
  ggplot(aes(x=lag,y=acf,group=quantile)) +
  geom_line(aes(color=quantile)) +
  ylim(-0.2,1) +
  facet_wrap(vars(model)) +
  scale_color_manual(values=colors) +
  theme_minimal()

# betas
acf.beta = getacf(models,lags,quantiles,burn,wchpar="beta")
acf.beta %>%
  mutate(quantile=as.character(quantile)) %>%
  ggplot(aes(x=lag,y=acf,group=quantile,color=quantile)) +
  geom_line(aes(color=quantile)) +
  ylim(-0.2,1) +
  # scale_colour_manual(palette="Blues") +
  scale_color_brewer(palette = "Blues") +
  facet_wrap(vars(model)) +
  theme_minimal()

