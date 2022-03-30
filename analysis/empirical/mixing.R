library(tidyverse)
library(reshape2)
library(coda)

source(here("src","shrinkage-functions.R"))

getacf = function(models,lags,quantiles,burn,wchpar){
  
  # omit sparse models when analyzing thetas
  if(all(wchpar!="beta") & any(str_detect(names(models),"sparse"))){
    models = models[-which(str_detect(names(models),"sparse"))]
  }

  # number of models
  m = length(models)
  
  # format model names
  model_names = names(models) %>% 
    str_remove(., "out.") %>%
    str_replace(.,"\\.","_") 

  # storage
  R = nrow(models[[1]]$betadraws[-(1:burn),])
  lags = lags[lags<R]
  out = array(dim=c(m,length(quantiles),length(lags)), 
              dimnames=list(model=model_names,quantile=quantiles,lag=lags))
  
  # start loop over models
  pb = txtProgressBar(min=0, max=m, initial=1,style=3) 
  for(i in 1:m){
    
    #  format draws (check if parameter is beta or indices for theta)
    if(all(wchpar=="beta")){
      nbeta = ncol(models[[i]]$betadraws)
      if(nbeta>1000){
        set.seed(1)
        cols = sample(1:nbeta,round(0.01*nbeta))
      }
      else{
        cols = 1:nbeta
      }
      draws = mcmc(models[[i]]$betadraws[-(1:burn),cols])
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
    mutate(quantile = factor(quantile,levels=rev(quantiles))) %>%
    # separate model into theta and beta
    mutate(theta = word(model,1,sep="_"),
           beta = word(model,2,sep="_")) %>%
    select(-model) %>%
    # formatting
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

figurepath = "/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/"

# load data
load(here("build","output","store_panel.RData"))

# build parameter index vector
treeindex = createindex(tree)
npar = treeindex$npar
npar_own = treeindex$npar_own

# load models
models = mget(ls()[str_detect(ls(),"out.")])

# specs
end = nrow(models[[1]]$betadraws)
burn = 0.5*end
quantiles = c(0.1,0.5,0.9)
lags = seq(0,50,by=1)
colors = RColorBrewer::brewer.pal(8,"Blues")[c(4,6,8)]

# thetas (upper level)
matplot(out.ridge.ridge$thetadraws[,1],type="l")
acf.theta.top = getacf(models,lags,quantiles,burn,wchpar=1:npar[1])
acf.theta.top %>%
  ggplot(aes(x=lag,y=acf,group=quantile)) +
  geom_line(aes(color=quantile)) +
  ylim(min(acf.theta.top$acf),1) +
  facet_grid(theta~beta, labeller=label_parsed) +
  scale_color_manual(values=colors) +
  theme_minimal()
ggsave(filename=paste0(figurepath,"mixing-acf-top.png"),
       height=3,width=6)

# thetas (middle level)
matplot(out.ridge.ridge$thetadraws[,npar[1]+1],type="l")
acf.theta.mid = getacf(models,lags,quantiles,burn,wchpar=(npar[1]+1):cumsum(npar)[2])
acf.theta.mid %>%  
  ggplot(aes(x=lag,y=acf,group=quantile)) +
  geom_line(aes(color=quantile)) +
  ylim(min(acf.theta.mid$acf),1) +
  facet_grid(theta~beta, labeller=label_parsed) +
  scale_color_manual(values=colors) +
  theme_minimal()
ggsave(filename=paste0(figurepath,"mixing-acf-mid.png"),
       height=3,width=6)

# betas
matplot(out.sparse.ridge$betadraws[,1],type="l")
acf.beta = getacf(models,lags,quantiles,burn,wchpar="beta")
acf.beta %>%
  ggplot(aes(x=lag,y=acf,group=quantile,color=quantile)) +
  geom_line(aes(color=quantile)) +
  ylim(min(acf.beta$acf),1) +
  scale_color_brewer(palette = "Blues") +
  facet_grid(theta~beta, labeller=label_parsed) +
  theme_minimal()
ggsave(filename=paste0(figurepath,"mixing-acf-beta.png"),
       height=4,width=6)
