library(coda)
library(GGally)
library(xtable)
library(reshape2)

# --------------------------------------------------------- #
# ggplot theme
# --------------------------------------------------------- #

theme_shrink = function(){
  
  theme_minimal() %+replace% 
  
  theme(panel.grid.minor = element_line(size = 0.2, color="white"), 
        panel.grid.major = element_line(size = 0.2, color="white"),
        panel.background = element_rect(fill="grey95", color="white"),
        axis.title.x = element_text(margin = margin(t=10, r=0, b=0, l=0)),
        axis.title.y = element_text(margin = margin(t=0, r=10, b=0, l=0), angle = 90))
  
}

# --------------------------------------------------------- #
# functions to summarize output
# --------------------------------------------------------- #

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
    pivot_wider(names_from=variable,values_from=value)
  
  return(out)
  
}

# function to compute summaries of price elasticities
summarize_elasticities_higher = function(out_list,index,p,wchpar){
  
  n = length(wchpar)
  mean = matrix(0,n,length(out_list))
  lower = matrix(0,n,length(out_list))
  upper = matrix(0,n,length(out_list))
  for(m in 1:length(out_list)){
    mean[,m] = apply(out_list[[m]]$thetadraws[index,wchpar],2,mean)
    lower[,m] = apply(out_list[[m]]$thetadraws[index,wchpar],2,function(x)quantile(x,0.025))
    upper[,m] = apply(out_list[[m]]$thetadraws[index,wchpar],2,function(x)quantile(x,0.975))
  }
  cols = str_remove(names(out_list),"out.")
  colnames(mean) = paste("mean",cols,sep="_")
  colnames(lower) = paste("lower",cols,sep="_")
  colnames(upper) = paste("upper",cols,sep="_")
  
  out = data.frame(cbind(mean,lower,upper)) %>%
    mutate(id = wchpar) %>%
    pivot_longer(-id, names_to = c("variable","theta","beta"), names_pattern="(.*)_(.*)[.](.*)") %>%
    mutate(model=paste(theta,beta,sep="_")) %>%
    pivot_wider(names_from=variable,values_from=value)
  
  return(out)
  
}

# function to compute summaries of promotion effects
summarize_promos = function(out_list,index,p){
  
  ids = NULL
  phinames = NULL
  for(i in 1:p){
    phinames = c(phinames,colnames(Clist[[i]]))
    ids = c(ids,rep(i,length(which(colnames(Clist[[i]]) %in% c("display","feature")))))
  }
  
  wch = which(phinames %in% c("display","feature"))
  estimates = matrix(nrow=length(wch),ncol=length(out_list))
  for(m in 1:length(out_list)){
    estimates[,m] = apply(out_list[[m]]$phidraws[index,wch],2,mean)
  }
  cols = str_remove(names(out_list),"out.")
  colnames(estimates) = paste(cols,sep="_")
  
  out = data.frame(estimates) %>%
    mutate(id=ids, variable = phinames[wch]) %>%
    pivot_longer(-c(id,variable),names_to = c("theta","beta"), names_pattern="(.*)[.](.*)") %>%
    mutate(model=paste(theta,beta,sep="_"))
  
  return(out)
  
}

# format model names
formatout = function(out, object, facet=NULL){
  
  if(object=="plot"){
    
    if(facet=="wrap"){
      out = out %>%
        mutate(theta = case_when(theta == "ridge" ~ "theta-Ridge",
                                 theta == "horseshoe" ~ "theta-Horseshoe",
                                 theta == "sparse" ~ "Standard~Shrinkage"),
               beta = case_when(beta == "ridge" ~ "beta-Ridge",
                                beta == "lasso" ~ "beta-Lasso",
                                beta == "horseshoe" ~ "beta-Horseshoe"),
               model = ifelse(theta == "Standard~Shrinkage", beta, paste0(theta,"~~",beta)),
               shrinkage = ifelse(theta == "Standard~Shrinkage", "Standard Shrinkage", "Hierarchical Shrinkage"),
               shrinkage = factor(shrinkage, levels=unique(shrinkage)),
               model = factor(model, levels=unique(model)),
               beta = factor(beta, levels=unique(beta)),
               theta = factor(theta, levels=unique(theta)))
      
    }
    if(facet=="grid"){
      out = out %>%
        mutate(theta = case_when(theta == "ridge" ~ "atop(textstyle(Hierarchical),theta-Ridge)",
                                 theta == "horseshoe" ~ "atop(textstyle(Hierarchical),theta-Horseshoe)",
                                 theta == "sparse" ~ "atop(textstyle(Standard),textstyle(Shrinkage))"),
               beta = case_when(beta == "ridge" ~ "beta-Ridge",
                                beta == "lasso" ~ "beta-Lasso",
                                beta == "horseshoe" ~ "beta-Horseshoe"),
               beta = factor(beta, levels=unique(beta)),
               theta = factor(theta, levels=unique(theta)))
    }
  }
  
  if(object=="table"){
    out = out %>%
      mutate(theta = case_when(theta == "ridge" ~ paste0("Hierarchical (","$\\theta$","-Ridge",")\\\\"),
                               theta == "horseshoe" ~ paste0("Hierarchical (","$\\theta$","-Horseshoe",")\\\\ "),
                               theta == "sparse" ~ paste0("Sparse\\\\ ")),
             theta = factor(theta, levels=unique(theta)),
             beta = case_when(beta == "ridge" ~ paste0("$\\quad\\beta$","-Ridge"),
                              beta == "lasso" ~ paste0("$\\quad\\beta$","-Lasso"),
                              beta == "horseshoe" ~ paste0("$\\quad\\beta$","-Horseshoe")),
             beta = factor(beta, levels=unique(beta)))
  }
  
  return(out)
  
}

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

# autocorrelation
getacf = function(models,lags,quantiles,burn,wchpar,own){
  
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

  # parameters
  if(all(wchpar=="beta")){
    nbeta = ncol(models[[1]]$betadraws)
    cols = sample(1:nbeta, min(1000,nbeta))
  }
  
  # start loop over models
  pb = txtProgressBar(min=0, max=m, initial=1,style=3) 
  for(i in 1:m){
    
    #  format draws (check if parameter is beta or indices for theta)
    if(all(wchpar=="beta") & !own){
      draws = mcmc(models[[i]]$betadraws[-(1:burn),cols])
    }
    else if(all(wchpar=="beta")){
      wchown = which(as.vector(diag(sqrt(nbeta)))==1)
      draws = mcmc(models[[i]]$betadraws[-(1:burn),wchown])
    }
    else if(!own){
      draws = mcmc(models[[i]]$thetadraws[-(1:burn),wchpar])
    }
    else{
      draws = mcmc(models[[i]]$thetaowndraws[-(1:burn),wchpar])
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
    rename(percentile=quantile)
  
  return(out)
  
}
