library(scales)

# ----------------------------------------------------------------- #
# exact
# ----------------------------------------------------------------- #

# function to compute summaries of price elasticities
summarize_shrinkage = function(out_list,index,prices){
  
  X = prices %>% 
    select(starts_with("PRICE_")) %>%
    as.matrix()
  p = ncol(X)
  n = nrow(X)
  XpX = crossprod(X)
  # XpX = n*diag(p)
  
  kappa = matrix(0,p^2,length(out_list))
  tausq = double(p)
  
  for(m in 1:length(out_list)){
    
    for(i in 1:p){
      
      kappatmp = matrix(0,length(index),p)
      j = 1
      for(r in index){
        
        Psi = double(p)
        if(is.null(out_list[[m]]$Psidraws)){
          Psi[i] = out_list[[m]]$lambdasqowndraws[r,i]
          Psi[-i] = out_list[[m]]$lambdasqdraws[r,((p-1)*i-(p-1)+1):((p-1)*i)]
        }
        else{
          Psi[i] = out_list[[m]]$Psiowndraws[r,i]
          Psi[-i] = out_list[[m]]$Psidraws[r,((p-1)*i-(p-1)+1):((p-1)*i)]
        }
        
        tausq[i] = out_list[[m]]$tausqowndraws[r,1]
        tausq[-i] = out_list[[m]]$tausqdraws[r,1]
        
        sigmasq = out_list[[m]]$sigmasqdraws[r,i]
        
        Lambdastar = Psi*tausq
        
        XpXi = chol2inv(chol(XpX/sigmasq + diag(1/Lambdastar)))
        kappatmp[j,] = diag(XpXi%*%diag(1/Lambdastar))
        j = j+1
      }
      
      kappa[(p*i-p+1):(p*i),m] = apply(kappatmp,2,median)
      
    }
    
  }
  
  # format output
  cols = str_remove(names(out_list),"out.")
  colnames(kappa) = paste("kappa",cols,sep="_")
  kappaout = data.frame(kappa) %>%
    mutate(id = rep(1:p,each=p),
           index = 1:p^2,
           own = index %in% c(diag(matrix(1:p^2,p,p)))) %>%
    pivot_longer(-c(id,index,own), names_to = c("variable","theta","beta"), 
                 values_to = "kappa", names_pattern="(.*)_(.*)[.](.*)") %>%
    mutate(model=paste(theta,beta,sep="_"))
  
  return(kappaout)
  
}


tmp = summarize_shrinkage(models,(0.95*end):end,prices)
tmp %>% 
  filter(own) %>%
  group_by(model,own) %>%
  summarise(mean(kappa))

# ----------------------------------------------------------------- #
# approximation
# ----------------------------------------------------------------- #

# function to compute summaries of price elasticities
summarize_shrinkage = function(out_list,index,prices,own=TRUE){

  tausqout = double(length(out_list))
  X = prices %>%
    select(starts_with("PRICE_")) %>%
    as.matrix()
  p = ncol(X)
  n = nrow(X)
  varx = apply(X,2,var)

  # own elasticities
  if(own){

    kappa = matrix(0,p,length(out_list))
    ssq = matrix(varx,length(index),p,byrow=TRUE)
    for(m in 1:length(out_list)){
      if(is.null(out_list[[m]]$Psiowndraws)){
        Psi = out_list[[m]]$lambdasqowndraws[index,1:p]
      }
      else{
        Psi = out_list[[m]]$Psiowndraws[index,1:p]
      }
      sigmasq = out_list[[m]]$sigmasqdraws[index,]
      tausq = matrix(out_list[[m]]$tausqowndraws[index,1],length(index),p)
      kappa[,m] = apply(1/(1+n*ssq*Psi*tausq/sigmasq),2,median)
      tausqout[m] = median(out_list[[m]]$tausqowndraws[index,1])
    }
    cols = str_remove(names(out_list),"out.")
    colnames(kappa) = paste("kappa",cols,sep="_")
    kappaout = data.frame(kappa) %>%
      mutate(id = 1:p) %>%
      pivot_longer(-c(id), names_to = c("variable","theta","beta"), names_pattern="(.*)_(.*)[.](.*)") %>%
      mutate(model=paste(theta,beta,sep="_")) %>%
      pivot_wider(names_from=variable,values_from=value)

  }
  
  # cross elasticities
  else{

    kappa = matrix(0,p^2-p,length(out_list))
    ssq = matrix(rep(varx,times=p-1),length(index),p^2-p,byrow=TRUE)
    for(m in 1:length(out_list)){
      if(is.null(out_list[[m]]$Psidraws)){
        Psi = out_list[[m]]$lambdasqdraws[index,1:(p^2-p)]
      }
      else{
        Psi = out_list[[m]]$Psidraws[index,1:(p^2-p)]
      }
      sigmasq = out_list[[m]]$sigmasqdraws[index,]
      sigmasq = matrix(rep(sigmasq,each=p-1),length(index),p^2-p,byrow=TRUE)
      tausq = matrix(out_list[[m]]$tausqdraws[index,1],length(index),p^2-p)
      kappa[,m] = apply(1/(1+n*ssq*Psi*tausq/sigmasq),2,median)
      tausqout[m] = median(out_list[[m]]$tausqdraws[index,1])
    }
    
    # format output
    cols = str_remove(names(out_list),"out.")
    colnames(kappa) = paste("kappa",cols,sep="_")
    kappaout = data.frame(kappa) %>%
      mutate(id = 1:(p^2-p)) %>%
      pivot_longer(-c(id), names_to = c("variable","theta","beta"), names_pattern="(.*)_(.*)[.](.*)") %>%
      mutate(model=paste(theta,beta,sep="_")) %>%
      pivot_wider(names_from=variable,values_from=value)

  }

  return(list(kappa=kappaout,
              tausq=data.frame(model=str_replace(str_remove(names(models),"out."),"\\.","_"),tausqout)))

}

debug(summarize_shrinkage)

own = summarize_shrinkage(models,burn:end,prices,own=TRUE)
tbl_own = own$kappa %>%
  group_by(model) %>%
  summarise(own_min=min(kappa),own_mean=mean(kappa),own_max=max(kappa)) %>%
  left_join(own$tausq) %>%
  rename("own_tausq"="tausqout") %>%
  relocate(own_tausq, .before=own_min)

cross = summarize_shrinkage(models,burn:end,prices,own=FALSE)
tbl_cross = cross$kappa %>%
  group_by(model) %>%
  summarise(cross_min=min(kappa),cross_mean=mean(kappa),cross_max=max(kappa)) %>%
  left_join(cross$tausq) %>%
  rename("cross_tausq"="tausqout") %>%
  relocate(cross_tausq, .before=cross_min)

tbl = tbl_own %>%
  left_join(tbl_cross) %>%
  arrange(-row_number())
print(xtable(tbl, digits=c(0,0,2,2,2,2,-2,2,2,2)), include.rownames = FALSE)

# by category
formatout(own$kappa, "plot", "wrap") %>%
  left_join(tree_names_clean,by=c("id"="PRODUCT")) %>%
  mutate(model = factor(model, levels=unique(model)),
         LARGE_CATEGORY = str_remove(LARGE_CATEGORY,"/ALCOHOLIC CID")) %>%
  ggplot(aes(x=kappa,color=model,linetype=model)) +
  stat_ecdf() +
  labs(x="shrinkage factor (kappa)",y="") +
  facet_wrap(vars(LARGE_CATEGORY)) +
  scale_color_manual(values=hue_pal()(9),labels = expression(beta-Ridge,
                                         beta-Lasso,
                                         beta-Horseshoe,
                                         theta-Ridge~~beta-Ridge,
                                         theta-Ridge~~beta-Lasso,
                                         theta-Ridge~~beta-Horseshoe,
                                         theta-Horseshoe~~beta-Ridge,
                                         theta-Horseshoe~~beta-Lasso,
                                         theta-Horseshoe~~beta-Horseshoe)) +
  scale_linetype_manual(values=c(2,2,2,1,1,1,1,1,1),
                        labels = expression(beta-Ridge,
                                            beta-Lasso,
                                            beta-Horseshoe,
                                            theta-Ridge~~beta-Ridge,
                                            theta-Ridge~~beta-Lasso,
                                            theta-Ridge~~beta-Horseshoe,
                                            theta-Horseshoe~~beta-Ridge,
                                            theta-Horseshoe~~beta-Lasso,
                                            theta-Horseshoe~~beta-Horseshoe)) +
  guides(color=guide_legend(nrow=3, title=""),
         linetype=guide_legend(nrow=3, title="")) +
  theme_shrink() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.text.align = 0,
        legend.text = element_text(size=8))
ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/shrinkage-category.png",
       height=7.5,width=6)

# # ----------------------------------------------------------------- #
# # by category
# # ----------------------------------------------------------------- #
# 
# # function to compute summaries of price elasticities
# summarize_shrinkage = function(out_list,index,p,tree_names){
#   
#   kappa = matrix(0,p,length(out_list))
#   for(m in 1:length(out_list)){
#     if(is.null(out_list[[m]]$Psiowndraws)){
#       Psi = out_list[[m]]$lambdasqowndraws[index,1:p]
#     }
#     else{
#       Psi = out_list[[m]]$Psiowndraws[index,1:p]
#     }
#     tausq = matrix(out_list[[m]]$tausqowndraws[index,1],length(index),p)
#     # kappa[,m] = apply(1/(1+Psi*tausq),2,median)
#     kappa[,m] = apply(1/(1+Psi),2,median)
#   }
#   cols = str_remove(names(out_list),"out.")
#   colnames(kappa) = paste("kappa",cols,sep="_")
# 
#   out = data.frame(kappa) %>%
#     mutate(id = 1:p,
#            LARGE_CATEGORY = tree_names$LARGE_CATEGORY) %>%
#     pivot_longer(-c(id,LARGE_CATEGORY), names_to = c("variable","theta","beta"), names_pattern="(.*)_(.*)[.](.*)") %>%
#     mutate(model=paste(theta,beta,sep="_")) %>%
#     pivot_wider(names_from=variable,values_from=value)
#   
#   # out = data.frame(kappa) %>%
#   #   mutate(id = rep(1:p,times=length(index)),
#   #          ind = rep(1:length(index),each=p),
#   #          LARGE_CATEGORY = rep(tree_names$LARGE_CATEGORY, times=length(index))) %>%
#   #   pivot_longer(-c(id,ind,LARGE_CATEGORY), names_to = c("variable","theta","beta"), names_pattern="(.*)_(.*)[.](.*)") %>%
#   #   mutate(model=paste(theta,beta,sep="_")) %>%
#   #   pivot_wider(names_from=variable,values_from=value)
#   
#   return(out)
#   
# }
# 
# debug(summarize_shrinkage)
# 
# 
# summarize_shrinkage(models,burn:end,p,tree_names_clean) %>%
#   group_by(model) %>%
#   summarise(mean(kappa))
# 
# 
# 
# # beta-lasso
# kappas.sparse.lasso = summarize_shrinkage(list(out.sparse.lasso=out.sparse.lasso),burn:end,p,tree_names_clean) %>%
#   group_by(model,LARGE_CATEGORY) %>%
#   summarise(kappa_min = min(kappa),
#             kappa_median = mean(kappa),
#             kappa_max = max(kappa)) %>%
#   pivot_wider(names_from="model",values_from = starts_with("kappa")) 
# 
# # beta-horseshoe
# kappas.sparse.horseshoe = summarize_shrinkage(list(out.sparse.horseshoe=out.sparse.horseshoe),burn:end,p,tree_names_clean) %>%
#   group_by(model,LARGE_CATEGORY) %>%
#   summarise(kappa_min = min(kappa),
#             kappa_median = mean(kappa),
#             kappa_max = max(kappa)) %>%
#   pivot_wider(names_from="model",values_from = starts_with("kappa")) 
# 
# # theta-ridge, beta-lasso
# kappas.ridge.lasso = summarize_shrinkage(list(out.ridge.lasso=out.ridge.lasso),burn:end,p,tree_names_clean) %>%
#   group_by(model,LARGE_CATEGORY) %>%
#   summarise(kappa_min = min(kappa),
#             kappa_mean = mean(kappa),
#             kappa_max = max(kappa)) %>%
#   pivot_wider(names_from="model",values_from = starts_with("kappa")) 
# 
# # theta-ridge, beta-horseshoe
# kappas.ridge.horseshoe = summarize_shrinkage(list(out.ridge.horseshoe=out.ridge.horseshoe),burn:end,p,tree_names_clean) %>%
#   group_by(model,LARGE_CATEGORY) %>%
#   summarise(kappa_min = min(kappa),
#             kappa_mean = mean(kappa),
#             kappa_max = max(kappa)) %>%
#   pivot_wider(names_from="model",values_from = starts_with("kappa")) 
# 
# # theta-horseshoe, beta-ridge
# kappas.horseshoe.ridge = summarize_shrinkage(list(out.horseshoe.ridge=out.horseshoe.ridge),burn:end,p,tree_names_clean) %>%
#   group_by(model,LARGE_CATEGORY) %>%
#   summarise(kappa_min = min(kappa),
#             kappa_mean = mean(kappa),
#             kappa_max = max(kappa)) %>%
#   pivot_wider(names_from="model",values_from = starts_with("kappa")) 
# 
# # theta-horseshoe, beta-horseshoe
# kappas.horseshoe.horseshoe = summarize_shrinkage(list(out.horseshoe.horseshoe=out.horseshoe.horseshoe),burn:end,p,tree_names_clean) %>%
#   group_by(model,LARGE_CATEGORY) %>%
#   summarise(kappa_min = min(kappa),
#             kappa_mean = mean(kappa),
#             kappa_max = max(kappa)) %>%
#   pivot_wider(names_from="model",values_from = starts_with("kappa")) 
# 
# # print
# kappas = kappas.sparse.horseshoe %>%
#   left_join(kappas.ridge.horseshoe) %>%
#   left_join(kappas.horseshoe.ridge)
# print(xtable(kappas), include.rownames = FALSE)
