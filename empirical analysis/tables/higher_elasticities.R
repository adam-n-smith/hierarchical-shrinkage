end = Mcmc$R/Mcmc$keep
burn = 0.5*end

out = out.ridge.ridge

# --------------------------------------------------------- #
# heatmap
# --------------------------------------------------------- #

# CATEGORY level
theta = matrix(apply(out$thetadraws[burn:end,1:npar[1]],2,mean),sqrt(npar[1]),sqrt(npar[1]))
colnames(theta) = unique(tree_names_clean[,1]) %>% pull()
rownames(theta) = unique(tree_names_clean[,1]) %>% pull()
reshape2::melt(theta) %>%
  mutate(Var1=as.factor(Var1),
         Var2=as.factor(Var2)) %>%
  ggplot(aes(x=Var1,y=fct_rev(Var2),fill=value)) +
  geom_tile(color="white") +
  labs(x="",y="") +
  scale_x_discrete(position = "top") +
  scale_fill_gradient2(low="#F8766D",mid="white",high="#00BFC4") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90,hjust=0))
ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/heatmap-cat.png",
       height=7,width=7)

# SUBCATEGORY level
theta = matrix(apply(out$thetadraws[burn:end,(cumsum(npar)[1]+1):cumsum(npar)[2]],2,mean),sqrt(npar[2]),sqrt(npar[2]))
colnames(theta) = unique(tree_names_clean[,2]) %>% pull() %>% to_any_case(case="title", parsing_option = 0)
rownames(theta) = unique(tree_names_clean[,2]) %>% pull() %>% to_any_case(case="title", parsing_option = 0)
reshape2::melt(theta) %>%
  mutate(Var1=as.factor(Var1),
         Var2=as.factor(Var2)) %>%
  ggplot(aes(x=Var1,y=fct_rev(Var2),fill=value)) +
  geom_tile(color="white") +
  labs(x="",y="") +
  scale_x_discrete(position = "top") +
  scale_fill_gradient2(low="#F8766D",mid="white",high="#00BFC4") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90,hjust=0))
ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/heatmap-subcat.png",
       height=7,width=7)


# --------------------------------------------------------- #
# own vs. cross category effects
# --------------------------------------------------------- #

# CATEGORY level
df = data.frame(theta = round(apply(out$thetadraw[burn:end,1:npar[1]],2,mean),3),
                SD = round(apply(out$thetadraw[burn:end,1:npar[1]],2,sd),3),
                i = rep(unique(product_table$LARGE_CATEGORY),n_distinct(product_table$LARGE_CATEGORY)),
                j = rep(unique(product_table$LARGE_CATEGORY),each=n_distinct(product_table$LARGE_CATEGORY))) %>%
  mutate(within = ifelse(i==j,"within category","across category"))
ggplot(df,aes(x=theta)) +
  geom_density(aes(color=within,linetype=within)) +
  labs(x="", y="density\n") +
  xlim(-0.2,0.2) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust=0.5, face="bold"))
ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/elasticities-cat.png",
       height=3,width=4)


# SUBCATEGORY level
df = data.frame(theta = round(apply(out$thetadraw[burn:end,(npar[1]+1):cumsum(npar)[2]],2,mean),3),
                lower = round(apply(out$thetadraw[burn:end,(npar[1]+1):cumsum(npar)[2]],2,function(x)quantile(x,0.025)),3),
                upper = round(apply(out$thetadraw[burn:end,(npar[1]+1):cumsum(npar)[2]],2,function(x)quantile(x,0.975)),3),
                i = rep(unique(product_table$SMALL_CATEGORY),n_distinct(product_table$SMALL_CATEGORY)),
                j = rep(unique(product_table$SMALL_CATEGORY),each=n_distinct(product_table$SMALL_CATEGORY))) %>%
  mutate(within = ifelse(i==j,"within subcategory","across subcategory"),
         significant = !(lower<0 & upper>0))
ggplot(df,aes(x=theta)) +
  geom_density(aes(color=within, linetype=within)) +
  labs(x="", y="density\n") +
  xlim(-0.2,0.2) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust=0.5, face="bold"))
ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/elasticities-subcat.png",
       height=3,width=4)

# --------------------------------------------------------- #
# thetas (own)
# --------------------------------------------------------- #

# df = data.frame(theta = round(apply(out$thetaowndraws[burn:end,1:npar_own[1]],2,mean),3),
#                 SD = round(apply(out$thetaowndraws[burn:end,1:npar_own[1]],2,sd),3))
# ggplot(df,aes(x=theta)) + 
#   geom_density(color="grey") + 
#   geom_density(data=df[df$focal,],aes(x=theta),color="#F8766D",size=1.2) +
#   labs(x="", y="density\n") +
#   theme_minimal() +
#   theme(legend.position="none",
#         plot.title = element_text(hjust=0.5))
# 
# df = data.frame(theta = round(apply(out$thetaowndraws[burn:end,1:npar_own[1]],2,mean),3),
#                 SD = round(apply(out$thetaowndraws[burn:end,1:npar_own[1]],2,sd),3))
# ggplot(df,aes(x=theta)) + 
#   geom_density(color="grey") + 
#   geom_density(data=df[df$focal,],aes(x=theta),color="#F8766D",size=1.2) +
#   labs(x="", y="density\n") +
#   theme_minimal() +
#   theme(legend.position="none",
#         plot.title = element_text(hjust=0.5))

# --------------------------------------------------------- #
# theta tables
# --------------------------------------------------------- #

end = Mcmc$R/Mcmc$keep
burn = 0.5*end

df = data.frame(mean = apply(out$thetadraw[burn:end,1:npar[1]],2,mean),
                sd = apply(out$thetadraw[burn:end,1:npar[1]],2,sd),
                lower = apply(out$thetadraw[burn:end,1:npar[1]],2,function(x)quantile(x,0.025)),
                upper = apply(out$thetadraw[burn:end,1:npar[1]],2,function(x)quantile(x,0.975)),
                ell = rep(unique(product_table$LARGE_CATEGORY),n_distinct(product_table$LARGE_CATEGORY)),
                k = rep(unique(product_table$LARGE_CATEGORY),each=n_distinct(product_table$LARGE_CATEGORY))) %>%
  mutate(focal = (k=="BEER/ALE/ALCOHOLIC CID"),
         within = k==ell)
tbl = df %>% 
  filter(focal) %>%
  arrange(-mean) %>%
  rownames_to_column() %>%
  mutate(mean = sprintf("%.3f",mean),
         sd = sprintf("%.4f",sd),
         # mean = ifelse(lower<0 & upper>0, mean, paste("\\bf",mean)),
         k = ifelse(rowname==1,k,"")) %>%
  select(k,ell,mean,sd)
print(xtable(tbl),include.rownames = F, sanitize.text.function = identity)

# subcategory
df = data.frame(mean = apply(out$thetadraw[burn:end,(npar[1]+1):cumsum(npar)[2]],2,mean),
                sd = apply(out$thetadraw[burn:end,(npar[1]+1):cumsum(npar)[2]],2,sd),
                lower = apply(out$thetadraw[burn:end,(npar[1]+1):cumsum(npar)[2]],2,function(x)quantile(x,0.025)),
                upper = apply(out$thetadraw[burn:end,(npar[1]+1):cumsum(npar)[2]],2,function(x)quantile(x,0.975)),
                i = rep(unique(tree_names_clean$SMALL_CATEGORY),n_distinct(tree_names_clean$SMALL_CATEGORY)),
                j = rep(unique(tree_names_clean$SMALL_CATEGORY),each=n_distinct(tree_names_clean$SMALL_CATEGORY))) %>%
  left_join(.,unique(tree_names_clean[,c("SMALL_CATEGORY","LARGE_CATEGORY")]),by=c("j"="SMALL_CATEGORY")) %>%
  mutate(focal_i = str_detect(j,"Imported"),
         focal_j = str_detect(j,"Domestic"),
         within =i==j)
tbl = df %>% 
  filter(focal_j) %>%
  arrange(-mean) %>%
  group_by(mean>0) %>%
  mutate(rank=order(abs(mean),decreasing=TRUE)) %>%
  filter(rank<=10) %>%
  ungroup() %>%
  rownames_to_column() %>%
  mutate(mean = sprintf("%.3f",mean),
         sd = sprintf("%.4f",round(sd,4)),
         # mean = ifelse(lower<0 & upper>0, mean, paste("\\bf",mean)),
         # i = to_any_case(i,case="title", parsing_option = 0),
         # j = to_any_case(j,case="title", parsing_option = 0),
         j = ifelse(rowname==1,j,"")) %>%
  select(j,i,mean,sd)
print(xtable(tbl),include.rownames = F, sanitize.text.function = identity)
