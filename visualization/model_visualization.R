library(tidyverse)

# --------------------------------------------------------- #
# visualization
# --------------------------------------------------------- #

# means
betahat.hier = apply(out_hierridge$betadraw[burn:end,],2,mean)
betahat = apply(out_ridge$betadraw[burn:end,],2,mean)
# sd
betasd.hier = apply(out_hierridge$betadraw[burn:end,],2,sd)
betasd = apply(out_ridge$betadraw[burn:end,],2,sd)
betasd = apply(out_lasso$betadraw[burn:end,],2,sd)
betasd = apply(out_horseshoe$betadraw[burn:end,],2,sd)

# compare betas
wchown = which(diag(p)==1)
wchcross = which(diag(p)==0)
plot(betahat[wchown],betahat.hier[wchown],
     xlab="horseshoe",ylab="hierarchical horseshoe")
abline(0,1)
plot(betahat[wchcross],betahat.hier[wchcross],
     xlab="horseshoe",ylab="hierarchical horseshoe")
abline(0,1)

# variance
plot(betasd[wchcross],betasd.hier[wchcross],
     xlab="horseshoe",ylab="hierarchical horseshoe")
abline(0,1)


# markups
mrkup.horse = betahat.horse[wchown]/(1+betahat.horse[wchown])
mrkup.hierhorse = betahat.hierhorse[wchown]/(1+betahat.hierhorse[wchown])
plot(mrkup.horse,mrkup.hierhorse,xlim=c(1,5),ylim=c(1,5))
abline(0,1)
hist(mrkup.horse - mrkup.hierhorse,xlim=c(-4,4),breaks=1000)
















# --------------------------------------------------------- #
# elasticity histograms
# --------------------------------------------------------- #

library(gridExtra)
wchcross = as.vector(diag(p)==0)
df = data.frame(ridge=apply(out_ridge$betadraw[burn:end,],2,mean),
                # lasso=apply(out_lasso$betadraw[burn:end,],2,mean),
                horseshoe=apply(out_horse$betadraw[burn:end,],2,mean),
                hierridge=apply(out_hierridge$betadraw[burn:end,],2,mean),
                hierhorse=apply(out_hierhorse$betadraw[burn:end,],2,mean),
                cross = ifelse(wchcross,"cross elasticity","own elasticity"))
labels = c(ridge = "ridge", lasso = "lasso",horseshoe = "horseshoe",
           hierridge = "hierarchical ridge",hierhorse = "hierarchical horseshoe")
a = df %>%
  pivot_longer(.,-cross,names_to="model") %>%
  mutate(model=factor(model,levels=c("ridge","lasso","horseshoe","hierridge","hierhorse")),
         cross=factor(cross,levels=c("own elasticity","cross elasticity"))) %>%
  filter(cross=="own elasticity") %>%
  ggplot(.,aes(x=value,y=..density..,group=model)) + 
  geom_histogram() + 
  facet_wrap(~model,
             scales="free",
             labeller = labeller(model = labels),
             nrow=1) + 
  labs(x="",y="density\n",title="Own-Price Elasticities") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust=0.5))
b = df %>%
  pivot_longer(.,-cross,names_to="model") %>%
  mutate(model=factor(model,levels=c("ridge","lasso","horseshoe","hierridge","hierhorse")),
         cross=factor(cross,levels=c("own elasticity","cross elasticity"))) %>%
  filter(abs(value)<10) %>%
  filter(cross=="cross elasticity") %>%
  ggplot(.,aes(x=value,y=..density..,group=model)) + 
  geom_histogram() + 
  facet_wrap(~model,
             scales="free",
             labeller = labeller(model = labels),
             nrow=1) + 
  labs(x="",y="density\n",title="Cross-Price Elasticities") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust=0.5))
grid.arrange(a,b)
# ggsave(file="figures/beta_histogram.png",
#        grid.arrange(a,b),
#        width = 11, height = 5)
df %>%
  pivot_longer(.,-cross,names_to="model") %>%
  mutate(model=factor(model,levels=c("ridge","lasso","horseshoe","hierridge","hierhorse")),
         cross=factor(cross,levels=c("own elasticity","cross elasticity"))) %>%
  filter(cross=="own elasticity") %>%
  group_by(model) %>%
  summarise(median(value))

# --------------------------------------------------------- #
# standard deviations
# --------------------------------------------------------- #

wchcross = as.vector(which(diag(p)==0))
df = data.frame(beta.ridge=apply(out_ridge$betadraw[burn:end,wchcross],2,sd),
                beta.lasso=apply(out_lasso$betadraw[burn:end,wchcross],2,sd),
                beta.horse=apply(out_horse$betadraw[burn:end,wchcross],2,sd),
                beta.hierridge=apply(out_hierridge$betadraw[burn:end,wchcross],2,sd))
a = ggplot(df,aes(x=beta.ridge,y=beta.hierridge)) + 
  geom_point(size=1, color="grey") + 
  geom_abline(intercept=0,slope=1) + 
  labs(x="\n ridge", y="hierarchical ridge\n") +
  theme_minimal()
b = ggplot(df,aes(x=beta.lasso,y=beta.hierridge)) + 
  geom_point(size=1, color="grey") + 
  geom_abline(intercept=0,slope=1) + 
  labs(x="\n lasso", y="hierarchical ridge\n") +
  theme_minimal()
c = ggplot(df,aes(x=beta.horse,y=beta.hierridge)) + 
  geom_point(size=1, color="grey") + 
  geom_abline(intercept=0,slope=1) + 
  labs(x="\n horseshoe", y="hierarchical ridge\n") +
  xlim(0,.5) +
  theme_minimal()
# grid.arrange(a,b,c, nrow=1)
ggsave(file="figures/beta_sd.png",
       grid.arrange(a,b,c, nrow=1),
       width = 8, height = 3, 
       device='png')

# --------------------------------------------------------- #
# thetas
# --------------------------------------------------------- #

# # thetahat1 = apply(out_hierhorse$thetadraw[burn:end,1:npar[1]],2,mean)
# # thetahat2 = apply(out_hierhorse$thetadraw[burn:end,(npar[1]+1):cumsum(npar)[2]],2,mean)
# thetahat1 = apply(out_hierridge$thetadraw[burn:end,1:npar[1]],2,mean)
# thetahat2 = apply(out_hierridge$thetadraw[burn:end,(npar[1]+1):cumsum(npar)[2]],2,mean)
# level = rep(c("top (categories)","middle (subcategories)"),npar[1:2])
# level = factor(level,levels=c("top (categories)","middle (subcategories)"))
# own = ifelse(c(as.vector(diag(10)==1),as.vector(diag(36)==1)),"within","across")
# own = factor(own,levels=c("within","across"))
# df = data.frame(theta = c(thetahat1,thetahat2),
#                 level = level,
#                 own = own)
# ggplot(df,aes(x=theta,color=level,fill=level)) + 
#   geom_density(alpha=.5) + 
#   facet_wrap(.~own) + 
#   # geom_vline(xintercept=0) +
#   labs(x="\n theta", y="density\n") +
#   theme_minimal() +
#   theme(legend.position = "bottom")
# dev.copy2pdf(file="figures/thetas.pdf",width=6,height=3)

# --------------------------------------------------------- #
# thetas
# --------------------------------------------------------- #

df = data.frame(theta = round(apply(out_hierridge$thetadraw[burn:end,1:npar[1]],2,mean),3),
                SD = round(apply(out_hierridge$thetadraw[burn:end,1:npar[1]],2,sd),3),
                i = rep(unique(product_table$LARGE_CATEGORY),n_distinct(product_table$LARGE_CATEGORY)),
                j = rep(unique(product_table$LARGE_CATEGORY),each=n_distinct(product_table$LARGE_CATEGORY))) %>%
  mutate(beer = (j=="BEER/ALE/ALCOHOLIC CID"))
ggplot(df,aes(x=theta,group=j)) + 
  geom_density(color="grey") + 
  geom_density(data=df[df$beer,],aes(x=theta),color="#F8766D",size=1.2) + 
  labs(x="", y="density\n") +
  theme_minimal() +
  theme(legend.position="none",
        plot.title = element_text(hjust=0.5))
# dev.copy2pdf(file="figures/thetas_beer.pdf",width=4,height=3)

df = data.frame(theta = apply(out_hierridge$thetadraw[burn:end,(npar[1]+1):cumsum(npar)[2]],2,mean),
                SD = round(apply(out_hierridge$thetadraw[burn:end,(npar[1]+1):cumsum(npar)[2]],2,sd),3),
                i = rep(unique(product_table$SMALL_CATEGORY),n_distinct(product_table$SMALL_CATEGORY)),
                j = rep(unique(product_table$SMALL_CATEGORY),each=n_distinct(product_table$SMALL_CATEGORY))) %>%
  left_join(.,unique(product_table[,c("SMALL_CATEGORY","LARGE_CATEGORY")]),by=c("j"="SMALL_CATEGORY")) %>%
  mutate(imported = (j=="IMPORTED BEER/ALE (INC NON-ALCOH"),
         domestic = (j=="DOMESTIC BEER/ALE (INC NON-ALCOH"))
ggplot(df,aes(x=theta,group=j)) + 
  geom_density(color="grey") + 
  geom_density(data=df[df$imported,],aes(x=theta),color="#F8766D",size=1.2) + 
  geom_density(data=df[df$domestic,],aes(x=theta),color="#00BFC4",linetype=2, size=1.2) + 
  labs(x="", y="density\n") +
  theme_minimal() +
  theme(legend.position="none",
        plot.title = element_text(hjust=0.5))
# dev.copy2pdf(file="figures/thetas_beer_imported.pdf",width=4,height=3)
