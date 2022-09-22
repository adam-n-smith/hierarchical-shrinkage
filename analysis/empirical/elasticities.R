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
# compute elasticities
# --------------------------------------------------------- #

# compute
end = nrow(models[[1]]$betadraws)
burn = 0.5*end
elasticities = summarize_elasticities(models,burn:end,p)
promotions = summarize_promos(models,burn:end,p)

# summary stats (own)
elasticities %>%
  filter(own) %>%
  mutate(sig = !(lower<0 & upper>0)) %>%
  group_by(model) %>%
  summarise(mean_all = mean(mean), 
            mean_sig = mean(mean[sig]),
            n_sig = mean(sig),
            share_neg = mean(mean<0))
elasticities %>%
  filter(own) %>%
  mutate(sig = !(lower<0 & upper>0)) %>%
  group_by(theta=="sparse") %>%
  summarise(mean_all = mean(mean), 
            median_all = median(mean), 
            n_sig = mean(sig),
            share_neg = mean(mean<0))

# summary stats (cross)
elasticities %>%
  filter(!own) %>%
  mutate(sig = !(lower<0 & upper>0)) %>%
  group_by(model) %>%
  summarise(mean_all = mean(mean), 
            n_sig = sum(sig))

# --------------------------------------------------------- #
# tables and figures for paper
# --------------------------------------------------------- #

# elasticity table
tbl = formatout(elasticities,"table") %>%
  mutate(own=ifelse(own==TRUE,"own","cross")) %>%
  group_by(theta,beta,own) %>%
  summarise(neg = 100*mean(mean<0),
            nsig = 100*mean(!(lower<0 & upper>0)),
            low = quantile(mean,0.1),
            mid = quantile(mean,0.5),
            avg = mean(mean),
            high = quantile(mean,0.9),.groups="drop") %>%
  pivot_wider(names_from=own,values_from=c(neg,nsig,avg,low,mid,high)) %>%
  mutate(beta = as.character(beta),
         beta = ifelse(beta=="$\\quad\\beta$-Ridge",paste0(theta,beta),beta)) %>%
  select(beta,contains("own"), avg_cross:high_cross)
print(xtable(tbl, digits=c(0,0,1,1,2,2,2,2,3,3,3,3)), 
      include.rownames=FALSE, 
      sanitize.text.function=function(x){x})

# effects of hierarchical shrinkage
qtlow = quantile(apply(prices[inweeks,-(1:2)],2,sd),0.25)
qthigh = quantile(apply(prices[inweeks,-(1:2)],2,sd),0.75)
sds = data.frame(id=1:p^2,
                 Var1=rep(1:p,each=p),
                 Var2=rep(1:p,times=p),
                 sd=rep(apply(prices[inweeks,-(1:2)],2,sd),times=p))
formatout(elasticities,"plot","wrap") %>%
  left_join(sds) %>%
  mutate(sd = ifelse(sd<qtlow,"Q1",ifelse(sd>qthigh,"Q3","(Q1,Q3)")),
         sd = relevel(factor(sd),"Q1")) %>% 
  filter(own) %>%
  select(id,beta,theta,mean,sd) %>%
  pivot_wider(names_from=theta,values_from=mean) %>%
  pivot_longer(-c(id,beta,sd,`Standard~Shrinkage`),names_to="theta",values_to="mean") %>%
  arrange(desc(sd)) %>%
  mutate(theta = factor(theta,levels=unique(theta))) %>%
  ggplot(aes(x=`Standard~Shrinkage`,y=mean)) +
  geom_point(aes(color=sd,shape=sd)) +
  geom_abline(slope=1,intercept=0) +
  scale_color_manual(values=c("#F8766D","grey75","#00BFC4")) +
  scale_shape_manual(values=c(17,16,15)) +
  xlim(-3,1) +
  ylim(-3,1) +
  labs(x="beta (standard shrinkage prior)",y="beta (hierarchical shrinkage prior)",shape="SD(price)",color="SD(price)") +
  facet_grid(theta~beta, labeller=label_parsed) +
  theme_shrink() +
  theme(legend.position = "bottom")
ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/elasticities-beta-scatter.png",
       height=4.5,width=5.5)

# own
formatout(elasticities,"plot","grid") %>%
  filter(own) %>%
  ggplot(aes(x=mean)) +
  geom_histogram(alpha=0.75) +
  xlim(-6,3) +
  facet_grid(theta~beta, labeller=label_parsed, scales="free") +
  labs(x="") +
  theme_shrink() +
  theme(axis.title.x = element_blank())
ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/elasticities-beta-own.png",
       height=4,width=5)

# cross
formatout(elasticities,"plot","grid") %>%
  filter(!own) %>%
  # mutate(significant = ifelse(lower<0 & upper>0,"no","yes")) %>%
  # ggplot(aes(x=mean, fill=significant)) + 
  ggplot(aes(x=mean)) +
  geom_histogram(alpha=0.75) +
  xlim(-0.25,0.25) +
  labs(x="") +
  facet_grid(theta~beta, labeller=label_parsed, scales="free") +
  theme_shrink() +
  theme(axis.title.x = element_blank())
ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/elasticities-beta-cross.png",
       height=4,width=5)

# own promos (display)
formatout(promotions,"plot","grid") %>%
  filter(variable=="display") %>%
  ggplot(aes(x=value)) +
  geom_histogram(alpha=0.75) +
  xlim(-2,5) +
  labs(x="") +
  facet_grid(theta~beta, labeller = label_parsed) +
  theme_shrink() +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        legend.title = element_blank())
ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/promotions-display.png",
       height=4.5,width=5)

# own promos (feature)
formatout(promotions,"plot","grid") %>%
  filter(variable=="feature") %>%
  ggplot(aes(x=value)) +
  geom_histogram(alpha=0.75) +
  xlim(-2,5) +
  labs(x="") +
  facet_grid(theta~beta, labeller = label_parsed) +
  theme_shrink() +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        legend.title = element_blank())
ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/promotions-feature.png",
       height=4.5,width=5)



