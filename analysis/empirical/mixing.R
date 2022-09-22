library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(coda)
library(here)

source(here("analysis","empirical","estimation-build.R"))
source(here("src","summary-functions.R"))

# load models
out_files = paste0("analysis/output/",dir(here("analysis","output")))
out_files = out_files[str_detect(out_files,"out-")]
lapply(out_files,load,.GlobalEnv)
models = mget(rev(ls()[str_detect(ls(),"out[.]")]))

# path
figurepath = "/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/"

# specs
end = nrow(models[[1]]$betadraws)
burn = 0.5*end
quantiles = c(0.1,0.5,0.9)
lags = seq(0,50,by=1)
colors = RColorBrewer::brewer.pal(8,"Blues")[c(4,6,8)]

# --------------------------------------------------------- #
# acf - cross elasticities
# --------------------------------------------------------- #

# index
partheta = cumsum(npar[-1])

# thetas (upper level)
acf.theta.top = getacf(models,lags,quantiles,burn,wchpar=(partheta[1]+1):partheta[2],own=FALSE)
formatout(acf.theta.top, "plot","grid") %>%
  ggplot(aes(x=lag,y=acf,group=percentile)) +
  geom_line(aes(color=percentile), show.legend=FALSE) +
  ylim(min(acf.theta.top$acf),1) +
  facet_grid(theta~beta, labeller=label_parsed) +
  scale_color_manual(values=colors) +
  theme_shrink()
ggsave(filename=paste0(figurepath,"mixing-acf-top.png"),
       height=3,width=5)

# thetas (middle level)
acf.theta.mid = getacf(models,lags,quantiles,burn,wchpar=1:partheta[1],own=FALSE)
formatout(acf.theta.mid, "plot","grid") %>%  
  ggplot(aes(x=lag,y=acf,group=percentile)) +
  geom_line(aes(color=percentile), show.legend=FALSE) +
  ylim(min(acf.theta.mid$acf),1) +
  facet_grid(theta~beta, labeller=label_parsed) +
  scale_color_manual(values=colors) +
  theme_shrink()
ggsave(filename=paste0(figurepath,"mixing-acf-mid.png"),
       height=3,width=5)

# betas
acf.beta = getacf(models,lags,quantiles,burn,wchpar="beta",own=FALSE)
formatout(acf.beta,"plot","grid") %>%
  ggplot(aes(x=lag,y=acf,group=percentile,color=percentile)) +
  geom_line(aes(color=percentile)) +
  ylim(min(acf.beta$acf),1) +
  scale_color_manual(values=colors) +
  facet_grid(theta~beta, labeller=label_parsed) +
  theme_shrink() +
  theme(legend.position = "bottom")
ggsave(filename=paste0(figurepath,"mixing-acf-beta.png"),
       height=4.5,width=5)

# --------------------------------------------------------- #
# acf - own elasticities
# --------------------------------------------------------- #

# index
partheta_own = cumsum(npar_own[-1])

# thetas (upper level)
acf.theta.top.own = getacf(models,lags,quantiles,burn,wchpar=(partheta_own[1]+1):partheta_own[2],own=TRUE)
formatout(acf.theta.top.own, "plot", "grid") %>%
  ggplot(aes(x=lag,y=acf,group=percentile)) +
  geom_line(aes(color=percentile), show.legend=FALSE) +
  ylim(min(acf.theta.top$acf),1) +
  facet_grid(theta~beta, labeller=label_parsed) +
  scale_color_manual(values=colors) +
  theme_shrink()
ggsave(filename=paste0(figurepath,"mixing-acf-top-own.png"),
       height=3,width=5)

# thetas (middle level)
acf.theta.mid.own = getacf(models,lags,quantiles,burn,wchpar=1:partheta_own[1],own=TRUE)
formatout(acf.theta.mid.own, "plot", "grid") %>%  
  ggplot(aes(x=lag,y=acf,group=percentile)) +
  geom_line(aes(color=percentile), show.legend=FALSE) +
  ylim(min(acf.theta.mid$acf),1) +
  facet_grid(theta~beta, labeller=label_parsed) +
  scale_color_manual(values=colors) +
  theme_shrink()
ggsave(filename=paste0(figurepath,"mixing-acf-mid-own.png"),
       height=3,width=5)

# betas
acf.beta.own = getacf(models,lags,quantiles,burn,wchpar="beta",own=TRUE)
formatout(acf.beta.own, "plot", "grid") %>%
  ggplot(aes(x=lag,y=acf,group=percentile,color=percentile)) +
  geom_line(aes(color=percentile)) +
  ylim(min(acf.beta.own$acf),1) +
  scale_color_manual(values=colors) +
  facet_grid(theta~beta, labeller=label_parsed) +
  theme_shrink() +
  theme(legend.position = "bottom")
ggsave(filename=paste0(figurepath,"mixing-acf-beta-own.png"),
       height=4.5,width=5)
