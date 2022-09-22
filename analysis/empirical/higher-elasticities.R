library(ggpubr)
library(snakecase)

source(here("src","summary-functions.R"))

end = nrow(models[[1]]$betadraws)
burn = 0.5*end

# model
out = out.ridge.ridge

# indices
partheta = cumsum(npar[-1])
partheta_own = cumsum(npar_own[-1])
wchtop = (partheta[1]+1):partheta[2]
wchmid = 1:partheta[1]
wchtop_own = (partheta_own[1]+1):partheta_own[2]
wchmid_own = 1:partheta_own[1]

# --------------------------------------------------------- #
# category-level elasticity matrix
# --------------------------------------------------------- #

df_top %>%
  ggplot(aes(i,forcats::fct_rev(j),fill=theta)) +
  geom_tile(color="white", lwd=1) +
  scale_fill_gradient2(low = "green", mid="white", high = "red") +
  labs(x="Change in Price of...", y="Change in Demand of...") +
  theme_shrink() +
  theme(axis.text.x = element_text(angle=90,hjust=1))

# --------------------------------------------------------- #
# theta own and theta across vs. within
# --------------------------------------------------------- #

df_top_own = data.frame(theta = round(apply(out$thetaowndraws[burn:end,wchtop_own],2,mean),3),
                        i = unique(product_table$LARGE_CATEGORY),
                        j = unique(product_table$LARGE_CATEGORY)) %>%
  mutate(within = "own",
         level = "category")
df_mid_own = data.frame(theta = round(apply(out$thetaowndraws[burn:end,wchmid_own],2,mean),3),
                        i = unique(product_table$SMALL_CATEGORY),
                        j = unique(product_table$SMALL_CATEGORY)) %>%
  mutate(within = "own",
         level = "subcategory")
df_top = data.frame(theta = round(apply(out$thetadraw[burn:end,wchtop],2,mean),3),
                    i = rep(unique(product_table$LARGE_CATEGORY),n_distinct(product_table$LARGE_CATEGORY)),
                    j = rep(unique(product_table$LARGE_CATEGORY),each=n_distinct(product_table$LARGE_CATEGORY))) %>%
  mutate(within = ifelse(i==j,"within","across"),
         level = "category")
df_mid = data.frame(theta = round(apply(out$thetadraw[burn:end,wchmid],2,mean),3),
                    i = rep(unique(product_table$SMALL_CATEGORY),n_distinct(product_table$SMALL_CATEGORY)),
                    j = rep(unique(product_table$SMALL_CATEGORY),each=n_distinct(product_table$SMALL_CATEGORY))) %>%
  mutate(within = ifelse(i==j,"within","across"),
         level = "subcategory")
df_own = rbind(df_top_own,df_mid_own)
plt_own = ggplot(df_own,aes(x=theta)) +
  geom_density(aes(linetype=level)) +
  labs(x="theta",title="Own-Price Elasticities") +
  theme_shrink() +
  xlim(-4,1) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust=0.5,size=11))
df = rbind(df_top,df_mid)
plt_cross = ggplot(df,aes(x=theta)) +
  geom_density(aes(linetype=level,color=within)) +
  labs(x="theta",title="Cross-Price Elasticities") +
  xlim(-0.06,0.06) +
  theme_shrink() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust=0.5,size=11)) +
  guides(color = guide_legend(order=2), 
         linetype = guide_legend(order=1))
ggarrange(plt_own,plt_cross,
          legend.grob = get_legend(plt_cross), legend="bottom")
ggsave(filename="/Users/adamsmith/Dropbox/Research/Hierarchical Shrinkage/paper/figures/elasticities-theta.png",
       height=3,width=6)

# --------------------------------------------------------- #
# example table
# --------------------------------------------------------- #

matrix(apply(out$thetadraw[burn:end,(partheta[1]+1):partheta[2]],2,mean),9,9)
round(matrix(apply(out$thetadraw[burn:end,1:partheta[1]],2,mean),28,28),3)

# top
dftop = data.frame(mean = apply(out$thetadraw[burn:end,wchtop],2,mean),
                   sd = apply(out$thetadraw[burn:end,wchtop],2,sd),
                   lower = apply(out$thetadraw[burn:end,wchtop],2,function(x)quantile(x,0.025)),
                   upper = apply(out$thetadraw[burn:end,wchtop],2,function(x)quantile(x,0.975)),
                   price_of = rep(unique(product_table$LARGE_CATEGORY),times=n_distinct(product_table$LARGE_CATEGORY)),
                   demand_of = rep(unique(product_table$LARGE_CATEGORY),each=n_distinct(product_table$LARGE_CATEGORY)))
dftop %>%
  group_by(demand_of) %>%
  summarise(max=price_of[which.max(mean)],
            maxval=max(mean),
            min=price_of[which.min(mean)],
            minval=min(mean)) %>%
  xtable(digits=3)      

# middle
dfmid = data.frame(mean = apply(out$thetadraw[burn:end,wchmid],2,mean),
                   sd = apply(out$thetadraw[burn:end,wchmid],2,sd),
                   lower = apply(out$thetadraw[burn:end,wchmid],2,function(x)quantile(x,0.025)),
                   upper = apply(out$thetadraw[burn:end,wchmid],2,function(x)quantile(x,0.975)),
                   price_of = rep(unique(product_table$SMALL_CATEGORY),times=n_distinct(product_table$SMALL_CATEGORY)),
                   demand_of = rep(unique(product_table$SMALL_CATEGORY),each=n_distinct(product_table$SMALL_CATEGORY)),
                   demand_of_cat = rep(unique(product_table[,2:3])$LARGE_CATEGORY,each=n_distinct(product_table[,2:3]))) %>%
  mutate(price_of = to_any_case(price_of,case="title",parsing_option = 0),
         price_of = str_remove(price_of,"\\([^(]*$"),
         demand_of = to_any_case(demand_of,case="title",parsing_option = 0),
         demand_of = str_remove(demand_of,"\\([^(]*$"),)

dfmid %>%
  group_by(demand_of,demand_of_cat) %>%
  summarise(max=price_of[which.max(mean)],
            maxval=max(mean),
            min=price_of[which.min(mean)],
            minval=min(mean)) %>%
  arrange(demand_of_cat) %>%
  select(-demand_of_cat) %>%
  xtable(digits=3)
