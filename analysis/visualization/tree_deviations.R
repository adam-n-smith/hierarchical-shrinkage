df = data.frame(lambda = apply(out.hierhorse$lambdadraw[burn:end,(cumsum(npar)[2]+1):cumsum(npar)[3]],2,mean),
                i = rep(product_table$BRAND,245-1),
                j = rep(product_table$BRAND,each=245-1))
df = data.frame(lambda = apply(out.hierhorse$lambdadraw[burn:end,(cumsum(npar)[2]+1):cumsum(npar)[3]],2,mean),
                i = rep(1:245,245-1),
                j = rep(1:245,each=245-1))

df %>%
  ggplot(aes(x=i,y=j)) + 
  geom_tile(aes(fill=log(lambda))) + 
  scale_fill_gradient2(low = "white", mid="pink", high = "blue", midpoint=8)

