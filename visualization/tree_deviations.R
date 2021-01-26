df = data.frame(lambda = apply(out_hierhorse$lambdadraw[burn:end,(cumsum(npar)[2]+1):cumsum(npar)[3]],2,mean),
                i = rep(product_table$BRAND,nrow(product_table)-1),
                j = rep(product_table$BRAND,each=nrow(product_table)-1))
df = data.frame(lambda = apply(out_hierhorse$lambdadraw[burn:end,(cumsum(npar)[2]+1):cumsum(npar)[3]],2,mean),
                i = rep(1:nrow(product_table),nrow(product_table)-1),
                j = rep(1:nrow(product_table),each=nrow(product_table)-1))

df %>%
  ggplot(aes(x=i,y=j)) + 
  geom_tile(aes(fill=1/(1+lambda))) + 
  scale_fill_gradient2(low = "white", mid="pink", high = "blue", midpoint=8)

