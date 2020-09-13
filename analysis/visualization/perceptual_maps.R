theta = apply(out.hierridge$thetadraw[burn:end,1:npar[1]],2,mean)
theta = matrix(theta,10,10)
theta = apply(out.hierridge$thetadraw[burn:end,(npar[1]+1):cumsum(npar)[2]],2,mean)
theta = matrix(theta,36,36)
pca = prcomp(theta, scale. = F, center = F)
plot(pca$x[, 1], pca$x[, 2], col=0)
text(pca$x[, 1], pca$x[, 2], labels = unique(product_table$SMALL_CATEGORY))

