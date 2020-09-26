library(gtools)
library(prodlim)

createindex = function(tree){
  
  # number of levels
  L = ncol(tree)
  
  # index positions of each element from current parameter vector
  K = max(tree[,1])
  indexpairlast = permutations(K,2,1:K,2,repeats.allowed=TRUE)
  
  # fill in parameter index
  out = NULL
  out[[1]] = rep(1,nrow(indexpairlast))
  
  # repeat for the remaining levels
  if(L>1){
    for(i in 2:L){
      subtree = unique(tree[,1:i])
      K = max(subtree)
      
      # group-level parameters assumed to be assumtric, use "combinations" to assume symmetry
      if(i<L){
        
        # index positions of each element from current parameter vector
        indexpair = permutations(K,2,1:K,repeats.allowed=TRUE)
        
        # index positions of each element from last level's parameter vector
        indexpairlast = permutations(max(tree[,i-1]),2,1:max(tree[,i-1]),repeats.allowed=TRUE)
        
        # match current positions with last positions
        index = matrix(subtree[indexpair,i-1],ncol=2)
        match = apply(index,1,function(x)row.match(x,indexpairlast))
        out[[i]] = match
      }
      
      # SKU-level parameters
      else{
        
        # index positions of each element from current parameter vector
        indexpair = permutations(K,2,1:K,repeats.allowed=TRUE)
        
        # index positions of each element from last level's parameter vector
        indexpairlast = permutations(max(tree[,i-1]),2,1:max(tree[,i-1]),repeats.allowed=TRUE)
        
        # match current positions with last positions
        index = matrix(subtree[indexpair,i-1],ncol=2)
        match = apply(index,1,function(x)row.match(x,indexpairlast))
        
        # ignore own elasticities
        own = apply(indexpair,1,function(x)x[1]==x[2])
        out[[i]] = match[!own]
      }
    }
    
  }
  
  # save parameter index matrix
  l = out[[L]]
  index = l
  if(L>1){
    for(i in (L-1):1){
      l = out[[i]][l]
      index = rbind(l,index)
    }
  }
  rownames(index) = NULL
  
  # number of parameters per level
  npar = apply(index,1,max)
  if(L>1){
    npar = c(npar[-1],ncol(index))
  }
  
  return(list(list=out,index=index,npar=npar))
  
}
