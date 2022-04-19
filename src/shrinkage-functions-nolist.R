library(gtools)
library(prodlim)


create_treeobjects = function(tree){
  
  # number of levels
  L = ncol(tree)
  p = nrow(tree)
  
  # number of parameters per level
  npar = apply(tree,2,max)^2
  npar[1] = npar[1] - p
  npar = c(npar,1)
  
  # number of parameters per level
  npar_own = apply(tree,2,max)
  npar_own = c(npar_own,1)
  
  # reverse
  tree = tree[,L:1]
  
  # index positions of each element from current parameter vector
  K = max(tree[,1])
  indexpairlast = permutations(K,2,1:K,2,repeats.allowed=TRUE)
  
  # fill in parameter index
  out = NULL
  out[[1]] = rep(1,nrow(indexpairlast))
  
  # repeat for the remaining levels
  if(L>1){
    
    # group-level parameters assumed to be asymmetric, use "combinations" to assume symmetry
    for(i in 2:L){
      
      subtree = unique(tree[,1:i])
      K = max(subtree)
      
      # index positions of each element from current parameter vector
      indexpair = permutations(K,2,1:K,repeats.allowed=TRUE)
      
      # index positions of each element from last level's parameter vector
      indexpairlast = permutations(max(tree[,i-1]),2,1:max(tree[,i-1]),repeats.allowed=TRUE)
      
      # match current positions with last positions
      index = matrix(subtree[indexpair,i-1],ncol=2)
      match = apply(index,1,function(x)row.match(x,data.frame(indexpairlast)))
      
      if(i<L){
        out[[i]] = match
      }
      else{
        own = apply(indexpair,1,function(x)x[1]==x[2])
        out[[i]] = match[!own]
      }
      
    }
    
  }
  
  # list for own parameters
  out_own = list(L)
  out_own[[1]] = rep(1,max(tree[,1]))
  for(l in 1:(L-1)){
    out_own[[l+1]] = unique(tree[,1:(l+1)])[,l]
  }
  
  return(list(npar=npar,npar_own=npar_own,parindextree=rev(out),parindextree_own=rev(out_own)))
  
}
