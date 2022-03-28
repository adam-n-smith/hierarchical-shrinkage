library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(here)

source(here("src","shrinkage-functions.R"))
sourceCpp(here("src","shrinkage-mcmc.cpp"))

# load data
load(here("build","output","store_panel.RData"))

# build tree objects
childrencounts = countchildren_cpp(tree)
treeindex = createindex(tree)
index = treeindex$index
list = treeindex$list
list_own = treeindex$list_own
npar = treeindex$npar
npar_own = treeindex$npar_own
L = ncol(tree)
p = nrow(tree)

# --------------------------------------------------------- #
# split data into training/test
# --------------------------------------------------------- #

nweeks = nrow(demand)
inweeks = demand %>%
  select(IRI_KEY,WEEK) %>%
  rowid_to_column() %>%
  filter(WEEK<=max(WEEK)-26) %>%
  pull(rowid)

Y = demand %>%
  slice(inweeks) %>%
  select(starts_with("UNITS")) %>%
  as.matrix()
Ytest = demand %>%
  slice(-inweeks) %>%
  select(starts_with("UNITS")) %>%
  as.matrix()
X = prices %>%
  slice(inweeks) %>%
  select(starts_with("PRICE")) %>%
  as.matrix()
Xtest = prices %>%
  slice(-inweeks) %>%
  select(starts_with("PRICE")) %>%
  as.matrix()

n = nrow(Y)
ntest = nrow(Ytest)

Clist = NULL
Clisttest = NULL
nphi = 0
nphivec = double(p)
for(i in 1:p){
  
  feati = feature[inweeks,2+i]
  dispi = display[inweeks,2+i]
  promo = feati + dispi
  
  # training
  C = cbind(matrix(1,nrow=n),as.matrix(seasonal[inweeks,c("SUMMER","HOLIDAY")]))
  if(any(promo!=0) & n_distinct(promo)>1){
    C = cbind(C,as.matrix(promo))
  }
  colnames(C) = NULL
  nphi = nphi + ncol(C)
  nphivec[i] = ncol(C)
  Clist[[i]] = C
  
  # test
  C = matrix(1,nrow=ntest)
  C = cbind(matrix(1,nrow=ntest),as.matrix(seasonal[-inweeks,c("SUMMER","HOLIDAY")]))
  if(any(promo!=0) & n_distinct(promo)>1){
    C =  cbind(C,as.matrix(feature[-inweeks,2+i] + display[-inweeks,2+i]))
  }
  colnames(C) = NULL
  Clisttest[[i]] = C
  
}
cumnphi = c(0,cumsum(nphivec))