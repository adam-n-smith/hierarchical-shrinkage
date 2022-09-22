library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(here)

source(here("src","shrinkage-functions.R"))
sourceCpp(here("src","shrinkage-mcmc.cpp"))

# load data
load(here("build","output","store_panel.RData"))

# build tree objects
if(!("objects" %in% ls())){
  objects = create_treeobjects(tree)
  parindextree = objects$parindextree
  parindextree_own = objects$parindextree_own
  npar = objects$npar
  npar_own = objects$npar_own
  L = ncol(tree)
  p = nrow(tree)
}

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

  # training
  C = as.matrix(cbind(rep(1,n),seasonal[inweeks,c("SUMMER","HOLIDAY")],dispi))
  colnames(C) = c("intercept","summer","holiday","display")
  if(any(feati!=0) & n_distinct(feati)>1){
    C = cbind(C,as.matrix(feati))
    colnames(C) = c("intercept","summer","holiday","display","feature")
  }
  nphi = nphi + ncol(C)
  nphivec[i] = ncol(C)
  Clist[[i]] = C

  # test
  C = as.matrix(cbind(rep(1,ntest),seasonal[-inweeks,c("SUMMER","HOLIDAY")],display[-inweeks,2+i]))
  colnames(C) = c("intercept","summer","holiday","display")
  if(any(feati!=0) & n_distinct(feati)>1){
    C = cbind(C,as.matrix(feature[-inweeks,2+i]))
    colnames(C) = c("intercept","summer","holiday","display","feature")
  }
  Clisttest[[i]] = C

}
cumnphi = c(0,cumsum(nphivec))
