// [[Rcpp::depends(RcppArmadillo)]]
// #define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// ---------------------------------------------------------------------- //
// data manipulation functions
// ---------------------------------------------------------------------- //

// unlist function
vec unlist(List const& list, int const& total_length, int const& list_end){
  
  // loop and fill
  NumericVector output = no_init(total_length);
  int index=0;
  for(int i=0;i<list_end;i++){
    NumericVector el = list[i];
    std::copy(el.begin(), el.end(), output.begin() + index);
    index += el.size();
  }
  
  return as<arma::vec>(output);
  
}

// rep_each function
vec rep_each(vec const& x, double const& each){
  
  double n = x.n_elem;
  
  // NumericVector output = no_init(n * each);
  std::vector<double> output(n*each);
  
  int index=0;
  for(int i=0;i<n;i++){
    std::fill(output.begin()+index, output.begin()+index+each, x(i));
    index += each;
  }
  
  
  return as<arma::vec>(wrap(output));
  
}

// replace function
//[[Rcpp::export]]
vec replace_cpp(vec const& x, vec const& y){
  
  int n = y.n_elem;
  vec output = x;
  for(int i=0;i<n;i++){
    output(find(x==i+1)).fill(y(i));
  }
  return output;
}

mat unique_rows(mat const& m) {
  
  uvec ulmt = zeros<arma::uvec>(m.n_rows);
  for(uword i = 0; i < m.n_rows; i++) {
    for(uword j = i + 1; j < m.n_rows; j++) {
      if (all(m.row(i)==m.row(j))) { ulmt(j) = 1; break; }
    }
  }
  
  return m.rows(find(ulmt == 0));
  
}

//[[Rcpp::export]]
List countchildren(mat const& tree, bool const& own){
  
  vec countsi;
  int L = tree.n_cols;
  List counts(L-1);
  for(int i=1; i<L; i++){
    
    int k = max(tree.col(i));
    vec tbl = zeros(k);
    
    // counting children for own elasticities
    if(own){

      countsi = zeros(k);
      for(int level=0; level<i; level++){
        mat subtree = unique_rows(tree.cols(span(level,i)));
        for(int j=0;j<k;j++){
          int lastcol = subtree.n_cols-1;
          uvec mtch = find(subtree.col(lastcol)==j+1);
          tbl(j) = mtch.size();
        }
        countsi += tbl;
      }
      
    }
    // counting children for cross elasticities 
    // (subtract off count of own elasticities on diagonals)
    else{
      
      countsi = zeros(k*k);
      for(int level=0; level<i; level++){
        mat subtree = unique_rows(tree.cols(span(level,i)));
        for(int j=0;j<k;j++){
          int lastcol = subtree.n_cols-1;
          uvec mtch = find(subtree.col(lastcol)==j+1);
          tbl(j) = mtch.size();
        }
        if(level==0) countsi += vectorise(tbl*trans(tbl) - diagmat(tbl));
        else countsi += vectorise(tbl*trans(tbl));
      }
      
    }
    counts[i-1] = countsi;
  }
  
  return counts;
}

// compute sums of x for each group
//[[Rcpp::export]]
vec groupsum_cpp(vec const& x, vec const& groups, int const& n){
  vec output = zeros(n);
  for(int j=0;j<n;j++){
    output(j) = sum(x(find(groups==j+1)));
  }
  return output;
}

// compute means of x for each group
//[[Rcpp::export]]
vec groupmean_cpp(vec const& x, vec const& groups, int const& n){
  vec output = zeros(n);
  for(int j=0;j<n;j++){
    uvec wch = find(groups==j+1);
    if(wch.n_elem>0)
      output(j) = mean(x(wch));
  }
  return output;
}

// print time
void print_time(int const& rep, int const& Rep, double const& ctime){
  // time
  if(rep==0){
    if(ctime/60.0*(Rep-rep-1)/(rep+1)<1){
      Rprintf("Estimated time: %.1f seconds \n",ctime*(Rep-rep-1)/(rep+1));
      Rprintf("  0%%   10   20   30   40   50   60   70   80   90   100%%\n");
      Rprintf("  |----|----|----|----|----|----|----|----|----|----|\n");
    }
    else{
      Rprintf("Estimated time: %.1f minutes \n",ctime/60.0*(Rep-rep-1)/(rep+1));
      Rprintf("  0%%   10   20   30   40   50   60   70   80   90   100%%\n");
      Rprintf("  |----|----|----|----|----|----|----|----|----|----|\n");
    }
  }
  if ((rep+1)%(Rep/50)==0){
    if(rep+1==Rep/50){
      Rprintf("  *");
    }
    else if(rep<Rep-1){
      Rprintf("*");
    }
    else Rprintf("*\n");
  }
}

// ---------------------------------------------------------------------- //
// sampling functions
// ---------------------------------------------------------------------- //

// generate inverse-gamma random variables
//[[Rcpp::export]]
vec randig(int const& n, vec const& shape, vec const& scale){
  // note that R::rgamma has (shape, scale) parameterization
  // if X ~ Gamma(shape,rate) then 1/X ~ InvGamma(shape,scale)
  vec output(n);
  for(int i=0;i<n;i++){
    output(i) = 1/R::rgamma(as_scalar(shape(i)),1/as_scalar(scale(i)));
  }
  
  // bound small values to prevent numerical overflow issues
  output(find(output<1.0e-6)).fill(1.0e-6);
  output(find(output>1.0e6)).fill(1.0e6);

  return output;
}

// generate inverse gaussian random variables
//[[Rcpp::export]]
vec randinvgaussian(int const& n, vec const& mu, double const& lambda){
  vec nu = randn(n);
  vec y = pow(nu,2.0);
  vec musq = pow(mu,2.0);
  vec output = mu + musq%y/2.0/lambda - mu/2.0/lambda%sqrt(4.0*lambda*mu%y + musq%pow(y,2.0));
  output(find(output<1.0e-6)).fill(1.0e-6);
  uvec wch = find(randu(n) > mu/(mu+output));
  output(wch) = musq(wch)/output(wch);
  return output;
}

// ---------------------------------------------------------------------- //
// MCMC functions
// ---------------------------------------------------------------------- //

// //[[Rcpp::export]]
// vec drawbeta(mat const& Y, mat const& X, bool fast){
//   
//   int n = Y.n_rows;
//   int p = Y.n_cols;
//   vec lamstar = ones(p*p);
//   vec beta = zeros(p*p);
//   vec betabar = zeros(p*p);
//   
//   if(fast){
//     
//     for(int i=0;i<p;i++){
// 
//       // precompute
//       mat LamXp = diagmat(lamstar(span(p*i,p*i+p-1))) * trans(X);
//       mat irootXt = solve(trimatu(chol(X*LamXp + eye(n,n))),eye(n,n));
// 
//       // Ystar
//       vec ytstar = vectorise(Y.col(i));
// 
//       // beta (conditional)
//       vec u = betabar(span(p*i,p*i+p-1)) + sqrt(lamstar(span(p*i,p*i+p-1)))%randn<vec>(p);
//       vec v = X*u + randn<vec>(n);
//       vec w = (irootXt*trans(irootXt))*(ytstar - v);
//       beta(span(p*i,p*i+p-1)) = u + LamXp*w;
// 
//     }
//     
//   }
//   else{
//     
//     for(int i=0;i<p;i++){
// 
//       // precompute
//       mat Lam = diagmat(lamstar(span(p*i,p*i+p-1)));
//       mat iroot = solve(trimatu(chol(trans(X)*X + Lam)),eye(p,p));
//       
//       // Ystar
//       vec ytstar = vectorise(Y.col(i));
//       
//       // beta (conditional)
//       vec betatilde = (iroot*trans(iroot))*(trans(X)*ytstar + Lam*betabar(span(p*i,p*i+p-1)));
//       beta(span(p*i,p*i+p-1)) = betatilde + iroot*randn<vec>(p);
// 
//     }
//     
//   }
//   
//   return beta;
// }
// 

//[[Rcpp::export]]
vec draw_theta(int const& ell, vec const& npar, List const& parindextree,
               vec const& theta, vec const& Psi, vec const& tausq){

  vec Psi_above, Psi_below, prec_above, prec_below, prec, theta_above, theta_below, thetatilde, theta_ell;
  
  vec endindex = cumsum(npar)-1;
  vec begindex = endindex-npar+1;
  int L = tausq.size();
  
  // above = prior
  // below = data

  // precision
  prec_above = 1 / Psi(span(begindex(ell),endindex(ell))) / tausq(ell);
  prec_above = prec_above % replace_cpp(parindextree[ell], 1 / Psi(span(begindex(ell+1),endindex(ell+1))));
  prec_below = 1 / Psi(span(begindex(ell-1),endindex(ell-1))) / tausq(ell-1);
  prec_below = groupsum_cpp(prec_below,parindextree[ell-1],npar(ell));
  prec = prec_above + prec_below;
  
  // mean
  theta_above = theta(span(begindex(ell+1),endindex(ell+1)));
  theta_above = replace_cpp(parindextree[ell],theta_above)/Psi(span(begindex(ell),endindex(ell)))/tausq(ell);
  theta_below = theta(span(begindex(ell-1),endindex(ell-1)))/Psi(span(begindex(ell-1),endindex(ell-1)))/tausq(ell-1);
  theta_below = groupsum_cpp(theta_below,parindextree[ell-1],npar(ell));
  thetatilde = 1/prec % (theta_above + theta_below);
  
  // draw theta
  theta_ell = thetatilde + sqrt(1/prec) % randn(npar(ell));

  return theta_ell;

}

//[[Rcpp::export]]
vec draw_lambdasq(int const& ell, vec const& npar, List const& parindextree, List const& childrencounts,
                  vec const& theta, vec const& Psi, vec const& lambdasq, vec const& xilambda, vec const& tausq,
                  std::string product_shrinkage, std::string group_shrinkage){

  vec counts, levelsums, shape, scale, mean, sums, xi, mutilde;
  vec out = ones(npar(ell));
  vec endindex = cumsum(npar)-1;
  vec begindex = endindex-npar+1;
  int L = tausq.n_elem;

  // product-lasso
  if(ell==0 && product_shrinkage=="lasso"){

    vec Psistar = replace_cpp(parindextree[ell],Psi(span(begindex(ell+1),endindex(ell+1))));

    // prior mean
    mean = replace_cpp(parindextree[ell],theta(span(begindex(ell+1),endindex(ell+1))));

    // beta - prior mean
    xi = theta(span(begindex(ell),endindex(ell))) - mean;

    // draw
    mutilde = sqrt(2.0)*sqrt(Psistar)*sqrt(tausq(ell))/abs(xi);
    out = 1/randinvgaussian(npar(ell),mutilde,2);

  }
  // product-horseshoe
  else if(ell==0 && product_shrinkage=="horseshoe"){

    vec Psistar = replace_cpp(parindextree[ell],Psi(span(begindex(ell+1),endindex(ell+1))));

    // shape
    shape = 0.5 + 0.5*ones(npar(ell));

    // scale: sums at level ell
    mean = replace_cpp(parindextree[ell],theta(span(begindex(ell+1),endindex(ell+1))));
    xi = theta(span(begindex(ell),endindex(ell))) - mean;
    levelsums = pow(xi,2.0)/Psistar/tausq(ell);
    scale = 1.0/xilambda(span(begindex(ell),endindex(ell))) + 0.5*levelsums;

    // draw
    out = randig(npar(ell),shape,scale);

  }
  // group-horseshoe
  else if(ell>0 && ell<L-1 && group_shrinkage=="horseshoe"){
    
    vec Psistar = ones(npar(ell));
    // vec Psistar = replace_cpp(parindextree[ell],Psi(span(begindex(ell+1),endindex(ell+1))));
    
    // shape: count children nodes at each level
    counts = as<vec>(childrencounts[ell-1]);
    shape = 0.5 + 0.5*(1 + counts);

    // scale (i): sums at level ell
    mean = replace_cpp(parindextree[ell],theta(span(begindex(ell+1),endindex(ell+1))));
    xi = theta(span(begindex(ell),endindex(ell))) - mean;
    levelsums = pow(xi,2.0)/Psistar/tausq(ell);

    // scale (ii): sums at all levels below ell
    for(int s=(ell-1);s>=0;s--){

      Psistar = replace_cpp(parindextree[s],Psistar) % lambdasq(span(begindex(s),endindex(s)));
      mean = replace_cpp(parindextree[s], theta(span(begindex(s+1),endindex(s+1))));
      xi = theta(span(begindex(s),endindex(s))) - mean;
      levelsums += groupsum_cpp(pow(xi,2.0)/Psistar/tausq(s),parindextree[s],npar(ell));
     
      // if(ell-s==1) levelsums += groupsum_cpp(pow(xi,2.0)/Psistar/tausq(s),parindextree[s],npar(ell));
      // else levelsums += groupsum_cpp(pow(xi,2.0)/Psistar/tausq(s), replace_cpp(parindextree[s],parindextree[s+1]),npar(ell));

    }

    // scale (iii)
    scale = 1.0/xilambda(span(begindex(ell),endindex(ell))) + 0.5*levelsums;

    // draw
    out = randig(npar(ell),shape,scale);

  }

  return out;

}

// vec lndnorm(vec const& y, vec const& mu, vec const& sigma){
//   return -log(sigma) - 0.5*((y-mu)/sigma) % ((y-mu)/sigma);
// 
// }
// //[[Rcpp::export]]
// List MHdraw_lambdasq(int const& ell, vec const& npar, List const& parindextree,
//                     vec const& theta, vec const& Psi, vec const& lambdasq, vec const& tausq){
// 
//   uvec wch;
//   vec llold, llnew;
//   vec accept = zeros(npar(ell));
//   
//   vec endindex = cumsum(npar)-1;
//   vec begindex = endindex-npar+1;
// 
//   vec loglambdasq = log(lambdasq(span(begindex(ell),endindex(ell))));
//   vec loglambdasqnew = loglambdasq + 1*randn(npar(ell));
// 
//   vec thetaell = theta(span(0,endindex(ell)));
//   vec thetabarell = theta(span(0,endindex(ell)));
//   vec Psinew = Psi(span(0,endindex(ell)));
// 
//   thetabarell(span(begindex(ell),endindex(ell))) = replace_cpp(parindextree[ell], theta(span(begindex(ell+1),endindex(ell+1))));
//   Psinew(span(begindex(ell),endindex(ell))) = exp(loglambdasqnew) % replace_cpp(parindextree[ell],Psi(span(begindex(ell+1),endindex(ell+1))));
//   vec logold = lndnorm(thetaell(span(begindex(ell),endindex(ell))), thetabarell(span(begindex(ell),endindex(ell))), sqrt(Psi(span(begindex(ell),endindex(ell)))*tausq(ell)));
//   vec lognew = lndnorm(thetaell(span(begindex(ell),endindex(ell))), thetabarell(span(begindex(ell),endindex(ell))), sqrt(Psinew(span(begindex(ell),endindex(ell)))*tausq(ell)));
//   for(int s=(ell-1);s>=0;s--){
//     thetabarell(span(begindex(s),endindex(s))) = replace_cpp(parindextree[s], theta(span(begindex(s+1),endindex(s+1))));
//     Psinew(span(begindex(s),endindex(s))) = lambdasq(span(begindex(s),endindex(s))) % replace_cpp(parindextree[s], Psi(span(begindex(s+1),endindex(s+1))));
//     llold = lndnorm(thetaell(span(begindex(s),endindex(s))), thetabarell(span(begindex(s),endindex(s))), sqrt(Psi(span(begindex(s),endindex(s)))*tausq(s)));
//     llnew = lndnorm(thetaell(span(begindex(s),endindex(s))), thetabarell(span(begindex(s),endindex(s))), sqrt(Psinew(span(begindex(s),endindex(s)))*tausq(s)));
//     logold += groupsum_cpp(llold,parindextree[s],npar(ell));
//     lognew += groupsum_cpp(llnew,parindextree[s],npar(ell));
//   }
// 
//   // lambda ~ C+(0,1) -> lambda^2 ~ inverted-beta
//   vec logoldprior = -0.5*loglambdasq - log(1+exp(loglambdasq));
//   vec lognewprior = -0.5*loglambdasqnew - log(1+exp(loglambdasqnew));
// 
//   vec ldiff = lognew + lognewprior - logold - logoldprior;
//   vec u = randu(npar(ell));
//   uvec wchacpt = find(u < exp(ldiff));
//   loglambdasq(wchacpt) = loglambdasqnew(wchacpt);
//   accept(wchacpt).fill(1);
//     
//   vec out = exp(loglambdasq);
//   // out(find(out<1e-8)).fill(1e-8);
//   // out(find(out>1e2)).fill(1e2);
//   
//   // return out;
//   return List::create(
//     Named("draw") = out,
//     Named("accept") = accept
//   );
// 
// }

//[[Rcpp::export]]
List rSURshrinkage(List Data, List Prior, List Mcmc, std::string Shrinkage, bool print){
  
  // data
  mat Y = Data["Y"];
  mat X = Data["X"];
  double p = X.n_cols;
  double n = X.n_rows;
  int npar = p*p-p;
  List Clist = Data["Clist"];
  
  // prior
  mat Aphi = Prior["Aphi"];
  vec phibar = Prior["phibar"];
  double a = Prior["a"];
  double b = Prior["b"];
  
  // mcmc
  int Rep = Mcmc["R"];
  int keep = Mcmc["keep"];
  
  // initialize
  int rep, mkeep;
  vec ytstar, phitilde, levelsums, diffsums, u, v, vari, w, sums, rate, scale, mutilde;
  mat Ystar, Xt, CtpCt, LamXtp, Xpy, XtLamXtp, Lamibetabar, irootX, irootXt, irootC, Ct, projmat, CIprojCp;
  mat XpX = trans(X)*X;
  vec nphi = zeros(p);
  mat Cstar;
  for(int i=0;i<p;i++){
    mat C = Clist[i];
    nphi(i) = C.n_cols;
    Cstar = join_rows(Cstar, C);
  }
  vec cumnphi = cumsum(nphi);
  mat CspCs = trans(Cstar)*Cstar;
  
  // initial values
  vec phi = zeros(sum(nphi));
  uvec wchown = find(eye<mat>(p,p)==1);
  uvec wchcross = find(eye<mat>(p,p)==0);
  vec beta = vectorise(inv(trans(X)*X+eye(p,p))*trans(X)*Y);
  vec betabar = zeros(p*p);
  vec lamstar = ones(p*p);
  vec tausq = ones(1);
  vec xitau = tausq;
  vec tausqown = ones(1);
  vec xitauown = tausqown;
  vec lambdasq = ones(npar);
  vec lambdasqown = ones(p);
  vec xilambda = ones(npar);
  vec xilambdaown = ones(p);
  vec sigmasq = ones(p);
  
  // storage matrices
  mat phidraws(Rep/keep,sum(nphi));
  mat betadraws(Rep/keep, p*p);
  mat lambdasqdraws(Rep/keep, npar);
  mat lambdasqowndraws(Rep/keep, p);
  mat tausqdraws(Rep/keep, 1);
  mat tausqowndraws(Rep/keep, 1);
  mat sigmasqdraws(Rep/keep, p);
  
  // print progress banner
  wall_clock timer;
  timer.tic();
  Datetime dt;
  if(print) Rprintf("MCMC Progress \n");
  
  // MCMC loop
  for (rep=0; rep<Rep; rep++){
    
    // -------------------------------------------------------------- //
    // phi, beta, sigmasq
    // -------------------------------------------------------------- //
    
    lamstar(wchcross) = lambdasq*as_scalar(tausq);
    lamstar(wchown) = lambdasqown*as_scalar(tausqown);
    
    for(int i=0;i<p;i++){
      
      // precompute
      Xt = X/sqrt(sigmasq(i));
      LamXtp = diagmat(lamstar(span(p*i,p*i+p-1))) *  trans(Xt);
      XtLamXtp = Xt * LamXtp;
      
      // phi (marginalizing over beta)
      // first use Woodbury to rewrite inverse as (Xt'Xt + Lam^-1)^-1 = Lam - LamXt'(I + XtLamXt')^-1XtLam
      // then projection matrix is can be written as:
      // I - Xt(Xt'Xt + Lam^-1)^-1Xt'
      // = I - Xt(Lam - LamXt'(I + XtLamXt')^-1XtLam)Xt'
      // = I - (XtLamXt' - XtLamXt'(I + XtLamXt')^-1XtLamXt'
      mat C = Clist[i];
      Ct = C/sqrt(sigmasq(i));
      irootXt = solve(trimatu(chol(symmatu(XtLamXtp + eye(n,n)))),eye(n,n));
      projmat = eye(n,n) - (XtLamXtp - XtLamXtp*irootXt*trans(irootXt)*XtLamXtp);
      irootC = solve(trimatu(chol(symmatu(trans(Ct)*projmat*Ct + Aphi(i,i)))), eye(nphi(i),nphi(i)));
      phitilde = (irootC*trans(irootC))*(trans(Ct)*projmat*Y.col(i)/sqrt(sigmasq(i)) + Aphi(i,i)*phibar(i));
      phi(span(cumnphi(i)-nphi(i),cumnphi(i)-1)) = phitilde + irootC*randn(nphi(i));
      
      // Ystar
      Ystar = Y.col(i) - C*phi(span(cumnphi(i)-nphi(i),cumnphi(i)-1));
      ytstar = vectorise(Ystar)/sqrt(sigmasq(i));
      
      // beta (conditional)
      u = betabar(span(p*i,p*i+p-1)) + sqrt(lamstar(span(p*i,p*i+p-1))) % randn<vec>(p);
      v = Xt * u + randn<vec>(n);
      w = (irootXt*trans(irootXt)) * (ytstar - v);
      beta(span(p*i,p*i+p-1)) = u + LamXtp * w;
      
      // sigmasq
      vec E = Y.col(i) - X*beta(span(p*i,p*i+p-1)) - C*phi(span(cumnphi(i)-nphi(i),cumnphi(i)-1));
      scale = a + 0.5*trans(E)*E;
      sigmasq(i) = 1/R::rgamma(b + 0.5*n,1/as_scalar(scale));
      
    }
    
    // tausq -------------------------------------------------------- //
    
    // cross
    levelsums = sum(pow(beta(wchcross),2.0)/lambdasq);
    tausq = 1/R::rgamma(0.5*(p*p-p+1), 1/as_scalar(1/xitau+0.5*levelsums));
    xitau = 1/R::rgamma(1,as_scalar(1/(1+1/tausq)));
    
    // own
    levelsums = sum(pow(beta(wchown),2.0));
    tausqown = 1/R::rgamma(0.5*(p+1), 1/as_scalar(1/xitauown+0.5*levelsums));
    xitauown = 1/R::rgamma(1,as_scalar(1/(1+1/tausqown)));
    
    // lambdasq  ---------------------------------------------------- //
    
    // lasso: Exp(1/2)=Gamma(1,1/2)
    if(Shrinkage=="lasso"){
      
      // cross
      mutilde = pow(2*tausq(0)/pow(beta(wchcross),2.0),0.5);
      lambdasq = 1/randinvgaussian(npar,mutilde,2);
      
      // own
      mutilde = pow(2*tausqown(0)/pow(beta(wchown),2.0),0.5);
      lambdasqown = 1/randinvgaussian(p,mutilde,2);
      
    }
    // horseshoe: C+(0,1)
    if(Shrinkage=="horseshoe"){
      
      // cross
      scale = 1/xilambda + 0.5*pow(beta(wchcross),2.0)/as_scalar(tausq);
      lambdasq = randig(npar,ones(npar),scale);
      xilambda = randig(npar,ones(npar),1+1/lambdasq);
      
      // own
      scale = 1/xilambdaown + 0.5*pow(beta(wchown),2.0)/as_scalar(tausqown);
      lambdasqown = randig(p,ones(p),scale);
      xilambdaown = randig(p,ones(p),1+1/lambdasqown);
      
    }
    
    // -------------------------------------------------------------- //
    // print time and store draws
    // -------------------------------------------------------------- //
    
    // time
    if(print) print_time(rep, Rep, timer.toc());
    
    // store draws
    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep - 1;
      betadraws(mkeep,span::all) = trans(beta);
      lambdasqdraws(mkeep,span::all) = trans(lambdasq);
      lambdasqowndraws(mkeep,span::all) = trans(lambdasqown);
      tausqdraws(mkeep,span::all) = tausq;
      tausqowndraws(mkeep,span::all) = tausqown;
      sigmasqdraws(mkeep,span::all) = trans(sigmasq);
      phidraws(mkeep,span::all) = trans(phi);
    }
    
  }
  
  // print total time elapsed
  if(print){
    if(timer.toc()/60.0<1) Rprintf("Total Time Elapsed: %.1f seconds \n",timer.toc());
    else Rprintf("Total Time Elapsed: %.1f minutes \n",timer.toc()/60.0);
  }
  
  return List::create(
    Named("betadraws") = betadraws,
    Named("lambdasqdraws") = lambdasqdraws,
    Named("lambdasqowndraws") = lambdasqowndraws,
    Named("tausqdraws") = tausqdraws,
    Named("tausqowndraws") = tausqowndraws,
    Named("sigmasqdraws") = sigmasqdraws,
    Named("phidraws") = phidraws
  );
  
}

//[[Rcpp::export]]
List rSURhiershrinkage(List const& Data, List const& Prior, List const& Mcmc, List const& Shrinkage, bool print){

  // data
  mat Y = Data["Y"];
  mat X = Data["X"];
  double p = X.n_cols;
  double n = X.n_rows;
  List Clist = Data["Clist"];
  vec npar = Data["npar"];
  vec npar_own = Data["npar_own"];
  List parindextree = Data["parindextree"];
  List parindextree_own = Data["parindextree_own"];
  mat tree = Data["tree"];
  int L = tree.n_cols;

  List childrencounts = countchildren(tree,false);
  List childrencounts_own = countchildren(tree,true);
  
  // prior
  double thetabar_cross = Prior["thetabar_cross"];
  double thetabar_own = Prior["thetabar_own"];
  mat Aphi = Prior["Aphi"];
  vec phibar = Prior["phibar"];
  double a = Prior["a"];
  double b = Prior["b"];

  // mcmc
  int Rep = Mcmc["R"];
  int keep = Mcmc["keep"];
  vec accept = zeros(sum(npar));
  
  // shrinkage
  std::string product_shrinkage = Shrinkage["product"];
  std::string group_shrinkage = Shrinkage["group"];

  // initialize
  int rep, mkeep;
  vec phitilde, ytstar, u, v, w, betadraw, rate, scale, mean, levelsums, res, Psistar, Psistar_own;
  mat Ystar, Xt, CtpCt, LamXtp, Xpy, XtLamXtp, Lamibetabar, irootX, irootXt, irootC, Ct, projmat, CIprojCp;

  // indices
  vec endindex = cumsum(npar)-1;
  vec begindex = endindex-npar+1;
  vec endindex_own = cumsum(npar_own)-1;
  vec begindex_own = endindex_own-npar_own+1;
  uvec wchown = find(eye<mat>(p,p)==1);
  uvec wchcross = find(eye<mat>(p,p)==0);
  
  // initial values: product elasticities
  vec beta = vectorise(inv(trans(X)*X+0.1*eye(p,p))*trans(X)*Y);
  
  // initial values: cross elasticities
  vec theta = zeros(sum(npar));
  theta(sum(npar)-1) = thetabar_cross;
  theta(span(0,npar(0)-1)) = beta(wchcross);
  vec lambdasq = ones(sum(npar));
  vec xilambda = ones(sum(npar));
  vec Psi = ones(sum(npar));
  vec tausq = ones(L);
  vec xitau = ones(L);
  
  // initial values: own elasticities
  vec theta_own = zeros(sum(npar_own));
  theta(sum(npar_own)-1) = thetabar_own;
  theta_own(span(0,npar_own(0)-1)) = beta(wchown);
  vec lambdasq_own = ones(sum(npar_own));
  vec xilambda_own = ones(sum(npar_own));
  vec Psi_own = ones(sum(npar_own));
  vec tausq_own = ones(L);
  vec xitau_own = ones(L);
  
  // initial values: everything else
  vec xi = zeros(p*p);
  vec betabar = zeros(p*p);
  vec Lambdastar = ones(p*p);
  vec sigmasq = ones(p);
  mat XpX = trans(X)*X;
  vec nphi = zeros(p);
  mat Cstar;
  for(int i=0;i<p;i++){
    mat C = Clist[i];
    nphi(i) = C.n_cols;
    Cstar = join_rows(Cstar, C);
  }
  vec cumnphi = cumsum(nphi);
  mat CspCs = trans(Cstar)*Cstar;
  vec phi = zeros(sum(nphi));
  
  // storage matrices
  // mat betadraws(Rep/keep, p*p);
  // mat thetadraws(Rep/keep, sum(npar(span(1,L-1))));
  // mat thetaowndraws(Rep/keep, sum(npar_own(span(1,L-1))));
  // mat lambdasqdraws(Rep/keep,sum(npar));
  // mat lambdasqowndraws(Rep/keep,sum(npar_own));
  // mat Psidraws(Rep/keep,sum(npar));
  // mat Psiowndraws(Rep/keep,sum(npar_own));
  // mat tausqdraws(Rep/keep, L);
  // mat tausqowndraws(Rep/keep, L);
  // mat sigmasqdraws(Rep/keep, p);
  // mat phidraws(Rep/keep,sum(nphi));
  mat betadraws(Rep/keep, p*p);
  mat lambdasqdraws(Rep/keep,p*p);
  mat thetadraws(Rep/keep, sum(npar(span(1,L-1))));
  mat thetaowndraws(Rep/keep, sum(npar_own(span(1,L-1))));
  mat psisqdraws(Rep/keep, sum(npar(span(1,L-1))));
  mat psisqowndraws(Rep/keep, sum(npar_own(span(1,L-1))));
  mat Psidraws(Rep/keep,sum(npar)-1);
  mat Psiowndraws(Rep/keep,sum(npar_own)-1);
  mat tausqdraws(Rep/keep, L);
  mat tausqowndraws(Rep/keep, L);
  mat sigmasqdraws(Rep/keep, p);
  mat phidraws(Rep/keep,sum(nphi));
  
  // print progress banner
  wall_clock timer;
  timer.tic();
  if(print) Rprintf("MCMC Progress \n");

  // MCMC loop
  for (rep=0; rep<Rep; rep++){

    // -------------------------------------------------------------- //
    // top -> bottom
    // -------------------------------------------------------------- //
    
    // HIGHER-LEVEL PARAMETERS -------------------------------------- //
    
    for(int ell=L-1; ell>0;ell--){
      
      // CROSS ELASTICITIES ----------------------------------------- //

      // theta
      theta(span(begindex(ell),endindex(ell))) = draw_theta(ell, npar, parindextree, theta, Psi, tausq);

      // tausq
      res = theta(span(begindex(ell),endindex(ell))) - replace_cpp(parindextree[ell],theta(span(begindex(ell+1),endindex(ell+1))));
      levelsums = sum(pow(res,2.0)/Psi(span(begindex(ell),endindex(ell))));
      tausq(ell) = 1.0/R::rgamma(0.5 + 0.5*npar(ell), 1.0/(1.0/xitau(ell)+0.5*as_scalar(levelsums)));
      if(tausq(ell)<1.0e-6) tausq(ell) = 1.0e-6;
      xitau(ell) = 1.0/R::rgamma(1,1.0/(1+1.0/tausq(ell)));

      // lambdasq
      lambdasq(span(begindex(ell),endindex(ell))) = draw_lambdasq(ell,npar,parindextree,childrencounts,theta,Psi,lambdasq,xilambda,tausq,product_shrinkage,group_shrinkage);
      xilambda(span(begindex(ell),endindex(ell))) = randig(npar(ell),ones(npar(ell)),1 + 1.0/lambdasq(span(begindex(ell),endindex(ell))));
      Psi(span(begindex(ell),endindex(ell))) = lambdasq(span(begindex(ell),endindex(ell))) % replace_cpp(parindextree[ell],Psi(span(begindex(ell+1),endindex(ell+1))));
      
      // OWN ELASTICITIES ------------------------------------------- //
      
      // theta
      theta_own(span(begindex_own(ell),endindex_own(ell))) = draw_theta(ell, npar_own, parindextree_own, theta_own, Psi_own, tausq_own);
      
      // tausq
      res = theta_own(span(begindex_own(ell),endindex_own(ell))) - replace_cpp(parindextree_own[ell],theta_own(span(begindex_own(ell+1),endindex_own(ell+1))));
      levelsums = sum(pow(res,2.0)/Psi_own(span(begindex_own(ell),endindex_own(ell))));
      tausq_own(ell) = 1.0/R::rgamma(0.5 + 0.5*npar_own(ell), 1.0/(1.0/xitau_own(ell)+0.5*as_scalar(levelsums)));
      xitau_own(ell) = 1.0/R::rgamma(1,1.0/(1+1.0/tausq_own(ell)));
      
      // lambdasq
      lambdasq_own(span(begindex_own(ell),endindex_own(ell))) = draw_lambdasq(ell,npar_own,parindextree_own,childrencounts_own,theta_own,Psi_own,lambdasq_own,xilambda_own,tausq_own,product_shrinkage,group_shrinkage);
      xilambda_own(span(begindex_own(ell),endindex_own(ell))) = randig(npar_own(ell),ones(npar_own(ell)),1 + 1.0/lambdasq_own(span(begindex_own(ell),endindex_own(ell))));
      Psi_own(span(begindex_own(ell),endindex_own(ell))) = lambdasq_own(span(begindex_own(ell),endindex_own(ell))) % replace_cpp(parindextree_own[ell],Psi_own(span(begindex_own(ell+1),endindex_own(ell+1))));
      
    } 
    
    // PRODUCT-LEVEL PARAMETERS ------------------------------------- //
    
    Lambdastar(wchcross) = Psi(span(begindex(0),endindex(0))) * tausq(0);
    Lambdastar(wchown) = Psi_own(span(begindex_own(0),endindex_own(0))) * tausq_own(0);
    betabar(wchcross) = replace_cpp(parindextree[0],theta(span(begindex(1),endindex(1))));
    betabar(wchown) = replace_cpp(parindextree_own[0],theta_own(span(begindex_own(1),endindex_own(1))));
    
    for(int i=0;i<p;i++){

      // precompute
      Xt = X/sqrt(sigmasq(i));
      LamXtp = diagmat(Lambdastar(span(p*i,p*i+p-1)))*trans(Xt);
      XtLamXtp = Xt*LamXtp;
      Ystar = Y.col(i) - X*betabar(span(p*i,p*i+p-1));

      // phi (marginalizing over beta)
      // first use Woodbury to rewrite inverse as (Xt'Xt + Lam^-1)^-1 = Lam - LamXt'(I + XtLamXt')^-1XtLam
      // then projection matrix is can be written as:
      // I - Xt(Xt'Xt + Lam^-1)^-1Xt'
      // = I - Xt(Lam - LamXt'(I + XtLamXt')^-1XtLam)Xt'
      // = I - (XtLamXt' - XtLamXt'(I + XtLamXt')^-1XtLamXt'
      mat C = Clist[i];
      Ct = C/sqrt(sigmasq(i));
      irootXt = solve(trimatu(chol(symmatu(XtLamXtp + eye(n,n)))),eye(n,n));
      projmat = eye(n,n) - (XtLamXtp - XtLamXtp*irootXt*trans(irootXt)*XtLamXtp);
      irootC = solve(trimatu(chol(symmatu(trans(Ct)*projmat*Ct + Aphi(i,i)))), eye(nphi(i),nphi(i)));
      phitilde = (irootC*trans(irootC))*(trans(Ct)*projmat*Ystar/sqrt(sigmasq(i)) + Aphi(i,i)*phibar(i));
      phi(span(cumnphi(i)-nphi(i),cumnphi(i)-1)) = phitilde + irootC*randn(nphi(i));

      // Ystar
      Ystar = Y.col(i) - C*phi(span(cumnphi(i)-nphi(i),cumnphi(i)-1)) - X*betabar(span(p*i,p*i+p-1));
      ytstar = vectorise(Ystar)/sqrt(sigmasq(i));

      // beta (conditional)
      u = sqrt(Lambdastar(span(p*i,p*i+p-1))) % randn<vec>(p);
      v = Xt * u + randn<vec>(n);
      w = (irootXt*trans(irootXt)) * (ytstar - v);
      xi(span(p*i,p*i+p-1)) = u + LamXtp * w;
      beta(span(p*i,p*i+p-1)) = betabar(span(p*i,p*i+p-1)) + xi(span(p*i,p*i+p-1));

      // sigmasq
      vec E = Y.col(i) - X*beta(span(p*i,p*i+p-1)) - C*phi(span(cumnphi(i)-nphi(i),cumnphi(i)-1));
      rate = a + 0.5*trans(E)*E;
      sigmasq(i) = 1.0/R::rgamma(b+0.5*n, 1.0/as_scalar(rate));

    }
    
    // fill in first set of thetas with last draw of beta
    theta(span(begindex(0),endindex(0))) = beta(wchcross);
    theta_own(span(0,npar_own(0)-1)) = beta(wchown);
    
    // tausq - CROSS
    levelsums = sum(pow(xi(wchcross),2.0)/Psi(span(begindex(0),endindex(0))));
    tausq(0) = 1.0/R::rgamma(0.5 + 0.5*npar(0), 1.0/(1.0/xitau(0)+0.5*as_scalar(levelsums)));
    xitau(0) = 1.0/R::rgamma(1,1.0/(1+1.0/tausq(0)));
    
    // lambdasq - CROSS
    lambdasq(span(begindex(0),endindex(0))) = draw_lambdasq(0,npar,parindextree,childrencounts,theta,Psi,lambdasq,xilambda,tausq,product_shrinkage,group_shrinkage);
    xilambda(span(begindex(0),endindex(0))) = randig(npar(0),ones(npar(0)),1 + 1.0/lambdasq(span(begindex(0),endindex(0))));
    Psi(span(begindex(0),endindex(0))) = lambdasq(span(begindex(0),endindex(0))) % replace_cpp(parindextree[0],Psi(span(begindex(1),endindex(1))));

    // tausq - OWN
    levelsums = sum(pow(xi(wchown),2.0)/Psi_own(span(begindex_own(0),endindex_own(0))));
    tausq_own(0) = 1.0/R::rgamma(0.5 + 0.5*npar_own(0), 1.0/(1.0/xitau_own(0)+0.5*as_scalar(levelsums)));
    xitau_own(0) = 1.0/R::rgamma(1,1.0/(1+1.0/tausq_own(0)));
    
    // lambdasq - OWN
    lambdasq_own(span(begindex_own(0),endindex_own(0))) = draw_lambdasq(0,npar_own,parindextree_own,childrencounts_own,theta_own,Psi_own,lambdasq_own,xilambda_own,tausq_own,product_shrinkage,group_shrinkage);
    xilambda_own(span(begindex_own(0),endindex_own(0))) = randig(npar_own(0),ones(npar_own(0)),1 + 1.0/lambdasq_own(span(begindex_own(0),endindex_own(0))));
    Psi_own(span(begindex_own(0),endindex_own(0))) = lambdasq_own(span(begindex_own(0),endindex_own(0))) % replace_cpp(parindextree_own[0],Psi_own(span(begindex_own(1),endindex_own(1))));
    
    // -------------------------------------------------------------- //
    // print time and store draws
    // -------------------------------------------------------------- //

    // time
    if(print) print_time(rep, Rep, timer.toc());
      
    // store draws
    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep - 1;
      betadraws(mkeep,span::all) = trans(beta);
      lambdasqdraws(mkeep*ones<uvec>(1),wchcross) = trans(lambdasq(span(begindex(0),endindex(0))));
      lambdasqdraws(mkeep*ones<uvec>(1),wchown) = trans(lambdasq_own(span(begindex_own(0),endindex_own(0))));
      thetadraws(mkeep,span::all) = trans(theta(span(begindex(1),endindex(L)-1)));
      thetaowndraws(mkeep,span::all) = trans(theta_own(span(begindex_own(1),endindex_own(L)-1)));
      // lambdasqdraws(mkeep,span::all) = trans(lambdasq);
      // lambdasqowndraws(mkeep,span::all) = trans(lambdasq_own);
      psisqdraws(mkeep,span::all) = trans(lambdasq(span(begindex(1),endindex(L)-1)));
      psisqowndraws(mkeep,span::all) = trans(lambdasq_own(span(begindex_own(1),endindex_own(L)-1)));
      Psidraws(mkeep,span::all) = trans(Psi(span(0,sum(npar)-2)));
      Psiowndraws(mkeep,span::all) = trans(Psi_own(span(0,sum(npar_own)-2)));
      tausqdraws(mkeep,span::all) = trans(tausq);
      tausqowndraws(mkeep,span::all) = trans(tausq_own);
      sigmasqdraws(mkeep,span::all) = trans(sigmasq);
      phidraws(mkeep,span::all) = trans(phi);
    }

  }

  // print total time elapsed
  if(print){
    if(timer.toc()/60.0<1) Rprintf("Total Time Elapsed: %.1f seconds \n",timer.toc());
    else Rprintf("Total Time Elapsed: %.1f minutes \n",timer.toc()/60.0);
  }

  return List::create(
    Named("betadraws") = betadraws,
    Named("lambdasqdraws") = lambdasqdraws,
    Named("thetadraws") = thetadraws,
    Named("thetaowndraws") = thetaowndraws,
    Named("psisqdraws") = psisqdraws,
    Named("psisqowndraws") = psisqowndraws,
    // Named("lambdasqdraws") = lambdasqdraws,
    // Named("lambdasqowndraws") = lambdasqowndraws,
    Named("Psidraws") = Psidraws,
    Named("Psiowndraws") = Psiowndraws,
    Named("tausqdraws") = tausqdraws,
    Named("tausqowndraws") = tausqowndraws,
    Named("sigmasqdraws") = sigmasqdraws,
    Named("phidraws") = phidraws
  );

}



