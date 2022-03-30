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

//[[Rcpp::export]]
List countchildren_cpp(mat const& tree){
  
  int L = tree.n_cols;
  List out(L-1);
  
  for(int ell=0;ell<L-1;ell++){
    
    vec uniq = unique(tree.col(ell));
    int n = uniq.size();
    mat counts = zeros(n,n);
    vec colvec = linspace(0,L-1,L);
    
    for(int i=0;i<n;i++){
      
      mat isubtree = tree(find(tree.col(ell)==i+1),find(colvec>=ell));

      for(int j=0;j<n;j++){
        
        mat jsubtree = tree(find(tree.col(ell)==j+1),find(colvec>=ell));

        vec counti = zeros(jsubtree.n_cols);
        vec countj = zeros(jsubtree.n_cols);
        for(int k=1;k<jsubtree.n_cols;k++){
          vec uniqi = unique(isubtree.col(k));
          vec uniqj = unique(jsubtree.col(k));
          counti(k) = uniqi.size();
          countj(k) = uniqj.size();
        }
        
        counts(i,j) = sum(counti % countj);

      }
      
    }
    
    out[ell] = vectorise(counts);
    
  }
  
  
  return out;
}

//[[Rcpp::export]]
List countchildrenown_cpp(mat const& tree){
  
  int L = tree.n_cols;
  List out(L-1);
  
  for(int ell=0;ell<L-1;ell++){
    
    vec uniq = unique(tree.col(ell));
    int n = uniq.size();
    mat counts = zeros(n);
    vec colvec = linspace(0,L-1,L);
    
    for(int i=0;i<n;i++){
      counts(i) = sum(tree.col(ell)==i+1);
    }
    
    out[ell] = vectorise(counts);
    
  }
  
  
  return out;
}


//[[Rcpp::export]]
vec create_Psi_ellmone_cpp(List const& lambdalist, List const& counts, List const& list, vec const& npar, int const& level){
  vec Psi_ellmone = ones(npar(0));
  if(level>1){
    Psi_ellmone = replace_cpp(list[1],Psi_ellmone) % replace_cpp(list[1],lambdalist[0]);
    if(level>2){
      for(int j=2;j<level;j++){
        Psi_ellmone = replace_cpp(list[j],Psi_ellmone) % replace_cpp(list[j],lambdalist[j-1]);
      }
    }
  }
  return Psi_ellmone;
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
  output(find(output<1.0e-8)).fill(1.0e-8);

  return output;
}

// generate inverse gaussian random variables
//[[Rcpp::export]]
vec randinvgaussian(int const& n, vec const& mu, double const& lambda){
  vec nu = randn(n);
  vec y = pow(nu,2.0);
  vec musq = pow(mu,2.0);
  vec output = mu + musq%y/2.0/lambda - mu/2.0/lambda%sqrt(4.0*lambda*mu%y + musq%pow(y,2.0));
  output(find(output<1.0e-8)).fill(1.0e-8);
  uvec wch = find(randu(n) > mu/(mu+output));
  output(wch) = musq(wch)/output(wch);
  return output;
}

// ---------------------------------------------------------------------- //
// MCMC functions
// ---------------------------------------------------------------------- //

//[[Rcpp::export]]
vec drawbeta(mat const& Y, mat const& X, bool fast){
  
  int n = Y.n_rows;
  int p = Y.n_cols;
  vec lamstar = ones(p*p);
  vec beta = zeros(p*p);
  vec betabar = zeros(p*p);
  
  if(fast){
    
    for(int i=0;i<p;i++){

      // precompute
      mat LamXp = diagmat(lamstar(span(p*i,p*i+p-1))) * trans(X);
      mat irootXt = solve(trimatu(chol(X*LamXp + eye(n,n))),eye(n,n));

      // Ystar
      vec ytstar = vectorise(Y.col(i));

      // beta (conditional)
      vec u = betabar(span(p*i,p*i+p-1)) + sqrt(lamstar(span(p*i,p*i+p-1)))%randn<vec>(p);
      vec v = X*u + randn<vec>(n);
      vec w = (irootXt*trans(irootXt))*(ytstar - v);
      beta(span(p*i,p*i+p-1)) = u + LamXp*w;

    }
    
  }
  else{
    
    for(int i=0;i<p;i++){

      // precompute
      mat Lam = diagmat(lamstar(span(p*i,p*i+p-1)));
      mat iroot = solve(trimatu(chol(trans(X)*X + Lam)),eye(p,p));
      
      // Ystar
      vec ytstar = vectorise(Y.col(i));
      
      // beta (conditional)
      vec betatilde = (iroot*trans(iroot))*(trans(X)*ytstar + Lam*betabar(span(p*i,p*i+p-1)));
      beta(span(p*i,p*i+p-1)) = betatilde + iroot*randn<vec>(p);

    }
    
  }
  
  return beta;
}

//[[Rcpp::export]]
List rSURshrinkage(List Data, List Prior, List Mcmc, std::string Shrinkage, bool print){
  
  // data
  mat Y = Data["Y"];
  mat X = Data["X"];
  double p = X.n_cols;
  double n = X.n_rows;
  // int npar = p*p;
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
  // betabar(wchown).fill(betabarii);
  // lamstar(wchown).fill(taubarii);
  vec tau = ones(1);
  vec xitau = tau;
  vec tauown = ones(1);
  vec xitauown = tauown;
  vec lambda = ones(npar);
  vec lambdaown = ones(p);
  vec xilambda = ones(npar);
  vec xilambdaown = ones(p);
  vec sigmasq = ones(p);
  
  // storage matrices
  mat phidraws(Rep/keep,sum(nphi));
  mat betadraws(Rep/keep, p*p);
  mat lambdadraws(Rep/keep, npar);
  mat lambdaowndraws(Rep/keep, p);
  mat taudraws(Rep/keep, 1);
  mat tauowndraws(Rep/keep, 1);
  mat sigmasqdraws(Rep/keep, p);
  
  // print progress banner
  wall_clock timer;
  timer.tic();
  Datetime dt;
  if(print){
    Rprintf("MCMC Progress \n");
  }
  
  // MCMC loop
  for (rep=0; rep<Rep; rep++){
    
    // -------------------------------------------------------------- //
    // phi, beta, sigmasq
    // -------------------------------------------------------------- //
    
    lamstar(wchcross) = lambda*as_scalar(tau);
    lamstar(wchown) = lambdaown*as_scalar(tauown);
    // lamstar(wchown) = ones(p)*as_scalar(tauown);
    // lamstar = lambda*as_scalar(tau);
    
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
    
    // tau ~ C+(0,1)//
    // cross elasticities
    levelsums = sum(pow(beta(wchcross),2.0)/lambda);
    tau = 1/R::rgamma(0.5*(p*p-p+1), 1/as_scalar(1/xitau+0.5*levelsums));
    xitau = 1/R::rgamma(1,as_scalar(1/(1+1/tau)));
    // own elasticities
    levelsums = sum(pow(beta(wchown),2.0));
    tauown = 1/R::rgamma(0.5*(p+1), 1/as_scalar(1/xitauown+0.5*levelsums));
    xitauown = 1/R::rgamma(1,as_scalar(1/(1+1/tauown)));
    
    // lambda (default to ridge)
    lambda = ones(npar); 
    lambdaown = ones(p); 

    // lasso: Exp(1/2)=Gamma(1,1/2)
    if(Shrinkage=="lasso"){
      
      // cross
      mutilde = pow(2*tau(0)/pow(beta(wchcross),2.0),0.5);
      // vec mutilde = pow(2*tau(0)/pow(beta,2.0),0.5);
      lambda = 1/randinvgaussian(npar,mutilde,2); 
      
      // own
      mutilde = pow(2*tauown(0)/pow(beta(wchown),2.0),0.5);
      // vec mutilde = pow(2*tau(0)/pow(beta,2.0),0.5);
      lambdaown = 1/randinvgaussian(p,mutilde,2); 
      
    }
    // horseshoe: C+(0,1)
    if(Shrinkage=="horseshoe"){
      
      // cross
      scale = 1/xilambda + 0.5*pow(beta(wchcross),2.0)/as_scalar(tau);
      lambda = randig(npar,ones(npar),scale);
      scale = 1+1/lambda;
      xilambda = randig(npar,ones(npar),scale);
      
      // own
      scale = 1/xilambdaown + 0.5*pow(beta(wchown),2.0)/as_scalar(tauown);
      lambdaown = randig(p,ones(p),scale);
      scale = 1+1/lambdaown;
      xilambdaown = randig(p,ones(p),scale);
      
    }
    
    // -------------------------------------------------------------- //
    // print time and store draws
    // -------------------------------------------------------------- //
    
    if(print){
      // time
      if(rep==0){
        if(timer.toc()/60.0*(Rep-rep-1)/(rep+1)<1){
          Rprintf("Estimated time: %.1f seconds \n",timer.toc()*(Rep-rep-1)/(rep+1));
          Rprintf("  0%%   10   20   30   40   50   60   70   80   90   100%%\n");
          Rprintf("  |----|----|----|----|----|----|----|----|----|----|\n");
        }
        else{
          Rprintf("Estimated time: %.1f minutes \n",timer.toc()/60.0*(Rep-rep-1)/(rep+1));
          Rprintf("  0%%   10   20   30   40   50   60   70   80   90   100%%\n");
          Rprintf("  |----|----|----|----|----|----|----|----|----|----|\n");
        }
      }
      if(Rep>50){
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
    }
    
    // store draws
    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      betadraws(mkeep-1,span::all) = trans(beta);
      lambdadraws(mkeep-1,span::all) = trans(lambda);
      lambdaowndraws(mkeep-1,span::all) = trans(lambdaown);
      taudraws(mkeep-1,span::all) = tau;
      tauowndraws(mkeep-1,span::all) = tauown;
      sigmasqdraws(mkeep-1,span::all) = trans(sigmasq);
      phidraws(mkeep-1,span::all) = trans(phi);
    }
    
  }
  
  // print total time elapsed
  if(print){
    if(timer.toc()/60.0<1){
      Rprintf("Total Time Elapsed: %.1f seconds \n",timer.toc());
    }
    else{
      Rprintf("Total Time Elapsed: %.1f minutes \n",timer.toc()/60.0);
    }
  }
  
  return List::create(
    Named("betadraws") = betadraws,
    Named("lambdadraws") = lambdadraws,
    Named("lambdaowndraws") = lambdaowndraws,
    Named("taudraws") = taudraws,
    Named("tauowndraws") = tauowndraws,
    Named("sigmasqdraws") = sigmasqdraws,
    Named("phidraws") = phidraws
  );
  
}

// function to sample upper-level thetas
vec drawtheta(int const& ell,List const& thetalist, double const& thetabar, List const& lambdalist,
              vec const& Psi_ell,vec const& Psi_ellmone,vec const& tau,vec const& npar, List const& list){
  
  vec denom, ellmone, ellpone_sums, v_kl, thetatilde, theta;
  
  // posterior variance
  denom = as<vec>(lambdalist[ell+1]) % Psi_ell * tau(ell+1);
  ellpone_sums = groupsum_cpp(1/denom,list[ell+1],npar(ell));
  if(ell>0){
    ellmone = 1/(as<vec>(lambdalist[ell])*tau(ell));
    ellmone = ellmone % replace_cpp(list[ell],1.0/Psi_ellmone);
  }
  else{
    ellmone = ones(npar(ell))*1/tau(ell);
  }
  v_kl = 1/(ellpone_sums + ellmone);
  
  // posterior mean
  ellpone_sums = groupsum_cpp(as<vec>(thetalist[ell+1])/denom,list[ell+1],npar(ell));
  if(ell>0){
    ellmone = replace_cpp(list[ell],thetalist[ell-1])/(as<vec>(lambdalist[ell])*tau(ell));
  }
  else{
    ellmone = ones(npar(ell))*thetabar/tau(ell);
  }
  thetatilde = v_kl % (ellpone_sums + ellmone);
  
  // draw theta
  theta = thetatilde + sqrt(v_kl) % randn(npar(ell));
  
  return theta;
  
}

// function to sample upper-level lambdas
vec drawlambda(int const& ell, List const& thetalist, double const& thetabar, List const& lambdalist, List const& xilambdalist,
               List const& childrencounts, vec const& Psi_ell, vec const& Psi_ellmone, vec const& tau,vec const& npar, List const& list,
               bool propagate, std::string const& group_shrinkage){
  
  vec counts, mean, levelsums, sums, diffsums, scale, Psiellmone;
  
  vec lambda = ones(npar(ell));
  vec xilambda = xilambdalist[ell];
  int L = thetalist.size();
  
  // shape: count children nodes at each level
  counts = ones(npar[ell]);
  if(propagate) counts = as<vec>(childrencounts[ell]);
  
  // rate: sums at level ell
  if(propagate && ell>0){
    Psiellmone = replace_cpp(list[ell],Psi_ellmone) % replace_cpp(list[ell],lambdalist[ell-1]);
  }
  else Psiellmone = ones(npar[ell]);
  if(ell>0) mean = replace_cpp(list[ell],thetalist[ell-1]);
  else mean = replace_cpp(list[ell],thetabar*ones(1));
  levelsums = pow(as<vec>(thetalist[ell])-mean,2.0)/(Psiellmone*tau(ell));
  
  if(ell>0){
    
    if(propagate){
      
      // sums at levels greater than ell
      for(int s=(ell+1);s<L;s++){
        
        // compute Psi minus one for level s
        Psiellmone = replace_cpp(list[s],Psiellmone) % replace_cpp(list[s],lambdalist[s-1]);

        // compute sums of scaled differences in theta
        mean = replace_cpp(list[s],thetalist[s-1]);
        sums = pow(as<vec>(thetalist[s])-mean,2.0)/(Psiellmone*tau(s));
        diffsums = groupsum_cpp(sums,replace_cpp(list[s],list[ell]),npar(ell));
        levelsums = levelsums + replace_cpp(list[ell],diffsums);
        
      }
    }
    
    // horseshoe: C+(0,1)
    if(group_shrinkage=="horseshoe"){
      scale = 1.0/xilambda + 0.5*levelsums;
      lambda = randig(npar(ell),0.5+0.5*counts,scale);
    }
    
  }
  
  return lambda;

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
  mat tree = Data["tree"];
  List childrencounts = Data["childrencounts"];
  List list = Data["list"];
  List list_own = Data["list_own"];
  int L = tree.n_cols;

  List childrencounts_own = countchildrenown_cpp(tree);
  
  // prior
  double thetabar_cross = Prior["thetabar_cross"];
  double thetabar_own = Prior["thetabar_own"];
  mat Aphi = Prior["Aphi"];
  vec phibar = Prior["phibar"];
  double a = Prior["a"];
  double b = Prior["b"];

  // mcmc
  int Rep = Mcmc["R"];
  int initial_run = Mcmc["initial_run"];
  int RepRun = Rep + initial_run;
  int keep = Mcmc["keep"];

  // shrinkage
  std::string product_shrinkage = Shrinkage["product"];
  std::string group_shrinkage = Shrinkage["group"];

  // initialize
  int rep, mkeep;
  vec ytstar, phitilde, Psi_Lmone, Psi_Lmone_own, Psi_ellmone, Psi_ellmone_own, Psi_ell, Psi_ell_own, ellpone_sums, ellmone, denom, counts,
  levelsums, diffsums, v_kl, thetatilde, theta, mean, u, v, vari, w, sums, rate, scale,
  lambda, xilambda, lambda_own, xilambda_own, mutilde;
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
  vec beta = vectorise(inv(trans(X)*X+0.1*eye(p,p))*trans(X)*Y);
  vec xi = zeros(p*p);
  vec betabar = zeros(p*p);
  vec lamstar = ones(p*p);
  // cross
  vec tau = ones(L);
  vec xitau = tau;
  List lambdalist = List(L);
  List xilambdalist = List(L);
  List thetalist = List(L);
  // own
  vec tau_own = ones(L);
  vec xitau_own = tau_own;
  List lambdalist_own = List(L);
  List xilambdalist_own = List(L);
  List thetalist_own = List(L);

  thetalist[L-1] = zeros(npar[L-1]);
  thetalist_own[L-1] = zeros(npar_own[L-1]);
  for(int ell=L-1;ell>=0;ell--){
    if(ell<L-1){
      thetalist[ell] = zeros(npar(ell));
      thetalist_own[ell] = zeros(npar_own(ell));
    }
    lambdalist[ell] = ones(npar(ell));
    xilambdalist[ell] = ones(npar(ell));
    lambdalist_own[ell] = ones(npar_own(ell));
    xilambdalist_own[ell] = ones(npar_own(ell));
  }
  vec sigmasq = ones(p);

  // storage matrices
  mat phidraws(Rep/keep,sum(nphi));
  mat betadraws(Rep/keep, p*p);
  mat thetadraws(Rep/keep,sum(npar(span(0,L-2))));
  mat thetaowndraws(Rep/keep,sum(npar_own(span(0,L-2))));
  mat lambdadraws(Rep/keep,sum(npar));
  mat lambdaowndraws(Rep/keep,sum(npar_own));
  mat taudraws(Rep/keep, L);
  mat tauowndraws(Rep/keep, L);
  mat sigmasqdraws(Rep/keep, p);

  // print progress banner
  wall_clock timer;
  timer.tic();
  if(print){
    Rprintf(" MCMC Progress \n");
  }

  bool propagate = TRUE;
  // bool propagate = FALSE;
  double numcutoff = 1.0e-4;

  // MCMC loop
  for (rep=0; rep<RepRun; rep++){

    // -------------------------------------------------------------- //
    // higher-level effects
    // -------------------------------------------------------------- //

    // for(int ell=0;ell<L-1;ell++){
    // 
    //   // theta (own) ---------------------------------------- //
    //   Psi_ellmone_own = ones(npar_own[ell]);
    //   if(propagate) Psi_ellmone_own = create_Psi_ellmone_cpp(lambdalist_own,childrencounts_own,list_own,npar_own,ell+1);
    //   Psi_ellmone_own(find(Psi_ellmone_own<numcutoff)).fill(numcutoff);
    //   Psi_ell_own = replace_cpp(list_own[ell+1],Psi_ellmone) % replace_cpp(list_own[ell+1],lambdalist_own[ell]);
    //   Psi_ell_own(find(Psi_ell_own<numcutoff)).fill(numcutoff);
    //   theta = drawtheta(ell,thetalist_own, thetabar_own, lambdalist_own, Psi_ell_own, Psi_ellmone_own, tau_own, npar_own, list_own);
    //   thetalist_own[ell] = theta;
    // 
    //   // tau (own) //
    //   if(ell>0) mean = replace_cpp(list_own[ell],thetalist_own[ell-1]);
    //   else mean = replace_cpp(list_own[ell],thetabar_own*ones(1));
    //   levelsums = sum(pow(as<vec>(thetalist_own[ell]) - mean,2.0)/(as<vec>(lambdalist_own[ell])));
    //   tau_own(ell) = 1.0/R::rgamma(0.5*(npar_own(ell)+1), 1.0/(1.0/xitau_own(ell)+0.5*as_scalar(levelsums)));
    //   if(tau_own(ell)<numcutoff) tau_own(ell)=numcutoff;
    //   xitau_own(ell) = 1.0/R::rgamma(1,1.0/(1+1.0/tau_own(ell)));
    //   if(xitau_own(ell)<numcutoff) xitau_own(ell)=numcutoff;
    // 
    //   // lambda (own)//
    //   if(group_shrinkage!="ridge" && rep>=initial_run){
    // 
    //     // draw lambda
    //     lambda = drawlambda(ell, thetalist_own, thetabar_own, lambdalist_own, xilambdalist_own,
    //                         childrencounts_own, Psi_ell_own, Psi_ellmone_own, tau_own, npar_own, list_own,
    //                         propagate, group_shrinkage="horseshoe");
    //     scale = 1 + 1.0/lambda;
    //     xilambda = randig(npar_own(ell),ones(npar_own(ell)),scale);
    // 
    //     // store draws
    //     lambdalist_own[ell] = lambda;
    //     xilambdalist_own[ell] = xilambda;
    // 
    //   }
    //   
    //   // theta (cross) ---------------------------------------- //
    //   Psi_ellmone = ones(npar[ell]); ////////////////////////////////////
    //   if(propagate) Psi_ellmone = create_Psi_ellmone_cpp(lambdalist,childrencounts,list,npar,ell+1);
    //   Psi_ellmone(find(Psi_ellmone<numcutoff)).fill(numcutoff);
    //   Psi_ell = replace_cpp(list[ell+1],Psi_ellmone) % replace_cpp(list[ell+1],lambdalist[ell]);
    //   Psi_ell(find(Psi_ell<numcutoff)).fill(numcutoff);
    //   theta = drawtheta(ell,thetalist, thetabar_cross, lambdalist, Psi_ell, Psi_ellmone, tau, npar, list);
    //   thetalist[ell] = theta;
    // 
    //   // tau (cross) //
    //   if(ell>0) mean = replace_cpp(list[ell],thetalist[ell-1]);
    //   else mean = replace_cpp(list[ell],thetabar_cross*ones(1));
    //   levelsums = sum(pow(as<vec>(thetalist[ell]) - mean,2.0)/(Psi_ellmone % as<vec>(lambdalist[ell])));
    //   tau(ell) = 1.0/R::rgamma(0.5*(npar(ell)+1), 1.0/(1.0/xitau(ell)+0.5*as_scalar(levelsums)));
    //   if(tau(ell)<numcutoff) tau(ell)=numcutoff;
    //   xitau(ell) = 1.0/R::rgamma(1,1.0/(1+1.0/tau(ell)));
    //   if(xitau(ell)<numcutoff) xitau(ell)=numcutoff;
    //   
    //   // lambda (cross)//
    //   if(group_shrinkage!="ridge" && rep>=initial_run){
    //     
    //     // draw lambda
    //     lambda = drawlambda(ell, thetalist, thetabar_cross, lambdalist, xilambdalist,
    //                         childrencounts, Psi_ell, Psi_ellmone, tau, npar,list,
    //                         propagate, group_shrinkage="horseshoe");
    //     scale = 1 + 1.0/lambda;
    //     xilambda = randig(npar(ell),ones(npar(ell)),scale);
    //     
    //     // store draws
    //     lambdalist[ell] = lambda;
    //     xilambdalist[ell] = xilambda;
    //     
    //   }
    //   
    // }

    // -------------------------------------------------------------- //
    // product-level elasticities + phi + sigmasq
    // -------------------------------------------------------------- //

    if(propagate){
      Psi_Lmone = create_Psi_ellmone_cpp(lambdalist,childrencounts,list,npar,L);
      Psi_Lmone(find(Psi_Lmone<numcutoff)).fill(numcutoff);
      Psi_Lmone_own = create_Psi_ellmone_cpp(lambdalist_own,childrencounts_own,list_own,npar_own,L);
      Psi_Lmone_own(find(Psi_Lmone_own<numcutoff)).fill(numcutoff);
    }
    else{
      Psi_Lmone = ones(npar[L-1]);
      Psi_Lmone_own = ones(npar_own[L-1]);
    }
    betabar(wchown) = replace_cpp(list_own[L-1],thetalist_own[L-2]);
    betabar(wchcross) = replace_cpp(list[L-1],thetalist[L-2]);
    lamstar(wchown) = as<vec>(lambdalist_own[L-1]) % Psi_Lmone_own * tau_own(L-1);
    lamstar(wchcross) = as<vec>(lambdalist[L-1]) % Psi_Lmone * tau(L-1);

    for(int i=0;i<p;i++){

      // precompute
      Xt = X/sqrt(sigmasq(i));
      LamXtp = diagmat(lamstar(span(p*i,p*i+p-1)))*trans(Xt);
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
      u = sqrt(lamstar(span(p*i,p*i+p-1))) % randn<vec>(p);
      v = Xt * u + randn<vec>(n);
      w = (irootXt*trans(irootXt)) * (ytstar - v);
      xi(span(p*i,p*i+p-1)) = u + LamXtp * w;
      beta(span(p*i,p*i+p-1)) = betabar(span(p*i,p*i+p-1)) + xi(span(p*i,p*i+p-1));

      // sigmasq
      vec E = Y.col(i) - X*beta(span(p*i,p*i+p-1)) - C*phi(span(cumnphi(i)-nphi(i),cumnphi(i)-1));
      rate = a + 0.5*trans(E)*E;
      sigmasq(i) = 1.0/R::rgamma(b+0.5*n, 1.0/as_scalar(rate));

    }

    // save elasticities
    thetalist[L-1] = vectorise(beta(wchcross));
    thetalist_own[L-1] = vectorise(beta(wchown));

    // tau (cross)//
    levelsums = sum(pow(xi(wchcross),2.0)/Psi_Lmone/as<vec>(lambdalist[L-1]));
    tau(L-1) = 1.0/R::rgamma(0.5*(npar(L-1)+1), 1.0/(1.0/xitau(L-1)+0.5*as_scalar(levelsums)));
    if(tau(L-1)<numcutoff) tau(L-1)=numcutoff;
    xitau(L-1) = 1.0/R::rgamma(1,1.0/(1+1.0/tau(L-1)));
    if(xitau(L-1)<numcutoff) xitau(L-1)=numcutoff;
    
    // tau (own)//
    levelsums = sum(pow(xi(wchown),2.0)/as<vec>(lambdalist_own[L-1]));
    tau_own(L-1) = 1.0/R::rgamma(0.5*(npar_own(L-1)+1), 1.0/(1.0/xitau_own(L-1)+0.5*as_scalar(levelsums)));
    if(tau_own(L-1)<numcutoff) tau_own(L-1)=numcutoff;
    xitau_own(L-1) = 1.0/R::rgamma(1,1.0/(1+1.0/tau_own(L-1)));
    if(xitau_own(L-1)<numcutoff) xitau_own(L-1)=numcutoff;
    
    // lambda (default to ridge)
    lambda = ones(npar(L-1));
    xilambda = as<vec>(xilambdalist[L-1]);
    lambda_own = ones(npar_own(L-1));
    xilambda_own = as<vec>(xilambdalist_own[L-1]);

    // lasso: Exp(1/2)=Gamma(1,1/2)
    if(product_shrinkage=="lasso" && rep>=initial_run){
      // cross
      mutilde = sqrt(2.0)*sqrt(Psi_Lmone)*sqrt(tau(L-1))/abs(xi(wchcross));
      lambda = 1/randinvgaussian(npar(L-1),mutilde,2);
      lambdalist[L-1] = lambda;

      // own
      mutilde = sqrt(2.0)*sqrt(tau_own(L-1))/abs(xi(wchown));
      lambda_own = 1/randinvgaussian(npar_own(L-1),mutilde,2);
      lambdalist_own[L-1] = lambda_own;
    }

    // horseshoe: C+(0,1)
    if(product_shrinkage=="horseshoe" && rep>=initial_run){

      // cross
      scale = 1.0/xilambda + 0.5*pow(xi(wchcross),2.0)/Psi_Lmone/as_scalar(tau(L-1));
      lambda = randig(npar(L-1),ones(npar(L-1)),scale);
      scale = 1 + 1.0/lambda;
      xilambda = randig(npar(L-1),ones(npar(L-1)),scale);
      lambdalist[L-1] = lambda;
      xilambdalist[L-1] = xilambda;

      // own
      scale = 1.0/xilambda_own + 0.5*pow(xi(wchown),2.0)/as_scalar(tau_own(L-1));
      lambda_own = randig(npar_own(L-1),ones(npar_own(L-1)),scale);
      scale = 1 + 1.0/lambda_own;
      xilambda = randig(npar_own(L-1),ones(npar_own(L-1)),scale);
      lambdalist_own[L-1] = lambda_own;
      xilambdalist_own[L-1] = xilambda_own;

    }

    // -------------------------------------------------------------- //
    // print time and store draws
    // -------------------------------------------------------------- //

    if(print){
      // time
      if(rep==0){
        if(timer.toc()/60.0*(RepRun-rep-1)/(rep+1)<1){
          Rprintf("Estimated time: %.1f seconds \n",timer.toc()*(RepRun-rep-1)/(rep+1));
          Rprintf("  0%%   10   20   30   40   50   60   70   80   90   100%%\n");
          Rprintf("  |----|----|----|----|----|----|----|----|----|----|\n");
        }
        else{
          Rprintf("Estimated time: %.1f minutes \n",timer.toc()/60.0*(RepRun-rep-1)/(rep+1));
          Rprintf("  0%%   10   20   30   40   50   60   70   80   90   100%%\n");
          Rprintf("  |----|----|----|----|----|----|----|----|----|----|\n");
        }
      }
      if(RepRun>50){
        if ((rep+1)%(RepRun/50)==0){
          if(rep+1==RepRun/50){
            Rprintf("  *");
          }
          else if(rep<RepRun-1){
            Rprintf("*");
          }
          else Rprintf("*\n");
        }
      }
    }

    // store draws
    if(rep>=initial_run && (rep+1)%keep==0){
      mkeep = (rep-initial_run+1)/keep;
      betadraws(mkeep-1,span::all) = trans(beta);
      thetadraws(mkeep-1,span::all) = trans(unlist(thetalist,sum(npar(span(0,L-2))),L-1));
      thetaowndraws(mkeep-1,span::all) = trans(unlist(thetalist_own,sum(npar_own(span(0,L-2))),L-1));
      lambdadraws(mkeep-1,span::all) = trans(unlist(lambdalist,sum(npar),lambdalist.size()));
      lambdaowndraws(mkeep-1,span::all) = trans(unlist(lambdalist_own,sum(npar_own),lambdalist_own.size()));
      taudraws(mkeep-1,span::all) = trans(tau);
      tauowndraws(mkeep-1,span::all) = trans(tau_own);
      sigmasqdraws(mkeep-1,span::all) = trans(sigmasq);
      phidraws(mkeep-1,span::all) = trans(phi);
    }

  }

  // print total time elapsed
  if(print){
    if(timer.toc()/60.0<1){
      Rprintf("Total Time Elapsed: %.1f seconds \n",timer.toc());
    }
    else{
      Rprintf("Total Time Elapsed: %.1f minutes \n",timer.toc()/60.0);
    }
  }

  return List::create(
    Named("betadraws") = betadraws,
    Named("thetadraws") = thetadraws,
    Named("thetaowndraws") = thetaowndraws,
    Named("lambdadraws") = lambdadraws,
    Named("lambdaowndraws") = lambdaowndraws,
    Named("taudraws") = taudraws,
    Named("tauowndraws") = tauowndraws,
    Named("sigmasqdraws") = sigmasqdraws,
    Named("phidraws") = phidraws
  );

}



// //[[Rcpp::export]]
// List rSURhiershrinkage(List const& Data, List const& Prior, List const& Mcmc, List const& Shrinkage, bool print){
//   
//   // data
//   mat Y = Data["Y"];
//   mat X = Data["X"];
//   double p = X.n_cols;
//   double n = X.n_rows;
//   List Clist = Data["Clist"];
//   vec npar = Data["npar"];
//   vec npar_own = Data["npar_own"];
//   mat tree = Data["tree"];
//   List childrencounts = Data["childrencounts"];
//   List list = Data["list"];
//   List list_own = Data["list_own"];
//   int L = tree.n_cols;
//   
//   List childrencounts_own = countchildrenown_cpp(tree);
//   
//   // prior
//   double thetabar_cross = Prior["thetabar_cross"];
//   double thetabar_own = Prior["thetabar_own"];
//   mat Aphi = Prior["Aphi"];
//   vec phibar = Prior["phibar"];
//   double a = Prior["a"];
//   double b = Prior["b"];
//   
//   // mcmc
//   int Rep = Mcmc["R"];
//   int initial_run = Mcmc["initial_run"];
//   int RepRun = Rep + initial_run;
//   int keep = Mcmc["keep"];
//   
//   // shrinkage
//   std::string product_shrinkage = Shrinkage["product"];
//   std::string group_shrinkage = Shrinkage["group"];
//   
//   // initialize
//   int rep, mkeep;
//   vec ytstar, phitilde, Psi_Lmone, Psi_ellmone, Psi_ellmone_own, Psi_ell, Psi_ell_own, ellpone_sums, ellmone, denom, counts,
//   levelsums, diffsums, v_kl, thetatilde, theta, mean, u, v, vari, w, sums, rate, scale,
//   lambda, xilambda, lambda_own, xilambda_own, mutilde;
//   mat Ystar, Xt, CtpCt, LamXtp, Xpy, XtLamXtp, Lamibetabar, irootX, irootXt, irootC, Ct, projmat, CIprojCp;
//   mat XpX = trans(X)*X;
//   vec nphi = zeros(p);
//   mat Cstar;
//   for(int i=0;i<p;i++){
//     mat C = Clist[i];
//     nphi(i) = C.n_cols;
//     Cstar = join_rows(Cstar, C);
//   }
//   vec cumnphi = cumsum(nphi);
//   mat CspCs = trans(Cstar)*Cstar;
//   
//   // initial values
//   vec phi = zeros(sum(nphi));
//   uvec wchown = find(eye<mat>(p,p)==1);
//   uvec wchcross = find(eye<mat>(p,p)==0);
//   vec beta = vectorise(inv(trans(X)*X+0.1*eye(p,p))*trans(X)*Y);
//   vec xi = zeros(p*p);
//   vec betabar = zeros(p*p);
//   vec lamstar = ones(p*p);
//   // cross
//   vec tau = ones(L);
//   vec xitau = tau;
//   List lambdalist = List(L);
//   List xilambdalist = List(L);
//   List thetalist = List(L);
//   // own
//   vec tau_own = ones(L);
//   vec xitau_own = tau_own;
//   List lambdalist_own = List(L);
//   List xilambdalist_own = List(L);
//   List thetalist_own = List(L);
//   
//   // thetalist[L-1] = vectorise(beta(wchcross));
//   // thetalist_own[L-1] = vectorise(beta(wchown));
//   thetalist[L-1] = zeros(npar[L-1]);
//   thetalist_own[L-1] = zeros(npar_own[L-1]);
//   for(int ell=L-1;ell>=0;ell--){
//     if(ell<L-1){
//       thetalist[ell] = zeros(npar(ell));
//       thetalist_own[ell] = zeros(npar_own(ell));
//       // if(group_shrinkage=="sparse"){
//       //   thetalist[ell] = zeros(npar(ell));
//       //   thetalist_own[ell] = zeros(npar_own(ell));
//       // }
//       // else{
//       //   thetalist[ell] = groupmean_cpp(thetalist[ell+1],list[ell+1],npar[ell]);
//       //   thetalist_own[ell] = groupmean_cpp(thetalist_own[ell+1],list_own[ell],npar_own(ell));
//       // }
//     }
//     lambdalist[ell] = ones(npar(ell));
//     xilambdalist[ell] = ones(npar(ell));
//     lambdalist_own[ell] = ones(npar_own(ell));
//     xilambdalist_own[ell] = ones(npar_own(ell));
//   }
//   vec sigmasq = ones(p);
//   
//   // storage matrices
//   mat phidraws(Rep/keep,sum(nphi));
//   mat betadraws(Rep/keep, p*p);
//   mat thetadraws(Rep/keep,sum(npar(span(0,L-2))));
//   mat thetaowndraws(Rep/keep,sum(npar_own(span(0,L-2))));
//   mat lambdadraws(Rep/keep,sum(npar));
//   mat lambdaowndraws(Rep/keep,sum(npar_own));
//   mat taudraws(Rep/keep, L);
//   mat tauowndraws(Rep/keep, L);
//   mat sigmasqdraws(Rep/keep, p);
//   
//   // print progress banner
//   wall_clock timer;
//   timer.tic();
//   if(print){
//     Rprintf(" MCMC Progress \n");
//   }
//   
//   bool propagate = TRUE;
//   // bool propagate = FALSE;
//   
//   // MCMC loop
//   for (rep=0; rep<RepRun; rep++){
//     
//     // -------------------------------------------------------------- //
//     // higher-level effects
//     // -------------------------------------------------------------- //
//     
//     // int ell=L-2;ell>=0;ell--
//     for(int ell=0;ell<L-1;ell++){
//       
//       // theta (own effects) ---------------------------------------- //
//       
//       // product of previous lambda parameters
//       Psi_ellmone_own = ones(npar_own[ell]); ////////////////////////////////////
//       if(propagate) Psi_ellmone_own = create_Psi_ellmone_cpp(lambdalist_own,childrencounts_own,list_own,npar_own,ell+1);
//       Psi_ell_own = replace_cpp(list_own[ell+1],Psi_ellmone) % replace_cpp(list_own[ell+1],lambdalist_own[ell]);
// 
//       // posterior variance
//       denom = as<vec>(lambdalist_own[ell+1]) % Psi_ell_own * tau_own(ell+1);
//       ellpone_sums = groupsum_cpp(1/denom,list_own[ell+1],npar_own(ell));
//       if(ell>0){
//         ellmone = 1/(as<vec>(lambdalist_own[ell])*tau_own(ell));
//         ellmone = ellmone % replace_cpp(list_own[ell],1.0/Psi_ellmone_own);
//       }
//       else{
//         ellmone = ones(npar_own(ell))*1/tau_own(ell);
//       }
//       v_kl = 1/(ellpone_sums + ellmone);
// 
//       // posterior mean
//       ellpone_sums = groupsum_cpp(as<vec>(thetalist_own[ell+1])/denom,list_own[ell+1],npar_own(ell));
//       if(ell>0){
//         ellmone = replace_cpp(list_own[ell],thetalist_own[ell-1])/(as<vec>(lambdalist_own[ell])*tau_own(ell));
//       }
//       else{
//         ellmone = ones(npar_own(ell))*thetabar_own/tau_own(ell);
//       }
//       thetatilde = v_kl % (ellpone_sums + ellmone);
// 
//       // draw theta
//       theta = thetatilde + sqrt(v_kl) % randn(npar_own(ell));
// 
//       // store draw
//       thetalist_own[ell] = theta;
// 
//       // tau //
//       if(ell>0) mean = replace_cpp(list_own[ell],thetalist_own[ell-1]);
//       else mean = replace_cpp(list_own[ell],thetabar_own*ones(1));
//       levelsums = sum(pow(as<vec>(thetalist_own[ell]) - mean,2.0)/(as<vec>(lambdalist_own[ell])));
//       tau_own(ell) = 1.0/R::rgamma(0.5*(npar_own(ell)+1), 1.0/(1.0/xitau_own(ell)+0.5*as_scalar(levelsums)));
//       xitau_own(ell) = 1.0/R::rgamma(1,1.0/(1+1.0/tau_own(ell)));
//       
//       // theta (cross effects) ---------------------------------------- //
//       
//       // product of previous lambda parameters
//       Psi_ellmone = ones(npar[ell]); ////////////////////////////////////
//       if(propagate) Psi_ellmone = create_Psi_ellmone_cpp(lambdalist,childrencounts,list,npar,ell+1);
//       Psi_ell = replace_cpp(list[ell+1],Psi_ellmone) % replace_cpp(list[ell+1],lambdalist[ell]);
//       
//       // posterior variance
//       denom = as<vec>(lambdalist[ell+1]) % Psi_ell * tau(ell+1);
//       ellpone_sums = groupsum_cpp(1/denom,list[ell+1],npar(ell));
//       if(ell>0){
//         ellmone = 1/(as<vec>(lambdalist[ell])*tau(ell));
//         ellmone = ellmone % replace_cpp(list[ell],1.0/Psi_ellmone);
//       }
//       else{
//         ellmone = ones(npar(ell)) * 1/tau(ell);
//       }
//       v_kl = 1/(ellpone_sums + ellmone);
//       
//       // posterior mean
//       ellpone_sums = groupsum_cpp(as<vec>(thetalist[ell+1])/denom,list[ell+1],npar(ell));
//       if(ell>0){
//         ellmone = replace_cpp(list[ell],thetalist[ell-1])/(as<vec>(lambdalist[ell])*tau(ell));
//         ellmone = ellmone % replace_cpp(list[ell],1.0/Psi_ellmone);
//       }
//       else{
//         ellmone = ones(npar(ell)) * thetabar_cross/tau(ell);
//       }
//       thetatilde = v_kl % (ellpone_sums + ellmone);
//       
//       // draw theta
//       theta = thetatilde + sqrt(v_kl) % randn(npar(ell));
//       
//       // store draw
//       thetalist[ell] = theta;
//       
//       // tau //
//       if(ell>0) mean = replace_cpp(list[ell],thetalist[ell-1]);
//       else mean = replace_cpp(list[ell],thetabar_cross*ones(1));
//       levelsums = sum(pow(as<vec>(thetalist[ell]) - mean,2.0)/(Psi_ellmone % as<vec>(lambdalist[ell])));
//       tau(ell) = 1.0/R::rgamma(0.5*(npar(ell)+1), 1.0/(1.0/xitau(ell)+0.5*as_scalar(levelsums)));
//       xitau(ell) = 1.0/R::rgamma(1,1.0/(1+1.0/tau(ell)));
//       
//       // lambda //
//       
//       if(group_shrinkage!="ridge"){
//         
//         // local shrinkage active only for ell>1 (otherwise fix lambda=1)
//         if(ell>0){
//           
//           // shape: count children nodes at each level
//           counts = ones(npar[ell]); ////////////////////////////////////
//           if(propagate) counts = as<vec>(childrencounts[ell]);
//           
//           // rate: sums at level ell
//           if(propagate && ell>0){
//             Psi_ellmone = replace_cpp(list[ell],Psi_ellmone) % replace_cpp(list[ell],lambdalist[ell-1]);
//           }
//           else Psi_ellmone = ones(npar[ell]); ////////////////////////////////////
//           if(ell>0) mean = replace_cpp(list[ell],thetalist[ell-1]);
//           else mean = replace_cpp(list[ell],thetabar_cross*ones(1));
//           levelsums = pow(as<vec>(thetalist[ell])-mean,2.0)/(Psi_ellmone*tau(ell));
//           
//           // mean = replace_cpp(list[ell],thetalist[ell-1]);
//           // levelsums = pow(as<vec>(thetalist[ell])-mean,2.0)/(Psi_ellmone*tau(ell));
//           
//           if(propagate){
//             
//             // sums at levels greater than ell ////////////////////////////////////
//             for(int s=(ell+1);s<L;s++){
//               
//               // compute Psi minus one for level s
//               Psi_ellmone = replace_cpp(list[s],Psi_ellmone) % replace_cpp(list[s],lambdalist[s-1]);
// 
//               // compute sums of scaled differences in theta
//               mean = replace_cpp(list[s],thetalist[s-1]);
//               sums = pow(as<vec>(thetalist[s])-mean,2.0)/(Psi_ellmone*tau(s));
//               diffsums = groupsum_cpp(sums,replace_cpp(list[s],list[ell]),npar(ell));
//               levelsums = levelsums + replace_cpp(list[ell],diffsums);
//               
//             }
//           }
//           
//           lambda = ones(npar(ell));
//           xilambda = as<vec>(xilambdalist[ell]);
//           
//           // horseshoe: C+(0,1)
//           if(group_shrinkage=="horseshoe" && rep>=initial_run){
//             scale = 1.0/xilambda + 0.5*levelsums;
//             lambda = randig(npar(ell),0.5+0.5*counts,scale);
//             scale = 1 + 1.0/lambda;
//             xilambda = randig(npar(ell),ones(npar(ell)),scale);
//           }
//           
//           // store draws
//           lambdalist[ell] = lambda;
//           xilambdalist[ell] = xilambda;
//           
//         }
//         
//       }
//       
//     }
//     
//     // -------------------------------------------------------------- //
//     // product-level elasticities + phi + sigmasq
//     // -------------------------------------------------------------- //
//     
//     if(propagate){
//       Psi_Lmone = create_Psi_ellmone_cpp(lambdalist,childrencounts,list,npar,L);
//     }
//     else Psi_Lmone = ones(npar[L-1]); ////////////////////////////////////
//     betabar(wchown) = replace_cpp(list_own[L-1],thetalist_own[L-2]);
//     lamstar(wchown) = as<vec>(lambdalist_own[L-1]) * tau_own(L-1);
//     betabar(wchcross) = replace_cpp(list[L-1],thetalist[L-2]);
//     lamstar(wchcross) = as<vec>(lambdalist[L-1]) % Psi_Lmone * tau(L-1);
//     
//     for(int i=0;i<p;i++){
//       
//       // precompute
//       Xt = X/sqrt(sigmasq(i));
//       LamXtp = diagmat(lamstar(span(p*i,p*i+p-1)))*trans(Xt);
//       XtLamXtp = Xt*LamXtp;
//       Ystar = Y.col(i) - X*betabar(span(p*i,p*i+p-1));
//       
//       // phi (marginalizing over beta)
//       // first use Woodbury to rewrite inverse as (Xt'Xt + Lam^-1)^-1 = Lam - LamXt'(I + XtLamXt')^-1XtLam
//       // then projection matrix is can be written as:
//       // I - Xt(Xt'Xt + Lam^-1)^-1Xt'
//       // = I - Xt(Lam - LamXt'(I + XtLamXt')^-1XtLam)Xt'
//       // = I - (XtLamXt' - XtLamXt'(I + XtLamXt')^-1XtLamXt'
//       mat C = Clist[i];
//       Ct = C/sqrt(sigmasq(i));
//       irootXt = solve(trimatu(chol(symmatu(XtLamXtp + eye(n,n)))),eye(n,n));
//       projmat = eye(n,n) - (XtLamXtp - XtLamXtp*irootXt*trans(irootXt)*XtLamXtp);
//       irootC = solve(trimatu(chol(symmatu(trans(Ct)*projmat*Ct + Aphi(i,i)))), eye(nphi(i),nphi(i)));
//       phitilde = (irootC*trans(irootC))*(trans(Ct)*projmat*Ystar/sqrt(sigmasq(i)) + Aphi(i,i)*phibar(i));
//       phi(span(cumnphi(i)-nphi(i),cumnphi(i)-1)) = phitilde + irootC*randn(nphi(i));
//       
//       // Ystar
//       Ystar = Y.col(i) - C*phi(span(cumnphi(i)-nphi(i),cumnphi(i)-1)) - X*betabar(span(p*i,p*i+p-1));
//       ytstar = vectorise(Ystar)/sqrt(sigmasq(i));
//       
//       // beta (conditional)
//       u = sqrt(lamstar(span(p*i,p*i+p-1))) % randn<vec>(p);
//       v = Xt * u + randn<vec>(n);
//       w = (irootXt*trans(irootXt)) * (ytstar - v);
//       xi(span(p*i,p*i+p-1)) = u + LamXtp * w;
//       beta(span(p*i,p*i+p-1)) = betabar(span(p*i,p*i+p-1)) + xi(span(p*i,p*i+p-1));
//       
//       // sigmasq
//       vec E = Y.col(i) - X*beta(span(p*i,p*i+p-1)) - C*phi(span(cumnphi(i)-nphi(i),cumnphi(i)-1));
//       rate = a + 0.5*trans(E)*E;
//       sigmasq(i) = 1.0/R::rgamma(b+0.5*n, 1.0/as_scalar(rate));
//       
//     }
//     
//     // save elasticities
//     thetalist[L-1] = vectorise(beta(wchcross));
//     thetalist_own[L-1] = vectorise(beta(wchown));
//     
//     // tau (cross)//
//     levelsums = sum(pow(xi(wchcross),2.0)/Psi_Lmone/as<vec>(lambdalist[L-1]));
//     tau(L-1) = 1.0/R::rgamma(0.5*(npar(L-1)+1), 1.0/(1.0/xitau(L-1)+0.5*as_scalar(levelsums)));
//     xitau(L-1) = 1.0/R::rgamma(1,1.0/(1+1.0/tau(L-1)));
//     
//     // tau (own)//
//     levelsums = sum(pow(xi(wchown),2.0)/as<vec>(lambdalist_own[L-1]));
//     tau_own(L-1) = 1.0/R::rgamma(0.5*(npar_own(L-1)+1), 1.0/(1.0/xitau_own(L-1)+0.5*as_scalar(levelsums)));
//     xitau_own(L-1) = 1.0/R::rgamma(1,1.0/(1+1.0/tau_own(L-1)));
//     
//     // lambda (default to ridge)
//     lambda = ones(npar(L-1));
//     xilambda = as<vec>(xilambdalist[L-1]);
//     lambda_own = ones(npar_own(L-1));
//     xilambda_own = as<vec>(xilambdalist_own[L-1]);
//     
//     // lasso: Exp(1/2)=Gamma(1,1/2)
//     if(product_shrinkage=="lasso" && rep>=initial_run){
//       // cross
//       mutilde = sqrt(2.0)*sqrt(Psi_Lmone)*sqrt(tau(L-1))/abs(xi(wchcross));
//       lambda = 1/randinvgaussian(npar(L-1),mutilde,2);
//       lambdalist[L-1] = lambda;
//       
//       // own
//       mutilde = sqrt(2.0)*sqrt(tau_own(L-1))/abs(xi(wchown));
//       lambda_own = 1/randinvgaussian(npar_own(L-1),mutilde,2);
//       lambdalist_own[L-1] = lambda_own;
//     }
//     
//     // horseshoe: C+(0,1)
//     if(product_shrinkage=="horseshoe" && rep>=initial_run){
//       
//       // cross
//       scale = 1.0/xilambda + 0.5*pow(xi(wchcross),2.0)/Psi_Lmone/as_scalar(tau(L-1));
//       lambda = randig(npar(L-1),ones(npar(L-1)),scale);
//       scale = 1 + 1.0/lambda;
//       xilambda = randig(npar(L-1),ones(npar(L-1)),scale);
//       lambdalist[L-1] = lambda;
//       xilambdalist[L-1] = xilambda;
//       
//       // own
//       scale = 1.0/xilambda_own + 0.5*pow(xi(wchown),2.0)/as_scalar(tau_own(L-1));
//       lambda_own = randig(npar_own(L-1),ones(npar_own(L-1)),scale);
//       scale = 1 + 1.0/lambda_own;
//       xilambda = randig(npar_own(L-1),ones(npar_own(L-1)),scale);
//       lambdalist_own[L-1] = lambda_own;
//       xilambdalist_own[L-1] = xilambda_own;
//       
//     }
//     
//     // -------------------------------------------------------------- //
//     // print time and store draws
//     // -------------------------------------------------------------- //
//     
//     if(print){
//       // time
//       if(rep==0){
//         if(timer.toc()/60.0*(RepRun-rep-1)/(rep+1)<1){
//           Rprintf("Estimated time: %.1f seconds \n",timer.toc()*(RepRun-rep-1)/(rep+1));
//           Rprintf("  0%%   10   20   30   40   50   60   70   80   90   100%%\n");
//           Rprintf("  |----|----|----|----|----|----|----|----|----|----|\n");
//         }
//         else{
//           Rprintf("Estimated time: %.1f minutes \n",timer.toc()/60.0*(RepRun-rep-1)/(rep+1));
//           Rprintf("  0%%   10   20   30   40   50   60   70   80   90   100%%\n");
//           Rprintf("  |----|----|----|----|----|----|----|----|----|----|\n");
//         }
//       }
//       if(RepRun>50){
//         if ((rep+1)%(RepRun/50)==0){
//           if(rep+1==RepRun/50){
//             Rprintf("  *");
//           }
//           else if(rep<RepRun-1){
//             Rprintf("*");
//           }
//           else Rprintf("*\n");
//         }
//       }
//     }
//     
//     // store draws
//     if(rep>=initial_run && (rep+1)%keep==0){
//       mkeep = (rep-initial_run+1)/keep;
//       betadraws(mkeep-1,span::all) = trans(beta);
//       thetadraws(mkeep-1,span::all) = trans(unlist(thetalist,sum(npar(span(0,L-2))),L-1));
//       thetaowndraws(mkeep-1,span::all) = trans(unlist(thetalist_own,sum(npar_own(span(0,L-2))),L-1));
//       lambdadraws(mkeep-1,span::all) = trans(unlist(lambdalist,sum(npar),lambdalist.size()));
//       lambdaowndraws(mkeep-1,span::all) = trans(unlist(lambdalist_own,sum(npar_own),lambdalist_own.size()));
//       taudraws(mkeep-1,span::all) = trans(tau);
//       tauowndraws(mkeep-1,span::all) = trans(tau_own);
//       sigmasqdraws(mkeep-1,span::all) = trans(sigmasq);
//       phidraws(mkeep-1,span::all) = trans(phi);
//     }
//     
//   }
//   
//   // print total time elapsed
//   if(print){
//     if(timer.toc()/60.0<1){
//       Rprintf("Total Time Elapsed: %.1f seconds \n",timer.toc());
//     }
//     else{
//       Rprintf("Total Time Elapsed: %.1f minutes \n",timer.toc()/60.0);
//     }
//   }
//   
//   return List::create(
//     Named("betadraws") = betadraws,
//     Named("thetadraws") = thetadraws,
//     Named("thetaowndraws") = thetaowndraws,
//     Named("lambdadraws") = lambdadraws,
//     Named("lambdaowndraws") = lambdaowndraws,
//     Named("taudraws") = taudraws,
//     Named("tauowndraws") = tauowndraws,
//     Named("sigmasqdraws") = sigmasqdraws,
//     Named("phidraws") = phidraws
//   );
//   
// }
