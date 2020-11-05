#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// HELPER FUNCTIONS ////////////////////////////////////////////////////////////

void vprint(arma::vec x, std::string msg){
  Rcout << msg << ": " << x << std::endl;
}
void mprint(arma::mat x, std::string msg){
  Rcout << msg << ": " << x << std::endl;
}
void iprint(arma::uword x, std::string msg){
  Rcout << msg << ": " << x << std::endl;
}
void dprint(double x, std::string msg){
  Rcout << msg << ": " << x << std::endl;
}
void nprint(std::string msg){
  Rcout << msg << std::endl;
}

arma::vec rep4k(arma::vec u){
  int n = u.n_elem;
  double tmp3;
  double ui;
  arma::vec out(n);
  for(int i=0; i<n; i++){
    ui = u(i);
    if(ui > 1.0){
      out(i) = 0.0;
    } else if(ui >= 0.0){
      tmp3 = 1.0 - ui;
      out(i) = tmp3*tmp3*tmp3;
    } else {
      out(i) = 0.0;
    }
  }
  return out;
}

arma::vec rep3k(arma::vec u){
  int n = u.n_elem;
  double ui;
  arma::vec out(n);
  for(int i=0; i<n; i++){
    ui = u(i);
    if(ui > 2.0){
      out(i) = 0.0;
    } else if(ui > 1.0){
      out(i) = 2.0 - 3.0*ui + 3.0/2.0*ui*ui - 1.0/4.0*ui*ui*ui;
    } else if(ui > 0.0){
      out(i) = 3.0*ui - 9.0/2.0*ui*ui + 7.0/4.0*ui*ui*ui;
    } else {
      out(i) = 0.0;
    }
  }
  return out;
}

arma::vec rep2k(arma::vec u){
  int n = u.n_elem;
  double ui;
  arma::vec out(n);
  for(int i=0; i<n; i++){
    ui = u(i);
    if(ui > 3.0){
      out(i) = 0.0;
    } else if(ui > 2.0){
      out(i) = 9.0/2.0 - 9.0/2.0*ui + 3.0/2.0*ui*ui - 1.0/6.0*ui*ui*ui;
    } else if(ui > 1.0){
      out(i) = -3.0/2.0 + 9.0/2.0*ui - 3.0*ui*ui + 7.0/12.0*ui*ui*ui;
    } else if(ui > 0.0){
      out(i) = 3.0/2.0*ui*ui - 11.0/12.0*ui*ui*ui;
    } else {
      out(i) = 0.0;
    }
  }
  return out;
}

arma::vec evenk(arma::vec u){
  int n = u.n_elem;
  double ui;
  arma::vec out(n);
  for(int i=0; i<n; i++){
    ui = u(i);
    if(ui > 4.0){
      out(i) = 0.0;
    } else if(ui > 3.0){
      out(i) = -1.0/6.0 * (ui-4)*(ui-4)*(ui-4);
    } else if(ui > 2.0){
      out(i) = -22.0/3.0 + 10.0*ui - 4.0*ui*ui + 1.0/2.0*ui*ui*ui;
    } else if(ui > 1.0){
      out(i) = 2.0/3.0 - 2.0*ui + 2.0*ui*ui - 1.0/2.0*ui*ui*ui;
    } else if(ui > 0.0){
      out(i) = 1.0/6.0*ui*ui*ui;
    } else {
      out(i) = 0.0;
    }
  }
  return out;
}

// [[Rcpp::export(".bs_even_C")]]
arma::mat bs_even(arma::vec time, int nk){
  arma::vec u;
  double dnk = nk - 1.0;
  u = dnk / time.max() * (time - time.min());
  int nbasis = nk + 2;
  int ni = nk - 2;
  int neven_basis = ni - 2;
  arma::mat out(time.n_elem, nbasis);
  out.col(0) = rep4k(u);
  out.col(1) = rep3k(u);
  out.col(2) = rep2k(u);
  arma::vec revu = 3 + neven_basis - u;
  out.col(nbasis - 1) = rep4k(revu);
  out.col(nbasis - 2) = rep3k(revu);
  out.col(nbasis - 3) = rep2k(revu);
  for(int i=0; i < neven_basis; i++){
    out.col(3 + i) = evenk(u - i);
  }
  return out;
}

arma::vec rmnorm(arma::vec mu, arma::mat covprec, bool is_prec) {//slower but more stable
  arma::vec z = rnorm(mu.n_elem);
  arma::vec evals;
  arma::mat evecs;
  if (is_prec) {
    return mu + solve(chol(covprec), z);
  } else {
    return mu + chol(covprec).t() * z;
    //arma::vec evals;
    //arma::mat evecs;
    //arma::eig_sym(evals, evecs, covprec);
    //arma::mat Dhalf = diagmat(sqrt(arma::clamp(evals, 0.0, evals.max())));
    //arma::mat Dhalf = diagmat(sqrt(evals));
    //return mu + evecs*Dhalf*evecs.t()*z;
  }
}

arma::mat penalty_mat(int p){
  arma::mat P(p,p);
  P.zeros();
  P(0,0) = 2;
  for(int i = 1; i<p-1; i++){
    P(i,i) = 2;
    P(i,i-1)=-1;
    P(i-1,i) =-1;
  }
  P(p-1, p-1) = 1;
  return P;
}


// Package Exports /////////////////////////////////////////////////////////////

// [[Rcpp::export(".g2g_fda")]]
List g2g_fda(arma::mat curves, arma::vec lasts, arma::vec time,
                  int p, int niter, int nburn){
  // SETUP //-------------------------------------------------------------------
  // Constants
  int ngrp = lasts.n_elem;
  int ncrv = curves.n_cols;
  int npts = curves.n_rows;
  arma::mat H = bs_even(time, p-2);
  arma::mat HtH = H.t()*H;
  arma::mat Hty = H.t()*curves;
  arma::vec crv_grp(ncrv);
  int counter = 0;
  for(int i=0; i<ncrv; i++){
    crv_grp(i) = counter;
    if(i == lasts(counter)){
      counter++;
    }
  }
  arma::mat P = penalty_mat(p);
  arma::mat Pi = arma::inv_sympd(P);

  // Parameters
  double sig2 = 1;
  double tau2 = 1;
  arma::mat beta_grp(p, ngrp);
  beta_grp.zeros();
  arma::mat beta_crv(p, ncrv);
  beta_crv.zeros();

  // Save Structure
  arma::vec sig2_chain(niter);
  sig2_chain.zeros();
  arma::vec tau2_chain(niter);
  tau2_chain.zeros();
  arma::cube beta_grp_chain(p, ngrp, niter);
  beta_grp_chain.zeros();
  arma::cube beta_crv_chain(p, ncrv, niter);
  beta_crv_chain.zeros();

  // Tmp variables
  arma::mat VVi;
  double sig2_sse;
  double tau2_sse;
  arma::vec bbdiff;

  // MCMC //--------------------------------------------------------------------

  for(int it = 0; it < nburn+niter; it++){

    // Update beta_crv
    VVi = arma::inv_sympd(1/tau2*P + 1/sig2*HtH);
    for(int i = 0; i < ncrv; i++){
      beta_crv.col(i) = rmnorm(VVi*(1/tau2*P*beta_grp.col(crv_grp(i)) + 1/sig2*Hty.col(i)), VVi, false);
    }

    // Update beta_grp
    for(int j = 0; j < ngrp; j++){
      if(j == 0){
        beta_grp.col(j) = rmnorm(arma::mean(beta_crv.cols(0, lasts(j)),1), 1/lasts(j)*tau2*Pi, false);
      } else {
        beta_grp.col(j) = rmnorm(arma::mean(beta_crv.cols(lasts(j-1)+1, lasts(j)), 1), 1/(lasts(j)-lasts(j-1))*tau2*Pi, false);
      }
    }

    // Update sig2
    sig2_sse = arma::accu(arma::square(curves - H*beta_crv));
    sig2 = 1/rgamma(1, 0.5*ncrv*npts, 1/(0.5*sig2_sse))(0);

    //Update tau2
    tau2_sse = 0;
    for(int i = 0; i < ncrv; i++){
      bbdiff = beta_crv.col(i) - beta_grp.col(crv_grp(i));
      tau2_sse += arma::as_scalar(bbdiff.t()*P*bbdiff);
    }
    tau2 = 1/rgamma(1, 0.5*ncrv*p, 1/(0.5*tau2_sse))(0);

    // Saves
    if(it >= nburn){
      Rcout << "\r Iteration: " << it-nburn+1 << " / " << niter;
      tau2_chain(it-nburn) = tau2;
      sig2_chain(it-nburn) = sig2;
      beta_grp_chain.slice(it-nburn) = beta_grp;
      beta_crv_chain.slice(it-nburn) = beta_crv;
    }

  }
  List chains = List::create(tau2_chain, sig2_chain, beta_grp_chain, beta_crv_chain);
  return chains;
}
