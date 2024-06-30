
// #define ARMA_64BIT_WORD 1
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>

// #define INT_MIN (-INT_MAX - 1)

using namespace Rcpp;
using namespace arma;
using namespace std;


/*
 * Auxiliary
 */
//' @keywords internal
//' @noRd
//' 
//' Data X has been centered.

//' Evaluate the log determinant of the covariance matrix and energy function values on each sample
//' point, by SVD decomposition to achieve efficeint computation.
//' @param X: a n x p  matrix, such as gene expression logcount matrix
//' @param Lam_vec0: a column vector, variance of each error component in factor model
//' @param W0: a p x q matrix, loading matrix
//' @param Ck: a q x q matrix,
//' @param Muk: a q x K matrix, gaussian mixture mean component
//' @param logdSk: a real number, the returned log determinant of covariance matrix
//' @param mSk: a colunmn vector, the returned energy funciton values. 
//' 
//' @return No return
//' 


// [[Rcpp::export]]
List multi_det_SkCpp2(const arma::mat& X, const arma::vec& Lam_vec0,
                      const arma::mat& W0, const arma::mat& Ck, 
                      const arma::rowvec Muk, const arma::mat& Sigmak) {
  int n = X.n_rows, q = W0.n_cols;
  
  arma::mat WC12, tmp2;
  arma::vec tmp1, s, tmp3;
  arma::mat U, V, X_tk;
  
  // Compute the SVD of Sigmak and calculate logdSk
  svd(U, s, V, Sigmak);
  WC12 = W0 * (U.each_row() % trans(sqrt(s))); 
  WC12 = WC12.each_col() % (1.0 / sqrt(Lam_vec0)); 
  arma::vec d = svd(WC12);
  double logdSk = -accu(log(1 + d % d)) - accu(log(Lam_vec0));
  
  // Compute the SVD of the inverse of Ck and calculate mSk
  //std::cout << "Ck:\n" << Ck.i() << std::endl;
  svd(U, s, V, Ck.i());
  WC12 = W0 * (U.each_row() % trans(sqrt(s))); 
  WC12 = WC12.each_col() % (1.0 / sqrt(Lam_vec0)); 
  X_tk = (X - repmat(Muk * W0.t(), n, 1)) % trans(repmat(1.0 / sqrt(Lam_vec0), 1, n)); 
  tmp1 = sum(X_tk % X_tk, 1); 
  tmp2 = X_tk * WC12; 
  tmp3 = sum(tmp2 % tmp2, 1); 
  arma::vec mSk = tmp1 - tmp3;
  
  // Convert Armadillo vectors to Rcpp vectors
  NumericVector mSkR = wrap(mSk);
  
  // Create a list to store the results
  List result;
  result["logdSk"] = logdSk;
  result["mSk"] = mSkR;
  
  return result;
}


sp_mat get_spNbs(ivec y, const sp_mat& Adj) {   
  // row is for pixel.
  //output a sparse matrix, i-th row contains labels of neighbor_i. 
  // Make const iterator
  arma::sp_mat::const_iterator start = Adj.begin(); 
  //arma::sp_mat::const_iterator end   = Adj.end();
  // Calculate number of nonzero points
  //int n = std::distance(start, end);
  int n = Adj.n_nonzero; 
  //cout << "n=" << n << endl;
  //cout << "n=" << Adj.n_nonzero << endl;
  
  sp_mat spNbs(y.n_elem, y.n_elem);    // neiborhood state matrix, matched with Adj.
  arma::sp_mat::const_iterator it = start; // Note spNbs is not a symmetric matrix, the nonzero in i-th row is the class label of sample i.
  for(int i = 0; i < n; ++i)
  {
    //temp(0) = it.row();
    //temp(1) = it.col();
    spNbs(it.row(), it.col()) = y(it.col()); 
    ++it; // increment
  }
  return spNbs.t(); // return the class label of neighbor matrix, i-th column is the neighbor label of sample i
}

// [[Rcpp::export]]
arma::mat calYenergy2D_sp(const arma::ivec& y, const arma::sp_mat& Adj, int K, const arma::vec alpha, const double beta)	{
  // Calculate the value of energy function of y, which is equal to negative logliklihood up to a constant
  int n = y.n_rows;
  arma::sp_mat spNbs_t = get_spNbs(y, Adj); // transform spNbs to iterate by column.
  arma::mat Uy(n, K);
  double n_sameS;
  int i, k, nn;
  for (k = 0; k < K; k++){
    for (i = 0; i < n; i++){
      arma::sp_mat col(spNbs_t.col(i)); // the class label of neighbors of i-th sample.
      n_sameS = 0;
      nn = col.n_nonzero; // the number of neighbors of i-th sample
      for (arma::sp_mat::iterator j = col.begin(); j != col.end(); ++j) {
        n_sameS += ((*j) == (k+1));
      }
      Uy(i, k) = alpha(k) + beta * (nn - n_sameS)/2;
   }
  }
  arma::mat C_mat = normalise(exp(-Uy), 1, 1); // pseudo likelihood of Y.
  Uy = -log(C_mat); // normalized Uy, this is the energy of y.
  return Uy;
  }


double obj_beta(const arma::ivec& y, const arma::mat& R, const arma::sp_mat& Adj, int K, const arma::vec alpha, const double beta)	{
  
  mat Uy = calYenergy2D_sp(y, Adj, K, alpha, beta); // Uy was normalized, so there is no need to normalized Uy. 
  //arma::mat C_mat = normalise(exp(-Uy), 1, 1); // set all rowSums to be ONE to get the likelihood
  //return accu(R % log(C_mat)); 
  return -accu(R % Uy);
}


// [[Rcpp::export]]
List runICM_sp (const arma::mat& X, const arma::mat& Ux,  arma::ivec& y, const arma::mat& Mu0,
                const arma::cube& Sigma0,const arma::sp_mat& Adj, const arma::vec& alpha, const arma::vec& beta_grid,
                double& beta, int maxIter_ICM)	{
  // Target: estimate Y, evaluate R, and update beta by using grid search.
  // basic info.
  int n = X.n_rows, K = Mu0.n_rows, q= Mu0.n_cols;
  int iter, k;
  
  // Estimate Y by ICM
  arma::vec Energy(maxIter_ICM);
  Energy(0) = INFINITY;
  arma::mat Uy(n, K);
  arma::mat U(n, K);
  arma::mat R(n, K);
  //--------------------------------------------------------------------------------	
  // ICM algrithm to estimate Y
  //--------------------------------------------------------------------------------
  // int Iteration = 1;
  for (iter = 1; iter < maxIter_ICM; iter ++ ) {
    Uy = calYenergy2D_sp(y, Adj, K, alpha, beta);
    U = Uy + Ux; // log likelihood of (x, y).
    arma::vec Umin = min(U, 1);
    arma::uvec y_u = index_min(U, 1);
    y = conv_to< ivec >::from(y_u) + 1;
    
    Energy(iter) = sum(Umin);
    if (Energy(iter) - Energy(iter - 1) > 1e-5) {
      //cout << "diff Energy = " << Energy(iter) - Energy(iter - 1)  << endl;
      break;
    }
    if (Energy(iter-1) - Energy(iter) < 1e-5){
      //cout << "ICM Converged at Iteration = " << iter  << endl;
      break;
    }
 }
  // calculate R and pseudo observed loglikelihood
  vec maxA1 = max(-U, 1);
  U = (-U - repmat(maxA1, 1, K));
  vec loglik_more_vec = sum(exp(U),1);
  double loglik;
  loglik = sum(log(loglik_more_vec) + maxA1); //  - n* p /2.0 * log(2* M_PI); 
  R = exp(U) / repmat(loglik_more_vec, 1, K);
  
  // vec energy = Energy.subvec(1, Iteration);
  // update beta: grid search.
  int ng_beta = beta_grid.n_elem;
  vec objBetaVec(ng_beta);
  for(k=0; k < ng_beta; ++k){
    objBetaVec(k) = obj_beta(y, R, Adj, K, alpha, beta_grid(k));
  }
  beta = beta_grid(index_max(objBetaVec));
  
  List output = List::create(
    Rcpp::Named("y") = y,
    Rcpp::Named("R") = R,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("hbeta") = beta_grid(index_max(objBetaVec))); // get the maximum beta.
  
  return output; 
}  



// [[Rcpp::export]]
List runICM_sp2 (const arma::mat& X, const arma::mat& W0, const arma::vec& Lam_vec0, const arma::mat& Mu0,
                const arma::cube& Sigma0,const arma::mat& Iden)	{
  // Target: estimate Y, evaluate R, Ez, Ck inverse, and update beta by using grid search.
  // basic info.
  int n = X.n_rows, K = Mu0.n_rows, q= Mu0.n_cols;
  int iter, k;
  
  // two cached objects used for parameters update.
  cube Cki_ara(q, q, K, fill::zeros), Ez(n,q,K, fill::zeros);
  double  logdSk;
  vec mSk(n);
  mat WtLW = W0.t() * (repmat(1.0/ Lam_vec0, 1, q) %  W0); //O(p^2 q) cache object 
  mat XLW = X * (repmat(1.0/ Lam_vec0, 1, q) % W0);
  // evaluate energy of x, Ux
  arma::mat Ux(n, K), Ck, inv_S;
  
  for (k = 0; k < K; k++)	{
    inv_S=inv_sympd(Sigma0.slice(k),inv_opts::allow_approx);
    //std::cout << "a:\n" << inv_S << std::endl;
    Ck = WtLW + inv_S;
    Cki_ara.slice(k) = inv_sympd(Ck,inv_opts::allow_approx);
    //std::cout << "Ck:\n" << Cki_ara.slice(k) << std::endl;
    multi_det_SkCpp2(X, Lam_vec0,W0, Ck, Mu0.row(k), Sigma0.slice(k));
    Ux.col(k) = -0.5*logdSk  + 0.5 * mSk; // calculate energy by column.
    Ez.slice(k) = (XLW + repmat(Mu0.row(k)* inv_S, n, 1)) * Ck.i();
  }
  
  // calculate R and pseudo observed loglikelihood
  vec maxA1 = max(-Ux, 1);
  Ux = (-Ux - repmat(maxA1, 1, K));
  vec loglik_more_vec = sum(exp(Ux),1);
  double loglik;
  loglik = sum(log(loglik_more_vec) + maxA1); //  - n* p /2.0 * log(2* M_PI); 
  
  List output = List::create(
    Rcpp::Named("Ez") = Ez,
    Rcpp::Named("Cki_ara") = Cki_ara,
    Rcpp::Named("loglik") = loglik);
  
  return output; 
} 



//' Evaluate the diagonal elements of three matrix multiplication such as W0*Cki*W0^T
//'  by SVD decomposition to achieve efficeint computation.
//' @param Cki: a q x q matrix,
//' @param W0: a p x q matrix
//' 
//' @return a column vector
//'   

vec decomp(const mat& Cki, const mat& W0){
  vec s, tmp1;
  mat U, V, WC12;
  // std::cout << "ak:\n" << Cki << std::endl;
  svd(U, s, V, Cki);
  WC12 = W0 * (U.each_row() % trans(sqrt(s)));
  tmp1 = sum(WC12 % WC12, 1);
  return tmp1;
}  


// [[Rcpp::export]] 
mat calculateX(const mat& R, const cube& Ez){
  int k, K= R.n_cols, q= Ez.n_cols, n= R.n_rows;
  mat Ezz(n, q, fill::zeros); // estimate Z, factor matrix
  for(k=0; k<K; k++){
    Ezz +=  Ez.slice(k) % repmat(R.col(k), 1, q);
  }
  return Ezz;
}



//' Evaluate the diagonal elements of three matrix multiplication such as W0*Cki*W0^T
//'  by SVD decomposition to achieve efficeint computation.
//' @param Cki: a q x q matrix,
//' @param W0: a p x q matrix
//' 
//' @return a column vector
//'  
// [[Rcpp::export]] 
mat update_W0(const mat& X, const mat& R, const cube& Ez, const cube& Ci,  const vec& N){
  int k, K= R.n_cols, q= Ez.n_cols, n= R.n_rows;
  mat tmpMat(n,q, fill::zeros), Ezzt(q,q, fill::zeros), tmpMat2;
  // vec N = sum(R.t(), 1);
  for(k=0; k<K; ++k){
    tmpMat2= repmat(R.col(k), 1, q) % Ez.slice(k);
    tmpMat += tmpMat2;
    Ezzt+= tmpMat2.t() * Ez.slice(k) + N(k) * Ci.slice(k);
  }
  return X.t() * tmpMat * Ezzt.i();
}

// update Sigma0
// [[Rcpp::export]]
cube update_Sigma0(const mat& R, const cube& Ez, const cube& Ci, const mat& Mu,  
                   const vec&  N){
  int k, K= R.n_cols, q= Mu.n_cols, n= R.n_rows;
  cube Sigma0(q,q,K);
  for(k = 0; k<K; ++k){
    Sigma0.slice(k) = (trans(Ez.slice(k) - repmat(Mu.row(k), n, 1)) % trans(repmat(R.col(k), 1, q))) * (Ez.slice(k) - repmat(Mu.row(k), n, 1)) + N(k)*  Ci.slice(k);
    Sigma0.slice(k) = Sigma0.slice(k) / N(k);
    
  }
  return(Sigma0);
}

// update Lambda
// [[Rcpp::export]]
vec update_Lam(const mat& R, const mat& X, const mat& W, const cube& Ez, const cube& Ci){
  int k, K= R.n_cols, p = X.n_cols;
  rowvec N = sum(R);
  vec Lsum(p,fill::zeros), Lam;
  mat tmpXk;
  for(k=0; k<K; ++k){
    tmpXk = (X - Ez.slice(k) * W.t() );
    Lsum += trans(sum(tmpXk % tmpXk % repmat(R.col(k), 1, p)));
    Lsum += N(k) * decomp(Ci.slice(k), W);
  }
  Lam = Lsum/(X.n_rows*1.0);
  return Lam; 
}

// // update sigma20: homo variance
// vec update_Lam2(const mat& R, const mat& X, const mat& W, const cube& Ez, const cube& Ci){
//   int k, K= R.n_cols, p = X.n_cols;
//   mat tmpXk;
//   double term1=0, term2=0;
//   for(k=0; k <K; ++k){
//     tmpXk = (X - Ez.slice(k) * W.t() );
//     term1 += accu(sum(tmpXk % tmpXk, 1) % R.col(k));
//     term2  += accu(decomp(Ci.slice(k), W))* accu(R.col(k));
//   }
//   double sigma20 = (term1+ term2)/ X.n_elem;
//   return(ones(p)* sigma20);
// }

//Calculate Q function
double Q_fun(const mat& X, const mat& R,  const cube& Ez, const arma::cube& Ci, 
             const mat& W0, const mat& Mu0, const cube& Sigma0, const vec& Pi0, 
             const vec& Lam_vec0){
  
  double Q = 0, tmp_scalar =0;
  int  K = Pi0.n_elem, n = X.n_rows;
  int q = Mu0.n_cols;
  mat tmpSig(q, q, fill::zeros), tmpMat;
  colvec Ezik;
  mat Ezzt(q,q, fill::zeros);
  for(int k=0; k<K; k++){
    tmpSig = Sigma0.slice(k);
    tmpMat = Ez.slice(k);
    for(int i = 0; i<n; i++){
      Ezik =  trans(tmpMat.row(i)); 
      Ezzt = Ci.slice(k) + Ezik * Ezik.t();
      
      tmp_scalar =  0.5* accu(log(Lam_vec0)) + 0.5* accu( (X.row(i) % X.row(i)) / Lam_vec0.t()) +
        0.5 * arma::as_scalar(accu(W0.t()%trans(repmat(1.0/Lam_vec0, 1, q))*W0 % Ezzt.t())- 2* (X.row(i)% (1.0/Lam_vec0))*W0*Ezik);
      Q +=  - R(i,k) * tmp_scalar + R(i,k)*(log(Pi0(k)) -  0.5* log(det(tmpSig))- 0.5* accu(inv_sympd(symmatu(tmpSig))%
        Ezzt.t())+ arma::as_scalar(Mu0.row(k) * solve(tmpSig, (Ezik- 0.5*trans(Mu0.row(k)))))); 
    }
  }
  return Q;
}
