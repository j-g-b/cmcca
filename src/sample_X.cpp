// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppDist)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <RcppDist.h>
#include <vector>
#include <random>
#include <typeinfo>
#include <algorithm>
#include <iostream>
#include <array>        // std::array
#include <chrono>
#include <math.h>
#include <set>
#include <ctime>
#include <fstream>

//
using namespace arma;
using namespace Rcpp;

//
Eigen::MatrixXd RcppToEigenMat(Rcpp::NumericMatrix & RMat){
  Eigen::MatrixXd EigMat = Eigen::MatrixXd::Zero(RMat.nrow(), RMat.ncol());
  for(int i = 0; i < RMat.nrow(); i = i + 1) {
    for(int j = 0; j < RMat.ncol(); j = j + 1){
      EigMat(i, j) = RMat(i, j);
    }
  }
  return(EigMat);
}

//
Eigen::VectorXd RcppToEigenVec(Rcpp::NumericVector & RVec){
  Eigen::VectorXd EigVec = Eigen::VectorXd::Zero(RVec.size());
  for(int i = 0; i < RVec.size(); i = i + 1) {
    EigVec[i] = RVec[i];
  }
  return(EigVec);
}

//
double AcceptanceRatioNorm(double x_new, Eigen::MatrixXd & X, Eigen::MatrixXd & A, Eigen::MatrixXd & B, int i, int k){
  double log_r;
  double r;
  log_r = (-0.5)*((1 + A(k, k))*(std::pow(x_new, 2) - std::pow(X(i, k), 2)) + 2*(x_new - X(i, k))*(A.row(k).dot(X.row(i)) - A(k, k)*X(i, k) + B(i, k))) + R::dnorm(X(i, k), 0, 1, TRUE) - R::dnorm(x_new, 0, 1, TRUE);
  r = std::exp(log_r);
  return(r);
}

//
double ProposeEpsHitRun(Eigen::MatrixXd & X, Eigen::MatrixXd & A, Eigen::MatrixXd & B, double i, Eigen::VectorXd & a, double & eps_lwr, double & eps_upr){
  double m;
  double v;
  double eps;
  v = 1.0/(1 + a.dot(A*a));
  m = -v*(a.dot(X.row(i).transpose() + A*X.row(i).transpose() + B.row(i).transpose()));
  eps = r_truncnorm(m, std::sqrt(v), eps_lwr, eps_upr);
  return(eps);
}

//
double AcceptanceRatioMVNorm(Eigen::VectorXd & x_new, Eigen::MatrixXd & X, Eigen::MatrixXd & A, Eigen::MatrixXd & B, int i, double eps){
  double log_r;
  double r;
  log_r = (-0.5)*(x_new.dot(x_new) - X.row(i).transpose().dot(X.row(i).transpose()) + x_new.dot(A*x_new) - X.row(i).transpose().dot(A*X.row(i).transpose()) + 2*(x_new - X.row(i).transpose()).dot(B.row(i).transpose())) - R::dnorm(eps, 0, 1, TRUE) + R::dnorm(0, 0, 1, TRUE);
  r = std::exp(log_r);
  return(r);
}

//
// [[Rcpp::export]]
Rcpp::List eigen_sample_X(Rcpp::NumericMatrix & YY, Rcpp::NumericMatrix & XX, Rcpp::NumericVector & WW,
                                   Rcpp::NumericMatrix & SS, Rcpp::NumericMatrix & AA, Rcpp::NumericMatrix & BB,
                                   Rcpp::IntegerVector & Sub)
{
  //
  int N = YY.nrow();
  int P = YY.ncol();
  int M = Sub.length();

  //
  // Define array containers for objects passed from R
  Eigen::MatrixXd Y = RcppToEigenMat(YY);
  Eigen::MatrixXd X = RcppToEigenMat(XX);
  Eigen::VectorXd W = RcppToEigenVec(WW);
  Eigen::MatrixXd S = RcppToEigenMat(SS);
  Eigen::MatrixXd A = RcppToEigenMat(AA);
  Eigen::MatrixXd B = RcppToEigenMat(BB);

  //
  // Define other arrays and variables to be used
  Eigen::VectorXd D = Eigen::VectorXd::Zero(N);
  Eigen::MatrixXd SigmaInv;
  Eigen::MatrixXd Sigma;
  Eigen::MatrixXd SigmaInvRev;
  Eigen::MatrixXd L;
  Eigen::VectorXd z = Eigen::VectorXd::Zero(P);
  Eigen::RowVectorXd Xprop;
  double eprob = 0.95;
  double q = R::qgamma(eprob, P / 2.0, 2.0, 1, 0);
  double log_r;
  double u;
  double obj = -1;
  double obj_new = -1;
  int count;
  bool is_cm;
  int i;

  //
  // Metropolis-Hastings
  for(int m = 0; m < M; m = m + 1){
    //
    i = Sub[m];
    //
    // Specify covariance matrix
    for(int j = 0; j < N; j = j + 1){
      if(j == i){
        D[j] = 1;
      } else {
        D[j] = 1.0 / ((X.row(i).transpose() - X.row(j).transpose()).dot(Y.row(i).transpose() - Y.row(j).transpose()));
        D[j] = std::pow(D[j], 2);
      }
    }
    SigmaInv = q*(Y.rowwise() - Y.row(i)).transpose()*D.asDiagonal()*(Y.rowwise() - Y.row(i));
    Sigma = SigmaInv.llt().solve(Eigen::MatrixXd::Identity(P, P));
    L = Sigma.llt().matrixL();

    //
    // Propose new Xi
    for(int k = 0; k < P; k = k + 1){
      z[k] = R::rnorm(0, 1);
    }
    Xprop = X.row(i) + (L*z).transpose();

    //
    // Compute covariance matrix for reverse proposal
    for(int j = 0; j < N; j = j + 1){
      if(j == i){
        D[j] = 1;
      } else {
        D[j] = 1.0 / ((Xprop.transpose() - X.row(j).transpose()).dot(Y.row(i).transpose() - Y.row(j).transpose()));
        D[j] = std::pow(D[j], 2);
      }
    }
    SigmaInvRev = q*(Y.rowwise() - Y.row(i)).transpose()*D.asDiagonal()*(Y.rowwise() - Y.row(i));

    //
    // Compute the log acceptance ratio
    log_r = -0.5*Xprop*A*Xprop.transpose();
    log_r = log_r + 0.5*X.row(i)*A*(X.row(i).transpose());
    log_r = log_r - (Xprop - X.row(i)).transpose().dot(B.row(i).transpose());
    log_r = log_r - 0.5*Xprop.squaredNorm() + 0.5*X.row(i).squaredNorm();
    log_r = log_r - 0.5*(X.row(i) - Xprop)*SigmaInvRev*(X.row(i) - Xprop).transpose();
    log_r = log_r + 0.5*std::log(SigmaInvRev.determinant());
    log_r = log_r + 0.5*(Xprop - X.row(i))*SigmaInv*(Xprop - X.row(i)).transpose();
    log_r = log_r - 0.5*std::log(SigmaInv.determinant());

    //
    // Check if new X is cm with Y
    count = 0;
    is_cm = true;
    for(int j = 0; j < N; j = j + 1){
      S(i, j) = S(i, j) + (Xprop - X.row(i)).transpose().dot((Y.row(i) - Y.row(j)).transpose());
    }
    obj = ((S.colwise() - W).rowwise() + W.transpose()).minCoeff();
    while(obj < -1e-12){
      W = ((-1.0*S).colwise() + W).colwise().maxCoeff();
      obj_new = ((S.colwise() - W).rowwise() + W.transpose()).minCoeff();
      if(std::abs(obj - obj_new) < 1e-12){
        count = count + 1;
      } else {
        count = 0;
      }
      if(count > 10 && obj_new < -1e-12){
        is_cm = false;
        break;
      }
      obj = obj_new;
    }

    // Accept / reject
    u = R::runif(0, 1);
    if(u < std::exp(log_r) && is_cm){
      X.row(i) = Xprop;
    } else {
      for(int j = 0; j < N; j = j + 1){
        S(i, j) = S(i, j) - (Xprop - X.row(i)).transpose().dot((Y.row(i) - Y.row(j)).transpose());
      }
    }
    //
    // Reset
    obj = -1;
  }
  //
  Rcpp::NumericMatrix rX(X.rows(), X.cols());
  //
  for(int i = 0; i < X.rows(); i = i + 1){
    for(int j = 0; j < X.cols(); j = j + 1){
      rX(i, j) = X(i, j);
    }
  }
  //
  Rcpp::NumericMatrix rS(S.rows(), S.cols());
  //
  for(int i = 0; i < S.rows(); i = i + 1){
    for(int j = 0; j < S.cols(); j = j + 1){
      rS(i, j) = S(i, j);
    }
  }
  //
  Rcpp::NumericVector rW(W.size());
  //
  for(int i = 0; i < W.size(); i = i + 1){
    rW[i] = W[i];
  }
  //
  return Rcpp::List::create(Rcpp::Named("X") = rX,
                            Rcpp::Named("S") = rS,
                            Rcpp::Named("W") = rW);
  //
}

