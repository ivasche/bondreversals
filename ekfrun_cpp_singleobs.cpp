// THIS IS AN ALTERNATIVE FORMULATION WITHOUT EXPLICIT IDENTIFICATION OF CREDIT RISK PREMIA Lambda
// kml = kappa minus lambda
// ktt = kappa times theta

#include <RcppArmadillo.h>
#include <boost/math/differentiation/finite_difference.hpp>  // included in BH  

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]    

// Functions to compute CIR2-implied yields

// A
// [[Rcpp::export]]
double A_tT_cpp(double mat, double kml, double ktt, double sigma) {
  
  double h = sqrt(kml * kml + 2.0 * sigma * sigma) ;
  return(pow((2.0*h*exp((kml+h)*mat/2.0))/(2.0*h+(kml+h)*(exp(mat*h)-1)), 2.0*ktt/(sigma*sigma)));
  
}

// B
// [[Rcpp::export]]
double B_tT_cpp(double mat, double kml, double ktt, double sigma) {
  
  double h = sqrt(kml * kml + 2.0 * sigma * sigma) ;
  return((2.0*(exp(mat*h)-1.0))/(2.0*h+(kml+h)*(exp(mat*h)-1.0))) ;
  
}

// vectorized A
// [[Rcpp::export]]
arma::vec vA_tT_cpp(arma::vec mat, double kml, double ktt, double sigma) {
  
  int n = mat.size();
  arma::vec z = arma::zeros(n, 1);
  double h = sqrt(kml * kml + 2.0 * sigma * sigma) ;
  double pwr = 2.0 * ktt/(sigma*sigma);
  
  for (int i=0; i<n; ++i) {
    z[i] = pow((2.0*h*exp((kml+h)*mat[i]/2.0))/(2.0*h+(kml+h)*(exp(mat[i]*h)-1.0)), pwr);
  }
  
  return(z) ;
}

// another vectorized A
// [[Rcpp::export]]
Rcpp::NumericVector vA_tT_cpp2(Rcpp::NumericVector mat, double kml, double ktt, double sigma) {
  
  Rcpp::NumericVector res = Rcpp::sapply(mat, [&](double x){return A_tT_cpp(x, kml, ktt, sigma);});
  return(res);
  
}

// vectorized B (returns a column vector = R matrix)
// [[Rcpp::export]]
arma::vec vB_tT_cpp(arma::vec mat, double kml, double ktt, double sigma) {
  
  int n = mat.size();
  arma::vec z = arma::zeros(n, 1);
  double h = sqrt(kml * kml + 2.0 * sigma * sigma) ;
  
  for (int i=0; i<n; ++i) {
    z[i] = (2.0*(exp(mat[i]*h)-1.0))/(2.0*h+(kml+h)*(exp(mat[i]*h)-1.0));
  }
  
  return(z) ;
}

// another vectorized B (returns a numerical vector)
// [[Rcpp::export]]
Rcpp::NumericVector vB_tT_cpp2(Rcpp::NumericVector mat, double kml, double ktt, double sigma) {
  
  Rcpp::NumericVector res = Rcpp::sapply(mat, [&](double x){return B_tT_cpp(x, kml, ktt, sigma);});
  return(res);
  
}

// price of a coupon bond
// [[Rcpp::export]]
double cb_price_cpp(double coup, arma::vec coup_dates, double coup_freq, double ytm, double notional) {
  
  int n = coup_dates.size();
  double z = 0;
  
  for (int i=0; i<n; ++i) {
    z += pow(1+ytm/coup_freq, -coup_freq*coup_dates[i])*coup/coup_freq*notional;
  }
  
  return(z + pow(1.0+ytm/coup_freq, -coup_freq*coup_dates.back())*notional);
  
}

// price of a coupon bond with coupons as a vector of (potentially distinct) cash flows
// [[Rcpp::export]]
double cb_price_full_cpp(arma::vec coup, arma::vec coup_dates, double coup_freq, double ytm, double notional) {
  
  int n = coup_dates.size();
  double z = 0;
  
  for (int i=0; i<n; ++i) {
    z += pow(1.0+ytm/coup_freq, -coup_freq*coup_dates[i])*coup[i]*notional;
  }
  
  return(z);
  
}

// price of a coupon bond with discount factors
// [[Rcpp::export]]
double cb_price_dfs_cpp(arma::vec coup, arma::vec dfs, double notional) {
  
  int n = coup.size();
  double z = 0;
  
  for (int i=0; i<n; ++i) {
    z += dfs[i]*coup[i]*notional;
  }
  
  return(z);
  
}

// price of a zero-coupon risk-free bond in the CIR2 model
// [[Rcpp::export]]
double PZC_CIR2_rf_cpp(double mat, double s1, double s2, double alpha, double kml1, double ktt1, 
                       double sigma1, double kml2, double ktt2, double sigma2) {
  
  return(exp(-alpha*mat + log(A_tT_cpp(mat, kml1, ktt1, sigma1)) + 
         log(A_tT_cpp(mat, kml2, ktt2, sigma2)) -
         s1*B_tT_cpp(mat, kml1, ktt1, sigma1) -
         s2*B_tT_cpp(mat, kml2, ktt2, sigma2))) ;
  
}

// price of a defaultable zero-coupon bond in Duffee (99) model
// [[Rcpp::export]]
double PZC_D99_def_cpp(double mat, double s1, double s2, double sh, double alpha, double kml1, double ktt1, 
                       double sigma1, double kml2, double ktt2, double sigma2, double alpha_def, 
                       double s1bar, double s2bar, double beta1, double beta2, double kmlh, double ktth,
                       double sigmah, double recov) {
  
  double V0 = exp(-(alpha+alpha_def-beta1*s1bar - beta2*s2bar)*mat + 
                  log(A_tT_cpp(mat, kml1, ktt1*(1.0+beta1), sigma1*sqrt(1.0+beta1))) +
                  log(A_tT_cpp(mat, kml2, ktt2*(1.0+beta2), sigma2*sqrt(1.0+beta2))) +
                  log(A_tT_cpp(mat, kmlh, ktth, sigmah))-
                  s1*B_tT_cpp(mat, kml1, ktt1*(1.0+beta1), sigma1*sqrt(1.0+beta1)) -
                  s2*B_tT_cpp(mat, kml2, ktt2*(1.0+beta2), sigma2*sqrt(1.0+beta2)) -
                  sh*B_tT_cpp(mat, kmlh, ktth, sigmah)) ;
  
  double P = PZC_CIR2_rf_cpp(mat, s1, s2, alpha, kml1, ktt1, sigma1, kml2, ktt2, sigma2) ;
  
  return(recov*P + (1.0-recov)*V0) ;
  
}

// vectorize price of a defaultable zero-coupon bond in Duffee (99) model
// [[Rcpp::export]]
Rcpp::NumericVector vPZC_D99_def_cpp(Rcpp::NumericVector mat, double s1, double s2, double sh, double alpha,
                                     double kml1, double ktt1, double sigma1,
                                     double kml2, double ktt2, double sigma2,
                                     double alpha_def, double s1bar, double s2bar, double beta1, double beta2, 
                                     double kmlh, double ktth, double sigmah, double recov) {
  
  Rcpp::NumericVector res = Rcpp::sapply(mat, [&](double x){
    
    return PZC_D99_def_cpp(x, s1, s2, sh, alpha, kml1, ktt1, sigma1, kml2, ktt2, sigma2, alpha_def, 
                           s1bar, s2bar, beta1, beta2, kmlh, ktth, sigmah, recov);
    
  });
  return(res);
  
}

// [A] pars_bonds: a list where each element (=bond) contains coupon, coupon frequency, coupon dates, notional and recovery
// 0. coup
// 1. coup_freq (not used!!)
// 2. c_dates
// 3. notional
// 4. recov

// [B] pars_s:
// 0. alpha
// 1. kml1
// 2. ktt1
// 3. sigma1
// 4. kml2
// 5. ktt2
// 6. sigma2
// 7. s1bar
// 8. s2bar

// [C] pars_h
// 0. alpha_def
// 1. beta1
// 2. beta2
// 3. kmlh
// 4. ktth
// 5. sigmah

// price of a single defaultable coupon bond in Duffee (99) model
// [[Rcpp::export]]

double Price_D99_def_single_cpp(double s_def, Rcpp::NumericVector s_1, Rcpp::NumericVector s_2,
                                Rcpp::List pars_s, Rcpp::List pars_h, Rcpp::List pars_bonds, int t_ind) {
  
  double coup = pars_bonds[0] ;
  double coup_freq = pars_bonds[1] ;
  Rcpp::List c_dates = pars_bonds[2] ;
  int ind = t_ind - 1;
  
  // order: vPZC_D99_def_cpp(mat, s1, s2, sh, alpha, kml1, ktt1, sigma1, kml2, ktt2, sigma2,
  //                         alpha_def, s1bar, s2bar, beta1, beta2, kmlh, ktth, sigmah, recov)
  Rcpp::NumericVector V = vPZC_D99_def_cpp(c_dates[ind],    // mat
                                           s_1[ind],        // s1
                                              s_2[ind],        // s2
                                                 s_def,           // sh
                                                 pars_s[0],       // alpha
                                                       pars_s[1],       // kml1
                                                             pars_s[2],       // ktt1
                                                                   pars_s[3],       // sigma1
                                                                         pars_s[4],       // kml2
                                                                               pars_s[5],       // ktt2
                                                                                     pars_s[6],       // sigma2
                                                                                           pars_h[0],       // alpha_def
                                                                                                 pars_s[7],       // s1bar
                                                                                                       pars_s[8],       // s2bar
                                                                                                             pars_h[1],       // beta1
                                                                                                                   pars_h[2],       // beta2
                                                                                                                         pars_h[3],       // kmlh
                                                                                                                               pars_h[4],       // ktth
                                                                                                                                     pars_h[5],       // sigmah
                                                                                                                                           pars_bonds[4]) ; // recov
  
  double P_coup = sum(V)*coup/coup_freq+V[V.length()-1] ;
  return(P_coup) ;
  
}

// price of a single defaultable coupon bond in Duffee (99) model
// [[Rcpp::export]]

double Price_D99_def_single_full_cpp(double s_def, Rcpp::NumericVector s_1, Rcpp::NumericVector s_2,
                                     Rcpp::List pars_s, Rcpp::List pars_h, Rcpp::List pars_bonds, int t_ind) {
  
  Rcpp::List coup = pars_bonds[0] ;
  Rcpp::List c_dates = pars_bonds[2] ;
  Rcpp::NumericVector aivec = pars_bonds[5] ;
  int ind = t_ind - 1;
  
  double ai = aivec[ind] ;
  
  Rcpp::NumericVector V = vPZC_D99_def_cpp(c_dates[ind],    // mat
                                           s_1[ind],        // s1
                                              s_2[ind],        // s2
                                                 s_def,           // sh
                                                 pars_s[0],       // alpha
                                                       pars_s[1],       // kml1
                                                             pars_s[2],       // ktt1
                                                                   pars_s[3],       // sigma1
                                                                         pars_s[4],       // kml2
                                                                               pars_s[5],       // ktt2
                                                                                     pars_s[6],       // sigma2
                                                                                           pars_h[0],       // alpha_def
                                                                                                 pars_s[7],       // s1bar
                                                                                                       pars_s[8],       // s2bar
                                                                                                             pars_h[1],       // beta1
                                                                                                                   pars_h[2],       // beta2
                                                                                                                         pars_h[3],       // kmlh
                                                                                                                               pars_h[4],       // ktth
                                                                                                                                     pars_h[5],       // sigmah
                                                                                                                                           pars_bonds[4]) ; // recov
  
  double P_coup = sum(Rcpp::as<arma::vec>(V) % Rcpp::as<arma::vec>(coup[ind])) - ai ;
  return(P_coup) ;
  
}


// vectorize price of a single defaultable coupon bond in Duffee (99) model
// [[Rcpp::export]]

Rcpp::NumericVector Price_D99_def_cpp(Rcpp::List pb, double s_def, Rcpp::NumericVector s_1, Rcpp::NumericVector s_2,
                                      Rcpp::List pars_s, Rcpp::List pars_h, int t_ind) {
  
  int n = pb.size() ;
  Rcpp::NumericVector res = Rcpp::NumericVector(n) ;
  
  for (int i=0; i<n; ++i) {
    //    res[i] = Price_D99_def_single_cpp(s_def, s_1, s_2, pars_s, pars_h, pb[i], t_ind) ;
    res[i] = Price_D99_def_single_full_cpp(s_def, s_1, s_2, pars_s, pars_h, pb[i], t_ind) ;
    
  }
  
  return(res);
  
}

// create a local function to differentiate Price_D99_def_cpp with the boost library finite_difference_derivative
// [[Rcpp::export]]

Rcpp::NumericVector Price_D99_def_linearized_cpp(Rcpp::List pb, double s_def, Rcpp::NumericVector s_1,
                                                 Rcpp::NumericVector s_2, Rcpp::List pars_s, Rcpp::List pars_h, int t_ind){
  
  int n = pb.size() ;
  Rcpp::NumericVector res = Rcpp::NumericVector(n) ;
  
  for (int i=0; i<n; ++i) {
    
    //    auto f = [&](double s_def) { return Price_D99_def_single_cpp(s_def, s_1, s_2, pars_s, pars_h, pb[i], t_ind) ; } ;
    auto f = [&](double s_def) { return Price_D99_def_single_full_cpp(s_def, s_1, s_2, pars_s, pars_h, pb[i], t_ind) ; } ;
    
    res[i] = boost::math::differentiation::finite_difference_derivative(f,s_def) ;
    
  }
  
  return(res);  
  
}

// A single run of an Extended Kalman Filter for Duffee (99)
// [[Rcpp::export]]

Rcpp::List EKF_run_cpp(arma::mat Rt, arma::cube a, arma::cube Phi, arma::mat U, arma::vec y0,
                       arma::cube w0, arma::cube w1, int Nobs, int Nfcs, int Nbonds, double dt,
                       Rcpp::List pars_b, Rcpp::List pars_s, Rcpp::List pars_h,
                       Rcpp::NumericVector s_1, Rcpp::NumericVector s_2) {
  
  // define a matrix to collect filtered states
  arma::mat y = arma::zeros(Nfcs, Nobs+1) ;
  y.col(0) = y0 ;
  
  // define an array to collect var-covar of filtered states
  arma::cube Sigma = arma::zeros(Nfcs, Nfcs, Nobs+1) ;
  
  // define an array to collect state-dependent Q
  arma::mat Q = arma::zeros(Nfcs, Nfcs) ;
  
  // initialize a vector to store data innovations within the loop
  arma::vec u = arma::zeros(Nbonds, 1) ; 
  
  // initialize a matrix to store var-covar of data innovation within the loop (+ its inverse)
  arma::mat Omega = arma::zeros(Nbonds, Nbonds) ;
  arma::mat invOmega = Omega ;
  
  // initialize Kalman gain
  arma::mat Kt = arma::zeros(Nfcs, Nbonds) ;
  
  // create an aux matrix to be used in the loop
  arma::mat It = arma::eye(Nfcs,Nfcs) ;
  
  // initialize a vector to store observations
  arma::vec obs = arma::zeros(Nbonds, 1);
  
  // initialize a vector to store gradient
  arma::vec B = arma::zeros(Nbonds, 1) ;
  
  // keep track of Log-Likelihood
  double LL = 0 ;
  
  // adjust starting values
  arma::vec linh = Price_D99_def_linearized_cpp(pars_b, arma::as_scalar(y.col(0)), s_1, s_2, pars_s, pars_h, 1) ; 
  arma::vec h = Price_D99_def_cpp(pars_b, arma::as_scalar(y.col(0)), s_1, s_2, pars_s, pars_h, 1) ;
  
  double linhcross = arma::as_scalar(linh.t()*linh) ;
  double hcrossinv = 1.0/arma::as_scalar(Phi.slice(0))/arma::as_scalar(Phi.slice(0)) ;
  
  y.col(0) = 1.0/arma::as_scalar(Phi.slice(0))*(linh.t()*(Rt.row(0).t() - h + 
    linh*arma::as_scalar(y.col(0)))/linhcross - a.slice(0)) ;
  Sigma.slice(0) = U*linhcross*hcrossinv + diagmat(w0.slice(0) + w1.slice(0) % y.col(0))*hcrossinv;
  
  // define a Kalman recursion
  for (int t = 1; t <= Nobs ; ++t) {
    
    // project the state
    y.col(t) = a.slice(t-1) + Phi.slice(t-1) * y.col(t-1) ;
    
    // calculate var-covar of state innovations that depend on y[,t-1]
    Q = diagmat(w0.slice(t-1) + w1.slice(t-1) % y.col(t-1)) ;
    
    // project var-covar of the state
    Sigma.slice(t) = Phi.slice(t-1) * Sigma.slice(t-1) * Phi.slice(t-1).t() + Q ;
    
    // calculate innovation in the observation relative to expectation
    obs = Rt.row(t-1).t();
    // u = obs - Rcpp::as<arma::vec>(FUN(y.col(t-1), t-1)) ;
    u = obs - Rcpp::as<arma::vec>(Price_D99_def_cpp(pars_b, arma::as_scalar(y.col(t)), s_1, s_2, pars_s, pars_h, t));
    
    // linearize the measurement
    //B = Rcpp::as<Rcpp::NumericVector>(D(FUN, y.col(t-1), t-1)) ;
    B = Price_D99_def_linearized_cpp(pars_b, arma::as_scalar(y.col(t)), s_1, s_2, pars_s, pars_h, t) ;
    
    // track the var-covar of this innovation
    Omega = B.t() * B * Sigma.slice(t) + U ;
    invOmega = arma::inv(Omega);
    
    // calculate Kalman gain
    Kt = arma::as_scalar(Sigma.slice(t) * invOmega) * B.t()  ;
    
    // update the state
    y.col(t) = y.col(t) + Kt * u ;
    
    // replace negative predcited state with zeroes
    arma::vec tmp = y.col(t);
    tmp.transform( [](double val) { return (val < 0.0) ? 0.0 : val; } );
    y.col(t) = tmp;
    
    // update the var-covar of the state
    Sigma.slice(t) = (It - Kt * B) * Sigma.slice(t) ;
    
    // update negative log-likelihood (skip  M/2*log(2*pi) part); Omega_det is a log determinant
    LL += 0.5*log(det(Omega)) + 0.5*arma::as_scalar(u.t()*u*invOmega) ;
    
  }
  
  Rcpp::List ret ;
  ret["y"] = y ;
  ret["Sigma"] = Sigma ;
  ret["LL"] = LL ;
  
  return(ret) ;
}

// High-level function of data and parameters that calls the stuff from above and runs an EKF Filter on Duffee (99)
// [[Rcpp::export]]

Rcpp::List D99_EKF_run_cpp(Rcpp::List pars_s, Rcpp::List pars_h, Rcpp::List pars_b, arma::mat data, arma::vec days_ind,
                           double delta, int Nobs, int Nfcs, int Nbonds, double me,
                           Rcpp::NumericVector s_1, Rcpp::NumericVector s_2) {
  
  // define SS matrices and other SS model parameters
  arma::cube a_in = arma::zeros(1, 1, Nobs) ;
  a_in.subcube(0,0,0,0,0,Nobs-1) = Rcpp::as<double>(pars_h[4])/Rcpp::as<double>(pars_h[3])*
    (1.0-exp(-Rcpp::as<double>(pars_h[3])*delta*days_ind)) ;
  
  arma::cube Phi_in = arma::zeros(1,1,Nobs) ;
  Phi_in.subcube(0,0,0,0,0,Nobs-1) = exp(-Rcpp::as<double>(pars_h[3])*delta*days_ind) ;
  
  arma::cube w0_in = arma::zeros(1,1,Nobs) ;  
  w0_in.subcube(0,0,0,0,0,Nobs-1) = Rcpp::as<double>(pars_h[5])*Rcpp::as<double>(pars_h[5])*Rcpp::as<double>(pars_h[4])/
    (2.0*Rcpp::as<double>(pars_h[3])*Rcpp::as<double>(pars_h[3]))*
    (1.0-exp(-Rcpp::as<double>(pars_h[3])*delta*days_ind)) %
    (1.0-exp(-Rcpp::as<double>(pars_h[3])*delta*days_ind)) ;
  
  arma::cube w1_in = arma::zeros(1,1,Nobs) ;    
  w1_in.subcube(0,0,0,0,0,Nobs-1) = Rcpp::as<double>(pars_h[5])*Rcpp::as<double>(pars_h[5])/Rcpp::as<double>(pars_h[3])*
      (exp(-Rcpp::as<double>(pars_h[3])*delta*days_ind)-
      exp(-2.0*Rcpp::as<double>(pars_h[3])*delta*days_ind)) ;
  
  arma::mat y0_in = arma::zeros(1,1) ;
  y0_in = Rcpp::as<double>(pars_h[4])/Rcpp::as<double>(pars_h[3]);  
  
  arma::mat U = arma::zeros(1,1) ;
  U(0,0) = me*me;
  
  // run the filter for given parameters
  Rcpp::List kfrun = EKF_run_cpp(data, a_in, Phi_in, U, y0_in, w0_in, w1_in, Nobs, Nfcs, Nbonds, delta,
                                 pars_b, pars_s, pars_h, s_1, s_2) ;
  return(kfrun) ;
  
}