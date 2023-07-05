#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix m_m1(NumericMatrix x, NumericMatrix y, int n){
    NumericMatrix out(n,n);
    for (int i = 0; i < n; ++i){
      for (int j = 0; j < n; ++j){
        out(i,j) = x(i,j) + y(i,j);
      }
    }
    return(out);
}

double loop_mdn(double x1, double x2, double u1, double u2, NumericMatrix s){
  double pi = 3.1415926535898;
  double x1c = x1 - u1;
  double x2c = x2 - u2;
  double rr = s(0,1)/s(0,0)/s(1,1);
  double rrr = 1 - rr*s(0,1);
  double out = exp(-(pow(x1c,2)/s(0,0) - 2*rr*x1c*x2c + pow(x2c,2)/s(1,1)) / 2 / rrr) / (2*pi*pow(rrr*s(0,0)*s(1,1),0.5));
  return(out);
}

// [[Rcpp::export]]
NumericVector f_est_k_prior_loop(NumericMatrix k_matr, List s_list, List ss_list,
  double mychisq, double p0_sp, double p0, double p_cut, double u1, double u2, NumericMatrix sigma){
    Environment pkg = Environment::namespace_env("MRprollim");
    Function f = pkg["my_pmvnEll"];
    int n = k_matr.nrow();
    NumericMatrix sl[n];
    for (int i = 0; i < n; ++i) {
      NumericMatrix z0 = as<NumericMatrix>(s_list[i]);
      sl[i]=z0;
    }
    NumericMatrix ssl[n];
    for (int i = 0; i < n; ++i) {
      NumericMatrix z0 = as<NumericMatrix>(ss_list[i]);
      ssl[i]=z0;
    }
    NumericVector mu(2);
    mu[0] = u1;
    mu[1] = u2;
    NumericVector x0(2);
    NumericVector p(n);

    if (p0_sp == -1) {
      for (int i = 0; i < n; ++i) {
        NumericMatrix s = m_m1(sigma, sl[i], 2);
        NumericVector p1 = f(mychisq, s, mu, ssl[i], x0, false);
        double p10 = p1[0];
        p[i] = loop_mdn(k_matr(i,0),k_matr(i,1),u1,u2,s) / p10;
      }
    } else {
      for (int i = 0; i < n; ++i) {
        NumericMatrix s = m_m1(sigma, sl[i], 2);
        NumericVector p1 = f(mychisq, s, mu, ssl[i], x0, false);
        double p10 = p1[0];
        double p00 = p0 * p_cut / (p0 * p_cut + (1 - p0) * p10);
        p[i] = p00 * loop_mdn(k_matr(i,0),k_matr(i,1),0,0,sl[i]) / p_cut + (1 - p00) * loop_mdn(k_matr(i,0),k_matr(i,1),u1,u2,s) / p10;
      }
    }
    return(p);
}