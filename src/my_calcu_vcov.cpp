#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix my_calcu_vcov(NumericMatrix x_c, IntegerVector y){
        int nr = x_c.nrow(), nc = x_c.ncol();
        NumericMatrix out(nc, nc);

        for (int i = 0; i < nc; ++i) {
          for (int k = i; k < nc; ++k) {
            NumericVector out0(nr);
            for (int j = 0; j < nr; ++j) {
              out0[j] = x_c(j, i) * x_c(j, k);
            }
            LogicalVector na = !is_na(out0);
            NumericVector out00 = out0[na];
            double s = sum(na);
            out(i, k) = sum(out00) / (s - 1) / y[i] / y[k] * s;
          }
        }
        for (int i = 1; i < nc; ++i) {
          for (int k = 0; k < i; ++k) {
            out(i, k) = out(k, i);
          }
        }
        return out;
      }