#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector my_allocate_cpp2(NumericMatrix x, NumericMatrix y){
    int ynrow = y.nrow(), yncol = y.ncol(), xncol = x.ncol();
    NumericVector out(xncol);
    std::vector<int> loc;

    for (int i = 0; i < ynrow; ++i) {
        NumericVector out0(xncol);
        for (int k = 0; k < xncol; ++k) {
            for (int j = 0; j < yncol; ++j) {
                out0[k] += pow(y(i, j) - x(j, k), 2);
            }
        }
        loc.clear();
        double m = out0[0] + 1;
        for (int q = 0; q < xncol; ++q) {
            if (out0[q] < m) {
                m = out0[q];
                loc.clear();
                loc.push_back(q);
            } else if (out0[q] == m) {
                loc.push_back(q);
            }
        }
        int n = loc.size();
        double v = 1;
        v = v / n;
        for (int q = 0; q < n; ++q) {
            int loc0 = loc[q];
            out[loc0] += v;
        }
    }
    return out;
}