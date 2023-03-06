#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List my_boot_adj(int n_rep, int n_beta, NumericMatrix p_matr, List der1_list, List der2_list){
           int nr = p_matr.nrow(), nc = p_matr.ncol(), n0 = nr * n_beta;
           List out1(n_rep);
           List out2(n_rep);
           List out(2);
           NumericMatrix d1_l_p[n_beta];
           for (int z = 0; z < n_beta; ++z) {NumericMatrix z0 = as<NumericMatrix>(der1_list[z]); d1_l_p[z]=z0;}
           NumericMatrix d2_l_p[n_beta * n_beta];
           for (int z = 0; z < (n_beta * n_beta); ++z) {NumericMatrix z0 = as<NumericMatrix>(der2_list[z]); d2_l_p[z]=z0;}
           
           double m;
           for (int i = 0; i < n_rep; ++i) {
             NumericVector p_b(nr);
             NumericMatrix der1_matr(nr, n_beta);
             NumericVector der2_array(n0 * n_beta);
             for (int j = 0; j < nr; ++j) {
               IntegerVector id = sample(nc, nc, true) - 1;
               m = 0;
               for (int z = 0; z < nc; ++z) m += p_matr(j, id[z]);
               p_b[j] = m / nc;
               for(int k = 0; k < n_beta; ++k){
                 m = 0;
                 for (int z = 0; z < nc; ++z) m += (d1_l_p[k])(j, id[z]);
                 der1_matr(j, k) = m / nc;
                 for(int l = 0; l < n_beta; ++l){
                   m = 0;
                   for (int z = 0; z < nc; ++z) m += (d2_l_p[n_beta * k + l])(j, id[z]);
                   der2_array[j * n_beta + k + l * n0] = m / nc;
                 }
               }
             }
             NumericVector der1(n_beta); 
             NumericMatrix der2(n_beta, n_beta);
             for (int j = 0; j < n_beta; ++j) {
               m = 0;
               for (int z = 0; z < nr; ++z) m += der1_matr(z, j) / p_b[z];
               der1[j] = - m / nr;
               for(int k = 0; k < n_beta; ++k){
                 m = 0;
                 for (int z = 0; z < nr; ++z) m += - der1_matr(z, j) * der1_matr(z, k) / p_b[z] / p_b[z] + der2_array(z * n_beta + j + k * n0) / p_b[z];
                 der2(j, k) = - m / nr;
               }
             }
             out1[i] = der1;
             out2[i] = der2;
           }
           out[0] = out1;
           out[1] = out2;
           return(out);
           }