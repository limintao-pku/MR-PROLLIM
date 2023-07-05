#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix m_m3(NumericMatrix x, NumericMatrix y, int n){
    NumericMatrix out(n,n);
    for (int i = 0; i < n; ++i){
      for (int j = 0; j < n; ++j){
        out(i,j) = x(i,j) + y(i,j);
      }
    }
    return(out);
}

double mydnorm2(double x1, double x2, double u1, double u2,
              double s12i, double s22i, double s1s2ri, double k1, double k2){
    double x1c = x1 - u1;
    double x2c = x2 - u2;
    double out = k1 * exp(k2 * (pow(x1c,2)*s12i - x1c*x2c*s1s2ri + pow(x2c,2)*s22i));
    return(out);
}

// [[Rcpp::export]]
NumericVector loop_p3_1_bi(double p1_sp, double p2_sp, double b1_t, double u1, double u2, double p1, double p2, double log_appr,
                NumericMatrix m_matrix, NumericMatrix s_h, 
                NumericMatrix post_sample_k1, NumericMatrix post_sample_k1_exp1,
                NumericMatrix post_sample_k2, NumericMatrix post_sample_k2_exp1, 
                NumericMatrix post_sample_p,
                NumericMatrix g1_matr, NumericMatrix g2_matr, List sigma_prime_list,
                IntegerVector sign_k1, IntegerVector sign_k2){
    double pi = 3.1415926535898;
    int n = sigma_prime_list.length();
    NumericMatrix spl[n];
    for (int i = 0; i < n; ++i) {
      NumericMatrix z0 = as<NumericMatrix>(sigma_prime_list[i]);
      spl[i]=z0;
    }
    int nr = m_matrix.nrow(), nc = post_sample_k1.ncol();  
    NumericVector pd(nr), pun(nr), pno(nr);
    double tu2 = 2 * u2;
    double s12i, s22i, s1s2ri, r2i, k1, k2, mu1, mu2, m;
    double z1[nc], z2[nc], z1_log[nc], z2_log[nc];
    double z11, z22, z1_m, z2_m, z1_m_log, z2_m_log, z1_z, z2_z;
    for (int i = 0; i < nr; ++i){
      if (log_appr == 2) {
        z11 = 0;
        z22 = 0;
        for (int j = 0; j < nc; ++j){
          z1[j] = post_sample_k1_exp1(i,j) * b1_t * post_sample_p(i,j) + 1;
          z2[j] = post_sample_k2_exp1(i,j) * b1_t * post_sample_p(i,j) + 1;
          z11 += z1[j]; 
          z22 += z2[j];
        }
        z1_m = z11 / nc;
        z2_m = z22 / nc;
        z1_m_log = log(z1_m);
        z2_m_log = log(z2_m);
        for (int j = 0; j < nc; ++j){
          z1_z = z1[j] / z1_m -1;
          z2_z = z2[j] / z2_m -1;
          z1_log[j] = z1_m_log + z1_z - pow(z1_z,2);
          z2_log[j] = z2_m_log + z2_z - pow(z2_z,2); 
        }  
      } else if (log_appr == 1) {
        z11 = 0;
        z22 = 0;
        for (int j = 0; j < nc; ++j){
          z1[j] = post_sample_k1_exp1(i,j) * b1_t * post_sample_p(i,j) + 1;
          z2[j] = post_sample_k2_exp1(i,j) * b1_t * post_sample_p(i,j) + 1;
          z11 += z1[j]; 
          z22 += z2[j];
        }
        z1_m = z11 / nc;
        z2_m = z22 / nc;
        z1_m_log = log(z1_m);
        z2_m_log = log(z2_m);
        for (int j = 0; j < nc; ++j){
          z1_log[j] = z1_m_log + z1[j] / z1_m -1;
          z2_log[j] = z2_m_log + z2[j] / z2_m -1; 
        } 
      } else {
        for (int j = 0; j < nc; ++j){
          z1_log[j] = log(post_sample_k1_exp1(i,j) * b1_t * post_sample_p(i,j) + 1);
          z2_log[j] = log(post_sample_k2_exp1(i,j) * b1_t * post_sample_p(i,j) + 1);
        } 
      }

      if (p2_sp == 0) {
        pd[i] = 0;
        NumericMatrix vcov = m_m3(s_h, spl[i], 2);
        s12i = 1/vcov(0,0);
        s22i = 1/vcov(1,1);
        s1s2ri = 2*vcov(0,1)*s12i*s22i;
        r2i = 1 - s1s2ri*vcov(0,1)/2;
        k1 = 1/(2*pi*pow(vcov(0,0)*vcov(1,1)*r2i,0.5));
        k2 = -1/(2*r2i);
      } else {
        NumericMatrix vcov = m_m3(s_h, spl[i], 2);
        s12i = 1/vcov(0,0);
        s22i = 1/vcov(1,1);
        s1s2ri = 2*vcov(0,1)*s12i*s22i;
        r2i = 1 - s1s2ri*vcov(0,1)/2;
        k1 = 1/(2*pi*pow(vcov(0,0)*vcov(1,1)*r2i,0.5));
        k2 = -1/(2*r2i);
        m = 0;
        for (int j = 0; j < nc; ++j){
          mu1 = u1 * post_sample_k1(i,j) + u2 * sign_k1[i] + g1_matr(i,j) + z1_log[j];
          mu2 = u1 * post_sample_k2(i,j) + tu2 * sign_k2[i] + g2_matr(i,j) + z2_log[j];
          m += mydnorm2(m_matrix(i,0), m_matrix(i,1), mu1, mu2, s12i, s22i, s1s2ri, k1, k2);
        }
        pd[i] = m / nc;
      }

      if (p2_sp == 1) {
        pun[i] = 0;
      } else {
        m = 0;
        for (int j = 0; j < nc; ++j){
          mu1 = u2 * sign_k1[i] + g1_matr(i,j) + z1_log[j];
          mu2 = tu2 * sign_k2[i] + g2_matr(i,j) + z2_log[j];
          m += mydnorm2(m_matrix(i,0), m_matrix(i,1), mu1, mu2, s12i, s22i, s1s2ri, k1, k2);
        }
        pun[i] = m / nc;
      }

      if (p1_sp == 1) {
        pno[i] = 0;
      } else {
        NumericMatrix vcov = spl[i];
        s12i = 1/vcov(0,0);
        s22i = 1/vcov(1,1);
        s1s2ri = 2*vcov(0,1)*s12i*s22i;
        r2i = 1 - s1s2ri*vcov(0,1)/2;
        k1 = 1/(2*pi*pow(vcov(0,0)*vcov(1,1)*r2i,0.5));
        k2 = -1/(2*r2i);
        m = 0;
        for (int j = 0; j < nc; ++j){
          mu1 = g1_matr(i,j) + z1_log[j];
          mu2 = g2_matr(i,j) + z2_log[j];
          m += mydnorm2(m_matrix(i,0), m_matrix(i,1), mu1, mu2, s12i, s22i, s1s2ri, k1, k2);
        }
        pno[i] = m / nc;
      }  
    }
    NumericVector out = p1*p2*pd + p1*(1-p2)*pun + (1-p1)*pno;
    return(out);
}

// [[Rcpp::export]]
NumericMatrix loop_p3_1_bi_indiv_p(double p1_sp, double p2_sp, double b1_t, double u1, double u2, double p1, double p2, double log_appr,
                NumericMatrix m_matrix, NumericMatrix s_h, 
                NumericMatrix post_sample_k1, NumericMatrix post_sample_k1_exp1,
                NumericMatrix post_sample_k2, NumericMatrix post_sample_k2_exp1, 
                NumericMatrix post_sample_p,
                NumericMatrix g1_matr, NumericMatrix g2_matr, List sigma_prime_list,
                IntegerVector sign_k1, IntegerVector sign_k2){
    double pi = 3.1415926535898;
    int n = sigma_prime_list.length();
    NumericMatrix spl[n];
    for (int i = 0; i < n; ++i) {
      NumericMatrix z0 = as<NumericMatrix>(sigma_prime_list[i]);
      spl[i]=z0;
    }
    int nr = m_matrix.nrow(), nc = post_sample_k1.ncol();  
    NumericMatrix pd(nr,nc), pun(nr,nc), pno(nr,nc);
    double tu2 = 2 * u2;
    double s12i, s22i, s1s2ri, r2i, k1, k2, mu1, mu2;
    double z1[nc], z2[nc], z1_log[nc], z2_log[nc];
    double z11, z22, z1_m, z2_m, z1_m_log, z2_m_log, z1_z, z2_z;

    for (int i = 0; i < nr; ++i){
      if (log_appr == 2) {
        z11 = 0;
        z22 = 0;
        for (int j = 0; j < nc; ++j){
          z1[j] = post_sample_k1_exp1(i,j) * b1_t * post_sample_p(i,j) + 1;
          z2[j] = post_sample_k2_exp1(i,j) * b1_t * post_sample_p(i,j) + 1;
          z11 += z1[j]; 
          z22 += z2[j];
        }
        z1_m = z11 / nc;
        z2_m = z22 / nc;
        z1_m_log = log(z1_m);
        z2_m_log = log(z2_m);
        for (int j = 0; j < nc; ++j){
          z1_z = z1[j] / z1_m -1;
          z2_z = z2[j] / z2_m -1;
          z1_log[j] = z1_m_log + z1_z - pow(z1_z,2);
          z2_log[j] = z2_m_log + z2_z - pow(z2_z,2); 
        }  
      } else if (log_appr == 1) {
        z11 = 0;
        z22 = 0;
        for (int j = 0; j < nc; ++j){
          z1[j] = post_sample_k1_exp1(i,j) * b1_t * post_sample_p(i,j) + 1;
          z2[j] = post_sample_k2_exp1(i,j) * b1_t * post_sample_p(i,j) + 1;
          z11 += z1[j]; 
          z22 += z2[j];
        }
        z1_m = z11 / nc;
        z2_m = z22 / nc;
        z1_m_log = log(z1_m);
        z2_m_log = log(z2_m);
        for (int j = 0; j < nc; ++j){
          z1_log[j] = z1_m_log + z1[j] / z1_m -1;
          z2_log[j] = z2_m_log + z2[j] / z2_m -1; 
        } 
      } else {
        for (int j = 0; j < nc; ++j){
          z1_log[j] = log(post_sample_k1_exp1(i,j) * b1_t * post_sample_p(i,j) + 1);
          z2_log[j] = log(post_sample_k2_exp1(i,j) * b1_t * post_sample_p(i,j) + 1);
        } 
      }

      if (p2_sp == 0) {
        for (int j = 0; j < nc; ++j){
          pd(i,j) = 0;
        }
        NumericMatrix vcov = m_m3(s_h, spl[i], 2);
        s12i = 1/vcov(0,0);
        s22i = 1/vcov(1,1);
        s1s2ri = 2*vcov(0,1)*s12i*s22i;
        r2i = 1 - s1s2ri*vcov(0,1)/2;
        k1 = 1/(2*pi*pow(vcov(0,0)*vcov(1,1)*r2i,0.5));
        k2 = -1/(2*r2i);
      } else {
        NumericMatrix vcov = m_m3(s_h, spl[i], 2);
        s12i = 1/vcov(0,0);
        s22i = 1/vcov(1,1);
        s1s2ri = 2*vcov(0,1)*s12i*s22i;
        r2i = 1 - s1s2ri*vcov(0,1)/2;
        k1 = 1/(2*pi*pow(vcov(0,0)*vcov(1,1)*r2i,0.5));
        k2 = -1/(2*r2i);
        for (int j = 0; j < nc; ++j){
          mu1 = u1 * post_sample_k1(i,j) + u2 * sign_k1[i] + g1_matr(i,j) + z1_log[j];
          mu2 = u1 * post_sample_k2(i,j) + tu2 * sign_k2[i] + g2_matr(i,j) + z2_log[j];
          pd(i,j) = mydnorm2(m_matrix(i,0), m_matrix(i,1), mu1, mu2, s12i, s22i, s1s2ri, k1, k2);
        }
      }

      if (p2_sp == 1) {
        for (int j = 0; j < nc; ++j){
          pun(i,j) = 0;
        }
      } else {
        for (int j = 0; j < nc; ++j){
          mu1 = u2 * sign_k1[i] + g1_matr(i,j) + z1_log[j];
          mu2 = tu2 * sign_k2[i] + g2_matr(i,j) + z2_log[j];
          pun(i,j) = mydnorm2(m_matrix(i,0), m_matrix(i,1), mu1, mu2, s12i, s22i, s1s2ri, k1, k2);
        }
      }

      if (p1_sp == 1) {
        for (int j = 0; j < nc; ++j){
          pno(i,j) = 0;
        }
      } else {
        NumericMatrix vcov = spl[i];
        s12i = 1/vcov(0,0);
        s22i = 1/vcov(1,1);
        s1s2ri = 2*vcov(0,1)*s12i*s22i;
        r2i = 1 - s1s2ri*vcov(0,1)/2;
        k1 = 1/(2*pi*pow(vcov(0,0)*vcov(1,1)*r2i,0.5));
        k2 = -1/(2*r2i);
        for (int j = 0; j < nc; ++j){
          mu1 = g1_matr(i,j) + z1_log[j];
          mu2 = g2_matr(i,j) + z2_log[j];
          pno(i,j) = mydnorm2(m_matrix(i,0), m_matrix(i,1), mu1, mu2, s12i, s22i, s1s2ri, k1, k2);
        }
      }  
    }

    NumericMatrix out(nr,nc);
    for (int i = 0; i < nr; ++i){
      for (int j = 0; j < nc; ++j){
        out(i,j) = p1*p2*pd(i,j) + p1*(1-p2)*pun(i,j) + (1-p1)*pno(i,j);
      }
    }
    return(out);
} 