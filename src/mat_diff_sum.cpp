#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(.mat_diff_sum)]]
double mat_diff_sum(Rcpp::NumericMatrix mat1, Rcpp::NumericMatrix mat2, double scale = 1, int delta_x = 0, int delta_y = 0){
  int nrow = mat1.nrow();
  int ncol = mat1.ncol();
  double mat_sum = 0;
  Rcpp::NumericMatrix m(nrow,ncol);
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      if(i - delta_x >= 0 && i - delta_x < nrow){
        if(j - delta_y >= 0 && j - delta_y < ncol){
          if(std::isfinite(mat1(i,j)) && std::isfinite(mat2(i - delta_x, j - delta_y))){
            mat_sum += pow((mat1(i, j) - mat2(i - delta_x, j - delta_y))/scale,2);
          }
        }
      }
    }
  }
  return mat_sum;
}
