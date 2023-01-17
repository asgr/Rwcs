#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".dilate_cpp")]]
IntegerMatrix dilate_cpp(IntegerMatrix segim, IntegerMatrix kern, IntegerVector expand=0){
  
  int srow = segim.nrow();
  int scol = segim.ncol();
  int krow = kern.nrow();
  int kcol = kern.ncol();
  int krow_off = ((krow - 1) / 2);
  int kcol_off = ((kcol - 1) / 2);
  bool checkseg = true;
  int max_segim = max(segim);
  IntegerMatrix segim_new(srow, scol);
  LogicalVector seglogic  (max_segim);
  
  // if(expand(0) == 0){
  //   expandall = true;
  // }else{
  //   expandall = false;
  //   
  // }
  
  if(expand(0) > 0){
    for(int k = 0; k < expand.length(); k++){
      if(expand(k) <= max_segim){
        seglogic(expand(k)-1) = true;
      }
    }
  }
  
  for (int i = 0; i < srow; i++) {
    for (int j = 0; j < scol; j++) {
      if(segim(i,j) > 0){
        if(expand(0) > 0){
          checkseg = seglogic(segim(i,j)-1);
        }
        if(checkseg){
          for (int m = std::max(0, krow_off - i); m < std::min(krow, krow_off - (i - srow)); m++) {
            for (int n = std::max(0, kcol_off - j); n < std::min(kcol, kcol_off - (j - scol)); n++) {
              if(kern(m,n) > 0){
                if(m != krow_off || n != kcol_off){
                  if(segim(i + m - krow_off,j + n - kcol_off) == 0){
                    int xloc = i + m - krow_off;
                    int yloc = j + n - kcol_off;
                    if(segim(i,j) < segim_new(xloc,yloc) || segim_new(xloc,yloc) == 0) {
                      segim_new(xloc,yloc) = segim(i,j);
                    }
                  }
                }else{
                  segim_new(i,j) = segim(i,j);
                }
              }
            }
          }
        }else{
          segim_new(i,j) = segim(i,j);
        }
      }
    }
  }
  return segim_new;
}

// [[Rcpp::export(".num_mat_add_cpp")]]
NumericMatrix num_mat_add(NumericMatrix base, NumericMatrix add, IntegerMatrix ind, IntegerVector offset){
  
  int row, col;
  //Rcout << ind.nrow() << std::endl;
  for(int i = 0; i < ind.nrow(); i++){
    //Rcout << i << std::endl;
    row = ind(i,0) - 1; // To convert from R to C index
    col = ind(i,1) - 1; // To convert from R to C index
    
    base(row + offset[0] - 1, col + offset[1] - 1) += add(row, col);
  }
  
  return base;
}

// [[Rcpp::export(".num_mat_add_mult_cpp")]]
NumericMatrix num_mat_add_mult(NumericMatrix base, NumericMatrix add, NumericMatrix mult, IntegerMatrix ind, IntegerVector offset){
  
  int row, col;
  //Rcout << ind.nrow() << std::endl;
  for(int i = 0; i < ind.nrow(); i++){
    //Rcout << i << std::endl;
    row = ind(i,0) - 1; // To convert from R to C index
    col = ind(i,1) - 1; // To convert from R to C index
    
    //Rcout << row << ' ' << col << ' ' << row + offset[0] - 1 << ' ' << col + offset[1] - 1 << std::endl;
    
    base(row + offset[0] - 1, col + offset[1] - 1) += add(row, col) * mult(row, col);
  }
  
  return base;
}

// [[Rcpp::export(".int_mat_add_cpp")]]
IntegerMatrix int_mat_add(IntegerMatrix base, IntegerMatrix add, IntegerMatrix ind, IntegerVector offset){
  
  int row, col;
  //Rcout << ind.nrow() << std::endl;
  for(int i = 0; i < ind.nrow(); i++){
    //Rcout << i << std::endl;
    row = ind(i,0) - 1; // To convert from R to C index
    col = ind(i,1) - 1; // To convert from R to C index
    
    base(row + offset[0] - 1, col + offset[1] - 1) += add(row, col);
  }
  
  return base;
}

// [[Rcpp::export(".int_mat_add_sin_cpp")]]
IntegerMatrix int_mat_add_sin(IntegerMatrix base, int add, IntegerMatrix ind, IntegerVector offset){
  
  int row, col;
  //Rcout << ind.nrow() << std::endl;
  for(int i = 0; i < ind.nrow(); i++){
    //Rcout << i << std::endl;
    row = ind(i,0) - 1; // To convert from R to C index
    col = ind(i,1) - 1; // To convert from R to C index
    
    base(row + offset[0] - 1, col + offset[1] - 1) += add;
  }
  
  return base;
}

// [[Rcpp::export(".stack_image_inVar_weight_mat_cpp")]]
SEXP stack_image_inVar_weight_mat(NumericMatrix post_image, NumericMatrix post_inVar, IntegerMatrix post_weight,
              NumericMatrix pre_image, NumericMatrix pre_inVar, NumericMatrix pre_weight, IntegerVector offset){
  int do_stack;
  
  for (int j = 0; j < pre_image.ncol(); j++) {
    for (int i = 0; i < pre_image.nrow(); i++) {
      do_stack = R_finite(pre_image(i,j) + pre_inVar(i,j)) && (pre_inVar(i,j) > 0);
      if(do_stack){
        post_image(i + offset[0] - 1, j + offset[1] - 1) += pre_image(i,j) * pre_inVar(i,j);
        post_inVar(i + offset[0] - 1, j + offset[1] - 1) += pre_inVar(i,j);
        post_weight(i + offset[0] - 1, j + offset[1] - 1) += pre_weight(i,j);
      }
    }
  }
  return R_NilValue;
}

// [[Rcpp::export(".stack_image_inVar_weight_int_cpp")]]
SEXP stack_image_inVar_weight_int(NumericMatrix post_image, NumericMatrix post_inVar, IntegerMatrix post_weight,
                            NumericMatrix pre_image, NumericMatrix pre_inVar, int pre_weight, IntegerVector offset){
  int do_stack;
  
  for (int j = 0; j < pre_image.ncol(); j++) {
    for (int i = 0; i < pre_image.nrow(); i++) {
      do_stack = R_finite(pre_image(i,j) + pre_inVar(i,j)) && (pre_inVar(i,j) > 0);
      if(do_stack){
        post_image(i + offset[0] - 1, j + offset[1] - 1) += pre_image(i,j) * pre_inVar(i,j);
        post_inVar(i + offset[0] - 1, j + offset[1] - 1) += pre_inVar(i,j);
        post_weight(i + offset[0] - 1, j + offset[1] - 1) += pre_weight;
      }
    }
  }
  return R_NilValue;
}

// [[Rcpp::export(".stack_image_weight_mat_cpp")]]
SEXP stack_image_weight_mat(NumericMatrix post_image, IntegerMatrix post_weight,
                            NumericMatrix pre_image, NumericMatrix pre_weight, IntegerVector offset){
  int do_stack;
  
  for (int j = 0; j < pre_image.ncol(); j++) {
    for (int i = 0; i < pre_image.nrow(); i++) {
      do_stack = R_finite(pre_image(i,j));
      if(do_stack){
        post_image(i + offset[0] - 1, j + offset[1] - 1) += pre_image(i,j) * pre_weight(i,j);
        post_weight(i + offset[0] - 1, j + offset[1] - 1) += pre_weight(i,j);
      }
    }
  }
  return R_NilValue;
}

// [[Rcpp::export(".stack_image_weight_int_cpp")]]
SEXP stack_image_weight_int(NumericMatrix post_image, IntegerMatrix post_weight,
                            NumericMatrix pre_image, int pre_weight, IntegerVector offset){
  int do_stack;
  
  for (int j = 0; j < pre_image.ncol(); j++) {
    for (int i = 0; i < pre_image.nrow(); i++) {
      do_stack = R_finite(pre_image(i,j));
      if(do_stack){
        post_image(i + offset[0] - 1, j + offset[1] - 1) += pre_image(i,j) * pre_weight;
        post_weight(i + offset[0] - 1, j + offset[1] - 1) += pre_weight;
      }
    }
  }
  return R_NilValue;
}

// [[Rcpp::export(".stack_exp_cpp")]]
SEXP stack_exp(NumericMatrix post_exp, NumericMatrix pre_exp, IntegerVector offset){
  int do_stack;
  
  for (int j = 0; j < pre_exp.ncol(); j++) {
    for (int i = 0; i < pre_exp.nrow(); i++) {
      do_stack = R_finite(pre_exp(i,j));
      if(do_stack){
        post_exp(i + offset[0] - 1, j + offset[1] - 1) += pre_exp(i,j);
      }
    }
  }
  return R_NilValue;
}

//masking versions

// [[Rcpp::export(".image_inVar_weight_mat_mask_cpp")]]
SEXP image_inVar_weight_mat_mask(NumericMatrix post_image, NumericMatrix post_inVar, IntegerMatrix post_weight,
                            NumericMatrix pre_image, NumericMatrix pre_inVar, NumericMatrix pre_weight, IntegerVector offset,
                            LogicalMatrix post_mask){
  int do_stack;
  
  for (int j = 0; j < pre_image.ncol(); j++) {
    for (int i = 0; i < pre_image.nrow(); i++) {
      do_stack = R_finite(pre_image(i,j) + pre_inVar(i,j)) && (pre_inVar(i,j) > 0) && (post_mask(i + offset[0] - 1, j + offset[1] - 1) == false);
      if(do_stack){
        post_image(i + offset[0] - 1, j + offset[1] - 1) += pre_image(i,j) * pre_inVar(i,j);
        post_inVar(i + offset[0] - 1, j + offset[1] - 1) += pre_inVar(i,j);
        post_weight(i + offset[0] - 1, j + offset[1] - 1) += pre_weight(i,j);
      }
    }
  }
  return R_NilValue;
}

// [[Rcpp::export(".stack_image_inVar_weight_int_mask_cpp")]]
SEXP stack_image_inVar_weight_int_mask(NumericMatrix post_image, NumericMatrix post_inVar, IntegerMatrix post_weight,
                                  NumericMatrix pre_image, NumericMatrix pre_inVar, int pre_weight, IntegerVector offset,
                                  LogicalMatrix post_mask){
  int do_stack;
  
  for (int j = 0; j < pre_image.ncol(); j++) {
    for (int i = 0; i < pre_image.nrow(); i++) {
      do_stack = R_finite(pre_image(i,j) + pre_inVar(i,j)) && (pre_inVar(i,j) > 0) && (post_mask(i + offset[0] - 1, j + offset[1] - 1) == false);
      if(do_stack){
        post_image(i + offset[0] - 1, j + offset[1] - 1) += pre_image(i,j) * pre_inVar(i,j);
        post_inVar(i + offset[0] - 1, j + offset[1] - 1) += pre_inVar(i,j);
        post_weight(i + offset[0] - 1, j + offset[1] - 1) += pre_weight;
      }
    }
  }
  return R_NilValue;
}

// [[Rcpp::export(".stack_image_weight_mat_mask_cpp")]]
SEXP stack_image_weight_mat_mask(NumericMatrix post_image, IntegerMatrix post_weight,
                            NumericMatrix pre_image, NumericMatrix pre_weight, IntegerVector offset,
                            LogicalMatrix post_mask){
  int do_stack;
  
  for (int j = 0; j < pre_image.ncol(); j++) {
    for (int i = 0; i < pre_image.nrow(); i++) {
      do_stack = R_finite(pre_image(i,j)) && (post_mask(i + offset[0] - 1, j + offset[1] - 1) == false);
      if(do_stack){
        post_image(i + offset[0] - 1, j + offset[1] - 1) += pre_image(i,j) * pre_weight(i,j);
        post_weight(i + offset[0] - 1, j + offset[1] - 1) += pre_weight(i,j);
      }
    }
  }
  return R_NilValue;
}

// [[Rcpp::export(".stack_image_weight_int_mask_cpp")]]
SEXP stack_image_weight_int_mask(NumericMatrix post_image, IntegerMatrix post_weight,
                            NumericMatrix pre_image, int pre_weight, IntegerVector offset,
                            LogicalMatrix post_mask){
  int do_stack;
  
  for (int j = 0; j < pre_image.ncol(); j++) {
    for (int i = 0; i < pre_image.nrow(); i++) {
      do_stack = R_finite(pre_image(i,j)) && (post_mask(i + offset[0] - 1, j + offset[1] - 1) == false);
      if(do_stack){
        post_image(i + offset[0] - 1, j + offset[1] - 1) += pre_image(i,j) * pre_weight;
        post_weight(i + offset[0] - 1, j + offset[1] - 1) += pre_weight;
      }
    }
  }
  return R_NilValue;
}

// [[Rcpp::export(".stack_exp_mask_cpp")]]
SEXP stack_exp_mask(NumericMatrix post_exp, NumericMatrix pre_exp, IntegerVector offset,
                    LogicalMatrix post_mask){
  int do_stack;
  
  for (int j = 0; j < pre_exp.ncol(); j++) {
    for (int i = 0; i < pre_exp.nrow(); i++) {
      do_stack = R_finite(pre_exp(i,j)) && (post_mask(i + offset[0] - 1, j + offset[1] - 1) == false);
      if(do_stack){
        post_exp(i + offset[0] - 1, j + offset[1] - 1) += pre_exp(i,j);
      }
    }
  }
  return R_NilValue;
}

