// While including system headrs, avoid clashes with the old wcsset POXIS
// function that is still exposed in Windows headers.
//
// This old function is *not* exposed if one defines the NO_OLDNAMES
// (MinGW-64bit) or _NO_OLDNAMES (MinGW-32bit) macros. In fact, that's what we
// do for the compilation of wcslib itself. However, these macros cause other
// names to be hidden, and make the C++ system headers internally incompatible,
// turning this into an unfeasible option for compiling this module.
#if defined(_WIN32) || defined(_MSC_VER) || defined(__MINGW32__) || defined (__MINGW64__)
#define wcsset wcsset_
#endif

#include <Rcpp.h>
#include <algorithm>
#include <utility>
#include <vector>

// End the hack above
#if defined(_WIN32) || defined(_MSC_VER) || defined(__MINGW32__) || defined (__MINGW64__)
#undef wcsset
#endif

#include <wcslib.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP Cwcs_s2p(Rcpp::NumericVector RA, Rcpp::NumericVector Dec,
              Rcpp::String CTYPE1 = "RA---TAN", Rcpp::String CTYPE2 = "DEC--TAN",
              double CRVAL1 = 0, double CRVAL2 = 0,
              double CRPIX1 = 0, double CRPIX2 = 0,
              double CD1_1 = 1, double CD1_2 = 0,
              double CD2_1 = 0, double CD2_2 = 1,
              Rcpp::String RADESYS = "ICRS", int EQUINOX = 2000,
              double PV1_1 = NA_REAL, double PV1_2 = NA_REAL, double PV1_3 = NA_REAL,
              double PV2_1 = NA_REAL, double PV2_2 = NA_REAL, double PV2_3 = NA_REAL,
              double PV2_4 = NA_REAL, double PV2_5 = NA_REAL
){
  int i;
  const int ncoord = RA.length();
  const int naxis = 2;
  
  //setup wcs
  struct wcsprm wcs;
  wcs.flag = -1;
  wcsini(1, naxis, &wcs);
  
  //insert wcs val
  wcs.crval[0] = CRVAL1;
  wcs.crval[1] = CRVAL2;
  
  //insert wcs pix
  wcs.crpix[0] = CRPIX1;
  wcs.crpix[1] = CRPIX2;
  
  //insert wcs pix
  wcs.crpix[0] = CRPIX1;
  wcs.crpix[1] = CRPIX2;
  
  //insert wcs cd matrix
  double m[2][2];
  m[0][0] = CD1_1;
  m[0][1] = CD1_2;
  m[1][0] = CD2_1;
  m[1][1] = CD2_2;
  wcs.pc = *m;
  
  //insert ctype
  strcpy(wcs.ctype[0], CTYPE1.get_cstring());
  strcpy(wcs.ctype[1], CTYPE2.get_cstring());
  
  //insert radesys and equinox
  strcpy(wcs.radesys, RADESYS.get_cstring());
  wcs.equinox = EQUINOX;
  
  //insert wcs pv
  wcs.npv = 0;
  if(!R_IsNA(PV1_1)){wcs.npv++;}
  if(!R_IsNA(PV1_2)){wcs.npv++;}
  if(!R_IsNA(PV1_3)){wcs.npv++;}
  if(!R_IsNA(PV2_1)){wcs.npv++;}
  if(!R_IsNA(PV2_2)){wcs.npv++;}
  if(!R_IsNA(PV2_3)){wcs.npv++;}
  if(!R_IsNA(PV2_4)){wcs.npv++;}
  if(!R_IsNA(PV2_5)){wcs.npv++;}
  
  int npvcount = 0;
  
  if(!R_IsNA(PV1_1)){
    wcs.pv[npvcount].i = 1;
    wcs.pv[npvcount].m = 1;
    wcs.pv[npvcount].value = PV1_1;
    npvcount++;
  }
  if(!R_IsNA(PV1_2)){
    wcs.pv[npvcount].i = 1;
    wcs.pv[npvcount].m = 2;
    wcs.pv[npvcount].value = PV1_2;
    npvcount++;
  }
  if(!R_IsNA(PV1_3)){
    wcs.pv[npvcount].i = 1;
    wcs.pv[npvcount].m = 3;
    wcs.pv[npvcount].value = PV1_3;
    npvcount++;
  }
  if(!R_IsNA(PV2_1)){
    wcs.pv[npvcount].i = 2;
    wcs.pv[npvcount].m = 1;
    wcs.pv[npvcount].value = PV2_1;
    npvcount++;
  }
  if(!R_IsNA(PV2_2)){
    wcs.pv[npvcount].i = 2;
    wcs.pv[npvcount].m = 2;
    wcs.pv[npvcount].value = PV2_2;
    npvcount++;
  }
  if(!R_IsNA(PV2_3)){
    wcs.pv[npvcount].i = 2;
    wcs.pv[npvcount].m = 3;
    wcs.pv[npvcount].value = PV2_3;
    npvcount++;
  }
  if(!R_IsNA(PV2_4)){
    wcs.pv[npvcount].i = 2;
    wcs.pv[npvcount].m = 4;
    wcs.pv[npvcount].value = PV2_4;
    npvcount++;
  }
  if(!R_IsNA(PV2_5)){
    wcs.pv[npvcount].i = 2;
    wcs.pv[npvcount].m = 5;
    wcs.pv[npvcount].value = PV2_5;
    npvcount++;
  }
  
  //set long and lat axis positions
  wcs.lng = 0;
  wcs.lat = 1;
  
  wcsset(&wcs);
  
  NumericMatrix world(naxis, ncoord);
  for (i = 0; i < ncoord; i++) {
    world(0, i) = RA[i];
    world(1, i) = Dec[i];
  }
  NumericVector phi(ncoord);
  NumericVector theta(ncoord);
  NumericMatrix img(naxis, ncoord);
  IntegerVector stat(ncoord);
  NumericMatrix pixel_matrix(naxis, ncoord);
  
  int status = wcss2p(&wcs, ncoord, naxis,
         &(world[0]), &(phi[0]), &(theta[0]), &(img[0]),
         &(pixel_matrix[0]), &(stat[0]));
         // TODO: check stat
  wcsfree(&wcs);
  if(status > 0 || stat[0] > 0){
    Rcout << "Failed s2p conversion!" << "\n";
    return(stat);
  }else{
    return(transpose(pixel_matrix));
  }
}

// [[Rcpp::export]]
SEXP Cwcs_p2s(Rcpp::NumericVector x, Rcpp::NumericVector y,
              Rcpp::String CTYPE1 = "RA---TAN", Rcpp::String CTYPE2 = "DEC--TAN",
              double CRVAL1 = 0, double CRVAL2 = 0,
              double CRPIX1 = 0, double CRPIX2 = 0,
              double CD1_1 = 1, double CD1_2 = 0,
              double CD2_1 = 0, double CD2_2 = 1,
              Rcpp::String RADESYS = "ICRS", int EQUINOX = 2000,
              double PV1_1 = NA_REAL, double PV1_2 = NA_REAL, double PV1_3 = NA_REAL,
              double PV2_1 = NA_REAL, double PV2_2 = NA_REAL, double PV2_3 = NA_REAL,
              double PV2_4 = NA_REAL, double PV2_5 = NA_REAL
){
  int i;
  const int ncoord = x.length();
  const int naxis = 2;
  
  //setup wcs
  struct wcsprm wcs;
  wcs.flag = -1;
  wcsini(1, naxis, &wcs);
  
  //insert wcs val
  wcs.crval[0] = CRVAL1;
  wcs.crval[1] = CRVAL2;
  
  //insert wcs pix
  wcs.crpix[0] = CRPIX1;
  wcs.crpix[1] = CRPIX2;
  
  //insert wcs pix
  wcs.crpix[0] = CRPIX1;
  wcs.crpix[1] = CRPIX2;
  
  //insert wcs cd matrix
  double m[2][2];
  m[0][0] = CD1_1;
  m[0][1] = CD1_2;
  m[1][0] = CD2_1;
  m[1][1] = CD2_2;
  wcs.pc = *m;
  
  //insert ctype
  strcpy(wcs.ctype[0], CTYPE1.get_cstring());
  strcpy(wcs.ctype[1], CTYPE2.get_cstring());
  
  //insert radesys and equinox
  strcpy(wcs.radesys, RADESYS.get_cstring());
  wcs.equinox = EQUINOX;
  
  //insert wcs pv
  wcs.npv = 0;
  if(!R_IsNA(PV1_1)){wcs.npv++;}
  if(!R_IsNA(PV1_2)){wcs.npv++;}
  if(!R_IsNA(PV1_3)){wcs.npv++;}
  if(!R_IsNA(PV2_1)){wcs.npv++;}
  if(!R_IsNA(PV2_2)){wcs.npv++;}
  if(!R_IsNA(PV2_3)){wcs.npv++;}
  if(!R_IsNA(PV2_4)){wcs.npv++;}
  if(!R_IsNA(PV2_5)){wcs.npv++;}
  
  int npvcount = 0;
  
  if(!R_IsNA(PV1_1)){
    wcs.pv[npvcount].i = 1;
    wcs.pv[npvcount].m = 1;
    wcs.pv[npvcount].value = PV1_1;
    npvcount++;
  }
  if(!R_IsNA(PV1_2)){
    wcs.pv[npvcount].i = 1;
    wcs.pv[npvcount].m = 2;
    wcs.pv[npvcount].value = PV1_2;
    npvcount++;
  }
  if(!R_IsNA(PV1_3)){
    wcs.pv[npvcount].i = 1;
    wcs.pv[npvcount].m = 3;
    wcs.pv[npvcount].value = PV1_3;
    npvcount++;
  }
  if(!R_IsNA(PV2_1)){
    wcs.pv[npvcount].i = 2;
    wcs.pv[npvcount].m = 1;
    wcs.pv[npvcount].value = PV2_1;
    npvcount++;
  }
  if(!R_IsNA(PV2_2)){
    wcs.pv[npvcount].i = 2;
    wcs.pv[npvcount].m = 2;
    wcs.pv[npvcount].value = PV2_2;
    npvcount++;
  }
  if(!R_IsNA(PV2_3)){
    wcs.pv[npvcount].i = 2;
    wcs.pv[npvcount].m = 3;
    wcs.pv[npvcount].value = PV2_3;
    npvcount++;
  }
  if(!R_IsNA(PV2_4)){
    wcs.pv[npvcount].i = 2;
    wcs.pv[npvcount].m = 4;
    wcs.pv[npvcount].value = PV2_4;
    npvcount++;
  }
  if(!R_IsNA(PV2_5)){
    wcs.pv[npvcount].i = 2;
    wcs.pv[npvcount].m = 5;
    wcs.pv[npvcount].value = PV2_5;
    npvcount++;
  }
  
  //set long and lat axis positions
  wcs.lng = 0;
  wcs.lat = 1;
  
  wcsset(&wcs);
  
  NumericMatrix pixel(naxis, ncoord);
  for (i = 0; i < ncoord; i++) {
    pixel(0, i) = x[i];
    pixel(1, i) = y[i];
  }
  NumericVector phi(ncoord);
  NumericVector theta(ncoord);
  NumericMatrix img(naxis, ncoord);
  IntegerVector stat(ncoord);
  NumericMatrix world_matrix(naxis, ncoord);
  
  int status = wcsp2s(&wcs, ncoord, naxis,
         &(pixel[0]), &(img[0]), &(phi[0]), &(theta[0]),
         &(world_matrix[0]), &(stat[0]));
  
  wcsfree(&wcs);
  if(status > 0 || stat[0] > 0){
    Rcout << "Failed p2s conversion!" << "\n";
    return(stat);
  }else{
    return(transpose(world_matrix));
  }
}

// [[Rcpp::export]]
SEXP Cwcs_head_p2s(Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::String header, 
                  int nkeyrec){
  int i;
  const int ncoord = x.length();
  const int naxis = 2;
  int nreject, nwcs;
  
  //setup wcs
  struct wcsprm* wcs;
  
  //pass header into wcs
  int status = wcspih((char *)header.get_cstring(), nkeyrec, WCSHDR_all, 2, &nreject, &nwcs, &wcs);
  
  if(status > 0){
    Rcout << "Failed WCS header read!" << "\n";
  }
  
  NumericMatrix pixel(naxis, ncoord);
  for (i = 0; i < ncoord; i++) {
    pixel(0, i) = x[i];
    pixel(1, i) = y[i];
  }
  NumericVector phi(ncoord);
  NumericVector theta(ncoord);
  NumericMatrix img(naxis, ncoord);
  IntegerVector stat(ncoord);
  NumericMatrix world_matrix(naxis, ncoord);
  
  status = wcsp2s(wcs, ncoord, naxis,
                  &(pixel[0]), &(img[0]), &(phi[0]), &(theta[0]),
                  &(world_matrix[0]), &(stat[0]));
  
  wcsvfree(&nwcs, &wcs);
  if(status > 0 || stat[0] > 0){
    Rcout << "Failed p2s conversion!" << "\n";
    //fprintf(stderr, "ERROR %d from wcsset(): %s.\n", status, wcs_errmsg[status]);
    return(stat);
  }else{
    return(transpose(world_matrix));
  }
}

// [[Rcpp::export]]
SEXP Cwcs_head_s2p(Rcpp::NumericVector RA, Rcpp::NumericVector Dec, Rcpp::String header, 
                   int nkeyrec){
  int i;
  const int ncoord = RA.length();
  const int naxis = 2;
  int nreject, nwcs;
  
  //setup wcs
  struct wcsprm* wcs;
  
  //pass header into wcs
  int status = wcspih((char *)header.get_cstring(), nkeyrec, WCSHDR_all, 2, &nreject, &nwcs, &wcs);
  
  if(status > 0){
    Rcout << "Failed WCS header read!" << "\n";
  }
  
  NumericMatrix world(naxis, ncoord);
  for (i = 0; i < ncoord; i++) {
    world(0, i) = RA[i];
    world(1, i) = Dec[i];
  }
  NumericVector phi(ncoord);
  NumericVector theta(ncoord);
  NumericMatrix img(naxis, ncoord);
  IntegerVector stat(ncoord);
  NumericMatrix pixel_matrix(naxis, ncoord);
  
  status = wcss2p(wcs, ncoord, naxis,
                  &(world[0]), &(phi[0]), &(theta[0]), &(img[0]),
                  &(pixel_matrix[0]), &(stat[0]));
  
  wcsvfree(&nwcs, &wcs);
  if(status > 0 || stat[0] > 0){
    Rcout << "Failed s2p conversion!" << "\n";
    //fprintf(stderr, "ERROR %d from wcsset(): %s.\n", status, wcs_errmsg[status]);
    return(stat);
  }else{
    return(transpose(pixel_matrix));
  }
}
