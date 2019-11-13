#include <Rcpp.h>
#include <algorithm>
#include <utility>
#include <vector>
#include <wcslib.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP Cwcs_s2p(Rcpp::NumericVector RA, Rcpp::NumericVector Dec,
              Rcpp::String CTYPE1, Rcpp::String CTYPE2,
              double CRVAL1 = 0, double CRVAL2 = 0,
              double CRPIX1 = 0, double CRPIX2 = 0,
              double CD1_1 = 1, double CD1_2 = 0,
              double CD2_1 = 0, double CD2_2 = 1,
              double PV1 = 0, double PV2 = 0
              ){

  int ncoord = RA.length();
  int nelem = 2;
  int naxis = 2;

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

  //insert wcs pv
  wcs.npv = 2;
  wcs.pv[0].value = PV1;
  wcs.pv[1].value = PV2;

  NumericMatrix world(nelem, ncoord);
  for (int i = 0; i < ncoord; i++) {
    world(0, i) = RA[i];
    world(1, i) = Dec[i];
  }
  NumericVector phi(ncoord);
  NumericVector theta(ncoord);
  NumericMatrix img(nelem, ncoord);
  IntegerVector stat(ncoord);
  NumericMatrix pixel(nelem, ncoord);

  wcss2p(&wcs, ncoord, nelem,
         &(world[0]), &(phi[0]), &(theta[0]), &(img[0]),
         &(pixel[0]), &(stat[0]));

  // TODO: check stat

  return(pixel);
}

// [[Rcpp::export]]
SEXP Cwcs_p2s(Rcpp::NumericVector x, Rcpp::NumericVector y,
              Rcpp::String CTYPE1, Rcpp::String CTYPE2,
              double CRVAL1 = 0, double CRVAL2 = 0,
              double CRPIX1 = 0, double CRPIX2 = 0,
              double CD1_1 = 1, double CD1_2 = 0,
              double CD2_1 = 0, double CD2_2 = 1,
              double PV1 = 0, double PV2 = 0
){
  
  int i,j;
  
  int ncoord = x.length();
  int nelem = 2;
  int naxis = 2;
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
  
  //insert wcs pv
  wcs.npv = 2;
  wcs.pv[0].value = PV1;
  wcs.pv[1].value = PV2;
  
  NumericMatrix pixel(nelem, ncoord);
  for (int i = 0; i < ncoord; i++) {
    pixel(0, i) = x[i];
    pixel(1, i) = y[i];
  }
  NumericVector phi(ncoord);
  NumericVector theta(ncoord);
  NumericMatrix img(nelem, ncoord);
  IntegerVector stat(ncoord);
  NumericMatrix world(nelem, ncoord);
  
  wcsp2s(&wcs, ncoord, nelem,
         &(pixel[0]), &(phi[0]), &(theta[0]), &(img[0]),
         &(world[0]), &(stat[0]));
  
  return(world);
}
