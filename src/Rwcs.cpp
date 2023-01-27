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

static const int naxis = 2;

static void enable_wcsperr()
{
  wcserr_enable(1);
  wcsprintf_set(nullptr);
}

static SEXP _wcss2p(struct wcsprm *wcs, NumericVector RA, NumericVector Dec)
{
  const int ncoord = RA.length();
  
  NumericMatrix world(naxis, ncoord);
  for (int i = 0; i < ncoord; i++) {
    world(0, i) = RA[i];
    world(1, i) = Dec[i];
  }
  NumericVector phi(ncoord);
  NumericVector theta(ncoord);
  NumericMatrix img(naxis, ncoord);
  IntegerVector stat(ncoord);
  NumericMatrix pixel_matrix(naxis, ncoord);
  auto status = wcss2p(wcs, ncoord, naxis,
                       &(world[0]), &(phi[0]), &(theta[0]), &(img[0]),
                       &(pixel_matrix[0]), &(stat[0]));

  if (status) {
    wcsperr(wcs, "");
    Rcerr << "Failed s2p conversion :(:\n" << wcsprintf_buf();
    return stat;
  }
  return transpose(pixel_matrix);
}

static SEXP _wcsp2s(struct wcsprm *wcs, NumericVector x, NumericVector y)
{
  const int ncoord = x.length();
  
  NumericMatrix pixel(naxis, ncoord);
  for (int i = 0; i < ncoord; i++) {
    pixel(0, i) = x[i];
    pixel(1, i) = y[i];
  }
  NumericVector phi(ncoord);
  NumericVector theta(ncoord);
  NumericMatrix img(naxis, ncoord);
  IntegerVector stat(ncoord);
  NumericMatrix world_matrix(naxis, ncoord);
  
  auto status = wcsp2s(wcs, ncoord, naxis,
                       &(pixel[0]), &(img[0]), &(phi[0]), &(theta[0]),
                       &(world_matrix[0]), &(stat[0]));

  if (status) {
    wcsperr(wcs, "");
    Rcerr << "Failed p2s conversion :(:\n" << wcsprintf_buf();
    return stat;
  }
  return transpose(world_matrix);
}

static void _wcsset(struct wcsprm* wcs,
                    Rcpp::String CTYPE1, Rcpp::String CTYPE2,
                    double CRVAL1, double CRVAL2, double CRPIX1, double CRPIX2,
                    double CD1_1, double CD1_2, double CD2_1, double CD2_2,
                    Rcpp::String RADESYS, int EQUINOX,
                    double PV1_0, double PV1_1, double PV1_2, double PV1_3, double PV1_4,
                    // double PV1_5, double PV1_6, double PV1_7, double PV1_8, double PV1_9, double PV1_10,
                    double PV2_0, double PV2_1, double PV2_2, double PV2_3, double PV2_4, double PV2_5
                    // double PV2_6, double PV2_7, double PV2_8, double PV2_9, double PV2_10
                    )
                    
{
  //setup wcs
  wcs->flag = -1;
  wcsini(1, naxis, wcs);

  //insert wcs val
  wcs->crval[0] = CRVAL1;
  wcs->crval[1] = CRVAL2;

  //insert wcs pix
  wcs->crpix[0] = CRPIX1;
  wcs->crpix[1] = CRPIX2;

  //insert wcs pix
  wcs->crpix[0] = CRPIX1;
  wcs->crpix[1] = CRPIX2;

  //insert wcs cd matrix
  wcs->pc[0] = CD1_1;
  wcs->pc[1] = CD1_2;
  wcs->pc[2] = CD2_1;
  wcs->pc[3] = CD2_2;

  //insert ctype
  strcpy(wcs->ctype[0], CTYPE1.get_cstring());
  strcpy(wcs->ctype[1], CTYPE2.get_cstring());

  //insert radesys and equinox
  strcpy(wcs->radesys, RADESYS.get_cstring());
  wcs->equinox = EQUINOX;

  //insert wcs pv
#define FILL_PV(WHICH, I, M) \
  if (!R_IsNA(WHICH)) {               \
    wcs->pv[wcs->npv].i = I;          \
    wcs->pv[wcs->npv].m = M;          \
    wcs->pv[wcs->npv].value = WHICH;  \
    wcs->npv++;                       \
  }

  wcs->npv = 0;
  FILL_PV(PV1_0, 1, 0);
  FILL_PV(PV1_1, 1, 1);
  FILL_PV(PV1_2, 1, 2);
  FILL_PV(PV1_3, 1, 3);
  FILL_PV(PV1_4, 1, 4);
  // FILL_PV(PV1_5, 1, 5);
  // FILL_PV(PV1_6, 1, 6);
  // FILL_PV(PV1_7, 1, 7);
  // FILL_PV(PV1_8, 1, 8);
  // FILL_PV(PV1_9, 1, 9);
  // FILL_PV(PV1_10, 1, 10);
  FILL_PV(PV2_0, 2, 0);
  FILL_PV(PV2_1, 2, 1);
  FILL_PV(PV2_2, 2, 2);
  FILL_PV(PV2_3, 2, 3);
  FILL_PV(PV2_4, 2, 4);
  FILL_PV(PV2_5, 2, 5);
  // FILL_PV(PV2_6, 2, 6);
  // FILL_PV(PV2_7, 2, 7);
  // FILL_PV(PV2_8, 2, 8);
  // FILL_PV(PV2_9, 2, 9);
  // FILL_PV(PV2_10, 2, 10);

  //set long and lat axis positions
  wcs->lng = 0;
  wcs->lat = 1;

  wcsset(wcs);
}

// [[Rcpp::export]]
SEXP Cwcs_s2p(Rcpp::NumericVector RA, Rcpp::NumericVector Dec,
              Rcpp::String CTYPE1 = "RA---TAN", Rcpp::String CTYPE2 = "DEC--TAN",
              double CRVAL1 = 0, double CRVAL2 = 0,
              double CRPIX1 = 0, double CRPIX2 = 0,
              double CD1_1 = 1, double CD1_2 = 0,
              double CD2_1 = 0, double CD2_2 = 1,
              Rcpp::String RADESYS = "ICRS", int EQUINOX = 2000,
              double PV1_0 = NA_REAL, double PV1_1 = NA_REAL, double PV1_2 = NA_REAL, double PV1_3 = NA_REAL, double PV1_4 = NA_REAL,
              double PV2_0 = NA_REAL, double PV2_1 = NA_REAL, double PV2_2 = NA_REAL, double PV2_3 = NA_REAL, double PV2_4 = NA_REAL, double PV2_5 = NA_REAL
              )
{
  enable_wcsperr();
  struct wcsprm wcs;
  _wcsset(&wcs, CTYPE1, CTYPE2, CRVAL1, CRVAL2, CRPIX1, CRPIX2, CD1_1, CD1_2, CD2_1, CD2_2, RADESYS, EQUINOX,
    PV1_0, PV1_1, PV1_2, PV1_3, PV1_4,
    // PV1_5, PV1_6, PV1_7, PV1_8, PV1_9, PV1_10,
    PV2_0, PV2_1, PV2_2, PV2_3, PV2_4, PV2_5
    // PV2_6, PV2_7, PV2_8, PV2_9, PV2_10
    );
  auto result = _wcss2p(&wcs, RA, Dec);
  wcsfree(&wcs);
  return result;
}

// [[Rcpp::export]]
SEXP Cwcs_p2s(Rcpp::NumericVector x, Rcpp::NumericVector y,
              Rcpp::String CTYPE1 = "RA---TAN", Rcpp::String CTYPE2 = "DEC--TAN",
              double CRVAL1 = 0, double CRVAL2 = 0,
              double CRPIX1 = 0, double CRPIX2 = 0,
              double CD1_1 = 1, double CD1_2 = 0,
              double CD2_1 = 0, double CD2_2 = 1,
              Rcpp::String RADESYS = "ICRS", int EQUINOX = 2000,
              double PV1_0 = NA_REAL, double PV1_1 = NA_REAL, double PV1_2 = NA_REAL, double PV1_3 = NA_REAL, double PV1_4 = NA_REAL,
              double PV2_0 = NA_REAL, double PV2_1 = NA_REAL, double PV2_2 = NA_REAL, double PV2_3 = NA_REAL, double PV2_4 = NA_REAL, double PV2_5 = NA_REAL
              )
{
  enable_wcsperr();
  struct wcsprm wcs;
  _wcsset(&wcs, CTYPE1, CTYPE2, CRVAL1, CRVAL2, CRPIX1, CRPIX2, CD1_1, CD1_2, CD2_1, CD2_2, RADESYS, EQUINOX,
    PV1_0, PV1_1, PV1_2, PV1_3, PV1_4, 
    //PV1_5, PV1_6, PV1_7, PV1_8, PV1_9, PV1_10,
    PV2_0, PV2_1, PV2_2, PV2_3, PV2_4, PV2_5
    //PV2_6, PV2_7, PV2_8, PV2_9, PV2_10
    );
  auto result = _wcsp2s(&wcs, x, y);
  wcsfree(&wcs);
  return result;
}

struct wcsprm* _read_from_header(int *nwcs, struct wcsprm** wcs, Rcpp::String header, int nkey, int WCSref, int ctrl)
{
  int nreject;
  int status = wcspih((char *)header.get_cstring(), nkey, WCSHDR_all, ctrl, &nreject, nwcs, wcs);
  
  if (status) {
    Rcerr << "Failed WCS header read :(\n";
    Rcerr << "ERROR " << status << " from wcspih(): " << wcs_errmsg[status] << '\n';
    return nullptr;
  }
  
  int alts[27]{};
  status = wcsidx(*nwcs, wcs, alts);
  if (status) {
    Rcerr << "ERROR " << status << " from wcsidx()(\n";
    return nullptr;
  }
  
  if (alts[WCSref] < 0) {
    Rcout << "Bad WCS projection selection!" << "\n";
    return nullptr;
  }
  
  return wcs[alts[WCSref]];
}

// [[Rcpp::export]]
SEXP Cwcs_head_p2s(Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::String header, 
                   int nkey, int WCSref=0, int ctrl=2)
{
  enable_wcsperr();
  
  int nwcs;
  struct wcsprm* wcs;
  auto wcs_at_ref = _read_from_header(&nwcs, &wcs, header, nkey, WCSref, ctrl);
  if (!wcs_at_ref) {
    wcsvfree(&nwcs, &wcs);
    return nullptr;
  }
  
  auto result = _wcsp2s(wcs_at_ref, x, y);
  wcsvfree(&nwcs, &wcs);
  return result;
}

// [[Rcpp::export]]
SEXP Cwcs_head_s2p(Rcpp::NumericVector RA, Rcpp::NumericVector Dec, Rcpp::String header, 
                   int nkey, int WCSref=0, int ctrl=2)
{
  enable_wcsperr();
  
  int nwcs;
  struct wcsprm* wcs;
  auto wcs_at_ref = _read_from_header(&nwcs, &wcs, header, nkey, WCSref, ctrl);
  if (!wcs_at_ref) {
    wcsvfree(&nwcs, &wcs);
    return nullptr;
  }
  
  auto result = _wcss2p(wcs_at_ref, RA, Dec);
  wcsvfree(&nwcs, &wcs);
  return result;
}
