#include <Rcpp.h>
#include <wcslib.h>

const double tol = 1.0e-10;

const int NELEM = 9;
const int NAXIS = 4;
const double CRPIX[4] =  {  513.0,  0.0,  0.0,  0.0};
const double PC[4][4] = {{    1.1,  0.0,  0.0,  0.0},
                         {    0.0,  1.0,  0.0,  0.1},
                         {    0.0,  0.0,  1.0,  0.0},
                         {    0.0,  0.2,  0.0,  1.0}};
const double CDELT[4] =  {-9.635265432e-6, 1.0, 0.1, -1.0};

char CTYPE[4][9] = {"WAVE-F2W", "XLAT-BON", "TIME-LOG", "XLON-BON"};

const double CRVAL[4] = {0.214982042, -30.0, 1.0, 150.0};
const double LONPOLE  = 150.0;
const double LATPOLE  = 999.0;
const double RESTFRQ  =   1.42040575e9;
const double RESTWAV  =   0.0;

int NPV = 3;
struct pvcard PV[3]; /* Projection parameters are set in main(). */

int check_error(struct wcsprm *wcs, int status, int exstatus, const char *exmsg)
{
  const char *errmsg = (status ? wcs->err->msg : "");

  if (status == exstatus && strcmp(errmsg, exmsg) == 0) {
    wcsperr(wcs, "OK: ");
    wcsprintf("...succeeded.\n");
  } else {
    wcsprintf("Expected error %d: '%s', got\n", exstatus, exmsg);
    wcsperr(wcs, "");
    wcsprintf("...failed.\n");
    return 1;
  }

  return 0;
}


using namespace Rcpp;

// [[Rcpp::export]]
int rcpp_hello_world() {

  int    i, k, lat, lng, nFail1 = 0, nFail2 = 0, stat[361], status;
  double img[361][NELEM], lat1, lng1, phi[361], pixel1[361][NELEM],
         pixel2[361][NELEM], r, resid, residmax, theta[361],
         world1[361][NELEM], world2[361][NELEM];
  struct wcsprm wcs;
  wcs.flag = -1;
  wcsini(1, NAXIS, &wcs);
    
  for (int j = 0; j < NAXIS; j++) {
    wcs.crpix[j] = CRPIX[j];
  }

  auto pcij = wcs.pc;
  for (i = 0; i < NAXIS; i++) {
    for (int j = 0; j < NAXIS; j++) {
      *(pcij++) = PC[i][j];
    }
  }

  for (i = 0; i < NAXIS; i++) {
    wcs.cdelt[i] = CDELT[i];
  }

  for (i = 0; i < NAXIS; i++) {
    strcpy(wcs.ctype[i], &CTYPE[i][0]);
  }

  for (i = 0; i < NAXIS; i++) {
    wcs.crval[i] = CRVAL[i];
  }

  wcs.lonpole = LONPOLE;
  wcs.latpole = LATPOLE;

  wcs.restfrq = RESTFRQ;
  wcs.restwav = RESTWAV;

  wcs.npv = NPV;
  for (i = 0; i < NPV; i++) {
    wcs.pv[i] = PV[i];
  }

  for (lat = 90; lat >= -90; lat--) {
    lat1 = (double)lat;

    for (lng = -180, k = 0; lng <= 180; lng++, k++) {
      lng1 = (double)lng;

      world1[k][wcs.lng] = lng1;
      world1[k][wcs.lat] = lat1;
    }

    if (wcss2p(&wcs, 361, NELEM, world1[0], phi, theta, img[0], pixel1[0],
               stat)) {
      printf("  At wcss2p#1 with lat1 == %f\n", lat1);
      wcsperr(&wcs, "  ");
      continue;
    }

    if (wcsp2s(&wcs, 361, NELEM, pixel1[0], img[0], phi, theta, world2[0],
               stat)) {
      printf("  At wcsp2s with lat1 == %f\n", lat1);
      wcsperr(&wcs, "  ");
      continue;
    }

    if (wcss2p(&wcs, 361, NELEM, world2[0], phi, theta, img[0], pixel2[0],
               stat)) {
      printf("  At wcss2p#2 with lat1 == %f\n", lat1);
      wcsperr(&wcs, "  ");
      continue;
    }

    for (k = 0; k < 361; k++) {
      resid = 0.0;
      for (i = 0; i < NAXIS; i++) {
        r = pixel2[k][i] - pixel1[k][i];
        resid += r*r;
      }

      resid = sqrt(resid);
      if (resid > residmax) residmax = resid;

      if (resid > tol) {
        nFail1++;
        printf("\nClosure error:\n"
               "world1:%18.12f%18.12f%18.12f%18.12f\n"
               "pixel1:%18.12f%18.12f%18.12f%18.12f\n"
               "world2:%18.12f%18.12f%18.12f%18.12f\n"
               "pixel2:%18.12f%18.12f%18.12f%18.12f\n",
          world1[k][0], world1[k][1], world1[k][2], world1[k][3],
          pixel1[k][0], pixel1[k][1], pixel1[k][2], pixel1[k][3],
          world2[k][0], world2[k][1], world2[k][2], world2[k][3],
          pixel2[k][0], pixel2[k][1], pixel2[k][2], pixel2[k][3]);
       }
    }
  }
  printf("wcsp2s/wcss2p: Maximum closure residual = %.1e pixel.\n", residmax);


  /* Test wcserr and wcsprintf() as well. */
  nFail2 = 0;
  wcsprintf_set(stdout);
  wcsprintf("\n\nIGNORE messages marked with 'OK', they test wcserr "
    "(and wcsprintf):\n");

  wcserr_enable(1);

  /* Test 1. */
  wcs.pv[2].value = UNDEFINED;
  status = wcsset(&wcs);
  nFail2 += check_error(&wcs, status, WCSERR_BAD_PARAM,
                        "Invalid parameter value");

  if (nFail1 || nFail2) {
    if (nFail1) {
      printf("\nFAIL: %d closure residuals exceed reporting tolerance.\n",
        nFail1);
    }

    if (nFail2) {
      printf("FAIL: %d error messages differ from that expected.\n", nFail2);
    }
  } else {
    printf("\nPASS: All closure residuals are within reporting tolerance.\n");
    printf("PASS: All error messages reported as expected.\n");
  }


  /* Clean up. */
  wcsfree(&wcs);

  return nFail1 + nFail2;
}
