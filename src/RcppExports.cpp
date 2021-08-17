// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Cwcs_s2p
SEXP Cwcs_s2p(Rcpp::NumericVector RA, Rcpp::NumericVector Dec, Rcpp::String CTYPE1, Rcpp::String CTYPE2, double CRVAL1, double CRVAL2, double CRPIX1, double CRPIX2, double CD1_1, double CD1_2, double CD2_1, double CD2_2, Rcpp::String RADESYS, int EQUINOX, double PV1_1, double PV1_2, double PV1_3, double PV2_1, double PV2_2, double PV2_3, double PV2_4, double PV2_5);
RcppExport SEXP _Rwcs_Cwcs_s2p(SEXP RASEXP, SEXP DecSEXP, SEXP CTYPE1SEXP, SEXP CTYPE2SEXP, SEXP CRVAL1SEXP, SEXP CRVAL2SEXP, SEXP CRPIX1SEXP, SEXP CRPIX2SEXP, SEXP CD1_1SEXP, SEXP CD1_2SEXP, SEXP CD2_1SEXP, SEXP CD2_2SEXP, SEXP RADESYSSEXP, SEXP EQUINOXSEXP, SEXP PV1_1SEXP, SEXP PV1_2SEXP, SEXP PV1_3SEXP, SEXP PV2_1SEXP, SEXP PV2_2SEXP, SEXP PV2_3SEXP, SEXP PV2_4SEXP, SEXP PV2_5SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type RA(RASEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Dec(DecSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type CTYPE1(CTYPE1SEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type CTYPE2(CTYPE2SEXP);
    Rcpp::traits::input_parameter< double >::type CRVAL1(CRVAL1SEXP);
    Rcpp::traits::input_parameter< double >::type CRVAL2(CRVAL2SEXP);
    Rcpp::traits::input_parameter< double >::type CRPIX1(CRPIX1SEXP);
    Rcpp::traits::input_parameter< double >::type CRPIX2(CRPIX2SEXP);
    Rcpp::traits::input_parameter< double >::type CD1_1(CD1_1SEXP);
    Rcpp::traits::input_parameter< double >::type CD1_2(CD1_2SEXP);
    Rcpp::traits::input_parameter< double >::type CD2_1(CD2_1SEXP);
    Rcpp::traits::input_parameter< double >::type CD2_2(CD2_2SEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type RADESYS(RADESYSSEXP);
    Rcpp::traits::input_parameter< int >::type EQUINOX(EQUINOXSEXP);
    Rcpp::traits::input_parameter< double >::type PV1_1(PV1_1SEXP);
    Rcpp::traits::input_parameter< double >::type PV1_2(PV1_2SEXP);
    Rcpp::traits::input_parameter< double >::type PV1_3(PV1_3SEXP);
    Rcpp::traits::input_parameter< double >::type PV2_1(PV2_1SEXP);
    Rcpp::traits::input_parameter< double >::type PV2_2(PV2_2SEXP);
    Rcpp::traits::input_parameter< double >::type PV2_3(PV2_3SEXP);
    Rcpp::traits::input_parameter< double >::type PV2_4(PV2_4SEXP);
    Rcpp::traits::input_parameter< double >::type PV2_5(PV2_5SEXP);
    rcpp_result_gen = Rcpp::wrap(Cwcs_s2p(RA, Dec, CTYPE1, CTYPE2, CRVAL1, CRVAL2, CRPIX1, CRPIX2, CD1_1, CD1_2, CD2_1, CD2_2, RADESYS, EQUINOX, PV1_1, PV1_2, PV1_3, PV2_1, PV2_2, PV2_3, PV2_4, PV2_5));
    return rcpp_result_gen;
END_RCPP
}
// Cwcs_p2s
SEXP Cwcs_p2s(Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::String CTYPE1, Rcpp::String CTYPE2, double CRVAL1, double CRVAL2, double CRPIX1, double CRPIX2, double CD1_1, double CD1_2, double CD2_1, double CD2_2, Rcpp::String RADESYS, int EQUINOX, double PV1_1, double PV1_2, double PV1_3, double PV2_1, double PV2_2, double PV2_3, double PV2_4, double PV2_5);
RcppExport SEXP _Rwcs_Cwcs_p2s(SEXP xSEXP, SEXP ySEXP, SEXP CTYPE1SEXP, SEXP CTYPE2SEXP, SEXP CRVAL1SEXP, SEXP CRVAL2SEXP, SEXP CRPIX1SEXP, SEXP CRPIX2SEXP, SEXP CD1_1SEXP, SEXP CD1_2SEXP, SEXP CD2_1SEXP, SEXP CD2_2SEXP, SEXP RADESYSSEXP, SEXP EQUINOXSEXP, SEXP PV1_1SEXP, SEXP PV1_2SEXP, SEXP PV1_3SEXP, SEXP PV2_1SEXP, SEXP PV2_2SEXP, SEXP PV2_3SEXP, SEXP PV2_4SEXP, SEXP PV2_5SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type CTYPE1(CTYPE1SEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type CTYPE2(CTYPE2SEXP);
    Rcpp::traits::input_parameter< double >::type CRVAL1(CRVAL1SEXP);
    Rcpp::traits::input_parameter< double >::type CRVAL2(CRVAL2SEXP);
    Rcpp::traits::input_parameter< double >::type CRPIX1(CRPIX1SEXP);
    Rcpp::traits::input_parameter< double >::type CRPIX2(CRPIX2SEXP);
    Rcpp::traits::input_parameter< double >::type CD1_1(CD1_1SEXP);
    Rcpp::traits::input_parameter< double >::type CD1_2(CD1_2SEXP);
    Rcpp::traits::input_parameter< double >::type CD2_1(CD2_1SEXP);
    Rcpp::traits::input_parameter< double >::type CD2_2(CD2_2SEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type RADESYS(RADESYSSEXP);
    Rcpp::traits::input_parameter< int >::type EQUINOX(EQUINOXSEXP);
    Rcpp::traits::input_parameter< double >::type PV1_1(PV1_1SEXP);
    Rcpp::traits::input_parameter< double >::type PV1_2(PV1_2SEXP);
    Rcpp::traits::input_parameter< double >::type PV1_3(PV1_3SEXP);
    Rcpp::traits::input_parameter< double >::type PV2_1(PV2_1SEXP);
    Rcpp::traits::input_parameter< double >::type PV2_2(PV2_2SEXP);
    Rcpp::traits::input_parameter< double >::type PV2_3(PV2_3SEXP);
    Rcpp::traits::input_parameter< double >::type PV2_4(PV2_4SEXP);
    Rcpp::traits::input_parameter< double >::type PV2_5(PV2_5SEXP);
    rcpp_result_gen = Rcpp::wrap(Cwcs_p2s(x, y, CTYPE1, CTYPE2, CRVAL1, CRVAL2, CRPIX1, CRPIX2, CD1_1, CD1_2, CD2_1, CD2_2, RADESYS, EQUINOX, PV1_1, PV1_2, PV1_3, PV2_1, PV2_2, PV2_3, PV2_4, PV2_5));
    return rcpp_result_gen;
END_RCPP
}
// Cwcs_head_p2s
SEXP Cwcs_head_p2s(Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::String header, int nkeyrec);
RcppExport SEXP _Rwcs_Cwcs_head_p2s(SEXP xSEXP, SEXP ySEXP, SEXP headerSEXP, SEXP nkeyrecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type header(headerSEXP);
    Rcpp::traits::input_parameter< int >::type nkeyrec(nkeyrecSEXP);
    rcpp_result_gen = Rcpp::wrap(Cwcs_head_p2s(x, y, header, nkeyrec));
    return rcpp_result_gen;
END_RCPP
}
// Cwcs_head_s2p
SEXP Cwcs_head_s2p(Rcpp::NumericVector RA, Rcpp::NumericVector Dec, Rcpp::String header, int nkeyrec);
RcppExport SEXP _Rwcs_Cwcs_head_s2p(SEXP RASEXP, SEXP DecSEXP, SEXP headerSEXP, SEXP nkeyrecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type RA(RASEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Dec(DecSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type header(headerSEXP);
    Rcpp::traits::input_parameter< int >::type nkeyrec(nkeyrecSEXP);
    rcpp_result_gen = Rcpp::wrap(Cwcs_head_s2p(RA, Dec, header, nkeyrec));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Rwcs_Cwcs_s2p", (DL_FUNC) &_Rwcs_Cwcs_s2p, 22},
    {"_Rwcs_Cwcs_p2s", (DL_FUNC) &_Rwcs_Cwcs_p2s, 22},
    {"_Rwcs_Cwcs_head_p2s", (DL_FUNC) &_Rwcs_Cwcs_head_p2s, 4},
    {"_Rwcs_Cwcs_head_s2p", (DL_FUNC) &_Rwcs_Cwcs_head_s2p, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_Rwcs(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
