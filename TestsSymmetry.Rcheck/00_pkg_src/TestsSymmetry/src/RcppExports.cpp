// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// getV
double getV(NumericVector x_, Nullable<NumericVector> y_);
RcppExport SEXP _TestsSymmetry_getV(SEXP x_SEXP, SEXP y_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x_(x_SEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type y_(y_SEXP);
    rcpp_result_gen = Rcpp::wrap(getV(x_, y_));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TestsSymmetry_getV", (DL_FUNC) &_TestsSymmetry_getV, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_TestsSymmetry(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}