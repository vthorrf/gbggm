// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// dpoly
Rcpp::List dpoly(const arma::mat& X, arma::mat R, std::vector<std::vector<double>> taus);
RcppExport SEXP _gbggm_dpoly(SEXP XSEXP, SEXP RSEXP, SEXP tausSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type R(RSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<double>> >::type taus(tausSEXP);
    rcpp_result_gen = Rcpp::wrap(dpoly(X, R, taus));
    return rcpp_result_gen;
END_RCPP
}
// euclidean
double euclidean(NumericVector x);
RcppExport SEXP _gbggm_euclidean(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(euclidean(x));
    return rcpp_result_gen;
END_RCPP
}
// pcor
arma::mat pcor(arma::mat R);
RcppExport SEXP _gbggm_pcor(SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(pcor(R));
    return rcpp_result_gen;
END_RCPP
}
// Modularity
double Modularity(NumericMatrix adjMatrix, IntegerVector community);
RcppExport SEXP _gbggm_Modularity(SEXP adjMatrixSEXP, SEXP communitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type adjMatrix(adjMatrixSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type community(communitySEXP);
    rcpp_result_gen = Rcpp::wrap(Modularity(adjMatrix, community));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gbggm_dpoly", (DL_FUNC) &_gbggm_dpoly, 3},
    {"_gbggm_euclidean", (DL_FUNC) &_gbggm_euclidean, 1},
    {"_gbggm_pcor", (DL_FUNC) &_gbggm_pcor, 1},
    {"_gbggm_Modularity", (DL_FUNC) &_gbggm_Modularity, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_gbggm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}