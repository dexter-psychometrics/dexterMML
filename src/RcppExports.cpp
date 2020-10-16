// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// mat_pre
Rcpp::List mat_pre(arma::imat& dat, const int max_score);
RcppExport SEXP _dexterMML_mat_pre(SEXP datSEXP, SEXP max_scoreSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::imat& >::type dat(datSEXP);
    Rcpp::traits::input_parameter< const int >::type max_score(max_scoreSEXP);
    rcpp_result_gen = Rcpp::wrap(mat_pre(dat, max_score));
    return rcpp_result_gen;
END_RCPP
}
// categorize
arma::imat categorize(const arma::ivec& inp, const arma::ivec& pni, const arma::ivec& icnp, const arma::ivec& pcni, const arma::ivec& ip, const arma::ivec& pi, const arma::imat& icat, const arma::ivec& imax, const int max_cat, arma::ivec& ix, arma::ivec& px);
RcppExport SEXP _dexterMML_categorize(SEXP inpSEXP, SEXP pniSEXP, SEXP icnpSEXP, SEXP pcniSEXP, SEXP ipSEXP, SEXP piSEXP, SEXP icatSEXP, SEXP imaxSEXP, SEXP max_catSEXP, SEXP ixSEXP, SEXP pxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::ivec& >::type inp(inpSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pni(pniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type icnp(icnpSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pcni(pcniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ip(ipSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type icat(icatSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type imax(imaxSEXP);
    Rcpp::traits::input_parameter< const int >::type max_cat(max_catSEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type ix(ixSEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type px(pxSEXP);
    rcpp_result_gen = Rcpp::wrap(categorize(inp, pni, icnp, pcni, ip, pi, icat, imax, max_cat, ix, px));
    return rcpp_result_gen;
END_RCPP
}
// estimate_2pl_dich_multigroup
Rcpp::List estimate_2pl_dich_multigroup(const arma::vec& a_start, const arma::vec& b_start, const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, arma::vec& theta, const arma::vec& mu_start, const arma::vec& sigma_start, const arma::ivec& gn, const arma::ivec& pgroup, const int ref_group);
RcppExport SEXP _dexterMML_estimate_2pl_dich_multigroup(SEXP a_startSEXP, SEXP b_startSEXP, SEXP pniSEXP, SEXP pcniSEXP, SEXP piSEXP, SEXP pxSEXP, SEXP thetaSEXP, SEXP mu_startSEXP, SEXP sigma_startSEXP, SEXP gnSEXP, SEXP pgroupSEXP, SEXP ref_groupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type a_start(a_startSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b_start(b_startSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pni(pniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pcni(pcniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type px(pxSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu_start(mu_startSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma_start(sigma_startSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type gn(gnSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pgroup(pgroupSEXP);
    Rcpp::traits::input_parameter< const int >::type ref_group(ref_groupSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_2pl_dich_multigroup(a_start, b_start, pni, pcni, pi, px, theta, mu_start, sigma_start, gn, pgroup, ref_group));
    return rcpp_result_gen;
END_RCPP
}
// estimate_nrm
Rcpp::List estimate_nrm(arma::imat& a, const arma::mat& b_start, const arma::ivec& ncat, const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, arma::vec& theta, const arma::vec& mu_start, const arma::vec& sigma_start, const arma::ivec& gn, const arma::ivec& pgroup, const int ref_group);
RcppExport SEXP _dexterMML_estimate_nrm(SEXP aSEXP, SEXP b_startSEXP, SEXP ncatSEXP, SEXP pniSEXP, SEXP pcniSEXP, SEXP piSEXP, SEXP pxSEXP, SEXP thetaSEXP, SEXP mu_startSEXP, SEXP sigma_startSEXP, SEXP gnSEXP, SEXP pgroupSEXP, SEXP ref_groupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::imat& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type b_start(b_startSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pni(pniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pcni(pcniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type px(pxSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu_start(mu_startSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma_start(sigma_startSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type gn(gnSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pgroup(pgroupSEXP);
    Rcpp::traits::input_parameter< const int >::type ref_group(ref_groupSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_nrm(a, b_start, ncat, pni, pcni, pi, px, theta, mu_start, sigma_start, gn, pgroup, ref_group));
    return rcpp_result_gen;
END_RCPP
}
// Oakes_2pl_dich
Rcpp::List Oakes_2pl_dich(const arma::vec& a, const arma::vec& b, arma::mat& r0, arma::mat& r1, const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, arma::vec& theta, const arma::vec& mu, const arma::vec& sigma, const arma::ivec& gn, const arma::ivec& pgroup, const int ref_group);
RcppExport SEXP _dexterMML_Oakes_2pl_dich(SEXP aSEXP, SEXP bSEXP, SEXP r0SEXP, SEXP r1SEXP, SEXP pniSEXP, SEXP pcniSEXP, SEXP piSEXP, SEXP pxSEXP, SEXP thetaSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP gnSEXP, SEXP pgroupSEXP, SEXP ref_groupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type r0(r0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type r1(r1SEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pni(pniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pcni(pcniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type px(pxSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type gn(gnSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pgroup(pgroupSEXP);
    Rcpp::traits::input_parameter< const int >::type ref_group(ref_groupSEXP);
    rcpp_result_gen = Rcpp::wrap(Oakes_2pl_dich(a, b, r0, r1, pni, pcni, pi, px, theta, mu, sigma, gn, pgroup, ref_group));
    return rcpp_result_gen;
END_RCPP
}
// prox_dich
Rcpp::List prox_dich(const arma::ivec& isum, const arma::ivec& psum, const arma::ivec& inp, const arma::ivec& pni, const arma::ivec& icnp, const arma::ivec& pcni, const arma::ivec& ip, const arma::ivec& pi, const int max_iter, const double min_change);
RcppExport SEXP _dexterMML_prox_dich(SEXP isumSEXP, SEXP psumSEXP, SEXP inpSEXP, SEXP pniSEXP, SEXP icnpSEXP, SEXP pcniSEXP, SEXP ipSEXP, SEXP piSEXP, SEXP max_iterSEXP, SEXP min_changeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::ivec& >::type isum(isumSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type psum(psumSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type inp(inpSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pni(pniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type icnp(icnpSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pcni(pcniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ip(ipSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< const double >::type min_change(min_changeSEXP);
    rcpp_result_gen = Rcpp::wrap(prox_dich(isum, psum, inp, pni, icnp, pcni, ip, pi, max_iter, min_change));
    return rcpp_result_gen;
END_RCPP
}
// start_lr
Rcpp::List start_lr(const arma::vec& theta, const arma::ivec& ip, const arma::ivec& ix, const arma::ivec& inp, const arma::ivec& icnp, const arma::vec& ibeta);
RcppExport SEXP _dexterMML_start_lr(SEXP thetaSEXP, SEXP ipSEXP, SEXP ixSEXP, SEXP inpSEXP, SEXP icnpSEXP, SEXP ibetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ip(ipSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ix(ixSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type inp(inpSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type icnp(icnpSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ibeta(ibetaSEXP);
    rcpp_result_gen = Rcpp::wrap(start_lr(theta, ip, ix, inp, icnp, ibeta));
    return rcpp_result_gen;
END_RCPP
}
// test_ll
double test_ll(arma::ivec& a, arma::vec& b, arma::vec& theta, arma::mat& r);
RcppExport SEXP _dexterMML_test_ll(SEXP aSEXP, SEXP bSEXP, SEXP thetaSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(test_ll(a, b, theta, r));
    return rcpp_result_gen;
END_RCPP
}
// test_df
arma::vec test_df(arma::ivec& a, arma::vec& b, arma::vec& theta, arma::mat& r);
RcppExport SEXP _dexterMML_test_df(SEXP aSEXP, SEXP bSEXP, SEXP thetaSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(test_df(a, b, theta, r));
    return rcpp_result_gen;
END_RCPP
}
// test_minimize
Rcpp::List test_minimize(arma::ivec& a, arma::vec& b, arma::vec& theta, arma::mat& r);
RcppExport SEXP _dexterMML_test_minimize(SEXP aSEXP, SEXP bSEXP, SEXP thetaSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(test_minimize(a, b, theta, r));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dexterMML_mat_pre", (DL_FUNC) &_dexterMML_mat_pre, 2},
    {"_dexterMML_categorize", (DL_FUNC) &_dexterMML_categorize, 11},
    {"_dexterMML_estimate_2pl_dich_multigroup", (DL_FUNC) &_dexterMML_estimate_2pl_dich_multigroup, 12},
    {"_dexterMML_estimate_nrm", (DL_FUNC) &_dexterMML_estimate_nrm, 13},
    {"_dexterMML_Oakes_2pl_dich", (DL_FUNC) &_dexterMML_Oakes_2pl_dich, 14},
    {"_dexterMML_prox_dich", (DL_FUNC) &_dexterMML_prox_dich, 10},
    {"_dexterMML_start_lr", (DL_FUNC) &_dexterMML_start_lr, 6},
    {"_dexterMML_test_ll", (DL_FUNC) &_dexterMML_test_ll, 4},
    {"_dexterMML_test_df", (DL_FUNC) &_dexterMML_test_df, 4},
    {"_dexterMML_test_minimize", (DL_FUNC) &_dexterMML_test_minimize, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_dexterMML(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
