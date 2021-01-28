// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// theta_2pl
arma::vec theta_2pl(const arma::imat& a, const arma::vec& A, const arma::mat& b, const arma::ivec& ncat, const arma::ivec& pni, const arma::ivec& pcni, arma::ivec& pi, const arma::ivec& px, const bool WLE, const bool USE_A);
RcppExport SEXP _dexterMML_theta_2pl(SEXP aSEXP, SEXP ASEXP, SEXP bSEXP, SEXP ncatSEXP, SEXP pniSEXP, SEXP pcniSEXP, SEXP piSEXP, SEXP pxSEXP, SEXP WLESEXP, SEXP USE_ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::imat& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pni(pniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pcni(pcniSEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type px(pxSEXP);
    Rcpp::traits::input_parameter< const bool >::type WLE(WLESEXP);
    Rcpp::traits::input_parameter< const bool >::type USE_A(USE_ASEXP);
    rcpp_result_gen = Rcpp::wrap(theta_2pl(a, A, b, ncat, pni, pcni, pi, px, WLE, USE_A));
    return rcpp_result_gen;
END_RCPP
}
// E_score
arma::vec E_score(const arma::vec& theta, const arma::vec& A, const arma::imat& a, const arma::mat& b, const arma::ivec& items, const arma::ivec& ncat);
RcppExport SEXP _dexterMML_E_score(SEXP thetaSEXP, SEXP ASEXP, SEXP aSEXP, SEXP bSEXP, SEXP itemsSEXP, SEXP ncatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type items(itemsSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ncat(ncatSEXP);
    rcpp_result_gen = Rcpp::wrap(E_score(theta, A, a, b, items, ncat));
    return rcpp_result_gen;
END_RCPP
}
// mat_pre
Rcpp::List mat_pre(const arma::imat& dat, const int max_score);
RcppExport SEXP _dexterMML_mat_pre(SEXP datSEXP, SEXP max_scoreSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::imat& >::type dat(datSEXP);
    Rcpp::traits::input_parameter< const int >::type max_score(max_scoreSEXP);
    rcpp_result_gen = Rcpp::wrap(mat_pre(dat, max_score));
    return rcpp_result_gen;
END_RCPP
}
// categorize
arma::imat categorize(const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::imat& icat, const arma::ivec& imax, const int max_cat, arma::ivec& px);
RcppExport SEXP _dexterMML_categorize(SEXP pniSEXP, SEXP pcniSEXP, SEXP piSEXP, SEXP icatSEXP, SEXP imaxSEXP, SEXP max_catSEXP, SEXP pxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::ivec& >::type pni(pniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pcni(pcniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type icat(icatSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type imax(imaxSEXP);
    Rcpp::traits::input_parameter< const int >::type max_cat(max_catSEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type px(pxSEXP);
    rcpp_result_gen = Rcpp::wrap(categorize(pni, pcni, pi, icat, imax, max_cat, px));
    return rcpp_result_gen;
END_RCPP
}
// plot_data
Rcpp::DataFrame plot_data(const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, const arma::ivec& inp, const arma::vec& theta, const arma::imat& a, const int item);
RcppExport SEXP _dexterMML_plot_data(SEXP pcniSEXP, SEXP piSEXP, SEXP pxSEXP, SEXP inpSEXP, SEXP thetaSEXP, SEXP aSEXP, SEXP itemSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::ivec& >::type pcni(pcniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type px(pxSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type inp(inpSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const int >::type item(itemSEXP);
    rcpp_result_gen = Rcpp::wrap(plot_data(pcni, pi, px, inp, theta, a, item));
    return rcpp_result_gen;
END_RCPP
}
// design_matrices
Rcpp::List design_matrices(const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& pg, const int nit, const int ng);
RcppExport SEXP _dexterMML_design_matrices(SEXP pniSEXP, SEXP pcniSEXP, SEXP piSEXP, SEXP pgSEXP, SEXP nitSEXP, SEXP ngSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::ivec& >::type pni(pniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pcni(pcniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pg(pgSEXP);
    Rcpp::traits::input_parameter< const int >::type nit(nitSEXP);
    Rcpp::traits::input_parameter< const int >::type ng(ngSEXP);
    rcpp_result_gen = Rcpp::wrap(design_matrices(pni, pcni, pi, pg, nit, ng));
    return rcpp_result_gen;
END_RCPP
}
// check_connected_c
int check_connected_c(const arma::imat& item, const arma::imat& group, const arma::ivec& item_fixed);
RcppExport SEXP _dexterMML_check_connected_c(SEXP itemSEXP, SEXP groupSEXP, SEXP item_fixedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::imat& >::type item(itemSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type group(groupSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type item_fixed(item_fixedSEXP);
    rcpp_result_gen = Rcpp::wrap(check_connected_c(item, group, item_fixed));
    return rcpp_result_gen;
END_RCPP
}
// estimate_nrm
Rcpp::List estimate_nrm(arma::imat& a, const arma::mat& b_start, const arma::ivec& ncat, const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, arma::vec& theta, const arma::vec& mu_start, const arma::vec& sigma_start, const arma::ivec& gn, const arma::ivec& pgroup, const arma::ivec& item_fixed, const int ref_group, const int max_iter, const int pgw);
RcppExport SEXP _dexterMML_estimate_nrm(SEXP aSEXP, SEXP b_startSEXP, SEXP ncatSEXP, SEXP pniSEXP, SEXP pcniSEXP, SEXP piSEXP, SEXP pxSEXP, SEXP thetaSEXP, SEXP mu_startSEXP, SEXP sigma_startSEXP, SEXP gnSEXP, SEXP pgroupSEXP, SEXP item_fixedSEXP, SEXP ref_groupSEXP, SEXP max_iterSEXP, SEXP pgwSEXP) {
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
    Rcpp::traits::input_parameter< const arma::ivec& >::type item_fixed(item_fixedSEXP);
    Rcpp::traits::input_parameter< const int >::type ref_group(ref_groupSEXP);
    Rcpp::traits::input_parameter< const int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< const int >::type pgw(pgwSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_nrm(a, b_start, ncat, pni, pcni, pi, px, theta, mu_start, sigma_start, gn, pgroup, item_fixed, ref_group, max_iter, pgw));
    return rcpp_result_gen;
END_RCPP
}
// test_ll_nrm
double test_ll_nrm(arma::ivec& a, arma::vec theta, arma::mat& r, const arma::vec& par);
RcppExport SEXP _dexterMML_test_ll_nrm(SEXP aSEXP, SEXP thetaSEXP, SEXP rSEXP, SEXP parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type r(rSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type par(parSEXP);
    rcpp_result_gen = Rcpp::wrap(test_ll_nrm(a, theta, r, par));
    return rcpp_result_gen;
END_RCPP
}
// test_gradient_nrm
arma::vec test_gradient_nrm(arma::ivec& a, arma::vec theta, arma::mat& r, const arma::vec& par);
RcppExport SEXP _dexterMML_test_gradient_nrm(SEXP aSEXP, SEXP thetaSEXP, SEXP rSEXP, SEXP parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type r(rSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type par(parSEXP);
    rcpp_result_gen = Rcpp::wrap(test_gradient_nrm(a, theta, r, par));
    return rcpp_result_gen;
END_RCPP
}
// test_hess_nrm
arma::mat test_hess_nrm(arma::ivec& a, arma::vec theta, arma::mat& r, const arma::vec& par);
RcppExport SEXP _dexterMML_test_hess_nrm(SEXP aSEXP, SEXP thetaSEXP, SEXP rSEXP, SEXP parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type r(rSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type par(parSEXP);
    rcpp_result_gen = Rcpp::wrap(test_hess_nrm(a, theta, r, par));
    return rcpp_result_gen;
END_RCPP
}
// Oakes_nrm
Rcpp::List Oakes_nrm(arma::imat& a, const arma::mat& b, const arma::ivec& ncat, arma::field<arma::mat>& r, const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, arma::vec& theta, const arma::vec& mu, const arma::vec& sigma, const arma::ivec& gn, const arma::ivec& pgroup, const arma::imat& dsg_ii, const arma::imat& dsg_gi, const arma::ivec& item_fixed, const int ref_group, const int pgw);
RcppExport SEXP _dexterMML_Oakes_nrm(SEXP aSEXP, SEXP bSEXP, SEXP ncatSEXP, SEXP rSEXP, SEXP pniSEXP, SEXP pcniSEXP, SEXP piSEXP, SEXP pxSEXP, SEXP thetaSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP gnSEXP, SEXP pgroupSEXP, SEXP dsg_iiSEXP, SEXP dsg_giSEXP, SEXP item_fixedSEXP, SEXP ref_groupSEXP, SEXP pgwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::imat& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::mat>& >::type r(rSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pni(pniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pcni(pcniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type px(pxSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type gn(gnSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pgroup(pgroupSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type dsg_ii(dsg_iiSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type dsg_gi(dsg_giSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type item_fixed(item_fixedSEXP);
    Rcpp::traits::input_parameter< const int >::type ref_group(ref_groupSEXP);
    Rcpp::traits::input_parameter< const int >::type pgw(pgwSEXP);
    rcpp_result_gen = Rcpp::wrap(Oakes_nrm(a, b, ncat, r, pni, pcni, pi, px, theta, mu, sigma, gn, pgroup, dsg_ii, dsg_gi, item_fixed, ref_group, pgw));
    return rcpp_result_gen;
END_RCPP
}
// loglikelihood_2pl
double loglikelihood_2pl(const arma::imat& a, const arma::vec& A, const arma::mat& b, const arma::ivec& ncat, const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, const arma::vec& theta, const arma::vec& mu, const arma::vec& sigma, const arma::ivec& pgroup);
RcppExport SEXP _dexterMML_loglikelihood_2pl(SEXP aSEXP, SEXP ASEXP, SEXP bSEXP, SEXP ncatSEXP, SEXP pniSEXP, SEXP pcniSEXP, SEXP piSEXP, SEXP pxSEXP, SEXP thetaSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP pgroupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::imat& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pni(pniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pcni(pcniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type px(pxSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pgroup(pgroupSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikelihood_2pl(a, A, b, ncat, pni, pcni, pi, px, theta, mu, sigma, pgroup));
    return rcpp_result_gen;
END_RCPP
}
// estimate_pl2
Rcpp::List estimate_pl2(arma::imat& a, const arma::vec& A_start, const arma::mat& b_start, const arma::ivec& ncat, const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, arma::vec& theta, const arma::vec& mu_start, const arma::vec& sigma_start, const arma::ivec& gn, const arma::ivec& pgroup, const arma::ivec& item_fixed, const arma::ivec ip, const arma::ivec& inp, const arma::ivec& icnp, const int ref_group, const int A_prior, const double A_mu, const double A_sigma, const int use_m2, const int max_iter, const int pgw);
RcppExport SEXP _dexterMML_estimate_pl2(SEXP aSEXP, SEXP A_startSEXP, SEXP b_startSEXP, SEXP ncatSEXP, SEXP pniSEXP, SEXP pcniSEXP, SEXP piSEXP, SEXP pxSEXP, SEXP thetaSEXP, SEXP mu_startSEXP, SEXP sigma_startSEXP, SEXP gnSEXP, SEXP pgroupSEXP, SEXP item_fixedSEXP, SEXP ipSEXP, SEXP inpSEXP, SEXP icnpSEXP, SEXP ref_groupSEXP, SEXP A_priorSEXP, SEXP A_muSEXP, SEXP A_sigmaSEXP, SEXP use_m2SEXP, SEXP max_iterSEXP, SEXP pgwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::imat& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type A_start(A_startSEXP);
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
    Rcpp::traits::input_parameter< const arma::ivec& >::type item_fixed(item_fixedSEXP);
    Rcpp::traits::input_parameter< const arma::ivec >::type ip(ipSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type inp(inpSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type icnp(icnpSEXP);
    Rcpp::traits::input_parameter< const int >::type ref_group(ref_groupSEXP);
    Rcpp::traits::input_parameter< const int >::type A_prior(A_priorSEXP);
    Rcpp::traits::input_parameter< const double >::type A_mu(A_muSEXP);
    Rcpp::traits::input_parameter< const double >::type A_sigma(A_sigmaSEXP);
    Rcpp::traits::input_parameter< const int >::type use_m2(use_m2SEXP);
    Rcpp::traits::input_parameter< const int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< const int >::type pgw(pgwSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_pl2(a, A_start, b_start, ncat, pni, pcni, pi, px, theta, mu_start, sigma_start, gn, pgroup, item_fixed, ip, inp, icnp, ref_group, A_prior, A_mu, A_sigma, use_m2, max_iter, pgw));
    return rcpp_result_gen;
END_RCPP
}
// test_ll_p2
double test_ll_p2(arma::ivec& a, arma::vec theta, arma::mat& r, const arma::vec& par, const int prior);
RcppExport SEXP _dexterMML_test_ll_p2(SEXP aSEXP, SEXP thetaSEXP, SEXP rSEXP, SEXP parSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type r(rSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const int >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(test_ll_p2(a, theta, r, par, prior));
    return rcpp_result_gen;
END_RCPP
}
// test_gradient_p2
arma::vec test_gradient_p2(arma::ivec& a, arma::vec theta, arma::mat& r, const arma::vec& par, const int prior);
RcppExport SEXP _dexterMML_test_gradient_p2(SEXP aSEXP, SEXP thetaSEXP, SEXP rSEXP, SEXP parSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type r(rSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const int >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(test_gradient_p2(a, theta, r, par, prior));
    return rcpp_result_gen;
END_RCPP
}
// test_hess_p2
arma::mat test_hess_p2(arma::ivec& a, arma::vec theta, arma::mat& r, const arma::vec& par, const int prior);
RcppExport SEXP _dexterMML_test_hess_p2(SEXP aSEXP, SEXP thetaSEXP, SEXP rSEXP, SEXP parSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type r(rSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const int >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(test_hess_p2(a, theta, r, par, prior));
    return rcpp_result_gen;
END_RCPP
}
// test_ll_v2
double test_ll_v2(const arma::imat& a, const arma::vec& A, const arma::mat& b, const arma::ivec& ncat, const arma::vec& theta, const arma::ivec& ip, const arma::ivec& pi, const arma::ivec& pcni, const arma::ivec& px, const arma::ivec& pgroup, const arma::ivec& inp, const arma::ivec& icnp, const arma::vec& mu, const arma::vec& sigma, const int item, const arma::vec& pars);
RcppExport SEXP _dexterMML_test_ll_v2(SEXP aSEXP, SEXP ASEXP, SEXP bSEXP, SEXP ncatSEXP, SEXP thetaSEXP, SEXP ipSEXP, SEXP piSEXP, SEXP pcniSEXP, SEXP pxSEXP, SEXP pgroupSEXP, SEXP inpSEXP, SEXP icnpSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP itemSEXP, SEXP parsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::imat& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ip(ipSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pcni(pcniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type px(pxSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pgroup(pgroupSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type inp(inpSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type icnp(icnpSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const int >::type item(itemSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pars(parsSEXP);
    rcpp_result_gen = Rcpp::wrap(test_ll_v2(a, A, b, ncat, theta, ip, pi, pcni, px, pgroup, inp, icnp, mu, sigma, item, pars));
    return rcpp_result_gen;
END_RCPP
}
// test_gr_v2
arma::vec test_gr_v2(const arma::imat& a, const arma::vec& A, const arma::mat& b, const arma::ivec& ncat, const arma::vec& theta, const arma::ivec& ip, const arma::ivec& pi, const arma::ivec& pcni, const arma::ivec& px, const arma::ivec& pgroup, const arma::ivec& inp, const arma::ivec& icnp, const arma::vec& mu, const arma::vec& sigma, const int item, const arma::vec& pars);
RcppExport SEXP _dexterMML_test_gr_v2(SEXP aSEXP, SEXP ASEXP, SEXP bSEXP, SEXP ncatSEXP, SEXP thetaSEXP, SEXP ipSEXP, SEXP piSEXP, SEXP pcniSEXP, SEXP pxSEXP, SEXP pgroupSEXP, SEXP inpSEXP, SEXP icnpSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP itemSEXP, SEXP parsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::imat& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ip(ipSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pcni(pcniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type px(pxSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pgroup(pgroupSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type inp(inpSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type icnp(icnpSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const int >::type item(itemSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pars(parsSEXP);
    rcpp_result_gen = Rcpp::wrap(test_gr_v2(a, A, b, ncat, theta, ip, pi, pcni, px, pgroup, inp, icnp, mu, sigma, item, pars));
    return rcpp_result_gen;
END_RCPP
}
// test_hess_v2
arma::mat test_hess_v2(const arma::imat& a, const arma::vec& A, const arma::mat& b, const arma::ivec& ncat, const arma::vec& theta, const arma::ivec& ip, const arma::ivec& pi, const arma::ivec& pcni, const arma::ivec& px, const arma::ivec& pgroup, const arma::ivec& inp, const arma::ivec& icnp, const arma::vec& mu, const arma::vec& sigma, const int item, const arma::vec& pars);
RcppExport SEXP _dexterMML_test_hess_v2(SEXP aSEXP, SEXP ASEXP, SEXP bSEXP, SEXP ncatSEXP, SEXP thetaSEXP, SEXP ipSEXP, SEXP piSEXP, SEXP pcniSEXP, SEXP pxSEXP, SEXP pgroupSEXP, SEXP inpSEXP, SEXP icnpSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP itemSEXP, SEXP parsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::imat& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ip(ipSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pcni(pcniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type px(pxSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pgroup(pgroupSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type inp(inpSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type icnp(icnpSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const int >::type item(itemSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pars(parsSEXP);
    rcpp_result_gen = Rcpp::wrap(test_hess_v2(a, A, b, ncat, theta, ip, pi, pcni, px, pgroup, inp, icnp, mu, sigma, item, pars));
    return rcpp_result_gen;
END_RCPP
}
// Oakes_pl2
Rcpp::List Oakes_pl2(arma::imat& a, const arma::vec& A, const arma::mat& b, const arma::ivec& ncat, arma::field<arma::mat>& r, const arma::ivec& pni, const arma::ivec& pcni, const arma::ivec& pi, const arma::ivec& px, arma::vec& theta, const arma::vec& mu, const arma::vec& sigma, const arma::ivec& gn, const arma::ivec& pgroup, const arma::imat& dsg_ii, const arma::imat& dsg_gi, const arma::ivec& item_fixed, const arma::ivec ip, const arma::ivec& inp, const arma::ivec& icnp, const int ref_group, const int A_prior, const double A_mu, const double A_sigma, const int pgw);
RcppExport SEXP _dexterMML_Oakes_pl2(SEXP aSEXP, SEXP ASEXP, SEXP bSEXP, SEXP ncatSEXP, SEXP rSEXP, SEXP pniSEXP, SEXP pcniSEXP, SEXP piSEXP, SEXP pxSEXP, SEXP thetaSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP gnSEXP, SEXP pgroupSEXP, SEXP dsg_iiSEXP, SEXP dsg_giSEXP, SEXP item_fixedSEXP, SEXP ipSEXP, SEXP inpSEXP, SEXP icnpSEXP, SEXP ref_groupSEXP, SEXP A_priorSEXP, SEXP A_muSEXP, SEXP A_sigmaSEXP, SEXP pgwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::imat& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::mat>& >::type r(rSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pni(pniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pcni(pcniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type px(pxSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type gn(gnSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pgroup(pgroupSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type dsg_ii(dsg_iiSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type dsg_gi(dsg_giSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type item_fixed(item_fixedSEXP);
    Rcpp::traits::input_parameter< const arma::ivec >::type ip(ipSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type inp(inpSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type icnp(icnpSEXP);
    Rcpp::traits::input_parameter< const int >::type ref_group(ref_groupSEXP);
    Rcpp::traits::input_parameter< const int >::type A_prior(A_priorSEXP);
    Rcpp::traits::input_parameter< const double >::type A_mu(A_muSEXP);
    Rcpp::traits::input_parameter< const double >::type A_sigma(A_sigmaSEXP);
    Rcpp::traits::input_parameter< const int >::type pgw(pgwSEXP);
    rcpp_result_gen = Rcpp::wrap(Oakes_pl2(a, A, b, ncat, r, pni, pcni, pi, px, theta, mu, sigma, gn, pgroup, dsg_ii, dsg_gi, item_fixed, ip, inp, icnp, ref_group, A_prior, A_mu, A_sigma, pgw));
    return rcpp_result_gen;
END_RCPP
}
// plausible_values_c
arma::mat plausible_values_c(const arma::vec& A, const arma::imat& a, const arma::mat& b, const arma::ivec& ncat, const arma::ivec& pni, const arma::ivec& pcni, arma::ivec& pi, const arma::ivec& px, const arma::ivec& pop, const arma::ivec& popn, const int npv, const arma::vec& starting_values, const int n_prior_updates, const int thin, const int pgw);
RcppExport SEXP _dexterMML_plausible_values_c(SEXP ASEXP, SEXP aSEXP, SEXP bSEXP, SEXP ncatSEXP, SEXP pniSEXP, SEXP pcniSEXP, SEXP piSEXP, SEXP pxSEXP, SEXP popSEXP, SEXP popnSEXP, SEXP npvSEXP, SEXP starting_valuesSEXP, SEXP n_prior_updatesSEXP, SEXP thinSEXP, SEXP pgwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pni(pniSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pcni(pcniSEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type px(pxSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pop(popSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type popn(popnSEXP);
    Rcpp::traits::input_parameter< const int >::type npv(npvSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type starting_values(starting_valuesSEXP);
    Rcpp::traits::input_parameter< const int >::type n_prior_updates(n_prior_updatesSEXP);
    Rcpp::traits::input_parameter< const int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const int >::type pgw(pgwSEXP);
    rcpp_result_gen = Rcpp::wrap(plausible_values_c(A, a, b, ncat, pni, pcni, pi, px, pop, popn, npv, starting_values, n_prior_updates, thin, pgw));
    return rcpp_result_gen;
END_RCPP
}
// sim_2plc
arma::imat sim_2plc(const arma::imat& a, const arma::vec& A, const arma::mat& b, const arma::ivec& ncat, const arma::vec& theta);
RcppExport SEXP _dexterMML_sim_2plc(SEXP aSEXP, SEXP ASEXP, SEXP bSEXP, SEXP ncatSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::imat& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_2plc(a, A, b, ncat, theta));
    return rcpp_result_gen;
END_RCPP
}
// test_nlm
void test_nlm(arma::ivec& a, arma::vec theta, arma::mat& r, const arma::vec& par_in);
RcppExport SEXP _dexterMML_test_nlm(SEXP aSEXP, SEXP thetaSEXP, SEXP rSEXP, SEXP par_inSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type r(rSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type par_in(par_inSEXP);
    test_nlm(a, theta, r, par_in);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dexterMML_theta_2pl", (DL_FUNC) &_dexterMML_theta_2pl, 10},
    {"_dexterMML_E_score", (DL_FUNC) &_dexterMML_E_score, 6},
    {"_dexterMML_mat_pre", (DL_FUNC) &_dexterMML_mat_pre, 2},
    {"_dexterMML_categorize", (DL_FUNC) &_dexterMML_categorize, 7},
    {"_dexterMML_plot_data", (DL_FUNC) &_dexterMML_plot_data, 7},
    {"_dexterMML_design_matrices", (DL_FUNC) &_dexterMML_design_matrices, 6},
    {"_dexterMML_check_connected_c", (DL_FUNC) &_dexterMML_check_connected_c, 3},
    {"_dexterMML_estimate_nrm", (DL_FUNC) &_dexterMML_estimate_nrm, 16},
    {"_dexterMML_test_ll_nrm", (DL_FUNC) &_dexterMML_test_ll_nrm, 4},
    {"_dexterMML_test_gradient_nrm", (DL_FUNC) &_dexterMML_test_gradient_nrm, 4},
    {"_dexterMML_test_hess_nrm", (DL_FUNC) &_dexterMML_test_hess_nrm, 4},
    {"_dexterMML_Oakes_nrm", (DL_FUNC) &_dexterMML_Oakes_nrm, 18},
    {"_dexterMML_loglikelihood_2pl", (DL_FUNC) &_dexterMML_loglikelihood_2pl, 12},
    {"_dexterMML_estimate_pl2", (DL_FUNC) &_dexterMML_estimate_pl2, 24},
    {"_dexterMML_test_ll_p2", (DL_FUNC) &_dexterMML_test_ll_p2, 5},
    {"_dexterMML_test_gradient_p2", (DL_FUNC) &_dexterMML_test_gradient_p2, 5},
    {"_dexterMML_test_hess_p2", (DL_FUNC) &_dexterMML_test_hess_p2, 5},
    {"_dexterMML_test_ll_v2", (DL_FUNC) &_dexterMML_test_ll_v2, 16},
    {"_dexterMML_test_gr_v2", (DL_FUNC) &_dexterMML_test_gr_v2, 16},
    {"_dexterMML_test_hess_v2", (DL_FUNC) &_dexterMML_test_hess_v2, 16},
    {"_dexterMML_Oakes_pl2", (DL_FUNC) &_dexterMML_Oakes_pl2, 25},
    {"_dexterMML_plausible_values_c", (DL_FUNC) &_dexterMML_plausible_values_c, 15},
    {"_dexterMML_sim_2plc", (DL_FUNC) &_dexterMML_sim_2plc, 5},
    {"_dexterMML_test_nlm", (DL_FUNC) &_dexterMML_test_nlm, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_dexterMML(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
