// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// seq_rcpp
Rcpp::NumericVector seq_rcpp(double from_, double to_, double by_);
RcppExport SEXP _GeoSim_seq_rcpp(SEXP from_SEXP, SEXP to_SEXP, SEXP by_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type from_(from_SEXP);
    Rcpp::traits::input_parameter< double >::type to_(to_SEXP);
    Rcpp::traits::input_parameter< double >::type by_(by_SEXP);
    rcpp_result_gen = Rcpp::wrap(seq_rcpp(from_, to_, by_));
    return rcpp_result_gen;
END_RCPP
}
// cokrige_cpp
arma::field<arma::mat> cokrige_cpp(arma::mat datacoord, arma::uvec index_missing, arma::mat coord, arma::mat model, arma::cube sill, arma::mat b, arma::mat sillnugget, arma::cube model_rotationmatrix);
RcppExport SEXP _GeoSim_cokrige_cpp(SEXP datacoordSEXP, SEXP index_missingSEXP, SEXP coordSEXP, SEXP modelSEXP, SEXP sillSEXP, SEXP bSEXP, SEXP sillnuggetSEXP, SEXP model_rotationmatrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type datacoord(datacoordSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type index_missing(index_missingSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type model(modelSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sill(sillSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sillnugget(sillnuggetSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type model_rotationmatrix(model_rotationmatrixSEXP);
    rcpp_result_gen = Rcpp::wrap(cokrige_cpp(datacoord, index_missing, coord, model, sill, b, sillnugget, model_rotationmatrix));
    return rcpp_result_gen;
END_RCPP
}
// cc_truncate_cpp
int cc_truncate_cpp(arma::mat y_simu, int nfield, arma::vec flag, arma::vec nthres, arma::vec thresholds);
RcppExport SEXP _GeoSim_cc_truncate_cpp(SEXP y_simuSEXP, SEXP nfieldSEXP, SEXP flagSEXP, SEXP nthresSEXP, SEXP thresholdsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y_simu(y_simuSEXP);
    Rcpp::traits::input_parameter< int >::type nfield(nfieldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type flag(flagSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nthres(nthresSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type thresholds(thresholdsSEXP);
    rcpp_result_gen = Rcpp::wrap(cc_truncate_cpp(y_simu, nfield, flag, nthres, thresholds));
    return rcpp_result_gen;
END_RCPP
}
// cc_search_cpp
arma::field<arma::mat> cc_search_cpp(arma::mat datacoord, arma::mat datavalue, arma::rowvec coord, arma::mat search_rotationmatrix, int octant, int ndata, int nxsup, int nysup, int nzsup, int xmnsup, int ymnsup, int zmnsup, int xsizsup, int ysizsup, int zsizsup, arma::rowvec ixsbtosr, arma::rowvec iysbtosr, arma::rowvec izsbtosr, arma::colvec nisb);
RcppExport SEXP _GeoSim_cc_search_cpp(SEXP datacoordSEXP, SEXP datavalueSEXP, SEXP coordSEXP, SEXP search_rotationmatrixSEXP, SEXP octantSEXP, SEXP ndataSEXP, SEXP nxsupSEXP, SEXP nysupSEXP, SEXP nzsupSEXP, SEXP xmnsupSEXP, SEXP ymnsupSEXP, SEXP zmnsupSEXP, SEXP xsizsupSEXP, SEXP ysizsupSEXP, SEXP zsizsupSEXP, SEXP ixsbtosrSEXP, SEXP iysbtosrSEXP, SEXP izsbtosrSEXP, SEXP nisbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type datacoord(datacoordSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type datavalue(datavalueSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type search_rotationmatrix(search_rotationmatrixSEXP);
    Rcpp::traits::input_parameter< int >::type octant(octantSEXP);
    Rcpp::traits::input_parameter< int >::type ndata(ndataSEXP);
    Rcpp::traits::input_parameter< int >::type nxsup(nxsupSEXP);
    Rcpp::traits::input_parameter< int >::type nysup(nysupSEXP);
    Rcpp::traits::input_parameter< int >::type nzsup(nzsupSEXP);
    Rcpp::traits::input_parameter< int >::type xmnsup(xmnsupSEXP);
    Rcpp::traits::input_parameter< int >::type ymnsup(ymnsupSEXP);
    Rcpp::traits::input_parameter< int >::type zmnsup(zmnsupSEXP);
    Rcpp::traits::input_parameter< int >::type xsizsup(xsizsupSEXP);
    Rcpp::traits::input_parameter< int >::type ysizsup(ysizsupSEXP);
    Rcpp::traits::input_parameter< int >::type zsizsup(zsizsupSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type ixsbtosr(ixsbtosrSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type iysbtosr(iysbtosrSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type izsbtosr(izsbtosrSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type nisb(nisbSEXP);
    rcpp_result_gen = Rcpp::wrap(cc_search_cpp(datacoord, datavalue, coord, search_rotationmatrix, octant, ndata, nxsup, nysup, nzsup, xmnsup, ymnsup, zmnsup, xsizsup, ysizsup, zsizsup, ixsbtosr, iysbtosr, izsbtosr, nisb));
    return rcpp_result_gen;
END_RCPP
}
// Gibbs_cosim_cpp
arma::mat Gibbs_cosim_cpp(arma::mat datacoord, arma::mat idata, arma::mat ydata, int nfield, arma::vec flag, arma::vec nthres, arma::vec thresholds, arma::mat model, arma::cube sill, arma::mat b, arma::mat sillnugget, int nrealiz, int niterations, arma::cube model_rotationmatrix, arma::mat search_rotationmatrix, int cc_unique, int octant, int ndata, int nxsup, int nysup, int nzsup, int xmnsup, int ymnsup, int zmnsup, int xsizsup, int ysizsup, int zsizsup, arma::rowvec ixsbtosr, arma::rowvec iysbtosr, arma::rowvec izsbtosr, arma::colvec nisb);
RcppExport SEXP _GeoSim_Gibbs_cosim_cpp(SEXP datacoordSEXP, SEXP idataSEXP, SEXP ydataSEXP, SEXP nfieldSEXP, SEXP flagSEXP, SEXP nthresSEXP, SEXP thresholdsSEXP, SEXP modelSEXP, SEXP sillSEXP, SEXP bSEXP, SEXP sillnuggetSEXP, SEXP nrealizSEXP, SEXP niterationsSEXP, SEXP model_rotationmatrixSEXP, SEXP search_rotationmatrixSEXP, SEXP cc_uniqueSEXP, SEXP octantSEXP, SEXP ndataSEXP, SEXP nxsupSEXP, SEXP nysupSEXP, SEXP nzsupSEXP, SEXP xmnsupSEXP, SEXP ymnsupSEXP, SEXP zmnsupSEXP, SEXP xsizsupSEXP, SEXP ysizsupSEXP, SEXP zsizsupSEXP, SEXP ixsbtosrSEXP, SEXP iysbtosrSEXP, SEXP izsbtosrSEXP, SEXP nisbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type datacoord(datacoordSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type idata(idataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ydata(ydataSEXP);
    Rcpp::traits::input_parameter< int >::type nfield(nfieldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type flag(flagSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nthres(nthresSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type thresholds(thresholdsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type model(modelSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sill(sillSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sillnugget(sillnuggetSEXP);
    Rcpp::traits::input_parameter< int >::type nrealiz(nrealizSEXP);
    Rcpp::traits::input_parameter< int >::type niterations(niterationsSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type model_rotationmatrix(model_rotationmatrixSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type search_rotationmatrix(search_rotationmatrixSEXP);
    Rcpp::traits::input_parameter< int >::type cc_unique(cc_uniqueSEXP);
    Rcpp::traits::input_parameter< int >::type octant(octantSEXP);
    Rcpp::traits::input_parameter< int >::type ndata(ndataSEXP);
    Rcpp::traits::input_parameter< int >::type nxsup(nxsupSEXP);
    Rcpp::traits::input_parameter< int >::type nysup(nysupSEXP);
    Rcpp::traits::input_parameter< int >::type nzsup(nzsupSEXP);
    Rcpp::traits::input_parameter< int >::type xmnsup(xmnsupSEXP);
    Rcpp::traits::input_parameter< int >::type ymnsup(ymnsupSEXP);
    Rcpp::traits::input_parameter< int >::type zmnsup(zmnsupSEXP);
    Rcpp::traits::input_parameter< int >::type xsizsup(xsizsupSEXP);
    Rcpp::traits::input_parameter< int >::type ysizsup(ysizsupSEXP);
    Rcpp::traits::input_parameter< int >::type zsizsup(zsizsupSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type ixsbtosr(ixsbtosrSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type iysbtosr(iysbtosrSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type izsbtosr(izsbtosrSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type nisb(nisbSEXP);
    rcpp_result_gen = Rcpp::wrap(Gibbs_cosim_cpp(datacoord, idata, ydata, nfield, flag, nthres, thresholds, model, sill, b, sillnugget, nrealiz, niterations, model_rotationmatrix, search_rotationmatrix, cc_unique, octant, ndata, nxsup, nysup, nzsup, xmnsup, ymnsup, zmnsup, xsizsup, ysizsup, zsizsup, ixsbtosr, iysbtosr, izsbtosr, nisb));
    return rcpp_result_gen;
END_RCPP
}
// cond_mooving_neigbor
arma::mat cond_mooving_neigbor(int m1, arma::mat simu, int nvar, int nrealiz, arma::mat datacoord, arma::mat cc_residuals, arma::mat coord, arma::mat search_rotationmatrix, int octant, int ndata, int nxsup, int nysup, int nzsup, int xmnsup, int ymnsup, int zmnsup, int xsizsup, int ysizsup, int zsizsup, arma::rowvec ixsbtosr, arma::rowvec iysbtosr, arma::rowvec izsbtosr, arma::colvec nisb, arma::mat model, arma::cube sill, arma::mat b, arma::mat sillnugget, arma::cube model_rotationmatrix);
RcppExport SEXP _GeoSim_cond_mooving_neigbor(SEXP m1SEXP, SEXP simuSEXP, SEXP nvarSEXP, SEXP nrealizSEXP, SEXP datacoordSEXP, SEXP cc_residualsSEXP, SEXP coordSEXP, SEXP search_rotationmatrixSEXP, SEXP octantSEXP, SEXP ndataSEXP, SEXP nxsupSEXP, SEXP nysupSEXP, SEXP nzsupSEXP, SEXP xmnsupSEXP, SEXP ymnsupSEXP, SEXP zmnsupSEXP, SEXP xsizsupSEXP, SEXP ysizsupSEXP, SEXP zsizsupSEXP, SEXP ixsbtosrSEXP, SEXP iysbtosrSEXP, SEXP izsbtosrSEXP, SEXP nisbSEXP, SEXP modelSEXP, SEXP sillSEXP, SEXP bSEXP, SEXP sillnuggetSEXP, SEXP model_rotationmatrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type simu(simuSEXP);
    Rcpp::traits::input_parameter< int >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< int >::type nrealiz(nrealizSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type datacoord(datacoordSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cc_residuals(cc_residualsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type search_rotationmatrix(search_rotationmatrixSEXP);
    Rcpp::traits::input_parameter< int >::type octant(octantSEXP);
    Rcpp::traits::input_parameter< int >::type ndata(ndataSEXP);
    Rcpp::traits::input_parameter< int >::type nxsup(nxsupSEXP);
    Rcpp::traits::input_parameter< int >::type nysup(nysupSEXP);
    Rcpp::traits::input_parameter< int >::type nzsup(nzsupSEXP);
    Rcpp::traits::input_parameter< int >::type xmnsup(xmnsupSEXP);
    Rcpp::traits::input_parameter< int >::type ymnsup(ymnsupSEXP);
    Rcpp::traits::input_parameter< int >::type zmnsup(zmnsupSEXP);
    Rcpp::traits::input_parameter< int >::type xsizsup(xsizsupSEXP);
    Rcpp::traits::input_parameter< int >::type ysizsup(ysizsupSEXP);
    Rcpp::traits::input_parameter< int >::type zsizsup(zsizsupSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type ixsbtosr(ixsbtosrSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type iysbtosr(iysbtosrSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type izsbtosr(izsbtosrSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type nisb(nisbSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type model(modelSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sill(sillSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sillnugget(sillnuggetSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type model_rotationmatrix(model_rotationmatrixSEXP);
    rcpp_result_gen = Rcpp::wrap(cond_mooving_neigbor(m1, simu, nvar, nrealiz, datacoord, cc_residuals, coord, search_rotationmatrix, octant, ndata, nxsup, nysup, nzsup, xmnsup, ymnsup, zmnsup, xsizsup, ysizsup, zsizsup, ixsbtosr, iysbtosr, izsbtosr, nisb, model, sill, b, sillnugget, model_rotationmatrix));
    return rcpp_result_gen;
END_RCPP
}
// Gibbs_plurisim_cpp
arma::mat Gibbs_plurisim_cpp(arma::mat datacoord, arma::mat idata, int nfield, arma::vec flag, arma::vec nthres, arma::vec thresholds, arma::mat model, arma::cube sill, arma::mat b, arma::mat sillnugget, int nrealiz, int niterations, arma::cube model_rotationmatrix, arma::mat search_rotationmatrix, int cc_unique, int octant, int ndata, int nxsup, int nysup, int nzsup, int xmnsup, int ymnsup, int zmnsup, int xsizsup, int ysizsup, int zsizsup, arma::rowvec ixsbtosr, arma::rowvec iysbtosr, arma::rowvec izsbtosr, arma::colvec nisb);
RcppExport SEXP _GeoSim_Gibbs_plurisim_cpp(SEXP datacoordSEXP, SEXP idataSEXP, SEXP nfieldSEXP, SEXP flagSEXP, SEXP nthresSEXP, SEXP thresholdsSEXP, SEXP modelSEXP, SEXP sillSEXP, SEXP bSEXP, SEXP sillnuggetSEXP, SEXP nrealizSEXP, SEXP niterationsSEXP, SEXP model_rotationmatrixSEXP, SEXP search_rotationmatrixSEXP, SEXP cc_uniqueSEXP, SEXP octantSEXP, SEXP ndataSEXP, SEXP nxsupSEXP, SEXP nysupSEXP, SEXP nzsupSEXP, SEXP xmnsupSEXP, SEXP ymnsupSEXP, SEXP zmnsupSEXP, SEXP xsizsupSEXP, SEXP ysizsupSEXP, SEXP zsizsupSEXP, SEXP ixsbtosrSEXP, SEXP iysbtosrSEXP, SEXP izsbtosrSEXP, SEXP nisbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type datacoord(datacoordSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type idata(idataSEXP);
    Rcpp::traits::input_parameter< int >::type nfield(nfieldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type flag(flagSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nthres(nthresSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type thresholds(thresholdsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type model(modelSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sill(sillSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sillnugget(sillnuggetSEXP);
    Rcpp::traits::input_parameter< int >::type nrealiz(nrealizSEXP);
    Rcpp::traits::input_parameter< int >::type niterations(niterationsSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type model_rotationmatrix(model_rotationmatrixSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type search_rotationmatrix(search_rotationmatrixSEXP);
    Rcpp::traits::input_parameter< int >::type cc_unique(cc_uniqueSEXP);
    Rcpp::traits::input_parameter< int >::type octant(octantSEXP);
    Rcpp::traits::input_parameter< int >::type ndata(ndataSEXP);
    Rcpp::traits::input_parameter< int >::type nxsup(nxsupSEXP);
    Rcpp::traits::input_parameter< int >::type nysup(nysupSEXP);
    Rcpp::traits::input_parameter< int >::type nzsup(nzsupSEXP);
    Rcpp::traits::input_parameter< int >::type xmnsup(xmnsupSEXP);
    Rcpp::traits::input_parameter< int >::type ymnsup(ymnsupSEXP);
    Rcpp::traits::input_parameter< int >::type zmnsup(zmnsupSEXP);
    Rcpp::traits::input_parameter< int >::type xsizsup(xsizsupSEXP);
    Rcpp::traits::input_parameter< int >::type ysizsup(ysizsupSEXP);
    Rcpp::traits::input_parameter< int >::type zsizsup(zsizsupSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type ixsbtosr(ixsbtosrSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type iysbtosr(iysbtosrSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type izsbtosr(izsbtosrSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type nisb(nisbSEXP);
    rcpp_result_gen = Rcpp::wrap(Gibbs_plurisim_cpp(datacoord, idata, nfield, flag, nthres, thresholds, model, sill, b, sillnugget, nrealiz, niterations, model_rotationmatrix, search_rotationmatrix, cc_unique, octant, ndata, nxsup, nysup, nzsup, xmnsup, ymnsup, zmnsup, xsizsup, ysizsup, zsizsup, ixsbtosr, iysbtosr, izsbtosr, nisb));
    return rcpp_result_gen;
END_RCPP
}
// set_seed
void set_seed(double seed);
RcppExport SEXP _GeoSim_set_seed(SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type seed(seedSEXP);
    set_seed(seed);
    return R_NilValue;
END_RCPP
}
// tbmain_cosimu_cpp
arma::mat tbmain_cosimu_cpp(arma::mat coord, arma::mat model, arma::vec cc_sigma, arma::cube A1, int nvar, arma::vec nlines, int nrealiz, arma::mat seed, arma::cube all_lines, arma::mat all_offset, arma::cube all_r, arma::cube all_phi, arma::cube all_theta, arma::mat valid_lines);
RcppExport SEXP _GeoSim_tbmain_cosimu_cpp(SEXP coordSEXP, SEXP modelSEXP, SEXP cc_sigmaSEXP, SEXP A1SEXP, SEXP nvarSEXP, SEXP nlinesSEXP, SEXP nrealizSEXP, SEXP seedSEXP, SEXP all_linesSEXP, SEXP all_offsetSEXP, SEXP all_rSEXP, SEXP all_phiSEXP, SEXP all_thetaSEXP, SEXP valid_linesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type model(modelSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cc_sigma(cc_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type A1(A1SEXP);
    Rcpp::traits::input_parameter< int >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nlines(nlinesSEXP);
    Rcpp::traits::input_parameter< int >::type nrealiz(nrealizSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type all_lines(all_linesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type all_offset(all_offsetSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type all_r(all_rSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type all_phi(all_phiSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type all_theta(all_thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type valid_lines(valid_linesSEXP);
    rcpp_result_gen = Rcpp::wrap(tbmain_cosimu_cpp(coord, model, cc_sigma, A1, nvar, nlines, nrealiz, seed, all_lines, all_offset, all_r, all_phi, all_theta, valid_lines));
    return rcpp_result_gen;
END_RCPP
}
// tbmain_simu_cpp
arma::mat tbmain_simu_cpp(arma::mat coord, arma::mat model, arma::vec cc_sigma, arma::cube A1, int nvar, arma::vec nlines, int nrealiz, arma::mat seed, arma::cube all_lines, arma::mat all_offset, arma::cube all_r, arma::cube all_phi, arma::cube all_theta, arma::mat valid_lines);
RcppExport SEXP _GeoSim_tbmain_simu_cpp(SEXP coordSEXP, SEXP modelSEXP, SEXP cc_sigmaSEXP, SEXP A1SEXP, SEXP nvarSEXP, SEXP nlinesSEXP, SEXP nrealizSEXP, SEXP seedSEXP, SEXP all_linesSEXP, SEXP all_offsetSEXP, SEXP all_rSEXP, SEXP all_phiSEXP, SEXP all_thetaSEXP, SEXP valid_linesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type model(modelSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cc_sigma(cc_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type A1(A1SEXP);
    Rcpp::traits::input_parameter< int >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nlines(nlinesSEXP);
    Rcpp::traits::input_parameter< int >::type nrealiz(nrealizSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type all_lines(all_linesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type all_offset(all_offsetSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type all_r(all_rSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type all_phi(all_phiSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type all_theta(all_thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type valid_lines(valid_linesSEXP);
    rcpp_result_gen = Rcpp::wrap(tbmain_simu_cpp(coord, model, cc_sigma, A1, nvar, nlines, nrealiz, seed, all_lines, all_offset, all_r, all_phi, all_theta, valid_lines));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GeoSim_seq_rcpp", (DL_FUNC) &_GeoSim_seq_rcpp, 3},
    {"_GeoSim_cokrige_cpp", (DL_FUNC) &_GeoSim_cokrige_cpp, 8},
    {"_GeoSim_cc_truncate_cpp", (DL_FUNC) &_GeoSim_cc_truncate_cpp, 5},
    {"_GeoSim_cc_search_cpp", (DL_FUNC) &_GeoSim_cc_search_cpp, 19},
    {"_GeoSim_Gibbs_cosim_cpp", (DL_FUNC) &_GeoSim_Gibbs_cosim_cpp, 31},
    {"_GeoSim_cond_mooving_neigbor", (DL_FUNC) &_GeoSim_cond_mooving_neigbor, 28},
    {"_GeoSim_Gibbs_plurisim_cpp", (DL_FUNC) &_GeoSim_Gibbs_plurisim_cpp, 30},
    {"_GeoSim_set_seed", (DL_FUNC) &_GeoSim_set_seed, 1},
    {"_GeoSim_tbmain_cosimu_cpp", (DL_FUNC) &_GeoSim_tbmain_cosimu_cpp, 14},
    {"_GeoSim_tbmain_simu_cpp", (DL_FUNC) &_GeoSim_tbmain_simu_cpp, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_GeoSim(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
