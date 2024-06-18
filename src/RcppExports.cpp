// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <RcppThread.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// select_yvine_cpp
Rcpp::List select_yvine_cpp(const Eigen::MatrixXd& data, std::vector<std::string> family_set, std::string par_method, std::string nonpar_method, double mult, std::string selcrit, const Eigen::VectorXd& weights, double psi0, bool preselect_families, size_t cores, const std::vector<std::string>& var_types);
RcppExport SEXP _bivinereg_select_yvine_cpp(SEXP dataSEXP, SEXP family_setSEXP, SEXP par_methodSEXP, SEXP nonpar_methodSEXP, SEXP multSEXP, SEXP selcritSEXP, SEXP weightsSEXP, SEXP psi0SEXP, SEXP preselect_familiesSEXP, SEXP coresSEXP, SEXP var_typesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type family_set(family_setSEXP);
    Rcpp::traits::input_parameter< std::string >::type par_method(par_methodSEXP);
    Rcpp::traits::input_parameter< std::string >::type nonpar_method(nonpar_methodSEXP);
    Rcpp::traits::input_parameter< double >::type mult(multSEXP);
    Rcpp::traits::input_parameter< std::string >::type selcrit(selcritSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< double >::type psi0(psi0SEXP);
    Rcpp::traits::input_parameter< bool >::type preselect_families(preselect_familiesSEXP);
    Rcpp::traits::input_parameter< size_t >::type cores(coresSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type var_types(var_typesSEXP);
    rcpp_result_gen = Rcpp::wrap(select_yvine_cpp(data, family_set, par_method, nonpar_method, mult, selcrit, weights, psi0, preselect_families, cores, var_types));
    return rcpp_result_gen;
END_RCPP
}
// cond_m_dist_cpp
Eigen::VectorXd cond_m_dist_cpp(const Eigen::MatrixXd& u, const Rcpp::List& vinecop_r, const int margin, size_t num_threads);
RcppExport SEXP _bivinereg_cond_m_dist_cpp(SEXP uSEXP, SEXP vinecop_rSEXP, SEXP marginSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type vinecop_r(vinecop_rSEXP);
    Rcpp::traits::input_parameter< const int >::type margin(marginSEXP);
    Rcpp::traits::input_parameter< size_t >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(cond_m_dist_cpp(u, vinecop_r, margin, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// cond_m2_dist_cpp
Eigen::VectorXd cond_m2_dist_cpp(const Eigen::MatrixXd& u, const Rcpp::List& vinecop_r, const int margin, size_t num_threads);
RcppExport SEXP _bivinereg_cond_m2_dist_cpp(SEXP uSEXP, SEXP vinecop_rSEXP, SEXP marginSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type vinecop_r(vinecop_rSEXP);
    Rcpp::traits::input_parameter< const int >::type margin(marginSEXP);
    Rcpp::traits::input_parameter< size_t >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(cond_m2_dist_cpp(u, vinecop_r, margin, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// cond_bi_dens_cpp
Eigen::VectorXd cond_bi_dens_cpp(const Eigen::MatrixXd& u, const Rcpp::List& vinecop_r, size_t num_threads);
RcppExport SEXP _bivinereg_cond_bi_dens_cpp(SEXP uSEXP, SEXP vinecop_rSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type vinecop_r(vinecop_rSEXP);
    Rcpp::traits::input_parameter< size_t >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(cond_bi_dens_cpp(u, vinecop_r, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// cond_m_dens_cpp
Eigen::VectorXd cond_m_dens_cpp(const Eigen::MatrixXd& u, const Rcpp::List& vinecop_r, const int margin, size_t num_threads);
RcppExport SEXP _bivinereg_cond_m_dens_cpp(SEXP uSEXP, SEXP vinecop_rSEXP, SEXP marginSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type vinecop_r(vinecop_rSEXP);
    Rcpp::traits::input_parameter< const int >::type margin(marginSEXP);
    Rcpp::traits::input_parameter< size_t >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(cond_m_dens_cpp(u, vinecop_r, margin, num_threads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bivinereg_select_yvine_cpp", (DL_FUNC) &_bivinereg_select_yvine_cpp, 11},
    {"_bivinereg_cond_m_dist_cpp", (DL_FUNC) &_bivinereg_cond_m_dist_cpp, 4},
    {"_bivinereg_cond_m2_dist_cpp", (DL_FUNC) &_bivinereg_cond_m2_dist_cpp, 4},
    {"_bivinereg_cond_bi_dens_cpp", (DL_FUNC) &_bivinereg_cond_bi_dens_cpp, 3},
    {"_bivinereg_cond_m_dens_cpp", (DL_FUNC) &_bivinereg_cond_m_dens_cpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_bivinereg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
