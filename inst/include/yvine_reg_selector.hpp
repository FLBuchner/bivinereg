#pragma once
#include <RcppThread.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/bicop/class.hpp>

namespace vinecopulib {
  class Bicop;

  namespace tools_stl {
    template<typename T>
    std::vector<T> rank(const std::vector<T>& x)
    {
      std::vector<size_t> r(x.size());
      auto order = tools_stl::get_order(x);
      for (auto i : order)
        r[order[i]] = static_cast<double>(i + 1);
      return r;
    }
  }
}

namespace bivinereg {
using namespace vinecopulib;

struct DVineFitTemporaries
{
  std::vector<Eigen::VectorXd> hfunc1;
  std::vector<Eigen::VectorXd> hfunc2;
  std::vector<Eigen::VectorXd> hfunc1_sub;
  std::vector<Eigen::VectorXd> hfunc2_sub;
  Eigen::VectorXd hfunc2_resp;
  std::vector<Bicop> pcs;
  std::vector<size_t> remaining_vars;
  std::vector<size_t> selected_vars;
  double crit;
};

class YVineRegSelector
{
public:
  YVineRegSelector(const Eigen::MatrixXd& data,
                   const std::vector<std::string>& var_types,
                   const FitControlsBicop& controls);

  void select_model();
  std::vector<size_t> get_selected_vars() const { return fit1_.selected_vars; }
  std::vector<std::vector<Bicop>> get_pcs() const { return pcs_; }

private:
  void extend_fit(DVineFitTemporaries& fit1,
                  DVineFitTemporaries& fit2,
                  size_t var) const;
  void initialize_var(DVineFitTemporaries& fit1,
                      DVineFitTemporaries& fit2,
                      size_t var) const;
  std::vector<std::string> get_edge_types(const DVineFitTemporaries& fit,
                                          size_t t) const;
  Eigen::MatrixXd get_edge_data(const DVineFitTemporaries& fit, size_t t) const;
  void fit_pair_copula(DVineFitTemporaries& fit,
                       size_t t,
                       const Eigen::MatrixXd& u_e) const;
  void fit_pair_copula(DVineFitTemporaries& fit1,
                       DVineFitTemporaries& fit2,
                       size_t t,
                       const Eigen::MatrixXd& u_e) const;
  void update_hfunc1(DVineFitTemporaries& fit,
                     size_t t,
                     const Eigen::MatrixXd& u_e) const;
  void update_hfunc1(DVineFitTemporaries& fit1,
                     DVineFitTemporaries& fit2,
                     size_t t,
                     const Eigen::MatrixXd& u_e) const;
  void update_hfunc2(DVineFitTemporaries& fit,
                     size_t t,
                     const Eigen::MatrixXd& u_e) const;
  void update_hfunc2(DVineFitTemporaries& fit1,
                     DVineFitTemporaries& fit2,
                     size_t t,
                     const Eigen::MatrixXd& u_e) const;
  void update_hfuncs(DVineFitTemporaries& fit,
                     size_t t,
                     const Eigen::MatrixXd& u_e) const;
  void update_hfuncs(DVineFitTemporaries& fit1,
                     DVineFitTemporaries& fit2,
                     size_t t,
                     const Eigen::MatrixXd& u_e) const;
  void update_selcrit(DVineFitTemporaries& fit1,
                      DVineFitTemporaries& fit2) const;
  void update_vars(DVineFitTemporaries& fit1,
                   DVineFitTemporaries& fit2,
                   size_t var) const;
  void update_status(DVineFitTemporaries& fit1,
                     DVineFitTemporaries& fit2,
                     size_t var) const;

  size_t p_;
  Eigen::MatrixXd data_;
  std::vector<std::string> var_types_;
  FitControlsBicop controls_;
  DVineFitTemporaries fit1_;
  DVineFitTemporaries fit2_;
  std::vector<std::vector<Bicop>> pcs_; // Becomes list of fitted pair copulas
};

inline YVineRegSelector::YVineRegSelector(
    const Eigen::MatrixXd& data,
    const std::vector<std::string>& var_types,
    const FitControlsBicop& controls)
  : p_(var_types.size() - 2)
  , data_(data)
  , var_types_(var_types)
  , controls_(controls)
{
  fit1_.hfunc1.resize(p_);
  fit1_.hfunc2.resize(p_);
  fit1_.hfunc1_sub.resize(p_); // sub only for discrete I think
  fit1_.hfunc2_sub.resize(p_);
  fit1_.pcs.resize(p_);
  fit1_.remaining_vars = tools_stl::seq_int(2, p_);
  fit1_.selected_vars.reserve(p_);
  fit1_.crit = 0.0;

  fit1_.hfunc2[0] = data_.col(0);
  fit1_.hfunc2_resp = data_.col(0);

  fit2_.hfunc1.resize(p_);
  fit2_.hfunc2.resize(p_);
  fit2_.hfunc1_sub.resize(p_);
  fit2_.hfunc2_sub.resize(p_);
  fit2_.pcs.resize(p_);
  fit2_.remaining_vars = tools_stl::seq_int(2, p_);
  fit2_.selected_vars.reserve(p_);
  fit2_.crit = 0.0;

  fit2_.hfunc2[0] = data_.col(1);
  fit2_.hfunc2_resp = data_.col(2);
}


inline void YVineRegSelector::select_model()
{
  std::mutex m; // required to synchronize write/reads to the selector
  auto num_threads = controls_.get_num_threads();
  RcppThread::ThreadPool pool(num_threads > 1 ? num_threads : 0);
  controls_.set_num_threads(0);

  while (fit1_.selected_vars.size() < p_) {
    auto old_fit1 = fit1_;  // fix current model (fit_ will be modified below)
    auto old_fit2 = fit2_;
    auto fit_replace_if_better = [&](size_t var) {
      DVineFitTemporaries new_fit1 = old_fit1;
      DVineFitTemporaries new_fit2 = old_fit2;
      this->extend_fit(new_fit1, new_fit2, var);
      // continue
      std::lock_guard<std::mutex> lk(m);  // synchronize
      if (new_fit1.crit + new_fit2.crit > fit1_.crit + fit2_.crit) {
        fit1_ = std::move(new_fit1);
        fit2_ = std::move(new_fit2);
      }
    };
    pool.map(fit_replace_if_better, old_fit1.remaining_vars); // is this even needed?
    pool.map(fit_replace_if_better, old_fit2.remaining_vars);
    pool.wait();

    if (fit1_.selected_vars == old_fit1.selected_vars)
      break;  // could not improve the selection criterion

    // model improved; store pair copulas for new variable
    auto p_sel = fit1_.selected_vars.size();
    // new tree with responses
    pcs_.push_back(std::vector<Bicop>{ fit1_.pcs[p_sel - 1] });
    pcs_[p_sel - 1].push_back(fit2_.pcs[p_sel - 1]);
    // add pair copulas to existing trees
    for (size_t t = 0; t < p_sel - 1; t++)
      pcs_[t].push_back(fit1_.pcs[t]);
  }
  pool.join();
  controls_.set_num_threads(num_threads);

  // need to fit pair copula between responses at the end
  auto p_sel = fit1_.selected_vars.size();
  Eigen::MatrixXd u_e(data_.rows(), 2);
  u_e.col(0) = fit1_.hfunc2_resp;
  u_e.col(1) = fit2_.hfunc2_resp;

  std::cout << fit1_.hfunc2_resp[0] << ", " << fit2_.hfunc2_resp[0] << "\n";

  std::vector<std::string> var_types;
  var_types.push_back("c");
  var_types.push_back("c");

  Bicop pc_final;

  pc_final.set_var_types(var_types);
  pc_final.select(u_e, controls_);
  pcs_.push_back(std::vector<Bicop>{ pc_final });
}

inline void YVineRegSelector::extend_fit(DVineFitTemporaries& fit1,
                                         DVineFitTemporaries& fit2,
                                         size_t var) const
{
  this->initialize_var(fit1, fit2, var);

  for (size_t t = 0; t < fit1.selected_vars.size() + 1; t++) {
    if (!(t == fit1.selected_vars.size())) { // Pair copulas between covariates
      auto u_e = this->get_edge_data(fit1, t);
      this->fit_pair_copula(fit1, fit2, t, u_e);
      this->update_hfuncs(fit1, fit2, t, u_e);
    } else { // Pair copula between responses and new covariate
      auto u_e1 = this->get_edge_data(fit1, t);
      auto u_e2 = this->get_edge_data(fit2, t);
      this->fit_pair_copula(fit1, t, u_e1);
      this->update_hfuncs(fit1, t, u_e1);
      this->fit_pair_copula(fit2, t, u_e2);
      this->update_hfuncs(fit2, t, u_e2);
    }
  }
  this->update_status(fit1, fit2, var);
}

inline void YVineRegSelector::initialize_var(DVineFitTemporaries& fit1,
                                             DVineFitTemporaries& fit2,
                                             size_t var) const
{
  fit1.hfunc1[0] = data_.col(var);
  fit1.hfunc1_sub[0] =
    (var_types_[var] == "d") ? data_.col(p_ + 1 + var) : Eigen::VectorXd();

  fit2.hfunc1[0] = data_.col(var);
  fit2.hfunc1_sub[0] =
    (var_types_[var] == "d") ? data_.col(p_ + 1 + var) : Eigen::VectorXd();
}

// obtain variable types for the new edge in tree t
inline std::vector<std::string> YVineRegSelector::get_edge_types(
    const DVineFitTemporaries& fit, size_t t) const
{
  // the variable type can be inferred from the existence of _sub data
  std::vector<std::string> var_types(2);
  var_types[0] = fit.hfunc2_sub[t].size() ? "d" : "c";
  var_types[1] = fit.hfunc1_sub[t].size() ? "d" : "c";
  return var_types;
}

// obtain data for the new edge in tree t
inline Eigen::MatrixXd YVineRegSelector::get_edge_data(
    const DVineFitTemporaries& fit, size_t t) const
{
  // (hfunc2 has been computed in previous fit, hfunc1 in previous tree)
  Eigen::MatrixXd u_e(data_.rows(), 2);
  u_e.col(0) = fit.hfunc2[t];
  u_e.col(1) = fit.hfunc1[t];

  if (fit.hfunc2_sub[t].size() | fit.hfunc1_sub[t].size()) {
    u_e.conservativeResize(u_e.rows(), 4);
    // use dummys for _sub data if variable is not discrete
    u_e.col(2) =
      fit.hfunc2_sub[t].size() ? fit.hfunc2_sub[t] : fit.hfunc2[t];
    u_e.col(3) =
      fit.hfunc1_sub[t].size() ? fit.hfunc1_sub[t] : fit.hfunc1[t];
  }

  return u_e;
}

inline void YVineRegSelector::fit_pair_copula(DVineFitTemporaries& fit,
                                              size_t t,
                                              const Eigen::MatrixXd& u_e) const
{
  auto var_types = this->get_edge_types(fit, t);
  fit.pcs[t].set_var_types(var_types);
  fit.pcs[t].select(u_e, controls_);
}

inline void YVineRegSelector::fit_pair_copula(DVineFitTemporaries& fit1,
                                              DVineFitTemporaries& fit2,
                                              size_t t,
                                              const Eigen::MatrixXd& u_e) const
{
  auto var_types = this->get_edge_types(fit1, t);
  fit1.pcs[t].set_var_types(var_types);
  // Fit pair copula in first D-vine
  fit1.pcs[t].select(u_e, controls_);
  // Copy fitted copula into second D-vine
  fit2.pcs[t] = fit1.pcs[t];
}

inline void YVineRegSelector::update_hfunc1(DVineFitTemporaries& fit,
                                            size_t t,
                                            const Eigen::MatrixXd& u_e) const
{
  fit.hfunc2_resp = fit.pcs[t].hfunc2(u_e);
  if (p_ == t + 1) { // selection is complete
    return;
  } else {
    fit.hfunc1[t + 1] = fit.pcs[t].hfunc1(u_e);
    if (fit.hfunc1_sub[t].size()) { // second variable is discrete
      auto u_e_sub = u_e;
      u_e_sub.col(1) = u_e.col(3);
      fit.hfunc1_sub[t + 1] = fit.pcs[t].hfunc1(u_e_sub);
    } else {
      fit.hfunc1_sub[t + 1] = Eigen::VectorXd();
    }
  }
}

inline void YVineRegSelector::update_hfunc1(DVineFitTemporaries& fit1,
                                            DVineFitTemporaries& fit2,
                                            size_t t,
                                            const Eigen::MatrixXd& u_e) const
{
  if (p_ == t + 1) // selection is complete
    return;
  fit1.hfunc1[t + 1] = fit1.pcs[t].hfunc1(u_e);
  fit2.hfunc1[t + 1] = fit1.hfunc1[t + 1];
  if (fit1.hfunc1_sub[t].size()) { // second variable is discrete
    auto u_e_sub = u_e;
    u_e_sub.col(1) = u_e.col(3);
    fit1.hfunc1_sub[t + 1] = fit1.pcs[t].hfunc1(u_e_sub);
    fit2.hfunc1_sub[t + 1] = fit1.hfunc1_sub[t + 1];
  } else {
    fit1.hfunc1_sub[t + 1] = Eigen::VectorXd();
    fit2.hfunc1_sub[t + 1] = fit1.hfunc1_sub[t + 1];
  }
}

inline void YVineRegSelector::update_hfunc2(DVineFitTemporaries& fit,
                                            size_t t,
                                            const Eigen::MatrixXd& u_e) const
{
  // use unneeded space to store hfunc2 (required when we add the next var);
  // will be shifted one up at the end
  fit.hfunc2[t] = fit.pcs[t].hfunc2(u_e);
  if (u_e.cols() > 2) {
    if (fit.hfunc2_sub[t].size()) { // first variable is discrete
      auto u_e_sub = u_e;
      u_e_sub.col(0) = u_e.col(2);
      fit.hfunc2_sub[t] = fit.pcs[t].hfunc2(u_e_sub);
    } else {
      fit.hfunc2_sub[t] = Eigen::VectorXd();
    }
  }

  if (t == fit.selected_vars.size()) { // all trees have been fit
    // shift hfunc2 entries into correct tree level
    std::rotate(
      fit.hfunc2.begin(), fit.hfunc2.end() - 1, fit.hfunc2.end());
    std::rotate(
      fit.hfunc2_sub.begin(), fit.hfunc2_sub.end() - 1, fit.hfunc2_sub.end());

    // fill first tree with actual observations
    fit.hfunc2[0] = fit.hfunc1[0];
    if (fit.hfunc1_sub[0].size()) {
      fit.hfunc2_sub[0] = fit.hfunc1_sub[0];
    } else {
      fit.hfunc2_sub[0] = Eigen::VectorXd();
    }
  }
}

inline void YVineRegSelector::update_hfunc2(DVineFitTemporaries& fit1,
                                            DVineFitTemporaries& fit2,
                                            size_t t,
                                            const Eigen::MatrixXd& u_e) const
{
  // use unneeded space to store hfunc2 (required when we add the next var);
  // will be shifted one up at the end
  fit1.hfunc2[t] = fit1.pcs[t].hfunc2(u_e);
  fit2.hfunc2[t] = fit1.hfunc2[t];
  if (u_e.cols() > 2) {
    if (fit1.hfunc2_sub[t].size()) { // first variable is discrete
      auto u_e_sub = u_e;
      u_e_sub.col(0) = u_e.col(2);
      fit1.hfunc2_sub[t] = fit1.pcs[t].hfunc2(u_e_sub);
      fit2.hfunc2_sub[t] = fit1.hfunc2_sub[t];
    } else {
      fit1.hfunc2_sub[t] = Eigen::VectorXd();
      fit2.hfunc2_sub[t] = Eigen::VectorXd();
    }
  }

  if (t == fit1.selected_vars.size()) { // all trees have been fit
    // shift hfunc2 entries into correct tree level
    std::rotate(
      fit1.hfunc2.begin(), fit1.hfunc2.end() - 1, fit1.hfunc2.end());
    std::rotate(
      fit1.hfunc2_sub.begin(), fit1.hfunc2_sub.end() - 1, fit1.hfunc2_sub.end());
    std::rotate(
      fit2.hfunc2.begin(), fit2.hfunc2.end() - 1, fit2.hfunc2.end());
    std::rotate(
      fit2.hfunc2_sub.begin(), fit2.hfunc2_sub.end() - 1, fit2.hfunc2_sub.end());

    // fill first tree with actual observations
    fit1.hfunc2[0] = fit1.hfunc1[0];
    fit2.hfunc2[0] = fit1.hfunc1[0];
    if (fit1.hfunc1_sub[0].size()) {
      fit1.hfunc2_sub[0] = fit1.hfunc1_sub[0];
      fit2.hfunc2_sub[0] = fit1.hfunc2_sub[0];
    } else {
      fit1.hfunc2_sub[0] = Eigen::VectorXd();
      fit2.hfunc2_sub[0] = Eigen::VectorXd();
    }
  }
}

void YVineRegSelector::update_hfuncs(DVineFitTemporaries& fit,
                                     size_t t,
                                     const Eigen::MatrixXd& u_e) const
{
  this->update_hfunc1(fit, t, u_e);
  this->update_hfunc2(fit, t, u_e);
}

void YVineRegSelector::update_hfuncs(DVineFitTemporaries& fit1,
                                     DVineFitTemporaries& fit2,
                                     size_t t,
                                     const Eigen::MatrixXd& u_e) const
{
  this->update_hfunc1(fit1, fit2, t, u_e);
  this->update_hfunc2(fit1, fit2, t, u_e);
}


// update value of the criterion used for variable and family selection
inline void YVineRegSelector::update_selcrit(DVineFitTemporaries& fit1,
                                             DVineFitTemporaries& fit2) const
{
  if (controls_.get_selection_criterion() == "loglik") {
    fit1.crit += fit1.pcs[fit1.selected_vars.size()].get_loglik();
    fit2.crit += fit2.pcs[fit2.selected_vars.size()].get_loglik();
  }
  if (controls_.get_selection_criterion() == "aic") {
    fit1.crit -= fit1.pcs[fit1.selected_vars.size()].get_aic();
    fit2.crit -= fit2.pcs[fit2.selected_vars.size()].get_aic();
  }
  if (controls_.get_selection_criterion() == "bic") {
    fit1.crit -= fit1.pcs[fit1.selected_vars.size()].get_bic();
    fit2.crit -= fit2.pcs[fit2.selected_vars.size()].get_bic();
  }
}

// remove var from remaining variables; add to selected variables
inline void YVineRegSelector::update_vars(DVineFitTemporaries& fit1,
                                          DVineFitTemporaries& fit2,
                                          size_t var)
  const
{
  fit1.remaining_vars.erase(
    std::remove(fit1.remaining_vars.begin(), fit1.remaining_vars.end(), var));
  fit1.selected_vars.push_back(var);

  fit2.remaining_vars = fit1.remaining_vars;
  fit2.selected_vars = fit1.selected_vars;
}

inline void YVineRegSelector::update_status(DVineFitTemporaries& fit1,
                                            DVineFitTemporaries& fit2,
                                            size_t var)
  const
{
  this->update_selcrit(fit1, fit2);
  this->update_vars(fit1, fit2, var);
}

}
