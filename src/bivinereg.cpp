#ifndef BOOST_NO_AUTO_PTR
#define BOOST_NO_AUTO_PTR
#endif

#ifndef BOOST_MATH_PROMOTE_DOUBLE_POLICY
#define BOOST_MATH_PROMOTE_DOUBLE_POLICY false
#else
#undef BOOST_MATH_PROMOTE_DOUBLE_POLICY
#define BOOST_MATH_PROMOTE_DOUBLE_POLICY false
#endif

#ifndef BOOST_ALLOW_DEPRECATED_HEADERS
#define BOOST_ALLOW_DEPRECATED_HEADERS
#endif

#include "yvine_reg_selector.hpp"
#include <RcppThread.h>
#include <vinecopulib-wrappers.hpp>
#include <vinecopulib/bicop/fit_controls.hpp>
#include <vinecopulib/misc/triangular_array.hpp>

using namespace vinecopulib;

// [[Rcpp::export]]
Rcpp::List
select_yvine_cpp(const Eigen::MatrixXd& data,
                 std::vector<std::string> family_set,
                 std::string par_method,
                 std::string nonpar_method,
                 double mult,
                 std::string selcrit,
                 const Eigen::VectorXd& weights,
                 double psi0,
                 bool preselect_families,
                 size_t cores,
                 const std::vector<std::string>& var_types)
{
  // set up the cpp fit controls from all the arguments ------
  std::vector<BicopFamily> fam_set(family_set.size());
  for (unsigned int fam = 0; fam < fam_set.size(); ++fam) {
    fam_set[fam] = to_cpp_family(family_set[fam]);
  }
  FitControlsBicop controls(fam_set,
                            par_method,
                            nonpar_method,
                            mult,
                            selcrit,
                            weights,
                            psi0,
                            preselect_families,
                            cores);

  // select the model -----------------------------------------
  bivinereg::YVineRegSelector selector(data, var_types, controls);
  selector.select_model();
  auto selected_vars = selector.get_selected_vars();
  auto pcs = selector.get_pcs();

  // make results R-compatible -------------------------------
  Rcpp::List vinecop_r;
  // rank ensures that vars are 1, ..., p_sel
  std::vector<size_t> response_vars = {0, 1};
  // auto order = tools_stl::cat(static_cast<size_t>(0), selected_vars);
  auto order = tools_stl::cat(response_vars, selected_vars);
  // Y-vine structure array

  auto rank_order = tools_stl::rank(order);

  std::vector<std::vector<size_t>> rows;

  for (size_t t = 0; t < order.size() - 2; t++) {
    std::vector<size_t> row = {rank_order[t + 2]};
    row.insert(row.end(), rank_order.begin() + 2 + t, rank_order.end());
    rows.push_back(row);
  }

  rows.push_back(std::vector<size_t> {2});

  TriangularArray<size_t> t_array(rows);

  auto new_struct = RVineStructure(rank_order, t_array, false, true); // Need RVineStructure !!!

  auto sv = selected_vars;
  std::sort(sv.begin(), sv.end());
  auto vt = std::vector<std::string>();
  vt.push_back(var_types[0]);
  vt.push_back(var_types[1]);
  for (auto v : sv)
    vt.push_back(var_types[v]);
  vinecop_r = Rcpp::List::create(
    Rcpp::Named("pair_copulas") = pair_copulas_wrap(pcs, order.size(), true),
    Rcpp::Named("structure") = rvine_structure_wrap(new_struct),
    Rcpp::Named("var_types") = vt,
    Rcpp::Named("npars") = Vinecop(new_struct, pcs).get_npars(),
    Rcpp::Named("nobs") = data.rows(),
    Rcpp::Named("loglik") = NAN,
    Rcpp::Named("threshold") = 0);
  vinecop_r.attr("class") =
    Rcpp::CharacterVector{ "vinecop", "vinecop_dist" };
  for (auto& v : selected_vars) // R indexing starts at 1
    v++;

  return Rcpp::List::create(Rcpp::Named("vine") = vinecop_r,
                            Rcpp::Named("selected_vars") = selected_vars);
}

// [[Rcpp::export]]
Eigen::VectorXd
cond_m_dist_cpp(const Eigen::MatrixXd& u,
                const Rcpp::List& vinecop_r,
                const int margin,
                size_t num_threads)
{
  tools_eigen::check_if_in_unit_cube(u);
  auto vinecop_cpp = vinecop_wrap(vinecop_r);
  auto vine_struct_ = vinecop_cpp.get_rvine_structure();
  auto d = vine_struct_.get_dim();
  if ((static_cast<size_t>(u.cols()) != d) &&
      (static_cast<size_t>(u.cols()) != 2 * d))
    throw std::runtime_error("data dimension is incompatible with model.");

  auto trunc_lvl = vine_struct_.get_trunc_lvl();
  auto order = vine_struct_.get_order();
  auto inverse_order = tools_stl::invert_permutation(order);

  Eigen::VectorXd p(u.rows());
  auto do_batch = [&](const tools_batch::Batch& b) {
    Eigen::MatrixXd hfunc1, hfunc2, u_e;
    hfunc1 = Eigen::MatrixXd::Zero(b.size, d);
    hfunc2 = Eigen::MatrixXd::Zero(b.size, d);

    // data have to be reordered to correspond to natural order
    for (size_t j = 0; j < d; ++j) {
      hfunc2.col(j) = u.block(b.begin, order[j] - 1, b.size, 1);
    }

    for (size_t tree = 0; tree < trunc_lvl; ++tree) {
      for (size_t edge = 0; edge < d - tree - 1; ++edge) {
        tools_interface::check_user_interrupt();
        Bicop edge_copula = vinecop_cpp.get_pair_copula(tree, edge);
        size_t m = vine_struct_.min_array(tree, edge);

        u_e = Eigen::MatrixXd(b.size, 2);
        u_e.col(0) = hfunc2.col(edge);
        if (m == vine_struct_.struct_array(tree, edge, true)) {
          u_e.col(1) = hfunc2.col(m - 1);
        } else {
          u_e.col(1) = hfunc1.col(m - 1);
        }

        // hfunc1 only needed in last tree
        hfunc1.col(edge) = edge_copula.hfunc1(u_e);
        hfunc2.col(edge) = edge_copula.hfunc2(u_e);
      }
    }
    if (margin == 1) {
      p.segment(b.begin, b.size) = hfunc2.col(0);
    } else {
      p.segment(b.begin, b.size) = hfunc1.col(0);
    }

  };

  RcppThread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);
  pool.map(do_batch, tools_batch::create_batches(u.rows(), num_threads));
  pool.join();

  return p;
}

// [[Rcpp::export]]
Eigen::VectorXd
cond_bi_dens_cpp(const Eigen::MatrixXd& u,
                 const Rcpp::List& vinecop_r,
                 size_t num_threads)
{
  tools_eigen::check_if_in_unit_cube(u);
  auto vinecop_cpp = vinecop_wrap(vinecop_r);
  auto rvine_structure_ = vinecop_cpp.get_rvine_structure();
  auto d_ = rvine_structure_.get_dim();
  auto pair_copulas_ = vinecop_cpp.get_all_pair_copulas();
  if ((static_cast<size_t>(u.cols()) != d_) &&
      (static_cast<size_t>(u.cols()) != 2 * d_))
    throw std::runtime_error("data dimension is incompatible with model.");

  // info about the vine structure (reverse rows (!) for more natural
  // indexing)
  size_t trunc_lvl = rvine_structure_.get_trunc_lvl();
  auto order = rvine_structure_.get_order();

  // initial value must be 1.0 for multiplication
  Eigen::VectorXd pdf = Eigen::VectorXd::Constant(u.rows(), 1.0);

  auto do_batch = [&](const tools_batch::Batch& b) {
    // temporary storage objects (all data must be in (0, 1))
    Eigen::MatrixXd hfunc1, hfunc2, u_e;
    hfunc1 = Eigen::MatrixXd::Zero(b.size, d_);
    hfunc2 = Eigen::MatrixXd::Zero(b.size, d_);

    // fill first row of hfunc2 matrix with evaluation points;
    // points have to be reordered to correspond to natural order
    for (size_t j = 0; j < d_; ++j) {
      hfunc2.col(j) = u.block(b.begin, order[j] - 1, b.size, 1);
    }

    for (size_t tree = 0; tree < trunc_lvl; ++tree) {
      tools_interface::check_user_interrupt(
        static_cast<double>(u.rows()) * static_cast<double>(d_) > 1e5);
      for (size_t edge = 0; edge < d_ - tree - 1; ++edge) {
        tools_interface::check_user_interrupt(edge % 100 == 0);
        // extract evaluation point from hfunction matrices (have been
        // computed in previous tree level)
        Bicop* edge_copula = &pair_copulas_[tree][edge];
        size_t m = rvine_structure_.min_array(tree, edge);

        u_e = Eigen::MatrixXd(b.size, 2);
        u_e.col(0) = hfunc2.col(edge);
        if (m == rvine_structure_.struct_array(tree, edge, true)) {
          u_e.col(1) = hfunc2.col(m - 1);
        } else {
          u_e.col(1) = hfunc1.col(m - 1);
        }

        if (edge == 0 || edge == 1) {
          pdf.segment(b.begin, b.size) =
            pdf.segment(b.begin, b.size)
               .cwiseProduct(edge_copula->pdf(u_e));
        }

        // h-functions are only evaluated if needed in next step
        if (rvine_structure_.needed_hfunc1(tree, edge)) {
          hfunc1.col(edge) = edge_copula->hfunc1(u_e);
        }
        if (rvine_structure_.needed_hfunc2(tree, edge)) {
          hfunc2.col(edge) = edge_copula->hfunc2(u_e);
        }
      }
    }
  };

  if (trunc_lvl > 0) {
    tools_thread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);
    pool.map(do_batch, tools_batch::create_batches(u.rows(), num_threads));
    pool.join();
  }

  return pdf;
}

// [[Rcpp::export]]
Eigen::VectorXd
cond_m_dens_cpp(const Eigen::MatrixXd& u,
                const Rcpp::List& vinecop_r,
                const int margin,
                size_t num_threads)
{
  tools_eigen::check_if_in_unit_cube(u);
  auto vinecop_cpp = vinecop_wrap(vinecop_r);
  auto rvine_structure_ = vinecop_cpp.get_rvine_structure();
  auto d_ = rvine_structure_.get_dim();
  auto pair_copulas_ = vinecop_cpp.get_all_pair_copulas();
  if ((static_cast<size_t>(u.cols()) != d_) &&
      (static_cast<size_t>(u.cols()) != 2 * d_))
    throw std::runtime_error("data dimension is incompatible with model.");

  // info about the vine structure (reverse rows (!) for more natural
  // indexing)
  size_t trunc_lvl = rvine_structure_.get_trunc_lvl();
  auto order = rvine_structure_.get_order();

  // initial value must be 1.0 for multiplication
  Eigen::VectorXd pdf = Eigen::VectorXd::Constant(u.rows(), 1.0);

  auto do_batch = [&](const tools_batch::Batch& b) {
    // temporary storage objects (all data must be in (0, 1))
    Eigen::MatrixXd hfunc1, hfunc2, u_e;
    hfunc1 = Eigen::MatrixXd::Zero(b.size, d_);
    hfunc2 = Eigen::MatrixXd::Zero(b.size, d_);

    // fill first row of hfunc2 matrix with evaluation points;
    // points have to be reordered to correspond to natural order
    for (size_t j = 0; j < d_; ++j) {
      hfunc2.col(j) = u.block(b.begin, order[j] - 1, b.size, 1);
    }

    for (size_t tree = 0; tree < trunc_lvl - 1; ++tree) {
      tools_interface::check_user_interrupt(
        static_cast<double>(u.rows()) * static_cast<double>(d_) > 1e5);
      for (size_t edge = 0; edge < d_ - tree - 1; ++edge) {
        tools_interface::check_user_interrupt(edge % 100 == 0);
        // extract evaluation point from hfunction matrices (have been
        // computed in previous tree level)
        Bicop* edge_copula = &pair_copulas_[tree][edge];
        size_t m = rvine_structure_.min_array(tree, edge);

        u_e = Eigen::MatrixXd(b.size, 2);
        u_e.col(0) = hfunc2.col(edge);
        if (m == rvine_structure_.struct_array(tree, edge, true)) {
          u_e.col(1) = hfunc2.col(m - 1);
        } else {
          u_e.col(1) = hfunc1.col(m - 1);
        }

        if (edge == margin) {
          pdf.segment(b.begin, b.size) =
            pdf.segment(b.begin, b.size)
               .cwiseProduct(edge_copula->pdf(u_e));
        }

        // h-functions are only evaluated if needed in next step
        if (rvine_structure_.needed_hfunc1(tree, edge)) {
          hfunc1.col(edge) = edge_copula->hfunc1(u_e);
        }
        if (rvine_structure_.needed_hfunc2(tree, edge)) {
          hfunc2.col(edge) = edge_copula->hfunc2(u_e);
        }
      }
    }
  };

  if (trunc_lvl > 0) {
    tools_thread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);
    pool.map(do_batch, tools_batch::create_batches(u.rows(), num_threads));
    pool.join();
  }

  return pdf;
}
