#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int do_it(const int i) {
  return i * i;
}
