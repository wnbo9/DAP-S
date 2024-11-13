#include <Rcpp.h>
using namespace Rcpp;

//' Square input using C++
//' @param x Input
//' @return Result
//' @export
// [[Rcpp::export]]
double square(double x) {
  return x * x;
}
