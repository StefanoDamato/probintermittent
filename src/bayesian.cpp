#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
List gammaDynamic(NumericVector y, double a0, double b0, double w) {

  int n = y.size();
  NumericVector a(n);
  NumericVector b(n);

  a[0] = w*a0;
  b[0] = w*b0;
  if (n >= 2){
    for (int t = 1; t < n; t++){
      a[t] = w*a[t-1] + y[t-1];
      b[t] = w*b[t-1] + 1;
    }
  }

  return List::create(
    _["a"] = a,
    _["b"] = b
  );
}


// [[Rcpp::export]]
List betaDynamic(NumericVector y, double v, double a0, double b0, double w) {

  int n = y.size();
  NumericVector a(n);
  NumericVector b(n);

  a[0] = w*a0 + (1-w);
  b[0] = w*b0;
  if (n >= 2){
    for (int t = 1; t < n; t++){
      a[t] = w*a[t-1] + (1-w) + v;
      b[t] = w*b[t-1] + y[t-1];
    }
  }

  return List::create(
    _["a"] = a,
    _["b"] = b
  );
}
