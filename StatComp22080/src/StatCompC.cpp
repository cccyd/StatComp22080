#include <Rcpp.h>
using namespace Rcpp;

//' @title A accuracy computation function
//' @description A accuracy computation function using Rcpp
//' @param N the number of samples
//' @param a matrix a
//' @param b matrix b
//' @return the accuracy ratio
//' @import Rcpp
//' @examples
//' \dontrun{
//' acc(3,as.matrix(c(1,2,3)),as.matrix(c(2,3,4)))
//' }
//' @export
// [[Rcpp::export]]
float acc(int N, NumericMatrix a, NumericMatrix b) {
  float result=0;
  int x=0, y=0;
  float ratio = 0;
  for(int i = 1; i < N; i++) {
    y = a[i];
    x = b[i];
    if (x==y){
      result ++;
    }
  }
  ratio = result / N;
  return(ratio);
}
