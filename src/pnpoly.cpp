#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export(.pnpoly)]]
LogicalVector pnpoly(NumericVector testx, NumericVector testy, NumericVector vertx, NumericVector verty)
{
  int i, j, k;
  int temp = 0;
  LogicalVector c(testx.size());
  
  for(k = 0; k < testx.size(); k++) {
    temp = 0;
    for (i = 0, j = vertx.size() - 1; i < vertx.size(); j = i++) {
      if ( ((verty[i] > testy[k]) != (verty[j] > testy[k])) &&
           (testx[k] < (vertx[j] - vertx[i]) * (testy[k] - verty[i]) / (verty[j] - verty[i]) + vertx[i]) )
        temp = !temp;
    }
    c[k] = temp;
  }
  
  return c;
}
