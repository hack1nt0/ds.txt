#include <Rcpp.h>
#include <vector>
#include <iostream>
using namespace Rcpp;
using namespace std;

// This is a simple function using Rcpp that creates an R list
// containing a character vector and a numeric vector.
//
// Learn more about how to use Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//
// and browse examples of code using Rcpp at:
//
//   http://gallery.rcpp.org/
//

struct A {
    A() { cout << "A()" << endl; }
};

// [[Rcpp::export]]
void rcpp_hello() {
    int nx = 3, np = 4, nc = 2;
    NumericMatrix P(nx, nc);
    NumericMatrix meanGrad(np, nc);
    NumericMatrix varGrad(np, nc);
    NumericMatrix maxVarGrad(np, nc);
    cout << P << endl;
}
