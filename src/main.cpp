#include "Rcwa.h"

#include <armadillo>

int main(){
  arma::cx_mat S;
  RCWA::initSMatrix(5, &S);
  std::cout << S;
  return 0;
}
