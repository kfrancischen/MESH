#include "Rcwa.h"
#include <armadillo>

void initSMatrix(
  size_t n,
  arma::cx_mat* S
)
{
  S->eye(n,n);
}

void getSMatrix(
  size_t numOfLayers,
  size_t startLayer,
  size_t n,
  const arma::vec* thicknessList,
  std::vector<arma::cx_mat*> MMatrices,
  std::vector<arma::cx_mat*> fMatrices
){

}

void populateQ(
  const std::complex<double> *q,
  arma::cx_mat* qR,
  arma::cx_mat* qL
){

}
