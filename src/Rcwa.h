#ifndef _RCWA_H
#define _RCWA_H
#define ARMA_DONT_USE_WRAPPER
#include <armadillo>
#include <complex>
#include <vector>

namespace RCWA{

  using namespace arma;
  typedef cx_mat RCWAMatrix;
  typedef std::vector< RCWAMatrix > RCWAMatrices;
  typedef fvec RCWAVector;
  typedef std::vector< int > SourceList;
  typedef std::complex<double> dcomplex;

  enum DIMENSION { NO, ONE, TWO };

  void initSMatrix(const int n, RCWAMatrix* S);


  RCWAMatrices getSMatrix(
    const int startLayer,
    const int Nx,
    const int Ny,
    const RCWAVector thicknessList,
    RCWAMatrices MMatrices,
    RCWAMatrices FMatrices
  );

  void populateQ(
    RCWAMatrix vL,
    RCWAMatrix vR,
    RCWAMatrix* qR,
    RCWAMatrix* qL
  );

  double poyntingFlux(
    const double omega,
    const RCWAVector thicknessList,
    double kx,
    double ky,
    const RCWAMatrices dielectricMatrixInverse,
    const RCWAMatrices dielectricMatrix,
    const RCWAMatrices dielectricImMatrix,
    const SourceList sourceList,
    const int targetLayer,
    const int nGx,
    const int nGy,
    const double *period,
    const DIMENSION d
  );
}
#endif
