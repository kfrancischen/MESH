#define ARMA_DONT_USE_WRAPPER
#include "Rcwa.h"
#include <armadillo>

void RCWA::initSMatrix(
  const int n,
  RCWAMatrix *S
)
{
  S->eye(n,n);
}

RCWA::RCWAMatrices RCWA::getSMatrix(
  const int startLayer,
  const int Nx,
  const int Ny,
  const RCWAVector thicknessList,
  RCWAMatrices MMatrices,
  RCWAMatrices FMatrices
){
  int numOfLayer = thicknessList.n_elem;
  int N = Nx * Ny;
  int r1 = 0, r2 = 2*N -1, r3 = 2*N, r4 = 4*N -1;

  RCWAMatrices SMatrices;
  for(size_t i = 0; i < numOfLayer; i++){
    RCWAMatrix S;
    RCWA::initSMatrix(N, &S);
    SMatrices.push_back(S);
  }

// propogating down
  for(size_t i = startLayer; i >=1; i--){
    RCWAMatrix I = solve(MMatrices[i], MMatrices[i-1]);
    RCWAMatrix leftTop = I(span(r3, r4), span(r3, r4));
    RCWAMatrix rightTop = I(span(r3, r4), span(r1, r2));
    RCWAMatrix leftBottom = I(span(r1, r2), span(r3, r4));
    RCWAMatrix rightBottom = I(span(r1, r2), span(r1, r2));
    SMatrices[i-1](span(r1, r2), span(r1, r2)) = solve(
      leftTop - FMatrices[i] * SMatrices[i](span(r1, r2), span(r3, r4)) * leftBottom,
      FMatrices[i] * SMatrices[i](span(r1, r2), span(r1,r2))
    );

    SMatrices[i-1](span(r1, r2), span(r3, r4)) = solve(
      leftTop - FMatrices[i] * SMatrices[i](span(r1, r2), span(r3, r4)) * leftBottom,
      FMatrices[i] * SMatrices[i](span(r1, r2), span(r3,r4)) * rightBottom - rightTop
    ) * FMatrices[i-1];

    SMatrices[i-1](span(r3,r4), span(r1,r2)) = SMatrices[i](span(r3,r4), span(r1,r2)) +
      SMatrices[i](span(r3,r4), span(r3,r4)) * leftBottom * SMatrices[i-1](span(r1,r2), span(r1,r2));

    SMatrices[i-1](span(r3,r4), span(r3,r4)) = SMatrices[i](span(r3,r4), span(r3,r4)) *
      (leftBottom * SMatrices[i-1](span(r1,r2), span(r3,r4)) + rightBottom * FMatrices[i-1]);
  }

// propogating up
  for(size_t i = startLayer; i < numOfLayer - 1; i++){
    RCWAMatrix I = solve(MMatrices[i], MMatrices[i-1]);
    RCWAMatrix leftTop = I(span(r1, r2), span(r1, r2));
    RCWAMatrix rightTop = I(span(r1, r2), span(r1, r2));
    RCWAMatrix leftBottom = I(span(r3, r4), span(r1, r2));
    RCWAMatrix rightBottom = I(span(r3, r4), span(r3, r4));
    SMatrices[i+1](span(r1, r2), span(r1, r2)) = solve(
      leftTop - FMatrices[i] * SMatrices[i](span(r1, r2), span(r3, r4)) * leftBottom,
      FMatrices[i] * SMatrices[i](span(r1, r2), span(r1,r2))
    );

    SMatrices[i+1](span(r1, r2), span(r3, r4)) = solve(
      leftTop - FMatrices[i] * SMatrices[i](span(r1, r2), span(r3, r4)) * leftBottom,
      FMatrices[i] * SMatrices[i](span(r1, r2), span(r3,r4)) * rightBottom - rightTop
    ) * FMatrices[i+1];

    SMatrices[i+1](span(r3,r4), span(r1,r2)) = SMatrices[i](span(r3,r4), span(r1,r2)) +
      SMatrices[i](span(r3,r4), span(r3,r4)) * leftBottom * SMatrices[i+1](span(r1,r2), span(r1,r2));

    SMatrices[i+1](span(r3,r4), span(r3,r4)) = SMatrices[i](span(r3,r4), span(r3,r4)) *
      (leftBottom * SMatrices[i+1](span(r1,r2), span(r3,r4)) + rightBottom * FMatrices[i+1]);
  }

  return SMatrices;
}


void RCWA::populateQ(
  RCWAMatrix vL,
  RCWAMatrix vR,
  RCWAMatrix *qR,
  RCWAMatrix *qL
){
  *qR = repmat(vL.st(), vR.n_rows, 1);
  *qL = repmat(vR, 1, vL.n_rows);
}

double RCWA::poyntingFlux(
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
){
  int Nx = 2 * nGx + 1, Ny = 2 * nGy + 1;
  int N = Nx * Ny;
  int r1 = 0, r2 = 2 * N -1, r3 = 2 * N, r4 = 4 * N -1;
  double Gx = 0, Gy = 0;
  dcomplex IMAG_I = dcomplex(0, 1);
  RCWAMatrix zeroPadding(N, N, fill::zeros);
  RCWAMatrix onePadding(N, N, fill::eye);

  switch (d) {
    case NO:
      break;
    case ONE:
      Gx = 2 * datum::pi / period[0];
      break;
    case TWO:
      Gx = 2 * datum::pi / period[0];
      Gy = 2 * datum::pi / period[1];
  }

  int numOfLayer = thicknessList.n_elem;
  RCWAMatrix Gx_list(1, Nx), Gy_list(1, Ny);
  for(int i = -nGx; i <= nGx; i++){
    Gx_list(0, i) = i * Gx;
  }
  for(int i = -nGy; i <= nGy; i++){
    Gy_list(0, i) = i * Gy;
  }
  RCWAMatrix Gx_mat, Gy_mat;
  populateQ(Gx_list, Gy_list, &Gx_mat, &Gy_mat);
  Gx_mat.reshape(N, 1);
  Gy_mat.reshape(N, 1);

  RCWAMatrix kxMat = diagmat(kx + Gx_mat);
  RCWAMatrix kyMat = diagmat(ky + Gy_mat);

  RCWAMatrices EMatrices(numOfLayer), HMatrices(numOfLayer), TMatrices(numOfLayer), MMatrices(numOfLayer);
  RCWAMatrices EigenValMatrices(numOfLayer), EigenVecMatrices(numOfLayer), FMatrices(numOfLayer);


  for(size_t i = 0; i < numOfLayer; i++){
    EMatrices[i] = join_vert(join_horiz(dielectricMatrix[i], zeroPadding), join_horiz(zeroPadding, dielectricMatrix[i]));
    HMatrices[i] = join_vert(join_horiz(dielectricMatrixInverse[i], zeroPadding), join_horiz(zeroPadding, dielectricMatrixInverse[i]));
  }


  RCWAMatrix KMatrix = join_vert(
    join_horiz(kxMat * kxMat, kxMat * kyMat),
    join_horiz(kyMat * kxMat, kyMat * kyMat)
  );



  for(size_t i = 0; i < numOfLayer; i++){

    TMatrices[i] = join_vert(
      join_horiz(kyMat * dielectricMatrixInverse[i] * kyMat, -kyMat * dielectricMatrixInverse[i] * kxMat),
      join_horiz(-kxMat * dielectricMatrixInverse[i] * kyMat, kyMat * dielectricMatrixInverse[i] * kyMat)
    );

    RCWAMatrix eigMatrix = EMatrices[i] * (omega * omega * onePadding - TMatrices[i]) - KMatrix;
    cx_vec eigVal;
    eig_gen(eigVal, EigenVecMatrices[i], eigMatrix);
    eigVal = sqrt(eigVal);
    eigVal = eigVal % sign(imag(eigVal));
    EigenVecMatrices[i] = diagmat(eigVal);

    if(i == 0 || i == numOfLayer - 1){
      FMatrices[i] = eye<RCWAMatrix>(2*N, 2*N);
    }
    else{
      FMatrices[i] = diagmat(exp(-IMAG_I * dcomplex(thicknessList(i),0) * eigVal));
    }

    MMatrices[i] = zeros<RCWAMatrix>(2 * N, 2 * N);
    MMatrices[i](span(r1, r2), span(r1, r2)) = (omega * eye<RCWAMatrix>(2*N, 2*N) - TMatrices[i] / omega) *
      EigenVecMatrices[i] * inv(EigenValMatrices[i]);
    MMatrices[i](span(r1, r2), span(r3, r4)) = MMatrices[i](span(r1, r2), span(r1, r2));
    MMatrices[i](span(r3, r4), span(r1,r2)) = EigenVecMatrices[i];
    MMatrices[i](span(r3, r4), span(r3,r4)) = MMatrices[i](span(r3, r4), span(r1,r2));

    MMatrices[i] = MMatrices[i] * inv(sqrt(diagmat(diagvec(MMatrices[i].t() * MMatrices[i]))));
  }

  // now computes the flux
  double flux = 0;

  for(size_t layerIdx = 0; layerIdx < targetLayer; layerIdx++){
    if(sourceList[layerIdx] == 0) continue;

    RCWAMatrices S_matrix_target = getSMatrix(targetLayer, Nx, Ny, thicknessList,
        MMatrices, FMatrices);

    RCWAMatrix grandImaginaryMatrix = join_vert(
      join_horiz(dielectricImMatrix[layerIdx], dielectricImMatrix[layerIdx]),
      join_horiz(dielectricImMatrix[layerIdx], dielectricImMatrix[layerIdx])
    );

    RCWAMatrix q_R, q_L;
    RCWAMatrix q(diagvec(EigenValMatrices[layerIdx]));
    populateQ(q, q, &q_R, &q_L);

    RCWAMatrix source = zeros<RCWAMatrix>(4*N, 3*N);
    source(span(0,N-1), span(2*N, 3*N-1)) = -kyMat * dielectricMatrixInverse[layerIdx] / omega;
    source(span(N, 2*N-1), span(2*N, 3*N-1)) = kxMat * dielectricMatrixInverse[layerIdx] / omega;
    source(span(2*N, 3*N-1), span(N, 2*N-1)) = onePadding;
    source(span(3*N, 4*N-1), span(0, N-1)) = -onePadding;

    RCWAMatrices NewFMatrices = FMatrices;
    NewFMatrices[layerIdx] = eye<RCWAMatrix>(2*N, 2*N);

    RCWAMatrices S_matrix = getSMatrix(layerIdx, Nx, Ny, thicknessList,
        MMatrices, NewFMatrices);

    RCWAMatrix targetFields = solve(MMatrices[layerIdx], source);

    RCWAMatrix P1 = solve(
      eye<RCWAMatrix>(2*N, 2*N) - S_matrix[targetLayer](span(r1, r2), span(r3, r4)) * S_matrix_target[numOfLayer-1](span(r3, r4), span(r1, r2)),
      S_matrix[targetLayer](span(r1, r2), span(r1, r2))
    );

    RCWAMatrix P2 = S_matrix_target[numOfLayer-1](span(r3, r4), span(r1, r2)) * P1;

    RCWAMatrix Q1 = eye<RCWAMatrix>(2*N, 2*N) - FMatrices[layerIdx] * S_matrix[0](span(r3, r4), span(r1, r2)) *
      FMatrices[layerIdx] * S_matrix[numOfLayer-1](span(r3, r4), span(r1, r2));

    RCWAMatrix Q2 = -FMatrices[layerIdx] * S_matrix[0](span(r3, r4), span(r1, r2));

    RCWAMatrix R = MMatrices[targetLayer] * join_vert(FMatrices[targetLayer] * P1, P2) *
      inv(Q1) * join_horiz(eye<RCWAMatrix>(2*N, 2*N), Q2);

    RCWAMatrix integralSelf, integralMutual;
    if(layerIdx == 0 || layerIdx == numOfLayer - 1){
      integralSelf = 1 / (IMAG_I * (q_L - conj(q_R)));
      integralMutual = zeros<RCWAMatrix>(2*N, 2*N);
    }
    else{
      integralSelf = (1 - exp(-IMAG_I * dcomplex(thicknessList[layerIdx], 0) * (q_L - conj(q_R)))) /
        (IMAG_I * (q_L - conj(q_R)));
      integralMutual = (exp(IMAG_I * dcomplex(thicknessList[layerIdx], 0) * conj(q_R)) - exp(-IMAG_I * dcomplex(thicknessList[layerIdx], 0) * q_R)) /
        (IMAG_I * (q_L + conj(q_R)));
    }

    RCWAMatrix integral = join_vert(
      join_horiz(integralSelf, integralMutual),
      join_horiz(integralMutual, integralSelf)
    );

    RCWAMatrix poyntingMat = (targetFields * grandImaginaryMatrix * targetFields.t()) % integral;

    poyntingMat = -R * poyntingMat * R.t();
    flux += real(trace(poyntingMat(span(r1, r2), span(r3, r4))));
  }
  return flux;

}


