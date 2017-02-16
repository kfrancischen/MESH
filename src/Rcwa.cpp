/* Copyright (C) 2016-2018, Stanford University
 * This file is part of MESH
 * Written by Kaifeng Chen (kfchen@stanford.edu)
 *
 * MESH is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * MESH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
#include "Rcwa.h"

/*============================================================
* Function similar to meshgrid in matlab
@arg:
 vL, vR: the input two vectors (here everything is a matrix)
 qL, qR: the input two vectors
 [qL, qR] = meshgrid(vl, vR)
==============================================================*/
void RCWA::meshGrid(
  const RCWAMatrix *vL,
  const RCWAMatrix *vR,
  RCWAMatrix *qL,
  RCWAMatrix *qR
){
  *qL = repmat((*vL).st(), vR->n_rows, 1);
  *qR = repmat(*vR, 1, vL->n_rows);
}

/*============================================================
* Function computing S matrix for each layers
@arg:
 startLayer: the starting layer for the propogation
 N: number of G total G
 numOfLayer: the number of layers
 MMatrices: matrix corresponding to the propogation in each layer
 FMatrices: matrix corresponding to the phase in each layer
 SMatrices: the S matrix for each layer
==============================================================*/
void RCWA::getSMatrices(
  const int startLayer,
  const int N,
  const int numOfLayer,
  RCWAMatrices* MMatrices,
  RCWAMatrices* FMatrices,
  RCWAMatrices* SMatrices
){
  int r1 = 0, r2 = 2*N -1, r3 = 2*N, r4 = 4*N -1;
// propogating down

  for(size_t i = startLayer; i >=1; i--){

    RCWAMatrix I = solve((*MMatrices)[i], (*MMatrices)[i-1]);
    RCWAMatrix leftTop = I(span(r3, r4), span(r3, r4));
    RCWAMatrix rightTop = I(span(r3, r4), span(r1, r2));
    RCWAMatrix leftBottom = I(span(r1, r2), span(r3, r4));
    RCWAMatrix rightBottom = I(span(r1, r2), span(r1, r2));

    (*SMatrices)[i-1](span(r1, r2), span(r1, r2)) = solve(
      leftTop - (*FMatrices)[i] * (*SMatrices)[i](span(r1, r2), span(r3, r4)) * leftBottom,
      (*FMatrices)[i] * (*SMatrices)[i](span(r1, r2), span(r1,r2))
    );

    (*SMatrices)[i-1](span(r1, r2), span(r3, r4)) = solve(
      leftTop - (*FMatrices)[i] * (*SMatrices)[i](span(r1, r2), span(r3, r4)) * leftBottom,
      (*FMatrices)[i] * (*SMatrices)[i](span(r1, r2), span(r3,r4)) * rightBottom - rightTop
    ) * (*FMatrices)[i-1];

    (*SMatrices)[i-1](span(r3,r4), span(r1,r2)) = (*SMatrices)[i](span(r3,r4), span(r1,r2)) +
      (*SMatrices)[i](span(r3,r4), span(r3,r4)) * leftBottom * (*SMatrices)[i-1](span(r1,r2), span(r1,r2));

    (*SMatrices)[i-1](span(r3,r4), span(r3,r4)) = (*SMatrices)[i](span(r3,r4), span(r3,r4)) *
      (leftBottom * (*SMatrices)[i-1](span(r1,r2), span(r3,r4)) + rightBottom * (*FMatrices)[i-1]);
  }
// propogating up
  for(size_t i = startLayer; i < numOfLayer - 1; i++){
    RCWAMatrix I = solve((*MMatrices)[i], (*MMatrices)[i+1]);

    RCWAMatrix leftTop = I(span(r1, r2), span(r1, r2));
    RCWAMatrix rightTop = I(span(r1, r2), span(r3, r4));
    RCWAMatrix leftBottom = I(span(r3, r4), span(r1, r2));
    RCWAMatrix rightBottom = I(span(r3, r4), span(r3, r4));

    (*SMatrices)[i+1](span(r1, r2), span(r1, r2)) = solve(
      leftTop - (*FMatrices)[i] * (*SMatrices)[i](span(r1, r2), span(r3, r4)) * leftBottom,
      (*FMatrices)[i] * (*SMatrices)[i](span(r1, r2), span(r1,r2))
    );

    (*SMatrices)[i+1](span(r1, r2), span(r3, r4)) = solve(
      leftTop - (*FMatrices)[i] * (*SMatrices)[i](span(r1, r2), span(r3, r4)) * leftBottom,
      (*FMatrices)[i] * (*SMatrices)[i](span(r1, r2), span(r3,r4)) * rightBottom - rightTop
    ) * (*FMatrices)[i+1];

    (*SMatrices)[i+1](span(r3,r4), span(r1,r2)) = (*SMatrices)[i](span(r3,r4), span(r1,r2)) +
      (*SMatrices)[i](span(r3,r4), span(r3,r4)) * leftBottom * (*SMatrices)[i+1](span(r1,r2), span(r1,r2));

    (*SMatrices)[i+1](span(r3,r4), span(r3,r4)) = (*SMatrices)[i](span(r3,r4), span(r3,r4)) *
      (leftBottom * (*SMatrices)[i+1](span(r1,r2), span(r3,r4)) + rightBottom * (*FMatrices)[i+1]);
  }
}

/*============================================================
* Function computing numbef of G for the system
@arg:
 nGx: positive G along x direction
 nGy: positive G along y direction
==============================================================*/
int RCWA::getN(
 const int nGx,
 const int nGy
){
  return (2*nGx + 1) * (2*nGy + 1);
}
/*============================================================
* Function computing the sinc function (sin(x) / x)
@arg:
 x: the input argument
 mask: the matrix denoting the zero positions of x
==============================================================*/
RCWA::RCWAMatrix RCWA::sinc(RCWAMatrix x, RCWAMatrix* mask){
  RCWAMatrix output = sin(x * datum::pi) / (x * datum::pi);
  if(mask != nullptr){
    RCWAMatrix allOnes = RCWAMatrix(size(*mask), fill::ones);
    output = output % (allOnes - *mask) + *mask;
  }
  return output;
}
/*============================================================
* Function computing G matrix for the system
@arg:
 startLayer: the starting layer for the propogation
 nGx: positive G along x direction
 nGy: positive G along y direction
 period: periodicity along two directions
 Gx_mat: the Gx matrix
 Gy_mat: the Gy matrix
 d: the dimension of the structure
==============================================================*/
void RCWA::getGMatrices(
  const int nGx,
  const int nGy,
  const double period[2],
  RCWAMatrix* Gx_mat,
  RCWAMatrix* Gy_mat,
  DIMENSION d
)
{

  int Nx = 2 * nGx + 1;
  int Ny = 2 * nGy + 1;
  int N = Nx * Ny;
  int Gx = 0, Gy = 0;
  switch (d) {
    case NO_:
      break;
    case ONE_:{
      Gx = 2 * datum::pi / period[0];
      break;
    }
    case TWO_:{
      Gx = 2 * datum::pi / period[0];
      Gy = 2 * datum::pi / period[1];
      break;
    }
  }
  RCWAMatrix Gx_list(1, Nx), Gy_list(1, Ny);
  for(int i = -nGx; i <= nGx; i++){
    Gx_list(0, i+nGx) = i * Gx;
  }
  for(int i = -nGy; i <= nGy; i++){
    Gy_list(0, i+nGy) = i * Gy;
  }

  meshGrid(&Gx_list, &Gy_list, Gx_mat, Gy_mat);

  Gx_mat->reshape(N, 1);
  Gy_mat->reshape(N, 1);
}

/*============================================================
* Function computing imaginary dielectric matrix for the system
@arg:
grandImaginaryMatrices: the matrices containing the imaginary part matrix
dielectricImMatrix: the imaginary part matrices
numOfLayer: the number of layer in the system
N: the number of G
==============================================================*/
void RCWA::getGrandImaginaryMatrices(
  RCWAMatrices* grandImaginaryMatrices,
  RCWAMatrices* dielectricImMatrix,
  int numOfLayer,
  int N
)
{
  for(size_t i = 0; i < numOfLayer; i++){
    RCWAMatrix grandImaginaryMatrix = zeros<RCWAMatrix>(3*N, 3*N);
    grandImaginaryMatrix(span(0, N-1), span(0, N-1)) = (*dielectricImMatrix)[i];
    grandImaginaryMatrix(span(N, 2*N-1), span(N, 2*N-1)) = (*dielectricImMatrix)[i];
    grandImaginaryMatrix(span(2*N, 3*N-1), span(2*N, 3*N-1)) = (*dielectricImMatrix)[i];
    grandImaginaryMatrices->push_back(grandImaginaryMatrix);
  }
}

/*============================================================
* Function computing E matrices for the system
@arg:
EMatrices: the Ematrices
dielectricMatrix: the dielectric matrices
numOfLayer: the number of layer in the system
N: the number of G
==============================================================*/
void RCWA::getEMatrices(
  RCWAMatrices* EMatrices,
  RCWAMatrices* dielectricMatrix,
  int numOfLayer,
  int N
){
  RCWAMatrix zeroPadding(N, N, fill::zeros);
  for(size_t i = 0; i < numOfLayer; i++){
    RCWAMatrix EMatrix = join_vert(
      join_horiz((*dielectricMatrix)[i], zeroPadding),
      join_horiz(zeroPadding, (*dielectricMatrix)[i])
    );
    EMatrices->push_back(EMatrix);
  }
}


/*============================================================
* Function computing the poynting vector at given (kx, ky)
@arg:
 omega: the angular frequency (normalized to c)
 thicknessList: the thickness for each layer
 kx: the k vector at x direction (normalized value)
 ky: the y vector at x direction (normalized value)
 EMatrices:  the E matrices for all layers
 grandImaginaryMatrices: collection of all imaginary matrices in all layers
 dielectricMatrixInverse: the inverse of dielectric matrix
 Gx_mat: the Gx matrix
 Gy_mat: the Gy matrix
 sourceList: list of 0 or 1 with the same size of thicknessList
 targetLayer: the targetLayer for the flux measurement
 N: total number of G
==============================================================*/
double RCWA::poyntingFlux(
  const double omega,
  const RCWAVector* thicknessList,
  double kx,
  double ky,
  const RCWAMatrices* EMatrices,
  const RCWAMatrices* grandImaginaryMatrices,
  const RCWAMatrices* dielectricMatrixInverse,
  const RCWAMatrix* Gx_mat,
  const RCWAMatrix* Gy_mat,
  const SourceList* sourceList,
  const int targetLayer,
  const int N
){

  /*======================================================
  this part initializes parameters
  =======================================================*/
  kx = kx * omega;
  ky = ky * omega;
  int r1 = 0, r2 = 2 * N -1, r3 = 2 * N, r4 = 4 * N -1;
  dcomplex IMAG_I = dcomplex(0, 1);
  RCWAMatrix onePadding4N(4*N, 4*N, fill::eye);
  RCWAMatrix onePadding2N(2*N, 2*N, fill::eye);
  RCWAMatrix onePadding1N(N, N, fill::eye);
  RCWAMatrix zeroPadding2N(2*N, 2*N, fill::zeros);
  RCWAMatrix zeroPadding4N(4*N, 4*N, fill::zeros);
  int numOfLayer = thicknessList->n_elem;

  // populate Gx and Gy matrices
  RCWAMatrix kxMat = diagmat(kx + *Gx_mat);
  RCWAMatrix kyMat = diagmat(ky + *Gy_mat);
  /*======================================================
  this part initializes structure matrices
  =======================================================*/
  RCWAMatrices TMatrices(numOfLayer), MMatrices(numOfLayer);
  RCWAMatrices EigenValMatrices(numOfLayer), EigenVecMatrices(numOfLayer), FMatrices(numOfLayer);


  // initialize K matrix
  RCWAMatrix KMatrix = join_vert(
    join_horiz(kxMat * kxMat, kxMat * kyMat),
    join_horiz(kyMat * kxMat, kyMat * kyMat)
  );
  /*======================================================
  This part solves RCWA
  e.g initialize M and F matrices, and compute the Eigen value problem
  =======================================================*/
  for(size_t i = 0; i < numOfLayer; i++){

    TMatrices[i] = join_vert(
      join_horiz(kyMat * (*dielectricMatrixInverse)[i] * kyMat, -kyMat * (*dielectricMatrixInverse)[i] * kxMat),
      join_horiz(-kxMat * (*dielectricMatrixInverse)[i] * kyMat, kxMat * (*dielectricMatrixInverse)[i] * kxMat)
    );
    RCWAMatrix eigMatrix = (*EMatrices)[i] * (POW2(omega) * onePadding2N - TMatrices[i]) - KMatrix;
    cx_vec eigVal;
    // here is the problem
    eig_gen(eigVal, EigenVecMatrices[i], eigMatrix);

    eigVal = sqrt(eigVal);
    eigVal = -eigVal % sign(imag(eigVal));
    EigenValMatrices[i] = diagmat(eigVal);
    if(i == 0 || i == numOfLayer - 1){
      FMatrices[i] = onePadding2N;
    }
    else{
      FMatrices[i] = diagmat(exp(-IMAG_I * dcomplex((*thicknessList)(i),0) * eigVal));
    }
    MMatrices[i] = zeroPadding4N;
    MMatrices[i](span(r1, r2), span(r1, r2)) = (omega * onePadding2N - TMatrices[i] / omega) *
      EigenVecMatrices[i] * (EigenValMatrices[i]).i();

    MMatrices[i](span(r1, r2), span(r3, r4)) = -MMatrices[i](span(r1, r2), span(r1, r2));
    MMatrices[i](span(r3, r4), span(r1, r2)) = EigenVecMatrices[i];
    MMatrices[i](span(r3, r4), span(r3, r4)) = MMatrices[i](span(r3, r4), span(r1, r2));
    // normalization
    MMatrices[i] = MMatrices[i] * (diagmat(sqrt(diagvec(MMatrices[i].t() * MMatrices[i])))).i();
  }

  /*======================================================
  This part initialize matrix for flux computation
  =======================================================*/

  double flux = 0;
  RCWAMatrices S_matrices_target;


  for(size_t i = 0; i < numOfLayer; i++){
    S_matrices_target.push_back(onePadding4N);
  }

  getSMatrices(targetLayer, N, numOfLayer,
      &MMatrices, &FMatrices, &S_matrices_target);

  RCWAMatrix q_R, q_L, targetFields, P1, P2, Q1, Q2, R;
  RCWAMatrix integralSelf, integralMutual, integral, poyntingMat;
  RCWAMatrices S_matrices(numOfLayer), NewFMatrices(numOfLayer);

  RCWAMatrix source = zeros<RCWAMatrix>(4*N, 3*N);
  source(span(2*N, 3*N-1), span(N, 2*N-1)) = onePadding1N;
  source(span(3*N, 4*N-1), span(0, N-1)) = -onePadding1N;
  /*======================================================
  This part compute flux by collecting emission from source layers
  =======================================================*/

  for(size_t layerIdx = 0; layerIdx < targetLayer; layerIdx++){

    // if is not source layer, then continue
    if((*sourceList)[layerIdx] == ISNOTSOURCE_) continue;

    // initial steps, propogate S matrix
    RCWAMatrix q(diagvec(EigenValMatrices[layerIdx]));
    meshGrid(&q, &q, &q_R, &q_L);

    // defining source
    source(span(0,N-1), span(2*N, 3*N-1)) = -kyMat * (*dielectricMatrixInverse)[layerIdx] / omega;
    source(span(N, 2*N-1), span(2*N, 3*N-1)) = kxMat * (*dielectricMatrixInverse)[layerIdx] / omega;
    source(span(2*N, 3*N-1), span(N, 2*N-1)) = onePadding1N;
    source(span(3*N, 4*N-1), span(0, N-1)) = -onePadding1N;

    // treat as if the source layer has no thickness
    NewFMatrices = FMatrices;
    NewFMatrices[layerIdx] = onePadding2N;

    for(size_t i = 0; i < numOfLayer; i++){
      S_matrices[i] = onePadding4N;
    }

    getSMatrices(layerIdx, N, numOfLayer,
        &MMatrices, &NewFMatrices, &S_matrices);


    // solve the source
    targetFields = solve(MMatrices[layerIdx], source);

    // calculating the P1 and P2
    P1 = solve(
      onePadding2N - S_matrices[targetLayer](span(r1, r2), span(r3, r4)) * S_matrices_target[numOfLayer-1](span(r3, r4), span(r1, r2)),
      S_matrices[targetLayer](span(r1, r2), span(r1, r2))
    );

    P2 = S_matrices_target[numOfLayer-1](span(r3, r4), span(r1, r2)) * P1;

    // calculating the Q1 and Q2
    Q1 = onePadding2N - FMatrices[layerIdx] * S_matrices[0](span(r3, r4), span(r1, r2)) *
      FMatrices[layerIdx] * S_matrices[numOfLayer-1](span(r3, r4), span(r1, r2));

    Q2 = -FMatrices[layerIdx] * S_matrices[0](span(r3, r4), span(r1, r2));

    // calculating R
    R = MMatrices[targetLayer] * join_vert(FMatrices[targetLayer] * P1, P2) *
      Q1.i() * join_horiz(onePadding2N, Q2);

    // calculating integrands
    if(layerIdx == 0 || layerIdx == numOfLayer - 1){
      integralSelf = 1 / (IMAG_I * (q_L - conj(q_R)));
      integralMutual = zeroPadding2N;
    }
    else{
      integralSelf = (1 - exp(-IMAG_I * dcomplex((*thicknessList)[layerIdx], 0) * (q_L - conj(q_R)))) /
        (IMAG_I * (q_L - conj(q_R)));
      integralMutual = (exp(IMAG_I * dcomplex((*thicknessList)[layerIdx], 0) * conj(q_R)) - exp(-IMAG_I * dcomplex((*thicknessList)[layerIdx], 0) * q_L)) /
        (IMAG_I * (q_L + conj(q_R)));
    }

    integral = join_vert(
      join_horiz(integralSelf, integralMutual),
      join_horiz(integralMutual, integralSelf)
    );

    // calculating kernel
    poyntingMat = (targetFields * (*grandImaginaryMatrices)[layerIdx] * targetFields.t()) % integral;

    poyntingMat = -R * poyntingMat * R.t();
    flux += real(trace(poyntingMat(span(r1, r2), span(r3, r4))));
  }
  return flux;

}


